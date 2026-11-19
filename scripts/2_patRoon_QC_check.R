#------------------------------------------------------------------------------#
# NTS to risk assessment - feature finding using patRoon and prioritization    #
# QC check of spiked samples                                                   #
#                                                                              #
# HILIC only                                                                   #
#                                                                              #
# Developed by Nienke Meekel, KWR Water Research Institute                     #
#                                                                              #
# R version 4.5.2                                                              #
# patRoon version 2.3.4                                                        #
# executed on VM CWG                                                           #
#                                                                              #
# Current version: 23 March 2026                                               #
# 404300/028 (From NTS to risk assessment)                                     #
# ---------------------------------------------------------------------------- #

# allocate additional memory
options(java.parameters = "-Xmx8000m")

# load required packages
library(patRoon) # for NTS feature finding and processing
library(dplyr) # for data handling
library(ggplot2) # for visualization

# check dependencies
verifyDependencies() # check if all dependencies are found

# set seed
seed <- 234
set.seed(seed)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Initialization ----
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

## Obtain analysis info
# read file with analysis info
anaInfo <- readxl::read_excel(path = "data\\analysisinfo.xlsx", sheet = "hilic_anonimized", col_names = TRUE)
anaInfo <- anaInfo[-c(1, 58), ] # skip the first and last row as these blanks are for stabilizing

# Load features and fGroups data
fList <- readRDS(file = "output\\data\\fList_hilic_original.RDS")
fGroups <- readRDS(file = "output\\data\\fGroups_hilic_original.RDS")


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Features - filtering ----
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# apply general filters (remove features outside retention range, with lower intensity than 10,000, and occurring in the blank)
fGroups <- patRoon::filter(fGroups,
                           preAbsMinIntensity = 100,
                           retentionRange = c(90, 1200), # RT between 2.4 and 20 minutes
                           mzRange = NULL,
                           blankThreshold = 5, # feature groups that are also present in blank analyses are filtered out unless their relative intensity is above this threshold
                           absMinIntensity = 10000, # minimum absolute intensity for features to be kept.
                           removeBlanks = TRUE
)

# apply replicate abundance filter
fGroups <- patRoon::filter(fGroups, relMinReplicateAbundance = 1)

# apply CV (or max replicate relative standard deviation of the intensity) filter (50%)
fGroups <- patRoon::filter(fGroups, maxReplicateIntRSD = 0.5)

# Calculate peak qualities
fGroups <- calculatePeakQualities(fGroups) # calculate peak qualities for all features and groups

# Remove feature groups with a GaussianSimilarity < 0.5 and totalScore < 5
fGroups <- patRoon::filter(fGroups, featQualityRange = list(GaussianSimilarity = c(0.5, Inf), totalScore = c(5, Inf)))

# Select QC standards
fGroups <- patRoon::filter(fGroups, rGroups = c("S0", "S1", "S2", "S3", "S4", "S5", "sw_addition", "dw_addition"))


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# QC spiked compounds ----
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# load suspect list with spiked compounds
suspects <- read.csv(file = "data\\libraries\\suspectlist_standards_hilic.csv", sep = ";")

# screen suspects and remove other feature groups that do not match
fGroupsQC <- screenSuspects(fGroups, suspects, onlyHits = TRUE, rtWindow = 6, mzWindow = 0.0015)

saveRDS(object = fGroupsQC, file = paste0("output\\QC\\", Sys.Date(), "_fGroups_hilic_QC.RDS"))


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# MS1 & MS2 data ----
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# Retrieve MS peak lists
set.seed(seed)
avgMSListParams <- getDefAvgPListParams(clusterMzWindow = 0.001)
mslists <- generateMSPeakLists(
  fGroups = fGroupsQC,
  algorithm = "mzr",
  maxMSRtWindow = 5,
  precursorMzWindow = 0.8, # +/- 0.8 Da = 1.6 Da isolation width
  avgFeatParams = avgMSListParams,
  avgFGroupParams = avgMSListParams,
  topMost = NULL
)

# Rule based filtering of MS peak lists (for screening MassBank etc.)
mslistsF <- patRoon::filter(
  mslists,
  withMSMS = TRUE,
  absMSIntThr = NULL,
  absMSMSIntThr = NULL,
  relMSIntThr = NULL,
  relMSMSIntThr = NULL,
  topMSPeaks = NULL,
  topMSMSPeaks = 25,
  deIsotopeMS = FALSE,
  deIsotopeMSMS = FALSE
)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Combine results ----
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# create data table with peak areas
result <- patRoon::as.data.table(fGroupsQC, area = TRUE)

# add column on whether MS2 was available or not
fGroups.ms2 <- names(which(lengths(mslists@averagedPeakLists) == 2)) # get names of feature groups with MS2 data
result$MS2 <- "no"
result$MS2[which(result$group %in% fGroups.ms2)] <- "yes"

# export results overview
xlsx::write.xlsx(
  x = as.data.frame(result), file = paste0("output\\QC\\", Sys.Date(), "_QC_results_hilic.xlsx"), sheetName = "general",
  col.names = TRUE, row.names = FALSE, append = FALSE
)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Plot chromatograms ----
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# load names of fGroups of which the peak shape needs to be inspected
peakshapes <- result$group
  
# plot chromatograms of the selected fGroups
for (i in peakshapes) {
  png(
    file = paste0("output\\QC\\plots\\chromatograms_20260323\\fGroup_", i, ".png"),
    width = 15, height = 15, units = "cm", res = 330
  )
  plotChroms(fGroupsQC[, i],
    retMin = TRUE,
    showPeakArea = TRUE,
    showFGroupRect = FALSE,
    colourBy = "rGroups",
    EICParams = getDefEICParams(topMost = NULL),
    showLegend = TRUE
  )
  dev.off()
  print(i)
}

# remove variables we don't need anymore
rm(peakshapes, i)
