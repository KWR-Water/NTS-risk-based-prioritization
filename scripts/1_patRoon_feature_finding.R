#------------------------------------------------------------------------------#
# NTS to risk assessment - feature finding using patRoon                       #
#                                                                              #
# HILIC only                                                                   #
#                                                                              #
# Developed by Nienke Meekel, KWR Water Research Institute                     #
#                                                                              #
# R version 4.5.2                                                              #
# patRoon version 2.3.4                                                        #
#                                                                              #
# Current version: 16 March 2026                                               #
# 404300/028 (From NTS to risk assessment)                                     #
# ---------------------------------------------------------------------------- #

# allocate additional memory
options(java.parameters = "-Xmx8000m")

# load required packages
library(patRoon) # for NTS feature finding and processing

# check dependencies
verifyDependencies() # check if all dependencies are found

# set seed
seed <- 234
set.seed(seed)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Initialization ----
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

## Convert .raw to mzML
# convert HILIC files from .raw to .mzML (only run once, set to 'TRUE')
do.data.pretreatment <- FALSE
if (do.data.pretreatment) {
  convertMSFiles(
    files = "",
    outPath = "data\\HILIC",
    algorithm = "pwiz",
    dirs = FALSE,
    from = "thermo",
    to = "mzML"
  )
}

## Obtain analysis info
# read file with analysis info
anaInfo <- readxl::read_excel(path = "data\\analysisinfo.xlsx", sheet = "hilic_anonimized", col_names = TRUE)
anaInfo <- anaInfo[-c(1, 58), ] # skip the first and last row as these blanks are for stabilizing


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Feature parameter optimization ----
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

do.feature.parameter.optimization <- FALSE

if (do.feature.parameter.optimization) {
  # get standards and addition files for feature parameter optimization
  indices <- which(anaInfo$group %in% c("S1", "S2", "S3", "S4", "S5", "sw_addition", "dw_addition"))
  anaInfo.fpo <- anaInfo[indices, ]
  row.names(anaInfo.fpo) <- NULL
  
  # define parameter sets for the OpenMS algorithm (with rather high noise threshold)
  set.seed(seed)
  pSet <- list(
    chromFWHM = c(5, 10),
    mzPPM = c(3, 5),
    minFWHM = c(4, 12),
    maxFWHM = c(35, 48),
    noiseThrInt = 150000
  )
  
  # optimize feature finding parameters
  ftOpt <- optimizeFeatureFinding(
    anaInfo = anaInfo.fpo, algorithm = "openms", pSet,
    isoIdent = "OpenMS", maxIterations = 30
  )
  saveRDS(object = ftOpt, file = paste0("output\\", Sys.Date(), "_featureparoptimization_openms_hilic.RDS"))
  
  # visualize results and export to .png file
  optimizedParameters(ftOpt)
  png(
    file = paste0("output\\feature parameter optimization\\", Sys.Date(), "_featureparoptimization_hilic.png"), width = 30, height = 20, units = "cm",
    res = 1200, bg = "transparent"
  )
  plot(ftOpt, paramSet = 1, DoEIteration = 6)
  dev.off()
  
  # clean environment
  rm(anaInfo.fpo, ftOpt, pSet, indices)
}

rm(do.data.pretreatment, do.feature.parameter.optimization)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Features - finding & grouping ----
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# Find all features with OpenMS (using settings from DoE, 29 October 2024)
fList <- findFeatures(anaInfo,
                      "openms",
                      noiseThrInt = 1000,
                      chromSNR = 5, 
                      chromFWHM = 6.14,
                      minFWHM = 1,
                      maxFWHM = 41.5,
                      mzPPM = 3
)
saveRDS(object = fList, file = paste0("output\\data\\", Sys.Date(), "_fList_hilic_original.RDS"))

# Group features with OpenMS
fGroups <- groupFeaturesOpenMS(fList,
                               rtalign = TRUE,
                               maxAlignRT = 30, # time in seconds
                               maxAlignMZ = 0.0015, # mz in Da
                               maxGroupRT = 6,
                               maxGroupMZ = 0.0015
)
saveRDS(object = fGroups, file = paste0("output\\data\\", Sys.Date(), "_fGroups_hilic_original.RDS"))


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Internal standard ----
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# screen for internal standard chlormequat-d9 (C5[2]H9H4ClN) or sotalol-d7 (7.49 min, C12H13H7[2]N2O3S m/z 280)
internal.standard <- data.frame(name = "chlormequat-d9", mz = 131.12941, rt = 527.64)
fGroups <- screenSuspects(
  fGroups = fGroups,
  suspects = internal.standard,
  rtWindow = 4,
  mzWindow = 0.0015,
  adduct = "[M+H]+"
)

# generate plot to visualize internal standard intensities and peak areas in all analyses
fGroups.is <- as.data.table(fGroups, areas = TRUE)
fGroups.is <- fGroups.is[which(!is.na(fGroups.is$susp_name)), ]

# prepare data frame for plotting
is.areas <- reshape2::melt(fGroups.is[, -c("group", "ret", "mz")])
colnames(is.areas) <- c("internal standard", "analysis", "area")

# generate histogram
plot.is <- ggplot(data = is.areas, aes(x = area)) +
  geom_histogram(binwidth = 500000, fill = "lightblue", alpha = 0.8) +
  theme_classic() +
  scale_y_continuous(expand = c(0, 0)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), legend.position = "none", legend.title = element_blank(), axis.title.x = element_blank())
plot.is
ggsave(plot = plot.is, file = paste0("output\\plots\\", Sys.Date(), "_internalstandards_hilic.png"), dpi = 330, background = "white")

# visually inspect peak shapes of internal standard
fGroups.is.only <- fGroups[, "M131_R527_13590"] #insert feature group of internal standard here
checkFeatures(fGroups.is.only)

# clean environment
rm(fGroups.is, is.areas, plot.is, fGroups.is.only)
