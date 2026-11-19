#------------------------------------------------------------------------------#
# NTS to risk assessment - feature finding using patRoon and prioritization    #
# using MS2Tox, MS2bioactivity, MS2Quant, structural alerts and trends         #
#                                                                              #
# HILIC only                                                                   #
#                                                                              #
# Developed by Nienke Meekel, KWR Water Research Institute                     #
#                                                                              #
# R version 4.5.2                                                              #
# patRoon version 2.3.4                                                        #
#                                                                              #
# Current version: 25 March 2026                                               #
# 404300/028 (From NTS to risk assessment)                                     #
# ---------------------------------------------------------------------------- #

# allocate additional memory
options(java.parameters = "-Xmx8000m")

# load required packages
library(patRoon) # for NTS feature finding and processing
library(dplyr) # for data handling
library(ggplot2) # for visualization
library(egg) # for extra ggplot themes

# check dependencies
options(patRoon.path.obabel = "C:/Program Files/OpenBabel-3.1.1") # specify location of openbabel
options(patRoon.path.SIRIUS = "C:/Program Files/sirius") # path to SIRIUS 5.8.5
verifyDependencies() # check if all dependencies are found

# set seed
seed <- 234
set.seed(seed)

# load structural alert functions
source("scripts\\structuralalertFunctions.R")


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Initialization ----
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

## Obtain analysis info
# read file with analysis info
anaInfo <- readxl::read_excel(path = "data\\analysisinfo.xlsx", sheet = "hilic_anonimized", col_names = TRUE)
anaInfo <- anaInfo[-c(1, 58), ] # skip the first and last row as these blanks are for stabilizing

## Load detected fGroups
fGroups <- readRDS(file = "output\\data\\fGroups_hilic_original.RDS")


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Features - filtering ----
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# Remove standards and other analyses that do not belong to samples - result NTS features only
fGroups <- patRoon::filter(fGroups, rGroups = c("S0", "S1", "S2", "S3", "S4", "S5", "blank_is", "blank_dw", "dw_addition", "blank_sw", "sw_addition"), negate = TRUE)

# apply general filters (remove features outside retention range, with lower intensity than 10,000, and occurring in the blank)
fGroups <- patRoon::filter(fGroups,
  preAbsMinIntensity = 100,
  retentionRange = c(144, 1200), # RT between 2.4 and 20 minutes
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

# inspect results and select feature groups that have an original peak area > 50k in any of the samples
fGroups.dt <- as.data.table(fGroups, areas = TRUE) # convert feature Groups object to data table
samples <- colnames(fGroups.dt)[-c(which(colnames(fGroups.dt) == "group" | colnames(fGroups.dt) == "mz" | colnames(fGroups.dt) == "ret"))] # select sample columns to filter on (i.e. do not filter on m/z, rt and group name)
fGroups.50k <- fGroups.dt %>% filter_at(samples, any_vars(. > 50000)) # apply 50k peak area filter
fGroups <- fGroups[, fGroups.50k$group] # select the corresponding feature groups

# remove variables we don't need anymore
rm(fGroups.dt, fGroups.50k, samples)

# save fGroups for non-target screening
saveRDS(object = fGroups, file = paste0("output\\data\\", Sys.Date(), "_fGroups_hilic_nts.RDS"))

# save workspace image
save.image(file = paste0("workspace images\\", Sys.Date(), "_fGroups_NTS.RData"))


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Componentization ----
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# cluster feature groups with similar chromatographic elution profiles and annotate by known chemical rules
components <- generateComponents(
  fGroups = fGroups,
  algorithm = "camera",
  ionization = "positive",
  onlyIsotopes = FALSE,
  minSize = 1,
  relMinReplicates = 0.3
)

# select feature groups with adducts [M+H]+, and assign [M+H]+ to those where the adduct is unknown
fGroups <- selectIons(fGroups, components, prefAdduct = "[M+H]+")

# find unique adducts and which feature groups belong to it
unique(adducts(fGroups))
which(adducts(fGroups) == "[M+Na+K]2+")
which(adducts(fGroups) == "[2M+Na+2K]3+")
which(adducts(fGroups) == "[M+H-HCOOH]+")

# remove feature groups with deviating adducts
fGroups <- fGroups[, -which(colnames(fGroups@groups) == "M447_R552_89981")] # [M+Na+K]2+
fGroups <- fGroups[, -which(colnames(fGroups@groups) == "M463_R553_92728")] # [2M+Na+2K]3+
fGroups <- fGroups[, -which(colnames(fGroups@groups) == "M579_R558_112330")] # [M+H-HCOOH]+


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# MS1 & MS2 data ----
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# Retrieve MS peak lists
set.seed(seed)
avgMSListParams <- getDefAvgPListParams(clusterMzWindow = 0.001)
mslists <- generateMSPeakLists(
  fGroups = fGroups,
  algorithm = "mzr",
  maxMSRtWindow = 5,
  precursorMzWindow = 0.8, # +/- 0.8 Da = 1.6 Da isolation width
  avgFeatParams = avgMSListParams,
  avgFGroupParams = avgMSListParams,
  topMost = NULL
)
saveRDS(object = mslists, file = paste0("output\\data\\", Sys.Date(), "_mslists_hilic.RDS")) # save mslists

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
# Suspect screening - MS2 database KWR ----
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# Load in-house library from .csv
suspects <- read.csv(file = "data\\libraries\\MS2 database referentiestoffen.csv", row.names = 1)

# screen fGroups against in-house library
fGroups <- screenSuspects(
  fGroups = fGroups,
  suspects = suspects,
  rtWindow = 12,
  mzWindow = 0.0015,
  adduct = "[M+H]+"
)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Annotation - formulas ----
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# Calculate candidate formulas using both MS1 and MS2 data
formulas.genform <- generateFormulas(
  fGroups = fGroups,
  MSPeakLists = mslists,
  algorithm = "genform",
  relMzDev = 5, # in ppm
  adduct = "[M+H]+",
  elements = "CHNOPSClBrF",
  hetero = TRUE,
  oc = TRUE,
  calculateFeatures = FALSE,
  featThresholdAnn = 1,
  absAlignMzDev = 0.0015, # in Da
  MSMode = "both", # include MS1 and MS2 data
  topMost = 5
)
saveRDS(object = formulas.genform, file = paste0("output\\data\\", Sys.Date(), "_formulas_hilic_genform.RDS"))

# Calculate candidate formulas using both MS1 and MS2 data
formulas.sirius <- generateFormulas( # use formulas for fingerprint predictions
  fGroups = fGroups,
  MSPeakLists = mslists,
  algorithm = "sirius",
  relMzDev = 5, # in ppm
  adduct = "[M+H]+",
  elements = "CHNOPSClBrF",
  profile = "orbitrap",
  database = "pubchem", # use PubChem as database for retrieval of formula candidates
  cores = 7,
  getFingerprints = TRUE,
  topMost = 5,
  calculateFeatures = FALSE,
  absAlignMzDev = 0.0015, # in Da
  extraOptsFormula = "--ppm-max-ms2=50",
  projectPath = "C:\\Users\\meekeni\\Documents\\sirius_hilic_formulas",
  verbose = TRUE,
  login = "interactive" # if no active login is present, username and password will be asked
)
saveRDS(object = formulas.sirius, file = paste0("output\\data\\", Sys.Date(), "_formulas_hilic_sirius.RDS"))


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Annotation - compounds ----
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# Calculate candidate structures and fingerprints using SIRIUS (version 5.8.5)
compounds.sirius <- generateCompounds(
  algorithm = "sirius",
  fGroups = fGroups,
  MSPeakLists = mslists,
  relMzDev = 5, # in ppm
  adduct = "[M+H]+",
  formulaDatabase = "pubchem",
  fingerIDDatabase = "pubchem",
  topMost = 5,
  topMostFormulas = 5,
  profile = "orbitrap",
  cores = 7,
  elements = "CHNOPSClBrF",
  extraOptsFormula = "--ppm-max-ms2=50", 
  projectPath = "C:\\Users\\meekeni\\Documents\\sirius_hilic_compounds",
  verbose = TRUE # ,
  # login = "interactive" # if no active login is present, username and password will be asked
)
saveRDS(object = compounds.sirius, file = paste0("output\\data\\", Sys.Date(), "_compounds_hilic_sirius.RDS"))

# Find candidate structures using MetFrag (version 2.5.0)
compounds.metfrag <-
  generateCompounds(
    fGroups = fGroups,
    MSPeakLists = mslistsF,
    algorithm = "metfrag",
    method = "CL",
    topMost = 5,
    dbRelMzDev = 5, # in ppm, for database search
    fragAbsMzDev = 0.05, # in Da, for fragment matching
    adduct = "[M+H]+",
    database = "pubchemlite",
    maxCandidatesToStop = 2500
  )
saveRDS(object = compounds.metfrag, file = paste0("output\\data\\", Sys.Date(), "_compounds_hilic_metfrag.RDS"))

# Load MassBank Europe and filter LC-HRMS data
mb <- loadMSLibrary("data\\libraries\\MassBank_NIST_2024_06.msp", "msp") # load MassBank.eu MSP library (downloaded from https://github.com/MassBank/MassBank-data/releases, release version 2024.06 (june 5), downloaded on 24th of july 2024)
unique(records(mb)[["Instrument_type"]]) # get all unique instrument types
mb <- patRoon::filter(mb, properties = list(Instrument_type = c("LC-ESI-ITFT", "ESI-ITFT", "ESI-QTOF", "LC-ESI-QTOF"))) # only select Orbitrap and QTOF data

# Find candidate structures using MassBank Europe (only look into [M+H]+ adducts)
compounds.mb <- generateCompounds(fGroups,
  MSPeakLists = mslistsF,
  algorithm = "library",
  MSLibrary = mb,
  minSim = 0.5,
  absMzDev = 0.05, # in Da
  spectrumType = "MS2",
  adduct = "[M+H]+",
  checkIons = "adduct"
)
saveRDS(object = compounds.mb, file = paste0("output\\data\\", Sys.Date(), "_compounds_hilic_massbank.RDS"))

# Find candidate structures using MassBank of North America (MoNA)
mona <- loadMSLibraryMSP( # load MoNA library (downloaded from https://mona.fiehnlab.ucdavis.edu/downloads, release version unknown, downloaded on 29th of October 2024)
  file = "data\\libraries\\MoNA-export-All_LC-MS-MS_Orbitrap.msp", parseComments = TRUE, prefCalcChemProps = TRUE,
  neutralChemProps = TRUE, potAdducts = TRUE, potAdductsLib = TRUE, absMzDev = 0.002, calcSPLASH = FALSE
)

# Find candidate structures using MassBank of North America (only look into [M+H]+ adducts)
compounds.mona <- generateCompounds(fGroups,
  MSPeakLists = mslistsF,
  algorithm = "library",
  MSLibrary = mona,
  minSim = 0.5,
  absMzDev = 0.05, # in Da
  spectrumType = "MS2",
  adduct = "[M+H]+",
  checkIons = "adduct"
) # Find candidates
saveRDS(object = compounds.mona, file = paste0("output\\data\\", Sys.Date(), "_compounds_hilic_mona.RDS"))

# save workspace image
save.image(file = paste0("workspace images\\", Sys.Date(), "_fGroups_NTS_formulas_compounds.RData"))


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Combine results ----
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# Incorporate formula information
compounds.mb <- addFormulaScoring(compounds = compounds.mb, formulas = formulas.sirius, updateScore = TRUE)
compounds.mona <- addFormulaScoring(compounds = compounds.mona, formulas = formulas.sirius, updateScore = TRUE)
compounds.metfrag <- addFormulaScoring(compounds = compounds.metfrag, formulas = formulas.sirius, updateScore = TRUE)

# select top-ranked candidates
compounds.mb <- patRoon::filter(compounds.mb, minExplainedPeaks = 2, topMost = 1)
compounds.mona <- patRoon::filter(compounds.mona, minExplainedPeaks = 2, topMost = 1)
compounds.sirius <- patRoon::filter(compounds.sirius, minExplainedPeaks = 2, topMost = 1)
compounds.metfrag <- patRoon::filter(compounds.metfrag, minExplainedPeaks = 2, topMost = 1)

# re-structure candidates obtained with SIRIUS
results.sirius <- patRoon::as.data.table(compounds.sirius, fGroups = fGroups) %>% # generate data table
  select(group, compoundName, score, SMILES, identifier, neutral_formula) %>% # select relevant variables
  dplyr::rename(
    compound.name.sirius = compoundName, # rename the variables for the final result table
    score.sirius = score,
    SMILES.sirius = SMILES,
    cid.sirius = identifier,
    formula.sirius = neutral_formula
  ) %>%
  mutate(cid.sirius = stringr::str_split_i(cid.sirius, ";", 1)) # split code so that CID identifier is extracted

# re-structure candidates obtained with MassBank
results.mb <- patRoon::as.data.table(compounds.mb, fGroups = fGroups) %>% # generate data table
  select(group, compoundName, score, libMatch, SMILES, identifier, neutral_formula) %>% # select relevant variables
  dplyr::rename(
    compound.name.mb = compoundName, # rename the variables for the final result table
    score.mb = score,
    libmatch.mb = libMatch,
    SMILES.mb = SMILES,
    identifier.mb = identifier,
    formula.mb = neutral_formula
  )

# re-structure candidates obtained with MoNA
results.mona <- patRoon::as.data.table(compounds.mona, fGroups = fGroups) %>% # generate data table
  select(group, compoundName, score, libMatch, SMILES, identifier, neutral_formula) %>% # select relevant variables
  dplyr::rename(
    compound.name.mona = compoundName, # rename the variables for the final result table
    score.mona = score,
    libmatch.mona = libMatch,
    SMILES.mona = SMILES,
    identifier.mona = identifier,
    formula.mona = neutral_formula
  )

# re-structure candidates obtained with MetFrag
results.metfrag <- patRoon::as.data.table(compounds.metfrag, fGroups = fGroups) %>% # generate data table
  select(group, compoundName, score, explainedPeaks, SMILES, identifier, neutral_formula, individualMoNAScore) %>% # select relevant variables
  dplyr::rename(
    compound.name.metfrag = compoundName, # rename the variables for the final result table
    score.metfrag = score,
    explainedpeaks.metfrag = explainedPeaks,
    SMILES.metfrag = SMILES,
    identifier.metfrag = identifier,
    formula.metfrag = neutral_formula,
    score.metfrag.mona = individualMoNAScore
  )

# re-structure formula candidates obtained with GenForm
results.genform <- patRoon::as.data.table(formulas.genform, fGroups = fGroups) %>% # generate data table
  select(group, neutral_formula, isoScore, explainedPeaks) %>% # elect relevant variables
  dplyr::rename(
    formula.genform = neutral_formula,
    isoScore.genform = isoScore,
    explainedPeaks.genform = explainedPeaks
  )

# Clean environment
rm(mb, mona)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Toxicity prediction - MS2Tox ----
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# predict toxicities
compounds.sirius <- predictTox(
  obj = compounds.sirius,
  LC50Mode = "static", # or "flow"
  type = "FP", # use SIRIUS+CSI:FingerID fingerprints for the prediction
  concUnit = "ugL"
) # determine LC50 values
View(as.data.table(compounds.sirius))

# assign toxicities to feature groups
fGroups.tox <- calculateTox(fGroups, compounds.sirius)
View(as.data.table(fGroups.tox))

# inspect predicted values (NOTE: check whether multiple predictions were made
# for 1 feature, the rest of the script assumes this is not the case as we selected the top 1 SIRIUS candidate)
toxicities(fGroups.tox)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Quantification with uncertainty - MS2Quant ----
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# load file with calibrants, concentrations and peak areas
calibrants <- readxl::read_excel(path = "data\\MS2Quant\\20241010_calibrants_xcalibur_hilic.xlsx", sheet = "MS2Quant_input", col_names = TRUE)

# remove entries where the concentration is 0, since MS2Quant cannot handle this
calibrants <- calibrants[-which(calibrants$conc == 0), ]
row.names(calibrants) <- NULL # reset row numbers

# remove all rows where the area is 0
calibrants <- calibrants[-which(calibrants$area == 0), ]
row.names(calibrants) <- NULL # reset row numbers

# remove all calibrants where we don't have enough data points, or calibrants that aren't linear
calibrants <- calibrants[-which(calibrants$name %in% c("tyrosine", "urotropine", "metamidophos", "diatrizoic acid", "tryptophane", "phenylalanine", "melamine")), ]
row.names(calibrants) <- NULL

# specify gradient
eluent <- data.frame(
  time = c(0, 4, 15, 26, 29.5), # check time stamp of last row, when does the gradient end?
  B = c(100, 90, 20, 20, 100)
)
eluent$time <- 60 * eluent$time # convert minutes to seconds

set.seed(seed)

# predict response factors
compounds.sirius <- predictRespFactors(
  compounds.sirius,
  fGroups = fGroups.tox,
  calibrants = calibrants,
  eluent = eluent,
  organicModifier = "MeCN",
  pHAq = 3,
  concUnit = "ugL",
  calibConcUnit = "ugL",
  type = "FP"
)

# inspect model specs
compounds.sirius@MS2QuantMeta

# save model plots
png(filename = paste0("output\\plots\\ms2quant\\", Sys.Date(), "_calibPlot.png"), res = 330, width = 15, height = 15, units = "cm")
compounds.sirius@MS2QuantMeta$calibPlot
dev.off()

png(filename = paste0("output\\plots\\ms2quant\\", Sys.Date(), "_individualPlots.png"), res = 330, width = 40, height = 20, units = "cm")
compounds.sirius@MS2QuantMeta$individualPlots
dev.off()

# predict concentrations
fGroups.haz <- calculateConcs(
  fGroups = fGroups.tox,
  featureAnn = compounds.sirius,
  areas = TRUE
)

# inspect predicted values (NOTE: check whether multiple predictions were made
# for 1 feature, the rest of the script assumes this is not the case as we selected the top 1 SIRIUS candidate)
View(concentrations(fGroups.haz))


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Structural alerts ----
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# obtain raw MS data from mslists object
ms.data <- unlist(mslists@peakLists, recursive = FALSE)

# retrieve fragments and neutral losses
ms.data.frag <- pbapply::pbsapply(ms.data, FUN = getFragments, simplify = FALSE, USE.NAMES = TRUE) # get fragments
ms.data.frag <- ms.data.frag[lapply(ms.data.frag, length) > 0] # remove feature groups without fragments
ms.data.nl <- pbapply::pbsapply(ms.data, FUN = getNeutralLosses, simplify = FALSE, USE.NAMES = TRUE) # get neutral losses
ms.data.nl <- ms.data.nl[lapply(ms.data.nl, length) > 0] # remove feature groups without neutral losses

# get predictions for the organophosphorus structural alerts
set.seed(seed)
result.organophosphorus <- getPredictionAlert(frags = ms.data.frag, nls = ms.data.nl, alert = "organophosphorus")
result.organophosphorus$name <- row.names(result.organophosphorus)
result.organophosphorus$group <- stringr::str_split_fixed(string = result.organophosphorus$name, pattern = "\\.", n = 2)[, 2]

# get predictions for the aromatic amine structural alerts
set.seed(seed)
result.aromaticamine <- getPredictionAlert(frags = ms.data.frag, nls = ms.data.nl, alert = "aromaticamine")
result.aromaticamine$name <- row.names(result.aromaticamine)
result.aromaticamine$group <- stringr::str_split_fixed(string = result.aromaticamine$name, pattern = "\\.", n = 2)[, 2]

# clean up data frames
result.aromaticamine$name <- NULL
fGroups.aromaticamine <- aggregate(alert ~ group, result.aromaticamine, median)
colnames(fGroups.aromaticamine) <- c("group", "aromatic.amine.alert")
result.organophosphorus$name <- NULL
fGroups.organophosphorus <- aggregate(alert ~ group, result.organophosphorus, median)
colnames(fGroups.organophosphorus) <- c("group", "organophosphorus.alert")

# clean environment
rm(ms.data, ms.data.frag, ms.data.nl, model.aromaticamine, model.organophosphorus, result.organophosphorus, result.aromaticamine)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Prioritization ----
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# Create data.table of fGroups and their peak areas, including hazard information
fGroups.haz <- patRoon::as.data.table(fGroups.haz,
  area = TRUE,
  toxAggrParams = getDefPredAggrParams(),
  concAggrParams = getDefPredAggrParams()
)

# combine results from feature groups, MassBank, SIRIUS, predicted toxicity and structural alerts
result <- fGroups.haz %>%
  rowwise() %>%
  mutate(maxConc = max(across(contains("_conc")))) %>%
  ungroup() %>%
  mutate(priorityScore = maxConc / (LC50)) %>%
  left_join(results.sirius, by = "group") %>%
  left_join(results.mb, by = "group") %>%
  left_join(results.mona, by = "group") %>%
  left_join(results.metfrag, by = "group") %>%
  left_join(fGroups.organophosphorus, by = "group") %>%
  left_join(fGroups.aromaticamine, by = "group") %>%
  mutate(ret = round(ret / 60, digits = 2)) # convert retention time from seconds to minutes

# add column on whether MS2 was available or not
fGroups.ms2 <- patRoon::filter(mslists, withMSMS = TRUE) # filter for feature groups with MS2 data
fGroups.ms2 <- names(fGroups.ms2@averagedPeakLists) # get names of feature groups with MS2 data
result$MS2 <- "no"
result$MS2[which(result$group %in% fGroups.ms2)] <- "yes"

# report average concentration per sample
result <- result %>%
  rowwise() %>%
  mutate(
    A_treated_ugL = mean(c_across(contains("pos-18_conc") | contains("pos-19_conc") | contains("pos-20_conc"))),
    A_untreated_ugL = mean(c_across(contains("pos-21_conc") | contains("pos-22_conc") | contains("pos-23_conc"))),
    A_P1_ugL = mean(c_across(contains("pos-25_conc") | contains("pos-26_conc") | contains("pos-27_conc"))),
    A_P2_ugL = mean(c_across(contains("pos-28_conc") | contains("pos-29_conc") | contains("pos-30_conc"))),
    A_P3_ugL = mean(c_across(contains("pos-31_conc") | contains("pos-32_conc") | contains("pos-33_conc"))),
    B_treated_ugL = mean(c_across(contains("pos-35_conc") | contains("pos-36_conc") | contains("pos-37_conc"))),
    B_untreated_ugL = mean(c_across(contains("pos-38_conc") | contains("pos-39_conc") | contains("pos-40_conc"))),
    B_P1_ugL = mean(c_across(contains("pos-41_conc") | contains("pos-42_conc") | contains("pos-43_conc"))),
    B_P2_ugL = mean(c_across(contains("pos-45_conc") | contains("pos-46_conc") | contains("pos-47_conc"))),
    B_P3_ugL = mean(c_across(contains("pos-48_conc") | contains("pos-49_conc") | contains("pos-50_conc")))
  )

# remove columns with _conc in colname, as we made new columns with average predicted concentration per sample
result <- result %>%
  select(-contains("_conc"))

# report average peak area per sample
result <- result %>%
  rowwise() %>%
  mutate(
    A_treated_area = mean(c_across(contains("pos-18") | contains("pos-19") | contains("pos-20"))),
    A_untreated_area = mean(c_across(contains("pos-21") | contains("pos-22") | contains("pos-23"))),
    A_P1_area = mean(c_across(contains("pos-25") | contains("pos-26") | contains("pos-27"))),
    A_P2_area = mean(c_across(contains("pos-28") | contains("pos-29") | contains("pos-30"))),
    A_P3_area = mean(c_across(contains("pos-31") | contains("pos-32") | contains("pos-33"))),
    B_treated_area = mean(c_across(contains("pos-35") | contains("pos-36") | contains("pos-37"))),
    B_untreated_area = mean(c_across(contains("pos-38") | contains("pos-39") | contains("pos-40"))),
    B_P1_area = mean(c_across(contains("pos-41") | contains("pos-42") | contains("pos-43"))),
    B_P2_area = mean(c_across(contains("pos-45") | contains("pos-46") | contains("pos-47"))),
    B_P3_area = mean(c_across(contains("pos-48") | contains("pos-49") | contains("pos-50")))
  )
result <- result[, -c(grep(pattern = "20190919HILIC_pos-", x = colnames(result)))] # remove measurement entries

# save workspace image
save.image(file = paste0("workspace images\\", Sys.Date(), "_fGroups_NTS_formulas_compounds_tox_quant_excl_mlinvitrotox_excl_ms2bioactivity.RData"))


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# MLinvitroTox ----
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# load MLinvitroTox results and create column with feature groups and formulas
mlinvitrotox <- read.csv(file = "data\\MLinvitroTox\\mlinvitrotox_model_predictions.csv")
mlinvitrotox$group <- stringr::str_split_i(mlinvitrotox$chem_id, "_", 2)
mlinvitrotox$formula <- stringr::str_split_i(mlinvitrotox$chem_id, "_", 3)
mlinvitrotox$group <- gsub(pattern = "-unknownCompound", replacement = "", x = mlinvitrotox$group)
mlinvitrotox$group <- gsub(pattern = "-", replacement = "_", x = mlinvitrotox$group)
mlinvitrotox$token <- paste0(mlinvitrotox$group, "|", mlinvitrotox$formula) # initialize token for matching

# initialize same token in result table
result$token.temp <- paste0(result$group, "|", result$formula.sirius)
mlinvitrotox.clean <- mlinvitrotox %>% filter(token %in% result$token.temp)

# combine MLinvitroTox results with result file
mlinvitrotox.active <- filter(
  mlinvitrotox.clean,
  similarity > 0.7 & probability > 0.7 # & endpoint_score > 0.5
)

# obtain statistics for the manuscript
unique(mlinvitrotox.active$MechanisticTarget) # get unique mechanistic target (note the 'not annotated' class, we do not include this in the total count)
length(unique(mlinvitrotox.active$aeid)) # get number of unique asssay id's

# add mlinvitrotox predictions to result file
result$mlinvitrotox.active <- NA # initialize column for mlinvitrotox results
result$mlinvitrotox.active[which(result$token.temp %in% mlinvitrotox.clean$token)] <- 0 # when a prediction is made by mlinvitrotox, turn NA into 0
result$mlinvitrotox.active[which(result$token.temp %in% mlinvitrotox.active$token)] <- 1 # when this prediction is active, turn 0 into 1

# remove token from result data frame
result$token.temp <- NULL


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# MS2Bioactivity ----
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# load MS2Bioactivity results
ms2bioactivity <- read.table(file = "data\\MS2Bioactivity\\hilic_compounds_predictions.tsv", header = TRUE)
colnames(ms2bioactivity) <- c("fgroup.ms2bioactivity", "formula.ms2bioactivity", "nr.ahr")

# initialize columns in results data frame for MS2bioactivity output
result$nr.ahr <- NA
result$formula.ms2bioactivity <- NA

# fill results data frame with MS2bioactivity results
for (i in 1:nrow(ms2bioactivity)) {
  fgroup <- ms2bioactivity$fgroup.ms2bioactivity[i]
  formula <- ms2bioactivity$formula.ms2bioactivity[i]
  x <- grep(fgroup, x = result$group)
  y <- grep(formula, x = result$formula.sirius)

  if (length(x) == 1) {
    result$nr.ahr[x] <- ms2bioactivity$nr.ahr[i]
    result$formula.ms2bioactivity[x] <- ms2bioactivity$formula.ms2bioactivity[i]
  } else { # in case an fGroup code from MS2bioactivity appears multiple times in the results dataframe, look into the formula predicted by SIRIUS
    z <- x[which(x %in% y)]
    result$nr.ahr[z] <- ms2bioactivity$nr.ahr[i]
    result$formula.ms2bioactivity[z] <- ms2bioactivity$formula.ms2bioactivity[i]
    rm(z)
  }
  rm(fgroup, x, y, formula)
}

# clean environment
rm(i)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Export ----
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# reorder columns for readability
column.order <- c(
  "group", "ret", "mz", "adduct", "neutralMass", "MS2", "LC50", "maxConc", "priorityScore", "mlinvitrotox.active",
  "nr.ahr", "aromatic.amine.alert", "organophosphorus.alert", "susp_name", "compound.name.sirius", "score.sirius",
  "SMILES.sirius", "cid.sirius", "formula.sirius", "compound.name.mb", "score.mb", "libmatch.mb", "SMILES.mb",
  "identifier.mb", "formula.mb", "compound.name.mona", "score.mona", "libmatch.mona", "SMILES.mona",
  "identifier.mona", "formula.mona", "compound.name.metfrag", "score.metfrag", "explainedpeaks.metfrag", "SMILES.metfrag",
  "identifier.metfrag", "formula.metfrag", "score.metfrag.mona", "A_treated_area", "A_untreated_area", "A_P1_area", "A_P2_area",
  "A_P3_area", "B_treated_area", "B_untreated_area", "B_P1_area", "B_P2_area", "B_P3_area", "A_treated_ugL", "A_untreated_ugL",
  "A_P1_ugL", "A_P2_ugL", "A_P3_ugL", "B_treated_ugL", "B_untreated_ugL", "B_P1_ugL", "B_P2_ugL", "B_P3_ugL", "formula.ms2bioactivity",
  "conc_types", "LC50_types"
)

result <- result[, column.order]

# export results overview
xlsx::write.xlsx(
  x = as.data.frame(result), file = paste0("output\\data\\", Sys.Date(), "_results_hilic.xlsx"), sheetName = "general",
  col.names = TRUE, row.names = FALSE, append = FALSE
)
xlsx::write.xlsx( # add mlinvitrotox results (filtered!) on a separate sheet
  x = as.data.frame(mlinvitrotox.active), file = paste0("output\\data\\", Sys.Date(), "_results_hilic.xlsx"), sheetName = "MLinvitroTox",
  col.names = TRUE, row.names = FALSE, append = TRUE
)

# export mlinvitrotox results to separate xlsx file
xlsx::write.xlsx( # add mlinvitrotox results (filtered!) on a separate sheet
  x = as.data.frame(mlinvitrotox.active), file = paste0("output\\data\\", Sys.Date(), "_results_hilic_mlinvitrotox.xlsx"), sheetName = "MLinvitroTox",
  col.names = TRUE, row.names = FALSE, append = FALSE
)

# save workspace image
save.image(file = paste0("workspace images\\", Sys.Date(), "_fGroups_NTS_formulas_compounds_tox_quant.RData"))


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Plot chromatograms ----
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# peakshapes of drinking water relevant features
peakshapes <- c(
  "M100_R542_2335", "M128_R495_12706", "M136_R227_14721", "M136_R216_14726",
  "M143_R543_16451", "M148_R424_17485", "M150_R196_17930", "M150_R227_18093",
  "M150_R202_18099", "M154_R535_18934", "M158_R276_19736", "M158_R443_19837",
  "M163_R589_20771", "M164_R217_21090", "M171_R555_23219", "M172_R339_23313",
  "M178_R202_27213", "M200_R216_33167", "M222_R630_37881", "M222_R221_38013",
  "M227_R314_39231", "M242_R416_42594", "M244_R379_42993", "M255_R444_45821",
  "M258_R516_46354", "M277_R311_50670", "M286_R455_55506", "M291_R282_56677",
  "M299_R433_59220", "M307_R310_60761", "M313_R406_62690", "M325_R412_65009",
  "M326_R238_65091", "M339_R400_69489", "M341_R225_70372", "M341_R403_70612",
  "M351_R373_72809", "M357_R441_74517", "M358_R596_74792", "M369_R373_76984",
  "M375_R369_78254", "M385_R454_79585", "M479_R557_95787", "M600_R159_115388"
)

# plot chromatograms of the selected fGroups
for (i in peakshapes) {
  png(
    file = paste0("output\\plots\\chromatograms\\fGroup_", i, ".png"),
    width = 15, height = 15, units = "cm", res = 330
  )
  plotChroms(fGroups[, i],
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


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Plot m/z rt distribution with toxquant data ----
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

library(ggplot2) # for visualization
library(egg) # for plot themes
library(scales) # for axis labels

# create plot for the m/z and retention time distribution
plot.mz.rt <- ggplot(data = result, aes(x = ret, y = mz, color = MS2)) +
  geom_point(size = 2, alpha = 0.6) +
  scale_colour_manual(values = c("yes" = "#6c88c4", "no" = "#de425b")) +
  xlab("retention time (min)") +
  ylab("m/z") +
  # ylim(0, 1000) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1300)) +
  scale_x_continuous(expand = c(0, 0), limits = c(0, 15)) +
  theme_article() +
  theme(legend.position = "none")
plot.mz.rt
ggsave(plot = plot.mz.rt, filename = paste0("output\\plots\\", Sys.Date(), "_fGroups_mz_rt.png"), width = 15, height = 10, units = "cm")

# Focus on fGroups with MS2
fGroups.ms2 <- fGroups[, fGroups.ms2]
data.plot.ms2 <- result %>% filter(MS2 == "yes")

# create plot with m/z and retention time distribution, and priority score
plot.mz.rt.priority <- ggplot(data = data.plot.ms2, aes(x = ret, y = mz, color = priorityScore)) +
  geom_point(data = data.plot.ms2 %>% filter(is.na(priorityScore)), size = 2, alpha = 0.8) + # plot NA in the first layer so that they do not overlap other data points
  geom_point(data = data.plot.ms2 %>% filter(!is.na(priorityScore)), size = 2, alpha = 0.8) +
  scale_color_gradient(low = "#6c88c4", high = "#de425b", na.value = "lightgrey") +
  xlab("retention time (min)") +
  ylab("m/z") +
  labs(colour = "priority score") +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1000)) +
  scale_x_continuous(expand = c(0, 0), limits = c(0, 15)) +
  theme_article() +
  theme(legend.position = c(.95, .95), legend.justification = c("right", "top"), legend.margin = margin(6, 6, 6, 6))
plot.mz.rt.priority
ggsave(plot = plot.mz.rt.priority, filename = paste0("output\\plots\\", Sys.Date(), "_fGroups_mz_rt_priorityscore.png"), width = 15, height = 10, units = "cm")

# create plot with m/z and retention time distribution, and predicted LC50 value
plot.mz.rt.lc50 <- ggplot(data = data.plot.ms2, aes(x = ret, y = mz, color = LC50)) +
  geom_point(data = data.plot.ms2 %>% filter(is.na(LC50)), size = 2, alpha = 0.8) + # plot NA in the first layer so that they do not overlap other data points
  geom_point(data = data.plot.ms2 %>% filter(!is.na(LC50)), size = 2, alpha = 0.8) +
  scale_color_gradient(low = "#de425b", high = "#6c88c4", na.value = "lightgrey") +
  xlab("retention time (min)") +
  ylab("m/z") +
  labs(colour = "LC50 (\U03BCg/L)") +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1000)) +
  scale_x_continuous(expand = c(0, 0), limits = c(0, 15)) +
  theme_article() +
  theme(legend.position = c(.95, .95), legend.justification = c("right", "top"), legend.margin = margin(6, 6, 6, 6))
plot.mz.rt.lc50
ggsave(plot = plot.mz.rt.lc50, filename = paste0("output\\plots\\", Sys.Date(), "_fGroups_mz_rt_LC50.png"), width = 15, height = 10, units = "cm")

# create plot with m/z and retention time distribution, and predicted concentration
plot.mz.rt.conc <- ggplot(data = data.plot.ms2, aes(x = ret, y = mz, color = maxConc)) +
  geom_point(data = data.plot.ms2 %>% filter(is.na(maxConc)), size = 2, alpha = 0.8) + # plot NA in the first layer so that they do not overlap other data points
  geom_point(data = data.plot.ms2 %>% filter(!is.na(maxConc)), size = 2, alpha = 0.8) +
  scale_color_gradient(low = "#6c88c4", high = "#de425b", na.value = "lightgrey") +
  xlab("retention time (min)") +
  ylab("m/z") +
  labs(colour = "concentration (\U03BCg/L)") +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1000)) +
  scale_x_continuous(expand = c(0, 0), limits = c(0, 15)) +
  theme_article() +
  theme(legend.position = c(.95, .95), legend.justification = c("right", "top"), legend.margin = margin(6, 6, 6, 6))
plot.mz.rt.conc
ggsave(plot = plot.mz.rt.conc, filename = paste0("output\\plots\\", Sys.Date(), "_fGroups_mz_rt_conc.png"), width = 15, height = 10, units = "cm")

# create plot with predicted concentration and predicted LC50
plot.conc.lc50 <- ggplot(data = data.plot.ms2, aes(x = maxConc, y = LC50, color = priorityScore)) +
  geom_point(data = data.plot.ms2 %>% filter(!is.na(maxConc)), size = 2, alpha = 0.8) +
  scale_color_gradient(low = "#6c88c4", high = "#de425b", na.value = "lightgrey") +
  xlab("concentration (\U03BCg/L)") +
  ylab("LC50 (\U03BCg/L)") +
  labs(colour = "priority score") +
  scale_y_continuous(transform = "log10", breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x))) +
  scale_x_continuous(transform = "log10", breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x))) +
  theme_article() +
  theme(legend.position = c(.95, .95), legend.justification = c("right", "top"), legend.margin = margin(6, 6, 6, 6))
plot.conc.lc50
ggsave(plot = plot.conc.lc50, filename = paste0("output\\plots\\", Sys.Date(), "_fGroups_conc_LC50.png"), width = 15, height = 10, units = "cm")

# create plot with m/z and retention time distribution, and MLinvitroTox activity
data.plot.ms2.mlinvitrotox <- data.plot.ms2 %>%
  mutate(
    mlinvitrotox.active = replace(mlinvitrotox.active, mlinvitrotox.active == 1, "active")
  ) %>%
  mutate(mlinvitrotox.active = replace(mlinvitrotox.active, mlinvitrotox.active == 0, "not active"))

plot.mz.rt.mlinvitrotox <- ggplot(data = data.plot.ms2.mlinvitrotox, aes(x = ret, y = mz, color = mlinvitrotox.active)) +
  geom_point(data = data.plot.ms2.mlinvitrotox %>% filter(is.na(mlinvitrotox.active)), size = 2, alpha = 0.8) + # plot NA in the first layer so that they do not overlap other data points
  geom_point(data = data.plot.ms2.mlinvitrotox %>% filter(!is.na(mlinvitrotox.active)), size = 2, alpha = 0.8) +
  scale_colour_manual(values = c("not active" = "#6c88c4", "active" = "#de425b"), na.value = "lightgrey") +
  xlab("retention time (min)") +
  ylab("m/z") +
  labs(colour = "MLinvitroTox") +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1000)) +
  scale_x_continuous(expand = c(0, 0), limits = c(0, 15)) +
  theme_article() +
  theme(legend.position = c(.95, .95), legend.justification = c("right", "top"), legend.margin = margin(6, 6, 6, 6))
plot.mz.rt.mlinvitrotox
ggsave(plot = plot.mz.rt.mlinvitrotox, filename = paste0("output\\plots\\", Sys.Date(), "_fGroups_mz_rt_mlinvitrotox.png"), width = 15, height = 10, units = "cm")

# create plot with m/z and retention time distribution, and mlinvitrotox activity combined with the priority score
data.plot.ms2$ps.ml <- data.plot.ms2$priorityScore * data.plot.ms2$mlinvitrotox.active
plot.mz.rt.mlinvitrotox.priority <- ggplot(data = data.plot.ms2, aes(x = ret, y = mz, color = ps.ml)) +
  geom_point(data = data.plot.ms2 %>% filter(is.na(ps.ml)), size = 2, alpha = 0.8) + # plot NA in the first layer so that they do not overlap other data points
  geom_point(data = data.plot.ms2 %>% filter(!is.na(ps.ml)), size = 2, alpha = 0.8) +
  scale_color_gradient(low = "#6c88c4", high = "#de425b", na.value = "lightgrey") +
  xlab("retention time (min)") +
  ylab("m/z") +
  labs(colour = "priority score") +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1000)) +
  scale_x_continuous(expand = c(0, 0), limits = c(0, 15)) +
  theme_article() +
  theme(legend.position = c(.95, .95), legend.justification = c("right", "top"), legend.margin = margin(6, 6, 6, 6))
plot.mz.rt.mlinvitrotox.priority
ggsave(plot = plot.mz.rt.mlinvitrotox.priority, filename = paste0("output\\plots\\", Sys.Date(), "_fGroups_mz_rt_mlinvitrotox_priority.png"), width = 15, height = 10, units = "cm")


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Plot drinking water relevant m/z rt distribution ----
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# plot priority scores for selected drinking water relevant mechanistic targets
data.plot.ms2$ps.ml.selected <- NA
mlinvitrotox.dw <- mlinvitrotox.clean %>% filter(MechanisticTarget %in% c(
  "Androgen Receptor Modulation", "Aryl Hydrocarbon Receptor Modulation", "Cell Growth", "Cell Morphology", "DNA Damage", "DNA Repair", "ER Alpha", "Estrogen Receptor Modulation",
  "Glucocorticoid Receptor Modulation", "Neurotransmission", "NR3C1", "Oxidative Stress", "p53 Modulation", "Peroxisome Proliferator Activated Receptor Modulation", "PPARG",
  "Pregnane X Receptor Modulation", "Progesterone Receptor Modulation", "Proliferation", "Thyroid Receptor Modulation", "TP53", "RAR-related Orphan Receptor Modulation",
  "Retinoic Acid Receptor Modulation", "Retinoid X Receptor Modulation", "Other Developmental Signaling Transcription Factors"
))

# select active features only, for the drinking water relevant mechanistic targets
mlinvitrotox.dw <- filter(
  mlinvitrotox.dw,
  similarity > 0.7 & probability > 0.7
)

# initialize same token as in result table
data.plot.ms2$token.temp <- paste0(data.plot.ms2$group, "|", data.plot.ms2$formula.sirius)

# add mlinvitrotox predictions to plot data frame
data.plot.ms2$ps.ml.selected <- 0 # initialize column for dw relevant mlinvitrotox results
data.plot.ms2$ps.ml.selected[which(data.plot.ms2$token.temp %in% mlinvitrotox.dw$token)] <- 1 # when this prediction is active, turn 0 into 1
data.plot.ms2$ps.ml.selected <- data.plot.ms2$ps.ml.selected * data.plot.ms2$ps.ml
data.plot.ms2$ps.ml.selected[which(data.plot.ms2$ps.ml.selected == 0)] <- NA

# include aromatic amine structural alert information in the plot
data.plot.ms2$alert.shape <- "no"
data.plot.ms2$alert.shape[which(data.plot.ms2$aromatic.amine.alert >= 0.5)] <- "yes"

plot.mz.rt.mlinvitrotox.priority.dw <- ggplot(data = data.plot.ms2, aes(x = ret, y = mz, color = ps.ml.selected)) +
  geom_point(data = data.plot.ms2 %>% filter(is.na(ps.ml.selected)), size = 2, alpha = 0.8, aes(shape = alert.shape)) + # plot NA in the first layer so that they do not overlap other data points
  geom_point(data = data.plot.ms2 %>% filter(!is.na(ps.ml.selected)), size = 2, alpha = 0.8, aes(shape = alert.shape)) +
  scale_color_gradient(low = "#6c88c4", high = "#de425b", na.value = "lightgrey") +
  xlab("retention time (min)") +
  ylab("m/z") +
  labs(colour = "priority score", shape = "aromatic amine") +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1000)) +
  scale_x_continuous(expand = c(0, 0), limits = c(0, 15)) +
  theme_article() +
  theme(legend.position = c(.95, .95), legend.justification = c("right", "top"), legend.margin = margin(6, 6, 6, 6))
plot.mz.rt.mlinvitrotox.priority.dw
ggsave(plot = plot.mz.rt.mlinvitrotox.priority.dw, filename = paste0("output\\plots\\", Sys.Date(), "_fGroups_mz_rt_mlinvitrotox_priority_drinking_water_relevant_alert.png"), width = 15, height = 10, units = "cm")

# export data frame for inspection
xlsx::write.xlsx(
  x = as.data.frame(data.plot.ms2), file = paste0("output\\data\\", Sys.Date(), "_results_hilic_prioritization_combined_scores.xlsx"), sheetName = "general",
  col.names = TRUE, row.names = FALSE, append = FALSE
)
