


# Script automatically generated on Fri Apr 21 09:45:34 2023
# and modified by GF


library(patRoon)
library(xcms)
library(tidyverse)

BiocParallel::register(BiocParallel::SerialParam(), default = TRUE)
memory.size(500000)


# -------------------------
# initialization
# -------------------------


# MetFrag CommandLine local database path
# MetFrag CommandLine downloaded from https://github.com/ipb-halle/MetFragRelaunched/releases/tag/v.2.5.0
options(patRoon.path.MetFragCL = "C:/MetFrag/MetFragCommandLine-2.5.0.jar")



# Load analysis table
anaInfo <- read.csv("analyses.csv")



# -------------------------
# features
# -------------------------

# Find all features
# NOTE: see the XCMS manual for many more options
fList <- findFeatures(anaInfo, "xcms3",
                      xcms::CentWaveParam(ppm = 7,
                                          peakwidth = c(4.92, 147.50),
                                          snthresh = 10,
                                          prefilter = c(3, 100),
                                          mzCenterFun = "wMean",
                                          integrate = 1,
                                          mzdiff = -0.0175,
                                          fitgauss = FALSE,
                                          noise = 0,
                                          firstBaselineCheck = TRUE),
                      verbose = FALSE)


# Group and align features between analyses
fGroups <- groupFeatures(fList, "xcms3", rtalign = TRUE,
                         groupParam = xcms::PeakDensityParam(sampleGroups = anaInfo$group,
                                                             bw = 24.76,
                                                             minFraction = 0,
                                                             minSamples = 1,
                                                             binSize = 0.045,
                                                             maxFeatures = 50),
                         
                         retAlignParam = xcms::ObiwarpParam(binSize = 0.045,
                                                            centerSample = integer(0), 
                                                            response = 1, 
                                                            distFun = "cor_opt",
                                                            gapInit = 1.16,
                                                            gapExtend = 3,
                                                            factorDiag = 2,
                                                            factorGap = 1,
                                                            localAlignment = FALSE,
                                                            initPenalty = 0,
                                                            subset = integer(0),
                                                            subsetAdjust = "average"))

# Basic rule based filtering
#fGroups <- filter(fGroups, preAbsMinIntensity = 100, absMinIntensity = 10000, relMinReplicateAbundance = 1,
#                  maxReplicateIntRSD = 0.75, blankThreshold = 5, removeBlanks = TRUE,
#                  retentionRange = NULL, mzRange = NULL)


# Try to export before the components
reportCSV(fGroups, path = "report_01justbeforecomponent", formulas = NULL, compounds = NULL, components = NULL)

# -------------------------
# componentization
# -------------------------

# Perform automatic generation of components
components <- generateComponents(fGroups, "camera", ionization = "positive")

# Try to export after the components
reportCSV(fGroups, path = "report_02justaftercomponent", formulas = NULL, compounds = NULL, components = components)



fGroups <- selectIons(fGroups, components, prefAdduct = "[M+H]+", onlyMonoIso = TRUE)

# Try to export after the selectIons
reportCSV(fGroups, path = "report_03justafterselectIons", formulas = NULL, compounds = NULL, components = components)

# -------------------------
# annotation
# -------------------------

# Retrieve MS peak lists
avgMSListParams <- getDefAvgPListParams(clusterMzWindow = 0.001, minIntensityPre = 1E3, minIntensityPost = 1E3)
mslists <- generateMSPeakLists(fGroups, "mzr", maxMSRtWindow = 20, precursorMzWindow = 1,
                               avgFeatParams = avgMSListParams,
                               avgFGroupParams = avgMSListParams)
# Rule based filtering of MS peak lists. You may want to tweak this. See the manual for more information.
#mslists <- filter(mslists, absMSIntThr = NULL, absMSMSIntThr = NULL, relMSIntThr = NULL, relMSMSIntThr = 0.05,
#                  topMSPeaks = NULL, topMSMSPeaks = 25)


## I skipped the formula!
# Calculate formula candidates
#formulas <- generateFormulas(fGroups, mslists, "genform", relMzDev = 5, adduct = "[M+H]+", elements = "CHNOClBrSP",
#                             calculateFeatures = FALSE, featThresholdAnn = 0.75)

# Calculate compound structure candidates with PubChemLite
# Pubchemlite local database path
## Pubchemlite database in csv format downloaded from https://zenodo.org/record/7684618
options(patRoon.path.MetFragPubChemLite = "C:/PubChemLite/PubChemLite_exposomics_20230224.csv")
myLocalDatabasePath <- "C:/PubChemLite/PubChemLite_exposomics_20230224.csv"

compounds_PCL <- generateCompounds(fGroups, mslists, "metfrag", method = "CL",
                                   dbRelMzDev = 5, fragRelMzDev = 5,
                                   fragAbsMzDev = 0.002, adduct = "[M+H]+",
                                   database = "csv", extraOpts = list(LocalDatabasePath = myLocalDatabasePath),
                                   scoreTypes = c("fragScore","score", "individualMoNAScore",
                                                  "AnnoTypeCount","Patent_Count", "PubMed_Count"),
                                   maxCandidatesToStop = 100)
compounds_PCL <- addFormulaScoring(compounds_PCL, formulas, updateScore = TRUE)

# -------------------------
# reporting - for PubChemLite
# -------------------------

reportCSV(fGroups, path = "report_PCL", formulas = NULL, compounds = compounds_PCL, components = components)

# Summary of MetFrag Results in a a Single Table 

MFsummary_PCL <- as.data.table(compounds_PCL)
write_tsv(MFsummary_PCL, "MFsummary_PCL.txt")


#########################
##
##
##

# Calculate compound structure candidates with WormJam
## The database in csv format was created flowing the RMarkdown prepared by Emma S.: https://gitlab.lcsb.uni.lu/eci/pubchem-docs/-/blob/main/taxonomy/Celegans/C_elegans_metabolites.Rmd

WormJamDatabasePath <- "C:/C_elegans_all_Struct_Info.csv"

compounds_WJ <- generateCompounds(fGroups, mslists, "metfrag", method = "CL",
                                  dbRelMzDev = 5, fragRelMzDev = 5,
                                  fragAbsMzDev = 0.002, adduct = "[M+H]+",
                                  database = "csv", extraOpts = list(LocalDatabasePath = WormJamDatabasePath),
                                  scoreTypes = c("LiteratureCount", "PatentCount"),
                                  maxCandidatesToStop = 100)
compounds_WJ <- addFormulaScoring(compounds_WJ, formulas, updateScore = TRUE)

# -------------------------
# reporting  - for wormjam
# -------------------------

reportCSV(fGroups, path = "report_WJ", formulas = NULL, compounds = compounds_WJ, components = components)

# Summary of MetFrag Results in a a Single Table 

MFsummary_WJ <- as.data.table(compounds_WJ)
write_tsv(MFsummary_WJ, "MFsummary_WJ.txt")






##
####
########
###################### .. and WITH the formula:

# Calculate formula candidates
formulas <- generateFormulas(fGroups, mslists, "genform", relMzDev = 5, adduct = "[M+H]+", elements = "CHNOClBrSP",
                             calculateFeatures = FALSE, featThresholdAnn = 0.75)

# Calculate compound structure candidates with PubChemLite
# Pubchemlite local database path
## Pubchemlite database in csv format downloaded from https://zenodo.org/record/7684618
options(patRoon.path.MetFragPubChemLite = "C:/PubChemLite/PubChemLite_exposomics_20230224.csv")
myLocalDatabasePath <- "C:/PubChemLite/PubChemLite_exposomics_20230224.csv"

compounds_PCL_with_formulas <- generateCompounds(fGroups, mslists, "metfrag", method = "CL",
                                   dbRelMzDev = 5, fragRelMzDev = 5,
                                   fragAbsMzDev = 0.002, adduct = "[M+H]+",
                                   database = "csv", extraOpts = list(LocalDatabasePath = myLocalDatabasePath),
                                   scoreTypes = c("fragScore","score", "individualMoNAScore",
                                                  "AnnoTypeCount","Patent_Count", "PubMed_Count"),
                                   maxCandidatesToStop = 100)
compounds_PCL_with_formulas  <- addFormulaScoring(compounds_PCL_with_formulas , formulas, updateScore = TRUE)

# -------------------------
# reporting - for PubChemLite
# -------------------------

reportCSV(fGroups, path = "report_PCL_with_formulas", formulas = formulas, compounds = compounds_PCL_with_formulas, components = components)

# Summary of MetFrag Results in a a Single Table 

MFsummary_PCL_with_formulas <- as.data.table(compounds_PCL_with_formulas)
write_tsv(MFsummary_PCL_with_formulas, "MFsummary_PCL_with_formulas.txt")


#########################
##
##
##

# Calculate compound structure candidates with WormJam
## The database in csv format was created flowing the RMarkdown prepared by Emma S.: https://gitlab.lcsb.uni.lu/eci/pubchem-docs/-/blob/main/taxonomy/Celegans/C_elegans_metabolites.Rmd

WormJamDatabasePath <- "C:/C_elegans_all_Struct_Info.csv"

compounds_WJ_with_formulas <- generateCompounds(fGroups, mslists, "metfrag", method = "CL",
                                  dbRelMzDev = 5, fragRelMzDev = 5,
                                  fragAbsMzDev = 0.002, adduct = "[M+H]+",
                                  database = "csv", extraOpts = list(LocalDatabasePath = WormJamDatabasePath),
                                  scoreTypes = c("LiteratureCount", "PatentCount"),
                                  maxCandidatesToStop = 100)
compounds_WJ_with_formulas <- addFormulaScoring(compounds_WJ_with_formulas, formulas, updateScore = TRUE)

# -------------------------
# reporting  - for wormjam
# -------------------------

reportCSV(fGroups, path = "report_WJ_with_formulas", formulas = formulas, compounds = compounds_WJ_with_formulas, components = components)

# Summary of MetFrag Results in a a Single Table 

MFsummary_WJ_with_formulas <- as.data.table(compounds_WJ_with_formulas)
write_tsv(MFsummary_WJ_with_formulas, "MFsummary_WJ_with_formulas.txt")






####### WJ bis!!!!  2023-08-10
##
##
##

# MetFrag CommandLine local database path
# MetFrag CommandLine downloaded from https://github.com/ipb-halle/MetFragRelaunched/releases/tag/v.2.5.0
options(patRoon.path.MetFragCL = "C:/MetFrag/MetFragCommandLine-2.5.0.jar")



# Calculate compound structure candidates with WorJam
## The database in csv format was created folowing the RMarkdown prepared by Emma S.: https://gitlab.lcsb.uni.lu/eci/pubchem-docs/-/blob/main/taxonomy/Celegans/C_elegans_metabolites.Rmd

WormJamDatabasePath <- "C:/C_elegans_all_Struct_Info.csv"

fGroups <- normInts(fGroups)


compounds_WJ_no_formula_bis <- generateCompounds(fGroups, mslists, "metfrag", method = "CL",
                                                 dbRelMzDev = 5, fragRelMzDev = 5,
                                                 fragAbsMzDev = 0.002, adduct = "[M+H]+",
                                                 database = "csv", extraOpts = list(LocalDatabasePath = WormJamDatabasePath),
                                                 scoreTypes = c("fragScore","score", "individualMoNAScore",
                                                                "LiteratureCount", "PatentCount"),
                                                 maxCandidatesToStop = 100)


# -------------------------
# reporting  - for wormjam
# -------------------------

reportCSV(fGroups, path = "report_WJ_no_formula_bis", formulas = NULL, compounds = compounds_WJ_no_formula_bis, components = components)

# Summary of MetFrag Results in a a Single Table 

MFsummary_WJ_no_formula_bis <- as.data.table(compounds_WJ_no_formula_bis)
write_tsv(MFsummary_WJ_no_formula_bis, "MFsummary_WJ_no_formula_bis.txt")




##### April 2024
### Suspect screening!

##fixing the fgroups errors:
fGroups@ISTDs <-  data.table::data.table()
fGroups@ISTDAssignments <-  data.table::data.table()
fGroups@concentrations <- fGroups@toxicities <- data.table::data.table()
fGroups@analysisInfo[["path"]][1:length(fGroups@analysisInfo[["path"]])] <- rep("C:/Users/frigerio.gianfranco/OneDrive - Ospedale San Raffaele/Documents/Lux/C Elegans/2022-09-20_elaborations/Patroon/HILIC_POS_Katie/PatRoon_only_new_2023_05_23/mzmLdata", length(fGroups@analysisInfo[["path"]]))
fGroups@features@analysisInfo[["path"]][1:length(fGroups@features@analysisInfo[["path"]])] <- rep("C:/Users/frigerio.gianfranco/OneDrive - Ospedale San Raffaele/Documents/Lux/C Elegans/2022-09-20_elaborations/Patroon/HILIC_POS_Katie/PatRoon_only_new_2023_05_23/mzmLdata", length(fGroups@analysisInfo[["path"]]))
fGroups@features@xdata@processingData@files <- gsub("C:\\Users\\gianfranco.frigerio\\OneDrive - University of Luxembourg\\Documents\\C Elegans\\2022-09-20_elaborations\\Patroon\\HILIC_POS_Katie\\PatRoon_only_new_2023_05_23\\mzmLdata\\",
                                                    "C:\\Users\\frigerio.gianfranco\\OneDrive - Ospedale San Raffaele\\Documents\\Lux\\C Elegans\\2022-09-20_elaborations\\Patroon\\HILIC_POS_Katie\\PatRoon_only_new_2023_05_23\\mzmLdata\\",
                                                    fGroups@features@xdata@processingData@files,
                                                    fixed=TRUE)




## loading the databases and preparing them for the suspect screening:

WJ_andPubChemRelMetabolites <- read_csv("C_elegans_all_Struct_Info.csv") %>%  ## database used during patRoon elaboration
  mutate(idname = paste0("WJplus_ID_", (1:length(Identifier))))

suspect_WJplus <- WJ_andPubChemRelMetabolites %>%
  select(name = idname, formula = MolecularFormula) %>%
  mutate(adduct = "[M+H]+")




WJ_zenodo <- read_tsv("WormJamMetabolites_20190909.txt") %>%
  mutate(idname = paste0("WJ_zenodo_ID_", (1:length(Name_neutral))))

suspect_WJ_zenodo <- WJ_zenodo %>%
  select(name = idname, formula = Formula_Neutral) %>%
  mutate(adduct = "[M+H]+")


WormJam_v0.1.0 <- read_tsv("WormJam_v0.1.0.txt")%>%
  mutate(idname = paste0("WJ_v0.1.0_ID_", (1:length(Name))))

suspect_WJ_v0.1.0 <- WormJam_v0.1.0 %>%
  select(name = idname, formula = Formula) %>%
  mutate(adduct = "[M+H]+")


unique_lipids_review <- read_tsv("unique_lipids.txt")

unique_metabolites_review <- read_tsv("unique_metabolites.txt")

unique_metabolipids_review <- rbind(unique_metabolites_review, unique_lipids_review) %>%
  mutate(idname = paste0("review_ID_", (1:length(Name))))


suspect_review <- unique_metabolipids_review %>%
  select(name = idname, formula = Formula) %>%
  mutate(adduct = "[M+H]+")

PubChemLite_database <- read_csv("C:/PL/PubChemLite_exposomics_20230224.csv") %>%   ## database used during patRoon elaboration
  mutate(idname = paste0("PCL_ID_", (1:length(Identifier))))

suspect_PCL <- PubChemLite_database %>%
  select(name = idname, formula = MolecularFormula) %>%
  mutate(adduct = "[M+H]+")

## running the suspect screening:

fGroupsSusp_WJplus <- screenSuspects(fGroups, suspect_WJplus) 

fGroupsSusp_WJ_zenodo <- screenSuspects(fGroups, suspect_WJ_zenodo)

fGroupsSusp_WJ_v0.1.0 <- screenSuspects(fGroups, suspect_WJ_v0.1.0)

fGroupsSusp_review <- screenSuspects(fGroups, suspect_review)

fGroupsSusp_PCL <- screenSuspects(fGroups, suspect_PCL)


## retrieving all the INFO from molecules:

suspscreen_WJplus <- as_tibble(fGroupsSusp_WJplus@screenInfo)
colnames(suspscreen_WJplus)[which(colnames(suspscreen_WJplus)=="name")] <- "idname"
if (all(is.na(pull(suspscreen_WJplus, "SMILES"))) & all(is.na(pull(suspscreen_WJplus, "InChI"))) & all(is.na(pull(suspscreen_WJplus, "InChIKey")))) {
  suspscreen_WJplus <- select(suspscreen_WJplus, -SMILES, -InChI, -InChIKey)
}

suspscreen_WJplus_joined <- left_join(x= suspscreen_WJplus, y = WJ_andPubChemRelMetabolites, by = "idname")




suspscreen_WJ_zenodo <- as_tibble(fGroupsSusp_WJ_zenodo@screenInfo)
colnames(suspscreen_WJ_zenodo)[which(colnames(suspscreen_WJ_zenodo)=="name")] <- "idname"
if (all(is.na(pull(suspscreen_WJ_zenodo, "SMILES"))) & all(is.na(pull(suspscreen_WJ_zenodo, "InChI"))) & all(is.na(pull(suspscreen_WJ_zenodo, "InChIKey")))) {
  suspscreen_WJ_zenodo <- select(suspscreen_WJ_zenodo, -SMILES, -InChI, -InChIKey)
}

suspscreen_WJ_zenodo_joined <- left_join(x= suspscreen_WJ_zenodo, y = WJ_zenodo, by = "idname")





suspscreen_WJ_v0.1.0 <- as_tibble(fGroupsSusp_WJ_v0.1.0@screenInfo)
colnames(suspscreen_WJ_v0.1.0)[which(colnames(suspscreen_WJ_v0.1.0)=="name")] <- "idname"
if (all(is.na(pull(suspscreen_WJ_v0.1.0, "SMILES"))) & all(is.na(pull(suspscreen_WJ_v0.1.0, "InChI"))) & all(is.na(pull(suspscreen_WJ_v0.1.0, "InChIKey")))) {
  suspscreen_WJ_v0.1.0 <- select(suspscreen_WJ_v0.1.0, -SMILES, -InChI, -InChIKey)
}

suspscreen_WJ_v0.1.0_joined <- left_join(x= suspscreen_WJ_v0.1.0, y = WormJam_v0.1.0, by = "idname")




suspscreen_review <- as_tibble(fGroupsSusp_review@screenInfo)
colnames(suspscreen_review)[which(colnames(suspscreen_review)=="name")] <- "idname"
if (all(is.na(pull(suspscreen_review, "SMILES"))) & all(is.na(pull(suspscreen_review, "InChI"))) & all(is.na(pull(suspscreen_review, "InChIKey")))) {
  suspscreen_review <- select(suspscreen_review, -SMILES, -InChI, -InChIKey)
}

suspscreen_review_joined <- left_join(x= suspscreen_review, y = unique_metabolipids_review, by = "idname")




suspscreen_PCL <- as_tibble(fGroupsSusp_PCL@screenInfo)
colnames(suspscreen_PCL)[which(colnames(suspscreen_PCL)=="name")] <- "idname"
if (all(is.na(pull(suspscreen_PCL, "SMILES"))) & all(is.na(pull(suspscreen_PCL, "InChI"))) & all(is.na(pull(suspscreen_PCL, "InChIKey")))) {
  suspscreen_PCL <- select(suspscreen_PCL, -SMILES, -InChI, -InChIKey)
}

suspscreen_PCL_joined <- left_join(x= suspscreen_PCL, y = PubChemLite_database, by = "idname")





## exporting the suspect screenings:

write_tsv(suspscreen_WJplus_joined, "suspscreen_HILICPOS_WJplus.txt")

write_tsv(suspscreen_WJ_zenodo_joined, "suspscreen_HILICPOS_WJ_zenodo.txt")

write_tsv(suspscreen_WJ_v0.1.0_joined, "suspscreen_HILICPOS_WJ_v0.1.0.txt")

write_tsv(suspscreen_review_joined, "suspscreen_HILICPOS_review.txt")

write_tsv(suspscreen_PCL_joined, "suspscreen_HILICPOS_PCL.txt")








##
#### 2024-04-10
######## Since there were several "Ignored following suspects for which no mass could be calculated"
#### updating adding more columns, if possible with the neutralMass


suspect_WJplus_updt <- WJ_andPubChemRelMetabolites %>%
  select(name = idname,
         neutralMass = MonoisotopicMass,
         formula = MolecularFormula,
         SMILES = SMILES,
         InChI = InChI) %>%
  mutate(adduct = "[M+H]+")



suspect_WJ_zenodo_updt <- WJ_zenodo %>%
  select(name = idname,
         #neutralMass = ,
         formula = Formula_Neutral,
         SMILES = SMILES_neutral,
         InChI = InChI_neutral) %>%
  mutate(adduct = "[M+H]+")


suspect_WJ_v0.1.0_updt <- WormJam_v0.1.0 %>%
  select(name = idname,
         #neutralMass = ,
         formula = Formula,
         SMILES = SMILES,
         InChI = InChI) %>%
  mutate(adduct = "[M+H]+")


suspect_review_updt <- unique_metabolipids_review %>%
  select(name = idname,
         #neutralMass = ,
         formula = Formula,
         SMILES = SMILES,
         InChI = InChI) %>%
  mutate(adduct = "[M+H]+")


suspect_PCL_updt <- PubChemLite_database %>%
  select(name = idname,
         neutralMass = MonoisotopicMass,
         formula = MolecularFormula,
         SMILES = SMILES,
         InChI = InChI) %>%
  mutate(adduct = "[M+H]+")



## running the suspect screening:

fGroupsSusp_WJplus_updt <- screenSuspects(fGroups, suspect_WJplus_updt) 

fGroupsSusp_WJ_zenodo_updt <- screenSuspects(fGroups, suspect_WJ_zenodo_updt)

fGroupsSusp_WJ_v0.1.0_updt <- screenSuspects(fGroups, suspect_WJ_v0.1.0_updt)

fGroupsSusp_review_updt <- screenSuspects(fGroups, suspect_review_updt)

fGroupsSusp_PCL_updt <- screenSuspects(fGroups, suspect_PCL_updt)


## retrieving all the INFO from molecules:

suspscreen_WJplus_updt <- as_tibble(fGroupsSusp_WJplus_updt@screenInfo)
colnames(suspscreen_WJplus_updt)[which(colnames(suspscreen_WJplus_updt)=="name")] <- "idname"


suspscreen_WJplus_joined_updt <- left_join(x= suspscreen_WJplus_updt, y = WJ_andPubChemRelMetabolites, by = "idname")




suspscreen_WJ_zenodo_updt <- as_tibble(fGroupsSusp_WJ_zenodo_updt@screenInfo)
colnames(suspscreen_WJ_zenodo_updt)[which(colnames(suspscreen_WJ_zenodo_updt)=="name")] <- "idname"

suspscreen_WJ_zenodo_joined_updt <- left_join(x= suspscreen_WJ_zenodo_updt, y = WJ_zenodo, by = "idname")





suspscreen_WJ_v0.1.0_updt <- as_tibble(fGroupsSusp_WJ_v0.1.0_updt@screenInfo)

colnames(suspscreen_WJ_v0.1.0_updt)[which(colnames(suspscreen_WJ_v0.1.0_updt)=="name")] <- "idname"


suspscreen_WJ_v0.1.0_joined_updt <- left_join(x= suspscreen_WJ_v0.1.0_updt, y = WormJam_v0.1.0, by = "idname")




suspscreen_review_updt <- as_tibble(fGroupsSusp_review_updt@screenInfo)

colnames(suspscreen_review_updt)[which(colnames(suspscreen_review_updt)=="name")] <- "idname"


suspscreen_review_joined_updt <- left_join(x= suspscreen_review_updt, y = unique_metabolipids_review, by = "idname")




suspscreen_PCL_updt <- as_tibble(fGroupsSusp_PCL_updt@screenInfo)

colnames(suspscreen_PCL_updt)[which(colnames(suspscreen_PCL_updt)=="name")] <- "idname"


suspscreen_PCL_joined_updt <- left_join(x= suspscreen_PCL_updt, y = PubChemLite_database, by = "idname")





## exporting the suspect screenings:

write_tsv(suspscreen_WJplus_joined_updt, "suspscreen_HILICPOS_WJplus_updt.txt")

write_tsv(suspscreen_WJ_zenodo_joined_updt, "suspscreen_HILICPOS_WJ_zenodo_updt.txt")

write_tsv(suspscreen_WJ_v0.1.0_joined_updt, "suspscreen_HILICPOS_WJ_v0.1.0_updt.txt")

write_tsv(suspscreen_review_joined_updt, "suspscreen_HILICPOS_review_updt.txt")

write_tsv(suspscreen_PCL_joined_updt, "suspscreen_HILICPOS_PCL_updt.txt")






## table summarising numbers:

suspscreen_summary <- tibble(database = c("WJplus",
                                          "WJ_zenodo",
                                          "WJ_v0.1.0",
                                          "review",
                                          "PCL"),
                             suspectsfound = c(length(pull(suspscreen_WJplus_joined, 1)),
                                               length(pull(suspscreen_WJ_zenodo_joined, 1)),
                                               length(pull(suspscreen_WJ_v0.1.0_joined, 1)),
                                               length(pull(suspscreen_review_joined, 1)),
                                               length(pull(suspscreen_PCL_joined, 1))),
                             suspectsfound_unique = c(length(unique(pull(suspscreen_WJplus_joined, 1))),
                                                      length(unique(pull(suspscreen_WJ_zenodo_joined, 1))),
                                                      length(unique(pull(suspscreen_WJ_v0.1.0_joined, 1))),
                                                      length(unique(pull(suspscreen_review_joined, 1))),
                                                      length(unique(pull(suspscreen_PCL_joined, 1)))),
                             total = c(length(pull(WJ_andPubChemRelMetabolites, 1)),
                                       length(pull(WJ_zenodo, 1)),
                                       length(pull(WormJam_v0.1.0, 1)),
                                       length(pull(unique_metabolipids_review, 1)),
                                       length(pull(PubChemLite_database, 1))))


write_tsv(suspscreen_summary, "suspscreen_summary_HILICPOS.txt")






## table summarising numbers:

suspscreen_summary_updt <- tibble(database = c("WJplus",
                                               "WJ_zenodo",
                                               "WJ_v0.1.0",
                                               "review",
                                               "PCL"),
                                  suspectsfound = c(length(pull(suspscreen_WJplus_joined_updt, 1)),
                                                    length(pull(suspscreen_WJ_zenodo_joined_updt, 1)),
                                                    length(pull(suspscreen_WJ_v0.1.0_joined_updt, 1)),
                                                    length(pull(suspscreen_review_joined_updt, 1)),
                                                    length(pull(suspscreen_PCL_joined_updt, 1))),
                                  suspectsfound_unique = c(length(unique(pull(suspscreen_WJplus_joined_updt, 1))),
                                                           length(unique(pull(suspscreen_WJ_zenodo_joined_updt, 1))),
                                                           length(unique(pull(suspscreen_WJ_v0.1.0_joined_updt, 1))),
                                                           length(unique(pull(suspscreen_review_joined_updt, 1))),
                                                           length(unique(pull(suspscreen_PCL_joined_updt, 1)))),
                                  total = c(length(pull(WJ_andPubChemRelMetabolites, 1)),
                                            length(pull(WJ_zenodo, 1)),
                                            length(pull(WormJam_v0.1.0, 1)),
                                            length(pull(unique_metabolipids_review, 1)),
                                            length(pull(PubChemLite_database, 1))))


write_tsv(suspscreen_summary_updt, "suspscreen_summary_HILICPOS_updt.txt")


## comparison:


suspscreen_summary_first_vs_updt <- tibble(database = c("WJplus",
                                                        "WJ_zenodo",
                                                        "WJ_v0.1.0",
                                                        "review",
                                                        "PCL"),
                                           unique_present_in_first_but_not_in_updt = c(sum(!(unique(pull(suspscreen_WJplus_joined, 1)) %in% unique(pull(suspscreen_WJplus_joined_updt, 1)))),
                                                                                       sum(!(unique(pull(suspscreen_WJ_zenodo_joined, 1)) %in% unique(pull(suspscreen_WJ_zenodo_joined_updt, 1)))),
                                                                                       sum(!(unique(pull(suspscreen_WJ_v0.1.0_joined, 1)) %in% unique(pull(suspscreen_WJ_v0.1.0_joined_updt, 1)))),
                                                                                       sum(!(unique(pull(suspscreen_review_joined, 1)) %in% unique(pull(suspscreen_review_joined_updt, 1)))),
                                                                                       sum(!(unique(pull(suspscreen_PCL_joined, 1)) %in% unique(pull(suspscreen_PCL_joined_updt, 1))))),
                                           unique_present_in_updt_but_not_in_first = c(sum(!(unique(pull(suspscreen_WJplus_joined_updt, 1)) %in% unique(pull(suspscreen_WJplus_joined, 1)))),
                                                                                       sum(!(unique(pull(suspscreen_WJ_zenodo_joined_updt, 1)) %in% unique(pull(suspscreen_WJ_zenodo_joined, 1)))),
                                                                                       sum(!(unique(pull(suspscreen_WJ_v0.1.0_joined_updt, 1)) %in% unique(pull(suspscreen_WJ_v0.1.0_joined, 1)))),
                                                                                       sum(!(unique(pull(suspscreen_review_joined_updt, 1)) %in% unique(pull(suspscreen_review_joined, 1)))),
                                                                                       sum(!(unique(pull(suspscreen_PCL_joined_updt, 1)) %in% unique(pull(suspscreen_PCL_joined, 1))))))

write_tsv(suspscreen_summary_first_vs_updt, "suspscreen_summary_first_vs_updt_HILICPOS.txt")





## 2024-11-28
### reperforming the annotation with the new WormJam integrated from Emma S.

WormJamDatabasePath_WJcombined <- "C:/WJ/20240829_C_elegans_all_Struct_Info.csv"

fGroups <- normInts(fGroups)


compounds_WJ_no_formula_combined <- generateCompounds(fGroups, mslists, "metfrag", method = "CL",
                                                                      dbRelMzDev = 5, fragRelMzDev = 5,
                                                                      fragAbsMzDev = 0.002, adduct = "[M+H]+",
                                                                      database = "csv", extraOpts = list(LocalDatabasePath = WormJamDatabasePath_WJcombined),
                                                                      scoreTypes = c("LiteratureCount", "PatentCount"),
                                                                      maxCandidatesToStop = 100)


# -------------------------
# reporting  - for wormjam_combined_fromEmmaAug2024
# -------------------------

reportCSV(fGroups, path = "report_WJcombined", formulas = NULL, compounds = compounds_WJ_no_formula_combined, components = components)

# Summary of MetFrag Results in a a Single Table 

MFsummary_WJ_no_formula_combined <- as.data.table(compounds_WJ_no_formula_combined)
write_tsv(MFsummary_WJ_no_formula_combined, "MFsummary_WJ_no_formula_combined.txt")







## 2024-11-29
### reperforming the annotation with the new WormJam integrated from Emma S.

WormJamDatabasePath_WJcombined <- "C:/WJ/20240829_C_elegans_all_Struct_Info.csv"

fGroups <- normInts(fGroups)


compounds_WJ_no_formula_combined_bis <- generateCompounds(fGroups, mslists, "metfrag", method = "CL",
                                                      dbRelMzDev = 5, fragRelMzDev = 5,
                                                      fragAbsMzDev = 0.002, adduct = "[M+H]+",
                                                      database = "csv", extraOpts = list(LocalDatabasePath = WormJamDatabasePath_WJcombined),
                                                      scoreTypes = c("fragScore","score", "individualMoNAScore",
                                                                     "LiteratureCount", "PatentCount"),
                                                      maxCandidatesToStop = 100)


# -------------------------
# reporting  - for wormjam_combined_fromEmmaAug2024
# -------------------------

reportCSV(fGroups, path = "report_WJcombined_bis", formulas = NULL, compounds = compounds_WJ_no_formula_combined_bis, components = components)

# Summary of MetFrag Results in a a Single Table 

MFsummary_WJ_no_formula_combined_bis <- as.data.table(compounds_WJ_no_formula_combined_bis)
write_tsv(MFsummary_WJ_no_formula_combined_bis, "MFsummary_WJ_no_formula_combined_bis.txt")


