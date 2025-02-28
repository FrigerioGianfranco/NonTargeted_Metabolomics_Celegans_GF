

library(tidyverse)

library(RChemMass)  # to get the Title of comounds
library(webchem)   # to get the CID from InChIKey from MS-DIAL data
library(classyfireR) # for the compound classification
library(networkD3)  # for the sankeyNetwork
library(UpSetR) # for the upset plot
library(gridExtra) # for tables close to box-plots
library(grid) # for tables close to box-plots
library(patchwork) # for combining the 2 plots
library(ggvenn) # for Eulero-venn with names
library(metaboliteIDmapping) # for converting from CID to KEGG or others
library(FELLA) # for enrichment analyses
library(igraph) # for printing FELLA graphs
library(venn) # for Eulero-venn with 5 sets

set.seed(1991)



samples_shipped_GF_mod <- read_tsv("samples_shipped_GF_mod.txt")


# changhing AnnoLevel with AnnotationLevel

AnnoLevel_to_AnnotationLevel <- function(a_df) {
  
  a_df_final <- a_df
  
  colnames(a_df_final)[which(colnames(a_df_final) == "AnnoLevel")] <- "AnnotationLevel"
  
  return(a_df_final)
}


# more info on the QC table

add_presence_in_scheme1or2_dry <- function(df_QC, df_QCsc1, df_QCsc2) {
  
  df_QC_output <- add_column(df_QC,
                             present_in_scheme1 = pull(df_QC, 1) %in% pull(df_QCsc1, 1),
                             present_in_scheme2 = pull(df_QC, 1) %in% pull(df_QCsc2, 1),
                             .after = "AnnotationLevel")
  
  if (!all(df_QC_output$present_in_scheme1 | df_QC_output$present_in_scheme2)) {stop("there are fearures not present not in scheme1 nor in scheme2")}
  
  return(df_QC_output)
  
}






# uploading all the tables with levels:

HILICPOS_MSDial_feat_info_annolevels <- read_tsv("HILICPOS_MSDial_feat_info_annolevels.txt") %>% AnnoLevel_to_AnnotationLevel()
HILICPOS_MSDial_feat_info_annolevels_QCfiltered <- read_tsv("HILICPOS_MSDial_feat_info_annolevels_QCfiltered.txt") %>% AnnoLevel_to_AnnotationLevel()
HILICPOS_MSDial_feat_info_annolevels_QCfiltered_scheme1_dry <- read_tsv("HILICPOS_MSDial_feat_info_annolevels_QCfiltered_scheme1_dry.txt") %>% AnnoLevel_to_AnnotationLevel()
HILICPOS_MSDial_feat_info_annolevels_QCfiltered_scheme2_dry <- read_tsv("HILICPOS_MSDial_feat_info_annolevels_QCfiltered_scheme2_dry.txt") %>% AnnoLevel_to_AnnotationLevel()
HILICPOS_MSDial_scheme1_dry_ANOVA_sign_annotat_withLevels <- read_tsv("HILICPOS_MSDial_scheme1_dry_ANOVA_sign_annotat_withLevels.txt") %>% AnnoLevel_to_AnnotationLevel()
HILICPOS_MSDial_scheme2_dry_ANOVA_sign_annotat_withLevels <- read_tsv("HILICPOS_MSDial_scheme2_dry_ANOVA_sign_annotat_withLevels.txt") %>% AnnoLevel_to_AnnotationLevel()
HILICPOS_MSDial_feat_info_annolevels_QCfiltered <- add_presence_in_scheme1or2_dry(HILICPOS_MSDial_feat_info_annolevels_QCfiltered, HILICPOS_MSDial_feat_info_annolevels_QCfiltered_scheme1_dry, HILICPOS_MSDial_feat_info_annolevels_QCfiltered_scheme2_dry)


HILICPOS_MSDial_WJ_feat_info_annolevels <- read_tsv("HILICPOS_MSDial_WJ_feat_info_annolevels.txt") %>% AnnoLevel_to_AnnotationLevel()
HILICPOS_MSDial_WJ_feat_info_annolevels_QCfiltered <- read_tsv("HILICPOS_MSDial_WJ_feat_info_annolevels_QCfiltered.txt") %>% AnnoLevel_to_AnnotationLevel()
HILICPOS_MSDial_WJ_feat_info_annolevels_QCfiltered_scheme1_dry <- read_tsv("HILICPOS_MSDial_WJ_feat_info_annolevels_QCfiltered_scheme1_dry.txt") %>% AnnoLevel_to_AnnotationLevel()
HILICPOS_MSDial_WJ_feat_info_annolevels_QCfiltered_scheme2_dry <- read_tsv("HILICPOS_MSDial_WJ_feat_info_annolevels_QCfiltered_scheme2_dry.txt") %>% AnnoLevel_to_AnnotationLevel()
HILICPOS_MSDial_WJ_scheme1_dry_ANOVA_sign_annotat_withLevels <- read_tsv("HILICPOS_MSDial_WJ_scheme1_dry_ANOVA_sign_annotat_withLevels.txt") %>% AnnoLevel_to_AnnotationLevel()
HILICPOS_MSDial_WJ_scheme2_dry_ANOVA_sign_annotat_withLevels <- read_tsv("HILICPOS_MSDial_WJ_scheme2_dry_ANOVA_sign_annotat_withLevels.txt") %>% AnnoLevel_to_AnnotationLevel()
HILICPOS_MSDial_WJ_feat_info_annolevels_QCfiltered <- add_presence_in_scheme1or2_dry(HILICPOS_MSDial_WJ_feat_info_annolevels_QCfiltered, HILICPOS_MSDial_WJ_feat_info_annolevels_QCfiltered_scheme1_dry, HILICPOS_MSDial_WJ_feat_info_annolevels_QCfiltered_scheme2_dry)



HILICPOS_patRoon_IPO_PCL_feat_info_annotation <- read_tsv("HILICPOS_patRoon_IPO_PCL_feat_info_annotation.txt") %>% AnnoLevel_to_AnnotationLevel()
HILICPOS_patRoon_IPO_PCL_feat_info_annotation_QCfiltered <- read_tsv("HILICPOS_patRoon_IPO_PCL_feat_info_annotation_QCfiltered.txt") %>% AnnoLevel_to_AnnotationLevel()
HILICPOS_patRoon_IPO_PCL_feat_info_annotation_QCfiltered_scheme1_dry <- read_tsv("HILICPOS_patRoon_IPO_PCL_feat_info_annotation_QCfiltered_scheme1_dry.txt") %>% AnnoLevel_to_AnnotationLevel()
HILICPOS_patRoon_IPO_PCL_feat_info_annotation_QCfiltered_scheme2_dry <- read_tsv("HILICPOS_patRoon_IPO_PCL_feat_info_annotation_QCfiltered_scheme2_dry.txt") %>% AnnoLevel_to_AnnotationLevel()
HILICPOS_patRoon_IPO_PCL_scheme1_dry_ANOVA_sign_annotat_withLevels <- read_tsv("HILICPOS_patRoon_IPO_PCL_scheme1_dry_ANOVA_sign_annotat_withLevels.txt") %>% AnnoLevel_to_AnnotationLevel()
HILICPOS_patRoon_IPO_PCL_scheme2_dry_ANOVA_sign_annotat_withLevels <- read_tsv("HILICPOS_patRoon_IPO_PCL_scheme2_dry_ANOVA_sign_annotat_withLevels.txt") %>% AnnoLevel_to_AnnotationLevel()
HILICPOS_patRoon_IPO_PCL_feat_info_annotation_QCfiltered <- add_presence_in_scheme1or2_dry(HILICPOS_patRoon_IPO_PCL_feat_info_annotation_QCfiltered, HILICPOS_patRoon_IPO_PCL_feat_info_annotation_QCfiltered_scheme1_dry, HILICPOS_patRoon_IPO_PCL_feat_info_annotation_QCfiltered_scheme2_dry)



HILICPOS_patRoon_IPO_WJ_feat_info_annotation <- read_tsv("HILICPOS_patRoon_IPO_WJ_feat_info_annotation.txt") %>% AnnoLevel_to_AnnotationLevel()
HILICPOS_patRoon_IPO_WJ_feat_info_annotation_QCfiltered <- read_tsv("HILICPOS_patRoon_IPO_WJ_feat_info_annotation_QCfiltered.txt") %>% AnnoLevel_to_AnnotationLevel()
HILICPOS_patRoon_IPO_WJ_feat_info_annotation_QCfiltered_scheme1_dry <- read_tsv("HILICPOS_patRoon_IPO_WJ_feat_info_annotation_QCfiltered_scheme1_dry.txt") %>% AnnoLevel_to_AnnotationLevel()
HILICPOS_patRoon_IPO_WJ_feat_info_annotation_QCfiltered_scheme2_dry <- read_tsv("HILICPOS_patRoon_IPO_WJ_feat_info_annotation_QCfiltered_scheme2_dry.txt") %>% AnnoLevel_to_AnnotationLevel()
HILICPOS_patRoon_IPO_WJ_scheme1_dry_ANOVA_sign_annotat_withLevels <- read_tsv("HILICPOS_patRoon_IPO_WJ_scheme1_dry_ANOVA_sign_annotat_withLevels.txt") %>% AnnoLevel_to_AnnotationLevel()
HILICPOS_patRoon_IPO_WJ_scheme2_dry_ANOVA_sign_annotat_withLevels <- read_tsv("HILICPOS_patRoon_IPO_WJ_scheme2_dry_ANOVA_sign_annotat_withLevels.txt") %>% AnnoLevel_to_AnnotationLevel()
HILICPOS_patRoon_IPO_WJ_feat_info_annotation_QCfiltered <- add_presence_in_scheme1or2_dry(HILICPOS_patRoon_IPO_WJ_feat_info_annotation_QCfiltered, HILICPOS_patRoon_IPO_WJ_feat_info_annotation_QCfiltered_scheme1_dry, HILICPOS_patRoon_IPO_WJ_feat_info_annotation_QCfiltered_scheme2_dry)






RPLCNEG_MSDial_feat_info_annolevels <- read_tsv("RPLCNEG_MSDial_feat_info_annolevels.txt") %>% AnnoLevel_to_AnnotationLevel()
RPLCNEG_MSDial_feat_info_annolevels_QCfiltered <- read_tsv("RPLCNEG_MSDial_feat_info_annolevels_QCfiltered.txt") %>% AnnoLevel_to_AnnotationLevel()
RPLCNEG_MSDial_feat_info_annolevels_QCfiltered_scheme1_dry <- read_tsv("RPLCNEG_MSDial_feat_info_annolevels_QCfiltered_scheme1_dry.txt") %>% AnnoLevel_to_AnnotationLevel()
RPLCNEG_MSDial_feat_info_annolevels_QCfiltered_scheme2_dry <- read_tsv("RPLCNEG_MSDial_feat_info_annolevels_QCfiltered_scheme2_dry.txt") %>% AnnoLevel_to_AnnotationLevel()
RPLCNEG_MSDial_scheme1_dry_ANOVA_sign_annotat_withLevels <- read_tsv("RPLCNEG_MSDial_scheme1_dry_ANOVA_sign_annotat_withLevels.txt") %>% AnnoLevel_to_AnnotationLevel()
RPLCNEG_MSDial_scheme2_dry_ANOVA_sign_annotat_withLevels <- read_tsv("RPLCNEG_MSDial_scheme2_dry_ANOVA_sign_annotat_withLevels.txt") %>% AnnoLevel_to_AnnotationLevel()
RPLCNEG_MSDial_feat_info_annolevels_QCfiltered <- add_presence_in_scheme1or2_dry(RPLCNEG_MSDial_feat_info_annolevels_QCfiltered, RPLCNEG_MSDial_feat_info_annolevels_QCfiltered_scheme1_dry, RPLCNEG_MSDial_feat_info_annolevels_QCfiltered_scheme2_dry)



RPLCNEG_MSDial_WJ_feat_info_annolevels <- read_tsv("RPLCNEG_MSDial_WJ_feat_info_annolevels.txt") %>% AnnoLevel_to_AnnotationLevel()
RPLCNEG_MSDial_WJ_feat_info_annolevels_QCfiltered <- read_tsv("RPLCNEG_MSDial_WJ_feat_info_annolevels_QCfiltered.txt") %>% AnnoLevel_to_AnnotationLevel()
RPLCNEG_MSDial_WJ_feat_info_annolevels_QCfiltered_scheme1_dry <- read_tsv("RPLCNEG_MSDial_WJ_feat_info_annolevels_QCfiltered_scheme1_dry.txt") %>% AnnoLevel_to_AnnotationLevel()
RPLCNEG_MSDial_WJ_feat_info_annolevels_QCfiltered_scheme2_dry <- read_tsv("RPLCNEG_MSDial_WJ_feat_info_annolevels_QCfiltered_scheme2_dry.txt") %>% AnnoLevel_to_AnnotationLevel()
RPLCNEG_MSDial_WJ_scheme1_dry_ANOVA_sign_annotat_withLevels <- read_tsv("RPLCNEG_MSDial_WJ_scheme1_dry_ANOVA_sign_annotat_withLevels.txt") %>% AnnoLevel_to_AnnotationLevel()
RPLCNEG_MSDial_WJ_scheme2_dry_ANOVA_sign_annotat_withLevels <- read_tsv("RPLCNEG_MSDial_WJ_scheme2_dry_ANOVA_sign_annotat_withLevels.txt") %>% AnnoLevel_to_AnnotationLevel()
RPLCNEG_MSDial_WJ_feat_info_annolevels_QCfiltered <- add_presence_in_scheme1or2_dry(RPLCNEG_MSDial_WJ_feat_info_annolevels_QCfiltered, RPLCNEG_MSDial_WJ_feat_info_annolevels_QCfiltered_scheme1_dry, RPLCNEG_MSDial_WJ_feat_info_annolevels_QCfiltered_scheme2_dry)




RPLCNEG_patRoon_IPO_PCL_feat_info_annotation <- read_tsv("RPLCNEG_patRoon_IPO_PCL_feat_info_annotation.txt") %>% AnnoLevel_to_AnnotationLevel()
RPLCNEG_patRoon_IPO_PCL_feat_info_annotation_QCfiltered <- read_tsv("RPLCNEG_patRoon_IPO_PCL_feat_info_annotation_QCfiltered.txt") %>% AnnoLevel_to_AnnotationLevel()
RPLCNEG_patRoon_IPO_PCL_feat_info_annotation_QCfiltered_scheme1_dry <- read_tsv("RPLCNEG_patRoon_IPO_PCL_feat_info_annotation_QCfiltered_scheme1_dry.txt") %>% AnnoLevel_to_AnnotationLevel()
RPLCNEG_patRoon_IPO_PCL_feat_info_annotation_QCfiltered_scheme2_dry <- read_tsv("RPLCNEG_patRoon_IPO_PCL_feat_info_annotation_QCfiltered_scheme2_dry.txt") %>% AnnoLevel_to_AnnotationLevel()
RPLCNEG_patRoon_IPO_PCL_scheme1_dry_ANOVA_sign_annot_Levels <- read_tsv("RPLCNEG_patRoon_IPO_PCL_scheme1_dry_ANOVA_sign_annot_Levels.txt") %>% AnnoLevel_to_AnnotationLevel()
RPLCNEG_patRoon_IPO_PCL_scheme2_dry_ANOVA_sign_annotat_withLevels <- read_tsv("RPLCNEG_patRoon_IPO_PCL_scheme2_dry_ANOVA_sign_annotat_withLevels.txt") %>% AnnoLevel_to_AnnotationLevel()
RPLCNEG_patRoon_IPO_PCL_feat_info_annotation_QCfiltered <- add_presence_in_scheme1or2_dry(RPLCNEG_patRoon_IPO_PCL_feat_info_annotation_QCfiltered, RPLCNEG_patRoon_IPO_PCL_feat_info_annotation_QCfiltered_scheme1_dry, RPLCNEG_patRoon_IPO_PCL_feat_info_annotation_QCfiltered_scheme2_dry)


RPLCNEG_patRoon_IPO_WJ_feat_info_annotation <- read_tsv("RPLCNEG_patRoon_IPO_WJ_feat_info_annotation.txt") %>% AnnoLevel_to_AnnotationLevel()
RPLCNEG_patRoon_IPO_WJ_feat_info_annotation_QCfiltered <- read_tsv("RPLCNEG_patRoon_IPO_WJ_feat_info_annotation_QCfiltered.txt") %>% AnnoLevel_to_AnnotationLevel()
RPLCNEG_patRoon_IPO_WJ_feat_info_annotation_QCfiltered_scheme1_dry <- read_tsv("RPLCNEG_patRoon_IPO_WJ_feat_info_annotation_QCfiltered_scheme1_dry.txt") %>% AnnoLevel_to_AnnotationLevel()
RPLCNEG_patRoon_IPO_WJ_feat_info_annotation_QCfiltered_scheme2_dry <- read_tsv("RPLCNEG_patRoon_IPO_WJ_feat_info_annotation_QCfiltered_scheme2_dry.txt") %>% AnnoLevel_to_AnnotationLevel()
RPLCNEG_patRoon_IPO_WJ_scheme1_dry_ANOVA_sign_annotat_withLevels <- read_tsv("RPLCNEG_patRoon_IPO_WJ_scheme1_dry_ANOVA_sign_annotat_withLevels.txt") %>% AnnoLevel_to_AnnotationLevel()
RPLCNEG_patRoon_IPO_WJ_scheme2_dry_ANOVA_sign_annotat_withLevels <- read_tsv("RPLCNEG_patRoon_IPO_WJ_scheme2_dry_ANOVA_sign_annotat_withLevels.txt") %>% AnnoLevel_to_AnnotationLevel()
RPLCNEG_patRoon_IPO_WJ_feat_info_annotation_QCfiltered <- add_presence_in_scheme1or2_dry(RPLCNEG_patRoon_IPO_WJ_feat_info_annotation_QCfiltered, RPLCNEG_patRoon_IPO_WJ_feat_info_annotation_QCfiltered_scheme1_dry, RPLCNEG_patRoon_IPO_WJ_feat_info_annotation_QCfiltered_scheme2_dry)





##data: 


MSDIAL_HLP4_scheme1_dry_transf <- read_tsv("MSDIAL_HLP4_shem1_dry_transf.txt")
MSDIAL_HLP4_scheme1_dry_transf$Sample_ID <- factor(MSDIAL_HLP4_scheme1_dry_transf$Sample_ID, levels = c("N2", "VC40", "VC1668", "UA57", "BR5270"))

MSDIAL_HLP4_scheme2_dry_transf <- read_tsv("MSDIAL_HLP4_shem2_dry_transf.txt")
MSDIAL_HLP4_scheme2_dry_transf$Sample_ID <- factor(MSDIAL_HLP4_scheme2_dry_transf$Sample_ID, levels = c("N2", "VC40", "VC1668", "UA57", "BR5270"))



MSDIAL_WJ_HLP4_scheme1_dry_transf <- read_tsv("MSDIAL_WJ_HLP4_shem1_dry_transf.txt")
MSDIAL_WJ_HLP4_scheme1_dry_transf$Sample_ID <- factor(MSDIAL_WJ_HLP4_scheme1_dry_transf$Sample_ID, levels = c("N2", "VC40", "VC1668", "UA57", "BR5270"))

MSDIAL_WJ_HLP4_scheme2_dry_transf <- read_tsv("MSDIAL_WJ_HLP4_shem2_dry_transf.txt")
MSDIAL_WJ_HLP4_scheme2_dry_transf$Sample_ID <- factor(MSDIAL_WJ_HLP4_scheme2_dry_transf$Sample_ID, levels = c("N2", "VC40", "VC1668", "UA57", "BR5270"))


patRoon_HLP4_shem1_dry_transf <- read_tsv("HLP4_shem1_dry_transf.txt")
patRoon_HLP4_shem1_dry_transf$Sample_ID <- factor(patRoon_HLP4_shem1_dry_transf$Sample_ID, levels = c("N2", "VC40", "VC1668", "UA57", "BR5270"))

patRoon_HLP4_shem2_dry_transf <- read_tsv("HLP4_shem2_dry_transf.txt")
patRoon_HLP4_shem2_dry_transf$Sample_ID <- factor(patRoon_HLP4_shem2_dry_transf$Sample_ID, levels = c("N2", "VC40", "VC1668", "UA57", "BR5270"))






MSDIAL_RPN4_scheme1_dry_transf <- read_tsv("MSDIAL_RPN4_scheme1_dry_transf.txt")
MSDIAL_RPN4_scheme1_dry_transf$Sample_ID <- factor(MSDIAL_RPN4_scheme1_dry_transf$Sample_ID, levels = c("N2", "VC40", "VC1668", "UA57", "BR5270"))

MSDIAL_RPN4_scheme2_dry_transf <- read_tsv("MSDIAL_RPN4_scheme2_dry_transf.txt")
MSDIAL_RPN4_scheme2_dry_transf$Sample_ID <- factor(MSDIAL_RPN4_scheme2_dry_transf$Sample_ID, levels = c("N2", "VC40", "VC1668", "UA57", "BR5270"))


MSDIAL_WJ_RPN4_scheme1_dry_transf <- read_tsv("MSDIAL_WJ_RPN4_scheme1_dry_transf.txt")
MSDIAL_WJ_RPN4_scheme1_dry_transf$Sample_ID <- factor(MSDIAL_WJ_RPN4_scheme1_dry_transf$Sample_ID, levels = c("N2", "VC40", "VC1668", "UA57", "BR5270"))

MSDIAL_WJ_RPN4_scheme2_dry_transf <- read_tsv("MSDIAL_WJ_RPN4_scheme2_dry_transf.txt")
MSDIAL_WJ_RPN4_scheme2_dry_transf$Sample_ID <- factor(MSDIAL_WJ_RPN4_scheme2_dry_transf$Sample_ID, levels = c("N2", "VC40", "VC1668", "UA57", "BR5270"))


patRoon_RPN4_scheme1_dry_transf <- read_tsv("RPN4_scheme1_dry_transf.txt")
patRoon_RPN4_scheme1_dry_transf$Sample_ID <- factor(patRoon_RPN4_scheme1_dry_transf$Sample_ID, levels = c("N2", "VC40", "VC1668", "UA57", "BR5270"))

patRoon_RPN4_scheme2_dry_transf <- read_tsv("RPN4_scheme2_dry_transf.txt")
patRoon_RPN4_scheme2_dry_transf$Sample_ID <- factor(patRoon_RPN4_scheme2_dry_transf$Sample_ID, levels = c("N2", "VC40", "VC1668", "UA57", "BR5270"))





# Getting all the summaries of annotation levels

all_data_list <- list(HILICPOS_MSDial_feat_info_annolevels = HILICPOS_MSDial_feat_info_annolevels,
                      HILICPOS_MSDial_feat_info_annolevels_QCfiltered = HILICPOS_MSDial_feat_info_annolevels_QCfiltered,
                      HILICPOS_MSDial_feat_info_annolevels_QCfiltered_scheme1_dry = HILICPOS_MSDial_feat_info_annolevels_QCfiltered_scheme1_dry,
                      HILICPOS_MSDial_feat_info_annolevels_QCfiltered_scheme2_dry = HILICPOS_MSDial_feat_info_annolevels_QCfiltered_scheme2_dry,
                      HILICPOS_MSDial_scheme1_dry_ANOVA_sign_annotat_withLevels = HILICPOS_MSDial_scheme1_dry_ANOVA_sign_annotat_withLevels,
                      HILICPOS_MSDial_scheme2_dry_ANOVA_sign_annotat_withLevels = HILICPOS_MSDial_scheme2_dry_ANOVA_sign_annotat_withLevels,
                      
                      HILICPOS_MSDial_WJ_feat_info_annolevels = HILICPOS_MSDial_WJ_feat_info_annolevels,
                      HILICPOS_MSDial_WJ_feat_info_annolevels_QCfiltered = HILICPOS_MSDial_WJ_feat_info_annolevels_QCfiltered,
                      HILICPOS_MSDial_WJ_feat_info_annolevels_QCfiltered_scheme1_dry = HILICPOS_MSDial_WJ_feat_info_annolevels_QCfiltered_scheme1_dry,
                      HILICPOS_MSDial_WJ_feat_info_annolevels_QCfiltered_scheme2_dry = HILICPOS_MSDial_WJ_feat_info_annolevels_QCfiltered_scheme2_dry,
                      HILICPOS_MSDial_WJ_scheme1_dry_ANOVA_sign_annotat_withLevels = HILICPOS_MSDial_WJ_scheme1_dry_ANOVA_sign_annotat_withLevels,
                      HILICPOS_MSDial_WJ_scheme2_dry_ANOVA_sign_annotat_withLevels = HILICPOS_MSDial_WJ_scheme2_dry_ANOVA_sign_annotat_withLevels,
                      
                      HILICPOS_patRoon_IPO_PCL_feat_info_annotation = HILICPOS_patRoon_IPO_PCL_feat_info_annotation,
                      HILICPOS_patRoon_IPO_PCL_feat_info_annotation_QCfiltered = HILICPOS_patRoon_IPO_PCL_feat_info_annotation_QCfiltered,
                      HILICPOS_patRoon_IPO_PCL_feat_info_annotation_QCfiltered_scheme1_dry = HILICPOS_patRoon_IPO_PCL_feat_info_annotation_QCfiltered_scheme1_dry,
                      HILICPOS_patRoon_IPO_PCL_feat_info_annotation_QCfiltered_scheme2_dry = HILICPOS_patRoon_IPO_PCL_feat_info_annotation_QCfiltered_scheme2_dry,
                      HILICPOS_patRoon_IPO_PCL_scheme1_dry_ANOVA_sign_annotat_withLevels = HILICPOS_patRoon_IPO_PCL_scheme1_dry_ANOVA_sign_annotat_withLevels,
                      HILICPOS_patRoon_IPO_PCL_scheme2_dry_ANOVA_sign_annotat_withLevels = HILICPOS_patRoon_IPO_PCL_scheme2_dry_ANOVA_sign_annotat_withLevels,
                      
                      HILICPOS_patRoon_IPO_WJ_feat_info_annotation = HILICPOS_patRoon_IPO_WJ_feat_info_annotation,
                      HILICPOS_patRoon_IPO_WJ_feat_info_annotation_QCfiltered = HILICPOS_patRoon_IPO_WJ_feat_info_annotation_QCfiltered,
                      HILICPOS_patRoon_IPO_WJ_feat_info_annotation_QCfiltered_scheme1_dry = HILICPOS_patRoon_IPO_WJ_feat_info_annotation_QCfiltered_scheme1_dry,
                      HILICPOS_patRoon_IPO_WJ_feat_info_annotation_QCfiltered_scheme2_dry = HILICPOS_patRoon_IPO_WJ_feat_info_annotation_QCfiltered_scheme2_dry,
                      HILICPOS_patRoon_IPO_WJ_scheme1_dry_ANOVA_sign_annotat_withLevels = HILICPOS_patRoon_IPO_WJ_scheme1_dry_ANOVA_sign_annotat_withLevels,
                      HILICPOS_patRoon_IPO_WJ_scheme2_dry_ANOVA_sign_annotat_withLevels = HILICPOS_patRoon_IPO_WJ_scheme2_dry_ANOVA_sign_annotat_withLevels,
                      
                      RPLCNEG_MSDial_feat_info_annolevels = RPLCNEG_MSDial_feat_info_annolevels,
                      RPLCNEG_MSDial_feat_info_annolevels_QCfiltered = RPLCNEG_MSDial_feat_info_annolevels_QCfiltered,
                      RPLCNEG_MSDial_feat_info_annolevels_QCfiltered_scheme1_dry = RPLCNEG_MSDial_feat_info_annolevels_QCfiltered_scheme1_dry,
                      RPLCNEG_MSDial_feat_info_annolevels_QCfiltered_scheme2_dry = RPLCNEG_MSDial_feat_info_annolevels_QCfiltered_scheme2_dry,
                      RPLCNEG_MSDial_scheme1_dry_ANOVA_sign_annotat_withLevels = RPLCNEG_MSDial_scheme1_dry_ANOVA_sign_annotat_withLevels,
                      RPLCNEG_MSDial_scheme2_dry_ANOVA_sign_annotat_withLevels = RPLCNEG_MSDial_scheme2_dry_ANOVA_sign_annotat_withLevels,
                      
                      RPLCNEG_MSDial_WJ_feat_info_annolevels = RPLCNEG_MSDial_WJ_feat_info_annolevels,
                      RPLCNEG_MSDial_WJ_feat_info_annolevels_QCfiltered = RPLCNEG_MSDial_WJ_feat_info_annolevels_QCfiltered,
                      RPLCNEG_MSDial_WJ_feat_info_annolevels_QCfiltered_scheme1_dry = RPLCNEG_MSDial_WJ_feat_info_annolevels_QCfiltered_scheme1_dry,
                      RPLCNEG_MSDial_WJ_feat_info_annolevels_QCfiltered_scheme2_dry = RPLCNEG_MSDial_WJ_feat_info_annolevels_QCfiltered_scheme2_dry,
                      RPLCNEG_MSDial_WJ_scheme1_dry_ANOVA_sign_annotat_withLevels = RPLCNEG_MSDial_WJ_scheme1_dry_ANOVA_sign_annotat_withLevels,
                      RPLCNEG_MSDial_WJ_scheme2_dry_ANOVA_sign_annotat_withLevels = RPLCNEG_MSDial_WJ_scheme2_dry_ANOVA_sign_annotat_withLevels,
                      
                      
                      RPLCNEG_patRoon_IPO_PCL_feat_info_annotation = RPLCNEG_patRoon_IPO_PCL_feat_info_annotation,
                      RPLCNEG_patRoon_IPO_PCL_feat_info_annotation_QCfiltered = RPLCNEG_patRoon_IPO_PCL_feat_info_annotation_QCfiltered,
                      RPLCNEG_patRoon_IPO_PCL_feat_info_annotation_QCfiltered_scheme1_dry = RPLCNEG_patRoon_IPO_PCL_feat_info_annotation_QCfiltered_scheme1_dry,
                      RPLCNEG_patRoon_IPO_PCL_feat_info_annotation_QCfiltered_scheme2_dry = RPLCNEG_patRoon_IPO_PCL_feat_info_annotation_QCfiltered_scheme2_dry,
                      RPLCNEG_patRoon_IPO_PCL_scheme1_dry_ANOVA_sign_annot_Levels = RPLCNEG_patRoon_IPO_PCL_scheme1_dry_ANOVA_sign_annot_Levels,
                      RPLCNEG_patRoon_IPO_PCL_scheme2_dry_ANOVA_sign_annotat_withLevels = RPLCNEG_patRoon_IPO_PCL_scheme2_dry_ANOVA_sign_annotat_withLevels,
                      
                      RPLCNEG_patRoon_IPO_WJ_feat_info_annotation = RPLCNEG_patRoon_IPO_WJ_feat_info_annotation,
                      RPLCNEG_patRoon_IPO_WJ_feat_info_annotation_QCfiltered = RPLCNEG_patRoon_IPO_WJ_feat_info_annotation_QCfiltered,
                      RPLCNEG_patRoon_IPO_WJ_feat_info_annotation_QCfiltered_scheme1_dry = RPLCNEG_patRoon_IPO_WJ_feat_info_annotation_QCfiltered_scheme1_dry,
                      RPLCNEG_patRoon_IPO_WJ_feat_info_annotation_QCfiltered_scheme2_dry = RPLCNEG_patRoon_IPO_WJ_feat_info_annotation_QCfiltered_scheme2_dry,
                      RPLCNEG_patRoon_IPO_WJ_scheme1_dry_ANOVA_sign_annotat_withLevels = RPLCNEG_patRoon_IPO_WJ_scheme1_dry_ANOVA_sign_annotat_withLevels,
                      RPLCNEG_patRoon_IPO_WJ_scheme2_dry_ANOVA_sign_annotat_withLevels = RPLCNEG_patRoon_IPO_WJ_scheme2_dry_ANOVA_sign_annotat_withLevels)


all_data_summary_list <- vector("list", length(all_data_list))

for (i in 1:length(all_data_list))  {
  
  DF <- all_data_list[[i]]
  DF_name <- names(all_data_list)[i]
  
  DF_summary <- DF %>%
    group_by(AnnotationLevel) %>%
    summarize(N=n()) %>%
    mutate(perc = round(N/sum(N)*100,1))
  
  DF_summary_totals <- tibble(AnnotationLevel = "TOTAL:",
                              N = sum(DF_summary$N),
                              perc = sum(DF_summary$perc))
  
  
  if (DF_summary_totals$N != length(pull(DF,colnames(DF)[1]))) {stop("Something wrong: total number of features are not matching!")}
  
  DF_summary_final <- rbind(DF_summary, DF_summary_totals)
  
  write_csv(DF_summary_final, paste0("AnnLEV_", DF_name, ".csv"))
  
  assign(paste0("AnnLEV_", DF_name), DF_summary_final, envir = .GlobalEnv)
  
  all_data_summary_list[[i]] <- DF_summary_final
  names(all_data_summary_list)[i] <- paste0("AnnLEV_", DF_name)
  
}







# Retrieving details for levels 2 and 3:

all_lev2_3_list_QC <- vector("list", length(names(all_data_list)[which(grepl("QCfiltered", names(all_data_list)) & !grepl("scheme", names(all_data_list)))]))

for (i in 1:length(all_data_list))  {

  DF <-  all_data_list[[i]]
  DF_name <- names(all_data_list)[i]
  
  if (grepl("QCfiltered", DF_name) & !grepl("scheme", DF_name)) {
    DF_2_3 <- filter(DF, AnnotationLevel %in% c("2a", "2b", "3a", "3b", "3c"))
    
    if (length(DF_2_3$AnnotationLevel) != 0) {
      if(grepl("MSDial", DF_name)) {
        
        DF_2_3_sel <- select(DF_2_3,
                             all_of(c("Alignment_ID",  "Average_Rt.min_", "Average_Mz", "AnnotationLevel", "present_in_scheme1", "present_in_scheme2",
                                      "Metabolite_name", "Formula", "INCHIKEY", "SMILES" )))
        
        colnames(DF_2_3_sel) <- str_replace_all(colnames(DF_2_3_sel), "Alignment_ID", "feature")
        colnames(DF_2_3_sel) <- str_replace_all(colnames(DF_2_3_sel), "Average_Rt.min_", "rt")
        colnames(DF_2_3_sel) <- str_replace_all(colnames(DF_2_3_sel), "Average_Mz", "mz")
        colnames(DF_2_3_sel) <- str_replace_all(colnames(DF_2_3_sel), "Metabolite_name", "compoundName")
        colnames(DF_2_3_sel) <- str_replace_all(colnames(DF_2_3_sel), "Formula", "neutral_formula")
        colnames(DF_2_3_sel) <- str_replace_all(colnames(DF_2_3_sel), "INCHIKEY", "InChIKey")
        
        TIBBLE_WITH_CID <- get_cid(DF_2_3_sel$InChIKey, from = "inchikey", match = "first")
        
        DF_2_3_sel <- add_column(DF_2_3_sel,
                                  .after = "compoundName",
                                  identifier = as.character(TIBBLE_WITH_CID$cid))
        
        DF_2_3_sel <- mutate(DF_2_3_sel, identifier = ifelse(is.na(identifier), "null", identifier))
        
        
        write_csv(DF_2_3_sel, paste0("AnnLEV2_3_", DF_name, ".csv"))
        
        assign(paste0("AnnLEV2_3_", DF_name), DF_2_3_sel, envir = .GlobalEnv)
        
        all_lev2_3_list_QC[[i]] <- DF_2_3_sel
        names(all_lev2_3_list_QC)[i] <- paste0("AnnLEV2_3_", DF_name)
        
        
      } else if (grepl("patRoon", DF_name)) {
        
        
        DF_2_3_sel <- select(DF_2_3,
                             all_of(c("group", "ret", "mz", "AnnotationLevel", "present_in_scheme1", "present_in_scheme2",
                                      "compoundName_topMONA", "identifier_topMONA", "neutral_formula_topMONA", "InChIKey_topMONA", "SMILES_topMONA")))
        
        
        colnames(DF_2_3_sel) <- str_replace_all(colnames(DF_2_3_sel), "group", "feature")
        colnames(DF_2_3_sel) <- str_replace_all(colnames(DF_2_3_sel), "ret", "rt")
        colnames(DF_2_3_sel) <- str_remove_all(colnames(DF_2_3_sel), "_topMONA")
        
        
        DF_2_3_sel$identifier <- as.character(DF_2_3_sel$identifier)
        
        DF_2_3_sel <- mutate(DF_2_3_sel, identifier = ifelse(is.na(identifier), "null", identifier))
        
        write_csv(DF_2_3_sel, paste0("AnnLEV2_3_", DF_name, ".csv"))
        
        
        assign(paste0("AnnLEV2_3_", DF_name), DF_2_3_sel, envir = .GlobalEnv)
        
        all_lev2_3_list_QC[[i]] <- DF_2_3_sel
        names(all_lev2_3_list_QC)[i] <- paste0("AnnLEV2_3_", DF_name)
        
        
        
      }
    }
  }
}




## creating a single table with all 2_3 tables

all_lev2_3_df_QC <- bind_rows(all_lev2_3_list_QC, .id = "type_analysis")

write_tsv(all_lev2_3_df_QC, "all_lev2_3_df_QC.txt")




unique_identifier_2_3_QC <- unique(pull(all_lev2_3_df_QC, "identifier"))


## creating separate table with info of the 2_3 compounds for each 

single_LEV2_3_list_QC <- vector("list", length(unique_identifier_2_3_QC))

for (a in unique_identifier_2_3_QC) {
 
  nametable <- paste0("singleLEV2_3_CID_", a)
  
  df_fil <- all_lev2_3_df_QC %>%
    filter(identifier == a)
  
  The_Title <- getPCdesc.title(query = a, from = "cid", timeout=10)$Title
  
  The_classification <- get_classification(df_fil$InChIKey[1])
  
  df_fil <- add_column(df_fil,
                       .after = "identifier",
                       Title = rep(The_Title, length(pull(df_fil, colnames(df_fil)[1])))) %>%
    add_column(.after = "Title",
               kingdom = ifelse(!is.null(The_classification),
                                rep(The_classification@classification[["Classification"]][1], length(pull(df_fil, colnames(df_fil)[1]))),
                                rep(NA, length(pull(df_fil, colnames(df_fil)[1])))),
               superclass = ifelse(!is.null(The_classification),
                                   rep(The_classification@classification[["Classification"]][2], length(pull(df_fil, colnames(df_fil)[1]))),
                                   rep(NA, length(pull(df_fil, colnames(df_fil)[1])))),
               class = ifelse(!is.null(The_classification),
                              rep(The_classification@classification[["Classification"]][3], length(pull(df_fil, colnames(df_fil)[1]))),
                              rep(NA, length(pull(df_fil, colnames(df_fil)[1])))),
               subclass = ifelse(!is.null(The_classification),
                                 rep(The_classification@classification[["Classification"]][4], length(pull(df_fil, colnames(df_fil)[1]))),
                                 rep(NA, length(pull(df_fil, colnames(df_fil)[1])))),
               level5 = ifelse(!is.null(The_classification),
                               rep(The_classification@classification[["Classification"]][5], length(pull(df_fil, colnames(df_fil)[1]))),
                               rep(NA, length(pull(df_fil, colnames(df_fil)[1])))),
               level6 = ifelse(!is.null(The_classification),
                               rep(The_classification@classification[["Classification"]][6], length(pull(df_fil, colnames(df_fil)[1]))),
                               rep(NA, length(pull(df_fil, colnames(df_fil)[1])))),
               level7 = ifelse(!is.null(The_classification),
                               rep(The_classification@classification[["Classification"]][7], length(pull(df_fil, colnames(df_fil)[1]))),
                               rep(NA, length(pull(df_fil, colnames(df_fil)[1])))))
  
  
  
  df_fil <- relocate(df_fil,
                     feature, rt, mz, type_analysis, AnnotationLevel,
                     present_in_scheme1, present_in_scheme2,
                     compoundName, identifier, neutral_formula, InChIKey, SMILES,
                     Title, kingdom, superclass, class, subclass, level5, level6, level7)
  
  if (length(unique(df_fil$AnnotationLevel))!=1) {
    warning(paste0("there are different levels for the identifier ", a))
  }
  
  write_csv(df_fil, paste0(nametable, ".csv"))
  
  assign(nametable, df_fil, envir = .GlobalEnv)
  
  single_LEV2_3_list_QC[[which(unique_identifier_2_3_QC==a)]] <- df_fil
  names(single_LEV2_3_list_QC)[which(unique_identifier_2_3_QC==a)] <- nametable
  
}



single_LEV2_3_df_QC <- bind_rows(single_LEV2_3_list_QC)


single_LEV2_3_df_QC$type_analysis <- str_replace_all(single_LEV2_3_df_QC$type_analysis, "AnnLEV2_3_HILICPOS_MSDial_feat_info_annolevels_QCfiltered", "HILICPOS_MSDial_PublicMSP")
single_LEV2_3_df_QC$type_analysis <- str_replace_all(single_LEV2_3_df_QC$type_analysis, "AnnLEV2_3_HILICPOS_MSDial_WJ_feat_info_annolevels_QCfiltered", "HILICPOS_MSDial_WormJamExpanded")
single_LEV2_3_df_QC$type_analysis <- str_replace_all(single_LEV2_3_df_QC$type_analysis, "AnnLEV2_3_HILICPOS_patRoon_IPO_PCL_feat_info_annotation_QCfiltered", "HILICPOS_patRoon_PubChemLite")
single_LEV2_3_df_QC$type_analysis <- str_replace_all(single_LEV2_3_df_QC$type_analysis, "AnnLEV2_3_HILICPOS_patRoon_IPO_WJ_feat_info_annotation_QCfiltered", "HILICPOS_patRoon_WormJamExpanded")
single_LEV2_3_df_QC$type_analysis <- str_replace_all(single_LEV2_3_df_QC$type_analysis, "AnnLEV2_3_RPLCNEG_MSDial_feat_info_annolevels_QCfiltered", "RPLCNEG_MSDial_PublicMSP")
single_LEV2_3_df_QC$type_analysis <- str_replace_all(single_LEV2_3_df_QC$type_analysis, "AnnLEV2_3_RPLCNEG_MSDial_WJ_feat_info_annolevels_QCfiltered", "RPLCNEG_MSDial_WormJamExpanded")
single_LEV2_3_df_QC$type_analysis <- str_replace_all(single_LEV2_3_df_QC$type_analysis, "AnnLEV2_3_RPLCNEG_patRoon_IPO_PCL_feat_info_annotation_QCfiltered", "RPLCNEG_patRoon_PubChemLite")
single_LEV2_3_df_QC$type_analysis <- str_replace_all(single_LEV2_3_df_QC$type_analysis, "AnnLEV2_3_RPLCNEG_patRoon_IPO_WJ_feat_info_annotation_QCfiltered", "RPLCNEG_patRoon_WormJamExpanded")


single_LEV2_3_df_QC <- mutate(single_LEV2_3_df_QC, identifiers_with_compoundName_instead_null = single_LEV2_3_df_QC$identifier)
single_LEV2_3_df_QC$identifiers_with_compoundName_instead_null[which(single_LEV2_3_df_QC$identifiers_with_compoundName_instead_null == "null")] <- single_LEV2_3_df_QC$compoundName[which(single_LEV2_3_df_QC$identifiers_with_compoundName_instead_null == "null")]


single_LEV2_3_df_QC <- single_LEV2_3_df_QC %>%
  add_column(type_chromatography = NA,
             type_processing_tool = NA,
             type_chemical_database = NA,
             type_annotation = NA,
             .after = "type_analysis") %>%
  mutate(type_chromatography = factor(ifelse(grepl("HILICPOS", type_analysis), "HILICPOS",
                                             ifelse(grepl("RPLCNEG", type_analysis), "RPLCNEG", NA)), levels = c("RPLCNEG", "HILICPOS")),
         type_processing_tool = factor(ifelse(grepl("patRoon", type_analysis), "patRoon",
                                              ifelse(grepl("MSDial", type_analysis), "MSDial", NA)), levels = c("patRoon", "MSDial")),
         type_chemical_database = factor(ifelse(grepl("PubChemLite", type_analysis), "PubChemLite",
                                              ifelse(grepl("PublicMSP", type_analysis), "PublicMSP",
                                                     ifelse(grepl("WormJamExpanded", type_analysis), "WormJamExpanded", NA))), levels = c("PubChemLite", "PublicMSP", "WormJamExpanded")),
         type_annotation = factor(ifelse(type_processing_tool == "patRoon" & type_chemical_database == "PubChemLite", "patRoon_PubChemLite",
                                         ifelse(type_processing_tool == "patRoon" & type_chemical_database == "WormJamExpanded", "patRoon_WormJamExpanded",
                                                ifelse(type_processing_tool == "MSDial" & type_chemical_database == "PublicMSP", "MSDial_PublicMSP",
                                                       ifelse(type_processing_tool == "MSDial" & type_chemical_database == "WormJamExpanded", "MSDial_WormJamExpanded", NA)))), levels = c("patRoon_PubChemLite", "patRoon_WormJamExpanded", "MSDial_PublicMSP", "MSDial_WormJamExpanded")))


# cleaning names!!
single_LEV2_3_df_QC <- mutate(single_LEV2_3_df_QC,
                              Name_to_use = as.character(NA))


for (i in 1:length(pull(single_LEV2_3_df_QC, 1))) {
  
  this_compoundName <- single_LEV2_3_df_QC$compoundName[i]
  this_title <- single_LEV2_3_df_QC$Title[i]
  
  if (is.na(this_title)) {
    this_name_to_use <- this_compoundName
  } else if (is.na(this_compoundName)) {
    this_name_to_use <- this_title
  } else if (grepl("CID", this_title)) {
    this_name_to_use <- this_compoundName
  } else if (nchar(this_title) < nchar(this_compoundName)) {
    this_name_to_use <- this_title
  } else {
    this_name_to_use <- this_compoundName
  }
  
  this_name_to_use <- sub("^L-", "", this_name_to_use)
  this_name_to_use <- sub("^D-", "", this_name_to_use)
  this_name_to_use <- sub("^DL-", "", this_name_to_use)
  
  
  if (nchar(gsub("[^A-Z]", "", this_name_to_use)) > 4) {
    this_name_to_use <- sub("^(\\w)(.*)", "\\1\\L\\2", this_name_to_use, perl = TRUE)
  }
  
  
  single_LEV2_3_df_QC$Name_to_use[i] <- this_name_to_use
}



single_LEV2_3_df_QC_summ <- single_LEV2_3_df_QC %>%
  group_by(type_annotation) %>%
  summarize(N=n())




# actually, I can add the type of analyses

single_schemsep_LEV2_3_df_QC_empty <- single_LEV2_3_df_QC[0,]
single_schemsep_LEV2_3_df_QC_empty <- mutate(single_schemsep_LEV2_3_df_QC_empty,
                                             type_scheme = factor(character(), levels = c("scheme1", "scheme2")))
single_schemsep_LEV2_3_df_QC <- single_schemsep_LEV2_3_df_QC_empty

for (i in 1:length(pull(single_LEV2_3_df_QC, 1))) {
  
  this_new_rows <- single_schemsep_LEV2_3_df_QC_empty
  
  if (single_LEV2_3_df_QC$present_in_scheme1[i]) {
    
    this_new_rows_a <- single_schemsep_LEV2_3_df_QC_empty
    this_new_rows_a[1, colnames(single_LEV2_3_df_QC)] <- single_LEV2_3_df_QC[i,]
    this_new_rows_a[1,"type_scheme"] <- factor("scheme1", levels = c("scheme1", "scheme2"))
    
    if (grepl("HILICPOS_", this_new_rows_a$type_analysis)) {
      this_new_rows_a$type_analysis <- gsub("HILICPOS_", "HILICPOS_scheme1_", this_new_rows_a$type_analysis)
    }
    if (grepl("RPLCNEG_", this_new_rows_a$type_analysis)) {
      this_new_rows_a$type_analysis <- gsub("RPLCNEG_", "RPLCNEG_scheme1_", this_new_rows_a$type_analysis)
    }
    
    this_new_rows <- bind_rows(this_new_rows, this_new_rows_a)
  }
  
  if (single_LEV2_3_df_QC$present_in_scheme2[i]) {
    
    this_new_rows_b <- single_schemsep_LEV2_3_df_QC_empty
    this_new_rows_b[1, colnames(single_LEV2_3_df_QC)] <- single_LEV2_3_df_QC[i,]
    this_new_rows_b[1,"type_scheme"] <- factor("scheme2", levels = c("scheme1", "scheme2"))
    
    if (grepl("HILICPOS_", this_new_rows_b$type_analysis)) {
      this_new_rows_b$type_analysis <- gsub("HILICPOS_", "HILICPOS_scheme2_", this_new_rows_b$type_analysis)
    }
    if (grepl("RPLCNEG_", this_new_rows_b$type_analysis)) {
      this_new_rows_b$type_analysis <- gsub("RPLCNEG_", "RPLCNEG_scheme2_", this_new_rows_b$type_analysis)
    }
    
    this_new_rows <- bind_rows(this_new_rows, this_new_rows_b)
  }
  
  single_schemsep_LEV2_3_df_QC <- bind_rows(single_schemsep_LEV2_3_df_QC,
                                            this_new_rows)
  
}

single_schemsep_LEV2_3_df_QC <- select(single_schemsep_LEV2_3_df_QC,
                                       feature, rt, mz, type_analysis, type_chromatography, type_processing_tool, type_chemical_database, type_scheme, AnnotationLevel,
                                       compoundName, identifier, neutral_formula, InChIKey, SMILES, Title,
                                       kingdom, superclass, class, subclass, level5, level6, level7, identifiers_with_compoundName_instead_null)



# actually, similar to present_in_scheme1 and present_in_scheme2, could do something similar for everything, even though loosing feature information


single_presencesLEV2_3_df_QC <- tibble(unique_identifiers_with_compoundName_instead_null = unique(single_LEV2_3_df_QC$identifiers_with_compoundName_instead_null),
                                       present_in_HILICPOS = as.logical(NA), present_in_RPLCNEG = as.logical(NA), present_in_scheme1 = as.logical(NA), present_in_scheme2 = as.logical(NA),
                                       present_in_MSDial = as.logical(NA), present_in_patRoon = as.logical(NA), present_in_PubChemLite = as.logical(NA), present_in_PublicMSP = as.logical(NA), present_in_WormJamExpanded = as.logical(NA),
                                       identifier = as.character(NA), compoundName = as.character(NA), neutral_formula = as.character(NA), InChIKey = as.character(NA),
                                       SMILES = as.character(NA), Title = as.character(NA), kingdom = as.character(NA), superclass = as.character(NA), class = as.character(NA),
                                       subclass = as.character(NA), level5 = as.character(NA), level6 = as.character(NA), level7 = as.character(NA))

compound_specification_columns <- c("identifier", "compoundName", "neutral_formula", "InChIKey", "SMILES", "Title", "kingdom", "superclass", "class", "subclass", "level5", "level6", "level7")
some_compound_specification_columns <- c("InChIKey", "Title", "kingdom", "superclass", "class", "subclass", "level5", "level6", "level7")

for (i in 1:length(single_presencesLEV2_3_df_QC$unique_identifiers_with_compoundName_instead_null)) {
  
  this_unique_identifiers_with_compoundName_instead_null <- single_presencesLEV2_3_df_QC$unique_identifiers_with_compoundName_instead_null[i]
  
  single_LEV2_3_df_QC_fill <- filter(single_LEV2_3_df_QC, identifiers_with_compoundName_instead_null == this_unique_identifiers_with_compoundName_instead_null)
  
  
  if (!all(map_lgl(single_LEV2_3_df_QC_fill[, some_compound_specification_columns],
                   ~ (all(.x == .x[1])) | any(is.na(.x))))) {warning(paste0("compound specifications different for identifier ", this_identifier))}
  
  
  single_presencesLEV2_3_df_QC[i, compound_specification_columns] <- single_LEV2_3_df_QC_fill[1, compound_specification_columns]
  
  if (any(grepl("HILICPOS", single_LEV2_3_df_QC_fill$type_analysis))) { single_presencesLEV2_3_df_QC$present_in_HILICPOS[i] <- TRUE } else { single_presencesLEV2_3_df_QC$present_in_HILICPOS[i] <- FALSE }
  if (any(grepl("RPLCNEG", single_LEV2_3_df_QC_fill$type_analysis))) { single_presencesLEV2_3_df_QC$present_in_RPLCNEG[i] <- TRUE } else { single_presencesLEV2_3_df_QC$present_in_RPLCNEG[i] <- FALSE }
  if (any(single_LEV2_3_df_QC_fill$present_in_scheme1)) { single_presencesLEV2_3_df_QC$present_in_scheme1[i] <- TRUE } else { single_presencesLEV2_3_df_QC$present_in_scheme1[i] <- FALSE }
  if (any(single_LEV2_3_df_QC_fill$present_in_scheme2)) { single_presencesLEV2_3_df_QC$present_in_scheme2[i] <- TRUE } else { single_presencesLEV2_3_df_QC$present_in_scheme2[i] <- FALSE }
  if (any(grepl("patRoon", single_LEV2_3_df_QC_fill$type_analysis))) { single_presencesLEV2_3_df_QC$present_in_patRoon[i] <- TRUE } else { single_presencesLEV2_3_df_QC$present_in_patRoon[i] <- FALSE }
  if (any(grepl("MSDial", single_LEV2_3_df_QC_fill$type_analysis))) { single_presencesLEV2_3_df_QC$present_in_MSDial[i] <- TRUE } else { single_presencesLEV2_3_df_QC$present_in_MSDial[i] <- FALSE }
  if (any(grepl("PubChemLite", single_LEV2_3_df_QC_fill$type_analysis))) { single_presencesLEV2_3_df_QC$present_in_PubChemLite[i] <- TRUE } else { single_presencesLEV2_3_df_QC$present_in_PubChemLite[i] <- FALSE }
  if (any(grepl("PublicMSP", single_LEV2_3_df_QC_fill$type_analysis))) { single_presencesLEV2_3_df_QC$present_in_PublicMSP[i] <- TRUE } else { single_presencesLEV2_3_df_QC$present_in_PublicMSP[i] <- FALSE }
  if (any(grepl("WormJamExpanded", single_LEV2_3_df_QC_fill$type_analysis))) { single_presencesLEV2_3_df_QC$present_in_WormJamExpanded[i] <- TRUE } else { single_presencesLEV2_3_df_QC$present_in_WormJamExpanded[i] <- FALSE }
  
  
}


# upset plot

# only annotation tools e chromatographic runs

single_presencesLEV2_3_df_QC_per_upset1 <- tibble(unique_identifiers_with_compoundName_instead_null = single_presencesLEV2_3_df_QC$unique_identifiers_with_compoundName_instead_null,
                                                  MSDial_PublicMSP = as.integer(NA),
                                                  MSDial_WormJamExpanded = as.integer(NA),
                                                  patRoon_PubChemLite = as.integer(NA),
                                                  patRoon_WormJamExpanded = as.integer(NA))


for (i in 1:length(pull(single_presencesLEV2_3_df_QC_per_upset1, 1))) {
  
  if (single_presencesLEV2_3_df_QC$present_in_MSDial[i] & single_presencesLEV2_3_df_QC$present_in_PublicMSP[i]) {
    single_presencesLEV2_3_df_QC_per_upset1$MSDial_PublicMSP[i] <- 1
  } else {
    single_presencesLEV2_3_df_QC_per_upset1$MSDial_PublicMSP[i] <- 0
  }
  
  if (single_presencesLEV2_3_df_QC$present_in_MSDial[i] & single_presencesLEV2_3_df_QC$present_in_WormJamExpanded[i]) {
    single_presencesLEV2_3_df_QC_per_upset1$MSDial_WormJamExpanded[i] <- 1
  } else {
    single_presencesLEV2_3_df_QC_per_upset1$MSDial_WormJamExpanded[i] <- 0
  }
  
  if (single_presencesLEV2_3_df_QC$present_in_patRoon [i] & single_presencesLEV2_3_df_QC$present_in_PubChemLite [i]) {
    single_presencesLEV2_3_df_QC_per_upset1$patRoon_PubChemLite[i] <- 1
  } else {
    single_presencesLEV2_3_df_QC_per_upset1$patRoon_PubChemLite[i] <- 0
  }
  
  if (single_presencesLEV2_3_df_QC$present_in_patRoon [i] & single_presencesLEV2_3_df_QC$present_in_WormJamExpanded[i]) {
    single_presencesLEV2_3_df_QC_per_upset1$patRoon_WormJamExpanded[i] <- 1
  } else {
    single_presencesLEV2_3_df_QC_per_upset1$patRoon_WormJamExpanded[i] <- 0
  }
  
  if (sum(single_presencesLEV2_3_df_QC_per_upset1$MSDial_PublicMSP[i],
          single_presencesLEV2_3_df_QC_per_upset1$MSDial_WormJamExpanded[i],
          single_presencesLEV2_3_df_QC_per_upset1$patRoon_PubChemLite[i],
          single_presencesLEV2_3_df_QC_per_upset1$patRoon_WormJamExpanded[i]) == 0) {stop(paste0(single_presencesLEV2_3_df_QC_per_upset1$unique_identifiers_with_compoundName_instead_null[i],
                                                                                                 " is not present in at least one annotation method"))}
  
}

single_presencesLEV2_3_df_QC_per_upset1 <- mutate_if(single_presencesLEV2_3_df_QC_per_upset1, is.double, as.integer)


#upset(single_presencesLEV2_3_df_QC_per_upset1)




# it throws an error that I don't know how to solve.. Error in xtfrm.data.frame(x) : cannot xtfrm data frames

# I'll try with the list:


single_schemsep_LEV2_3_df_QC_list_for_upset1 <- list(MSDial_PublicMSP = unique(filter(single_schemsep_LEV2_3_df_QC, type_processing_tool == "MSDial" & type_chemical_database == "PublicMSP")$identifiers_with_compoundName_instead_null),
                                                     MSDial_WormJamExpanded = unique(filter(single_schemsep_LEV2_3_df_QC, type_processing_tool == "MSDial" & type_chemical_database == "WormJamExpanded")$identifiers_with_compoundName_instead_null),
                                                     patRoon_PubChemLite = unique(filter(single_schemsep_LEV2_3_df_QC, type_processing_tool == "patRoon" & type_chemical_database == "PubChemLite")$identifiers_with_compoundName_instead_null),
                                                     patRoon_WormJamExpanded = unique(filter(single_schemsep_LEV2_3_df_QC, type_processing_tool == "patRoon" & type_chemical_database == "WormJamExpanded")$identifiers_with_compoundName_instead_null))

upset(fromList(single_schemsep_LEV2_3_df_QC_list_for_upset1))


png("upset1.png", width = 800, height = 600)
upset(fromList(single_schemsep_LEV2_3_df_QC_list_for_upset1))
dev.off() 



## with also the chromatographic:

single_schemsep_LEV2_3_df_QC_list_for_upset2 <- list(HILICPOS_MSDial_PublicMSP = unique(filter(single_schemsep_LEV2_3_df_QC, type_chromatography == "HILICPOS" & type_processing_tool == "MSDial" & type_chemical_database == "PublicMSP")$identifiers_with_compoundName_instead_null),
                                                     HILICPOS_MSDial_WormJamExpanded = unique(filter(single_schemsep_LEV2_3_df_QC, type_chromatography == "HILICPOS" & type_processing_tool == "MSDial" & type_chemical_database == "WormJamExpanded")$identifiers_with_compoundName_instead_null),
                                                     HILICPOS_patRoon_PubChemLite = unique(filter(single_schemsep_LEV2_3_df_QC, type_chromatography == "HILICPOS" & type_processing_tool == "patRoon" & type_chemical_database == "PubChemLite")$identifiers_with_compoundName_instead_null),
                                                     HILICPOS_patRoon_WormJamExpanded = unique(filter(single_schemsep_LEV2_3_df_QC, type_chromatography == "HILICPOS" & type_processing_tool == "patRoon" & type_chemical_database == "WormJamExpanded")$identifiers_with_compoundName_instead_null),
                                                     RPLCNEG_MSDial_PublicMSP = unique(filter(single_schemsep_LEV2_3_df_QC, type_chromatography == "RPLCNEG" & type_processing_tool == "MSDial" & type_chemical_database == "PublicMSP")$identifiers_with_compoundName_instead_null),
                                                     RPLCNEG_MSDial_WormJamExpanded = unique(filter(single_schemsep_LEV2_3_df_QC, type_chromatography == "RPLCNEG" & type_processing_tool == "MSDial" & type_chemical_database == "WormJamExpanded")$identifiers_with_compoundName_instead_null),
                                                     RPLCNEG_patRoon_PubChemLite = unique(filter(single_schemsep_LEV2_3_df_QC, type_chromatography == "RPLCNEG" & type_processing_tool == "patRoon" & type_chemical_database == "PubChemLite")$identifiers_with_compoundName_instead_null),
                                                     RPLCNEG_patRoon_WormJamExpanded = unique(filter(single_schemsep_LEV2_3_df_QC, type_chromatography == "RPLCNEG" & type_processing_tool == "patRoon" & type_chemical_database == "WormJamExpanded")$identifiers_with_compoundName_instead_null))

upset(fromList(single_schemsep_LEV2_3_df_QC_list_for_upset2))

png("upset2.png", width = 800, height = 600)
upset(fromList(single_schemsep_LEV2_3_df_QC_list_for_upset2))
dev.off() 



# Eulero-Venn:


new_euler_graph01 <- ggvenn(single_schemsep_LEV2_3_df_QC_list_for_upset1,
                            c("MSDial_PublicMSP", "MSDial_WormJamExpanded", "patRoon_PubChemLite", "patRoon_WormJamExpanded"),
                            show_elements = FALSE, label_sep = "\n", text_size = 6) + 
  theme(panel.background = element_rect(fill = "white"))
  
ggsave(filename = "new_euler_graph01.png", plot = new_euler_graph01,  width = 25, height = 25, units = "in")




###
######
############ -> I am doing everything again using the variable "Name_to_use", instead:
######
###


# add the type of analyses

single_schemsep_LEV2_3_df_QC_empty_bis <- single_LEV2_3_df_QC[0,]
single_schemsep_LEV2_3_df_QC_empty_bis <- mutate(single_schemsep_LEV2_3_df_QC_empty_bis,
                                             type_scheme = factor(character(), levels = c("scheme1", "scheme2")))
single_schemsep_LEV2_3_df_QC_bis <- single_schemsep_LEV2_3_df_QC_empty_bis

for (i in 1:length(pull(single_LEV2_3_df_QC, 1))) {
  
  this_new_rows <- single_schemsep_LEV2_3_df_QC_empty_bis
  
  if (single_LEV2_3_df_QC$present_in_scheme1[i]) {
    
    this_new_rows_a <- single_schemsep_LEV2_3_df_QC_empty_bis
    this_new_rows_a[1, colnames(single_LEV2_3_df_QC)] <- single_LEV2_3_df_QC[i,]
    this_new_rows_a[1,"type_scheme"] <- factor("scheme1", levels = c("scheme1", "scheme2"))
    
    if (grepl("HILICPOS_", this_new_rows_a$type_analysis)) {
      this_new_rows_a$type_analysis <- gsub("HILICPOS_", "HILICPOS_scheme1_", this_new_rows_a$type_analysis)
    }
    if (grepl("RPLCNEG_", this_new_rows_a$type_analysis)) {
      this_new_rows_a$type_analysis <- gsub("RPLCNEG_", "RPLCNEG_scheme1_", this_new_rows_a$type_analysis)
    }
    
    this_new_rows <- bind_rows(this_new_rows, this_new_rows_a)
    
  }
  
  if (single_LEV2_3_df_QC$present_in_scheme2[i]) {
    
    this_new_rows_b <- single_schemsep_LEV2_3_df_QC_empty_bis
    this_new_rows_b[1, colnames(single_LEV2_3_df_QC)] <- single_LEV2_3_df_QC[i,]
    this_new_rows_b[1,"type_scheme"] <- factor("scheme2", levels = c("scheme1", "scheme2"))
    
    if (grepl("HILICPOS_", this_new_rows_b$type_analysis)) {
      this_new_rows_b$type_analysis <- gsub("HILICPOS_", "HILICPOS_scheme2_", this_new_rows_b$type_analysis)
    }
    if (grepl("RPLCNEG_", this_new_rows_b$type_analysis)) {
      this_new_rows_b$type_analysis <- gsub("RPLCNEG_", "RPLCNEG_scheme2_", this_new_rows_b$type_analysis)
    }
    
    this_new_rows <- bind_rows(this_new_rows, this_new_rows_b)
  }
  
  single_schemsep_LEV2_3_df_QC_bis <- bind_rows(single_schemsep_LEV2_3_df_QC_bis,
                                            this_new_rows)
  
}

single_schemsep_LEV2_3_df_QC_bis <- select(single_schemsep_LEV2_3_df_QC_bis,
                                       feature, rt, mz, type_analysis, type_chromatography, type_processing_tool, type_chemical_database, type_scheme, AnnotationLevel,
                                       Name_to_use, compoundName, identifier, neutral_formula, InChIKey, SMILES, Title,
                                       kingdom, superclass, class, subclass, level5, level6, level7, identifiers_with_compoundName_instead_null)



# actually, similar to present_in_scheme1 and present_in_scheme2, could do something similar for everything, even though loosing feature information


single_presencesLEV2_3_df_QC_bis <- tibble(Name_to_use = unique(single_LEV2_3_df_QC$Name_to_use),
                                       present_in_HILICPOS = as.logical(NA), present_in_RPLCNEG = as.logical(NA), present_in_scheme1 = as.logical(NA), present_in_scheme2 = as.logical(NA),
                                       present_in_MSDial = as.logical(NA), present_in_patRoon = as.logical(NA), present_in_PubChemLite = as.logical(NA), present_in_PublicMSP = as.logical(NA), present_in_WormJamExpanded = as.logical(NA),
                                       identifier = as.character(NA), compoundName = as.character(NA), neutral_formula = as.character(NA), InChIKey = as.character(NA),
                                       SMILES = as.character(NA), Title = as.character(NA), identifiers_with_compoundName_instead_null = as.character(NA),
                                       kingdom = as.character(NA), superclass = as.character(NA), class = as.character(NA),
                                       subclass = as.character(NA), level5 = as.character(NA), level6 = as.character(NA), level7 = as.character(NA))

compound_specification_columns_bis <- c("identifier", "compoundName", "neutral_formula", "InChIKey", "SMILES", "Title", "identifiers_with_compoundName_instead_null", "kingdom", "superclass", "class", "subclass", "level5", "level6", "level7")



for (i in 1:length(single_presencesLEV2_3_df_QC_bis$Name_to_use)) {
  
  this_Name_to_use <- single_presencesLEV2_3_df_QC_bis$Name_to_use[i]
  
  single_LEV2_3_df_QC_fill <- filter(single_LEV2_3_df_QC, Name_to_use == this_Name_to_use)
  
  
  if (!all(map_lgl(single_LEV2_3_df_QC_fill[, "neutral_formula"],
                   ~ (all(.x == .x[1])) | any(is.na(.x))))) {warning(paste0("compound specifications different for ", this_Name_to_use))}
  
  
  single_presencesLEV2_3_df_QC_bis[i, compound_specification_columns] <- single_LEV2_3_df_QC_fill[1, compound_specification_columns]
  
  if (any(grepl("HILICPOS", single_LEV2_3_df_QC_fill$type_analysis))) { single_presencesLEV2_3_df_QC_bis$present_in_HILICPOS[i] <- TRUE } else { single_presencesLEV2_3_df_QC_bis$present_in_HILICPOS[i] <- FALSE }
  if (any(grepl("RPLCNEG", single_LEV2_3_df_QC_fill$type_analysis))) { single_presencesLEV2_3_df_QC_bis$present_in_RPLCNEG[i] <- TRUE } else { single_presencesLEV2_3_df_QC_bis$present_in_RPLCNEG[i] <- FALSE }
  if (any(single_LEV2_3_df_QC_fill$present_in_scheme1)) { single_presencesLEV2_3_df_QC_bis$present_in_scheme1[i] <- TRUE } else { single_presencesLEV2_3_df_QC_bis$present_in_scheme1[i] <- FALSE }
  if (any(single_LEV2_3_df_QC_fill$present_in_scheme2)) { single_presencesLEV2_3_df_QC_bis$present_in_scheme2[i] <- TRUE } else { single_presencesLEV2_3_df_QC_bis$present_in_scheme2[i] <- FALSE }
  if (any(grepl("patRoon", single_LEV2_3_df_QC_fill$type_analysis))) { single_presencesLEV2_3_df_QC_bis$present_in_patRoon[i] <- TRUE } else { single_presencesLEV2_3_df_QC_bis$present_in_patRoon[i] <- FALSE }
  if (any(grepl("MSDial", single_LEV2_3_df_QC_fill$type_analysis))) { single_presencesLEV2_3_df_QC_bis$present_in_MSDial[i] <- TRUE } else { single_presencesLEV2_3_df_QC_bis$present_in_MSDial[i] <- FALSE }
  if (any(grepl("PubChemLite", single_LEV2_3_df_QC_fill$type_analysis))) { single_presencesLEV2_3_df_QC_bis$present_in_PubChemLite[i] <- TRUE } else { single_presencesLEV2_3_df_QC_bis$present_in_PubChemLite[i] <- FALSE }
  if (any(grepl("PublicMSP", single_LEV2_3_df_QC_fill$type_analysis))) { single_presencesLEV2_3_df_QC_bis$present_in_PublicMSP[i] <- TRUE } else { single_presencesLEV2_3_df_QC_bis$present_in_PublicMSP[i] <- FALSE }
  if (any(grepl("WormJamExpanded", single_LEV2_3_df_QC_fill$type_analysis))) { single_presencesLEV2_3_df_QC_bis$present_in_WormJamExpanded[i] <- TRUE } else { single_presencesLEV2_3_df_QC_bis$present_in_WormJamExpanded[i] <- FALSE }
  
}


write_tsv(single_presencesLEV2_3_df_QC_bis, "single_presencesLEV2_3_df_QC_bis.txt")



# upset plot

# only annotation tools e chromatographic runs

single_presencesLEV2_3_df_QC_per_upsetb1 <- tibble(Name_to_use = single_presencesLEV2_3_df_QC_bis$Name_to_use,
                                                   MSDial_PublicMSP = as.integer(NA),
                                                   MSDial_WormJamExpanded = as.integer(NA),
                                                   patRoon_PubChemLite = as.integer(NA),
                                                   patRoon_WormJamExpanded = as.integer(NA))


for (i in 1:length(pull(single_presencesLEV2_3_df_QC_per_upsetb1, 1))) {
  
  if (single_presencesLEV2_3_df_QC_bis$present_in_MSDial[i] & single_presencesLEV2_3_df_QC_bis$present_in_PublicMSP[i]) {
    single_presencesLEV2_3_df_QC_per_upsetb1$MSDial_PublicMSP[i] <- 1
  } else {
    single_presencesLEV2_3_df_QC_per_upsetb1$MSDial_PublicMSP[i] <- 0
  }
  
  if (single_presencesLEV2_3_df_QC_bis$present_in_MSDial[i] & single_presencesLEV2_3_df_QC_bis$present_in_WormJamExpanded[i]) {
    single_presencesLEV2_3_df_QC_per_upsetb1$MSDial_WormJamExpanded[i] <- 1
  } else {
    single_presencesLEV2_3_df_QC_per_upsetb1$MSDial_WormJamExpanded[i] <- 0
  }
  
  if (single_presencesLEV2_3_df_QC_bis$present_in_patRoon [i] & single_presencesLEV2_3_df_QC_bis$present_in_PubChemLite [i]) {
    single_presencesLEV2_3_df_QC_per_upsetb1$patRoon_PubChemLite[i] <- 1
  } else {
    single_presencesLEV2_3_df_QC_per_upsetb1$patRoon_PubChemLite[i] <- 0
  }
  
  if (single_presencesLEV2_3_df_QC_bis$present_in_patRoon [i] & single_presencesLEV2_3_df_QC_bis$present_in_WormJamExpanded[i]) {
    single_presencesLEV2_3_df_QC_per_upsetb1$patRoon_WormJamExpanded[i] <- 1
  } else {
    single_presencesLEV2_3_df_QC_per_upsetb1$patRoon_WormJamExpanded[i] <- 0
  }
  
  if (sum(single_presencesLEV2_3_df_QC_per_upsetb1$MSDial_PublicMSP[i],
          single_presencesLEV2_3_df_QC_per_upsetb1$MSDial_WormJamExpanded[i],
          single_presencesLEV2_3_df_QC_per_upsetb1$patRoon_PubChemLite[i],
          single_presencesLEV2_3_df_QC_per_upsetb1$patRoon_WormJamExpanded[i]) == 0) {stop(paste0(single_presencesLEV2_3_df_QC_per_upsetb1$unique_identifiers_with_compoundName_instead_null[i],
                                                                                                 " is not present in at least one annotation method"))}
  
}

single_presencesLEV2_3_df_QC_per_upsetb1 <- mutate_if(single_presencesLEV2_3_df_QC_per_upsetb1, is.double, as.integer)
glimpse(single_presencesLEV2_3_df_QC_per_upsetb1)

#upset(single_presencesLEV2_3_df_QC_per_upsetb1)




# it throws an error that I don't know how to solve.. Error in xtfrm.data.frame(x) : cannot xtfrm data frames

# I'll try with the list:


single_schemsep_LEV2_3_df_QC_list_for_upsetb1 <- list(MSDial_PublicMSP = filter(single_presencesLEV2_3_df_QC_bis, present_in_MSDial == TRUE & present_in_PublicMSP == TRUE)$Name_to_use,
                                                      MSDial_WormJamExpanded = filter(single_presencesLEV2_3_df_QC_bis, present_in_MSDial == TRUE & present_in_WormJamExpanded == TRUE)$Name_to_use,
                                                      patRoon_PubChemLite = filter(single_presencesLEV2_3_df_QC_bis, present_in_patRoon == TRUE & present_in_PubChemLite == TRUE)$Name_to_use,
                                                      patRoon_WormJamExpanded = filter(single_presencesLEV2_3_df_QC_bis, present_in_patRoon == TRUE & present_in_WormJamExpanded == TRUE)$Name_to_use)


map_dbl(single_schemsep_LEV2_3_df_QC_list_for_upsetb1, length)


png("upsetb1.png", width = 800, height = 600)
upset(fromList(single_schemsep_LEV2_3_df_QC_list_for_upsetb1), text.scale = 1.8)
dev.off()






# Eulero-Venn:


new_euler_graph_b01 <- ggvenn(single_schemsep_LEV2_3_df_QC_list_for_upsetb1,
                            c("MSDial_PublicMSP", "MSDial_WormJamExpanded", "patRoon_PubChemLite", "patRoon_WormJamExpanded"),
                            show_elements = FALSE, show_percentage = TRUE, label_sep = "\n", text_size = 6) + 
  theme(panel.background = element_rect(fill = "white"))
ggsave(filename = "new_euler_graph_b01.png", plot = new_euler_graph_b01,  width = 25, height = 25, units = "in")



new_euler_graph_b02 <- ggvenn(single_schemsep_LEV2_3_df_QC_list_for_upsetb1,
                              c("MSDial_PublicMSP", "MSDial_WormJamExpanded", "patRoon_PubChemLite", "patRoon_WormJamExpanded"),
                              show_elements = FALSE, show_percentage = FALSE, label_sep = "\n", text_size = 6) + 
  theme(panel.background = element_rect(fill = "white"))
ggsave(filename = "new_euler_graph_b02.png", plot = new_euler_graph_b02,  width = 25, height = 25, units = "in")


### with shorter names:

single_schemsep_LEV2_3_df_QC_list_for_upsetb1_shrt_nms <- single_schemsep_LEV2_3_df_QC_list_for_upsetb1
names(single_schemsep_LEV2_3_df_QC_list_for_upsetb1_shrt_nms)[which(names(single_schemsep_LEV2_3_df_QC_list_for_upsetb1_shrt_nms)=="MSDial_PublicMSP")] <- "MS-DIAL\nPublicMSP"
names(single_schemsep_LEV2_3_df_QC_list_for_upsetb1_shrt_nms)[which(names(single_schemsep_LEV2_3_df_QC_list_for_upsetb1_shrt_nms)=="MSDial_WormJamExpanded")] <- "MS-DIAL\nWormJamExpanded"
names(single_schemsep_LEV2_3_df_QC_list_for_upsetb1_shrt_nms)[which(names(single_schemsep_LEV2_3_df_QC_list_for_upsetb1_shrt_nms)=="patRoon_PubChemLite")] <- "patRoon\nPubChemLite"
names(single_schemsep_LEV2_3_df_QC_list_for_upsetb1_shrt_nms)[which(names(single_schemsep_LEV2_3_df_QC_list_for_upsetb1_shrt_nms)=="patRoon_WormJamExpanded")] <- "patRoon\nWormJamExpanded"

new_euler_graph_b02_shrt_nms <- ggvenn(single_schemsep_LEV2_3_df_QC_list_for_upsetb1_shrt_nms,
                              c("MS-DIAL\nPublicMSP", "MS-DIAL\nWormJamExpanded", "patRoon\nPubChemLite", "patRoon\nWormJamExpanded"),
                              show_elements = FALSE, show_percentage = FALSE, label_sep = "\n", text_size = 25, set_name_size = 9) + 
  theme(panel.background = element_rect(fill = "white"))
ggsave(filename = "new_euler_graph_b02_shrt_nms.png", plot = new_euler_graph_b02_shrt_nms,  width = 25, height = 20, units = "in")



## to add compounds detail to that picture:


single_pres2_3_QC_A_only_MSDIAL_publicMSP <- filter(single_presencesLEV2_3_df_QC_bis, present_in_MSDial==TRUE & present_in_PublicMSP == TRUE & present_in_patRoon ==FALSE & present_in_PubChemLite == FALSE & present_in_WormJamExpanded ==FALSE)
write_tsv(single_pres2_3_QC_A_only_MSDIAL_publicMSP, "single_pres2_3_QC_A_only_MSDIAL_publicMSP.txt")

single_pres2_3_QC_B_MSDIAL_publicMSP_WJ <- filter(single_presencesLEV2_3_df_QC_bis, present_in_MSDial==TRUE & present_in_PublicMSP == TRUE & present_in_patRoon ==FALSE & present_in_PubChemLite == FALSE & present_in_WormJamExpanded ==TRUE)
write_tsv(single_pres2_3_QC_B_MSDIAL_publicMSP_WJ, "single_pres2_3_QC_B_MSDIAL_publicMSP_WJ.txt")

single_pres2_3_QC_C_only_MSDIAL_WJ <- filter(single_presencesLEV2_3_df_QC_bis, present_in_MSDial==TRUE & present_in_PublicMSP == FALSE & present_in_patRoon ==FALSE & present_in_PubChemLite == FALSE & present_in_WormJamExpanded ==TRUE)
write_tsv(single_pres2_3_QC_C_only_MSDIAL_WJ, "single_pres2_3_QC_C_only_MSDIAL_WJ.txt")

single_pres2_3_QC_D_only_patRoonPCL <- filter(single_presencesLEV2_3_df_QC_bis, present_in_MSDial==FALSE & present_in_PublicMSP == FALSE & present_in_patRoon ==TRUE & present_in_PubChemLite == TRUE & present_in_WormJamExpanded ==FALSE)
write_tsv(single_pres2_3_QC_D_only_patRoonPCL, "single_pres2_3_QC_D_only_patRoonPCL.txt")

single_pres2_3_QC_E_patRoon_PCL_WJ <- filter(single_presencesLEV2_3_df_QC_bis, present_in_MSDial==FALSE & present_in_PublicMSP == FALSE & present_in_patRoon ==TRUE & present_in_PubChemLite == TRUE & present_in_WormJamExpanded ==TRUE)
write_tsv(single_pres2_3_QC_E_patRoon_PCL_WJ, "single_pres2_3_QC_E_patRoon_PCL_WJL.txt")

single_pres2_3_QC_F_only_patRoon_WJ <- filter(single_presencesLEV2_3_df_QC_bis, present_in_MSDial==FALSE & present_in_PublicMSP == FALSE & present_in_patRoon ==TRUE & present_in_PubChemLite == FALSE & present_in_WormJamExpanded ==TRUE)
write_tsv(single_pres2_3_QC_F_only_patRoon_WJ, "single_pres2_3_QC_F_only_patRoon_WJ.txt")

single_pres2_3_QC_G_MSDIALpublic_patRoonPCL <- filter(single_presencesLEV2_3_df_QC_bis, present_in_MSDial==TRUE & present_in_PublicMSP == TRUE & present_in_patRoon ==TRUE & present_in_PubChemLite == TRUE & present_in_WormJamExpanded ==FALSE)
write_tsv(single_pres2_3_QC_G_MSDIALpublic_patRoonPCL, "single_pres2_3_QC_G_MSDIALpublic_patRoonPCL.txt")

single_pres2_3_QC_H_all <- filter(single_presencesLEV2_3_df_QC_bis, present_in_MSDial==TRUE & present_in_PublicMSP == TRUE & present_in_patRoon ==TRUE & present_in_PubChemLite == TRUE & present_in_WormJamExpanded ==TRUE)
write_tsv(single_pres2_3_QC_H_all, "single_pres2_3_QC_H_all.txt")

single_pres2_3_QC_I_MSDIALpublicWJ_patRoon_WJ <- filter(single_presencesLEV2_3_df_QC_bis, present_in_MSDial==TRUE & present_in_PublicMSP == TRUE & present_in_patRoon ==TRUE & present_in_PubChemLite == FALSE & present_in_WormJamExpanded ==TRUE)
write_tsv(single_pres2_3_QC_I_MSDIALpublicWJ_patRoon_WJ, "single_pres2_3_QC_I_MSDIALpublicWJ_patRoon_WJ.txt")



### biologic significance:
# grouping:

single_presencesLEV2_3_df_QC_bis_sum_kingdom <- single_presencesLEV2_3_df_QC_bis %>% 
  group_by(kingdom) %>%
  summarize(N = n()) %>%
  mutate(Perc = round(N/sum(N)*100, digits = 1))

single_presencesLEV2_3_df_QC_bis_sum_superclass <- single_presencesLEV2_3_df_QC_bis %>% 
  group_by(superclass) %>%
  summarize(N = n()) %>%
  mutate(Perc = round(N/sum(N)*100, digits = 1))

single_presencesLEV2_3_df_QC_bis_sum_class <- single_presencesLEV2_3_df_QC_bis %>% 
  group_by(class) %>%
  summarize(N = n()) %>%
  mutate(Perc = round(N/sum(N)*100, digits = 1))

single_presencesLEV2_3_df_QC_bis_sum_subclass <- single_presencesLEV2_3_df_QC_bis %>% 
  group_by(subclass) %>%
  summarize(N = n()) %>%
  mutate(Perc = round(N/sum(N)*100, digits = 1))

single_presencesLEV2_3_df_QC_bis_sum_level5 <- single_presencesLEV2_3_df_QC_bis %>% 
  group_by(level5) %>%
  summarize(N = n()) %>%
  mutate(Perc = round(N/sum(N)*100, digits = 1))

single_presencesLEV2_3_df_QC_bis_sum_level6 <- single_presencesLEV2_3_df_QC_bis %>% 
  group_by(level6) %>%
  summarize(N = n()) %>%
  mutate(Perc = round(N/sum(N)*100, digits = 1))

single_presencesLEV2_3_df_QC_bis_sum_level7 <- single_presencesLEV2_3_df_QC_bis %>% 
  group_by(level7) %>%
  summarize(N = n()) %>%
  mutate(Perc = round(N/sum(N)*100, digits = 1))




######## 
## Sankey Diagram with categories of molecules:

# it's 0-indexed!
# that's why I add this_source - 1  and  this_target - 1

single_presencesLEV2_3_df_QC_bis_noNAkd <- filter(single_presencesLEV2_3_df_QC_bis, !is.na(kingdom))

all_names_QC <- c(unique(single_presencesLEV2_3_df_QC_bis_noNAkd$kingdom),
                  unique(single_presencesLEV2_3_df_QC_bis_noNAkd$superclass),
                  unique(single_presencesLEV2_3_df_QC_bis_noNAkd$class),
                  unique(single_presencesLEV2_3_df_QC_bis_noNAkd$subclass),
                  unique(single_presencesLEV2_3_df_QC_bis_noNAkd$level5),
                  unique(single_presencesLEV2_3_df_QC_bis_noNAkd$level6),
                  unique(single_presencesLEV2_3_df_QC_bis_noNAkd$level7))


all_names_QC <- all_names_QC[!is.na(all_names_QC)]

##a check: this MUST be FALSE (no duplicated in the names of categories!)
any(duplicated(all_names_QC))

all_links_QC <- data.frame(source = integer(length = 0),
                        target = integer(length = 0),
                        value = double(length = 0))


for (i in 1:length(pull(single_presencesLEV2_3_df_QC_bis_noNAkd, 1))) {
  this_source <- which(all_names_QC == single_presencesLEV2_3_df_QC_bis_noNAkd$kingdom[i])
  
  this_target <- which(all_names_QC == single_presencesLEV2_3_df_QC_bis_noNAkd$superclass[i])
  
  this_value <- 1
  
  if (length(this_source) == 1 & length(this_target) == 1) {
    
    this_source <- this_source - 1
    
    this_target <- this_target - 1
    
    this_df <- data.frame(source = as.integer(this_source),
                          target = as.integer(this_target),
                          value = as.double(this_value))
    
    all_links_QC <- rbind(all_links_QC, this_df)
    
  }
}


for (i in 1:length(pull(single_presencesLEV2_3_df_QC_bis_noNAkd, 1))) {
  this_source <- which(all_names_QC == single_presencesLEV2_3_df_QC_bis_noNAkd$superclass[i])
  
  this_target <-which(all_names_QC == single_presencesLEV2_3_df_QC_bis_noNAkd$class[i])
  
  this_value <- 1
  
  if (length(this_source) == 1 & length(this_target) == 1) {
    
    this_source <- this_source - 1
    
    this_target <- this_target - 1
    
    this_df <- data.frame(source = as.integer(this_source),
                          target = as.integer(this_target),
                          value = as.double(this_value))
    
    all_links_QC <- rbind(all_links_QC, this_df)
    
  }
}

for (i in 1:length(pull(single_presencesLEV2_3_df_QC_bis_noNAkd, 1))) {
  this_source <- which(all_names_QC == single_presencesLEV2_3_df_QC_bis_noNAkd$class[i])
  
  this_target <-which(all_names_QC == single_presencesLEV2_3_df_QC_bis_noNAkd$subclass[i])
  
  this_value <- 1
  
  if (length(this_source) == 1 & length(this_target) == 1) {
    
    this_source <- this_source - 1
    
    this_target <- this_target - 1
    
    this_df <- data.frame(source = as.integer(this_source),
                          target = as.integer(this_target),
                          value = as.double(this_value))
    
    all_links_QC <- rbind(all_links_QC, this_df)
    
  }
}


for (i in 1:length(pull(single_presencesLEV2_3_df_QC_bis_noNAkd, 1))) {
  this_source <- which(all_names_QC == single_presencesLEV2_3_df_QC_bis_noNAkd$subclass[i])
  
  this_target <-which(all_names_QC == single_presencesLEV2_3_df_QC_bis_noNAkd$level5[i])
  
  this_value <- 1
  
  if (length(this_source) == 1 & length(this_target) == 1) {
    
    this_source <- this_source - 1
    
    this_target <- this_target - 1
    
    this_df <- data.frame(source = as.integer(this_source),
                          target = as.integer(this_target),
                          value = as.double(this_value))
    
    all_links_QC <- rbind(all_links_QC, this_df)
    
  }
}


for (i in 1:length(pull(single_presencesLEV2_3_df_QC_bis_noNAkd, 1))) {
  this_source <- which(all_names_QC == single_presencesLEV2_3_df_QC_bis_noNAkd$level5[i])
  
  this_target <-which(all_names_QC == single_presencesLEV2_3_df_QC_bis_noNAkd$level6[i])
  
  this_value <- 1
  
  if (length(this_source) == 1 & length(this_target) == 1) {
    
    this_source <- this_source - 1
    
    this_target <- this_target - 1
    
    this_df <- data.frame(source = as.integer(this_source),
                          target = as.integer(this_target),
                          value = as.double(this_value))
    
    all_links_QC <- rbind(all_links_QC, this_df)
    
  }
}


for (i in 1:length(pull(single_presencesLEV2_3_df_QC_bis_noNAkd, 1))) {
  this_source <- which(all_names_QC == single_presencesLEV2_3_df_QC_bis_noNAkd$level6[i])
  
  this_target <-which(all_names_QC == single_presencesLEV2_3_df_QC_bis_noNAkd$level7[i])
  
  this_value <- 1
  
  if (length(this_source) == 1 & length(this_target) == 1) {
    
    this_source <- this_source - 1
    
    this_target <- this_target - 1
    
    this_df <- data.frame(source = as.integer(this_source),
                          target = as.integer(this_target),
                          value = as.double(this_value))
    
    all_links_QC <- rbind(all_links_QC, this_df)
    
  }
}





all_links_QC <- mutate(all_links_QC, from_to = paste("from", source, "to", target))

Unique_Duplicated_QC <- unique(all_links_QC$from_to[which(duplicated(all_links_QC$from_to))])

for (d in Unique_Duplicated_QC) {
  df_fil <- filter(all_links_QC, from_to == d)
  
  this_added_value <- length(df_fil$from_to)
  
  all_links_QC[which(all_links_QC$from_to==d)[1],"value"] <- this_added_value
}

all_links_QC <- filter(all_links_QC, !duplicated(from_to))




List_for_Sankey_Diagram_QC <- list(nodes = data.frame(name = all_names_QC),
                                links = all_links_QC)


My_sankeyNetwork_QC <- sankeyNetwork(Links = List_for_Sankey_Diagram_QC$links, Nodes = List_for_Sankey_Diagram_QC$nodes,
                                     Source = "source",
                                     Target = "target",
                                     Value = "value",
                                     NodeID = "name",
                                     fontSize = 12, nodeWidth = 30,
                                     sinksRight = FALSE)
saveNetwork(My_sankeyNetwork_QC, "My_sankeyNetwork_QC.html")
webshot::webshot("My_sankeyNetwork_QC.html","My_sankeyNetwork_QC.png", vwidth = 1600, vheight = 3900)




### similarly but without the several contained in MS-dial publicMSP

single_presencesLEV2_3_df_QC_bis_NO_MSDPMSP <- single_presencesLEV2_3_df_QC_bis[which(!(single_presencesLEV2_3_df_QC_bis$present_in_MSDial == TRUE &
                                                                                          single_presencesLEV2_3_df_QC_bis$present_in_PublicMSP == TRUE &
                                                                                          single_presencesLEV2_3_df_QC_bis$present_in_WormJamExpanded == FALSE &
                                                                                          single_presencesLEV2_3_df_QC_bis$present_in_patRoon == FALSE &
                                                                                          single_presencesLEV2_3_df_QC_bis$present_in_PubChemLite == FALSE)),]



### biologic significance:
# grouping:

single_presencesLEV2_3_df_QC_bis_NO_MSDPMSP_sum_kingdom <- single_presencesLEV2_3_df_QC_bis_NO_MSDPMSP %>% 
  group_by(kingdom) %>%
  summarize(N = n()) %>%
  mutate(Perc = round(N/sum(N)*100, digits = 1))

single_presencesLEV2_3_df_QC_bis_NO_MSDPMSP_sum_superclass <- single_presencesLEV2_3_df_QC_bis_NO_MSDPMSP %>% 
  group_by(superclass) %>%
  summarize(N = n()) %>%
  mutate(Perc = round(N/sum(N)*100, digits = 1))

single_presencesLEV2_3_df_QC_bis_NO_MSDPMSP_sum_class <- single_presencesLEV2_3_df_QC_bis_NO_MSDPMSP %>% 
  group_by(class) %>%
  summarize(N = n()) %>%
  mutate(Perc = round(N/sum(N)*100, digits = 1))

single_presencesLEV2_3_df_QC_bis_NO_MSDPMSP_sum_subclass <- single_presencesLEV2_3_df_QC_bis_NO_MSDPMSP %>% 
  group_by(subclass) %>%
  summarize(N = n()) %>%
  mutate(Perc = round(N/sum(N)*100, digits = 1))

single_presencesLEV2_3_df_QC_bis_NO_MSDPMSP_sum_level5 <- single_presencesLEV2_3_df_QC_bis_NO_MSDPMSP %>% 
  group_by(level5) %>%
  summarize(N = n()) %>%
  mutate(Perc = round(N/sum(N)*100, digits = 1))

single_presencesLEV2_3_df_QC_bis_NO_MSDPMSP_sum_level6 <- single_presencesLEV2_3_df_QC_bis_NO_MSDPMSP %>% 
  group_by(level6) %>%
  summarize(N = n()) %>%
  mutate(Perc = round(N/sum(N)*100, digits = 1))

single_presencesLEV2_3_df_QC_bis_NO_MSDPMSP_sum_level7 <- single_presencesLEV2_3_df_QC_bis_NO_MSDPMSP %>% 
  group_by(level7) %>%
  summarize(N = n()) %>%
  mutate(Perc = round(N/sum(N)*100, digits = 1))




######## 
## Sankey Diagram with categories of molecules:

# it's 0-indexed!
# that's why I add this_source - 1  and  this_target - 1

single_presencesLEV2_3_df_QC_bis_NO_MSDPMSP_noNAkd <- filter(single_presencesLEV2_3_df_QC_bis_NO_MSDPMSP, !is.na(kingdom))

all_names_QC_NO_MSDPMSP_QC <- c(unique(single_presencesLEV2_3_df_QC_bis_NO_MSDPMSP_noNAkd$kingdom),
                                unique(single_presencesLEV2_3_df_QC_bis_NO_MSDPMSP_noNAkd$superclass),
                                unique(single_presencesLEV2_3_df_QC_bis_NO_MSDPMSP_noNAkd$class),
                                unique(single_presencesLEV2_3_df_QC_bis_NO_MSDPMSP_noNAkd$subclass),
                                unique(single_presencesLEV2_3_df_QC_bis_NO_MSDPMSP_noNAkd$level5),
                                unique(single_presencesLEV2_3_df_QC_bis_NO_MSDPMSP_noNAkd$level6),
                                unique(single_presencesLEV2_3_df_QC_bis_NO_MSDPMSP_noNAkd$level7))


all_names_QC_NO_MSDPMSP_QC <- all_names_QC_NO_MSDPMSP_QC[!is.na(all_names_QC_NO_MSDPMSP_QC)]

##a check: this MUST be FALSE (no duplicated in the names of categories!)
any(duplicated(all_names_QC_NO_MSDPMSP_QC))

all_links_QC_NO_MSDPMSP_QC <- data.frame(source = integer(length = 0),
                                         target = integer(length = 0),
                                         value = double(length = 0))


for (i in 1:length(pull(single_presencesLEV2_3_df_QC_bis_NO_MSDPMSP_noNAkd, 1))) {
  this_source <- which(all_names_QC_NO_MSDPMSP_QC == single_presencesLEV2_3_df_QC_bis_NO_MSDPMSP_noNAkd$kingdom[i])
  
  this_target <- which(all_names_QC_NO_MSDPMSP_QC == single_presencesLEV2_3_df_QC_bis_NO_MSDPMSP_noNAkd$superclass[i])
  
  this_value <- 1
  
  if (length(this_source) == 1 & length(this_target) == 1) {
    
    this_source <- this_source - 1
    
    this_target <- this_target - 1
    
    this_df <- data.frame(source = as.integer(this_source),
                          target = as.integer(this_target),
                          value = as.double(this_value))
    
    all_links_QC_NO_MSDPMSP_QC <- rbind(all_links_QC_NO_MSDPMSP_QC, this_df)
    
  }
}


for (i in 1:length(pull(single_presencesLEV2_3_df_QC_bis_NO_MSDPMSP_noNAkd, 1))) {
  this_source <- which(all_names_QC_NO_MSDPMSP_QC == single_presencesLEV2_3_df_QC_bis_NO_MSDPMSP_noNAkd$superclass[i])
  
  this_target <-which(all_names_QC_NO_MSDPMSP_QC == single_presencesLEV2_3_df_QC_bis_NO_MSDPMSP_noNAkd$class[i])
  
  this_value <- 1
  
  if (length(this_source) == 1 & length(this_target) == 1) {
    
    this_source <- this_source - 1
    
    this_target <- this_target - 1
    
    this_df <- data.frame(source = as.integer(this_source),
                          target = as.integer(this_target),
                          value = as.double(this_value))
    
    all_links_QC_NO_MSDPMSP_QC <- rbind(all_links_QC_NO_MSDPMSP_QC, this_df)
    
  }
}

for (i in 1:length(pull(single_presencesLEV2_3_df_QC_bis_NO_MSDPMSP_noNAkd, 1))) {
  this_source <- which(all_names_QC_NO_MSDPMSP_QC == single_presencesLEV2_3_df_QC_bis_NO_MSDPMSP_noNAkd$class[i])
  
  this_target <-which(all_names_QC_NO_MSDPMSP_QC == single_presencesLEV2_3_df_QC_bis_NO_MSDPMSP_noNAkd$subclass[i])
  
  this_value <- 1
  
  if (length(this_source) == 1 & length(this_target) == 1) {
    
    this_source <- this_source - 1
    
    this_target <- this_target - 1
    
    this_df <- data.frame(source = as.integer(this_source),
                          target = as.integer(this_target),
                          value = as.double(this_value))
    
    all_links_QC_NO_MSDPMSP_QC <- rbind(all_links_QC_NO_MSDPMSP_QC, this_df)
    
  }
}


for (i in 1:length(pull(single_presencesLEV2_3_df_QC_bis_NO_MSDPMSP_noNAkd, 1))) {
  this_source <- which(all_names_QC_NO_MSDPMSP_QC == single_presencesLEV2_3_df_QC_bis_NO_MSDPMSP_noNAkd$subclass[i])
  
  this_target <-which(all_names_QC_NO_MSDPMSP_QC == single_presencesLEV2_3_df_QC_bis_NO_MSDPMSP_noNAkd$level5[i])
  
  this_value <- 1
  
  if (length(this_source) == 1 & length(this_target) == 1) {
    
    this_source <- this_source - 1
    
    this_target <- this_target - 1
    
    this_df <- data.frame(source = as.integer(this_source),
                          target = as.integer(this_target),
                          value = as.double(this_value))
    
    all_links_QC_NO_MSDPMSP_QC <- rbind(all_links_QC_NO_MSDPMSP_QC, this_df)
    
  }
}


for (i in 1:length(pull(single_presencesLEV2_3_df_QC_bis_NO_MSDPMSP_noNAkd, 1))) {
  this_source <- which(all_names_QC_NO_MSDPMSP_QC == single_presencesLEV2_3_df_QC_bis_NO_MSDPMSP_noNAkd$level5[i])
  
  this_target <-which(all_names_QC_NO_MSDPMSP_QC == single_presencesLEV2_3_df_QC_bis_NO_MSDPMSP_noNAkd$level6[i])
  
  this_value <- 1
  
  if (length(this_source) == 1 & length(this_target) == 1) {
    
    this_source <- this_source - 1
    
    this_target <- this_target - 1
    
    this_df <- data.frame(source = as.integer(this_source),
                          target = as.integer(this_target),
                          value = as.double(this_value))
    
    all_links_QC_NO_MSDPMSP_QC <- rbind(all_links_QC_NO_MSDPMSP_QC, this_df)
    
  }
}


for (i in 1:length(pull(single_presencesLEV2_3_df_QC_bis_NO_MSDPMSP_noNAkd, 1))) {
  this_source <- which(all_names_QC_NO_MSDPMSP_QC == single_presencesLEV2_3_df_QC_bis_NO_MSDPMSP_noNAkd$level6[i])
  
  this_target <-which(all_names_QC_NO_MSDPMSP_QC == single_presencesLEV2_3_df_QC_bis_NO_MSDPMSP_noNAkd$level7[i])
  
  this_value <- 1
  
  if (length(this_source) == 1 & length(this_target) == 1) {
    
    this_source <- this_source - 1
    
    this_target <- this_target - 1
    
    this_df <- data.frame(source = as.integer(this_source),
                          target = as.integer(this_target),
                          value = as.double(this_value))
    
    all_links_QC_NO_MSDPMSP_QC <- rbind(all_links_QC_NO_MSDPMSP_QC, this_df)
    
  }
}





all_links_QC_NO_MSDPMSP_QC <- mutate(all_links_QC_NO_MSDPMSP_QC, from_to = paste("from", source, "to", target))

Unique_Duplicated_QC_NO_MSDPMSP <- unique(all_links_QC_NO_MSDPMSP_QC$from_to[which(duplicated(all_links_QC_NO_MSDPMSP_QC$from_to))])

for (d in Unique_Duplicated_QC_NO_MSDPMSP) {
  df_fil <- filter(all_links_QC_NO_MSDPMSP_QC, from_to == d)
  
  this_added_value <- length(df_fil$from_to)
  
  all_links_QC_NO_MSDPMSP_QC[which(all_links_QC_NO_MSDPMSP_QC$from_to==d)[1],"value"] <- this_added_value
}

all_links_QC_NO_MSDPMSP_QC <- filter(all_links_QC_NO_MSDPMSP_QC, !duplicated(from_to))




List_for_Sankey_Diagram_QC_NO_MSDPMSP <- list(nodes = data.frame(name = all_names_QC_NO_MSDPMSP_QC),
                                              links = all_links_QC_NO_MSDPMSP_QC)


My_sankeyNetwork_QC_NO_MSDPMSP <- sankeyNetwork(Links = List_for_Sankey_Diagram_QC_NO_MSDPMSP$links, Nodes = List_for_Sankey_Diagram_QC_NO_MSDPMSP$nodes,
                                                Source = "source",
                                                Target = "target",
                                                Value = "value",
                                                NodeID = "name",
                                                fontSize = 12, nodeWidth = 30,
                                                sinksRight = FALSE)
saveNetwork(My_sankeyNetwork_QC_NO_MSDPMSP, "My_sankeyNetwork_QC_NO_MSDPMSP.html")
webshot::webshot("My_sankeyNetwork_QC_NO_MSDPMSP.html","My_sankeyNetwork_QC_NO_MSDPMSP.png", vwidth = 1900, vheight = 1400)



### similarly but without the several contained in MS-dial

single_presencesLEV2_3_df_QC_bis_NO_MSD <- single_presencesLEV2_3_df_QC_bis[which(single_presencesLEV2_3_df_QC_bis$present_in_patRoon == TRUE),]



### biologic significance:
# grouping:

single_presencesLEV2_3_df_QC_bis_NO_MSD_sum_kingdom <- single_presencesLEV2_3_df_QC_bis_NO_MSD %>% 
  group_by(kingdom) %>%
  summarize(N = n()) %>%
  mutate(Perc = round(N/sum(N)*100, digits = 1))

single_presencesLEV2_3_df_QC_bis_NO_MSD_sum_superclass <- single_presencesLEV2_3_df_QC_bis_NO_MSD %>% 
  group_by(superclass) %>%
  summarize(N = n()) %>%
  mutate(Perc = round(N/sum(N)*100, digits = 1))

single_presencesLEV2_3_df_QC_bis_NO_MSD_sum_class <- single_presencesLEV2_3_df_QC_bis_NO_MSD %>% 
  group_by(class) %>%
  summarize(N = n()) %>%
  mutate(Perc = round(N/sum(N)*100, digits = 1))

single_presencesLEV2_3_df_QC_bis_NO_MSD_sum_subclass <- single_presencesLEV2_3_df_QC_bis_NO_MSD %>% 
  group_by(subclass) %>%
  summarize(N = n()) %>%
  mutate(Perc = round(N/sum(N)*100, digits = 1))

single_presencesLEV2_3_df_QC_bis_NO_MSD_sum_level5 <- single_presencesLEV2_3_df_QC_bis_NO_MSD %>% 
  group_by(level5) %>%
  summarize(N = n()) %>%
  mutate(Perc = round(N/sum(N)*100, digits = 1))

single_presencesLEV2_3_df_QC_bis_NO_MSD_sum_level6 <- single_presencesLEV2_3_df_QC_bis_NO_MSD %>% 
  group_by(level6) %>%
  summarize(N = n()) %>%
  mutate(Perc = round(N/sum(N)*100, digits = 1))

single_presencesLEV2_3_df_QC_bis_NO_MSD_sum_level7 <- single_presencesLEV2_3_df_QC_bis_NO_MSD %>% 
  group_by(level7) %>%
  summarize(N = n()) %>%
  mutate(Perc = round(N/sum(N)*100, digits = 1))




######## 
## Sankey Diagram with categories of molecules:

# it's 0-indexed!
# that's why I add this_source - 1  and  this_target - 1

single_presencesLEV2_3_df_QC_bis_NO_MSD_noNAkd <- filter(single_presencesLEV2_3_df_QC_bis_NO_MSD, !is.na(kingdom))

all_names_QC_NO_MSD_QC <- c(unique(single_presencesLEV2_3_df_QC_bis_NO_MSD_noNAkd$kingdom),
                            unique(single_presencesLEV2_3_df_QC_bis_NO_MSD_noNAkd$superclass),
                            unique(single_presencesLEV2_3_df_QC_bis_NO_MSD_noNAkd$class),
                            unique(single_presencesLEV2_3_df_QC_bis_NO_MSD_noNAkd$subclass),
                            unique(single_presencesLEV2_3_df_QC_bis_NO_MSD_noNAkd$level5),
                            unique(single_presencesLEV2_3_df_QC_bis_NO_MSD_noNAkd$level6),
                            unique(single_presencesLEV2_3_df_QC_bis_NO_MSD_noNAkd$level7))


all_names_QC_NO_MSD_QC <- all_names_QC_NO_MSD_QC[!is.na(all_names_QC_NO_MSD_QC)]

##a check: this MUST be FALSE (no duplicated in the names of categories!)
any(duplicated(all_names_QC_NO_MSD_QC))

all_links_QC_NO_MSD_QC <- data.frame(source = integer(length = 0),
                                     target = integer(length = 0),
                                     value = double(length = 0))


for (i in 1:length(pull(single_presencesLEV2_3_df_QC_bis_NO_MSD_noNAkd, 1))) {
  this_source <- which(all_names_QC_NO_MSD_QC == single_presencesLEV2_3_df_QC_bis_NO_MSD_noNAkd$kingdom[i])
  
  this_target <- which(all_names_QC_NO_MSD_QC == single_presencesLEV2_3_df_QC_bis_NO_MSD_noNAkd$superclass[i])
  
  this_value <- 1
  
  if (length(this_source) == 1 & length(this_target) == 1) {
    
    this_source <- this_source - 1
    
    this_target <- this_target - 1
    
    this_df <- data.frame(source = as.integer(this_source),
                          target = as.integer(this_target),
                          value = as.double(this_value))
    
    all_links_QC_NO_MSD_QC <- rbind(all_links_QC_NO_MSD_QC, this_df)
    
  }
}


for (i in 1:length(pull(single_presencesLEV2_3_df_QC_bis_NO_MSD_noNAkd, 1))) {
  this_source <- which(all_names_QC_NO_MSD_QC == single_presencesLEV2_3_df_QC_bis_NO_MSD_noNAkd$superclass[i])
  
  this_target <-which(all_names_QC_NO_MSD_QC == single_presencesLEV2_3_df_QC_bis_NO_MSD_noNAkd$class[i])
  
  this_value <- 1
  
  if (length(this_source) == 1 & length(this_target) == 1) {
    
    this_source <- this_source - 1
    
    this_target <- this_target - 1
    
    this_df <- data.frame(source = as.integer(this_source),
                          target = as.integer(this_target),
                          value = as.double(this_value))
    
    all_links_QC_NO_MSD_QC <- rbind(all_links_QC_NO_MSD_QC, this_df)
    
  }
}

for (i in 1:length(pull(single_presencesLEV2_3_df_QC_bis_NO_MSD_noNAkd, 1))) {
  this_source <- which(all_names_QC_NO_MSD_QC == single_presencesLEV2_3_df_QC_bis_NO_MSD_noNAkd$class[i])
  
  this_target <-which(all_names_QC_NO_MSD_QC == single_presencesLEV2_3_df_QC_bis_NO_MSD_noNAkd$subclass[i])
  
  this_value <- 1
  
  if (length(this_source) == 1 & length(this_target) == 1) {
    
    this_source <- this_source - 1
    
    this_target <- this_target - 1
    
    this_df <- data.frame(source = as.integer(this_source),
                          target = as.integer(this_target),
                          value = as.double(this_value))
    
    all_links_QC_NO_MSD_QC <- rbind(all_links_QC_NO_MSD_QC, this_df)
    
  }
}


for (i in 1:length(pull(single_presencesLEV2_3_df_QC_bis_NO_MSD_noNAkd, 1))) {
  this_source <- which(all_names_QC_NO_MSD_QC == single_presencesLEV2_3_df_QC_bis_NO_MSD_noNAkd$subclass[i])
  
  this_target <-which(all_names_QC_NO_MSD_QC == single_presencesLEV2_3_df_QC_bis_NO_MSD_noNAkd$level5[i])
  
  this_value <- 1
  
  if (length(this_source) == 1 & length(this_target) == 1) {
    
    this_source <- this_source - 1
    
    this_target <- this_target - 1
    
    this_df <- data.frame(source = as.integer(this_source),
                          target = as.integer(this_target),
                          value = as.double(this_value))
    
    all_links_QC_NO_MSD_QC <- rbind(all_links_QC_NO_MSD_QC, this_df)
    
  }
}


for (i in 1:length(pull(single_presencesLEV2_3_df_QC_bis_NO_MSD_noNAkd, 1))) {
  this_source <- which(all_names_QC_NO_MSD_QC == single_presencesLEV2_3_df_QC_bis_NO_MSD_noNAkd$level5[i])
  
  this_target <-which(all_names_QC_NO_MSD_QC == single_presencesLEV2_3_df_QC_bis_NO_MSD_noNAkd$level6[i])
  
  this_value <- 1
  
  if (length(this_source) == 1 & length(this_target) == 1) {
    
    this_source <- this_source - 1
    
    this_target <- this_target - 1
    
    this_df <- data.frame(source = as.integer(this_source),
                          target = as.integer(this_target),
                          value = as.double(this_value))
    
    all_links_QC_NO_MSD_QC <- rbind(all_links_QC_NO_MSD_QC, this_df)
    
  }
}


for (i in 1:length(pull(single_presencesLEV2_3_df_QC_bis_NO_MSD_noNAkd, 1))) {
  this_source <- which(all_names_QC_NO_MSD_QC == single_presencesLEV2_3_df_QC_bis_NO_MSD_noNAkd$level6[i])
  
  this_target <-which(all_names_QC_NO_MSD_QC == single_presencesLEV2_3_df_QC_bis_NO_MSD_noNAkd$level7[i])
  
  this_value <- 1
  
  if (length(this_source) == 1 & length(this_target) == 1) {
    
    this_source <- this_source - 1
    
    this_target <- this_target - 1
    
    this_df <- data.frame(source = as.integer(this_source),
                          target = as.integer(this_target),
                          value = as.double(this_value))
    
    all_links_QC_NO_MSD_QC <- rbind(all_links_QC_NO_MSD_QC, this_df)
    
  }
}





all_links_QC_NO_MSD_QC <- mutate(all_links_QC_NO_MSD_QC, from_to = paste("from", source, "to", target))

Unique_Duplicated_QC_NO_MSD <- unique(all_links_QC_NO_MSD_QC$from_to[which(duplicated(all_links_QC_NO_MSD_QC$from_to))])

for (d in Unique_Duplicated_QC_NO_MSD) {
  df_fil <- filter(all_links_QC_NO_MSD_QC, from_to == d)
  
  this_added_value <- length(df_fil$from_to)
  
  all_links_QC_NO_MSD_QC[which(all_links_QC_NO_MSD_QC$from_to==d)[1],"value"] <- this_added_value
}

all_links_QC_NO_MSD_QC <- filter(all_links_QC_NO_MSD_QC, !duplicated(from_to))




List_for_Sankey_Diagram_QC_NO_MSD <- list(nodes = data.frame(name = all_names_QC_NO_MSD_QC),
                                          links = all_links_QC_NO_MSD_QC)


My_sankeyNetwork_QC_NO_MSD <- sankeyNetwork(Links = List_for_Sankey_Diagram_QC_NO_MSD$links, Nodes = List_for_Sankey_Diagram_QC_NO_MSD$nodes,
                                            Source = "source",
                                            Target = "target",
                                            Value = "value",
                                            NodeID = "name",
                                            fontSize = 12, nodeWidth = 30,
                                            sinksRight = FALSE)
saveNetwork(My_sankeyNetwork_QC_NO_MSD, "My_sankeyNetwork_QC_NO_MSD.html")
webshot::webshot("My_sankeyNetwork_QC_NO_MSD.html","My_sankeyNetwork_QC_NO_MSD.png", vwidth = 1600, vheight = 800)





### similarly but only those contained in MS-dial

single_presencesLEV2_3_df_QC_bis_ONLY_MSD <- single_presencesLEV2_3_df_QC_bis[which(single_presencesLEV2_3_df_QC_bis$present_in_MSDial == TRUE),]



### biologic significance:
# grouping:

single_presencesLEV2_3_df_QC_bis_ONLY_MSD_sum_kingdom <- single_presencesLEV2_3_df_QC_bis_ONLY_MSD %>% 
  group_by(kingdom) %>%
  summarize(N = n()) %>%
  mutate(Perc = round(N/sum(N)*100, digits = 1))

single_presencesLEV2_3_df_QC_bis_ONLY_MSD_sum_superclass <- single_presencesLEV2_3_df_QC_bis_ONLY_MSD %>% 
  group_by(superclass) %>%
  summarize(N = n()) %>%
  mutate(Perc = round(N/sum(N)*100, digits = 1))

single_presencesLEV2_3_df_QC_bis_ONLY_MSD_sum_class <- single_presencesLEV2_3_df_QC_bis_ONLY_MSD %>% 
  group_by(class) %>%
  summarize(N = n()) %>%
  mutate(Perc = round(N/sum(N)*100, digits = 1))

single_presencesLEV2_3_df_QC_bis_ONLY_MSD_sum_subclass <- single_presencesLEV2_3_df_QC_bis_ONLY_MSD %>% 
  group_by(subclass) %>%
  summarize(N = n()) %>%
  mutate(Perc = round(N/sum(N)*100, digits = 1))

single_presencesLEV2_3_df_QC_bis_ONLY_MSD_sum_level5 <- single_presencesLEV2_3_df_QC_bis_ONLY_MSD %>% 
  group_by(level5) %>%
  summarize(N = n()) %>%
  mutate(Perc = round(N/sum(N)*100, digits = 1))

single_presencesLEV2_3_df_QC_bis_ONLY_MSD_sum_level6 <- single_presencesLEV2_3_df_QC_bis_ONLY_MSD %>% 
  group_by(level6) %>%
  summarize(N = n()) %>%
  mutate(Perc = round(N/sum(N)*100, digits = 1))

single_presencesLEV2_3_df_QC_bis_ONLY_MSD_sum_level7 <- single_presencesLEV2_3_df_QC_bis_ONLY_MSD %>% 
  group_by(level7) %>%
  summarize(N = n()) %>%
  mutate(Perc = round(N/sum(N)*100, digits = 1))





######## 
## Sankey Diagram with categories of molecules:

# it's 0-indexed!
# that's why I add this_source - 1  and  this_target - 1

single_presencesLEV2_3_df_QC_bis_ONLY_MSD_noNAkd <- filter(single_presencesLEV2_3_df_QC_bis_ONLY_MSD, !is.na(kingdom))

all_names_QC_ONLY_MSD_QC <- c(unique(single_presencesLEV2_3_df_QC_bis_ONLY_MSD_noNAkd$kingdom),
                              unique(single_presencesLEV2_3_df_QC_bis_ONLY_MSD_noNAkd$superclass),
                              unique(single_presencesLEV2_3_df_QC_bis_ONLY_MSD_noNAkd$class),
                              unique(single_presencesLEV2_3_df_QC_bis_ONLY_MSD_noNAkd$subclass),
                              unique(single_presencesLEV2_3_df_QC_bis_ONLY_MSD_noNAkd$level5),
                              unique(single_presencesLEV2_3_df_QC_bis_ONLY_MSD_noNAkd$level6),
                              unique(single_presencesLEV2_3_df_QC_bis_ONLY_MSD_noNAkd$level7))


all_names_QC_ONLY_MSD_QC <- all_names_QC_ONLY_MSD_QC[!is.na(all_names_QC_ONLY_MSD_QC)]

##a check: this MUST be FALSE (no duplicated in the names of categories!)
any(duplicated(all_names_QC_ONLY_MSD_QC))

all_links_QC_ONLY_MSD_QC <- data.frame(source = integer(length = 0),
                                       target = integer(length = 0),
                                       value = double(length = 0))


for (i in 1:length(pull(single_presencesLEV2_3_df_QC_bis_ONLY_MSD_noNAkd, 1))) {
  this_source <- which(all_names_QC_ONLY_MSD_QC == single_presencesLEV2_3_df_QC_bis_ONLY_MSD_noNAkd$kingdom[i])
  
  this_target <- which(all_names_QC_ONLY_MSD_QC == single_presencesLEV2_3_df_QC_bis_ONLY_MSD_noNAkd$superclass[i])
  
  this_value <- 1
  
  if (length(this_source) == 1 & length(this_target) == 1) {
    
    this_source <- this_source - 1
    
    this_target <- this_target - 1
    
    this_df <- data.frame(source = as.integer(this_source),
                          target = as.integer(this_target),
                          value = as.double(this_value))
    
    all_links_QC_ONLY_MSD_QC <- rbind(all_links_QC_ONLY_MSD_QC, this_df)
    
  }
}


for (i in 1:length(pull(single_presencesLEV2_3_df_QC_bis_ONLY_MSD_noNAkd, 1))) {
  this_source <- which(all_names_QC_ONLY_MSD_QC == single_presencesLEV2_3_df_QC_bis_ONLY_MSD_noNAkd$superclass[i])
  
  this_target <-which(all_names_QC_ONLY_MSD_QC == single_presencesLEV2_3_df_QC_bis_ONLY_MSD_noNAkd$class[i])
  
  this_value <- 1
  
  if (length(this_source) == 1 & length(this_target) == 1) {
    
    this_source <- this_source - 1
    
    this_target <- this_target - 1
    
    this_df <- data.frame(source = as.integer(this_source),
                          target = as.integer(this_target),
                          value = as.double(this_value))
    
    all_links_QC_ONLY_MSD_QC <- rbind(all_links_QC_ONLY_MSD_QC, this_df)
    
  }
}

for (i in 1:length(pull(single_presencesLEV2_3_df_QC_bis_ONLY_MSD_noNAkd, 1))) {
  this_source <- which(all_names_QC_ONLY_MSD_QC == single_presencesLEV2_3_df_QC_bis_ONLY_MSD_noNAkd$class[i])
  
  this_target <-which(all_names_QC_ONLY_MSD_QC == single_presencesLEV2_3_df_QC_bis_ONLY_MSD_noNAkd$subclass[i])
  
  this_value <- 1
  
  if (length(this_source) == 1 & length(this_target) == 1) {
    
    this_source <- this_source - 1
    
    this_target <- this_target - 1
    
    this_df <- data.frame(source = as.integer(this_source),
                          target = as.integer(this_target),
                          value = as.double(this_value))
    
    all_links_QC_ONLY_MSD_QC <- rbind(all_links_QC_ONLY_MSD_QC, this_df)
    
  }
}


for (i in 1:length(pull(single_presencesLEV2_3_df_QC_bis_ONLY_MSD_noNAkd, 1))) {
  this_source <- which(all_names_QC_ONLY_MSD_QC == single_presencesLEV2_3_df_QC_bis_ONLY_MSD_noNAkd$subclass[i])
  
  this_target <-which(all_names_QC_ONLY_MSD_QC == single_presencesLEV2_3_df_QC_bis_ONLY_MSD_noNAkd$level5[i])
  
  this_value <- 1
  
  if (length(this_source) == 1 & length(this_target) == 1) {
    
    this_source <- this_source - 1
    
    this_target <- this_target - 1
    
    this_df <- data.frame(source = as.integer(this_source),
                          target = as.integer(this_target),
                          value = as.double(this_value))
    
    all_links_QC_ONLY_MSD_QC <- rbind(all_links_QC_ONLY_MSD_QC, this_df)
    
  }
}


for (i in 1:length(pull(single_presencesLEV2_3_df_QC_bis_ONLY_MSD_noNAkd, 1))) {
  this_source <- which(all_names_QC_ONLY_MSD_QC == single_presencesLEV2_3_df_QC_bis_ONLY_MSD_noNAkd$level5[i])
  
  this_target <-which(all_names_QC_ONLY_MSD_QC == single_presencesLEV2_3_df_QC_bis_ONLY_MSD_noNAkd$level6[i])
  
  this_value <- 1
  
  if (length(this_source) == 1 & length(this_target) == 1) {
    
    this_source <- this_source - 1
    
    this_target <- this_target - 1
    
    this_df <- data.frame(source = as.integer(this_source),
                          target = as.integer(this_target),
                          value = as.double(this_value))
    
    all_links_QC_ONLY_MSD_QC <- rbind(all_links_QC_ONLY_MSD_QC, this_df)
    
  }
}


for (i in 1:length(pull(single_presencesLEV2_3_df_QC_bis_ONLY_MSD_noNAkd, 1))) {
  this_source <- which(all_names_QC_ONLY_MSD_QC == single_presencesLEV2_3_df_QC_bis_ONLY_MSD_noNAkd$level6[i])
  
  this_target <-which(all_names_QC_ONLY_MSD_QC == single_presencesLEV2_3_df_QC_bis_ONLY_MSD_noNAkd$level7[i])
  
  this_value <- 1
  
  if (length(this_source) == 1 & length(this_target) == 1) {
    
    this_source <- this_source - 1
    
    this_target <- this_target - 1
    
    this_df <- data.frame(source = as.integer(this_source),
                          target = as.integer(this_target),
                          value = as.double(this_value))
    
    all_links_QC_ONLY_MSD_QC <- rbind(all_links_QC_ONLY_MSD_QC, this_df)
    
  }
}





all_links_QC_ONLY_MSD_QC <- mutate(all_links_QC_ONLY_MSD_QC, from_to = paste("from", source, "to", target))

Unique_Duplicated_QC_ONLY_MSD <- unique(all_links_QC_ONLY_MSD_QC$from_to[which(duplicated(all_links_QC_ONLY_MSD_QC$from_to))])

for (d in Unique_Duplicated_QC_ONLY_MSD) {
  df_fil <- filter(all_links_QC_ONLY_MSD_QC, from_to == d)
  
  this_added_value <- length(df_fil$from_to)
  
  all_links_QC_ONLY_MSD_QC[which(all_links_QC_ONLY_MSD_QC$from_to==d)[1],"value"] <- this_added_value
}

all_links_QC_ONLY_MSD_QC <- filter(all_links_QC_ONLY_MSD_QC, !duplicated(from_to))




List_for_Sankey_Diagram_QC_ONLY_MSD <- list(nodes = data.frame(name = all_names_QC_ONLY_MSD_QC),
                                            links = all_links_QC_ONLY_MSD_QC)


My_sankeyNetwork_QC_ONLY_MSD <- sankeyNetwork(Links = List_for_Sankey_Diagram_QC_ONLY_MSD$links, Nodes = List_for_Sankey_Diagram_QC_ONLY_MSD$nodes,
                                              Source = "source",
                                              Target = "target",
                                              Value = "value",
                                              NodeID = "name",
                                              fontSize = 12, nodeWidth = 30,
                                              sinksRight = FALSE)
saveNetwork(My_sankeyNetwork_QC_ONLY_MSD, "My_sankeyNetwork_QC_ONLY_MSD.html")
webshot::webshot("My_sankeyNetwork_QC_ONLY_MSD.html","My_sankeyNetwork_QC_ONLY_MSD.png", vwidth = 3400, vheight = 2600)









#####
#######
### doing everything again considering only 2a levels:


single_LEV2a_df_QC <- filter(single_LEV2_3_df_QC, AnnotationLevel == "2a")



single_presencesLEV2a_df_QC <- tibble(Name_to_use = unique(single_LEV2a_df_QC$Name_to_use),
                                      present_in_HILICPOS = as.logical(NA), present_in_RPLCNEG = as.logical(NA), present_in_scheme1 = as.logical(NA), present_in_scheme2 = as.logical(NA),
                                      present_in_MSDial = as.logical(NA), present_in_patRoon = as.logical(NA), present_in_PubChemLite = as.logical(NA), present_in_PublicMSP = as.logical(NA), present_in_WormJamExpanded = as.logical(NA),
                                      identifier = as.character(NA), compoundName = as.character(NA), neutral_formula = as.character(NA), InChIKey = as.character(NA),
                                      SMILES = as.character(NA), Title = as.character(NA), identifiers_with_compoundName_instead_null = as.character(NA),
                                      kingdom = as.character(NA), superclass = as.character(NA), class = as.character(NA),
                                      subclass = as.character(NA), level5 = as.character(NA), level6 = as.character(NA), level7 = as.character(NA))

compound_specification_columns <- c("identifier", "compoundName", "neutral_formula", "InChIKey", "SMILES", "Title", "identifiers_with_compoundName_instead_null", "kingdom", "superclass", "class", "subclass", "level5", "level6", "level7")



for (i in 1:length(single_presencesLEV2a_df_QC$Name_to_use)) {
  
  this_Name_to_use <- single_presencesLEV2a_df_QC$Name_to_use[i]
  
  single_LEV2a_df_QC_fill <- filter(single_LEV2a_df_QC, Name_to_use == this_Name_to_use)
  
  
  if (!all(map_lgl(single_LEV2a_df_QC_fill[, "neutral_formula"],
                   ~ (all(.x == .x[1])) | any(is.na(.x))))) {warning(paste0("compound specifications different for ", this_Name_to_use))}
  
  
  single_presencesLEV2a_df_QC[i, compound_specification_columns] <- single_LEV2a_df_QC_fill[1, compound_specification_columns]
  
  if (any(grepl("HILICPOS", single_LEV2a_df_QC_fill$type_analysis))) { single_presencesLEV2a_df_QC$present_in_HILICPOS[i] <- TRUE } else { single_presencesLEV2a_df_QC$present_in_HILICPOS[i] <- FALSE }
  if (any(grepl("RPLCNEG", single_LEV2a_df_QC_fill$type_analysis))) { single_presencesLEV2a_df_QC$present_in_RPLCNEG[i] <- TRUE } else { single_presencesLEV2a_df_QC$present_in_RPLCNEG[i] <- FALSE }
  if (any(single_LEV2a_df_QC_fill$present_in_scheme1)) { single_presencesLEV2a_df_QC$present_in_scheme1[i] <- TRUE } else { single_presencesLEV2a_df_QC$present_in_scheme1[i] <- FALSE }
  if (any(single_LEV2a_df_QC_fill$present_in_scheme2)) { single_presencesLEV2a_df_QC$present_in_scheme2[i] <- TRUE } else { single_presencesLEV2a_df_QC$present_in_scheme2[i] <- FALSE }
  if (any(grepl("patRoon", single_LEV2a_df_QC_fill$type_analysis))) { single_presencesLEV2a_df_QC$present_in_patRoon[i] <- TRUE } else { single_presencesLEV2a_df_QC$present_in_patRoon[i] <- FALSE }
  if (any(grepl("MSDial", single_LEV2a_df_QC_fill$type_analysis))) { single_presencesLEV2a_df_QC$present_in_MSDial[i] <- TRUE } else { single_presencesLEV2a_df_QC$present_in_MSDial[i] <- FALSE }
  if (any(grepl("PubChemLite", single_LEV2a_df_QC_fill$type_analysis))) { single_presencesLEV2a_df_QC$present_in_PubChemLite[i] <- TRUE } else { single_presencesLEV2a_df_QC$present_in_PubChemLite[i] <- FALSE }
  if (any(grepl("PublicMSP", single_LEV2a_df_QC_fill$type_analysis))) { single_presencesLEV2a_df_QC$present_in_PublicMSP[i] <- TRUE } else { single_presencesLEV2a_df_QC$present_in_PublicMSP[i] <- FALSE }
  if (any(grepl("WormJamExpanded", single_LEV2a_df_QC_fill$type_analysis))) { single_presencesLEV2a_df_QC$present_in_WormJamExpanded[i] <- TRUE } else { single_presencesLEV2a_df_QC$present_in_WormJamExpanded[i] <- FALSE }
  
}


# upset plot

# only annotation tools e chromatographic runs

single_presencesLEV2a_df_QC_per_upsetb1 <- tibble(Name_to_use = single_presencesLEV2a_df_QC$Name_to_use,
                                                  MSDial_PublicMSP = as.integer(NA),
                                                  MSDial_WormJamExpanded = as.integer(NA),
                                                  patRoon_PubChemLite = as.integer(NA),
                                                  patRoon_WormJamExpanded = as.integer(NA))


for (i in 1:length(pull(single_presencesLEV2a_df_QC_per_upsetb1, 1))) {
  
  if (single_presencesLEV2a_df_QC$present_in_MSDial[i] & single_presencesLEV2a_df_QC$present_in_PublicMSP[i]) {
    single_presencesLEV2a_df_QC_per_upsetb1$MSDial_PublicMSP[i] <- 1
  } else {
    single_presencesLEV2a_df_QC_per_upsetb1$MSDial_PublicMSP[i] <- 0
  }
  
  if (single_presencesLEV2a_df_QC$present_in_MSDial[i] & single_presencesLEV2a_df_QC$present_in_WormJamExpanded[i]) {
    single_presencesLEV2a_df_QC_per_upsetb1$MSDial_WormJamExpanded[i] <- 1
  } else {
    single_presencesLEV2a_df_QC_per_upsetb1$MSDial_WormJamExpanded[i] <- 0
  }
  
  if (single_presencesLEV2a_df_QC$present_in_patRoon [i] & single_presencesLEV2a_df_QC$present_in_PubChemLite [i]) {
    single_presencesLEV2a_df_QC_per_upsetb1$patRoon_PubChemLite[i] <- 1
  } else {
    single_presencesLEV2a_df_QC_per_upsetb1$patRoon_PubChemLite[i] <- 0
  }
  
  if (single_presencesLEV2a_df_QC$present_in_patRoon [i] & single_presencesLEV2a_df_QC$present_in_WormJamExpanded[i]) {
    single_presencesLEV2a_df_QC_per_upsetb1$patRoon_WormJamExpanded[i] <- 1
  } else {
    single_presencesLEV2a_df_QC_per_upsetb1$patRoon_WormJamExpanded[i] <- 0
  }
  
  if (sum(single_presencesLEV2a_df_QC_per_upsetb1$MSDial_PublicMSP[i],
          single_presencesLEV2a_df_QC_per_upsetb1$MSDial_WormJamExpanded[i],
          single_presencesLEV2a_df_QC_per_upsetb1$patRoon_PubChemLite[i],
          single_presencesLEV2a_df_QC_per_upsetb1$patRoon_WormJamExpanded[i]) == 0) {stop(paste0(single_presencesLEV2a_df_QC_per_upsetb1$unique_identifiers_with_compoundName_instead_null[i],
                                                                                                 " is not present in at least one annotation method"))}
  
}

single_presencesLEV2a_df_QC_per_upsetb1 <- mutate_if(single_presencesLEV2a_df_QC_per_upsetb1, is.double, as.integer)
glimpse(single_presencesLEV2a_df_QC_per_upsetb1)

#upset(single_presencesLEV2a_df_QC_per_upsetb1)




# it throws an error that I don't know how to solve.. Error in xtfrm.data.frame(x) : cannot xtfrm data frames

# I'll try with the list:


single_schemsep_LEV2a_df_QC_list_for_upsetb1 <- list(MSDial_PublicMSP = filter(single_presencesLEV2a_df_QC, present_in_MSDial == TRUE & present_in_PublicMSP == TRUE)$Name_to_use,
                                                     MSDial_WormJamExpanded = filter(single_presencesLEV2a_df_QC, present_in_MSDial == TRUE & present_in_WormJamExpanded == TRUE)$Name_to_use,
                                                     patRoon_PubChemLite = filter(single_presencesLEV2a_df_QC, present_in_patRoon == TRUE & present_in_PubChemLite == TRUE)$Name_to_use,
                                                     patRoon_WormJamExpanded = filter(single_presencesLEV2a_df_QC, present_in_patRoon == TRUE & present_in_WormJamExpanded == TRUE)$Name_to_use)


map_dbl(single_schemsep_LEV2a_df_QC_list_for_upsetb1, length)


png("upsetb1_only2a.png", width = 800, height = 600)
upset(fromList(single_schemsep_LEV2a_df_QC_list_for_upsetb1), text.scale = 1.8)
dev.off()






# Eulero-Venn:


new_euler_graph_b01_only2a <- ggvenn(single_schemsep_LEV2a_df_QC_list_for_upsetb1,
                                     c("MSDial_PublicMSP", "MSDial_WormJamExpanded", "patRoon_PubChemLite", "patRoon_WormJamExpanded"),
                                     show_elements = FALSE, show_percentage = TRUE, label_sep = "\n", text_size = 6) + 
  theme(panel.background = element_rect(fill = "white"))
ggsave(filename = "new_euler_graph_b01_only2a.png", plot = new_euler_graph_b01_only2a,  width = 25, height = 25, units = "in")



new_euler_graph_b02_only2a <- ggvenn(single_schemsep_LEV2a_df_QC_list_for_upsetb1,
                                     c("MSDial_PublicMSP", "MSDial_WormJamExpanded", "patRoon_PubChemLite", "patRoon_WormJamExpanded"),
                                     show_elements = FALSE, show_percentage = FALSE, label_sep = "\n", text_size = 6) + 
  theme(panel.background = element_rect(fill = "white"))
ggsave(filename = "new_euler_graph_b02_only2a.png", plot = new_euler_graph_b02_only2a,  width = 25, height = 25, units = "in")








### biologic significance:
# grouping:

single_presencesLEV2a_df_QC_sum_kingdom <- single_presencesLEV2a_df_QC %>% 
  group_by(kingdom) %>%
  summarize(N = n()) %>%
  mutate(Perc = round(N/sum(N)*100, digits = 1))

single_presencesLEV2a_df_QC_sum_superclass <- single_presencesLEV2a_df_QC %>% 
  group_by(superclass) %>%
  summarize(N = n()) %>%
  mutate(Perc = round(N/sum(N)*100, digits = 1))

single_presencesLEV2a_df_QC_sum_class <- single_presencesLEV2a_df_QC %>% 
  group_by(class) %>%
  summarize(N = n()) %>%
  mutate(Perc = round(N/sum(N)*100, digits = 1))

single_presencesLEV2a_df_QC_sum_subclass <- single_presencesLEV2a_df_QC %>% 
  group_by(subclass) %>%
  summarize(N = n()) %>%
  mutate(Perc = round(N/sum(N)*100, digits = 1))

single_presencesLEV2a_df_QC_sum_level5 <- single_presencesLEV2a_df_QC %>% 
  group_by(level5) %>%
  summarize(N = n()) %>%
  mutate(Perc = round(N/sum(N)*100, digits = 1))

single_presencesLEV2a_df_QC_sum_level6 <- single_presencesLEV2a_df_QC %>% 
  group_by(level6) %>%
  summarize(N = n()) %>%
  mutate(Perc = round(N/sum(N)*100, digits = 1))

single_presencesLEV2a_df_QC_sum_level7 <- single_presencesLEV2a_df_QC %>% 
  group_by(level7) %>%
  summarize(N = n()) %>%
  mutate(Perc = round(N/sum(N)*100, digits = 1))





######## 
## Sankey Diagram with categories of molecules:

# it's 0-indexed!
# that's why I add this_source - 1  and  this_target - 1


single_presencesLEV2a_df_QC_noNAkd <- filter(single_presencesLEV2a_df_QC, !is.na(kingdom))


all_names_2QC <- c(unique(single_presencesLEV2a_df_QC_noNAkd$kingdom),
                   unique(single_presencesLEV2a_df_QC_noNAkd$superclass),
                   unique(single_presencesLEV2a_df_QC_noNAkd$class),
                   unique(single_presencesLEV2a_df_QC_noNAkd$subclass),
                   unique(single_presencesLEV2a_df_QC_noNAkd$level5),
                   unique(single_presencesLEV2a_df_QC_noNAkd$level6),
                   unique(single_presencesLEV2a_df_QC_noNAkd$level7))

all_names_2QC <- all_names_2QC[!is.na(all_names_2QC)]

##a check: this MUST be FALSE (no duplicated in the names of categories!)
any(duplicated(all_names_2QC))

all_links_2QC <- data.frame(source = integer(length = 0),
                            target = integer(length = 0),
                            value = double(length = 0))


for (i in 1:length(pull(single_presencesLEV2a_df_QC_noNAkd, 1))) {
  this_source <- which(all_names_2QC == single_presencesLEV2a_df_QC_noNAkd$kingdom[i])
  
  this_target <-which(all_names_2QC == single_presencesLEV2a_df_QC_noNAkd$superclass[i])
  
  this_value <- 1
  
  if (length(this_source) == 1 & length(this_target) == 1) {
    
    this_source <- this_source - 1
    
    this_target <- this_target - 1
    
    this_df <- data.frame(source = as.integer(this_source),
                          target = as.integer(this_target),
                          value = as.double(this_value))
    
    all_links_2QC <- rbind(all_links_2QC, this_df)
    
  }
}


for (i in 1:length(pull(single_presencesLEV2a_df_QC_noNAkd, 1))) {
  this_source <- which(all_names_2QC == single_presencesLEV2a_df_QC_noNAkd$superclass[i])
  
  this_target <-which(all_names_2QC == single_presencesLEV2a_df_QC_noNAkd$class[i])
  
  this_value <- 1
  
  if (length(this_source) == 1 & length(this_target) == 1) {
    
    this_source <- this_source - 1
    
    this_target <- this_target - 1
    
    this_df <- data.frame(source = as.integer(this_source),
                          target = as.integer(this_target),
                          value = as.double(this_value))
    
    all_links_2QC <- rbind(all_links_2QC, this_df)
    
  }
}

for (i in 1:length(pull(single_presencesLEV2a_df_QC_noNAkd, 1))) {
  this_source <- which(all_names_2QC == single_presencesLEV2a_df_QC_noNAkd$class[i])
  
  this_target <-which(all_names_2QC == single_presencesLEV2a_df_QC_noNAkd$subclass[i])
  
  this_value <- 1
  
  if (length(this_source) == 1 & length(this_target) == 1) {
    
    this_source <- this_source - 1
    
    this_target <- this_target - 1
    
    this_df <- data.frame(source = as.integer(this_source),
                          target = as.integer(this_target),
                          value = as.double(this_value))
    
    all_links_2QC <- rbind(all_links_2QC, this_df)
    
  }
}


for (i in 1:length(pull(single_presencesLEV2a_df_QC_noNAkd, 1))) {
  this_source <- which(all_names_2QC == single_presencesLEV2a_df_QC_noNAkd$subclass[i])
  
  this_target <-which(all_names_2QC == single_presencesLEV2a_df_QC_noNAkd$level5[i])
  
  this_value <- 1
  
  if (length(this_source) == 1 & length(this_target) == 1) {
    
    this_source <- this_source - 1
    
    this_target <- this_target - 1
    
    this_df <- data.frame(source = as.integer(this_source),
                          target = as.integer(this_target),
                          value = as.double(this_value))
    
    all_links_2QC <- rbind(all_links_2QC, this_df)
    
  }
}


for (i in 1:length(pull(single_presencesLEV2a_df_QC_noNAkd, 1))) {
  this_source <- which(all_names_2QC == single_presencesLEV2a_df_QC_noNAkd$level5[i])
  
  this_target <-which(all_names_2QC == single_presencesLEV2a_df_QC_noNAkd$level6[i])
  
  this_value <- 1
  
  if (length(this_source) == 1 & length(this_target) == 1) {
    
    this_source <- this_source - 1
    
    this_target <- this_target - 1
    
    this_df <- data.frame(source = as.integer(this_source),
                          target = as.integer(this_target),
                          value = as.double(this_value))
    
    all_links_2QC <- rbind(all_links_2QC, this_df)
    
  }
}


for (i in 1:length(pull(single_presencesLEV2a_df_QC_noNAkd, 1))) {
  this_source <- which(all_names_2QC == single_presencesLEV2a_df_QC_noNAkd$level6[i])
  
  this_target <-which(all_names_2QC == single_presencesLEV2a_df_QC_noNAkd$level7[i])
  
  this_value <- 1
  
  if (length(this_source) == 1 & length(this_target) == 1) {
    
    this_source <- this_source - 1
    
    this_target <- this_target - 1
    
    this_df <- data.frame(source = as.integer(this_source),
                          target = as.integer(this_target),
                          value = as.double(this_value))
    
    all_links_2QC <- rbind(all_links_2QC, this_df)
    
  }
}





all_links_2QC <- mutate(all_links_2QC, from_to = paste("from", source, "to", target))

Unique_Duplicated_2QC <- unique(all_links_2QC$from_to[which(duplicated(all_links_2QC$from_to))])

for (d in Unique_Duplicated_2QC) {
  df_fil <- filter(all_links_2QC, from_to == d)
  
  this_added_value <- length(df_fil$from_to)
  
  all_links_2QC[which(all_links_2QC$from_to==d)[1],"value"] <- this_added_value
}

all_links_2QC <- filter(all_links_2QC, !duplicated(from_to))




List_for_Sankey_Diagram_2QC <- list(nodes = data.frame(name = all_names_2QC),
                                    links = all_links_2QC)


My_sankeyNetwork_2QC <- sankeyNetwork(Links = List_for_Sankey_Diagram_2QC$links, Nodes = List_for_Sankey_Diagram_2QC$nodes,
                                      Source = "source",
                                      Target = "target",
                                      Value = "value",
                                      NodeID = "name",
                                      fontSize = 12, nodeWidth = 50,
                                      sinksRight = FALSE)
saveNetwork(My_sankeyNetwork_2QC, "My_sankeyNetwork_2QC.html")
webshot::webshot("My_sankeyNetwork_2QC.html","My_sankeyNetwork_2QC.png", vwidth = 1600, vheight = 2200)















############
########################
####################################
################################################
############################################################
################################################
####################################
########################
############




# Retrieving details for levels 2 and 3 SINGIFICANT OF THE ANOVA:

all_lev2_3_list <- vector("list", length(all_data_list))

for (i in 1:length(all_data_list))  {
  
  DF <-  all_data_list[[i]]
  DF_name <- names(all_data_list)[i]
  
  if (grepl("ANOVA", DF_name)) {
    DF_2_3 <- filter(DF, AnnotationLevel %in% c("2a", "2b", "3a", "3b", "3c"))
    
    if (length(DF_2_3$AnnotationLevel) != 0) {
      if(grepl("MSDial", DF_name)) {
        
        DF_2_3_fill <- select(DF_2_3,
                              c("feature", "ANOVA_p_value", "ANOVA_p_value_FDR",                   
                                "N2_vs_VC40_pval", "N2_vs_VC1668_pval", "N2_vs_UA57_pval",                     
                                "N2_vs_BR5270_pval", "VC40_vs_VC1668_pval", "VC40_vs_UA57_pval",                   
                                "VC40_vs_BR5270_pval", "VC1668_vs_UA57_pval", "VC1668_vs_BR5270_pval",               
                                "UA57_vs_BR5270_pval", "N2_vs_VC40_higlow", "N2_vs_VC1668_higlow",                 
                                "N2_vs_UA57_higlow", "N2_vs_BR5270_higlow", "VC40_vs_VC1668_higlow",               
                                "VC40_vs_UA57_higlow", "VC40_vs_BR5270_higlow", "VC1668_vs_UA57_higlow",               
                                "VC1668_vs_BR5270_higlow", "UA57_vs_BR5270_higlow", "VC40_vs_BR5270",
                                "AnnotationLevel", "Metabolite_name", "INCHIKEY"))
        
        
        colnames(DF_2_3_fill) <- str_replace_all(colnames(DF_2_3_fill), "Metabolite_name", "compoundName")
        colnames(DF_2_3_fill) <- str_replace_all(colnames(DF_2_3_fill), "INCHIKEY", "InChIKey")
        
        TIBBLE_WITH_CID <- get_cid(DF_2_3_fill$InChIKey, from = "inchikey", match = "first")
        
        DF_2_3_fill <- add_column(DF_2_3_fill,
                                  .after = "InChIKey",
                                  identifier = as.character(TIBBLE_WITH_CID$cid))
        
        DF_2_3_fill <- mutate(DF_2_3_fill, identifier = ifelse(is.na(identifier), "null", identifier))
        
        
        write_csv(DF_2_3_fill, paste0("AnnLEV2_3_", DF_name, ".csv"))
        
        assign(paste0("AnnLEV2_3_", DF_name), DF_2_3_fill, envir = .GlobalEnv)
        
        all_lev2_3_list[[i]] <- DF_2_3_fill
        names(all_lev2_3_list)[i] <- paste0("AnnLEV2_3_", DF_name)
        
        
      } else if (grepl("patRoon", DF_name)) {
        
        DF_2_3_fill <- select(DF_2_3,
                              c("feature", "ANOVA_p_value", "ANOVA_p_value_FDR",                   
                                "N2_vs_VC40_pval", "N2_vs_VC1668_pval", "N2_vs_UA57_pval",                     
                                "N2_vs_BR5270_pval", "VC40_vs_VC1668_pval", "VC40_vs_UA57_pval",                   
                                "VC40_vs_BR5270_pval", "VC1668_vs_UA57_pval", "VC1668_vs_BR5270_pval",               
                                "UA57_vs_BR5270_pval", "N2_vs_VC40_higlow", "N2_vs_VC1668_higlow",                 
                                "N2_vs_UA57_higlow", "N2_vs_BR5270_higlow", "VC40_vs_VC1668_higlow",               
                                "VC40_vs_UA57_higlow", "VC40_vs_BR5270_higlow", "VC1668_vs_UA57_higlow",               
                                "VC1668_vs_BR5270_higlow", "UA57_vs_BR5270_higlow", "VC40_vs_BR5270",
                                "AnnotationLevel", "compoundName_topMONA", "InChIKey_topMONA", "identifier_topMONA"))
        
        colnames(DF_2_3_fill) <- str_remove_all(colnames(DF_2_3_fill), "_topMONA")
        
        DF_2_3_fill$identifier <- as.character(DF_2_3_fill$identifier)
        
        DF_2_3_fill <- mutate(DF_2_3_fill, identifier = ifelse(is.na(identifier), "null", identifier))
        
        write_csv(DF_2_3_fill, paste0("AnnLEV2_3_", DF_name, ".csv"))
        
        
        assign(paste0("AnnLEV2_3_", DF_name), DF_2_3_fill, envir = .GlobalEnv)
        
        all_lev2_3_list[[i]] <- DF_2_3_fill
        names(all_lev2_3_list)[i] <- paste0("AnnLEV2_3_", DF_name)
        
        
        
      }
    }
  }
}





## creating a single table with all 2_3 tables

all_lev2_3_df <- bind_rows(all_lev2_3_list, .id = "type_analysis")

write_tsv(all_lev2_3_df, "all_lev2_3_df.txt")


unique_identifier_2_3 <- unique(pull(all_lev2_3_df, "identifier"))

minFDR_2_3 <- vector("double", length(unique_identifier_2_3))

for (a in unique_identifier_2_3) {
  df_fil <- filter(all_lev2_3_df, identifier == a)
  minFDR_a <- min(pull(df_fil, "ANOVA_p_value_FDR"))
  
  minFDR_2_3[which(unique_identifier_2_3==a)] <- minFDR_a
}

minFDR_rank_2_3 <- rank(minFDR_2_3, ties.method = "first")


## creating separate table with info of the 2_3 compounds for each 

single_LEV2_3_list <- vector("list", length(unique_identifier_2_3))

for (a in unique_identifier_2_3) {
  The_rank_number <- ifelse(minFDR_rank_2_3[which(unique_identifier_2_3==a)]<10,
                            paste0("00", as.character(minFDR_rank_2_3[which(unique_identifier_2_3==a)])),
                            ifelse(minFDR_rank_2_3[which(unique_identifier_2_3==a)]<100,
                                   paste0("0", as.character(minFDR_rank_2_3[which(unique_identifier_2_3==a)])),
                                   as.character(minFDR_rank_2_3[which(unique_identifier_2_3==a)])))
  nametable <- paste0("singleLEV2_3_",
                      The_rank_number,
                      "_CID",
                      a)
  
  df_fil <- all_lev2_3_df %>%
    filter(identifier == a) %>%
    arrange(ANOVA_p_value_FDR)
  
  The_Title <- getPCdesc.title(query = a, from = "cid", timeout=10)$Title
  
  The_classification <- get_classification(df_fil$InChIKey[1])
  
  df_fil <- add_column(df_fil,
                       .after = "identifier",
                       Title = rep(The_Title, length(pull(df_fil, colnames(df_fil)[1])))) %>%
    add_column(.before = "identifier", rank_order = The_rank_number) %>%
    add_column(.after = "Title",
               kingdom = ifelse(!is.null(The_classification),
                                rep(The_classification@classification[["Classification"]][1], length(pull(df_fil, colnames(df_fil)[1]))),
                                rep(NA, length(pull(df_fil, colnames(df_fil)[1])))),
               superclass = ifelse(!is.null(The_classification),
                                   rep(The_classification@classification[["Classification"]][2], length(pull(df_fil, colnames(df_fil)[1]))),
                                   rep(NA, length(pull(df_fil, colnames(df_fil)[1])))),
               class = ifelse(!is.null(The_classification),
                              rep(The_classification@classification[["Classification"]][3], length(pull(df_fil, colnames(df_fil)[1]))),
                              rep(NA, length(pull(df_fil, colnames(df_fil)[1])))),
               subclass = ifelse(!is.null(The_classification),
                                 rep(The_classification@classification[["Classification"]][4], length(pull(df_fil, colnames(df_fil)[1]))),
                                 rep(NA, length(pull(df_fil, colnames(df_fil)[1])))),
               level5 = ifelse(!is.null(The_classification),
                               rep(The_classification@classification[["Classification"]][5], length(pull(df_fil, colnames(df_fil)[1]))),
                               rep(NA, length(pull(df_fil, colnames(df_fil)[1])))),
               level6 = ifelse(!is.null(The_classification),
                               rep(The_classification@classification[["Classification"]][6], length(pull(df_fil, colnames(df_fil)[1]))),
                               rep(NA, length(pull(df_fil, colnames(df_fil)[1])))),
               level7 = ifelse(!is.null(The_classification),
                               rep(The_classification@classification[["Classification"]][7], length(pull(df_fil, colnames(df_fil)[1]))),
                               rep(NA, length(pull(df_fil, colnames(df_fil)[1])))))
  
  
  
  df_fil <- relocate(df_fil,
                     feature, AnnotationLevel, rank_order,
                     identifier, InChIKey, Title, compoundName,
                     kingdom, superclass, class, subclass, level5, level6, level7,
                     type_analysis, ANOVA_p_value, ANOVA_p_value_FDR,
                     N2_vs_VC40_pval, N2_vs_VC1668_pval, N2_vs_UA57_pval, N2_vs_BR5270_pval, 
                     VC40_vs_VC1668_pval, VC40_vs_UA57_pval, VC40_vs_BR5270_pval, VC1668_vs_UA57_pval, VC1668_vs_BR5270_pval, UA57_vs_BR5270_pval, N2_vs_VC40_higlow, N2_vs_VC1668_higlow, 
                     N2_vs_UA57_higlow, N2_vs_BR5270_higlow, VC40_vs_VC1668_higlow, VC40_vs_UA57_higlow, VC40_vs_BR5270_higlow, VC1668_vs_UA57_higlow, VC1668_vs_BR5270_higlow, UA57_vs_BR5270_higlow, VC40_vs_BR5270)
  
  if (length(unique(df_fil$AnnotationLevel))!=1) {
    warning(paste0("there are different levels for the identifier ", The_rank_number, "_CID", a))
  }
  
  write_csv(df_fil, paste0(nametable, ".csv"))
  
  assign(nametable, df_fil, envir = .GlobalEnv)
  
  single_LEV2_3_list[[minFDR_rank_2_3[which(unique_identifier_2_3==a)]]] <- df_fil
  names(single_LEV2_3_list)[minFDR_rank_2_3[which(unique_identifier_2_3==a)]] <- nametable
  
}




all_lev2_3_df_arranged <- bind_rows(single_LEV2_3_list)

all_lev2_3_df_arranged$type_analysis <- str_remove_all(all_lev2_3_df_arranged$type_analysis, "AnnLEV2_3_")
all_lev2_3_df_arranged$type_analysis <- str_remove_all(all_lev2_3_df_arranged$type_analysis, "_ANOVA_sign_annotat_withLevels")
all_lev2_3_df_arranged$type_analysis <- str_remove_all(all_lev2_3_df_arranged$type_analysis, "_ANOVA_sign_annot_Levels")
all_lev2_3_df_arranged$type_analysis <- str_replace_all(all_lev2_3_df_arranged$type_analysis, "PCL", "PubChemLite")
all_lev2_3_df_arranged$type_analysis <- str_replace_all(all_lev2_3_df_arranged$type_analysis, "WJ", "WormJamExpanded")
all_lev2_3_df_arranged$type_analysis <- str_replace_all(all_lev2_3_df_arranged$type_analysis, "RPLCNEG_MSDial_scheme1_dry", "RPLCNEG_scheme1_MSDial_PublicMSP")
all_lev2_3_df_arranged$type_analysis <- str_replace_all(all_lev2_3_df_arranged$type_analysis, "HILICPOS_MSDial_scheme2_dry", "HILICPOS_scheme2_MSDial_PublicMSP")
all_lev2_3_df_arranged$type_analysis <- str_replace_all(all_lev2_3_df_arranged$type_analysis, "RPLCNEG_MSDial_scheme2_dry", "RPLCNEG_scheme2_MSDial_PublicMSP")
all_lev2_3_df_arranged$type_analysis <- str_replace_all(all_lev2_3_df_arranged$type_analysis, "HILICPOS_MSDial_scheme1_dry", "HILICPOS_scheme1_MSDial_PublicMSP")
all_lev2_3_df_arranged$type_analysis <- str_replace_all(all_lev2_3_df_arranged$type_analysis, "HILICPOS_patRoon_IPO_PubChemLite_scheme2_dry", "HILICPOS_scheme2_patRoon_PubChemLite")
all_lev2_3_df_arranged$type_analysis <- str_replace_all(all_lev2_3_df_arranged$type_analysis, "HILICPOS_patRoon_IPO_WormJamExpanded_scheme2_dry", "HILICPOS_scheme2_patRoon_WormJamExpanded")
all_lev2_3_df_arranged$type_analysis <- str_replace_all(all_lev2_3_df_arranged$type_analysis, "RPLCNEG_patRoon_IPO_PubChemLite_scheme2_dry", "RPLCNEG_scheme2_patRoon_PubChemLite")
all_lev2_3_df_arranged$type_analysis <- str_replace_all(all_lev2_3_df_arranged$type_analysis, "RPLCNEG_patRoon_IPO_WormJamExpanded_scheme2_dry", "RPLCNEG_scheme2_patRoon_WormJamExpanded")
all_lev2_3_df_arranged$type_analysis <- str_replace_all(all_lev2_3_df_arranged$type_analysis, "RPLCNEG_patRoon_IPO_PubChemLite_scheme1_dry", "RPLCNEG_scheme1_patRoon_PubChemLite")
all_lev2_3_df_arranged$type_analysis <- str_replace_all(all_lev2_3_df_arranged$type_analysis, "RPLCNEG_patRoon_IPO_WormJamExpanded_scheme1_dry", "RPLCNEG_scheme1_patRoon_WormJamExpanded")
all_lev2_3_df_arranged$type_analysis <- str_replace_all(all_lev2_3_df_arranged$type_analysis, "HILICPOS_patRoon_IPO_PubChemLite_scheme1_dry", "HILICPOS_scheme1_patRoon_PubChemLite")
all_lev2_3_df_arranged$type_analysis <- str_replace_all(all_lev2_3_df_arranged$type_analysis, "HILICPOS_patRoon_IPO_WormJamExpanded_scheme1_dry", "HILICPOS_scheme1_patRoon_WormJamExpanded")
all_lev2_3_df_arranged$type_analysis <- str_replace_all(all_lev2_3_df_arranged$type_analysis, "RPLCNEG_MSDial_WormJamExpanded_scheme1_dry", "RPLCNEG_scheme1_MSDial_WormJamExpanded")
all_lev2_3_df_arranged$type_analysis <- str_replace_all(all_lev2_3_df_arranged$type_analysis, "HILICPOS_MSDial_WormJamExpanded_scheme2_dry", "HILICPOS_scheme2_MSDial_WormJamExpanded")
all_lev2_3_df_arranged$type_analysis <- str_replace_all(all_lev2_3_df_arranged$type_analysis, "RPLCNEG_MSDial_WormJamExpanded_scheme2_dry", "RPLCNEG_scheme2_MSDial_WormJamExpanded")
all_lev2_3_df_arranged$type_analysis <- str_replace_all(all_lev2_3_df_arranged$type_analysis, "HILICPOS_MSDial_WormJamExpanded_scheme1_dry", "HILICPOS_scheme1_MSDial_WormJamExpanded")




write_tsv(all_lev2_3_df_arranged, "all_lev2_3_df_arranged.txt")



# cleaning names!!
all_lev2_3_df_arranged <- mutate(all_lev2_3_df_arranged,
                                 Name_to_use = as.character(NA))

for (i in 1:length(pull(all_lev2_3_df_arranged, 1))) {
  
  this_compoundName <- all_lev2_3_df_arranged$compoundName[i]
  this_title <- all_lev2_3_df_arranged$Title[i]
  
  if (is.na(this_title)) {
    this_name_to_use <- this_compoundName
  } else if (is.na(this_compoundName)) {
    this_name_to_use <- this_title
  } else if (grepl("CID", this_title)) {
    this_name_to_use <- this_compoundName
  } else if (nchar(this_title) < nchar(this_compoundName)) {
    this_name_to_use <- this_title
  } else {
    this_name_to_use <- this_compoundName
  }
  
  this_name_to_use <- sub("^L-", "", this_name_to_use)
  this_name_to_use <- sub("^D-", "", this_name_to_use)
  this_name_to_use <- sub("^DL-", "", this_name_to_use)
  
  if (nchar(gsub("[^A-Z]", "", this_name_to_use)) > 4) {
    this_name_to_use <- sub("^(\\w)(.*)", "\\1\\L\\2", this_name_to_use, perl = TRUE)
  }
  
  
  all_lev2_3_df_arranged$Name_to_use[i] <- this_name_to_use
}



all_lev2_3_df_arranged_perlev <- all_lev2_3_df_arranged
all_lev2_3_df_arranged_perlev$AnnotationLevel <- factor(all_lev2_3_df_arranged_perlev$AnnotationLevel, levels = c("2a", "2b", "3a", "3b", "3c"))

all_lev2_3_df_arranged_perlev <- arrange(all_lev2_3_df_arranged_perlev, AnnotationLevel)
write_tsv(all_lev2_3_df_arranged_perlev, "all_lev2_3_df_arranged_perlev.txt")






#  great, filtering just one line per compound, just to be able to make some general graphs later

all_lev2_3_df_one_compound <- all_lev2_3_df_arranged_perlev %>%
  mutate(NON_DUPLICATED = !duplicated(Name_to_use)) %>%
  filter(NON_DUPLICATED) %>%
  select(-NON_DUPLICATED)

write_tsv(all_lev2_3_df_one_compound, "all_lev2_3_df_one_compound.txt")



all_lev2_3_df_one_compound %>% group_by(AnnotationLevel) %>% summarize(N=n())




### biologic significance:
# grouping:

all_lev2_3_df_one_compound_sum_kingdom <- all_lev2_3_df_one_compound %>% 
  group_by(kingdom) %>%
  summarize(N = n()) %>%
  mutate(Perc = round(N/sum(N)*100, digits = 1))

all_lev2_3_df_one_compound_sum_superclass <- all_lev2_3_df_one_compound %>% 
  group_by(superclass) %>%
  summarize(N = n()) %>%
  mutate(Perc = round(N/sum(N)*100, digits = 1))

all_lev2_3_df_one_compound_sum_class <- all_lev2_3_df_one_compound %>% 
  group_by(class) %>%
  summarize(N = n()) %>%
  mutate(Perc = round(N/sum(N)*100, digits = 1))

all_lev2_3_df_one_compound_sum_subclass <- all_lev2_3_df_one_compound %>% 
  group_by(subclass) %>%
  summarize(N = n()) %>%
  mutate(Perc = round(N/sum(N)*100, digits = 1))

all_lev2_3_df_one_compound_sum_level5 <- all_lev2_3_df_one_compound %>% 
  group_by(level5) %>%
  summarize(N = n()) %>%
  mutate(Perc = round(N/sum(N)*100, digits = 1))

all_lev2_3_df_one_compound_sum_level6 <- all_lev2_3_df_one_compound %>% 
  group_by(level6) %>%
  summarize(N = n()) %>%
  mutate(Perc = round(N/sum(N)*100, digits = 1))

all_lev2_3_df_one_compound_sum_level7 <- all_lev2_3_df_one_compound %>% 
  group_by(level7) %>%
  summarize(N = n()) %>%
  mutate(Perc = round(N/sum(N)*100, digits = 1))


## without N2 significant:

all_lev2_3_df_without_N2 <- filter(all_lev2_3_df_one_compound, is.na(N2_vs_VC40_higlow) & is.na(N2_vs_VC1668_higlow) & is.na(N2_vs_UA57_higlow) & is.na(N2_vs_BR5270_higlow))

## VC40

# N2 > VC40

all_lev2_3_df_N2_VC40 <- filter(all_lev2_3_df_one_compound, grepl("N2 > VC40", N2_vs_VC40_higlow))

# N2 < VC40

all_lev2_3_df_VC40_N2 <- filter(all_lev2_3_df_one_compound, grepl("VC40 > N2", N2_vs_VC40_higlow))


## VC1668

# N2 > VC1668

all_lev2_3_df_N2_VC1668 <- filter(all_lev2_3_df_one_compound, grepl("N2 > VC1668", N2_vs_VC1668_higlow))

# N2 < VC1668

all_lev2_3_df_VC1668_N2 <- filter(all_lev2_3_df_one_compound, grepl("VC1668 > N2", N2_vs_VC1668_higlow))


## UA57

# N2 > UA57

all_lev2_3_df_N2_UA57 <- filter(all_lev2_3_df_one_compound, grepl("N2 > UA57", N2_vs_UA57_higlow))

# N2 < UA57

all_lev2_3_df_UA57_N2 <- filter(all_lev2_3_df_one_compound, grepl("UA57 > N2", N2_vs_UA57_higlow))


## BR5270

# N2 > BR5270

all_lev2_3_df_N2_BR5270 <- filter(all_lev2_3_df_one_compound, grepl("N2 > BR5270", N2_vs_BR5270_higlow))

# N2 < BR5270

all_lev2_3_df_BR5270_N2 <- filter(all_lev2_3_df_one_compound, grepl("BR5270 > N2", N2_vs_BR5270_higlow))




###
####
#########
######## 
## Sankey Diagram with categories of molecules:

# it's 0-indexed!
# that's why I add this_source - 1  and  this_target - 1




all_names <- c(unique(all_lev2_3_df_one_compound$kingdom),
               unique(all_lev2_3_df_one_compound$superclass),
               unique(all_lev2_3_df_one_compound$class),
               unique(all_lev2_3_df_one_compound$subclass),
               unique(all_lev2_3_df_one_compound$level5),
               unique(all_lev2_3_df_one_compound$level6),
               unique(all_lev2_3_df_one_compound$level7))

all_names <- all_names[!is.na(all_names)]

##a check: this MUST be FALSE (no duplicated in the names of categories!)
any(duplicated(all_names))

all_links <- data.frame(source = integer(length = 0),
                        target = integer(length = 0),
                        value = double(length = 0))


for (i in 1:length(pull(all_lev2_3_df_one_compound, colnames(all_lev2_3_df_one_compound)[1]))) {
  this_source <- which(all_names == all_lev2_3_df_one_compound$kingdom[i])
  
  this_target <-which(all_names == all_lev2_3_df_one_compound$superclass[i])
  
  this_value <- 1
  
  if (length(this_source) == 1 & length(this_target) == 1) {
    
    this_source <- this_source - 1
    
    this_target <- this_target - 1
    
    this_df <- data.frame(source = as.integer(this_source),
                          target = as.integer(this_target),
                          value = as.double(this_value))
    
    all_links <- rbind(all_links, this_df)
    
  }
}


for (i in 1:length(pull(all_lev2_3_df_one_compound, colnames(all_lev2_3_df_one_compound)[1]))) {
  this_source <- which(all_names == all_lev2_3_df_one_compound$superclass[i])
  
  this_target <-which(all_names == all_lev2_3_df_one_compound$class[i])
  
  this_value <- 1
  
  if (length(this_source) == 1 & length(this_target) == 1) {
    
    this_source <- this_source - 1
    
    this_target <- this_target - 1
    
    this_df <- data.frame(source = as.integer(this_source),
                          target = as.integer(this_target),
                          value = as.double(this_value))
    
    all_links <- rbind(all_links, this_df)
    
  }
}

for (i in 1:length(pull(all_lev2_3_df_one_compound, colnames(all_lev2_3_df_one_compound)[1]))) {
  this_source <- which(all_names == all_lev2_3_df_one_compound$class[i])
  
  this_target <-which(all_names == all_lev2_3_df_one_compound$subclass[i])
  
  this_value <- 1
  
  if (length(this_source) == 1 & length(this_target) == 1) {
    
    this_source <- this_source - 1
    
    this_target <- this_target - 1
    
    this_df <- data.frame(source = as.integer(this_source),
                          target = as.integer(this_target),
                          value = as.double(this_value))
    
    all_links <- rbind(all_links, this_df)
    
  }
}


for (i in 1:length(pull(all_lev2_3_df_one_compound, colnames(all_lev2_3_df_one_compound)[1]))) {
  this_source <- which(all_names == all_lev2_3_df_one_compound$subclass[i])
  
  this_target <-which(all_names == all_lev2_3_df_one_compound$level5[i])
  
  this_value <- 1
  
  if (length(this_source) == 1 & length(this_target) == 1) {
    
    this_source <- this_source - 1
    
    this_target <- this_target - 1
    
    this_df <- data.frame(source = as.integer(this_source),
                          target = as.integer(this_target),
                          value = as.double(this_value))
    
    all_links <- rbind(all_links, this_df)
    
  }
}


for (i in 1:length(pull(all_lev2_3_df_one_compound, colnames(all_lev2_3_df_one_compound)[1]))) {
  this_source <- which(all_names == all_lev2_3_df_one_compound$level5[i])
  
  this_target <-which(all_names == all_lev2_3_df_one_compound$level6[i])
  
  this_value <- 1
  
  if (length(this_source) == 1 & length(this_target) == 1) {
    
    this_source <- this_source - 1
    
    this_target <- this_target - 1
    
    this_df <- data.frame(source = as.integer(this_source),
                          target = as.integer(this_target),
                          value = as.double(this_value))
    
    all_links <- rbind(all_links, this_df)
    
  }
}


for (i in 1:length(pull(all_lev2_3_df_one_compound, colnames(all_lev2_3_df_one_compound)[1]))) {
  this_source <- which(all_names == all_lev2_3_df_one_compound$level6[i])
  
  this_target <-which(all_names == all_lev2_3_df_one_compound$level7[i])
  
  this_value <- 1
  
  if (length(this_source) == 1 & length(this_target) == 1) {
    
    this_source <- this_source - 1
    
    this_target <- this_target - 1
    
    this_df <- data.frame(source = as.integer(this_source),
                          target = as.integer(this_target),
                          value = as.double(this_value))
    
    all_links <- rbind(all_links, this_df)
    
  }
}





all_links <- mutate(all_links, from_to = paste("from", source, "to", target))

Unique_Duplicated <- unique(all_links$from_to[which(duplicated(all_links$from_to))])

for (d in Unique_Duplicated) {
  df_fil <- filter(all_links, from_to == d)
  
  this_added_value <- length(df_fil$from_to)
  
  all_links[which(all_links$from_to==d)[1],"value"] <- this_added_value
}

all_links <- filter(all_links, !duplicated(from_to))




List_for_Sankey_Diagram <- list(nodes = data.frame(name = all_names),
                                links = all_links)


My_sankeyNetwork <- sankeyNetwork(Links = List_for_Sankey_Diagram$links, Nodes = List_for_Sankey_Diagram$nodes,
                                  Source = "source",
                                  Target = "target",
                                  Value = "value",
                                  NodeID = "name",
                                  fontSize = 12, nodeWidth = 30,
                                  sinksRight = FALSE)
saveNetwork(My_sankeyNetwork, "My_sankeyNetwork.html")
webshot::webshot("My_sankeyNetwork.html","My_sankeyNetwork.png", vwidth = 1600, vheight = 900)






### upset plot!

List_for_Upset_plot <- list(RPLCNEG_scheme1_patRoon_PubChemLite = all_lev2_3_df_arranged$Name_to_use[all_lev2_3_df_arranged$type_analysis == "RPLCNEG_scheme1_patRoon_PubChemLite"],
                            RPLCNEG_scheme2_patRoon_PubChemLite = all_lev2_3_df_arranged$Name_to_use[all_lev2_3_df_arranged$type_analysis == "RPLCNEG_scheme2_patRoon_PubChemLite"],
                            RPLCNEG_scheme1_patRoon_WormJamExpanded = all_lev2_3_df_arranged$Name_to_use[all_lev2_3_df_arranged$type_analysis == "RPLCNEG_scheme1_patRoon_WormJamExpanded"],
                            RPLCNEG_scheme2_patRoon_WormJamExpanded = all_lev2_3_df_arranged$Name_to_use[all_lev2_3_df_arranged$type_analysis == "RPLCNEG_scheme2_patRoon_WormJamExpanded"],
                            
                            RPLCNEG_scheme1_MSDial_PublicMSP = all_lev2_3_df_arranged$Name_to_use[all_lev2_3_df_arranged$type_analysis == "RPLCNEG_scheme1_MSDial_PublicMSP"],
                            RPLCNEG_scheme2_MSDial_PublicMSP = all_lev2_3_df_arranged$Name_to_use[all_lev2_3_df_arranged$type_analysis == "RPLCNEG_scheme2_MSDial_PublicMSP"],
                            RPLCNEG_scheme1_MSDial_WormJamExpanded = all_lev2_3_df_arranged$Name_to_use[all_lev2_3_df_arranged$type_analysis == "RPLCNEG_scheme1_MSDial_WormJamExpanded"],
                            RPLCNEG_scheme2_MSDial_WormJamExpanded = all_lev2_3_df_arranged$Name_to_use[all_lev2_3_df_arranged$type_analysis == "RPLCNEG_scheme2_MSDial_WormJamExpanded"],
                            
                            HILICPOS_scheme1_patRoon_PubChemLite = all_lev2_3_df_arranged$Name_to_use[all_lev2_3_df_arranged$type_analysis == "HILICPOS_scheme1_patRoon_PubChemLite"],
                            HILICPOS_scheme2_patRoon_PubChemLite = all_lev2_3_df_arranged$Name_to_use[all_lev2_3_df_arranged$type_analysis == "HILICPOS_scheme2_patRoon_PubChemLite"],
                            HILICPOS_scheme1_patRoon_WormJamExpanded = all_lev2_3_df_arranged$Name_to_use[all_lev2_3_df_arranged$type_analysis == "HILICPOS_scheme1_patRoon_WormJamExpanded"],
                            HILICPOS_scheme2_patRoon_WormJamExpanded = all_lev2_3_df_arranged$Name_to_use[all_lev2_3_df_arranged$type_analysis == "HILICPOS_scheme2_patRoon_WormJamExpanded"],
                            
                            HILICPOS_scheme1_MSDial_PublicMSP = all_lev2_3_df_arranged$Name_to_use[all_lev2_3_df_arranged$type_analysis == "HILICPOS_scheme1_MSDial_PublicMSP"],
                            HILICPOS_scheme2_MSDial_PublicMSP = all_lev2_3_df_arranged$Name_to_use[all_lev2_3_df_arranged$type_analysis == "HILICPOS_scheme2_MSDial_PublicMSP"],
                            HILICPOS_scheme1_MSDial_WormJamExpanded = all_lev2_3_df_arranged$Name_to_use[all_lev2_3_df_arranged$type_analysis == "HILICPOS_scheme1_MSDial_WormJamExpanded"],
                            HILICPOS_scheme2_MSDial_WormJamExpanded = all_lev2_3_df_arranged$Name_to_use[all_lev2_3_df_arranged$type_analysis == "HILICPOS_scheme2_MSDial_WormJamExpanded"])



upset(fromList(List_for_Upset_plot),
      nsets = 1000000,
      nintersects = NA)


png("my_upset.png", width = 1400, height = 1200, res = 100)

upset(fromList(List_for_Upset_plot),
      nsets = 1000000,
      nintersects = NA)

dev.off()




####### Creation of a personalised graph!!
##### This will be used as Fig. 4!


all_lev2_3_df_arranged_perlev_forGFgraph <- mutate(all_lev2_3_df_arranged_perlev, molec = Name_to_use)


all_lev2_3_df_arranged_perlev_forGFgraph <- mutate(all_lev2_3_df_arranged_perlev_forGFgraph, molec_abbr = ifelse(nchar(molec)>70, paste0(substr(molec, 1, 70), "..."), molec))

all_lev2_3_df_arranged_perlev_forGFgraph$identifier <- factor(all_lev2_3_df_arranged_perlev_forGFgraph$identifier, levels = unique(all_lev2_3_df_arranged_perlev_forGFgraph$identifier))
all_lev2_3_df_arranged_perlev_forGFgraph$Title <- factor(all_lev2_3_df_arranged_perlev_forGFgraph$Title, levels = unique(all_lev2_3_df_arranged_perlev_forGFgraph$Title))
all_lev2_3_df_arranged_perlev_forGFgraph$Name_to_use <- factor(all_lev2_3_df_arranged_perlev_forGFgraph$Name_to_use, levels = unique(all_lev2_3_df_arranged_perlev_forGFgraph$Name_to_use))
all_lev2_3_df_arranged_perlev_forGFgraph$molec <- factor(all_lev2_3_df_arranged_perlev_forGFgraph$molec, levels = unique(all_lev2_3_df_arranged_perlev_forGFgraph$molec))
all_lev2_3_df_arranged_perlev_forGFgraph$molec_abbr <- factor(all_lev2_3_df_arranged_perlev_forGFgraph$molec_abbr, levels = unique(all_lev2_3_df_arranged_perlev_forGFgraph$molec_abbr))


all_lev2_3_df_arranged_perlev_forGFgraph$type_analysis <- factor(all_lev2_3_df_arranged_perlev_forGFgraph$type_analysis, levels = sort(unique(all_lev2_3_df_arranged_perlev_forGFgraph$type_analysis), decreasing = TRUE))



colours_named_vector <- c(RPLCNEG_scheme2_patRoon_WormJamExpanded = "dodgerblue1",
                          RPLCNEG_scheme1_patRoon_WormJamExpanded = "dodgerblue1",
                          RPLCNEG_scheme2_patRoon_PubChemLite = "deepskyblue1",
                          RPLCNEG_scheme1_patRoon_PubChemLite = "deepskyblue1",
                          
                          RPLCNEG_scheme2_MSDial_WormJamExpanded = "darkorchid4",
                          RPLCNEG_scheme1_MSDial_WormJamExpanded = "darkorchid4",
                          RPLCNEG_scheme2_MSDial_PublicMSP = "darkorchid1",
                          RPLCNEG_scheme1_MSDial_PublicMSP = "darkorchid1",
                          
                          
                          HILICPOS_scheme2_patRoon_WormJamExpanded ="red1",
                          HILICPOS_scheme1_patRoon_WormJamExpanded ="red1",
                          HILICPOS_scheme2_patRoon_PubChemLite = "orangered1",
                          HILICPOS_scheme1_patRoon_PubChemLite = "orangered1",
                          
                          HILICPOS_scheme2_MSDial_WormJamExpanded = "goldenrod4",
                          HILICPOS_scheme1_MSDial_WormJamExpanded = "goldenrod4",
                          HILICPOS_scheme2_MSDial_PublicMSP = "goldenrod1",
                          HILICPOS_scheme1_MSDial_PublicMSP = "goldenrod1")


colours_level_annotation <- c(`2a` = "olivedrab3",
                              `3a` = "gold2",
                              `3b` = "darkorange1",
                              `3c` = "brown2")

names_named_vector <- c(RPLCNEG_scheme2dry_patRoon_WormJamExpanded = "RPLCNEG scheme2 patRoon WormJamExpanded",
                        RPLCNEG_scheme1dry_patRoon_WormJamExpanded = "RPLCNEG scheme1 patRoon WormJamExpanded",
                        RPLCNEG_scheme2dry_patRoon_PubChemLite = "RPLCNEG scheme2 patRoon PubChemLite",
                        RPLCNEG_scheme1dry_patRoon_PubChemLite = "RPLCNEG scheme1 patRoon PubChemLite",
                        
                        RPLCNEG_scheme2dry_MSDial_WormJamExpanded = "RPLCNEG scheme2 MSDial WormJamExpanded",
                        RPLCNEG_scheme1dry_MSDial_WormJamExpanded = "RPLCNEG scheme1 MSDial WormJamExpanded",
                        RPLCNEG_scheme2dry_MSDial_PublicMSP = "RPLCNEG scheme2 MSDial_PublicMSP",
                        RPLCNEG_scheme1dry_MSDial_PublicMSP = "RPLCNEG scheme1 MSDial_PublicMSP",
                        
                        
                        HILICPOS_scheme2dry_patRoon_WormJamExpanded = "HILICPOS scheme2 patRoon WormJamExpanded",
                        HILICPOS_scheme1dry_patRoon_WormJamExpanded = "HILICPOS scheme1 patRoon WormJamExpanded",
                        HILICPOS_scheme2dry_patRoon_PubChemLite = "HILICPOS scheme2 patRoon PubChemLite",
                        HILICPOS_scheme1dry_patRoon_PubChemLite = "HILICPOS scheme1 patRoon PubChemLite",
                        
                        HILICPOS_scheme2dry_MSDial_WormJamExpanded = "HILICPOS scheme2 MSDial WormJamExpanded",
                        HILICPOS_scheme1dry_MSDial_WormJamExpanded = "HILICPOS scheme1 MSDial WormJamExpanded",
                        HILICPOS_scheme2dry_MSDial_PublicMSP = "HILICPOS scheme2 MSDial Public MSP",
                        HILICPOS_scheme1dry_MSDial_PublicMSP = "HILICPOS scheme1 MSDial Public MSP")


ggplot(all_lev2_3_df_arranged_perlev_forGFgraph, aes(x = molec_abbr, y = type_analysis)) +
  geom_point(aes(col = AnnotationLevel), shape = 15) + 
  scale_color_manual(values = colours_level_annotation) +
  ylab("Analyses type") +
  xlab("Annotated metabolites") +
  scale_x_discrete(position = "top") +
  scale_y_discrete(labels = names_named_vector) + 
  theme_light() +
  theme(axis.text.x = element_text(angle = 90, hjust = 0)) + 
  theme(legend.position = "none")
ggsave("my_GFgraph.png", width = 20.5, height = 7, units = "in")



#######again my graph, but transposed:
all_lev2_3_df_arranged_perlev_forGFgrapht <- all_lev2_3_df_arranged_perlev_forGFgraph
all_lev2_3_df_arranged_perlev_forGFgrapht$molec_abbr <- droplevels(all_lev2_3_df_arranged_perlev_forGFgrapht$molec_abbr)
all_lev2_3_df_arranged_perlev_forGFgrapht$molec_abbr <- factor(all_lev2_3_df_arranged_perlev_forGFgrapht$molec_abbr, levels = rev(levels(all_lev2_3_df_arranged_perlev_forGFgrapht$molec_abbr)))
all_lev2_3_df_arranged_perlev_forGFgrapht$type_analysis <- factor(all_lev2_3_df_arranged_perlev_forGFgrapht$type_analysis, levels = rev(levels(all_lev2_3_df_arranged_perlev_forGFgrapht$type_analysis)))



mygraph1t <- ggplot(all_lev2_3_df_arranged_perlev_forGFgrapht, aes(y = molec_abbr, x = type_analysis)) +
  geom_point(aes(col = AnnotationLevel), shape = 15, size = 2) + 
  scale_color_manual(values = colours_level_annotation) +
  ylab("Annotated metabolites") +
  xlab("Analyses type") +
  scale_x_discrete(position = "top", labels = names_named_vector) +
  scale_y_discrete(position = "left") +
  coord_cartesian(xlim = c(1, length(levels(all_lev2_3_df_arranged_perlev_forGFgrapht$type_analysis))), clip = "off") +
  theme_light() +
  theme(axis.text.x = element_text(angle = 90, hjust = 0),
        #legend.position = "none",
        axis.title.y = element_text(margin = margin(t = 0, r = 25, b = 0, l = 0)))

ggsave(filename = "my_GFgrapht.png", plot = mygraph1t,  width = 6.5, height = 17, units = "in")




# my graph to make understandable the statistically significance across groups

all_lev2_3_df_arranged_perlev_forGFgraph2 <- filter(all_lev2_3_df_arranged_perlev_forGFgrapht, !duplicated(molec_abbr))

all_lev2_3_df_arranged_perlev_forGFgraph2 <- all_lev2_3_df_arranged_perlev_forGFgraph2 %>%
  mutate(N2_vs_others = factor(rep(NA, length(pull(all_lev2_3_df_arranged_perlev_forGFgraph2, colnames(all_lev2_3_df_arranged_perlev_forGFgraph2)[1]))),
                               levels = c("N2 vs VC40", "N2 vs VC1668", "N2 vs UA57", "N2 vs BR5270"))) %>%
  mutate(higher_lower = factor(rep(NA, length(pull(all_lev2_3_df_arranged_perlev_forGFgraph2, colnames(all_lev2_3_df_arranged_perlev_forGFgraph2)[1]))),
                               levels = c("higer", "none", "lower")))

### a check: MUST be FALSE:
any(duplicated(all_lev2_3_df_arranged_perlev_forGFgraph2$molec_abbr))


for (a in all_lev2_3_df_arranged_perlev_forGFgraph2$molec_abbr) {
  
  df_fill <- filter(all_lev2_3_df_arranged_perlev_forGFgraph2, molec_abbr == a)
  
  df_between <- df_fill
  df_between[2,] <- df_fill[1,]
  df_between[3,] <- df_fill[1,]
  df_between[4,] <- df_fill[1,]
  
  df_between[1, "N2_vs_others"] <- "N2 vs VC40"
  df_between[2, "N2_vs_others"] <- "N2 vs VC1668"
  df_between[3, "N2_vs_others"] <- "N2 vs UA57"
  df_between[4, "N2_vs_others"] <- "N2 vs BR5270"
  
  df_between[1, "higher_lower"] <- ifelse(grepl("N2 > VC40", df_between$N2_vs_VC40_higlow[1]), "higer",
                                          ifelse(grepl("VC40 > N2", df_between$N2_vs_VC40_higlow[1]), "lower", "none"))
  df_between[2, "higher_lower"] <- ifelse(grepl("N2 > VC1668", df_between$N2_vs_VC1668_higlow[1]), "higer",
                                          ifelse(grepl("VC1668 > N2", df_between$N2_vs_VC1668_higlow[1]), "lower", "none"))
  df_between[3, "higher_lower"] <- ifelse(grepl("N2 > UA57", df_between$N2_vs_UA57_higlow[1]), "higer",
                                          ifelse(grepl("UA57 > N2", df_between$N2_vs_UA57_higlow[1]), "lower", "none"))
  df_between[4, "higher_lower"] <- ifelse(grepl("N2 > BR5270", df_between$N2_vs_BR5270_higlow[1]), "higer",
                                          ifelse(grepl("BR5270 > N2", df_between$N2_vs_BR5270_higlow[1]), "lower", "none"))
  
  if ((which(all_lev2_3_df_arranged_perlev_forGFgraph2$molec_abbr == a) - 1)  ==0) {
    indexup <- 0
  } else {
    indexup <- 1:(which(all_lev2_3_df_arranged_perlev_forGFgraph2$molec_abbr == a) - 1)
  }
  
  
  
  if ((which(all_lev2_3_df_arranged_perlev_forGFgraph2$molec_abbr == a) + 1)  == (length(all_lev2_3_df_arranged_perlev_forGFgraph2$molec_abbr)+1)) {
    indexdown <- 0
  } else {
    indexdown <- (which(all_lev2_3_df_arranged_perlev_forGFgraph2$molec_abbr == a) + 1):(length(all_lev2_3_df_arranged_perlev_forGFgraph2$molec_abbr))
  }
  
  
  df_part1 <- all_lev2_3_df_arranged_perlev_forGFgraph2[indexup, ]
  
  df_part2 <- all_lev2_3_df_arranged_perlev_forGFgraph2[indexdown, ]
  
  
  all_lev2_3_df_arranged_perlev_forGFgraph2 <- rbind(df_part1,
                                                     df_between,
                                                     df_part2)
  
  
}

all_lev2_3_df_arranged_perlev_forGFgraph2 <- mutate(all_lev2_3_df_arranged_perlev_forGFgraph2,
                                                    Value_for_color = vector(mode = "double", length = length(pull(all_lev2_3_df_arranged_perlev_forGFgraph2, colnames(all_lev2_3_df_arranged_perlev_forGFgraph2)[molec_abbr]))))



for (i in 1:length(all_lev2_3_df_arranged_perlev_forGFgraph2$Value_for_color)) {
  
  if (all_lev2_3_df_arranged_perlev_forGFgraph2$N2_vs_others[i] == "N2 vs VC40") {
    if (all_lev2_3_df_arranged_perlev_forGFgraph2$higher_lower[i] == "none" & is.na(all_lev2_3_df_arranged_perlev_forGFgraph2$N2_vs_VC40_higlow[i])) {
      all_lev2_3_df_arranged_perlev_forGFgraph2$Value_for_color[i] <- 0
    } else if (all_lev2_3_df_arranged_perlev_forGFgraph2$higher_lower[i] == "higer" & all_lev2_3_df_arranged_perlev_forGFgraph2$N2_vs_VC40_higlow[i] == "N2 > VC40") {
      all_lev2_3_df_arranged_perlev_forGFgraph2$Value_for_color[i] <- -(log10(all_lev2_3_df_arranged_perlev_forGFgraph2$N2_vs_VC40_pval[i]))
    } else if (all_lev2_3_df_arranged_perlev_forGFgraph2$higher_lower[i] == "lower" & all_lev2_3_df_arranged_perlev_forGFgraph2$N2_vs_VC40_higlow[i] == "VC40 > N2") {
      all_lev2_3_df_arranged_perlev_forGFgraph2$Value_for_color[i] <- log10(all_lev2_3_df_arranged_perlev_forGFgraph2$N2_vs_VC40_pval[i])
    } else {
      stop("something wrong!")
    }
  } else if (all_lev2_3_df_arranged_perlev_forGFgraph2$N2_vs_others[i] == "N2 vs VC1668") {
    if (all_lev2_3_df_arranged_perlev_forGFgraph2$higher_lower[i] == "none" & is.na(all_lev2_3_df_arranged_perlev_forGFgraph2$N2_vs_VC1668_higlow[i])) {
      all_lev2_3_df_arranged_perlev_forGFgraph2$Value_for_color[i] <- 0
    } else if (all_lev2_3_df_arranged_perlev_forGFgraph2$higher_lower[i] == "higer" & all_lev2_3_df_arranged_perlev_forGFgraph2$N2_vs_VC1668_higlow[i] == "N2 > VC1668") {
      all_lev2_3_df_arranged_perlev_forGFgraph2$Value_for_color[i] <- -(log10(all_lev2_3_df_arranged_perlev_forGFgraph2$N2_vs_VC1668_pval[i]))
    } else if (all_lev2_3_df_arranged_perlev_forGFgraph2$higher_lower[i] == "lower" & all_lev2_3_df_arranged_perlev_forGFgraph2$N2_vs_VC1668_higlow[i] == "VC1668 > N2") {
      all_lev2_3_df_arranged_perlev_forGFgraph2$Value_for_color[i] <- log10(all_lev2_3_df_arranged_perlev_forGFgraph2$N2_vs_VC1668_pval[i])
    } else {
      stop("something wrong!")
    }
  } else if (all_lev2_3_df_arranged_perlev_forGFgraph2$N2_vs_others[i] == "N2 vs UA57") {
    if (all_lev2_3_df_arranged_perlev_forGFgraph2$higher_lower[i] == "none" & is.na(all_lev2_3_df_arranged_perlev_forGFgraph2$N2_vs_UA57_higlow[i])) {
      all_lev2_3_df_arranged_perlev_forGFgraph2$Value_for_color[i] <- 0
    } else if (all_lev2_3_df_arranged_perlev_forGFgraph2$higher_lower[i] == "higer" & all_lev2_3_df_arranged_perlev_forGFgraph2$N2_vs_UA57_higlow[i] == "N2 > UA57") {
      all_lev2_3_df_arranged_perlev_forGFgraph2$Value_for_color[i] <- -(log10(all_lev2_3_df_arranged_perlev_forGFgraph2$N2_vs_UA57_pval[i]))
    } else if (all_lev2_3_df_arranged_perlev_forGFgraph2$higher_lower[i] == "lower" & all_lev2_3_df_arranged_perlev_forGFgraph2$N2_vs_UA57_higlow[i] == "UA57 > N2") {
      all_lev2_3_df_arranged_perlev_forGFgraph2$Value_for_color[i] <- log10(all_lev2_3_df_arranged_perlev_forGFgraph2$N2_vs_UA57_pval[i])
    } else {
      stop("something wrong!")
    }
  } else if (all_lev2_3_df_arranged_perlev_forGFgraph2$N2_vs_others[i] == "N2 vs BR5270") {
    if (all_lev2_3_df_arranged_perlev_forGFgraph2$higher_lower[i] == "none" & is.na(all_lev2_3_df_arranged_perlev_forGFgraph2$N2_vs_BR5270_higlow[i])) {
      all_lev2_3_df_arranged_perlev_forGFgraph2$Value_for_color[i] <- 0
    } else if (all_lev2_3_df_arranged_perlev_forGFgraph2$higher_lower[i] == "higer" & all_lev2_3_df_arranged_perlev_forGFgraph2$N2_vs_BR5270_higlow[i] == "N2 > BR5270") {
      all_lev2_3_df_arranged_perlev_forGFgraph2$Value_for_color[i] <- -(log10(all_lev2_3_df_arranged_perlev_forGFgraph2$N2_vs_BR5270_pval[i]))
    } else if (all_lev2_3_df_arranged_perlev_forGFgraph2$higher_lower[i] == "lower" & all_lev2_3_df_arranged_perlev_forGFgraph2$N2_vs_BR5270_higlow[i] == "BR5270 > N2") {
      all_lev2_3_df_arranged_perlev_forGFgraph2$Value_for_color[i] <- log10(all_lev2_3_df_arranged_perlev_forGFgraph2$N2_vs_BR5270_pval[i])
    } else {
      stop("something wrong!")
    }
  } else {
    stop("something wrong!")
  } 
}




write_tsv(all_lev2_3_df_arranged_perlev_forGFgraph2, "all_lev2_3_df_arranged_perlev_forGFgraph2.txt")


colours_named_vector_higher_lower <- c(higer = "red",
                                       none = "white",
                                       lower = "blue")



## onlt higher, lower 

ggplot(all_lev2_3_df_arranged_perlev_forGFgraph2, aes(y = molec_abbr, x = N2_vs_others)) +
  geom_tile(aes(fill = higher_lower), colour = "grey50") +
  scale_fill_manual(values = colours_named_vector_higher_lower) +
  ylab("Annotated metabolites") +
  scale_x_discrete(position = "top") +
  scale_y_discrete(position = "left") +
  theme_light() +
  theme(axis.text.x = element_text(angle = 90, hjust = 0)) + 
  theme(legend.position = "none")
ggsave("my_GFgraph2discr.png", width = 4.5, height = 19, units = "in")




all_lev2_3_df_arranged_perlev_forGFgraph2_ZEROasNA <- mutate(all_lev2_3_df_arranged_perlev_forGFgraph2,
                                                             Value_for_color = ifelse(Value_for_color==0, NA, Value_for_color))

mygraph2 <- ggplot(all_lev2_3_df_arranged_perlev_forGFgraph2_ZEROasNA, aes(y = molec_abbr, x = N2_vs_others)) +
  geom_tile(aes(fill = Value_for_color), colour = "grey50") +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", na.value = "white") +
  ylab("Annotated metabolites") +
  xlab("N2 vs others") +
  scale_x_discrete(position = "top") +
  scale_y_discrete(position = "left") +
  theme_light() +
  theme(axis.text.x = element_text(angle = 90, hjust = 0),
        legend.position = "none",
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank())
ggsave(filename = "my_GFgraph2cont_ZEROasNA.png", plot = mygraph2, width = 4.5, height = 19, units = "in")


### combining the two plots:

mygraph1t_tocombine <- ggplot(all_lev2_3_df_arranged_perlev_forGFgrapht, aes(y = molec_abbr, x = type_analysis)) +
  geom_point(aes(col = AnnotationLevel), shape = 15, size = 2) + 
  scale_color_manual(values = colours_level_annotation) +
  ylab("Annotated metabolites") +
  xlab("Analyses type") +
  scale_x_discrete(position = "top", labels = names_named_vector) +
  scale_y_discrete(position = "left") +
  coord_cartesian(xlim = c(1, length(levels(all_lev2_3_df_arranged_perlev_forGFgrapht$type_analysis))), clip = "off") +
  theme_light() +
  theme(legend.position = c(-1.3, 0.85),
        axis.text.x = element_text(angle = 90, hjust = 0),
        axis.title.y = element_text(margin = margin(t = 0, r = 25, b = 0, l = 0)))

mygraph2_tocombine <- ggplot(all_lev2_3_df_arranged_perlev_forGFgraph2_ZEROasNA, aes(y = molec_abbr, x = N2_vs_others)) +
  geom_tile(aes(fill = Value_for_color), colour = "grey50") +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", na.value = "white") +
  xlab("N2 vs others") +
  scale_x_discrete(position = "top") +
  theme_light() +
  theme(axis.text.x = element_text(angle = 90, hjust = 0),
        legend.position = "none",
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank())


mygraph1t2_combined <- mygraph1t_tocombine | mygraph2_tocombine + plot_layout(widths = c(6, 2))

ggsave(filename = "my_GFgraph1t2_combined.png", plot = mygraph1t2_combined, width = 8, height = 12, units = "in")




#### same plots with shorter names in the Analyses type:

all_lev2_3_df_arranged_perlev_forGFgrapht_shrt_nms <- mutate(all_lev2_3_df_arranged_perlev_forGFgrapht,
                                                             type_analysis_shrt_nms = factor(NA, levels = c("HP s1 MD_p", "HP s1 MD_WJ", "HP s1 pR_PCL", "HP s1 pR_WJ",
                                                                                                            "HP s2 MD_p", "HP s2 MD_WJ", "HP s2 pR_PCL", "HP s2 pR_WJ",
                                                                                                            "RN s1 MD_p", "RN s1 MD_WJ", "RN s1 pR_PCL", "RN s1 pR_WJ",
                                                                                                            "RN s2 MD_p", "RN s2 MD_WJ", "RN s2 pR_PCL", "RN s2 pR_WJ")))

all_lev2_3_df_arranged_perlev_forGFgrapht_shrt_nms$type_analysis_shrt_nms[which(all_lev2_3_df_arranged_perlev_forGFgrapht_shrt_nms$type_analysis == "HILICPOS_scheme1_MSDial_PublicMSP")] <- "HP s1 MD_p"
all_lev2_3_df_arranged_perlev_forGFgrapht_shrt_nms$type_analysis_shrt_nms[which(all_lev2_3_df_arranged_perlev_forGFgrapht_shrt_nms$type_analysis == "HILICPOS_scheme1_MSDial_WormJamExpanded")] <- "HP s1 MD_WJ"
all_lev2_3_df_arranged_perlev_forGFgrapht_shrt_nms$type_analysis_shrt_nms[which(all_lev2_3_df_arranged_perlev_forGFgrapht_shrt_nms$type_analysis == "HILICPOS_scheme1_patRoon_PubChemLite")] <- "HP s1 pR_PCL"
all_lev2_3_df_arranged_perlev_forGFgrapht_shrt_nms$type_analysis_shrt_nms[which(all_lev2_3_df_arranged_perlev_forGFgrapht_shrt_nms$type_analysis == "HILICPOS_scheme1_patRoon_WormJamExpanded")] <- "HP s1 pR_WJ"
all_lev2_3_df_arranged_perlev_forGFgrapht_shrt_nms$type_analysis_shrt_nms[which(all_lev2_3_df_arranged_perlev_forGFgrapht_shrt_nms$type_analysis == "HILICPOS_scheme2_MSDial_PublicMSP")] <- "HP s2 MD_p"
all_lev2_3_df_arranged_perlev_forGFgrapht_shrt_nms$type_analysis_shrt_nms[which(all_lev2_3_df_arranged_perlev_forGFgrapht_shrt_nms$type_analysis == "HILICPOS_scheme2_MSDial_WormJamExpanded")] <- "HP s2 MD_WJ"
all_lev2_3_df_arranged_perlev_forGFgrapht_shrt_nms$type_analysis_shrt_nms[which(all_lev2_3_df_arranged_perlev_forGFgrapht_shrt_nms$type_analysis == "HILICPOS_scheme2_patRoon_PubChemLite")] <- "HP s2 pR_PCL"
all_lev2_3_df_arranged_perlev_forGFgrapht_shrt_nms$type_analysis_shrt_nms[which(all_lev2_3_df_arranged_perlev_forGFgrapht_shrt_nms$type_analysis == "HILICPOS_scheme2_patRoon_WormJamExpanded")] <- "HP s2 pR_WJ"
all_lev2_3_df_arranged_perlev_forGFgrapht_shrt_nms$type_analysis_shrt_nms[which(all_lev2_3_df_arranged_perlev_forGFgrapht_shrt_nms$type_analysis == "RPLCNEG_scheme1_MSDial_PublicMSP")] <- "RN s1 MD_p"
all_lev2_3_df_arranged_perlev_forGFgrapht_shrt_nms$type_analysis_shrt_nms[which(all_lev2_3_df_arranged_perlev_forGFgrapht_shrt_nms$type_analysis == "RPLCNEG_scheme1_MSDial_WormJamExpanded")] <- "RN s1 MD_WJ"
all_lev2_3_df_arranged_perlev_forGFgrapht_shrt_nms$type_analysis_shrt_nms[which(all_lev2_3_df_arranged_perlev_forGFgrapht_shrt_nms$type_analysis == "RPLCNEG_scheme1_patRoon_PubChemLite")] <- "RN s1 pR_PCL"
all_lev2_3_df_arranged_perlev_forGFgrapht_shrt_nms$type_analysis_shrt_nms[which(all_lev2_3_df_arranged_perlev_forGFgrapht_shrt_nms$type_analysis == "RPLCNEG_scheme1_patRoon_WormJamExpanded")] <- "RN s1 pR_WJ"
all_lev2_3_df_arranged_perlev_forGFgrapht_shrt_nms$type_analysis_shrt_nms[which(all_lev2_3_df_arranged_perlev_forGFgrapht_shrt_nms$type_analysis == "RPLCNEG_scheme2_MSDial_PublicMSP")] <- "RN s2 MD_p"
all_lev2_3_df_arranged_perlev_forGFgrapht_shrt_nms$type_analysis_shrt_nms[which(all_lev2_3_df_arranged_perlev_forGFgrapht_shrt_nms$type_analysis == "RPLCNEG_scheme2_MSDial_WormJamExpanded")] <- "RN s2 MD_WJ"
all_lev2_3_df_arranged_perlev_forGFgrapht_shrt_nms$type_analysis_shrt_nms[which(all_lev2_3_df_arranged_perlev_forGFgrapht_shrt_nms$type_analysis == "RPLCNEG_scheme2_patRoon_PubChemLite")] <- "RN s2 pR_PCL"
all_lev2_3_df_arranged_perlev_forGFgrapht_shrt_nms$type_analysis_shrt_nms[which(all_lev2_3_df_arranged_perlev_forGFgrapht_shrt_nms$type_analysis == "RPLCNEG_scheme2_patRoon_WormJamExpanded")] <- "RN s2 pR_WJ"



mygraph1t_tocombine_shrt_nms <- ggplot(all_lev2_3_df_arranged_perlev_forGFgrapht_shrt_nms, aes(y = molec_abbr, x = type_analysis_shrt_nms)) +
  geom_point(aes(col = AnnotationLevel), shape = 15, size = 2) + 
  scale_color_manual(values = colours_level_annotation) +
  ylab("Annotated metabolites") +
  xlab("Analyses type") +
  scale_x_discrete(position = "top") +
  scale_y_discrete(position = "left") +
  coord_cartesian(xlim = c(1, length(levels(all_lev2_3_df_arranged_perlev_forGFgrapht_shrt_nms$type_analysis_shrt_nms))), clip = "off") +
  theme_light() +
  theme(legend.position = c(-1.3, 0.85),
        axis.text.x = element_text(angle = 90, hjust = 0),
        axis.title.y = element_text(margin = margin(t = 0, r = 25, b = 0, l = 0)))
  

mygraph2_tocombine_shrt_nms <- ggplot(all_lev2_3_df_arranged_perlev_forGFgraph2_ZEROasNA, aes(y = molec_abbr, x = N2_vs_others)) +
  geom_tile(aes(fill = Value_for_color), colour = "grey50") +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", na.value = "white") +
  xlab("N2 vs others") +
  scale_x_discrete(position = "top") +
  theme_light() +
  theme(axis.text.x = element_text(angle = 90, hjust = 0),
        legend.position = "none",
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank())


mygraph1t2_combined_shrt_nms <- mygraph1t_tocombine_shrt_nms | mygraph2_tocombine_shrt_nms + plot_layout(widths = c(6, 2))

ggsave(filename = "my_GFgraph1t2_combined_shrt_nms.png", plot = mygraph1t2_combined_shrt_nms, width = 8, height = 12, units = "in", dpi = 300)









#####
# only 2a:
all_lev2_3_df_arranged_perlev_forGFgrapht_only2a <- all_lev2_3_df_arranged_perlev_forGFgrapht %>%
  filter(AnnotationLevel == "2a") %>%
  mutate(AnnotationLevel = droplevels(AnnotationLevel),
         molec_abbr = droplevels(molec_abbr),
         type_analysis = droplevels(type_analysis))


mygraph1t_tocombine_only2a <- ggplot(all_lev2_3_df_arranged_perlev_forGFgrapht_only2a, aes(y = molec_abbr, x = type_analysis)) +
  geom_point(aes(col = type_analysis), shape = 15, size = 2) + 
  scale_color_manual(values = colours_named_vector) +
  ylab("Annotated metabolites") +
  xlab("Analyses type") +
  scale_x_discrete(position = "top", labels = names_named_vector) +
  scale_y_discrete(position = "left") +
  coord_cartesian(xlim = c(1, length(levels(all_lev2_3_df_arranged_perlev_forGFgrapht_only2a$type_analysis))), clip = "off") +
  theme_light() +
  theme(axis.text.x = element_text(angle = 90, hjust = 0),
        legend.position =  "none",
        axis.title.y = element_text(margin = margin(t = 0, r = 25, b = 0, l = 0)))

all_lev2_3_df_arranged_perlev_forGFgraph2_ZEROasNA_only2a <- all_lev2_3_df_arranged_perlev_forGFgraph2_ZEROasNA %>%
  filter(AnnotationLevel == "2a") %>%
  mutate(AnnotationLevel = droplevels(AnnotationLevel))

mygraph2_tocombine_only2a <- ggplot(all_lev2_3_df_arranged_perlev_forGFgraph2_ZEROasNA_only2a, aes(y = molec_abbr, x = N2_vs_others)) +
  geom_tile(aes(fill = Value_for_color), colour = "grey50") +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", na.value = "white") +
  xlab("N2 vs others") +
  scale_x_discrete(position = "top") +
  theme_light() +
  theme(axis.text.x = element_text(angle = 90, hjust = 0),
        legend.position = "none",
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank())


mygraph1t2_combined_only2a <- mygraph1t_tocombine_only2a | mygraph2_tocombine_only2a + plot_layout(widths = c(6, 2))

ggsave(filename = "my_GFgraph1t2_combined_only2a.png", plot = mygraph1t2_combined_only2a, width = 7, height = 8, units = "in")





####

####box plots!!!!

all_pictures_list <- vector("list", length(all_lev2_3_df_one_compound$Name_to_use))


for (i in 1:length(all_lev2_3_df_one_compound$feature)) {
  
  this_feat_name <- all_lev2_3_df_one_compound$feature[i]
  
  this_identifier <- all_lev2_3_df_one_compound$identifier[i]
  
  this_InChIKey <- all_lev2_3_df_one_compound$InChIKey[i]
  
  this_Title <- all_lev2_3_df_one_compound$Title[i]
  
  this_compoundName <- all_lev2_3_df_one_compound$compoundName[i]
  
  this_AnnotationLevel <- as.character(all_lev2_3_df_one_compound$AnnotationLevel[i])
  
  this_ANOVA_FDR_Pvalue <- all_lev2_3_df_one_compound$ANOVA_p_value_FDR[i]
  
  this_type_analysis <- all_lev2_3_df_one_compound$type_analysis[i]
  
  this_molec <- as.character(all_lev2_3_df_arranged_perlev_forGFgraph$molec[which(all_lev2_3_df_arranged_perlev_forGFgraph$identifier == this_identifier)[1]])
  
  this_molec_abbr <- as.character(all_lev2_3_df_arranged_perlev_forGFgraph$molec_abbr[which(all_lev2_3_df_arranged_perlev_forGFgraph$identifier == this_identifier)[1]])
  
  this_Name_to_use <- all_lev2_3_df_one_compound$Name_to_use[i]
  
  
  data_to_use_for_this_boxplot <- tibble()
  
  if (this_type_analysis == "HILICPOS_scheme1_MSDial_PublicMSP") {data_to_use_for_this_boxplot <- MSDIAL_HLP4_scheme1_dry_transf }
  if (this_type_analysis == "HILICPOS_scheme2_MSDial_PublicMSP") {data_to_use_for_this_boxplot <- MSDIAL_HLP4_scheme2_dry_transf }
  
  if (this_type_analysis == "HILICPOS_scheme1_MSDial_WormJamExpanded") {data_to_use_for_this_boxplot <- MSDIAL_WJ_HLP4_scheme1_dry_transf }
  if (this_type_analysis == "HILICPOS_scheme2_MSDial_WormJamExpanded") {data_to_use_for_this_boxplot <- MSDIAL_WJ_HLP4_scheme2_dry_transf }
  
  if (this_type_analysis == "HILICPOS_scheme1_patRoon_PubChemLite" | this_type_analysis == "HILICPOS_scheme1_patRoon_WormJamExpanded") {data_to_use_for_this_boxplot <- patRoon_HLP4_shem1_dry_transf }
  if (this_type_analysis == "HILICPOS_scheme2_patRoon_PubChemLite" | this_type_analysis == "HILICPOS_scheme2_patRoon_WormJamExpanded") {data_to_use_for_this_boxplot <- patRoon_HLP4_shem2_dry_transf }
  
  if (this_type_analysis == "RPLCNEG_scheme1_MSDial_PublicMSP") {data_to_use_for_this_boxplot <- MSDIAL_RPN4_scheme1_dry_transf }
  if (this_type_analysis == "RPLCNEG_scheme2_MSDial_PublicMSP") {data_to_use_for_this_boxplot <- MSDIAL_RPN4_scheme2_dry_transf }
  
  if (this_type_analysis == "RPLCNEG_scheme1_MSDial_WormJamExpanded") {data_to_use_for_this_boxplot <- MSDIAL_WJ_RPN4_scheme1_dry_transf }
  if (this_type_analysis == "RPLCNEG_scheme2_MSDial_WormJamExpanded") {data_to_use_for_this_boxplot <- MSDIAL_WJ_RPN4_scheme2_dry_transf }
  
  if (this_type_analysis == "RPLCNEG_scheme1_patRoon_PubChemLite" | this_type_analysis == "RPLCNEG_scheme1_patRoon_WormJamExpanded") {data_to_use_for_this_boxplot <- patRoon_RPN4_scheme1_dry_transf }
  if (this_type_analysis == "RPLCNEG_scheme2_patRoon_PubChemLite" | this_type_analysis == "RPLCNEG_scheme2_patRoon_WormJamExpanded") {data_to_use_for_this_boxplot <- patRoon_RPN4_scheme2_dry_transf }
  
  if (identical(data_to_use_for_this_boxplot, tibble())) {stop("tibble is empty..!")}
  
  if (!this_feat_name %in% colnames(data_to_use_for_this_boxplot)) {stop(paste0(this_feat_name, " not present in the data table ", this_type_analysis, "\n"))} 
  
  
  data_to_use_for_this_boxplot_fill <- select(data_to_use_for_this_boxplot,
                                              all_of(c("ANALYSIS", "Sample_ID", this_feat_name)))
  colnames(data_to_use_for_this_boxplot_fill)[which(colnames(data_to_use_for_this_boxplot_fill)=="Sample_ID")] <- "SampleGroup"
  
  
  this_box_plot <- ggplot(data_to_use_for_this_boxplot_fill, aes(x = SampleGroup, y = !!sym(this_feat_name), fill = SampleGroup)) +
    geom_boxplot() + 
    ylab(this_molec_abbr) +
    theme_classic() +
    theme(panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank(),
          axis.text.y = element_blank(),
          axis.title.y = element_text(size = 18),
          axis.ticks.y = element_blank(),
          axis.text.x = element_text(size = 16),
          axis.title.x = element_blank(),
          legend.position = "none")
  
  This_complementary_Table_first <- tibble(Characteristic = c("Pubchem CID:",
                                                              "InChIKey:",
                                                              "Title:",
                                                              "Compound name:",
                                                              "Name used:",
                                                              "Annotation level:",
                                                              "ANOVA FDR PValue:"),
                                           Value = c(this_identifier,
                                                     this_InChIKey,
                                                     this_Title,
                                                     this_compoundName,
                                                     this_Name_to_use,
                                                     this_AnnotationLevel,
                                                     ifelse(this_ANOVA_FDR_Pvalue<0.001, sprintf("%.3e", this_ANOVA_FDR_Pvalue), sprintf("%.4f", this_ANOVA_FDR_Pvalue))))
  
  
  
  LSD_info <- tibble(Characteristic = character(),
                     Value = character())
  
  
  vector_with_all_cases <- c("N2_vs_VC40_higlow", "N2_vs_VC1668_higlow", "N2_vs_UA57_higlow", "N2_vs_BR5270_higlow", "VC40_vs_VC1668_higlow", "VC40_vs_UA57_higlow", "VC40_vs_BR5270", "VC1668_vs_UA57_higlow", "VC1668_vs_BR5270_higlow", "UA57_vs_BR5270_higlow")
  vector_with_sign_pairwise_diff <- vector_with_all_cases
  
  for (e in vector_with_all_cases) {
    if (is.na(pull(all_lev2_3_df_one_compound, e)[i])) {
      vector_with_sign_pairwise_diff <- vector_with_sign_pairwise_diff[-which(vector_with_sign_pairwise_diff==e)]
    }
  }
  
  if (length(vector_with_sign_pairwise_diff) == 0) {
    LSD_info <- tibble(Characteristic = "Fisher's LSD pairwise comparison:",
                       Value = "No pairwise significant differences")
  } else {
    
    this_complete_text_full <- character()
    
    for (v in vector_with_sign_pairwise_diff) {
      
      Pval_col_name <- ifelse(v == "N2_vs_VC40_higlow", "N2_vs_VC40_pval",
                              ifelse(v == "N2_vs_VC1668_higlow", "N2_vs_VC1668_pval",
                                     ifelse(v == "N2_vs_UA57_higlow", "N2_vs_UA57_pval",
                                            ifelse(v == "N2_vs_BR5270_higlow", "N2_vs_BR5270_pval",
                                                   ifelse(v == "VC40_vs_VC1668_higlow", "VC40_vs_VC1668_pval",
                                                          ifelse(v == "VC40_vs_UA57_higlow", "VC40_vs_UA57_pval",
                                                                 ifelse(v == "VC40_vs_BR5270", "VC40_vs_BR5270_pval",
                                                                        ifelse(v == "VC1668_vs_UA57_higlow", "VC1668_vs_UA57_pval",
                                                                               ifelse(v == "VC1668_vs_BR5270_higlow", "VC1668_vs_BR5270_pval",
                                                                                      ifelse(v == "UA57_vs_BR5270_higlow", "UA57_vs_BR5270_pval"))))))))))
      
      This_Pvalue <- pull(all_lev2_3_df_one_compound, Pval_col_name)[i]
      
      this_complete_text_single <- character()
      
      if (which(vector_with_sign_pairwise_diff == v) == 1) {
        this_complete_text_single <- paste0(pull(all_lev2_3_df_one_compound, v)[i], " - Pvalue: ", ifelse(This_Pvalue<0.001, sprintf("%.3e", This_Pvalue), sprintf("%.4f", This_Pvalue)))
      } else {
        this_complete_text_single <- paste0("\n", pull(all_lev2_3_df_one_compound, v)[i], " - Pvalue: ", ifelse(This_Pvalue<0.001, sprintf("%.3e", This_Pvalue), sprintf("%.4f", This_Pvalue)))
      }
      
      this_complete_text_full <- paste0(this_complete_text_full, this_complete_text_single)
      
    }
    
    LSD_info <- tibble(Characteristic = "Fisher's LSD pairwise comparison:",
                       Value = this_complete_text_full)
    
  }
  
  
  
  type_of_analyses_info <- tibble(Characteristic = "Type of analysis:",
                                  Value = str_replace_all(this_type_analysis, "_", " "))
  
  all_lev2_3_df_arranged_perlev_fill <- filter(all_lev2_3_df_arranged_perlev,
                                               Name_to_use == this_Name_to_use,
                                               type_analysis != this_type_analysis)
  all_lev2_3_df_arranged_perlev_fill <- filter(all_lev2_3_df_arranged_perlev_fill,
                                               !duplicated(type_analysis))
  
  if (length(all_lev2_3_df_arranged_perlev_fill$Name_to_use) > 0) {
    
    this_other_analysis_full <- character()
    
    for (u in 1:length(all_lev2_3_df_arranged_perlev_fill$Name_to_use)) {
      if (u == 1) {
        this_other_analysis_single <- str_replace_all(all_lev2_3_df_arranged_perlev_fill$type_analysis[u], "_", " ")
      } else {
        this_other_analysis_single <- paste0("\n", str_replace_all(all_lev2_3_df_arranged_perlev_fill$type_analysis[u], "_", " "))
      }
      
      this_other_analysis_full <- paste0(this_other_analysis_full, this_other_analysis_single)
      
    }
    
    type_of_analyses_info <- rbind(type_of_analyses_info,  tibble(Characteristic = "Other analyses in which significative:",
                                                                  Value = this_other_analysis_full))
    
  }
  
  This_complementary_Table <- rbind(This_complementary_Table_first, LSD_info, type_of_analyses_info)
  
  this_grid_table <- tableGrob(This_complementary_Table,
                               rows = NULL,
                               cols = NULL,
                               theme = ttheme_minimal(core=list(fg_params=list(hjust=0, x=0.1, y = 1, vjust = 1))))
  
  this_grid_table$widths  <- unit(rep(1/ncol(this_grid_table), ncol(this_grid_table)), "npc")
  
  
  
  this_combined_picture <- grid.arrange(this_box_plot, this_grid_table, ncol = 2)
  
  the_order_number <- ifelse(i<10,
                             paste0("00", as.character(i)),
                             ifelse(i<100,
                                    paste0("0", as.character(i)),
                                    as.character(i)))
  
  
  
  ggsave(filename = paste0("mys", the_order_number, "_", this_identifier, ".png"),
         plot = this_combined_picture,
         device = "png",
         width = 35,
         height = 14,
         units = "cm",
         dpi = 300)
  
  
  
  all_pictures_list[[i]] <- this_combined_picture
  
}





#### other statistics for the paper

all_lev2_3_df_arranged_perlev_analyses_type_stat_table <- mutate(all_lev2_3_df_arranged_perlev,
                                                                 this_row_number = 1:length(all_lev2_3_df_arranged_perlev$Name_to_use))

rows_to_remove <- numeric()


for (a in unique(all_lev2_3_df_arranged_perlev_analyses_type_stat_table$Name_to_use)) {
  df_fill <- filter(all_lev2_3_df_arranged_perlev_analyses_type_stat_table, Name_to_use == a)
  this_rows_to_remove <- df_fill$this_row_number[duplicated(df_fill$type_analysis)]
  rows_to_remove <- c(rows_to_remove,this_rows_to_remove)
}


all_lev2_3_df_arranged_perlev_analyses_type_stat_table_fill <- all_lev2_3_df_arranged_perlev_analyses_type_stat_table[-rows_to_remove,]


all_lev2_3_df_arranged_perlev_analyses_type_stat <- all_lev2_3_df_arranged_perlev_analyses_type_stat_table_fill %>%
  group_by(type_analysis) %>%
  summarize(N = n())

##MS-DIAL
Name_to_use_with_MSDIAL <- numeric()

for (a in unique(all_lev2_3_df_arranged_perlev_analyses_type_stat_table_fill$Name_to_use)) {
  df_fill <- filter(all_lev2_3_df_arranged_perlev_analyses_type_stat_table_fill, Name_to_use == a)
  
  if (any(grepl("MSDial", df_fill$type_analysis))) {
    Name_to_use_with_MSDIAL <- c(Name_to_use_with_MSDIAL, a)
  }
}

##patROON
Name_to_use_with_patRoon <- numeric()

for (a in unique(all_lev2_3_df_arranged_perlev_analyses_type_stat_table_fill$Name_to_use)) {
  df_fill <- filter(all_lev2_3_df_arranged_perlev_analyses_type_stat_table_fill, Name_to_use == a)
  
  if (any(grepl("patRoon", df_fill$type_analysis))) {
    Name_to_use_with_patRoon <- c(Name_to_use_with_patRoon, a)
  }
}


##HILIC POS

Name_to_use_with_HILICPOS <- numeric()

for (a in unique(all_lev2_3_df_arranged_perlev_analyses_type_stat_table_fill$Name_to_use)) {
  df_fill <- filter(all_lev2_3_df_arranged_perlev_analyses_type_stat_table_fill, Name_to_use == a)
  
  if (any(grepl("HILICPOS", df_fill$type_analysis))) {
    Name_to_use_with_HILICPOS <- c(Name_to_use_with_HILICPOS, a)
  }
}

##RPLC NEG

Name_to_use_with_RPLCNEG <- numeric()

for (a in unique(all_lev2_3_df_arranged_perlev_analyses_type_stat_table_fill$Name_to_use)) {
  df_fill <- filter(all_lev2_3_df_arranged_perlev_analyses_type_stat_table_fill, Name_to_use == a)
  
  if (any(grepl("RPLCNEG", df_fill$type_analysis))) {
    Name_to_use_with_RPLCNEG <- c(Name_to_use_with_RPLCNEG, a)
  }
}

##scheme1

Name_to_use_with_scheme1 <- numeric()

for (a in unique(all_lev2_3_df_arranged_perlev_analyses_type_stat_table_fill$Name_to_use)) {
  df_fill <- filter(all_lev2_3_df_arranged_perlev_analyses_type_stat_table_fill, Name_to_use == a)
  
  if (any(grepl("scheme1", df_fill$type_analysis))) {
    Name_to_use_with_scheme1 <- c(Name_to_use_with_scheme1, a)
  }
}

##scheme2

Name_to_use_with_scheme2 <- numeric()

for (a in unique(all_lev2_3_df_arranged_perlev_analyses_type_stat_table_fill$Name_to_use)) {
  df_fill <- filter(all_lev2_3_df_arranged_perlev_analyses_type_stat_table_fill, Name_to_use == a)
  
  if (any(grepl("scheme2", df_fill$type_analysis))) {
    Name_to_use_with_scheme2 <- c(Name_to_use_with_scheme2, a)
  }
}




length(Name_to_use_with_MSDIAL)
length(Name_to_use_with_patRoon)
length(Name_to_use_with_HILICPOS)
length(Name_to_use_with_RPLCNEG)
length(Name_to_use_with_scheme1)
length(Name_to_use_with_scheme2)






## N2 vs VC40

all_lev2_3_df_arranged_perlev_forGFgraph2_ZEROasNA_N2_VC40_summ <- all_lev2_3_df_arranged_perlev_forGFgraph2_ZEROasNA %>%
  filter(N2_vs_others == "N2 vs VC40") %>%
  group_by(higher_lower) %>%
  summarize(N= n())

## N2 vs VC1668

all_lev2_3_df_arranged_perlev_forGFgraph2_ZEROasNA_N2_VC1668_summ <- all_lev2_3_df_arranged_perlev_forGFgraph2_ZEROasNA %>%
  filter(N2_vs_others == "N2 vs VC1668") %>%
  group_by(higher_lower) %>%
  summarize(N= n())

## N2 vs UA57

all_lev2_3_df_arranged_perlev_forGFgraph2_ZEROasNA_N2_UA57_summ <- all_lev2_3_df_arranged_perlev_forGFgraph2_ZEROasNA %>%
  filter(N2_vs_others == "N2 vs UA57") %>%
  group_by(higher_lower) %>%
  summarize(N= n())

## N2 vs BR5270

all_lev2_3_df_arranged_perlev_forGFgraph2_ZEROasNA_N2_BR5270_summ <- all_lev2_3_df_arranged_perlev_forGFgraph2_ZEROasNA %>%
  filter(N2_vs_others == "N2 vs BR5270") %>%
  group_by(higher_lower) %>%
  summarize(N= n())






#### preparing dataset for pathway analyses? - first try with metaboAnalyst - but actually failed as idea...


analysis_and_sampleID <- tibble(bind_rows(select(MSDIAL_HLP4_scheme1_dry_transf, 1:2),
                                          select(MSDIAL_HLP4_scheme2_dry_transf, 1:2),
                                          select(MSDIAL_WJ_HLP4_scheme1_dry_transf, 1:2),
                                          select(MSDIAL_WJ_HLP4_scheme2_dry_transf, 1:2),
                                          select(patRoon_HLP4_shem1_dry_transf, 1:2),
                                          select(patRoon_HLP4_shem2_dry_transf, 1:2),
                                          select(MSDIAL_RPN4_scheme1_dry_transf, 1:2),
                                          select(MSDIAL_RPN4_scheme2_dry_transf, 1:2),
                                          select(MSDIAL_WJ_RPN4_scheme1_dry_transf, 1:2),
                                          select(MSDIAL_WJ_RPN4_scheme2_dry_transf, 1:2),
                                          select(patRoon_RPN4_scheme1_dry_transf, 1:2),
                                          select(patRoon_RPN4_scheme2_dry_transf, 1:2))) %>%
  filter(!duplicated(ANALYSIS)) %>%
  arrange(ANALYSIS)   %>%  arrange(Sample_ID)




all_lev2_3_df_one_compound_w_data <- all_lev2_3_df_one_compound


all_lev2_3_df_one_compound_w_data[, analysis_and_sampleID$ANALYSIS] <- NA



for (i in 1:length(pull(all_lev2_3_df_one_compound_w_data, 1))) {
  
  this_type_analysis <- all_lev2_3_df_one_compound_w_data$type_analysis[i]
  
  
  if (this_type_analysis == "HILICPOS_scheme1_MSDial_PublicMSP") {this_data_to_use <- MSDIAL_HLP4_scheme1_dry_transf }
  if (this_type_analysis == "HILICPOS_scheme2_MSDial_PublicMSP") {this_data_to_use <- MSDIAL_HLP4_scheme2_dry_transf }
  
  if (this_type_analysis == "HILICPOS_scheme1_MSDial_WormJamExpanded") {this_data_to_use <- MSDIAL_WJ_HLP4_scheme1_dry_transf }
  if (this_type_analysis == "HILICPOS_scheme2_MSDial_WormJamExpanded") {this_data_to_use <- MSDIAL_WJ_HLP4_scheme2_dry_transf }
  
  if (this_type_analysis == "HILICPOS_scheme1_patRoon_PubChemLite" | this_type_analysis == "HILICPOS_scheme1_patRoon_WormJamExpanded") {this_data_to_use <- patRoon_HLP4_shem1_dry_transf }
  if (this_type_analysis == "HILICPOS_scheme2_patRoon_PubChemLite" | this_type_analysis == "HILICPOS_scheme2_patRoon_WormJamExpanded") {this_data_to_use <- patRoon_HLP4_shem2_dry_transf }
  
  if (this_type_analysis == "RPLCNEG_scheme1_MSDial_PublicMSP") {this_data_to_use <- MSDIAL_RPN4_scheme1_dry_transf }
  if (this_type_analysis == "RPLCNEG_scheme2_MSDial_PublicMSP") {this_data_to_use <- MSDIAL_RPN4_scheme2_dry_transf }
  
  if (this_type_analysis == "RPLCNEG_scheme1_MSDial_WormJamExpanded") {this_data_to_use <- MSDIAL_WJ_RPN4_scheme1_dry_transf }
  if (this_type_analysis == "RPLCNEG_scheme2_MSDial_WormJamExpanded") {this_data_to_use <- MSDIAL_WJ_RPN4_scheme2_dry_transf }
  
  if (this_type_analysis == "RPLCNEG_scheme1_patRoon_PubChemLite" | this_type_analysis == "RPLCNEG_scheme1_patRoon_WormJamExpanded") {this_data_to_use <- patRoon_RPN4_scheme1_dry_transf }
  if (this_type_analysis == "RPLCNEG_scheme2_patRoon_PubChemLite" | this_type_analysis == "RPLCNEG_scheme2_patRoon_WormJamExpanded") {this_data_to_use <- patRoon_RPN4_scheme2_dry_transf }
  
  
  
  
  
  dataset_fil_vector <- pull(this_data_to_use, all_lev2_3_df_one_compound_w_data$feature[i])
  names(dataset_fil_vector) <- pull(this_data_to_use, 1)
  
  
  
  for (sa in names(dataset_fil_vector)) {
    all_lev2_3_df_one_compound_w_data[i, sa] <- dataset_fil_vector[sa]
  }
  
  
}



##### Trying the FELLA pacakge
### dataset preparation:

all_lev2_3_df_one_compound_w_data_more_codes <- add_column(all_lev2_3_df_one_compound_w_data,
                                                           HMDBcode = as.character(rep(NA, length(pull(all_lev2_3_df_one_compound_w_data, 1)))),
                                                           KEGGcode = as.character(rep(NA, length(pull(all_lev2_3_df_one_compound_w_data, 1)))),
                                                           .after = "compoundName")

for (i in 1:length(pull(all_lev2_3_df_one_compound_w_data_more_codes, 1))) {
  
  if (pull(all_lev2_3_df_one_compound_w_data_more_codes, "identifier")[i] %in% metabolitesMapping$CID) {
    
    all_lev2_3_df_one_compound_w_data_more_codes[i, "HMDBcode"] <- pull(filter(metabolitesMapping, CID == pull(all_lev2_3_df_one_compound_w_data, "identifier")[i]), "HMDB")[1]
    all_lev2_3_df_one_compound_w_data_more_codes[i, "KEGGcode"] <- pull(filter(metabolitesMapping, CID == pull(all_lev2_3_df_one_compound_w_data, "identifier")[i]), "KEGG")[1]
    
  } else {
    
    print(paste0("CID ", pull(all_lev2_3_df_one_compound_w_data_more_codes, "identifier")[i], ", ", pull(all_lev2_3_df_one_compound_w_data, "compoundName")[i]  ,", is not present in the metabolitesMapping library"))
    
  }
}

# 

## manually adding a few missing:


all_lev2_3_df_one_compound_w_data_more_codes[which(all_lev2_3_df_one_compound_w_data_more_codes$identifier=="92136"), "KEGGcode"] <- "C00956"
all_lev2_3_df_one_compound_w_data_more_codes[which(all_lev2_3_df_one_compound_w_data_more_codes$identifier=="1000"), "KEGGcode"] <- "C02735"
all_lev2_3_df_one_compound_w_data_more_codes[which(all_lev2_3_df_one_compound_w_data_more_codes$identifier=="75619"), "KEGGcode"] <- "C02997"
all_lev2_3_df_one_compound_w_data_more_codes[which(all_lev2_3_df_one_compound_w_data_more_codes$identifier=="6267"), "KEGGcode"] <- "C00152"
all_lev2_3_df_one_compound_w_data_more_codes[which(all_lev2_3_df_one_compound_w_data_more_codes$identifier=="135398635"), "KEGGcode"] <- "C00387"
all_lev2_3_df_one_compound_w_data_more_codes[which(all_lev2_3_df_one_compound_w_data_more_codes$identifier=="617"), "KEGGcode"] <- "C00716"
all_lev2_3_df_one_compound_w_data_more_codes[which(all_lev2_3_df_one_compound_w_data_more_codes$identifier=="135398641"), "KEGGcode"] <- "C00294"
all_lev2_3_df_one_compound_w_data_more_codes[which(all_lev2_3_df_one_compound_w_data_more_codes$identifier=="107738"), "KEGGcode"] <- "C03017"
all_lev2_3_df_one_compound_w_data_more_codes[which(all_lev2_3_df_one_compound_w_data_more_codes$identifier=="92832"), "KEGGcode"] <- "C02727"
all_lev2_3_df_one_compound_w_data_more_codes[which(all_lev2_3_df_one_compound_w_data_more_codes$identifier=="440311"), "KEGGcode"] <- "C04368"
all_lev2_3_df_one_compound_w_data_more_codes[which(all_lev2_3_df_one_compound_w_data_more_codes$identifier=="440120"), "KEGGcode"] <- "C03793"
all_lev2_3_df_one_compound_w_data_more_codes[which(all_lev2_3_df_one_compound_w_data_more_codes$identifier=="41211"), "KEGGcode"] <- "C01367"
all_lev2_3_df_one_compound_w_data_more_codes[which(all_lev2_3_df_one_compound_w_data_more_codes$identifier=="135398631"), "KEGGcode"] <- "C00144"
all_lev2_3_df_one_compound_w_data_more_codes[which(all_lev2_3_df_one_compound_w_data_more_codes$identifier=="439258"), "KEGGcode"] <- "C00542"
all_lev2_3_df_one_compound_w_data_more_codes[which(all_lev2_3_df_one_compound_w_data_more_codes$identifier=="439756"), "KEGGcode"] <- "C02571"
all_lev2_3_df_one_compound_w_data_more_codes[which(all_lev2_3_df_one_compound_w_data_more_codes$identifier=="11953814"), "KEGGcode"] <- "C02838"
all_lev2_3_df_one_compound_w_data_more_codes[which(all_lev2_3_df_one_compound_w_data_more_codes$identifier=="439224"), "KEGGcode"] <- "C00386"
all_lev2_3_df_one_compound_w_data_more_codes[which(all_lev2_3_df_one_compound_w_data_more_codes$identifier=="94244"), "KEGGcode"] <- "C11332"
all_lev2_3_df_one_compound_w_data_more_codes[which(all_lev2_3_df_one_compound_w_data_more_codes$identifier=="3208"), "KEGGcode"] <- "C16525"
all_lev2_3_df_one_compound_w_data_more_codes[which(all_lev2_3_df_one_compound_w_data_more_codes$identifier=="3210"), "KEGGcode"] <- "C16522"
all_lev2_3_df_one_compound_w_data_more_codes[which(all_lev2_3_df_one_compound_w_data_more_codes$identifier=="65085"), "KEGGcode"] <- "C08432"
all_lev2_3_df_one_compound_w_data_more_codes[which(all_lev2_3_df_one_compound_w_data_more_codes$identifier=="7018721"), "KEGGcode"] <- "C21016"
all_lev2_3_df_one_compound_w_data_more_codes[which(all_lev2_3_df_one_compound_w_data_more_codes$identifier=="6287"), "KEGGcode"] <- "C00183"
all_lev2_3_df_one_compound_w_data_more_codes[which(all_lev2_3_df_one_compound_w_data_more_codes$identifier=="6450191"), "KEGGcode"] <- "C17425"
all_lev2_3_df_one_compound_w_data_more_codes[which(all_lev2_3_df_one_compound_w_data_more_codes$identifier=="27476"), "KEGGcode"] <- "C02494"
all_lev2_3_df_one_compound_w_data_more_codes[which(all_lev2_3_df_one_compound_w_data_more_codes$identifier=="5950"), "KEGGcode"] <- "C00041"
all_lev2_3_df_one_compound_w_data_more_codes[which(all_lev2_3_df_one_compound_w_data_more_codes$identifier=="5288725"), "KEGGcode"] <- "C02721"
all_lev2_3_df_one_compound_w_data_more_codes[which(all_lev2_3_df_one_compound_w_data_more_codes$identifier=="439227"), "KEGGcode"] <- "C00408"
all_lev2_3_df_one_compound_w_data_more_codes[which(all_lev2_3_df_one_compound_w_data_more_codes$identifier=="9547068"), "KEGGcode"] <- "C21484"
all_lev2_3_df_one_compound_w_data_more_codes[which(all_lev2_3_df_one_compound_w_data_more_codes$identifier=="965"), "KEGGcode"] <- "C00712"
all_lev2_3_df_one_compound_w_data_more_codes[which(all_lev2_3_df_one_compound_w_data_more_codes$identifier=="4668"), "KEGGcode"] <- "C08362"
all_lev2_3_df_one_compound_w_data_more_codes[which(all_lev2_3_df_one_compound_w_data_more_codes$identifier=="3931"), "KEGGcode"] <- "C01595"
all_lev2_3_df_one_compound_w_data_more_codes[which(all_lev2_3_df_one_compound_w_data_more_codes$identifier=="138392130"), "KEGGcode"] <- "C08265"
all_lev2_3_df_one_compound_w_data_more_codes[which(all_lev2_3_df_one_compound_w_data_more_codes$identifier=="163841"), "KEGGcode"] <- "C16300"
all_lev2_3_df_one_compound_w_data_more_codes[which(all_lev2_3_df_one_compound_w_data_more_codes$identifier=="8216"), "KEGGcode"] <- "C08316"



write_tsv(all_lev2_3_df_one_compound_w_data_more_codes, "all_lev2_3_df_one_compound_w_data_more_codes.txt")


print(paste0("there are ", sum(!is.na(all_lev2_3_df_one_compound_w_data_more_codes$KEGGcode)), " valid KEGG codes out of ", length(all_lev2_3_df_one_compound_w_data_more_codes$KEGGcode), " compounds"))
## there are 46 valid KEGG codes out of 85 compounds !!


### actually using the FELLA package now:



tot_c_eleg <- buildGraphFromKEGGREST(organism = "cel", filter.path = NULL)

buildDataFromGraph(keggdata.graph = tot_c_eleg,
                   databaseDir = "C:/databases/FELLA/c_elegans",
                   internalDir = FALSE,
                   matrices = "diffusion",
                   normality = "diffusion")

fella.data_c_elegans <- loadKEGGdata(databaseDir = "C:/databases/FELLA/c_elegans",
                                     internalDir = FALSE,
                                     loadMatrix = "diffusion")

Significative_molecules_KEGG_codes <- all_lev2_3_df_one_compound_w_data_more_codes$KEGGcode[which(!is.na(all_lev2_3_df_one_compound_w_data_more_codes$KEGGcode))]

enrich_analysis_sign_molecules <- defineCompounds(compounds = Significative_molecules_KEGG_codes,
                                                  data = fella.data_c_elegans)

print(paste0(length(getInput(enrich_analysis_sign_molecules)), " codes used; ",
             length(getExcluded(enrich_analysis_sign_molecules)), " codes excluded; ",
             length(Significative_molecules_KEGG_codes), " total codes"))


enrich_analysis_sign_molecules_run <- runDiffusion(object = enrich_analysis_sign_molecules,
                                                   data = fella.data_c_elegans,
                                                   approx = "normality")

enrich_analysis_sign_molecules_run_resultsGraph <- generateResultsGraph(object = enrich_analysis_sign_molecules_run,
                                                                        method = "diffusion",
                                                                        data = fella.data_c_elegans,
                                                                        nlimit = 250)

exportResults(format = "csv", file = "enrich_analysis_sign_molecules_run_resultsGraph.csv",
              method = "diffusion", object = enrich_analysis_sign_molecules_run, data = fella.data_c_elegans)


generated_table <- generateResultsTable(method = "diffusion",
                                        nlimit = 1000, 
                                        #threshold = 0.05, plimit = 15, LabelLengthAtPlot = 45,
                                        #capPscores = 1e-06,
                                        object = enrich_analysis_sign_molecules_run, data = fella.data_c_elegans)

plotGraph(graph = enrich_analysis_sign_molecules_run_resultsGraph,
          #layout = FALSE, graph.layout = NULL,
          #plotLegend = TRUE, plot.fun = "plot.igraph", NamesAsLabels = TRUE,
          vertex.label.cex = 0.5)                                                                     



write_graph(enrich_analysis_sign_molecules_run_resultsGraph, "enrich_analysis_graph_to_export.png")

png("enrich_analysis_sign_molecules_run_resultsGraph.png", width = 1000, height = 1000)
plotGraph(graph = enrich_analysis_sign_molecules_run_resultsGraph,
          #layout = FALSE, graph.layout = NULL,
          #plotLegend = TRUE, plot.fun = "plot.igraph", NamesAsLabels = TRUE,
          vertex.label.cex = 1)
dev.off()





#### only those significant in some:

write_csv(all_lev2_3_df_without_N2, "onlysomesign_00_all_lev2_3_df_without_N2.csv")
write_csv(all_lev2_3_df_N2_VC40, "onlysomesign_01_all_lev2_3_df_N2_VC40.csv")
write_csv(all_lev2_3_df_VC40_N2, "onlysomesign_02_all_lev2_3_df_VC40_N2.csv")
write_csv(all_lev2_3_df_N2_VC1668, "onlysomesign_03_all_lev2_3_df_N2_VC1668.csv")
write_csv(all_lev2_3_df_VC1668_N2, "onlysomesign_04_all_lev2_3_df_VC1668_N2.csv")
write_csv(all_lev2_3_df_N2_UA57, "onlysomesign_05_all_lev2_3_df_N2_UA57.csv")
write_csv(all_lev2_3_df_UA57_N2, "onlysomesign_06_all_lev2_3_df_UA57_N2.csv")
write_csv(all_lev2_3_df_N2_BR5270, "onlysomesign_07_all_lev2_3_df_N2_BR5270.csv")
write_csv(all_lev2_3_df_BR5270_N2, "onlysomesign_08_all_lev2_3_df_BR5270_N2.csv")


### preparing tables again, with the codes to do the FELLA pathway


## without N2 significant:

all_lev2_3_df_one_compound_w_data_more_codes_without_N2 <- filter(all_lev2_3_df_one_compound_w_data_more_codes, is.na(N2_vs_VC40_higlow) & is.na(N2_vs_VC1668_higlow) & is.na(N2_vs_UA57_higlow) & is.na(N2_vs_BR5270_higlow))

## VC40

# N2 > VC40

all_lev2_3_df_one_compound_w_data_more_codes_N2_VC40 <- filter(all_lev2_3_df_one_compound_w_data_more_codes, grepl("N2 > VC40", N2_vs_VC40_higlow))

# N2 < VC40

all_lev2_3_df_one_compound_w_data_more_codes_VC40_N2 <- filter(all_lev2_3_df_one_compound_w_data_more_codes, grepl("VC40 > N2", N2_vs_VC40_higlow))

# N2 or NC40

all_lev2_3_df_one_compound_w_data_more_codes_N2orVC40 <- filter(all_lev2_3_df_one_compound_w_data_more_codes, grepl("N2 > VC40", N2_vs_VC40_higlow) | grepl("VC40 > N2", N2_vs_VC40_higlow))




## VC1668

# N2 > VC1668

all_lev2_3_df_one_compound_w_data_more_codes_N2_VC1668 <- filter(all_lev2_3_df_one_compound_w_data_more_codes, grepl("N2 > VC1668", N2_vs_VC1668_higlow))

# N2 < VC1668

all_lev2_3_df_one_compound_w_data_more_codes_VC1668_N2 <- filter(all_lev2_3_df_one_compound_w_data_more_codes, grepl("VC1668 > N2", N2_vs_VC1668_higlow))


# N2 or VC1668

all_lev2_3_df_one_compound_w_data_more_codes_N2orVC1668 <- filter(all_lev2_3_df_one_compound_w_data_more_codes, grepl("N2 > VC1668", N2_vs_VC1668_higlow) | grepl("VC1668 > N2", N2_vs_VC1668_higlow))





## UA57

# N2 > UA57

all_lev2_3_df_one_compound_w_data_more_codes_N2_UA57 <- filter(all_lev2_3_df_one_compound_w_data_more_codes, grepl("N2 > UA57", N2_vs_UA57_higlow))

# N2 < UA57

all_lev2_3_df_one_compound_w_data_more_codes_UA57_N2 <- filter(all_lev2_3_df_one_compound_w_data_more_codes, grepl("UA57 > N2", N2_vs_UA57_higlow))


# N2 or UA57

all_lev2_3_df_one_compound_w_data_more_codes_N2orUA57 <- filter(all_lev2_3_df_one_compound_w_data_more_codes, grepl("N2 > UA57", N2_vs_UA57_higlow) | grepl("UA57 > N2", N2_vs_UA57_higlow))





## BR5270

# N2 > BR5270

all_lev2_3_df_one_compound_w_data_more_codes_N2_BR5270 <- filter(all_lev2_3_df_one_compound_w_data_more_codes, grepl("N2 > BR5270", N2_vs_BR5270_higlow))

# N2 < BR5270

all_lev2_3_df_one_compound_w_data_more_codes_BR5270_N2 <- filter(all_lev2_3_df_one_compound_w_data_more_codes, grepl("BR5270 > N2", N2_vs_BR5270_higlow))


# N2 or BR5270

all_lev2_3_df_one_compound_w_data_more_codes_N2orBR5270 <- filter(all_lev2_3_df_one_compound_w_data_more_codes, grepl("N2 > BR5270", N2_vs_BR5270_higlow) | grepl("BR5270 > N2", N2_vs_BR5270_higlow))


## export:


#### only those significant in some with the code:

write_csv(all_lev2_3_df_one_compound_w_data_more_codes_without_N2, "onlysigninsomewcode_00_all_lev2_3_df_one_compound_w_data_more_codes_without_N2.csv")
write_csv(all_lev2_3_df_one_compound_w_data_more_codes_N2_VC40, "onlysigninsomewcode_01_all_lev2_3_df_one_compound_w_data_more_codes_N2_VC40.csv")
write_csv(all_lev2_3_df_one_compound_w_data_more_codes_VC40_N2, "onlysigninsomewcode_02_all_lev2_3_df_one_compound_w_data_more_codes_VC40_N2.csv")
write_csv(all_lev2_3_df_one_compound_w_data_more_codes_N2orVC40, "onlysigninsomewcode_03_all_lev2_3_df_one_compound_w_data_more_codes_N2orVC40.csv")
write_csv(all_lev2_3_df_one_compound_w_data_more_codes_N2_VC1668, "onlysigninsomewcode_04_all_lev2_3_df_one_compound_w_data_more_codes_N2_VC1668.csv")
write_csv(all_lev2_3_df_one_compound_w_data_more_codes_VC1668_N2, "onlysigninsomewcode_05_all_lev2_3_df_one_compound_w_data_more_codes_VC1668_N2.csv")
write_csv(all_lev2_3_df_one_compound_w_data_more_codes_N2orVC1668, "onlysigninsomewcode_06_all_lev2_3_df_one_compound_w_data_more_codes_N2orVC1668.csv")
write_csv(all_lev2_3_df_one_compound_w_data_more_codes_N2_UA57, "onlysigninsomewcode_07_all_lev2_3_df_one_compound_w_data_more_codes_N2_UA57.csv")
write_csv(all_lev2_3_df_one_compound_w_data_more_codes_UA57_N2, "onlysigninsomewcode_08_all_lev2_3_df_one_compound_w_data_more_codes_UA57_N2.csv")
write_csv(all_lev2_3_df_one_compound_w_data_more_codes_N2orUA57, "onlysigninsomewcode_09_all_lev2_3_df_one_compound_w_data_more_codes_N2orUA57.csv")
write_csv(all_lev2_3_df_one_compound_w_data_more_codes_N2_BR5270, "onlysigninsomewcode_10_all_lev2_3_df_one_compound_w_data_more_codes_N2_BR5270.csv")
write_csv(all_lev2_3_df_one_compound_w_data_more_codes_BR5270_N2, "onlysigninsomewcode_11_all_lev2_3_df_one_compound_w_data_more_codes_BR5270_N2.csv")
write_csv(all_lev2_3_df_one_compound_w_data_more_codes_N2orBR5270, "onlysigninsomewcode_12_all_lev2_3_df_one_compound_w_data_more_codes_N2orBR5270.csv")




## Doing FELLA now:



### actually using the FELLA package now (only relevant part for now):

## N2orVC40

Significative_molecules_KEGG_codes_N2orVC40 <- all_lev2_3_df_one_compound_w_data_more_codes_N2orVC40$KEGGcode[which(!is.na(all_lev2_3_df_one_compound_w_data_more_codes_N2orVC40$KEGGcode))]

enrich_analysis_sign_molecules_N2orVC40 <- defineCompounds(compounds = Significative_molecules_KEGG_codes_N2orVC40,
                                                           data = fella.data_c_elegans)

print(paste0(length(getInput(enrich_analysis_sign_molecules_N2orVC40)), " codes used; ",
             length(getExcluded(enrich_analysis_sign_molecules_N2orVC40)), " codes excluded; ",
             length(Significative_molecules_KEGG_codes_N2orVC40), " total codes"))


enrich_analysis_sign_molecules_N2orVC40_run <- runDiffusion(object = enrich_analysis_sign_molecules_N2orVC40,
                                                            data = fella.data_c_elegans,
                                                            approx = "normality")

enrich_analysis_sign_molecules_N2orVC40_run_resultsGraph <- generateResultsGraph(object = enrich_analysis_sign_molecules_N2orVC40_run,
                                                                                 method = "diffusion",
                                                                                 data = fella.data_c_elegans,
                                                                                 nlimit = 250)

exportResults(format = "csv", file = "onlysigninsomewcodeFELLA_enrich_analysis_sign_molecules_N2orVC40_run_resultsGraph.csv",
              method = "diffusion", object = enrich_analysis_sign_molecules_N2orVC40_run, data = fella.data_c_elegans)


generated_table_N2orVC40 <- generateResultsTable(method = "diffusion",
                                                 nlimit = 1000, 
                                                 #threshold = 0.05, plimit = 15, LabelLengthAtPlot = 45,
                                                 #capPscores = 1e-06,
                                                 object = enrich_analysis_sign_molecules_N2orVC40_run, data = fella.data_c_elegans)

plotGraph(graph = enrich_analysis_sign_molecules_N2orVC40_run_resultsGraph,
          #layout = FALSE, graph.layout = NULL,
          #plotLegend = TRUE, plot.fun = "plot.igraph", NamesAsLabels = TRUE,
          vertex.label.cex = 0.5)                                                                     



write_graph(enrich_analysis_sign_molecules_N2orVC40_run_resultsGraph, "onlysigninsomewcodeFELLA_enrich_analysis_graph_to_export_N2orVC40.png")

png("onlysigninsomewcodeFELLA_enrich_analysis_sign_molecules_N2orVC40_run_resultsGraph.png", width = 1000, height = 1000)
plotGraph(graph = enrich_analysis_sign_molecules_N2orVC40_run_resultsGraph,
          #layout = FALSE, graph.layout = NULL,
          #plotLegend = TRUE, plot.fun = "plot.igraph", NamesAsLabels = TRUE,
          vertex.label.cex = 1)
dev.off()





## N2orVC1668

Significative_molecules_KEGG_codes_N2orVC1668 <- all_lev2_3_df_one_compound_w_data_more_codes_N2orVC1668$KEGGcode[which(!is.na(all_lev2_3_df_one_compound_w_data_more_codes_N2orVC1668$KEGGcode))]

enrich_analysis_sign_molecules_N2orVC1668 <- defineCompounds(compounds = Significative_molecules_KEGG_codes_N2orVC1668,
                                                             data = fella.data_c_elegans)

print(paste0(length(getInput(enrich_analysis_sign_molecules_N2orVC1668)), " codes used; ",
             length(getExcluded(enrich_analysis_sign_molecules_N2orVC1668)), " codes excluded; ",
             length(Significative_molecules_KEGG_codes_N2orVC1668), " total codes"))


enrich_analysis_sign_molecules_N2orVC1668_run <- runDiffusion(object = enrich_analysis_sign_molecules_N2orVC1668,
                                                              data = fella.data_c_elegans,
                                                              approx = "normality")

enrich_analysis_sign_molecules_N2orVC1668_run_resultsGraph <- generateResultsGraph(object = enrich_analysis_sign_molecules_N2orVC1668_run,
                                                                                   method = "diffusion",
                                                                                   data = fella.data_c_elegans,
                                                                                   nlimit = 250)

exportResults(format = "csv", file = "onlysigninsomewcodeFELLA_enrich_analysis_sign_molecules_N2orVC1668_run_resultsGraph.csv",
              method = "diffusion", object = enrich_analysis_sign_molecules_N2orVC1668_run, data = fella.data_c_elegans)


generated_table_N2orVC1668 <- generateResultsTable(method = "diffusion",
                                                   nlimit = 1000, 
                                                   #threshold = 0.05, plimit = 15, LabelLengthAtPlot = 45,
                                                   #capPscores = 1e-06,
                                                   object = enrich_analysis_sign_molecules_N2orVC1668_run, data = fella.data_c_elegans)

plotGraph(graph = enrich_analysis_sign_molecules_N2orVC1668_run_resultsGraph,
          #layout = FALSE, graph.layout = NULL,
          #plotLegend = TRUE, plot.fun = "plot.igraph", NamesAsLabels = TRUE,
          vertex.label.cex = 0.5)                                                                     



write_graph(enrich_analysis_sign_molecules_N2orVC1668_run_resultsGraph, "onlysigninsomewcodeFELLA_enrich_analysis_graph_to_export_N2orVC1668.png")

png("onlysigninsomewcodeFELLA_enrich_analysis_sign_molecules_N2orVC1668_run_resultsGraph.png", width = 1000, height = 1000)
plotGraph(graph = enrich_analysis_sign_molecules_N2orVC1668_run_resultsGraph,
          #layout = FALSE, graph.layout = NULL,
          #plotLegend = TRUE, plot.fun = "plot.igraph", NamesAsLabels = TRUE,
          vertex.label.cex = 1)
dev.off()




## N2orUA57

Significative_molecules_KEGG_codes_N2orUA57 <- all_lev2_3_df_one_compound_w_data_more_codes_N2orUA57$KEGGcode[which(!is.na(all_lev2_3_df_one_compound_w_data_more_codes_N2orUA57$KEGGcode))]

enrich_analysis_sign_molecules_N2orUA57 <- defineCompounds(compounds = Significative_molecules_KEGG_codes_N2orUA57,
                                                           data = fella.data_c_elegans)

print(paste0(length(getInput(enrich_analysis_sign_molecules_N2orUA57)), " codes used; ",
             length(getExcluded(enrich_analysis_sign_molecules_N2orUA57)), " codes excluded; ",
             length(Significative_molecules_KEGG_codes_N2orUA57), " total codes"))


enrich_analysis_sign_molecules_N2orUA57_run <- runDiffusion(object = enrich_analysis_sign_molecules_N2orUA57,
                                                            data = fella.data_c_elegans,
                                                            approx = "normality")

enrich_analysis_sign_molecules_N2orUA57_run_resultsGraph <- generateResultsGraph(object = enrich_analysis_sign_molecules_N2orUA57_run,
                                                                                 method = "diffusion",
                                                                                 data = fella.data_c_elegans,
                                                                                 nlimit = 250)

exportResults(format = "csv", file = "onlysigninsomewcodeFELLA_enrich_analysis_sign_molecules_N2orUA57_run_resultsGraph.csv",
              method = "diffusion", object = enrich_analysis_sign_molecules_N2orUA57_run, data = fella.data_c_elegans)


generated_table_N2orUA57 <- generateResultsTable(method = "diffusion",
                                                 nlimit = 1000, 
                                                 #threshold = 0.05, plimit = 15, LabelLengthAtPlot = 45,
                                                 #capPscores = 1e-06,
                                                 object = enrich_analysis_sign_molecules_N2orUA57_run, data = fella.data_c_elegans)

plotGraph(graph = enrich_analysis_sign_molecules_N2orUA57_run_resultsGraph,
          #layout = FALSE, graph.layout = NULL,
          #plotLegend = TRUE, plot.fun = "plot.igraph", NamesAsLabels = TRUE,
          vertex.label.cex = 0.5)                                                                     



write_graph(enrich_analysis_sign_molecules_N2orUA57_run_resultsGraph, "onlysigninsomewcodeFELLA_enrich_analysis_graph_to_export_N2orUA57.png")

png("onlysigninsomewcodeFELLA_enrich_analysis_sign_molecules_N2orUA57_run_resultsGraph.png", width = 1000, height = 1000)
plotGraph(graph = enrich_analysis_sign_molecules_N2orUA57_run_resultsGraph,
          #layout = FALSE, graph.layout = NULL,
          #plotLegend = TRUE, plot.fun = "plot.igraph", NamesAsLabels = TRUE,
          vertex.label.cex = 1)
dev.off()



## N2orBR5270

Significative_molecules_KEGG_codes_N2orBR5270 <- all_lev2_3_df_one_compound_w_data_more_codes_N2orBR5270$KEGGcode[which(!is.na(all_lev2_3_df_one_compound_w_data_more_codes_N2orBR5270$KEGGcode))]

enrich_analysis_sign_molecules_N2orBR5270 <- defineCompounds(compounds = Significative_molecules_KEGG_codes_N2orBR5270,
                                                             data = fella.data_c_elegans)

print(paste0(length(getInput(enrich_analysis_sign_molecules_N2orBR5270)), " codes used; ",
             length(getExcluded(enrich_analysis_sign_molecules_N2orBR5270)), " codes excluded; ",
             length(Significative_molecules_KEGG_codes_N2orBR5270), " total codes"))


enrich_analysis_sign_molecules_N2orBR5270_run <- runDiffusion(object = enrich_analysis_sign_molecules_N2orBR5270,
                                                              data = fella.data_c_elegans,
                                                              approx = "normality")

enrich_analysis_sign_molecules_N2orBR5270_run_resultsGraph <- generateResultsGraph(object = enrich_analysis_sign_molecules_N2orBR5270_run,
                                                                                   method = "diffusion",
                                                                                   data = fella.data_c_elegans,
                                                                                   nlimit = 250)

exportResults(format = "csv", file = "onlysigninsomewcodeFELLA_enrich_analysis_sign_molecules_N2orBR5270_run_resultsGraph.csv",
              method = "diffusion", object = enrich_analysis_sign_molecules_N2orBR5270_run, data = fella.data_c_elegans)


generated_table_N2orBR5270 <- generateResultsTable(method = "diffusion",
                                                   nlimit = 1000, 
                                                   #threshold = 0.05, plimit = 15, LabelLengthAtPlot = 45,
                                                   #capPscores = 1e-06,
                                                   object = enrich_analysis_sign_molecules_N2orBR5270_run, data = fella.data_c_elegans)

plotGraph(graph = enrich_analysis_sign_molecules_N2orBR5270_run_resultsGraph,
          #layout = FALSE, graph.layout = NULL,
          #plotLegend = TRUE, plot.fun = "plot.igraph", NamesAsLabels = TRUE,
          vertex.label.cex = 0.5)                                                                     



write_graph(enrich_analysis_sign_molecules_N2orBR5270_run_resultsGraph, "onlysigninsomewcodeFELLA_enrich_analysis_graph_to_export_N2orBR5270.png")

png("onlysigninsomewcodeFELLA_enrich_analysis_sign_molecules_N2orBR5270_run_resultsGraph.png", width = 1000, height = 1000)
plotGraph(graph = enrich_analysis_sign_molecules_N2orBR5270_run_resultsGraph,
          #layout = FALSE, graph.layout = NULL,
          #plotLegend = TRUE, plot.fun = "plot.igraph", NamesAsLabels = TRUE,
          vertex.label.cex = 1)
dev.off()




####



