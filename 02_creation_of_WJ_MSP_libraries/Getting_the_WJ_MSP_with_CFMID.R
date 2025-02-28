

##### Getting MSP with also fragments:

library(tidyverse)
library(webchem)
library(chem)
library(RChemMass)


WormJam_combined <- read_csv("20240829_C_elegans_all_Struct_Info.csv")

WormJam_combined <- mutate(WormJam_combined, Title = as.character(NA))
for (i in 1:length(pull(WormJam_combined, 1))) {
  WormJam_combined[i, "Title"] <- getPCdesc.title(query = WormJam_combined$Identifier[i], from = "cid", timeout=10)$Title
}


## some checks (should be all FALSE):

any(duplicated(WormJam_combined$Identifier))
any(duplicated(WormJam_combined$SMILES))
any(duplicated(WormJam_combined$InChI))
any(duplicated(WormJam_combined$InChIKey))
any(duplicated(WormJam_combined$CompoundName))
any(duplicated(WormJam_combined$Title))

any(is.na(WormJam_combined$Identifier))
any(is.na(WormJam_combined$SMILES))
any(is.na(WormJam_combined$InChI))
any(is.na(WormJam_combined$InChIKey))
any(is.na(WormJam_combined$CompoundName))
any(is.na(WormJam_combined$Title))

## there are some duplicates in the InChI and in the InChIKey. It shuld be fine as I don't use them.

## extracting info to run CFMID:

WormJam_combined_SMILES <- select(WormJam_combined, Identifier, SMILES)
write_tsv(WormJam_combined_SMILES, "WormJam_combined_SMILES.txt")

## actually that is too long, so reduced in subpart:


WormJam_combined_SMILES_part1of9 <- WormJam_combined_SMILES[1:500,]
write_tsv(WormJam_combined_SMILES_part1of9, "WormJam_combined_SMILES_part1of9.txt")

WormJam_combined_SMILES_part2of9 <- WormJam_combined_SMILES[501:1000,]
write_tsv(WormJam_combined_SMILES_part2of9, "WormJam_combined_SMILES_part2of9.txt")

WormJam_combined_SMILES_part3of9 <- WormJam_combined_SMILES[1001:1500,]
write_tsv(WormJam_combined_SMILES_part3of9, "WormJam_combined_SMILES_part3of9.txt")

WormJam_combined_SMILES_part4of9 <- WormJam_combined_SMILES[1501:2000,]
write_tsv(WormJam_combined_SMILES_part4of9, "WormJam_combined_SMILES_part4of9.txt")

WormJam_combined_SMILES_part5of9 <- WormJam_combined_SMILES[2001:2500,]
write_tsv(WormJam_combined_SMILES_part5of9, "WormJam_combined_SMILES_part5of9.txt")

WormJam_combined_SMILES_part6of9 <- WormJam_combined_SMILES[2501:3000,]
write_tsv(WormJam_combined_SMILES_part6of9, "WormJam_combined_SMILES_part6of9.txt")

WormJam_combined_SMILES_part7of9 <- WormJam_combined_SMILES[3001:3500,]
write_tsv(WormJam_combined_SMILES_part7of9, "WormJam_combined_SMILES_part7of9.txt")

WormJam_combined_SMILES_part8of9 <- WormJam_combined_SMILES[3501:4000,]
write_tsv(WormJam_combined_SMILES_part8of9, "WormJam_combined_SMILES_part8of9.txt")

WormJam_combined_SMILES_part9of9 <- WormJam_combined_SMILES[4000:length(pull(WormJam_combined_SMILES, 1)),]
write_tsv(WormJam_combined_SMILES_part9of9, "WormJam_combined_SMILES_part9of9.txt")



## before running the code, eliminate the first row (the col names)

#######
########## fragmentation simulation performed on CFM-ID command line
#######

# Code used for this project:
  
# POS:

## total - NOT COMPLEATED:  
# docker run --rm=true -v ${pwd}:/cfmid/public/ -i wishartlab/cfmid:latest sh -c "cd /cfmid/public/; cfm-predict 'WormJam_combined_SMILES.txt' 0.001 /trained_models_cfmid4.0/[M+H]+/param_output.log /trained_models_cfmid4.0/[M+H]+/param_config.txt 1 WormJam_combined_SMILES_predicted_output_POS.txt"


## 1of9 - DONE!

# docker run --rm=true -v ${pwd}:/cfmid/public/ -i wishartlab/cfmid:latest sh -c "cd /cfmid/public/; cfm-predict 'WormJam_combined_SMILES_part1of9.txt' 0.001 /trained_models_cfmid4.0/[M+H]+/param_output.log /trained_models_cfmid4.0/[M+H]+/param_config.txt 1 WormJam_combined_SMILES_predicted_output_POS_part1of9.txt"


## 2of9 - DONE!

# docker run --rm=true -v ${pwd}:/cfmid/public/ -i wishartlab/cfmid:latest sh -c "cd /cfmid/public/; cfm-predict 'WormJam_combined_SMILES_part2of9.txt' 0.001 /trained_models_cfmid4.0/[M+H]+/param_output.log /trained_models_cfmid4.0/[M+H]+/param_config.txt 1 WormJam_combined_SMILES_predicted_output_POS_part2of9.txt"


## 3of9 - DONE!

# docker run --rm=true -v ${pwd}:/cfmid/public/ -i wishartlab/cfmid:latest sh -c "cd /cfmid/public/; cfm-predict 'WormJam_combined_SMILES_part3of9.txt' 0.001 /trained_models_cfmid4.0/[M+H]+/param_output.log /trained_models_cfmid4.0/[M+H]+/param_config.txt 1 WormJam_combined_SMILES_predicted_output_POS_part3of9.txt"


## 4of9 - DONE!

# docker run --rm=true -v ${pwd}:/cfmid/public/ -i wishartlab/cfmid:latest sh -c "cd /cfmid/public/; cfm-predict 'WormJam_combined_SMILES_part4of9.txt' 0.001 /trained_models_cfmid4.0/[M+H]+/param_output.log /trained_models_cfmid4.0/[M+H]+/param_config.txt 1 WormJam_combined_SMILES_predicted_output_POS_part4of9.txt"

# the laptop shut down, when 70 molecules where left. So, run a "_02":

# docker run --rm=true -v ${pwd}:/cfmid/public/ -i wishartlab/cfmid:latest sh -c "cd /cfmid/public/; cfm-predict 'WormJam_combined_SMILES_part4of9_02.txt' 0.001 /trained_models_cfmid4.0/[M+H]+/param_output.log /trained_models_cfmid4.0/[M+H]+/param_config.txt 1 WormJam_combined_SMILES_predicted_output_POS_part4of9_02.txt"



## 5of9 - DONE!

# docker run --rm=true -v ${pwd}:/cfmid/public/ -i wishartlab/cfmid:latest sh -c "cd /cfmid/public/; cfm-predict 'WormJam_combined_SMILES_part5of9.txt' 0.001 /trained_models_cfmid4.0/[M+H]+/param_output.log /trained_models_cfmid4.0/[M+H]+/param_config.txt 1 WormJam_combined_SMILES_predicted_output_POS_part5of9.txt"


## 6of9 - DONE!

# docker run --rm=true -v ${pwd}:/cfmid/public/ -i wishartlab/cfmid:latest sh -c "cd /cfmid/public/; cfm-predict 'WormJam_combined_SMILES_part6of9.txt' 0.001 /trained_models_cfmid4.0/[M+H]+/param_output.log /trained_models_cfmid4.0/[M+H]+/param_config.txt 1 WormJam_combined_SMILES_predicted_output_POS_part6of9.txt"


## 7of9 - DONE!

# docker run --rm=true -v ${pwd}:/cfmid/public/ -i wishartlab/cfmid:latest sh -c "cd /cfmid/public/; cfm-predict 'WormJam_combined_SMILES_part7of9.txt' 0.001 /trained_models_cfmid4.0/[M+H]+/param_output.log /trained_models_cfmid4.0/[M+H]+/param_config.txt 1 WormJam_combined_SMILES_predicted_output_POS_part7of9.txt"


## 8of9 - DONE!

# docker run --rm=true -v ${pwd}:/cfmid/public/ -i wishartlab/cfmid:latest sh -c "cd /cfmid/public/; cfm-predict 'WormJam_combined_SMILES_part8of9.txt' 0.001 /trained_models_cfmid4.0/[M+H]+/param_output.log /trained_models_cfmid4.0/[M+H]+/param_config.txt 1 WormJam_combined_SMILES_predicted_output_POS_part8of9.txt"


## 9of9 - DONE!

# docker run --rm=true -v ${pwd}:/cfmid/public/ -i wishartlab/cfmid:latest sh -c "cd /cfmid/public/; cfm-predict 'WormJam_combined_SMILES_part9of9.txt' 0.001 /trained_models_cfmid4.0/[M+H]+/param_output.log /trained_models_cfmid4.0/[M+H]+/param_config.txt 1 WormJam_combined_SMILES_predicted_output_POS_part9of9.txt"





# NEG:

## total - NOT DONE AT ALL:
# docker run --rm=true -v ${pwd}:/cfmid/public/ -i wishartlab/cfmid:latest sh -c "cd /cfmid/public/; cfm-predict 'WormJam_combined_SMILES.txt' 0.001 /trained_models_cfmid4.0/[M-H]-/param_output.log /trained_models_cfmid4.0/[M-H]-/param_config.txt 1 WormJam_combined_SMILES_predicted_output_NEG.txt"


## 1of9 - DONE!

# docker run --rm=true -v ${pwd}:/cfmid/public/ -i wishartlab/cfmid:latest sh -c "cd /cfmid/public/; cfm-predict 'WormJam_combined_SMILES_part1of9.txt' 0.001 /trained_models_cfmid4.0/[M-H]-/param_output.log /trained_models_cfmid4.0/[M-H]-/param_config.txt 1 WormJam_combined_SMILES_predicted_output_NEG_part1of9.txt"


## 2of9 - DONE!

# docker run --rm=true -v ${pwd}:/cfmid/public/ -i wishartlab/cfmid:latest sh -c "cd /cfmid/public/; cfm-predict 'WormJam_combined_SMILES_part2of9.txt' 0.001 /trained_models_cfmid4.0/[M-H]-/param_output.log /trained_models_cfmid4.0/[M-H]-/param_config.txt 1 WormJam_combined_SMILES_predicted_output_NEG_part2of9.txt"


## 3of9 - DONE!

# docker run --rm=true -v ${pwd}:/cfmid/public/ -i wishartlab/cfmid:latest sh -c "cd /cfmid/public/; cfm-predict 'WormJam_combined_SMILES_part3of9.txt' 0.001 /trained_models_cfmid4.0/[M-H]-/param_output.log /trained_models_cfmid4.0/[M-H]-/param_config.txt 1 WormJam_combined_SMILES_predicted_output_NEG_part3of9.txt"


## 4of9 - DONE!

# docker run --rm=true -v ${pwd}:/cfmid/public/ -i wishartlab/cfmid:latest sh -c "cd /cfmid/public/; cfm-predict 'WormJam_combined_SMILES_part4of9.txt' 0.001 /trained_models_cfmid4.0/[M-H]-/param_output.log /trained_models_cfmid4.0/[M-H]-/param_config.txt 1 WormJam_combined_SMILES_predicted_output_NEG_part4of9.txt"


## 5of9 - DONE! 

# docker run --rm=true -v ${pwd}:/cfmid/public/ -i wishartlab/cfmid:latest sh -c "cd /cfmid/public/; cfm-predict 'WormJam_combined_SMILES_part5of9.txt' 0.001 /trained_models_cfmid4.0/[M-H]-/param_output.log /trained_models_cfmid4.0/[M-H]-/param_config.txt 1 WormJam_combined_SMILES_predicted_output_NEG_part5of9.txt"


## 6of9 - DONE!

# docker run --rm=true -v ${pwd}:/cfmid/public/ -i wishartlab/cfmid:latest sh -c "cd /cfmid/public/; cfm-predict 'WormJam_combined_SMILES_part6of9.txt' 0.001 /trained_models_cfmid4.0/[M-H]-/param_output.log /trained_models_cfmid4.0/[M-H]-/param_config.txt 1 WormJam_combined_SMILES_predicted_output_NEG_part6of9.txt"


## 7of9 - DONE!

# docker run --rm=true -v ${pwd}:/cfmid/public/ -i wishartlab/cfmid:latest sh -c "cd /cfmid/public/; cfm-predict 'WormJam_combined_SMILES_part7of9.txt' 0.001 /trained_models_cfmid4.0/[M-H]-/param_output.log /trained_models_cfmid4.0/[M-H]-/param_config.txt 1 WormJam_combined_SMILES_predicted_output_NEG_part7of9.txt"


## 8of9 - DONE!

# docker run --rm=true -v ${pwd}:/cfmid/public/ -i wishartlab/cfmid:latest sh -c "cd /cfmid/public/; cfm-predict 'WormJam_combined_SMILES_part8of9.txt' 0.001 /trained_models_cfmid4.0/[M-H]-/param_output.log /trained_models_cfmid4.0/[M-H]-/param_config.txt 1 WormJam_combined_SMILES_predicted_output_NEG_part8of9.txt"


## 9of9 - running:

# docker run --rm=true -v ${pwd}:/cfmid/public/ -i wishartlab/cfmid:latest sh -c "cd /cfmid/public/; cfm-predict 'WormJam_combined_SMILES_part9of9.txt' 0.001 /trained_models_cfmid4.0/[M-H]-/param_output.log /trained_models_cfmid4.0/[M-H]-/param_config.txt 1 WormJam_combined_SMILES_predicted_output_NEG_part9of9.txt"








#####



creating_in_silico_fragmentation_list <- function(readLines_char_vector) {
  
  output_list <- vector(mode = "list", length = length(which(grepl("#In-silico", readLines_char_vector))))
  
  for (i in 1:length(output_list)) {
    
    if (i != length(output_list)) {
      this_molecule_vector <- readLines_char_vector[which(grepl("#In-silico", readLines_char_vector))[i]:((which(grepl("#In-silico", readLines_char_vector))[i+1])-1)]
    } else {
      this_molecule_vector <- readLines_char_vector[which(grepl("#In-silico", readLines_char_vector))[i]:length(readLines_char_vector)]
    }
    
    output_list[[i]] <- list(ID = str_remove_all(this_molecule_vector[which(grepl("#ID=", this_molecule_vector))], "#ID="),
                             SMILES = str_remove_all(this_molecule_vector[which(grepl("#SMILES=", this_molecule_vector))], "#SMILES="),
                             InChI = str_remove_all(this_molecule_vector[which(grepl("#InChI=", this_molecule_vector))], "#InChI="),
                             InChiKey = str_remove_all(this_molecule_vector[which(grepl("#InChiKey=", this_molecule_vector))], "#InChiKey="),
                             Formula = str_remove_all(this_molecule_vector[which(grepl("#Formula=", this_molecule_vector))], "#Formula="),
                             PMass = str_remove_all(this_molecule_vector[which(grepl("#PMass=", this_molecule_vector))], "#PMass="),
                             all_else = this_molecule_vector[which(!(grepl("#In-silico", this_molecule_vector) |
                                                                       grepl("#PREDICTED", this_molecule_vector) |
                                                                       grepl("#ID=", this_molecule_vector) |
                                                                       grepl("#SMILES=", this_molecule_vector) |
                                                                       grepl("#InChI=", this_molecule_vector) |  
                                                                       grepl("#InChiKey=", this_molecule_vector) |
                                                                       grepl("#Formula=", this_molecule_vector) |
                                                                       grepl("#PMass=", this_molecule_vector)))])
    
    names(output_list)[i] <- str_remove_all(this_molecule_vector[which(grepl("#ID=", this_molecule_vector))], "#ID=")
    
  }
  
  return(output_list)
  
}


getting_summary_fragment_matrix <- function(all_else_from_CFMID, weight_energy0 = 1/3, weight_energy1 = 1/3, weight_energy2 = 1/3) {
  
  
  energy0_data <- all_else_from_CFMID[(which(all_else_from_CFMID == "energy0")+1):(which(all_else_from_CFMID == "energy1")-1)]
  energy0_fragm <- as.numeric(sub(" .*", "", energy0_data))
  energy0_intensity <- as.numeric(sub("^[^ ]+ ([^ ]+).*", "\\1", energy0_data))
  
  matrix_energy0 <- matrix(data =c(energy0_fragm, energy0_intensity),
                           nrow = length(energy0_data),
                           ncol = 2,
                           byrow = FALSE,
                           dimnames = list(NULL, c("fragment", "intensity")))
  
  
  energy1_data <- all_else_from_CFMID[(which(all_else_from_CFMID == "energy1")+1):(which(all_else_from_CFMID == "energy2")-1)]
  energy1_fragm <- as.numeric(sub(" .*", "", energy1_data))
  energy1_intensity <- as.numeric(sub("^[^ ]+ ([^ ]+).*", "\\1", energy1_data))
  
  matrix_energy1 <- matrix(data =c(energy1_fragm, energy1_intensity),
                           nrow = length(energy1_data),
                           ncol = 2,
                           byrow = FALSE,
                           dimnames = list(NULL, c("fragment", "intensity")))
  
  
  energy2_data <- all_else_from_CFMID[(which(all_else_from_CFMID == "energy2")+1):(which(all_else_from_CFMID == "")[1]-1)]
  energy2_fragm <- as.numeric(sub(" .*", "", energy2_data))
  energy2_intensity <- as.numeric(sub("^[^ ]+ ([^ ]+).*", "\\1", energy2_data))
  
  matrix_energy2 <- matrix(data =c(energy2_fragm, energy2_intensity),
                           nrow = length(energy2_data),
                           ncol = 2,
                           byrow = FALSE,
                           dimnames = list(NULL, c("fragment", "intensity")))
  
  
  
  matrix_energy0_weighted <- matrix_energy0
  matrix_energy0_weighted[,2] <- matrix_energy0_weighted[,2]*weight_energy0
  
  matrix_energy1_weighted <- matrix_energy1
  matrix_energy1_weighted[,2] <- matrix_energy1_weighted[,2]*weight_energy1
  
  matrix_energy2_weighted <- matrix_energy2
  matrix_energy2_weighted[,2] <- matrix_energy2_weighted[,2]*weight_energy2
  
  
  matrix_energy_summed <- rbind(matrix_energy0_weighted, matrix_energy1_weighted, matrix_energy2_weighted)
  matrix_energy_summed_cleaned <- matrix_energy_summed
  
  if (any(duplicated(matrix_energy_summed[,1]))) {
    
    row_to_remove <- numeric()
    
    for (d in unique(matrix_energy_summed[,1][duplicated(matrix_energy_summed[,1])])) {
      
      i_rows_this_dupl <- which(matrix_energy_summed[,1] == d)
      
      this_dupl_intensity_sum <- sum(matrix_energy_summed[,2][i_rows_this_dupl])
      
      matrix_energy_summed_cleaned[i_rows_this_dupl[1],2] <- this_dupl_intensity_sum
      
      row_to_remove <- c(row_to_remove, i_rows_this_dupl[2:length(i_rows_this_dupl)])
    }
    
    matrix_energy_summed_cleaned <- matrix_energy_summed_cleaned[-(row_to_remove),]
  }
  
  return(matrix_energy_summed_cleaned)
  
}





creating_MSP_text <- function(compound_table, id_col = NULL, name_col = "NAME", CID_col = NULL, INCHIKEY_col = NULL, SMILES_col = NULL, InChI_col = NULL, MOLECULAR_FORMULA_col = NULL, NEUTRAL_MASS_col = "MONOISOTOPIC_MASS", Ontology_col = NULL,
                              polarity = c("NEG", "POS"), CFMID_output_from_SMILES_list = NULL, CFMID_output_from_InChI_list = NULL, ID_used = NULL,
                              weight_energy0 = 1/3, weight_energy1 = 1/3, weight_energy2 = 1/3) {

  
  if (!identical(tolower(polarity), c("neg", "pos"))) {
    if (length(polarity) != 1) {stop('polarity must be one of "NEG", "POS"')}
    if (is.na(polarity)) {stop('polarity must be one of "NEG", "POS"')}
  }
  polarity <- tolower(polarity)
  polarity <- match.arg(polarity, c("neg", "pos"))
  
  final_MSP_text <- character()
  
  if (any(duplicated(pull(compound_table, name_col)))) {warning("There are duplicated names in the compound_table, this is not ideal...!!!!")}
  
  for (i in 1:length(pull(compound_table, 1))) {
    
    this_tibble <- tibble(NAME = pull(compound_table, name_col)[i],
                          PRECURSORMZ = ifelse(polarity == "neg", as.numeric(pull(compound_table, NEUTRAL_MASS_col)[i]) - 1.0078, as.numeric(pull(compound_table, NEUTRAL_MASS_col)[i]) + 1.0078),
                          PRECURSORTYPE = ifelse(polarity == "neg", "[M-H]-", "[M+H]+"),
                          IONMODE = ifelse(polarity == "neg", "Negative", "Positive"),
                          RETENTIONTIME = 1,
                          CCS = NA,
                          FORMULA = ifelse(is.null(MOLECULAR_FORMULA_col), NA, pull(compound_table, MOLECULAR_FORMULA_col)[i]),
                          ONTOLOGY = ifelse(is.null(Ontology_col), NA, pull(compound_table, Ontology_col)[i]),
                          SMILES = ifelse(is.null(SMILES_col), NA, pull(compound_table, SMILES_col)[i]),
                          INCHIKEY = ifelse(is.null(INCHIKEY_col), NA, pull(compound_table, INCHIKEY_col)[i]),
                          INSTRUMENTTYPE = NA,
                          COLLISIONENERGY = NA,
                          COMMENT_ID = ifelse(is.null(id_col), NA, pull(compound_table, id_col)[i]),
                          COMMENT_CID = ifelse(is.null(CID_col), NA, pull(compound_table, CID_col)[i]),
                          COMMENT_InChI = ifelse(is.null(InChI_col), NA, pull(compound_table, InChI_col)[i]),
                          FRAGMENTS = FALSE)
    
    ID_used_vector <- as.character(pull(compound_table, ID_used))
    
    if (!is.null(CFMID_output_from_SMILES_list)) {
      if (ID_used_vector[i] %in% names(CFMID_output_from_SMILES_list)) {
        present_in_CFMID_output_from_SMILES_list <- TRUE
        this_frag_data_fromSMILES_list <- CFMID_output_from_SMILES_list[[ID_used_vector[i]]]
      } else {
        present_in_CFMID_output_from_SMILES_list <- FALSE
      }
    } else {
      present_in_CFMID_output_from_SMILES_list <- FALSE
    }
    
    if (!is.null(CFMID_output_from_InChI_list)) {
      if (ID_used_vector[i] %in% names(CFMID_output_from_InChI_list)) {
        present_in_CFMID_output_from_InChI_list <- TRUE
        this_frag_data_fromSInChI_list <- CFMID_output_from_InChI_list[[ID_used_vector[i]]]
      } else {
        present_in_CFMID_output_from_InChI_list <- FALSE
      }
    } else {
      present_in_CFMID_output_from_InChI_list <- FALSE
    }
    
    
    if (present_in_CFMID_output_from_SMILES_list) {
      if (length(this_frag_data_fromSMILES_list$ID) == 1) {this_tibble$COMMENT_ID <- this_frag_data_fromSMILES_list$ID}
      if (length(this_frag_data_fromSMILES_list$SMILES) == 1) {this_tibble$SMILES <- this_frag_data_fromSMILES_list$SMILES}
      if (length(this_frag_data_fromSMILES_list$InChI) == 1) {this_tibble$COMMENT_InChI <- this_frag_data_fromSMILES_list$InChI}
      if (length(this_frag_data_fromSMILES_list$InChiKey) == 1) {this_tibble$INCHIKEY <- this_frag_data_fromSMILES_list$InChiKey}
      if (length(this_frag_data_fromSMILES_list$Formula) == 1) {this_tibble$FORMULA <- this_frag_data_fromSMILES_list$Formula}
      if (length(this_frag_data_fromSMILES_list$PMass) == 1) {this_tibble$PRECURSORMZ <- this_frag_data_fromSMILES_list$PMass}
      
      this_tibble$FRAGMENTS <- TRUE
      
      fragment_matrix <- getting_summary_fragment_matrix(all_else_from_CFMID = this_frag_data_fromSMILES_list$all_else, weight_energy0 = weight_energy0, weight_energy1 = weight_energy1, weight_energy2 = weight_energy2)
      
      
    } else if (present_in_CFMID_output_from_InChI_list) {
      if (length(this_frag_data_fromSInChI_list$ID) == 1) {this_tibble$COMMENT_ID <- this_frag_data_fromSInChI_list$ID}
      if (length(this_frag_data_fromSInChI_list$SMILES) == 1) {this_tibble$SMILES <- this_frag_data_fromSInChI_list$SMILES}
      if (length(this_frag_data_fromSInChI_list$InChI) == 1) {this_tibble$COMMENT_InChI <- this_frag_data_fromSInChI_list$InChI}
      if (length(this_frag_data_fromSInChI_list$InChiKey) == 1) {this_tibble$INCHIKEY <- this_frag_data_fromSInChI_list$InChiKey}
      if (length(this_frag_data_fromSInChI_list$Formula) == 1) {this_tibble$FORMULA <- this_frag_data_fromSInChI_list$Formula}
      if (length(this_frag_data_fromSInChI_list$PMass) == 1) {this_tibble$PRECURSORMZ <- this_frag_data_fromSInChI_list$PMass}
      
      this_tibble$FRAGMENTS <- TRUE
      
      fragment_matrix <- getting_summary_fragment_matrix(all_else_from_CFMID = this_frag_data_fromSMILES_list$all_else, weight_energy0 = weight_energy0, weight_energy1 = weight_energy1, weight_energy2 = weight_energy2)
      
    }
    
    
    
    this_MSP_text <- paste0("NAME: ", this_tibble$NAME, "\n",
                            "PRECURSORMZ: ", this_tibble$PRECURSORMZ, "\n",
                            "PRECURSORTYPE: ", this_tibble$PRECURSORTYPE, "\n",
                            "IONMODE: ", this_tibble$IONMODE, "\n",
                            "RETENTIONTIME: ", this_tibble$RETENTIONTIME, "\n",
                            "CCS: ", this_tibble$CCS, "\n",
                            "FORMULA: ", this_tibble$FORMULA, "\n",
                            "ONTOLOGY: ", this_tibble$ONTOLOGY, "\n",
                            "SMILES: ", this_tibble$SMILES, "\n",
                            "INCHIKEY: ", this_tibble$INCHIKEY, "\n",
                            "INSTRUMENTTYPE: ", this_tibble$INSTRUMENTTYPE, "\n",
                            "COLLISIONENERGY: ", this_tibble$COLLISIONENERGY, "\n",
                            "COMMENT: ", ifelse(this_tibble$FRAGMENTS, paste0("Spectra generated by CFMID, weight energy0: ", weight_energy0,", energy1: ", weight_energy1 ,", energy2: ", weight_energy2), "NO in-silico fragment"), "; ID=", this_tibble$COMMENT_ID, "; CID=", this_tibble$COMMENT_CID, "; InChI=", this_tibble$COMMENT_InChI, "\n")
    
    
    if (grepl("\\n{2,}$", this_MSP_text)) {this_MSP_text <- gsub("\\n{2,}$", "\n", this_MSP_text)}
    
    if (this_tibble$FRAGMENTS) {
      
      this_MSP_text <- paste0(this_MSP_text,
                              "Num Peaks: ", ifelse(is.matrix(fragment_matrix), nrow(fragment_matrix), 1), "\n")
      
      
      if (is.matrix(fragment_matrix)) {
        for (f in 1:nrow(fragment_matrix)) {
          this_MSP_text <- paste0(this_MSP_text,
                                  fragment_matrix[f, 1], " ", fragment_matrix[f, 2], ifelse(f == nrow(fragment_matrix), "\n\n", "\n"))
        }
      } else {
        this_MSP_text <- paste0(this_MSP_text,
                                fragment_matrix[1], " ", fragment_matrix[2], "\n\n")
      }
      
      
    } else {
      this_MSP_text <- paste0(this_MSP_text,
                              "Num Peaks: 1\n",
                              this_tibble$PRECURSORMZ, " 100\n\n")
    }
                            
    if (grepl("\\n{3,}$", this_MSP_text)) {
      
      this_MSP_text <- gsub("\\n{3,}$", "\n\n", this_MSP_text)
      
    }
    
    
    final_MSP_text <- paste0(final_MSP_text, this_MSP_text)
    
  }
  
  return(final_MSP_text)
  
}



## the check:

checking_frag_output <- function(df_molecules, list_obtained) {
  
  if (any(duplicated(names(list_obtained)))) {
    stop("There are some duplicated in the names of the list")
  }
  
  for (a in names(list_obtained)) {
    if (length(list_obtained[[a]]) != 7) {stop("There are no exatly 7 element for each list")}
    if (a != list_obtained[[a]][[1]]) {stop("The first element is not the same as the name of the list")}
  }
  
  checking_table <- mutate(df_molecules,
                           present = NA)
  checking_table[,1] <- as.character(pull(checking_table, 1))
  
  not_present_character <- character()
  
  for (i in 1:length(checking_table$present)) {
    if (pull(checking_table, 1)[i] %in% names(list_obtained)) {
      checking_table[i, "present"] <- TRUE
    } else {
      checking_table[i, "present"] <- FALSE
      not_present_character <- c(not_present_character, pull(checking_table, 1)[i])
    }
  }
  
  
  if (length(not_present_character)>0) {
    cat(paste0("Not present are the following: ", paste0(not_present_character, collapse = "; "), "\n"))
  } else {
    cat("All present!! Yee!!!")
  }
  
  
  return(checking_table)
  
}


merge_list_parts <- function(list_of_all, if_duplicated_consider_from_higher = TRUE) {
  
  if (is.null(names(list_of_all))) {stop("put some names to the the list_of_all!")}
  
  all_molecules_df <- tibble(molecules = character(),
                             from_part = numeric(),
                             name_of_that_part = character())
  
  for (p in 1:length(list_of_all)) {
    this_list <- list_of_all[[p]]
    this_part_name <- names(list_of_all)[p]
    
    if (any(duplicated(names(this_list)))) {stop(paste0("There are some duplicated in the list ", this_part_name))}
    
    this_df <- tibble(molecules = names(this_list),
                      from_part = rep(p, length(names(this_list))),
                      name_of_that_part = rep(this_part_name, length(names(this_list))))
    
    
    all_molecules_df <- bind_rows(all_molecules_df, this_df)
    
  }
  
  
  if (any(duplicated(all_molecules_df$molecules))) {
    
    rows_to_remove <- numeric()
    for (d in unique(all_molecules_df$molecules[which(duplicated(all_molecules_df$molecules))])) {
      all_molecules_df_fil <- filter(all_molecules_df, molecules == d)
      
      if (if_duplicated_consider_from_higher) {
        part_to_keep <- max(all_molecules_df_fil$from_part)
      } else {
        part_to_keep <- min(all_molecules_df_fil$from_part)
      }
      
      for (sd in 1:length(all_molecules_df_fil$molecules)) {
        if (all_molecules_df_fil$from_part[sd] != part_to_keep) {
          this_row_to_remove <- which(all_molecules_df$molecules == all_molecules_df_fil$molecules[sd] & all_molecules_df$from_part == all_molecules_df_fil$from_part[sd])
          if (length(this_row_to_remove) !=1) {stop("something wrong...")}
          
          rows_to_remove <- c(rows_to_remove, this_row_to_remove)
          
          cat(paste0("Removed molecule ", d, " from part ", all_molecules_df_fil$name_of_that_part[sd], "\n"))
          
        } else {
          
          cat(paste0("Kept ", d, " from part ", all_molecules_df_fil$name_of_that_part[sd], "\n"))
          
        }
      }
      
    }
    
    all_molecules_df <- all_molecules_df[-rows_to_remove, ]
    
  }
  
  #just to check it worked:
  if (any(duplicated(all_molecules_df$molecules))) {stop("something wired...!")}
  
  
  list_output <- vector(mode = "list", length = length(all_molecules_df$molecules))
  names(list_output) <- all_molecules_df$molecules
  
  for (m in 1:length(all_molecules_df$molecules)) {
    list_output[[all_molecules_df$molecules[m]]] <- list_of_all[[all_molecules_df$from_part[m]]][[all_molecules_df$molecules[m]]]
  }
  
  return(list_output)
  
}



doing_everything_to_merge_the_list <- function(character_vector_with_output_part_file_names) {
  
  the_list_of_all <- vector(mode = "list", length = length(character_vector_with_output_part_file_names))
  names(the_list_of_all) <- str_remove_all(character_vector_with_output_part_file_names, ".txt")
  
  for (tx in names(the_list_of_all)) {
    
    this_lines <- readLines(paste0(tx, ".txt"))
    
    the_list_of_all[[tx]] <- creating_in_silico_fragmentation_list(this_lines)
    
  }
  
  the_merged_list <- merge_list_parts(the_list_of_all)
  
  return(the_merged_list)
}



##########


## POS:


WormJam_combined_SMILES_predicted_output_POS_list <- doing_everything_to_merge_the_list(c("WormJam_combined_SMILES_predicted_output_POS_part1of9.txt",
                                                                                          "WormJam_combined_SMILES_predicted_output_POS_part2of9.txt",
                                                                                          "WormJam_combined_SMILES_predicted_output_POS_part3of9.txt",
                                                                                          "WormJam_combined_SMILES_predicted_output_POS_part4of9_01.txt",
                                                                                          "WormJam_combined_SMILES_predicted_output_POS_part4of9_02.txt",
                                                                                          "WormJam_combined_SMILES_predicted_output_POS_part5of9.txt",
                                                                                          "WormJam_combined_SMILES_predicted_output_POS_part6of9.txt",
                                                                                          "WormJam_combined_SMILES_predicted_output_POS_part7of9.txt",
                                                                                          "WormJam_combined_SMILES_predicted_output_POS_part8of9.txt",
                                                                                          "WormJam_combined_SMILES_predicted_output_POS_part9of9.txt"))

# Removed molecule 118703056 from part WormJam_combined_SMILES_predicted_output_POS_part4of9_01
# Kept 118703056 from part WormJam_combined_SMILES_predicted_output_POS_part4of9_02
# Removed molecule 16738686 from part WormJam_combined_SMILES_predicted_output_POS_part8of9
# Kept 16738686 from part WormJam_combined_SMILES_predicted_output_POS_part9of9



output_POS_checking_tab <- checking_frag_output(df_molecules = WormJam_combined_SMILES,
                                                list_obtained = WormJam_combined_SMILES_predicted_output_POS_list)


## Not present are the following: 15939; 281; 444255; 121489512

# -> SO I WANTO TO RUN AGAIN CFM-ID ON THOSE:

WormJam_combined_SMILES_try_again_POS <- filter(WormJam_combined_SMILES, Identifier %in% c(15939, 281, 444255, 121489512))

write_tsv(WormJam_combined_SMILES_try_again_POS, "WormJam_combined_SMILES_try_again_POS.txt")


### So running it on those 4:

# docker run --rm=true -v ${pwd}:/cfmid/public/ -i wishartlab/cfmid:latest sh -c "cd /cfmid/public/; cfm-predict 'WormJam_combined_SMILES_try_again_POS.txt' 0.001 /trained_models_cfmid4.0/[M+H]+/param_output.log /trained_models_cfmid4.0/[M+H]+/param_config.txt 1 WormJam_combined_SMILES_predicted_output_POS_part_try_again.txt"

# could not ionize those 4!!




##

# NEG:

WormJam_combined_SMILES_predicted_output_NEG_list <- doing_everything_to_merge_the_list(c("WormJam_combined_SMILES_predicted_output_NEG_part1of9.txt",
                                                                                          "WormJam_combined_SMILES_predicted_output_NEG_part2of9.txt",
                                                                                          "WormJam_combined_SMILES_predicted_output_NEG_part3of9.txt",
                                                                                          "WormJam_combined_SMILES_predicted_output_NEG_part4of9.txt",
                                                                                          "WormJam_combined_SMILES_predicted_output_NEG_part5of9.txt",
                                                                                          "WormJam_combined_SMILES_predicted_output_NEG_part6of9.txt",
                                                                                          "WormJam_combined_SMILES_predicted_output_NEG_part7of9.txt",
                                                                                          "WormJam_combined_SMILES_predicted_output_NEG_part8of9.txt",
                                                                                          "WormJam_combined_SMILES_predicted_output_NEG_part9of9.txt"))

# Removed molecule 16738686 from part WormJam_combined_SMILES_predicted_output_NEG_part8of9
# Kept 16738686 from part WormJam_combined_SMILES_predicted_output_NEG_part9of9



output_NEG_checking_tab <- checking_frag_output(df_molecules = WormJam_combined_SMILES,
                                                list_obtained = WormJam_combined_SMILES_predicted_output_NEG_list)

# Not present are the following: 2; 133; 187; 248; 249; 280; 85; 305; 457; 134; 1014; 1130; 3380; 3855; 4510; 5886; 5893; 5943; 10918; 13597; 13805; 14181; 3624; 15939; 34756; 61106; 71921; 75834; 79076; 94191; 107739; 121992; 129658; 131151; 161234; 165491; 439415; 439756; 439924; 160339; 460603; 439285; 2733536; 5283588; 446962; 446961; 3037043; 5313082; 5283485; 5313476; 5313355; 5313164; 440121; 5283486; 167760; 1132; 5313387; 5313597; 23724625; 23724626; 23724627; 5313088; 5313092; 5313094; 5313098; 5313099; 5313107; 5313108; 5313110; 5313143; 449006; 5313173; 5313193; 5313196; 5313208; 5313214; 5313219; 5313225; 5313276; 5313279; 5313280; 5313282; 24778772; 5313337; 5313348; 5313366; 5313381; 5313384; 5313390; 5313426; 5313429; 5313479; 5313482; 5313486; 5313488; 5313489; 5313491; 5313526; 5313539; 5313605; 5313634; 24892761; 25203331; 97535; 44229245; 3006797; 52922221; 52922231; 52922449; 52922451; 52922453; 5313287; 52922464; 52922466; 52922656; 52922658; 52922686; 52922698; 52922730; 52922732; 52922740; 52922746; 52922748; 52922788; 52922790; 52922798; 52922804; 52922806; 52922808; 52922846; 52922856; 52922862; 52922864; 52922866; 52923078; 52923080; 52923254; 52923256; 52923258; 52923316; 52923376; 53477722; 53477802; 53478606; 53478608; 53478670; 53478680; 53478682; 53478702; 53478724; 53478726; 53478728; 53478732; 53478738; 53478740; 53478742; 53478744; 53478766; 53478768; 53478788; 53478790; 53478810; 53478812; 53478834; 53478836; 53478880; 53478882; 53478974; 53478976; 53478978; 53478980; 53478982; 53479014; 53479038; 53479094; 53479096; 53479492; 132849242; 135398701; 135445750; 136165274; 281; 168382; 157837; 213145; 439460; 192884; 160907; 3246939; 6912389; 5280613; 449186; 6604726; 46907934; 49787001; 6604733; 53481696; 70679041; 90659871; 135405047; 162393723; 162393725; 448301; 439830; 439773; 70679178; 70679210; 70679214; 70679244; 70679246; 118987308; 118987311; 118987314; 118987316; 118987318; 118987320; 118987322; 118987324; 118987326; 462; 555; 8263; 823; 68867; 86555; 123702; 3036899; 161930; 168380; 440568; 440832; 444675; 497300; 1131; 444255; 6426852; 3082079; 6426897; 6426902; 6426904; 10177003; 177510; 5313990; 16035470; 5313988; 22833597; 5313975; 5313977; 5313981; 5313983; 5314001; 9546668; 9546674; 52924052; 52924054; 53477824; 53480468; 53481654; 70679122; 54348826; 19838751; 71768170; 118796871; 121489512; 122164833; 131770061



# -> SO I WANTO TO RUN AGAIN CFM-ID ON THOSE:

WormJam_combined_SMILES_try_again_NEG <- filter(WormJam_combined_SMILES, Identifier %in% c(2, 133, 187, 248, 249, 280, 85, 305, 457, 134, 1014, 1130, 3380, 3855, 4510, 5886, 5893, 5943,
                                                                                           10918, 13597, 13805, 14181, 3624, 15939, 34756, 61106, 71921, 75834, 79076, 94191, 107739, 121992,
                                                                                           129658, 131151, 161234, 165491, 439415, 439756, 439924, 160339, 460603, 439285, 2733536, 5283588, 446962,
                                                                                           446961, 3037043, 5313082, 5283485, 5313476, 5313355, 5313164, 440121, 5283486, 167760, 1132, 5313387, 5313597,
                                                                                           23724625, 23724626, 23724627, 5313088, 5313092, 5313094, 5313098, 5313099, 5313107, 5313108, 5313110, 5313143,
                                                                                           449006, 5313173, 5313193, 5313196, 5313208, 5313214, 5313219, 5313225, 5313276, 5313279, 5313280, 5313282, 24778772,
                                                                                           5313337, 5313348, 5313366, 5313381, 5313384, 5313390, 5313426, 5313429, 5313479, 5313482, 5313486, 5313488, 5313489,
                                                                                           5313491, 5313526, 5313539, 5313605, 5313634, 24892761, 25203331, 97535, 44229245, 3006797, 52922221, 52922231, 52922449,
                                                                                           52922451, 52922453, 5313287, 52922464, 52922466, 52922656, 52922658, 52922686, 52922698, 52922730, 52922732, 52922740, 52922746,
                                                                                           52922748, 52922788, 52922790, 52922798, 52922804, 52922806, 52922808, 52922846, 52922856, 52922862, 52922864, 52922866, 52923078, 52923080,
                                                                                           52923254, 52923256, 52923258, 52923316, 52923376, 53477722, 53477802, 53478606, 53478608, 53478670, 53478680, 53478682, 53478702, 53478724,
                                                                                           53478726, 53478728, 53478732, 53478738, 53478740, 53478742, 53478744, 53478766, 53478768, 53478788, 53478790, 53478810, 53478812, 53478834,
                                                                                           53478836, 53478880, 53478882, 53478974, 53478976, 53478978, 53478980, 53478982, 53479014, 53479038, 53479094, 53479096, 53479492, 132849242,
                                                                                           135398701, 135445750, 136165274, 281, 168382, 157837, 213145, 439460, 192884, 160907, 3246939, 6912389, 5280613, 449186, 6604726, 46907934,
                                                                                           49787001, 6604733, 53481696, 70679041, 90659871, 135405047, 162393723, 162393725, 448301, 439830, 439773, 70679178, 70679210, 70679214,
                                                                                           70679244, 70679246, 118987308, 118987311, 118987314, 118987316, 118987318, 118987320, 118987322, 118987324, 118987326, 462, 555, 8263,
                                                                                           823, 68867, 86555, 123702, 3036899, 161930, 168380, 440568, 440832, 444675, 497300, 1131, 444255, 6426852, 3082079, 6426897, 6426902,
                                                                                           6426904, 10177003, 177510, 5313990, 16035470, 5313988, 22833597, 5313975, 5313977, 5313981, 5313983, 5314001, 9546668, 9546674, 52924052,
                                                                                           52924054, 53477824, 53480468, 53481654, 70679122, 54348826, 19838751, 71768170, 118796871, 121489512, 122164833, 131770061))
# length: 266

write_tsv(WormJam_combined_SMILES_try_again_NEG, "WormJam_combined_SMILES_try_again_NEG.txt")


### So running it on those 266:

# docker run --rm=true -v ${pwd}:/cfmid/public/ -i wishartlab/cfmid:latest sh -c "cd /cfmid/public/; cfm-predict 'WormJam_combined_SMILES_try_again_NEG.txt' 0.001 /trained_models_cfmid4.0/[M-H]-/param_output.log /trained_models_cfmid4.0/[M-H]-/param_config.txt 1 WormJam_combined_SMILES_predicted_output_NEG_part_try_again.txt"


# Still could not ionize!!



####################
## SO! Creating the MSP file!


# in POS the collision energy was 20, in NEG 30.

WormJam_combined_MSPtext_POS <- creating_MSP_text(compound_table = WormJam_combined, id_col = NULL, name_col = "Title", CID_col = "Identifier", INCHIKEY_col = "InChIKey", SMILES_col = "SMILES", InChI_col = "InChI", MOLECULAR_FORMULA_col = "MolecularFormula", NEUTRAL_MASS_col = "MonoisotopicMass", Ontology_col = NULL,
                                                  polarity = "POS", CFMID_output_from_SMILES_list = WormJam_combined_SMILES_predicted_output_POS_list, CFMID_output_from_InChI_list = NULL, ID_used = "Identifier",
                                                  weight_energy0 = 0.1, weight_energy1 = 0.8, weight_energy2 = 0.1)
write(WormJam_combined_MSPtext_POS, file = "WormJam_combined_MSPtext_POS.msp")



WormJam_combined_MSPtext_NEG <- creating_MSP_text(compound_table = WormJam_combined, id_col = NULL, name_col = "Title", CID_col = "Identifier", INCHIKEY_col = "InChIKey", SMILES_col = "SMILES", InChI_col = "InChI", MOLECULAR_FORMULA_col = "MolecularFormula", NEUTRAL_MASS_col = "MonoisotopicMass", Ontology_col = NULL,
                                                  polarity = "NEG", CFMID_output_from_SMILES_list = WormJam_combined_SMILES_predicted_output_NEG_list, CFMID_output_from_InChI_list = NULL, ID_used = "Identifier",
                                                  weight_energy0 = 0.1, weight_energy1 = 0.45, weight_energy2 = 0.45)
write(WormJam_combined_MSPtext_NEG, file = "WormJam_combined_MSPtext_NEG.msp")



