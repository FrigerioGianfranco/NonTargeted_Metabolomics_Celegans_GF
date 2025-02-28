

library(tidyverse)
library(eulerr)
library(agricolae) #for the Fisher LSD
set.seed(1991)



#### function to export a vector with only features with no NAs:

get_noNAs_features <- function(dataset) {
  
  return_vector <- character()
  
  for (i in 1:length(pull(dataset, 1))) {
    
    if (!all(is.na(as.numeric(dataset[i, 2:length(colnames(dataset))])))) {
      return_vector <- c(return_vector, pull(dataset, 1)[i])
    }
  }
  
  cat(paste0(length(return_vector), " out of ", length(pull(dataset, 1)), " where kept.\n"))
  
  return(return_vector)
}


sample_info <- read_tsv("samples_shipped_GF_mod.txt")
sample_info$Scheme <- factor(sample_info$Scheme, levels = unique(sample_info$Scheme))
sample_info$Solvent <- factor(sample_info$Solvent, levels = unique(sample_info$Solvent))
sample_info$Schemesolvent <- factor(sample_info$Schemesolvent, levels = unique(sample_info$Schemesolvent))
sample_info$rawdried <- factor(sample_info$rawdried, levels = unique(sample_info$rawdried))
sample_info$Sample_ID <- factor(sample_info$Sample_ID, levels = unique(sample_info$Sample_ID))
sample_info$Content <- factor(sample_info$Content, levels = unique(sample_info$Content))
sample_info$Collect_batch <- factor(sample_info$Collect_batch, levels = unique(sample_info$Collect_batch))
sample_info$Solvent_num <- factor(sample_info$Solvent_num, levels = unique(sample_info$Solvent_num))
sample_info$Comments <- factor(sample_info$Comments, levels = unique(sample_info$Comments))

sample_info_dry <- filter(sample_info, rawdried == "DRIED")
dry_samples_codes <- sample_info_dry$Analyses_name







#######

HLPnrf <- read_tsv("C_elegans_HILIC_POS_Area_0_2025_01_21_18_49_19.txt")
HLP_feat_info <- read_tsv("C_elegans_HILIC_POS_Area_featINFO_0_2025_01_21_18_49_19.txt")


## adding an "X" at the beginning of each sample that starts with a number and replacing parenthesis with dots

colnames(HLPnrf)[grepl("^[[:digit:]]+", colnames(HLPnrf))] <- paste0("X", colnames(HLPnrf)[grepl("^[[:digit:]]+", colnames(HLPnrf))])
colnames(HLPnrf) <- str_replace_all(colnames(HLPnrf), "[//(//)]", ".")
colnames(HLPnrf) <- str_replace_all(colnames(HLPnrf), "[//)//)]", ".")
colnames(HLPnrf) <- str_replace_all(colnames(HLPnrf), "[// //)]", "_")

colnames(HLP_feat_info) <- str_replace_all(colnames(HLP_feat_info), "[// //)]", "_")
colnames(HLP_feat_info) <- str_replace_all(colnames(HLP_feat_info), "[/////)]", "")
colnames(HLP_feat_info) <- str_replace_all(colnames(HLP_feat_info), "[//(//)]", ".")
colnames(HLP_feat_info) <- str_replace_all(colnames(HLP_feat_info), "[//)//)]", ".")

## replacing zero with NA

HLPnrf[HLPnrf == 0] <- NA
HLPnrf[1,1] <- 0

## adding a "feat" at the beginning of feature number

add_feat <- function(x) {
  
  x_final <- as.character(rep(NA, length(x)))
  
  for (i in 1:length(x)) {
    zeros <- if (x[i] > 1000000) {
      stop("too many features!! XD")
    } else if (x[i] >= 100000) {
      ""
    } else if (x[i] >= 10000) {
      "0"
    } else if (x[i] >= 1000) {
      "00"
    } else if (x[i] >= 100) {
      "000"
    } else if (x[i] >= 10) {
      "0000"
    } else {
      "00000"
    }
    
    x_final[i] <- paste0("feat", zeros, x[i])
    
  }
  
  return(x_final)
}

HLPnrf$Alignment_ID <- add_feat(HLPnrf$Alignment_ID)
HLP_feat_info$Alignment_ID <- add_feat(HLP_feat_info$Alignment_ID)


all_featuresHLPnrf <- get_noNAs_features(HLPnrf)


analyses_dry_samples <- colnames(HLPnrf)[which(grepl(paste(dry_samples_codes, collapse = "|"), colnames(HLPnrf)))]


HLPnf <- select(HLPnrf, all_of(c("Alignment_ID", analyses_dry_samples)))

all_featuresHLPnf <- get_noNAs_features(HLPnf)

HLP <- filter(HLPnf, Alignment_ID %in% all_featuresHLPnf)

write_tsv(HLP, "HLP.txt")





## removing the first analyses of the run, for each category of samples, and creating the feature numbers

HLP0nf <- select(HLP,
                 
                 -X136_Pool_N2_ACN_pool,
                 -X137_Pool_VC40_ACN_pool,
                 -X138_Pool_VC1668_ACN_pool,
                 -X139_Pool_BR5270_ACN_pool,
                 -X140_Pool_UA57_ACN_pool,
                 -X143_Pool_VC40_ACN_pool,
                 
                 -X424_.dry._Pool_N2_2B_pool,
                 -X425_.dry._Pool_VC40_2B_pool,
                 -X426_.dry._Pool_VC1668_2B_pool,
                 -X427_.dry._Pool_BR5270_2B_pool,
                 -X428_.dry._Pool_UA57_2B_pool,
                 -X431_.dry._Pool_VC40_2B_pool)




all_featuresHLP0nf <- get_noNAs_features(HLP0nf)

HLP0 <- filter(HLP0nf, Alignment_ID %in% all_featuresHLP0nf)

write_tsv(HLP0, "HLP0.txt")



HLP_feat_info <- filter(HLP_feat_info, Alignment_ID %in% all_featuresHLP0nf)
write_tsv(HLP_feat_info, "HLP_feat_info.txt")








## Doing all the separated QC cleaning



## the dry samples


# dry_N2_ACN

HLP0_dry_N2_ACN <- HLP0 %>%
  mutate(X083_.dry._N2_blank_ACN_1 = ifelse(is.na(X083_.dry._N2_blank_ACN_1), 0, X083_.dry._N2_blank_ACN_1),
         X084_.dry._N2_blank_ACN_2 = ifelse(is.na(X084_.dry._N2_blank_ACN_2), 0, X084_.dry._N2_blank_ACN_2)) %>%
  rowwise() %>%
  mutate(MEDIA = mean(c(X136_Pool_N2_ACN_pool_20220714203441,
                        X136_Pool_N2_ACN_pool_20220715004714,
                        X136_Pool_N2_ACN_pool_20220715072327), na.rm=TRUE),
         SD = sd(c(X136_Pool_N2_ACN_pool_20220714203441,
                   X136_Pool_N2_ACN_pool_20220715004714,
                   X136_Pool_N2_ACN_pool_20220715072327), na.rm=TRUE),
         RSD = SD/MEDIA*100,
         COUNT = mean(!is.na(c(X136_Pool_N2_ACN_pool_20220714203441,
                               X136_Pool_N2_ACN_pool_20220715004714,
                               X136_Pool_N2_ACN_pool_20220715072327)))*100,
         BLANK_CONTR = mean(mean(X083_.dry._N2_blank_ACN_1, X084_.dry._N2_blank_ACN_2)/MEDIA*100)
  )

HLP0_dry_N2_ACNf <- HLP0_dry_N2_ACN %>%
  filter(COUNT>50) %>%
  filter(RSD<50) %>%
  filter(BLANK_CONTR<50)




# dry_VC40_ACN

HLP0_dry_VC40_ACN <- HLP0 %>%
  mutate(X089_.dry._VC40_blank_ACN_1 = ifelse(is.na(X089_.dry._VC40_blank_ACN_1), 0, X089_.dry._VC40_blank_ACN_1),
         X090_.dry._VC40_blank_ACN_2 = ifelse(is.na(X090_.dry._VC40_blank_ACN_2), 0, X090_.dry._VC40_blank_ACN_2)) %>%
  rowwise() %>%
  mutate(MEDIA = mean(c(X137_Pool_VC40_ACN_pool_20220714205244,
                        X137_Pool_VC40_ACN_pool_20220715010516,
                        X137_Pool_VC40_ACN_pool_20220715074131), na.rm=TRUE),
         SD = sd(c(X137_Pool_VC40_ACN_pool_20220714205244,
                   X137_Pool_VC40_ACN_pool_20220715010516,
                   X137_Pool_VC40_ACN_pool_20220715074131), na.rm=TRUE),
         RSD = SD/MEDIA*100,
         COUNT = mean(!is.na(c(X137_Pool_VC40_ACN_pool_20220714205244,
                               X137_Pool_VC40_ACN_pool_20220715010516,
                               X137_Pool_VC40_ACN_pool_20220715074131)))*100,
         BLANK_CONTR = mean(mean(X089_.dry._VC40_blank_ACN_1, X090_.dry._VC40_blank_ACN_2)/MEDIA*100)
  )

HLP0_dry_VC40_ACNf <- HLP0_dry_VC40_ACN %>%
  filter(COUNT>50) %>%
  filter(RSD<50) %>%
  filter(BLANK_CONTR<50)



# dry_VC1668_ACN

HLP0_dry_VC1668_ACN <- HLP0 %>%
  mutate(X095_.dry._VC1668_blank_ACN_1 = ifelse(is.na(X095_.dry._VC1668_blank_ACN_1), 0, X095_.dry._VC1668_blank_ACN_1),
         X096_.dry._VC1668_blank_ACN_2 = ifelse(is.na(X096_.dry._VC1668_blank_ACN_2), 0, X096_.dry._VC1668_blank_ACN_2)) %>%
  rowwise() %>%
  mutate(MEDIA = mean(c(X138_Pool_VC1668_ACN_pool_20220714211048,
                        X138_Pool_VC1668_ACN_pool_20220715012316,
                        X138_Pool_VC1668_ACN_pool_20220715075934), na.rm=TRUE),
         SD = sd(c(X138_Pool_VC1668_ACN_pool_20220714211048,
                   X138_Pool_VC1668_ACN_pool_20220715012316,
                   X138_Pool_VC1668_ACN_pool_20220715075934), na.rm=TRUE),
         RSD = SD/MEDIA*100,
         COUNT = mean(!is.na(c(X138_Pool_VC1668_ACN_pool_20220714211048,
                               X138_Pool_VC1668_ACN_pool_20220715012316,
                               X138_Pool_VC1668_ACN_pool_20220715075934)))*100,
         BLANK_CONTR = mean(mean(X095_.dry._VC1668_blank_ACN_1, X096_.dry._VC1668_blank_ACN_2)/MEDIA*100)
  )

HLP0_dry_VC1668_ACNf <- HLP0_dry_VC1668_ACN %>%
  filter(COUNT>50) %>%
  filter(RSD<50) %>%
  filter(BLANK_CONTR<50)



# dry_BR5270_ACN

HLP0_dry_BR5270_ACN <- HLP0 %>%
  mutate(X101_.dry._BR5270_blank_ACN_1 = ifelse(is.na(X101_.dry._BR5270_blank_ACN_1), 0, X101_.dry._BR5270_blank_ACN_1),
         X102_.dry._BR5270_blank_ACN_2 = ifelse(is.na(X102_.dry._BR5270_blank_ACN_2), 0, X102_.dry._BR5270_blank_ACN_2)) %>%
  rowwise() %>%
  mutate(MEDIA = mean(c(X139_Pool_BR5270_ACN_pool_20220714212851,
                        X139_Pool_BR5270_ACN_pool_20220715014117,
                        X139_Pool_BR5270_ACN_pool_20220715081739), na.rm=TRUE),
         SD = sd(c(X139_Pool_BR5270_ACN_pool_20220714212851,
                   X139_Pool_BR5270_ACN_pool_20220715014117,
                   X139_Pool_BR5270_ACN_pool_20220715081739), na.rm=TRUE),
         RSD = SD/MEDIA*100,
         COUNT = mean(!is.na(c(X139_Pool_BR5270_ACN_pool_20220714212851,
                               X139_Pool_BR5270_ACN_pool_20220715014117,
                               X139_Pool_BR5270_ACN_pool_20220715081739)))*100,
         BLANK_CONTR = mean(mean(X101_.dry._BR5270_blank_ACN_1, X102_.dry._BR5270_blank_ACN_2)/MEDIA*100)
  )

HLP0_dry_BR5270_ACNf <- HLP0_dry_BR5270_ACN %>%
  filter(COUNT>50) %>%
  filter(RSD<50) %>%
  filter(BLANK_CONTR<50)



# dry_UA57_ACN

HLP0_dry_UA57_ACN <- HLP0 %>%
  mutate(X107_.dry._UA57_blank_ACN_1 = ifelse(is.na(X107_.dry._UA57_blank_ACN_1), 0, X107_.dry._UA57_blank_ACN_1),
         X108_.dry._UA57_blank_ACN_2 = ifelse(is.na(X108_.dry._UA57_blank_ACN_2), 0, X108_.dry._UA57_blank_ACN_2)) %>%
  rowwise() %>%
  mutate(MEDIA = mean(c(X140_Pool_UA57_ACN_pool_20220714214654,
                        X140_Pool_UA57_ACN_pool_20220715015918,
                        X140_Pool_UA57_ACN_pool_20220715083543), na.rm=TRUE),
         SD = sd(c(X140_Pool_UA57_ACN_pool_20220714214654,
                   X140_Pool_UA57_ACN_pool_20220715015918,
                   X140_Pool_UA57_ACN_pool_20220715083543), na.rm=TRUE),
         RSD = SD/MEDIA*100,
         COUNT = mean(!is.na(c(X140_Pool_UA57_ACN_pool_20220714214654,
                               X140_Pool_UA57_ACN_pool_20220715015918,
                               X140_Pool_UA57_ACN_pool_20220715083543)))*100,
         BLANK_CONTR = mean(mean(X107_.dry._UA57_blank_ACN_1, X108_.dry._UA57_blank_ACN_2)/MEDIA*100)
  )

HLP0_dry_UA57_ACNf <- HLP0_dry_UA57_ACN %>%
  filter(COUNT>50) %>%
  filter(RSD<50) %>%
  filter(BLANK_CONTR<50)



# dry_VC40_ACN_bis


HLP0_dry_VC40_ACN_bis <- HLP0 %>%
  mutate(X125_.dry._VC40_blank_ACN_1 = ifelse(is.na(X125_.dry._VC40_blank_ACN_1), 0, X125_.dry._VC40_blank_ACN_1),
         X126_.dry._VC40_blank_ACN_2 = ifelse(is.na(X126_.dry._VC40_blank_ACN_2), 0, X126_.dry._VC40_blank_ACN_2)) %>%
  rowwise() %>%
  mutate(MEDIA = mean(c(X143_Pool_VC40_ACN_pool_20220714220458,
                        X143_Pool_VC40_ACN_pool_20220715021719,
                        X143_Pool_VC40_ACN_pool_20220715085345), na.rm=TRUE),
         SD = sd(c(X143_Pool_VC40_ACN_pool_20220714220458,
                   X143_Pool_VC40_ACN_pool_20220715021719,
                   X143_Pool_VC40_ACN_pool_20220715085345), na.rm=TRUE),
         RSD = SD/MEDIA*100,
         COUNT = mean(!is.na(c(X143_Pool_VC40_ACN_pool_20220714220458,
                               X143_Pool_VC40_ACN_pool_20220715021719,
                               X143_Pool_VC40_ACN_pool_20220715085345)))*100,
         BLANK_CONTR = mean(mean(X125_.dry._VC40_blank_ACN_1, X126_.dry._VC40_blank_ACN_2)/MEDIA*100)
  )

HLP0_dry_VC40_ACN_bisf <- HLP0_dry_VC40_ACN_bis %>%
  filter(COUNT>50) %>%
  filter(RSD<50) %>%
  filter(BLANK_CONTR<50)









# dry_N2_2B

HLP0_dry_N2_2B <- HLP0 %>%
  mutate(X371_.dry._N2_blank_2B_1 = ifelse(is.na(X371_.dry._N2_blank_2B_1), 0, X371_.dry._N2_blank_2B_1),
         X372_.dry._N2_blank_2B_2 = ifelse(is.na(X372_.dry._N2_blank_2B_2), 0, X372_.dry._N2_blank_2B_2)) %>%
  rowwise() %>%
  mutate(MEDIA = mean(c(X424_.dry._Pool_N2_2B_pool_20220716105522,
                        X424_.dry._Pool_N2_2B_pool_20220716150725,
                        X424_.dry._Pool_N2_2B_pool_20220716214339), na.rm=TRUE),
         SD = sd(c(X424_.dry._Pool_N2_2B_pool_20220716105522,
                   X424_.dry._Pool_N2_2B_pool_20220716150725,
                   X424_.dry._Pool_N2_2B_pool_20220716214339), na.rm=TRUE),
         RSD = SD/MEDIA*100,
         COUNT = mean(!is.na(c(X424_.dry._Pool_N2_2B_pool_20220716105522,
                               X424_.dry._Pool_N2_2B_pool_20220716150725,
                               X424_.dry._Pool_N2_2B_pool_20220716214339)))*100,
         BLANK_CONTR = mean(mean(X371_.dry._N2_blank_2B_1, X372_.dry._N2_blank_2B_2)/MEDIA*100)
  )

HLP0_dry_N2_2Bf <- HLP0_dry_N2_2B %>%
  filter(COUNT>50) %>%
  filter(RSD<50) %>%
  filter(BLANK_CONTR<50)



# dry_VC40_2B

HLP0_dry_VC40_2B <- HLP0 %>%
  mutate(X377_.dry._VC40_blank_2B_1 = ifelse(is.na(X377_.dry._VC40_blank_2B_1), 0, X377_.dry._VC40_blank_2B_1),
         X378_.dry._VC40_blank_2B_2 = ifelse(is.na(X378_.dry._VC40_blank_2B_2), 0, X378_.dry._VC40_blank_2B_2)) %>%
  rowwise() %>%
  mutate(MEDIA = mean(c(X425_.dry._Pool_VC40_2B_pool_20220716111323,
                        X425_.dry._Pool_VC40_2B_pool_20220716152528,
                        X425_.dry._Pool_VC40_2B_pool_20220716220141), na.rm=TRUE),
         SD = sd(c(X425_.dry._Pool_VC40_2B_pool_20220716111323,
                   X425_.dry._Pool_VC40_2B_pool_20220716152528,
                   X425_.dry._Pool_VC40_2B_pool_20220716220141), na.rm=TRUE),
         RSD = SD/MEDIA*100,
         COUNT = mean(!is.na(c(X425_.dry._Pool_VC40_2B_pool_20220716111323,
                               X425_.dry._Pool_VC40_2B_pool_20220716152528,
                               X425_.dry._Pool_VC40_2B_pool_20220716220141)))*100,
         BLANK_CONTR = mean(mean(X377_.dry._VC40_blank_2B_1, X378_.dry._VC40_blank_2B_2)/MEDIA*100)
  )

HLP0_dry_VC40_2Bf <- HLP0_dry_VC40_2B %>%
  filter(COUNT>50) %>%
  filter(RSD<50) %>%
  filter(BLANK_CONTR<50)



# dry_VC1668_2B

HLP0_dry_VC1668_2B <- HLP0 %>%
  mutate(X383_.dry._VC1668_blank_2B_1 = ifelse(is.na(X383_.dry._VC1668_blank_2B_1), 0, X383_.dry._VC1668_blank_2B_1),
         X384_.dry._VC1668_blank_2B_2 = ifelse(is.na(X384_.dry._VC1668_blank_2B_2), 0, X384_.dry._VC1668_blank_2B_2)) %>%
  rowwise() %>%
  mutate(MEDIA = mean(c(X426_.dry._Pool_VC1668_2B_pool_20220716113126,
                        X426_.dry._Pool_VC1668_2B_pool_20220716154332,
                        X426_.dry._Pool_VC1668_2B_pool_20220716221942), na.rm=TRUE),
         SD = sd(c(X426_.dry._Pool_VC1668_2B_pool_20220716113126,
                   X426_.dry._Pool_VC1668_2B_pool_20220716154332,
                   X426_.dry._Pool_VC1668_2B_pool_20220716221942), na.rm=TRUE),
         RSD = SD/MEDIA*100,
         COUNT = mean(!is.na(c(X426_.dry._Pool_VC1668_2B_pool_20220716113126,
                               X426_.dry._Pool_VC1668_2B_pool_20220716154332,
                               X426_.dry._Pool_VC1668_2B_pool_20220716221942)))*100,
         BLANK_CONTR = mean(mean(X383_.dry._VC1668_blank_2B_1, X384_.dry._VC1668_blank_2B_2)/MEDIA*100)
  )

HLP0_dry_VC1668_2Bf <- HLP0_dry_VC1668_2B %>%
  filter(COUNT>50) %>%
  filter(RSD<50) %>%
  filter(BLANK_CONTR<50)



# dry_BR5270_2B

HLP0_dry_BR5270_2B <- HLP0 %>%
  mutate(X389_.dry._BR5270_blank_2B_1 = ifelse(is.na(X389_.dry._BR5270_blank_2B_1), 0, X389_.dry._BR5270_blank_2B_1),
         X390_.dry._BR5270_blank_2B_2 = ifelse(is.na(X390_.dry._BR5270_blank_2B_2), 0, X390_.dry._BR5270_blank_2B_2)) %>%
  rowwise() %>%
  mutate(MEDIA = mean(c(X427_.dry._Pool_BR5270_2B_pool_20220716114922,
                        X427_.dry._Pool_BR5270_2B_pool_20220716160136,
                        X427_.dry._Pool_BR5270_2B_pool_20220716223744), na.rm=TRUE),
         SD = sd(c(X427_.dry._Pool_BR5270_2B_pool_20220716114922,
                   X427_.dry._Pool_BR5270_2B_pool_20220716160136,
                   X427_.dry._Pool_BR5270_2B_pool_20220716223744), na.rm=TRUE),
         RSD = SD/MEDIA*100,
         COUNT = mean(!is.na(c(X427_.dry._Pool_BR5270_2B_pool_20220716114922,
                               X427_.dry._Pool_BR5270_2B_pool_20220716160136,
                               X427_.dry._Pool_BR5270_2B_pool_20220716223744)))*100,
         BLANK_CONTR = mean(mean(X389_.dry._BR5270_blank_2B_1, X390_.dry._BR5270_blank_2B_2)/MEDIA*100)
  )

HLP0_dry_BR5270_2Bf <- HLP0_dry_BR5270_2B %>%
  filter(COUNT>50) %>%
  filter(RSD<50) %>%
  filter(BLANK_CONTR<50)



# dry_UA57_2B

HLP0_dry_UA57_2B <- HLP0 %>%
  mutate(X395_.dry._UA57_blank_2B_1 = ifelse(is.na(X395_.dry._UA57_blank_2B_1), 0, X395_.dry._UA57_blank_2B_1),
         X396_.dry._UA57_blank_2B_2 = ifelse(is.na(X396_.dry._UA57_blank_2B_2), 0, X396_.dry._UA57_blank_2B_2)) %>%
  rowwise() %>%
  mutate(MEDIA = mean(c(X428_.dry._Pool_UA57_2B_pool_20220716120723,
                        X428_.dry._Pool_UA57_2B_pool_20220716161939,
                        X428_.dry._Pool_UA57_2B_pool_20220716225547), na.rm=TRUE),
         SD = sd(c(X428_.dry._Pool_UA57_2B_pool_20220716120723,
                   X428_.dry._Pool_UA57_2B_pool_20220716161939,
                   X428_.dry._Pool_UA57_2B_pool_20220716225547), na.rm=TRUE),
         RSD = SD/MEDIA*100,
         COUNT = mean(!is.na(c(X428_.dry._Pool_UA57_2B_pool_20220716120723,
                               X428_.dry._Pool_UA57_2B_pool_20220716161939,
                               X428_.dry._Pool_UA57_2B_pool_20220716225547)))*100,
         BLANK_CONTR = mean(mean(X395_.dry._UA57_blank_2B_1, X396_.dry._UA57_blank_2B_2)/MEDIA*100)
  )

HLP0_dry_UA57_2Bf <- HLP0_dry_UA57_2B %>%
  filter(COUNT>50) %>%
  filter(RSD<50) %>%
  filter(BLANK_CONTR<50)



# dry_VC40_2B_bis

HLP0_dry_VC40_2B_bis <- HLP0 %>%
  mutate(X413_.dry._VC40_blank_2B_1 = ifelse(is.na(X413_.dry._VC40_blank_2B_1), 0, X413_.dry._VC40_blank_2B_1),
         X414_.dry._VC40_blank_2B_2 = ifelse(is.na(X414_.dry._VC40_blank_2B_2), 0, X414_.dry._VC40_blank_2B_2)) %>%
  rowwise() %>%
  mutate(MEDIA = mean(c(X431_.dry._Pool_VC40_2B_pool_20220716122523,
                        X431_.dry._Pool_VC40_2B_pool_20220716163744,
                        X431_.dry._Pool_VC40_2B_pool_20220716231347), na.rm=TRUE),
         SD = sd(c(X431_.dry._Pool_VC40_2B_pool_20220716122523,
                   X431_.dry._Pool_VC40_2B_pool_20220716163744,
                   X431_.dry._Pool_VC40_2B_pool_20220716231347), na.rm=TRUE),
         RSD = SD/MEDIA*100,
         COUNT = mean(!is.na(c(X431_.dry._Pool_VC40_2B_pool_20220716122523,
                               X431_.dry._Pool_VC40_2B_pool_20220716163744,
                               X431_.dry._Pool_VC40_2B_pool_20220716231347)))*100,
         BLANK_CONTR = mean(mean(X413_.dry._VC40_blank_2B_1, X414_.dry._VC40_blank_2B_2)/MEDIA*100)
  )

HLP0_dry_VC40_2B_bisf <- HLP0_dry_VC40_2B_bis %>%
  filter(COUNT>50) %>%
  filter(RSD<50) %>%
  filter(BLANK_CONTR<50)





### merge!!

HLP_unique_feat_vector <- unique(c(HLP0_dry_N2_ACNf$Alignment_ID, HLP0_dry_VC40_ACNf$Alignment_ID, HLP0_dry_VC1668_ACNf$Alignment_ID, HLP0_dry_BR5270_ACNf$Alignment_ID, HLP0_dry_UA57_ACNf$Alignment_ID, HLP0_dry_VC40_ACN_bisf$Alignment_ID,
                                   HLP0_dry_N2_2Bf$Alignment_ID, HLP0_dry_VC40_2Bf$Alignment_ID, HLP0_dry_VC1668_2Bf$Alignment_ID, HLP0_dry_BR5270_2Bf$Alignment_ID, HLP0_dry_UA57_2Bf$Alignment_ID, HLP0_dry_VC40_2B_bisf$Alignment_ID))


HLP1 <- filter(HLP0, Alignment_ID %in% HLP_unique_feat_vector)

write_tsv(HLP1, "HLP1.txt")


## some visualizations of the results
### to me it seems more convenient to put samples in rows, and features in columns!!

HLP1_t <- HLP1 %>%
  t()

HLP1_t_first_col <- matrix(data = rownames(HLP1_t), nrow = nrow(HLP1_t), ncol = 1)
HLP1_t <- cbind(HLP1_t_first_col, HLP1_t)
colnames(HLP1_t) <- c("ANALYSIS", as.vector(HLP1_t[1,])[-1])

HLP1_tc <- HLP1_t[-1,]
rownames(HLP1_tc) <- NULL

HLP2 <- as_tibble(HLP1_tc) %>%
  mutate_at(colnames(HLP1_tc)[-1], as.numeric)

write_tsv(HLP2, "HLP2.txt")


#### 




HLP3 <- HLP2 %>%
  mutate(CORRESP_FOUND = ifelse(strsplit(ANALYSIS, split = "_")[[1]][1] %in% sample_info$Analyses_name, TRUE, FALSE),
         .after = ANALYSIS)




if (mean(HLP3$CORRESP_FOUND) == 1) {
  for (a in colnames(sample_info)) {
    HLP3 <- mutate(HLP3, NEWCOLOUMN = rep(NA, length(HLP3$ANALYSIS)), .before = CORRESP_FOUND)
    colnames(HLP3)[colnames(HLP3) == "NEWCOLOUMN"] <- a
    for (ra in 1:length(HLP3$ANALYSIS)) {
      for (rs in 1:length(sample_info$Analyses_name)) {
        if (sample_info$Analyses_name[rs] == strsplit(HLP3$ANALYSIS, split = "_")[[ra]][1]) {
          HLP3[ra, a] <- sample_info[rs, a]
        }
      }
    }
  }
} else {
  stop("There are some analyses of which name is not included in the list of shipped samples!")
}

write_tsv(HLP3, "HLP3.txt")



## a check
#samples not included in the sequence of analyses:
samples_not_in_the_HLP_sequence <- filter(sample_info_dry, Analyses_name %in% sample_info_dry$Analyses_name[!sample_info_dry$Analyses_name %in% HLP3$Analyses_name])




### then, I would do Eulero-venn

                            
HLP_unique_feat_vector_dry_schem1 <- unique(c(HLP0_dry_N2_ACNf$Alignment_ID, HLP0_dry_VC40_ACNf$Alignment_ID, HLP0_dry_VC1668_ACNf$Alignment_ID, HLP0_dry_BR5270_ACNf$Alignment_ID, HLP0_dry_UA57_ACNf$Alignment_ID, HLP0_dry_VC40_ACN_bisf$Alignment_ID))

HLP_unique_feat_vector_dry_schem2 <- unique(c(HLP0_dry_N2_2Bf$Alignment_ID, HLP0_dry_VC40_2Bf$Alignment_ID, HLP0_dry_VC1668_2Bf$Alignment_ID, HLP0_dry_BR5270_2Bf$Alignment_ID, HLP0_dry_UA57_2Bf$Alignment_ID, HLP0_dry_VC40_2B_bisf$Alignment_ID))


#Eulero-venn combined:

HLP_EV_df <- tibble(feature = HLP_unique_feat_vector,
                    dry_scheme1 = HLP_unique_feat_vector %in% HLP_unique_feat_vector_dry_schem1,
                    dry_scheme2bottom = HLP_unique_feat_vector %in% HLP_unique_feat_vector_dry_schem2)

# better:

HLP_EV_df_better <- HLP_EV_df
colnames(HLP_EV_df_better) <- c("feature", "scheme1", "scheme2")


jpeg(file="HILIC_POS_EuleroVenn.jpeg")
plot(euler(HLP_EV_df_better[, 2:3], shape = "ellipse"), fill = c("cadetblue1", "indianred1"), quantities = list(cex = 2.5), labels = list(cex = 2))
dev.off()


#
###
#################### tables for the ANOVA::

HLP4 <- filter(HLP3, Sample_ID != "methBK")

HLP4$Sample_ID <- droplevels(HLP4$Sample_ID)


HLP4_shem1_dry <- filter(HLP4, Scheme == "scheme1", rawdried == "DRIED", Content == "worm") %>%
  #keeping only: ANALYSIS and Sample_ID
  select(-Scheme,             
         -Analyses_name,
         -Tube_num,
         -Label_top,
         -Label_side,
         -Solvent,
         -Schemesolvent,
         -rawdried,
         -Content,
         -Collect_batch,
         -Solvent_num,
         -Comments,
         -Side_reLABELING_white_tape,
         -CORRESP_FOUND)


HLP4_shem2_dry <- filter(HLP4, Scheme == "scheme2", rawdried == "DRIED", Content == "worm") %>%
  #keeping only: ANALYSIS and Sample_ID
  select(-Scheme,             
         -Analyses_name,
         -Tube_num,
         -Label_top,
         -Label_side,
         -Solvent,
         -Schemesolvent,
         -rawdried,
         -Content,
         -Collect_batch,
         -Solvent_num,
         -Comments,
         -Side_reLABELING_white_tape,
         -CORRESP_FOUND)

###


### function to transform the feature intensities:

transf_feat_intensities <- function(feat_int) {
  
  DfStAn <- feat_int
  
  
  pareto_scale <- function(x, na.rm = TRUE) {(x - mean(x, na.rm = na.rm)) / sqrt(sd(x, na.rm))}
  
  ## missing value
  DfStAn_mr <- DfStAn
  
  
  for (a in colnames(DfStAn)[3:length(colnames(DfStAn))]) {
    if (mean(is.na(pull(DfStAn, a))) == 1) {
      DfStAn_mr <- DfStAn_mr[, colnames(DfStAn_mr)[-which(colnames(DfStAn_mr)==a)]]
    } else {
      for (i in 1:length(pull(DfStAn, a))) {
        if (is.na(pull(DfStAn, a)[i])) {
          DfStAn_mr[i,a] <-  min(pull(DfStAn, a), na.rm = TRUE)/5
        }
      }
    }
  }
  
  
  # transf
  
  DfStAn_mrn <- DfStAn_mr %>%
    mutate_at(colnames(DfStAn_mr)[3:length(colnames(DfStAn_mr))],
              log10) %>%
    mutate_at(colnames(DfStAn_mr)[3:length(colnames(DfStAn_mr))],
              pareto_scale)
  
  return(DfStAn_mrn)
  
}


### function to obtain the ANOVA table:

getting_ANOVA_results <- function(feat_int_transf) {
  
  
  DfStAn_mrn <- feat_int_transf
  
  
  
  ANOVA_results <- tibble(feature = as.character(colnames(DfStAn_mrn)[3:length(colnames(DfStAn_mrn))]),
                          ANOVA_p_value = as.numeric(rep(NA,length(colnames(DfStAn_mrn)[3:length(colnames(DfStAn_mrn))]))),
                          
                          N2_vs_VC40_pval = as.numeric(rep(NA,length(colnames(DfStAn_mrn)[3:length(colnames(DfStAn_mrn))]))),
                          N2_vs_VC1668_pval = as.numeric(rep(NA,length(colnames(DfStAn_mrn)[3:length(colnames(DfStAn_mrn))]))),
                          N2_vs_UA57_pval = as.numeric(rep(NA,length(colnames(DfStAn_mrn)[3:length(colnames(DfStAn_mrn))]))),
                          N2_vs_BR5270_pval = as.numeric(rep(NA,length(colnames(DfStAn_mrn)[3:length(colnames(DfStAn_mrn))]))),
                          VC40_vs_VC1668_pval = as.numeric(rep(NA,length(colnames(DfStAn_mrn)[3:length(colnames(DfStAn_mrn))]))),
                          VC40_vs_UA57_pval = as.numeric(rep(NA,length(colnames(DfStAn_mrn)[3:length(colnames(DfStAn_mrn))]))),
                          VC40_vs_BR5270_pval = as.numeric(rep(NA,length(colnames(DfStAn_mrn)[3:length(colnames(DfStAn_mrn))]))),
                          VC1668_vs_UA57_pval = as.numeric(rep(NA,length(colnames(DfStAn_mrn)[3:length(colnames(DfStAn_mrn))]))),
                          VC1668_vs_BR5270_pval = as.numeric(rep(NA,length(colnames(DfStAn_mrn)[3:length(colnames(DfStAn_mrn))]))),
                          UA57_vs_BR5270_pval = as.numeric(rep(NA,length(colnames(DfStAn_mrn)[3:length(colnames(DfStAn_mrn))]))),
                          
                          N2_vs_VC40_higlow = as.character(rep(NA,length(colnames(DfStAn_mrn)[3:length(colnames(DfStAn_mrn))]))),
                          N2_vs_VC1668_higlow = as.character(rep(NA,length(colnames(DfStAn_mrn)[3:length(colnames(DfStAn_mrn))]))),
                          N2_vs_UA57_higlow = as.character(rep(NA,length(colnames(DfStAn_mrn)[3:length(colnames(DfStAn_mrn))]))),
                          N2_vs_BR5270_higlow = as.character(rep(NA,length(colnames(DfStAn_mrn)[3:length(colnames(DfStAn_mrn))]))),
                          VC40_vs_VC1668_higlow = as.character(rep(NA,length(colnames(DfStAn_mrn)[3:length(colnames(DfStAn_mrn))]))),
                          VC40_vs_UA57_higlow = as.character(rep(NA,length(colnames(DfStAn_mrn)[3:length(colnames(DfStAn_mrn))]))),
                          VC40_vs_BR5270_higlow = as.character(rep(NA,length(colnames(DfStAn_mrn)[3:length(colnames(DfStAn_mrn))]))),
                          VC1668_vs_UA57_higlow = as.character(rep(NA,length(colnames(DfStAn_mrn)[3:length(colnames(DfStAn_mrn))]))),
                          VC1668_vs_BR5270_higlow = as.character(rep(NA,length(colnames(DfStAn_mrn)[3:length(colnames(DfStAn_mrn))]))),
                          UA57_vs_BR5270_higlow = as.character(rep(NA,length(colnames(DfStAn_mrn)[3:length(colnames(DfStAn_mrn))]))))
  
  
  for (a in colnames(DfStAn_mrn)[3:length(colnames(DfStAn_mrn))]) {
    
    this_index <- which(ANOVA_results$feature == a)
    
    if (length(this_index) != 1) {stop("there is something wrong!")}
    
    ANOVA_this_feat <- aov(formula = as.formula(paste(a, "~ Sample_ID")), data = DfStAn_mrn)
    
    LSD_this_feat <- LSD.test(ANOVA_this_feat,"Sample_ID")
    
    PAIRWISE_PVALUE_this_feat <- pairwise.t.test(pull(DfStAn_mrn, a), pull(DfStAn_mrn, "Sample_ID"), p.adjust.method = "none")
    
    
    ANOVA_results[this_index, "ANOVA_p_value"] <- summary(ANOVA_this_feat)[[1]][["Pr(>F)"]][1]
    
    if (length(colnames(PAIRWISE_PVALUE_this_feat[["p.value"]])) == 4 & length(rownames(PAIRWISE_PVALUE_this_feat[["p.value"]])) == 4) {
      
      if (mean(colnames(PAIRWISE_PVALUE_this_feat[["p.value"]]) == c("N2", "VC40", "VC1668", "BR5270")) == 1 & mean(rownames(PAIRWISE_PVALUE_this_feat[["p.value"]]) == c("VC40", "VC1668", "BR5270", "UA57")) == 1) {
        ANOVA_results[this_index, "N2_vs_VC40_pval"] <- PAIRWISE_PVALUE_this_feat[["p.value"]]["VC40", "N2"]
        ANOVA_results[this_index, "N2_vs_VC1668_pval"] <- PAIRWISE_PVALUE_this_feat[["p.value"]]["VC1668", "N2"]
        ANOVA_results[this_index, "N2_vs_UA57_pval"] <- PAIRWISE_PVALUE_this_feat[["p.value"]]["UA57", "N2"]
        ANOVA_results[this_index, "N2_vs_BR5270_pval"] <- PAIRWISE_PVALUE_this_feat[["p.value"]]["BR5270", "N2"]
        ANOVA_results[this_index, "VC40_vs_VC1668_pval"] <- PAIRWISE_PVALUE_this_feat[["p.value"]]["VC1668", "VC40"]
        ANOVA_results[this_index, "VC40_vs_UA57_pval"] <- PAIRWISE_PVALUE_this_feat[["p.value"]]["UA57", "VC40"]
        ANOVA_results[this_index, "VC40_vs_BR5270_pval"] <- PAIRWISE_PVALUE_this_feat[["p.value"]]["BR5270", "VC40"]
        ANOVA_results[this_index, "VC1668_vs_UA57_pval"] <- PAIRWISE_PVALUE_this_feat[["p.value"]]["UA57", "VC1668"]
        ANOVA_results[this_index, "VC1668_vs_BR5270_pval"] <- PAIRWISE_PVALUE_this_feat[["p.value"]]["BR5270", "VC1668"]
        ANOVA_results[this_index, "UA57_vs_BR5270_pval"] <- PAIRWISE_PVALUE_this_feat[["p.value"]]["UA57", "BR5270"]
        
        if(ANOVA_results[this_index, "N2_vs_VC40_pval"] > 0.05) {
          ANOVA_results[this_index, "N2_vs_VC40_higlow"] <- NA
        } else {
          if (LSD_this_feat[["means"]]["N2",1] > LSD_this_feat[["means"]]["VC40", 1]) {
            ANOVA_results[this_index, "N2_vs_VC40_higlow"] <- "N2 > VC40"
          } else {
            ANOVA_results[this_index, "N2_vs_VC40_higlow"] <- "VC40 > N2"
          }
        }
        
        if(ANOVA_results[this_index, "N2_vs_VC1668_pval"] > 0.05) {
          ANOVA_results[this_index, "N2_vs_VC1668_higlow"] <- NA
        } else {
          if (LSD_this_feat[["means"]]["N2",1] > LSD_this_feat[["means"]]["VC1668", 1]) {
            ANOVA_results[this_index, "N2_vs_VC1668_higlow"] <- "N2 > VC1668"
          } else {
            ANOVA_results[this_index, "N2_vs_VC1668_higlow"] <- "VC1668 > N2"
          }
        }
        
        if(ANOVA_results[this_index, "N2_vs_UA57_pval"] > 0.05) {
          ANOVA_results[this_index, "N2_vs_UA57_higlow"] <- NA
        } else {
          if (LSD_this_feat[["means"]]["N2",1] > LSD_this_feat[["means"]]["UA57", 1]) {
            ANOVA_results[this_index, "N2_vs_UA57_higlow"] <- "N2 > UA57"
          } else {
            ANOVA_results[this_index, "N2_vs_UA57_higlow"] <- "UA57 > N2"
          }
        }
        
        if(ANOVA_results[this_index, "N2_vs_BR5270_pval"] > 0.05) {
          ANOVA_results[this_index, "N2_vs_BR5270_higlow"] <- NA
        } else {
          if (LSD_this_feat[["means"]]["N2",1] > LSD_this_feat[["means"]]["BR5270", 1]) {
            ANOVA_results[this_index, "N2_vs_BR5270_higlow"] <- "N2 > BR5270"
          } else {
            ANOVA_results[this_index, "N2_vs_BR5270_higlow"] <- "BR5270 > N2"
          }
        }
        
        if(ANOVA_results[this_index, "VC40_vs_VC1668_pval"] > 0.05) {
          ANOVA_results[this_index, "VC40_vs_VC1668_higlow"] <- NA
        } else {
          if (LSD_this_feat[["means"]]["VC40",1] > LSD_this_feat[["means"]]["VC1668", 1]) {
            ANOVA_results[this_index, "VC40_vs_VC1668_higlow"] <- "VC40 > VC1668"
          } else {
            ANOVA_results[this_index, "VC40_vs_VC1668_higlow"] <- "VC1668 > VC40"
          }
        }
        
        if(ANOVA_results[this_index, "VC40_vs_UA57_pval"] > 0.05) {
          ANOVA_results[this_index, "VC40_vs_UA57_higlow"] <- NA
        } else {
          if (LSD_this_feat[["means"]]["VC40",1] > LSD_this_feat[["means"]]["UA57", 1]) {
            ANOVA_results[this_index, "VC40_vs_UA57_higlow"] <- "VC40 > UA57"
          } else {
            ANOVA_results[this_index, "VC40_vs_UA57_higlow"] <- "UA57 > VC40"
          }
        }
        
        if(ANOVA_results[this_index, "VC40_vs_BR5270_pval"] > 0.05) {
          ANOVA_results[this_index, "VC40_vs_BR5270"] <- NA
        } else {
          if (LSD_this_feat[["means"]]["VC40",1] > LSD_this_feat[["means"]]["BR5270", 1]) {
            ANOVA_results[this_index, "VC40_vs_BR5270"] <- "VC40 > BR5270"
          } else {
            ANOVA_results[this_index, "VC40_vs_BR5270"] <- "BR5270 > VC40"
          }
        }
        
        if(ANOVA_results[this_index, "VC1668_vs_UA57_pval"] > 0.05) {
          ANOVA_results[this_index, "VC1668_vs_UA57_higlow"] <- NA
        } else {
          if (LSD_this_feat[["means"]]["VC1668",1] > LSD_this_feat[["means"]]["UA57", 1]) {
            ANOVA_results[this_index, "VC1668_vs_UA57_higlow"] <- "VC1668 > UA57"
          } else {
            ANOVA_results[this_index, "VC1668_vs_UA57_higlow"] <- "UA57 > VC1668"
          }
        }
        
        if(ANOVA_results[this_index, "VC1668_vs_BR5270_pval"] > 0.05) {
          ANOVA_results[this_index, "VC1668_vs_BR5270_higlow"] <- NA
        } else {
          if (LSD_this_feat[["means"]]["VC1668",1] > LSD_this_feat[["means"]]["BR5270", 1]) {
            ANOVA_results[this_index, "VC1668_vs_BR5270_higlow"] <- "VC1668 > BR5270"
          } else {
            ANOVA_results[this_index, "VC1668_vs_BR5270_higlow"] <- "BR5270 > VC1668"
          }
        }
        
        if(ANOVA_results[this_index, "UA57_vs_BR5270_pval"] > 0.05) {
          ANOVA_results[this_index, "UA57_vs_BR5270_higlow"] <- NA
        } else {
          if (LSD_this_feat[["means"]]["UA57",1] > LSD_this_feat[["means"]]["BR5270", 1]) {
            ANOVA_results[this_index, "UA57_vs_BR5270_higlow"] <- "UA57 > BR5270"
          } else {
            ANOVA_results[this_index, "UA57_vs_BR5270_higlow"] <- "BR5270 > UA57"
          }
        }
      } else {
        stop("something wrong, try to check")
      }
      
      
    } else if (length(colnames(PAIRWISE_PVALUE_this_feat[["p.value"]])) == 3 & length(rownames(PAIRWISE_PVALUE_this_feat[["p.value"]])) == 3) {
      if (mean(colnames(PAIRWISE_PVALUE_this_feat[["p.value"]]) == c("VC40", "VC1668", "BR5270")) == 1 & mean(rownames(PAIRWISE_PVALUE_this_feat[["p.value"]]) == c("VC1668", "BR5270", "UA57")) == 1) {
        
        ANOVA_results[this_index, "VC40_vs_VC1668_pval"] <- PAIRWISE_PVALUE_this_feat[["p.value"]]["VC1668", "VC40"]
        ANOVA_results[this_index, "VC40_vs_UA57_pval"] <- PAIRWISE_PVALUE_this_feat[["p.value"]]["UA57", "VC40"]
        ANOVA_results[this_index, "VC40_vs_BR5270_pval"] <- PAIRWISE_PVALUE_this_feat[["p.value"]]["BR5270", "VC40"]
        ANOVA_results[this_index, "VC1668_vs_UA57_pval"] <- PAIRWISE_PVALUE_this_feat[["p.value"]]["UA57", "VC1668"]
        ANOVA_results[this_index, "VC1668_vs_BR5270_pval"] <- PAIRWISE_PVALUE_this_feat[["p.value"]]["BR5270", "VC1668"]
        ANOVA_results[this_index, "UA57_vs_BR5270_pval"] <- PAIRWISE_PVALUE_this_feat[["p.value"]]["UA57", "BR5270"]
        
        
        if(ANOVA_results[this_index, "VC40_vs_VC1668_pval"] > 0.05) {
          ANOVA_results[this_index, "VC40_vs_VC1668_higlow"] <- NA
        } else {
          if (LSD_this_feat[["means"]]["VC40",1] > LSD_this_feat[["means"]]["VC1668", 1]) {
            ANOVA_results[this_index, "VC40_vs_VC1668_higlow"] <- "VC40 > VC1668"
          } else {
            ANOVA_results[this_index, "VC40_vs_VC1668_higlow"] <- "VC1668 > VC40"
          }
        }
        
        if(ANOVA_results[this_index, "VC40_vs_UA57_pval"] > 0.05) {
          ANOVA_results[this_index, "VC40_vs_UA57_higlow"] <- NA
        } else {
          if (LSD_this_feat[["means"]]["VC40",1] > LSD_this_feat[["means"]]["UA57", 1]) {
            ANOVA_results[this_index, "VC40_vs_UA57_higlow"] <- "VC40 > UA57"
          } else {
            ANOVA_results[this_index, "VC40_vs_UA57_higlow"] <- "UA57 > VC40"
          }
        }
        
        if(ANOVA_results[this_index, "VC40_vs_BR5270_pval"] > 0.05) {
          ANOVA_results[this_index, "VC40_vs_BR5270"] <- NA
        } else {
          if (LSD_this_feat[["means"]]["VC40",1] > LSD_this_feat[["means"]]["BR5270", 1]) {
            ANOVA_results[this_index, "VC40_vs_BR5270"] <- "VC40 > BR5270"
          } else {
            ANOVA_results[this_index, "VC40_vs_BR5270"] <- "BR5270 > VC40"
          }
        }
        
        if(ANOVA_results[this_index, "VC1668_vs_UA57_pval"] > 0.05) {
          ANOVA_results[this_index, "VC1668_vs_UA57_higlow"] <- NA
        } else {
          if (LSD_this_feat[["means"]]["VC1668",1] > LSD_this_feat[["means"]]["UA57", 1]) {
            ANOVA_results[this_index, "VC1668_vs_UA57_higlow"] <- "VC1668 > UA57"
          } else {
            ANOVA_results[this_index, "VC1668_vs_UA57_higlow"] <- "UA57 > VC1668"
          }
        }
        
        if(ANOVA_results[this_index, "VC1668_vs_BR5270_pval"] > 0.05) {
          ANOVA_results[this_index, "VC1668_vs_BR5270_higlow"] <- NA
        } else {
          if (LSD_this_feat[["means"]]["VC1668",1] > LSD_this_feat[["means"]]["BR5270", 1]) {
            ANOVA_results[this_index, "VC1668_vs_BR5270_higlow"] <- "VC1668 > BR5270"
          } else {
            ANOVA_results[this_index, "VC1668_vs_BR5270_higlow"] <- "BR5270 > VC1668"
          }
        }
        
        if(ANOVA_results[this_index, "UA57_vs_BR5270_pval"] > 0.05) {
          ANOVA_results[this_index, "UA57_vs_BR5270_higlow"] <- NA
        } else {
          if (LSD_this_feat[["means"]]["UA57",1] > LSD_this_feat[["means"]]["BR5270", 1]) {
            ANOVA_results[this_index, "UA57_vs_BR5270_higlow"] <- "UA57 > BR5270"
          } else {
            ANOVA_results[this_index, "UA57_vs_BR5270_higlow"] <- "BR5270 > UA57"
          }
        }
      } else {stop("something wrong, try to check")}
    } else {stop("something wrong, try to check")}
  }
  
  
  ANOVA_results <- add_column(ANOVA_results,
                              .after = "ANOVA_p_value",
                              ANOVA_p_value_FDR = p.adjust(ANOVA_results$ANOVA_p_value, method = "fdr"))
  
  return(ANOVA_results)
}




##


HLP4_shem1_dry_transf <- transf_feat_intensities(HLP4_shem1_dry)
HLP4_shem1_dry_ANOVA <- getting_ANOVA_results(HLP4_shem1_dry_transf)


HLP4_shem2_dry_transf <- transf_feat_intensities(HLP4_shem2_dry) 
HLP4_shem2_dry_ANOVA <- getting_ANOVA_results(HLP4_shem2_dry_transf)



write_tsv(HLP4_shem1_dry_transf, "MSDIAL_WJ_HLP4_shem1_dry_transf.txt")
write_tsv(HLP4_shem1_dry_ANOVA, "MSDIAL_WJ_HLP4_shem1_dry_ANOVA.txt")

write_tsv(HLP4_shem2_dry_transf, "MSDIAL_WJ_HLP4_shem2_dry_transf.txt")
write_tsv(HLP4_shem2_dry_ANOVA, "MSDIAL_WJ_HLP4_shem2_dry_ANOVA.txt")




## only ANOVA significant:


HLP4_shem1_dry_ANOVA_sign <- filter(HLP4_shem1_dry_ANOVA, ANOVA_p_value_FDR < 0.05)

HLP4_shem2_dry_ANOVA_sign <- filter(HLP4_shem2_dry_ANOVA, ANOVA_p_value_FDR < 0.05)



write_tsv(HLP4_shem1_dry_ANOVA_sign, "MSDIAL_HLP4_shem1_dry_ANOVA_sign.txt")

write_tsv(HLP4_shem2_dry_ANOVA_sign, "MSDIAL_HLP4_shem2_dry_ANOVA_sign.txt")







## Eulero-Venn with significant features:

HLP_EV_df_QCfil_sign <- tibble(feature = HLP_unique_feat_vector,
                               dry_scheme1 = HLP_unique_feat_vector %in% HLP4_shem1_dry_ANOVA_sign$feature,
                               dry_scheme2bottom = HLP_unique_feat_vector %in% HLP4_shem2_dry_ANOVA_sign$feature)



# better:

HLP_EV_df_QCfil_sign_better <- HLP_EV_df_QCfil_sign
colnames(HLP_EV_df_QCfil_sign_better) <- c("feature", "scheme1", "scheme2")


jpeg(file="HILIC_POS_EuleroVenn_ANOVA.jpeg")
plot(euler(HLP_EV_df_QCfil_sign_better[, 2:3], shape = "ellipse"), fill = c("cadetblue1", "indianred1"), quantities = list(cex = 2.5), labels = list(cex = 2))
dev.off()





## Actually, I added Eulero-Venn with features also before the filtering with QCs

all_feat <- HLP0$Alignment_ID



# dried scheme 1

feat_HLP0_dry_N2_ACN <- HLP0$Alignment_ID[!is.na(HLP0$X136_Pool_N2_ACN_pool_20220714203441) | !is.na(HLP0$X136_Pool_N2_ACN_pool_20220715004714) | !is.na(HLP0$X136_Pool_N2_ACN_pool_20220715072327)]
feat_HLP0_dry_N2_ACN_pres <- HLP0$Alignment_ID %in% feat_HLP0_dry_N2_ACN

feat_HLP0_dry_VC40_ACN <- HLP0$Alignment_ID[!is.na(HLP0$X137_Pool_VC40_ACN_pool_20220714205244) | !is.na(HLP0$X137_Pool_VC40_ACN_pool_20220715010516) | !is.na(HLP0$X137_Pool_VC40_ACN_pool_20220715074131)]
feat_HLP0_dry_VC40_ACN_pres <- HLP0$Alignment_ID %in% feat_HLP0_dry_VC40_ACN

feat_HLP0_dry_VC1668_ACN <- HLP0$Alignment_ID[!is.na(HLP0$X138_Pool_VC1668_ACN_pool_20220714211048) | !is.na(HLP0$X138_Pool_VC1668_ACN_pool_20220715012316) | !is.na(HLP0$X138_Pool_VC1668_ACN_pool_20220715075934)]
feat_HLP0_dry_VC1668_ACN_pres <- HLP0$Alignment_ID %in% feat_HLP0_dry_VC1668_ACN

feat_HLP0_dry_BR5270_ACN <- HLP0$Alignment_ID[!is.na(HLP0$X139_Pool_BR5270_ACN_pool_20220714212851) | !is.na(HLP0$X139_Pool_BR5270_ACN_pool_20220715014117) | !is.na(HLP0$X139_Pool_BR5270_ACN_pool_20220715081739)]
feat_HLP0_dry_BR5270_ACN_pres <- HLP0$Alignment_ID %in% feat_HLP0_dry_BR5270_ACN

feat_HLP0_dry_UA57_ACN <- HLP0$Alignment_ID[!is.na(HLP0$X140_Pool_UA57_ACN_pool_20220714214654) | !is.na(HLP0$X140_Pool_UA57_ACN_pool_20220715015918) | !is.na(HLP0$X140_Pool_UA57_ACN_pool_20220715083543)]
feat_HLP0_dry_UA57_ACN_pres <- HLP0$Alignment_ID %in% feat_HLP0_dry_UA57_ACN

feat_HLP0_dry_VC40_ACN_bis <- HLP0$Alignment_ID[!is.na(HLP0$X143_Pool_VC40_ACN_pool_20220714220458) | !is.na(HLP0$X143_Pool_VC40_ACN_pool_20220715021719) | !is.na(HLP0$X143_Pool_VC40_ACN_pool_20220715085345)]
feat_HLP0_dry_VC40_ACN_bis_pres <- HLP0$Alignment_ID %in% feat_HLP0_dry_VC40_ACN_bis


#  dry_2B

feat_HLP0_dry_N2_2B <- HLP0$Alignment_ID[!is.na(HLP0$X424_.dry._Pool_N2_2B_pool_20220716105522) | !is.na(HLP0$X424_.dry._Pool_N2_2B_pool_20220716150725) | !is.na(HLP0$X424_.dry._Pool_N2_2B_pool_20220716214339)]
feat_HLP0_dry_N2_2B_pres <- HLP0$Alignment_ID %in% feat_HLP0_dry_N2_2B

feat_HLP0_dry_VC40_2B <- HLP0$Alignment_ID[!is.na(HLP0$X425_.dry._Pool_VC40_2B_pool_20220716111323) | !is.na(HLP0$X425_.dry._Pool_VC40_2B_pool_20220716152528) | !is.na(HLP0$X425_.dry._Pool_VC40_2B_pool_20220716220141)]
feat_HLP0_dry_VC40_2B_pres <- HLP0$Alignment_ID %in% feat_HLP0_dry_VC40_2B

feat_HLP0_dry_VC1668_2B <- HLP0$Alignment_ID[!is.na(HLP0$X426_.dry._Pool_VC1668_2B_pool_20220716113126) | !is.na(HLP0$X426_.dry._Pool_VC1668_2B_pool_20220716154332) | !is.na(HLP0$X426_.dry._Pool_VC1668_2B_pool_20220716221942)]
feat_HLP0_dry_VC1668_2B_pres <- HLP0$Alignment_ID %in% feat_HLP0_dry_VC1668_2B

feat_HLP0_dry_BR5270_2B <- HLP0$Alignment_ID[!is.na(HLP0$X427_.dry._Pool_BR5270_2B_pool_20220716114922) | !is.na(HLP0$X427_.dry._Pool_BR5270_2B_pool_20220716160136) | !is.na(HLP0$X427_.dry._Pool_BR5270_2B_pool_20220716223744)]
feat_HLP0_dry_BR5270_2B_pres <- HLP0$Alignment_ID %in% feat_HLP0_dry_BR5270_2B

feat_HLP0_dry_UA57_2B <- HLP0$Alignment_ID[!is.na(HLP0$X428_.dry._Pool_UA57_2B_pool_20220716120723) | !is.na(HLP0$X428_.dry._Pool_UA57_2B_pool_20220716161939) | !is.na(HLP0$X428_.dry._Pool_UA57_2B_pool_20220716225547)]
feat_HLP0_dry_UA57_2B_pres <- HLP0$Alignment_ID %in% feat_HLP0_dry_UA57_2B

feat_HLP0_dry_VC40_2B_bis <- HLP0$Alignment_ID[!is.na(HLP0$X431_.dry._Pool_VC40_2B_pool_20220716122523) | !is.na(HLP0$X431_.dry._Pool_VC40_2B_pool_20220716163744) | !is.na(HLP0$X431_.dry._Pool_VC40_2B_pool_20220716231347)]
feat_HLP0_dry_VC40_2B_bis_pres <- HLP0$Alignment_ID %in% feat_HLP0_dry_VC40_2B_bis


HLP_EV_df_tot <- tibble(feature = all_feat,
                        dry_scheme1 = feat_HLP0_dry_N2_ACN_pres | feat_HLP0_dry_VC40_ACN_pres | feat_HLP0_dry_VC1668_ACN_pres | feat_HLP0_dry_BR5270_ACN_pres | feat_HLP0_dry_UA57_ACN_pres | feat_HLP0_dry_VC40_ACN_bis_pres,
                        dry_scheme2bottom = feat_HLP0_dry_N2_2B_pres | feat_HLP0_dry_VC40_2B_pres | feat_HLP0_dry_VC1668_2B_pres | feat_HLP0_dry_BR5270_2B_pres | feat_HLP0_dry_UA57_2B_pres | feat_HLP0_dry_VC40_2B_bis_pres)


HLP_EV_df_tot_better <- HLP_EV_df_tot
colnames(HLP_EV_df_tot_better) <- c("feature", "scheme1", "scheme2")


jpeg(file="HILIC_POS_EuleroVenn_tot.jpeg")
plot(euler(HLP_EV_df_tot_better[, 2:3], shape = "ellipse"), fill = c("cadetblue1", "indianred1"), quantities = list(cex = 2.5), labels = list(cex = 2))
dev.off()




###
# combining significant features with top candidate in annotation




HLP4_shem1_dry_ANOVA_sign_annotation <- HLP4_shem1_dry_ANOVA_sign

for (a in colnames(HLP_feat_info)) {
  HLP4_shem1_dry_ANOVA_sign_annotation[,a] <- rep(NA, length(HLP4_shem1_dry_ANOVA_sign_annotation$feature))
  
  HLP4_shem1_dry_ANOVA_sign_annotation[[a]] <- as(HLP4_shem1_dry_ANOVA_sign_annotation[[a]], Class = class(HLP_feat_info[[a]]))
  
}


for (i in 1:length(HLP4_shem1_dry_ANOVA_sign_annotation$feature)) {
  for (u in 1:length(HLP_feat_info$Alignment_ID)) {
    if (HLP4_shem1_dry_ANOVA_sign_annotation$feature[i] == HLP_feat_info$Alignment_ID[u]) {
      HLP4_shem1_dry_ANOVA_sign_annotation[i,colnames(HLP_feat_info)] <- HLP_feat_info[u,]
    }
  }
}

write_tsv(HLP4_shem1_dry_ANOVA_sign_annotation, "MSDIAL_HLP4_shem1_dry_ANOVA_sign_annotation.txt")






HLP4_shem2_dry_ANOVA_sign_annotation <- HLP4_shem2_dry_ANOVA_sign

for (a in colnames(HLP_feat_info)) {
  HLP4_shem2_dry_ANOVA_sign_annotation[,a] <- rep(NA, length(HLP4_shem2_dry_ANOVA_sign_annotation$feature))
  
  HLP4_shem2_dry_ANOVA_sign_annotation[[a]] <- as(HLP4_shem2_dry_ANOVA_sign_annotation[[a]], Class = class(HLP_feat_info[[a]]))
  
}


for (i in 1:length(HLP4_shem2_dry_ANOVA_sign_annotation$feature)) {
  for (u in 1:length(HLP_feat_info$Alignment_ID)) {
    if (HLP4_shem2_dry_ANOVA_sign_annotation$feature[i] == HLP_feat_info$Alignment_ID[u]) {
      HLP4_shem2_dry_ANOVA_sign_annotation[i,colnames(HLP_feat_info)] <- HLP_feat_info[u,]
    }
  }
}

write_tsv(HLP4_shem2_dry_ANOVA_sign_annotation, "MSDIAL_HLP4_shem2_dry_ANOVA_sign_annotation.txt")






#
##
###### annotation levels, considering the full list of features:


HLP_feat_info_annolevels <- add_column(HLP_feat_info, AnnoLevel = NA, .after = 3) %>%
  mutate(AnnoLevel = ifelse(MS_MS_assigned == TRUE & Matched_peaks_count != "null" & Matched_peaks_count >= 3 & Weighted_dot_product != "null" & Weighted_dot_product >= 0.7 & Matched_peaks_percentage != "null" & Matched_peaks_percentage >= 0.5 & INCHIKEY != "null", "2a",
                            ifelse(MS_MS_assigned == TRUE & Matched_peaks_count != "null" & Matched_peaks_count >= 3 & Weighted_dot_product != "null" & Weighted_dot_product >= 0.7 & Matched_peaks_percentage != "null" & Matched_peaks_percentage >= 0.5 & INCHIKEY == "null", "2b",
                                   ifelse(MS_MS_assigned == TRUE & Matched_peaks_count != "null" & Matched_peaks_count >= 3 & Weighted_dot_product != "null" & Weighted_dot_product >= 0.5 & Weighted_dot_product < 0.7 & Matched_peaks_percentage != "null" & Matched_peaks_percentage >= 0.5, "3a",
                                          ifelse(MS_MS_assigned == TRUE & Matched_peaks_count != "null" & Matched_peaks_count < 3 & Weighted_dot_product != "null" & Weighted_dot_product >= 0.5 & Matched_peaks_percentage != "null" & Matched_peaks_percentage >= 0.5 & INCHIKEY != "null", "3b",
                                                 ifelse(MS_MS_assigned == TRUE & Matched_peaks_count != "null" & Matched_peaks_count < 3 & Weighted_dot_product != "null" & Weighted_dot_product >= 0.5 & Matched_peaks_percentage != "null" & Matched_peaks_percentage >= 0.5 & INCHIKEY == "null", "3c",
                                                        ifelse(MS_MS_assigned == TRUE & Weighted_dot_product != "null" & Weighted_dot_product < 0.5 & Matched_peaks_percentage != "null" & Matched_peaks_percentage < 0.5, "4a",
                                                               ifelse(MS_MS_assigned == FALSE & Formula != "null", "4a", "5"))))))))


write_tsv(HLP_feat_info_annolevels, "HILICPOS_MSDial_WJ_feat_info_annolevels.txt")


# QC filtered:

HLP_feat_info_annolevels_QCfiltered <- filter(HLP_feat_info_annolevels, Alignment_ID %in% HLP_unique_feat_vector)

write_tsv(HLP_feat_info_annolevels_QCfiltered, "HILICPOS_MSDial_WJ_feat_info_annolevels_QCfiltered.txt")




# QC filtered, per sample type:

HLP_feat_info_annolevels_QCfiltered_scheme1_dry <- filter(HLP_feat_info_annolevels, Alignment_ID %in% HLP_unique_feat_vector_dry_schem1)

write_tsv(HLP_feat_info_annolevels_QCfiltered_scheme1_dry, "HILICPOS_MSDial_WJ_feat_info_annolevels_QCfiltered_scheme1_dry.txt")



HLP_feat_info_annolevels_QCfiltered_scheme2_dry <- filter(HLP_feat_info_annolevels, Alignment_ID %in% HLP_unique_feat_vector_dry_schem2)

write_tsv(HLP_feat_info_annolevels_QCfiltered_scheme2_dry, "HILICPOS_MSDial_WJ_feat_info_annolevels_QCfiltered_scheme2_dry.txt")






## annotation levels, only for statistically significant:


HLP_feat_info_annolevel_empty_dried_scheme1 <- HLP_feat_info_annolevels[0,]
HLP_feat_info_annolevel_empty_dried_scheme1[1:length(HLP4_shem1_dry_ANOVA_sign$feature),] <- NA


HLP4_shem1_dry_ANOVA_sign_annotat_withLevels <- cbind(HLP4_shem1_dry_ANOVA_sign, HLP_feat_info_annolevel_empty_dried_scheme1)

for (i in 1:length(HLP4_shem1_dry_ANOVA_sign$feature)) {
  indices <- which(HLP_feat_info_annolevels$Alignment_ID == HLP4_shem1_dry_ANOVA_sign$feature[i])
  
  HLP4_shem1_dry_ANOVA_sign_annotat_withLevels[i,colnames(HLP_feat_info_annolevel_empty_dried_scheme1)] <- HLP_feat_info_annolevels[indices,]
  
}

write_tsv(HLP4_shem1_dry_ANOVA_sign_annotat_withLevels, "HILICPOS_MSDial_WJ_scheme1_dry_ANOVA_sign_annotat_withLevels.txt")







HLP_feat_info_annolevel_empty_dried_scheme2bottom <- HLP_feat_info_annolevels[0,]
HLP_feat_info_annolevel_empty_dried_scheme2bottom[1:length(HLP4_shem2_dry_ANOVA_sign$feature),] <- NA


HLP4_shem2_dry_ANOVA_sign_annotat_withLevels <- cbind(HLP4_shem2_dry_ANOVA_sign, HLP_feat_info_annolevel_empty_dried_scheme2bottom)

for (i in 1:length(HLP4_shem2_dry_ANOVA_sign$feature)) {
  indices <- which(HLP_feat_info_annolevels$Alignment_ID == HLP4_shem2_dry_ANOVA_sign$feature[i])
  
  HLP4_shem2_dry_ANOVA_sign_annotat_withLevels[i,colnames(HLP_feat_info_annolevel_empty_dried_scheme2bottom)] <- HLP_feat_info_annolevels[indices,]
  
}

write_tsv(HLP4_shem2_dry_ANOVA_sign_annotat_withLevels, "HILICPOS_MSDial_WJ_scheme2_dry_ANOVA_sign_annotat_withLevels.txt")



