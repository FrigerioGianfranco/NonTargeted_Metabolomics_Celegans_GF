

library(tidyverse)
library(eulerr)
library(agricolae) #for the Fisher LSD
set.seed(1991)

#### function to export a vector with only features with no NAs:

get_noNAs_features <- function(dataset) {
  
  return_vector <- character()
  
  for (i in 1:length(pull(dataset, 1))) {
    
    if (!all(is.na(as.numeric(dataset[i, 4:length(colnames(dataset))])))) {
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

RPNnrf <- read_csv("featureGroupsXCMS3_updt.csv") %>%
  select(-adduct, -neutralMass)

## adding an "X" at the beginning of each sample that starts with a number and replacing parenthesis with dots

colnames(RPNnrf)[grepl("^[[:digit:]]+", colnames(RPNnrf))] <- paste0("X", colnames(RPNnrf)[grepl("^[[:digit:]]+", colnames(RPNnrf))])
colnames(RPNnrf) <- str_replace_all(colnames(RPNnrf), "[//(//)]", ".")
colnames(RPNnrf) <- str_replace_all(colnames(RPNnrf), "[//)//)]", ".")


## replacing zero with NA

RPNnrf[RPNnrf == 0] <- NA



all_featuresRPNnrf <- get_noNAs_features(RPNnrf)

analyses_dry_samples <- colnames(RPNnrf)[which(grepl(paste(dry_samples_codes, collapse = "|"), colnames(RPNnrf)))]


RPNnf <- select(RPNnrf, all_of(c("group", "ret", "mz", analyses_dry_samples)))

all_featuresRPNnf <- get_noNAs_features(RPNnf)

RPN <- filter(RPNnf, group %in% all_featuresRPNnf)

write_tsv(RPN, "RPN.txt")


## removing the first analyses of the run, for each category of samples, and creating the feature numbers

RPN0nf <- select(RPN,
                 
                 -X136_Pool_N2_ACN_pool,
                 -X137_Pool_VC40_ACN_pool,
                 -X138_Pool_VC1668_ACN_pool,
                 -X139_Pool_BR5270_ACN_pool,
                 -X140_Pool_UA57_ACN_pool,
                 -X143_Pool_VC40_ACN_pool,
                 
                 -X280_.dry._Pool_N2_2U_pool,
                 -X281_.dry._Pool_VC40_2U_pool,
                 -X282_.dry._Pool_VC1668_2U_pool,
                 -X283_.dry._Pool_BR5270_2U_pool,
                 -X284_.dry._Pool_UA57_2U_pool,
                 -X287_.dry._Pool_VC40_2U_pool)



all_featuresRPN0nf <- get_noNAs_features(RPN0nf)

RPN0 <- filter(RPN0nf, group %in% all_featuresRPN0nf)

write_tsv(RPN0, "RPN0.txt")

RPN_feat_info <- select(RPN0, group, ret, mz)
write_tsv(RPN_feat_info, "RPN_feat_info.txt")






## Doing all the separated QC cleaning


## the dry samples


# dry_N2_ACN

RPN0_dry_N2_ACN <- RPN0 %>%
  mutate(X083_.dry._N2_blank_ACN_1 = ifelse(is.na(X083_.dry._N2_blank_ACN_1), 0, X083_.dry._N2_blank_ACN_1),
         X084_.dry._N2_blank_ACN_2 = ifelse(is.na(X084_.dry._N2_blank_ACN_2), 0, X084_.dry._N2_blank_ACN_2)) %>%
  rowwise() %>%
  mutate(MEDIA = mean(c(X136_Pool_N2_ACN_pool_20220712144651,
                        X136_Pool_N2_ACN_pool_20220712182851,
                        X136_Pool_N2_ACN_pool_20220713001646), na.rm=TRUE),
         SD = sd(c(X136_Pool_N2_ACN_pool_20220712144651,
                   X136_Pool_N2_ACN_pool_20220712182851,
                   X136_Pool_N2_ACN_pool_20220713001646), na.rm=TRUE),
         RSD = SD/MEDIA*100,
         COUNT = mean(!is.na(c(X136_Pool_N2_ACN_pool_20220712144651,
                               X136_Pool_N2_ACN_pool_20220712182851,
                               X136_Pool_N2_ACN_pool_20220713001646)))*100,
         BLANK_CONTR = mean(mean(X083_.dry._N2_blank_ACN_1, X084_.dry._N2_blank_ACN_2)/MEDIA*100)
  )

RPN0_dry_N2_ACNf <- RPN0_dry_N2_ACN %>%
  filter(COUNT>50) %>%
  filter(RSD<50) %>%
  filter(BLANK_CONTR<50)




# dry_VC40_ACN

RPN0_dry_VC40_ACN <- RPN0 %>%
  mutate(X089_.dry._VC40_blank_ACN_1 = ifelse(is.na(X089_.dry._VC40_blank_ACN_1), 0, X089_.dry._VC40_blank_ACN_1),
         X090_.dry._VC40_blank_ACN_2 = ifelse(is.na(X090_.dry._VC40_blank_ACN_2), 0, X090_.dry._VC40_blank_ACN_2)) %>%
  rowwise() %>%
  mutate(MEDIA = mean(c(X137_Pool_VC40_ACN_pool_20220712150243,
                        X137_Pool_VC40_ACN_pool_20220712184443,
                        X137_Pool_VC40_ACN_pool_20220713003236), na.rm=TRUE),
  SD = sd(c(X137_Pool_VC40_ACN_pool_20220712150243,
            X137_Pool_VC40_ACN_pool_20220712184443,
            X137_Pool_VC40_ACN_pool_20220713003236), na.rm=TRUE),
  RSD = SD/MEDIA*100,
  COUNT = mean(!is.na(c(X137_Pool_VC40_ACN_pool_20220712150243,
                        X137_Pool_VC40_ACN_pool_20220712184443,
                        X137_Pool_VC40_ACN_pool_20220713003236)))*100,
  BLANK_CONTR = mean(mean(X089_.dry._VC40_blank_ACN_1, X090_.dry._VC40_blank_ACN_2)/MEDIA*100)
  )

RPN0_dry_VC40_ACNf <- RPN0_dry_VC40_ACN %>%
  filter(COUNT>50) %>%
  filter(RSD<50) %>%
  filter(BLANK_CONTR<50)



# dry_VC1668_ACN

RPN0_dry_VC1668_ACN <- RPN0 %>%
  mutate(X095_.dry._VC1668_blank_ACN_1 = ifelse(is.na(X095_.dry._VC1668_blank_ACN_1), 0, X095_.dry._VC1668_blank_ACN_1),
         X096_.dry._VC1668_blank_ACN_2 = ifelse(is.na(X096_.dry._VC1668_blank_ACN_2), 0, X096_.dry._VC1668_blank_ACN_2)) %>%
  rowwise() %>%
  mutate(MEDIA = mean(c(X138_Pool_VC1668_ACN_pool_20220712151834,
                        X138_Pool_VC1668_ACN_pool_20220712190033,
                        X138_Pool_VC1668_ACN_pool_20220713004824), na.rm=TRUE),
  SD = sd(c(X138_Pool_VC1668_ACN_pool_20220712151834,
            X138_Pool_VC1668_ACN_pool_20220712190033,
            X138_Pool_VC1668_ACN_pool_20220713004824), na.rm=TRUE),
  RSD = SD/MEDIA*100,
  COUNT = mean(!is.na(c(X138_Pool_VC1668_ACN_pool_20220712151834,
                        X138_Pool_VC1668_ACN_pool_20220712190033,
                        X138_Pool_VC1668_ACN_pool_20220713004824)))*100,
  BLANK_CONTR = mean(mean(X095_.dry._VC1668_blank_ACN_1, X096_.dry._VC1668_blank_ACN_2)/MEDIA*100)
  )

RPN0_dry_VC1668_ACNf <- RPN0_dry_VC1668_ACN %>%
  filter(COUNT>50) %>%
  filter(RSD<50) %>%
  filter(BLANK_CONTR<50)



# dry_BR5270_ACN

RPN0_dry_BR5270_ACN <- RPN0 %>%
  mutate(X101_.dry._BR5270_blank_ACN_1 = ifelse(is.na(X101_.dry._BR5270_blank_ACN_1), 0, X101_.dry._BR5270_blank_ACN_1),
         X102_.dry._BR5270_blank_ACN_2 = ifelse(is.na(X102_.dry._BR5270_blank_ACN_2), 0, X102_.dry._BR5270_blank_ACN_2)) %>%
  rowwise() %>%
  mutate(MEDIA = mean(c(X139_Pool_BR5270_ACN_pool_20220712153427,
                        X139_Pool_BR5270_ACN_pool_20220712191626,
                        X139_Pool_BR5270_ACN_pool_20220713010412), na.rm=TRUE),
  SD = sd(c(X139_Pool_BR5270_ACN_pool_20220712153427,
            X139_Pool_BR5270_ACN_pool_20220712191626,
            X139_Pool_BR5270_ACN_pool_20220713010412), na.rm=TRUE),
  RSD = SD/MEDIA*100,
  COUNT = mean(!is.na(c(X139_Pool_BR5270_ACN_pool_20220712153427,
                        X139_Pool_BR5270_ACN_pool_20220712191626,
                        X139_Pool_BR5270_ACN_pool_20220713010412)))*100,
  BLANK_CONTR = mean(mean(X101_.dry._BR5270_blank_ACN_1, X102_.dry._BR5270_blank_ACN_2)/MEDIA*100)
  )

RPN0_dry_BR5270_ACNf <- RPN0_dry_BR5270_ACN %>%
  filter(COUNT>50) %>%
  filter(RSD<50) %>%
  filter(BLANK_CONTR<50)



# dry_UA57_ACN

RPN0_dry_UA57_ACN <- RPN0 %>%
  mutate(X107_.dry._UA57_blank_ACN_1 = ifelse(is.na(X107_.dry._UA57_blank_ACN_1), 0, X107_.dry._UA57_blank_ACN_1),
         X108_.dry._UA57_blank_ACN_2 = ifelse(is.na(X108_.dry._UA57_blank_ACN_2), 0, X108_.dry._UA57_blank_ACN_2)) %>%
  rowwise() %>%
  mutate(MEDIA = mean(c(X140_Pool_UA57_ACN_pool_20220712155017,
                        X140_Pool_UA57_ACN_pool_20220712193216,
                        X140_Pool_UA57_ACN_pool_20220713012002), na.rm=TRUE),
  SD = sd(c(X140_Pool_UA57_ACN_pool_20220712155017,
            X140_Pool_UA57_ACN_pool_20220712193216,
            X140_Pool_UA57_ACN_pool_20220713012002), na.rm=TRUE),
  RSD = SD/MEDIA*100,
  COUNT = mean(!is.na(c(X140_Pool_UA57_ACN_pool_20220712155017,
                        X140_Pool_UA57_ACN_pool_20220712193216,
                        X140_Pool_UA57_ACN_pool_20220713012002)))*100,
  BLANK_CONTR = mean(mean(X107_.dry._UA57_blank_ACN_1, X108_.dry._UA57_blank_ACN_2)/MEDIA*100)
  )

RPN0_dry_UA57_ACNf <- RPN0_dry_UA57_ACN %>%
  filter(COUNT>50) %>%
  filter(RSD<50) %>%
  filter(BLANK_CONTR<50)



# dry_VC40_ACN_bis


RPN0_dry_VC40_ACN_bis <- RPN0 %>%
  mutate(X125_.dry._VC40_blank_ACN_1 = ifelse(is.na(X125_.dry._VC40_blank_ACN_1), 0, X125_.dry._VC40_blank_ACN_1),
         X126_.dry._VC40_blank_ACN_2 = ifelse(is.na(X126_.dry._VC40_blank_ACN_2), 0, X126_.dry._VC40_blank_ACN_2)) %>%
  rowwise() %>%
  mutate(MEDIA = mean(c(X143_Pool_VC40_ACN_pool_20220712160609,
                        X143_Pool_VC40_ACN_pool_20220712194808,
                        X143_Pool_VC40_ACN_pool_20220713013553), na.rm=TRUE),
  SD = sd(c(X143_Pool_VC40_ACN_pool_20220712160609,
            X143_Pool_VC40_ACN_pool_20220712194808,
            X143_Pool_VC40_ACN_pool_20220713013553), na.rm=TRUE),
  RSD = SD/MEDIA*100,
  COUNT = mean(!is.na(c(X143_Pool_VC40_ACN_pool_20220712160609,
                        X143_Pool_VC40_ACN_pool_20220712194808,
                        X143_Pool_VC40_ACN_pool_20220713013553)))*100,
  BLANK_CONTR = mean(mean(X125_.dry._VC40_blank_ACN_1, X126_.dry._VC40_blank_ACN_2)/MEDIA*100)
  )

RPN0_dry_VC40_ACN_bisf <- RPN0_dry_VC40_ACN_bis %>%
  filter(COUNT>50) %>%
  filter(RSD<50) %>%
  filter(BLANK_CONTR<50)


#  dry_2U!!!


# dry_N2_2U

RPN0_dry_N2_2U <- RPN0 %>%
  mutate(X227_.dry._N2_blank_2U_1 = ifelse(is.na(X227_.dry._N2_blank_2U_1), 0, X227_.dry._N2_blank_2U_1),
         X228_.dry._N2_blank_2U_2 = ifelse(is.na(X228_.dry._N2_blank_2U_2), 0, X228_.dry._N2_blank_2U_2)) %>%
  rowwise() %>%
  mutate(MEDIA = mean(c(X280_.dry._Pool_N2_2U_pool_20220713072400,
                        X280_.dry._Pool_N2_2U_pool_20220713110534,
                        X280_.dry._Pool_N2_2U_pool_20220713165311), na.rm=TRUE),
  SD = sd(c(X280_.dry._Pool_N2_2U_pool_20220713072400,
            X280_.dry._Pool_N2_2U_pool_20220713110534,
            X280_.dry._Pool_N2_2U_pool_20220713165311), na.rm=TRUE),
  RSD = SD/MEDIA*100,
  COUNT = mean(!is.na(c(X280_.dry._Pool_N2_2U_pool_20220713072400,
                        X280_.dry._Pool_N2_2U_pool_20220713110534,
                        X280_.dry._Pool_N2_2U_pool_20220713165311)))*100,
  BLANK_CONTR = mean(mean(X227_.dry._N2_blank_2U_1, X228_.dry._N2_blank_2U_2)/MEDIA*100)
  )

RPN0_dry_N2_2Uf <- RPN0_dry_N2_2U %>%
  filter(COUNT>50) %>%
  filter(RSD<50) %>%
  filter(BLANK_CONTR<50)



# dry_VC40_2U

RPN0_dry_VC40_2U <- RPN0 %>%
  mutate(X233_.dry._VC40_blank_2U_1 = ifelse(is.na(X233_.dry._VC40_blank_2U_1), 0, X233_.dry._VC40_blank_2U_1),
         X234_.dry._VC40_blank_2U_2 = ifelse(is.na(X234_.dry._VC40_blank_2U_2), 0, X234_.dry._VC40_blank_2U_2)) %>%
  rowwise() %>%
  mutate(MEDIA = mean(c(X281_.dry._Pool_VC40_2U_pool_20220713073952,
                        X281_.dry._Pool_VC40_2U_pool_20220713112126,
                        X281_.dry._Pool_VC40_2U_pool_20220713170903), na.rm=TRUE),
  SD = sd(c(X281_.dry._Pool_VC40_2U_pool_20220713073952,
            X281_.dry._Pool_VC40_2U_pool_20220713112126,
            X281_.dry._Pool_VC40_2U_pool_20220713170903), na.rm=TRUE),
  RSD = SD/MEDIA*100,
  COUNT = mean(!is.na(c(X281_.dry._Pool_VC40_2U_pool_20220713073952,
                        X281_.dry._Pool_VC40_2U_pool_20220713112126,
                        X281_.dry._Pool_VC40_2U_pool_20220713170903)))*100,
  BLANK_CONTR = mean(mean(X233_.dry._VC40_blank_2U_1, X234_.dry._VC40_blank_2U_2)/MEDIA*100)
  )

RPN0_dry_VC40_2Uf <- RPN0_dry_VC40_2U %>%
  filter(COUNT>50) %>%
  filter(RSD<50) %>%
  filter(BLANK_CONTR<50)



# dry_VC1668_2U

RPN0_dry_VC1668_2U <- RPN0 %>%
  mutate(X239_.dry._VC1668_blank_2U_1 = ifelse(is.na(X239_.dry._VC1668_blank_2U_1), 0, X239_.dry._VC1668_blank_2U_1),
         X240_.dry._VC1668_blank_2U_2 = ifelse(is.na(X240_.dry._VC1668_blank_2U_2), 0, X240_.dry._VC1668_blank_2U_2)) %>%
  rowwise() %>%
  mutate(MEDIA = mean(c(X282_.dry._Pool_VC1668_2U_pool_20220713075544,
                        X282_.dry._Pool_VC1668_2U_pool_20220713113716,
                        X282_.dry._Pool_VC1668_2U_pool_20220713172455), na.rm=TRUE),
  SD = sd(c(X282_.dry._Pool_VC1668_2U_pool_20220713075544,
            X282_.dry._Pool_VC1668_2U_pool_20220713113716,
            X282_.dry._Pool_VC1668_2U_pool_20220713172455), na.rm=TRUE),
  RSD = SD/MEDIA*100,
  COUNT = mean(!is.na(c(X282_.dry._Pool_VC1668_2U_pool_20220713075544,
                        X282_.dry._Pool_VC1668_2U_pool_20220713113716,
                        X282_.dry._Pool_VC1668_2U_pool_20220713172455)))*100,
  BLANK_CONTR = mean(mean(X239_.dry._VC1668_blank_2U_1, X240_.dry._VC1668_blank_2U_2)/MEDIA*100)
  )

RPN0_dry_VC1668_2Uf <- RPN0_dry_VC1668_2U %>%
  filter(COUNT>50) %>%
  filter(RSD<50) %>%
  filter(BLANK_CONTR<50)



# dry_BR5270_2U

RPN0_dry_BR5270_2U <- RPN0 %>%
  mutate(X245_.dry._BR5270_blank_2U_1 = ifelse(is.na(X245_.dry._BR5270_blank_2U_1), 0, X245_.dry._BR5270_blank_2U_1),
         X246_.dry._BR5270_blank_2U_2 = ifelse(is.na(X246_.dry._BR5270_blank_2U_2), 0, X246_.dry._BR5270_blank_2U_2)) %>%
  rowwise() %>%
  mutate(MEDIA = mean(c(X283_.dry._Pool_BR5270_2U_pool_20220713081134,
                        X283_.dry._Pool_BR5270_2U_pool_20220713115308,
                        X283_.dry._Pool_BR5270_2U_pool_20220713174046), na.rm=TRUE),
  SD = sd(c(X283_.dry._Pool_BR5270_2U_pool_20220713081134,
            X283_.dry._Pool_BR5270_2U_pool_20220713115308,
            X283_.dry._Pool_BR5270_2U_pool_20220713174046), na.rm=TRUE),
  RSD = SD/MEDIA*100,
  COUNT = mean(!is.na(c(X283_.dry._Pool_BR5270_2U_pool_20220713081134,
                        X283_.dry._Pool_BR5270_2U_pool_20220713115308,
                        X283_.dry._Pool_BR5270_2U_pool_20220713174046)))*100,
  BLANK_CONTR = mean(mean(X245_.dry._BR5270_blank_2U_1, X246_.dry._BR5270_blank_2U_2)/MEDIA*100)
  )

RPN0_dry_BR5270_2Uf <- RPN0_dry_BR5270_2U %>%
  filter(COUNT>50) %>%
  filter(RSD<50) %>%
  filter(BLANK_CONTR<50)



# dry_UA57_2U

RPN0_dry_UA57_2U <- RPN0 %>%
  mutate(X251_.dry._UA57_blank_2U_1 = ifelse(is.na(X251_.dry._UA57_blank_2U_1), 0, X251_.dry._UA57_blank_2U_1),
         X252_.dry._UA57_blank_2U_2 = ifelse(is.na(X252_.dry._UA57_blank_2U_2), 0, X252_.dry._UA57_blank_2U_2)) %>%
  rowwise() %>%
  mutate(MEDIA = mean(c(X284_.dry._Pool_UA57_2U_pool_20220713082725,
                        X284_.dry._Pool_UA57_2U_pool_20220713120900,
                        X284_.dry._Pool_UA57_2U_pool_20220713175638), na.rm=TRUE),
  SD = sd(c(X284_.dry._Pool_UA57_2U_pool_20220713082725,
            X284_.dry._Pool_UA57_2U_pool_20220713120900,
            X284_.dry._Pool_UA57_2U_pool_20220713175638), na.rm=TRUE),
  RSD = SD/MEDIA*100,
  COUNT = mean(!is.na(c(X284_.dry._Pool_UA57_2U_pool_20220713082725,
                        X284_.dry._Pool_UA57_2U_pool_20220713120900,
                        X284_.dry._Pool_UA57_2U_pool_20220713175638)))*100,
  BLANK_CONTR = mean(mean(X251_.dry._UA57_blank_2U_1, X252_.dry._UA57_blank_2U_2)/MEDIA*100)
  )

RPN0_dry_UA57_2Uf <- RPN0_dry_UA57_2U %>%
  filter(COUNT>50) %>%
  filter(RSD<50) %>%
  filter(BLANK_CONTR<50)



# dry_VC40_2U_bis

RPN0_dry_VC40_2U_bis <- RPN0 %>%
  mutate(X269_.dry._VC40_blank_2U_1 = ifelse(is.na(X269_.dry._VC40_blank_2U_1), 0, X269_.dry._VC40_blank_2U_1),
         X270_.dry._VC40_blank_2U_2 = ifelse(is.na(X270_.dry._VC40_blank_2U_2), 0, X270_.dry._VC40_blank_2U_2)) %>%
  rowwise() %>%
  mutate(MEDIA = mean(c(X287_.dry._Pool_VC40_2U_pool_20220713084317,
                        X287_.dry._Pool_VC40_2U_pool_20220713122450,
                        X287_.dry._Pool_VC40_2U_pool_20220713181230), na.rm=TRUE),
  SD = sd(c(X287_.dry._Pool_VC40_2U_pool_20220713084317,
            X287_.dry._Pool_VC40_2U_pool_20220713122450,
            X287_.dry._Pool_VC40_2U_pool_20220713181230), na.rm=TRUE),
  RSD = SD/MEDIA*100,
  COUNT = mean(!is.na(c(X287_.dry._Pool_VC40_2U_pool_20220713084317,
                        X287_.dry._Pool_VC40_2U_pool_20220713122450,
                        X287_.dry._Pool_VC40_2U_pool_20220713181230)))*100,
  BLANK_CONTR = mean(mean(X269_.dry._VC40_blank_2U_1, X270_.dry._VC40_blank_2U_2)/MEDIA*100)
  )

RPN0_dry_VC40_2U_bisf <- RPN0_dry_VC40_2U_bis %>%
  filter(COUNT>50) %>%
  filter(RSD<50) %>%
  filter(BLANK_CONTR<50)




### merge!!

RPN_unique_feat_vector <- unique(c(RPN0_dry_N2_ACNf$group, RPN0_dry_VC40_ACNf$group, RPN0_dry_VC1668_ACNf$group, RPN0_dry_BR5270_ACNf$group, RPN0_dry_UA57_ACNf$group, RPN0_dry_VC40_ACN_bisf$group,
                                   RPN0_dry_N2_2Uf$group, RPN0_dry_VC40_2Uf$group, RPN0_dry_VC1668_2Uf$group, RPN0_dry_BR5270_2Uf$group, RPN0_dry_UA57_2Uf$group, RPN0_dry_VC40_2U_bisf$group))


RPN1 <- filter(RPN0, group %in% RPN_unique_feat_vector)

write_tsv(RPN1, "RPN1.txt")



#2022-10-04
## some visualizations of the results
### to me it seems more convenient to put samples in rows, and features in columns!!



RPN1_t <- RPN1 %>%
  select(-ret, -mz) %>%
  t()

RPN1_t_first_col <- matrix(data = rownames(RPN1_t), nrow = nrow(RPN1_t), ncol = 1)
RPN1_t <- cbind(RPN1_t_first_col, RPN1_t)
colnames(RPN1_t) <- c("ANALYSIS", as.vector(RPN1_t[1,])[-1])

RPN1_tc <- RPN1_t[-1,]
rownames(RPN1_tc) <- NULL

RPN2 <- as_tibble(RPN1_tc) %>%
  mutate_at(colnames(RPN1_tc)[-1], as.numeric)

write_tsv(RPN2, "RPN2.txt")








RPN3 <- RPN2 %>%
  mutate(CORRESP_FOUND = ifelse(strsplit(ANALYSIS, split = "_")[[1]][1] %in% sample_info$Analyses_name, TRUE, FALSE),
         .after = ANALYSIS)




if (mean(RPN3$CORRESP_FOUND) == 1) {
  for (a in colnames(sample_info)) {
    RPN3 <- mutate(RPN3, NEWCOLOUMN = rep(NA, length(RPN3$ANALYSIS)), .before = CORRESP_FOUND)
    colnames(RPN3)[colnames(RPN3) == "NEWCOLOUMN"] <- a
    for (ra in 1:length(RPN3$ANALYSIS)) {
      for (rs in 1:length(sample_info$Analyses_name)) {
        if (sample_info$Analyses_name[rs] == strsplit(RPN3$ANALYSIS, split = "_")[[ra]][1]) {
          RPN3[ra, a] <- sample_info[rs, a]
        }
      }
    }
  }
} else {
  stop("There are some analyses of which name is not included in the list of shipped samples!")
}

write_tsv(RPN3, "RPN3.txt")



## a check
#samples not included in the sequence of analyses:
samples_not_in_the_RPN_sequence <- filter(sample_info_dry, Analyses_name %in% sample_info_dry$Analyses_name[!sample_info_dry$Analyses_name %in% RPN3$Analyses_name])




### then, I would do Eulero-venn
### 2023-03-31 Eulero-Venn

#Eulero-venn among the dry, comparing scheme 1 vs 2:

RPN_feat_dry_ACN <- unique(c(RPN0_dry_N2_ACNf$group, RPN0_dry_VC40_ACNf$group, RPN0_dry_VC1668_ACNf$group, RPN0_dry_BR5270_ACNf$group, RPN0_dry_UA57_ACNf$group, RPN0_dry_VC40_ACN_bisf$group))
RPN_feat_dry_2U <- unique(c(RPN0_dry_N2_2Uf$group, RPN0_dry_VC40_2Uf$group, RPN0_dry_VC1668_2Uf$group, RPN0_dry_BR5270_2Uf$group, RPN0_dry_UA57_2Uf$group, RPN0_dry_VC40_2U_bisf$group))


RPN_EV_only_dry_df <- tibble(feature = RPN_unique_feat_vector,
                             scheme1 = RPN_unique_feat_vector %in% RPN_feat_dry_ACN,
                             scheme2upper = RPN_unique_feat_vector %in% RPN_feat_dry_2U)


jpeg(file="RPLC_NEG_EuleroVenn_only_dry.jpeg")
plot(euler(RPN_EV_only_dry_df[, 2:3], shape = "ellipse"), quantities = TRUE)
dev.off()




#Eulero-venn combined:

RPN_EV_df_QCfil <- tibble(feature = RPN_unique_feat_vector,
                          dry_scheme1 = RPN_unique_feat_vector %in% RPN_feat_dry_ACN,
                          dry_scheme2upper = RPN_unique_feat_vector %in% RPN_feat_dry_2U)

jpeg(file="RPLC_NEG_EuleroVenn.jpeg")
plot(euler(RPN_EV_df_QCfil[, 2:3], shape = "ellipse"), quantities = TRUE)
dev.off()


## better:
RPN_EV_df_QCfil_better <- RPN_EV_df_QCfil
colnames(RPN_EV_df_QCfil_better) <- c("feature", "scheme1", "scheme2")

jpeg(file="RPLC_NEG_EuleroVenn_only_dry_bis.jpeg")
plot(euler(RPN_EV_df_QCfil_better[, 2:3], shape = "ellipse"), fill = c("cadetblue1", "indianred1"), quantities = list(cex = 2.5), labels = list(cex = 2))
dev.off()


#
###
#################### tables for the ANOVA::

RPN4 <- filter(RPN3, Sample_ID != "methBK")

RPN4$Sample_ID <- droplevels(RPN4$Sample_ID)



RPN4_scheme1_dry <- filter(RPN4, Scheme == "scheme1", rawdried == "DRIED", Content == "worm") %>%
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


RPN4_scheme2_dry <- filter(RPN4, Scheme == "scheme2", Content == "worm") %>%
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




RPN4_scheme1_dry_transf <- transf_feat_intensities(RPN4_scheme1_dry)

RPN4_scheme1_dry_ANOVA <- getting_ANOVA_results(RPN4_scheme1_dry_transf)



RPN4_scheme2_dry_transf <- transf_feat_intensities(RPN4_scheme2_dry)

RPN4_scheme2_dry_ANOVA <- getting_ANOVA_results(RPN4_scheme2_dry_transf)



write_tsv(RPN4_scheme1_dry_transf, "RPN4_scheme1_dry_transf.txt")
write_tsv(RPN4_scheme1_dry_ANOVA, "RPN4_scheme1_dry_ANOVA.txt")
write_tsv(RPN4_scheme2_dry_transf, "RPN4_scheme2_dry_transf.txt")
write_tsv(RPN4_scheme2_dry_ANOVA, "RPN4_scheme2_dry_ANOVA.txt")


## only ANOVA significant:


RPN4_scheme1_dry_ANOVA_sign <- filter(RPN4_scheme1_dry_ANOVA, ANOVA_p_value_FDR < 0.05)
RPN4_scheme2_dry_ANOVA_sign <- filter(RPN4_scheme2_dry_ANOVA, ANOVA_p_value_FDR < 0.05)


write_tsv(RPN4_scheme1_dry_ANOVA_sign, "RPN4_scheme1_dry_ANOVA_sign.txt")
write_tsv(RPN4_scheme2_dry_ANOVA_sign, "RPN4_scheme2_dry_ANOVA_sign.txt")



## Eulero-Venn with significant features:

RPN_EV_df_QCfil_sign <- tibble(feature = RPN_unique_feat_vector,
                               dry_scheme1 = RPN_unique_feat_vector %in% RPN4_scheme1_dry_ANOVA_sign$feature,
                               dry_scheme2upper = RPN_unique_feat_vector %in% RPN4_scheme2_dry_ANOVA_sign$feature)

RPN_EV_df_QCfil_sign_better <- RPN_EV_df_QCfil_sign
colnames(RPN_EV_df_QCfil_sign_better) <- c("feature", "scheme1", "scheme2")


jpeg(file="RPLC_NEG_EuleroVenn_ANOVA.jpeg")
plot(euler(RPN_EV_df_QCfil_sign_better[, 2:3], shape = "ellipse"), fill = c("cadetblue1", "indianred1"), quantities = list(cex = 2.5), labels = list(cex = 2))
dev.off()



###
# combining significant features with top candidate in annotation

MFsummary <- read_tsv("MFsummary_PCL_no_formula_updt.txt")


MFsummary_colnames <- setNames(rep("", length(colnames(MFsummary))), colnames(MFsummary))
MFsummary_empty <- bind_rows(MFsummary_colnames)[0, ]




MFsummary_empty_for_RPN4_scheme1_dry_ANOVA_sign_annotation <- MFsummary_empty
MFsummary_empty_for_RPN4_scheme1_dry_ANOVA_sign_annotation[1:length(RPN4_scheme1_dry_ANOVA_sign$feature),] <- NA

RPN4_scheme1_dry_ANOVA_sign_annotation <- cbind(RPN4_scheme1_dry_ANOVA_sign, MFsummary_empty_for_RPN4_scheme1_dry_ANOVA_sign_annotation)

for (s in 1:length(RPN4_scheme1_dry_ANOVA_sign_annotation$feature)) {
  if (RPN4_scheme1_dry_ANOVA_sign_annotation$feature[s] %in% MFsummary$group) {
    RPN4_scheme1_dry_ANOVA_sign_annotation[s,7:39] <- filter(MFsummary, group == RPN4_scheme1_dry_ANOVA_sign_annotation$feature[s])[1,]
  }
}

write_tsv(RPN4_scheme1_dry_ANOVA_sign_annotation, "RPN4_scheme1_dry_ANOVA_sign_annotation.txt")




MFsummary_empty_for_RPN4_scheme2_dry_ANOVA_sign_annotation <- MFsummary_empty
MFsummary_empty_for_RPN4_scheme2_dry_ANOVA_sign_annotation[1:length(RPN4_scheme2_dry_ANOVA_sign$feature),] <- NA

RPN4_scheme2_dry_ANOVA_sign_annotation <- cbind(RPN4_scheme2_dry_ANOVA_sign, MFsummary_empty_for_RPN4_scheme2_dry_ANOVA_sign_annotation)

for (s in 1:length(RPN4_scheme2_dry_ANOVA_sign_annotation$feature)) {
  if (RPN4_scheme2_dry_ANOVA_sign_annotation$feature[s] %in% MFsummary$group) {
    RPN4_scheme2_dry_ANOVA_sign_annotation[s,7:39] <- filter(MFsummary, group == RPN4_scheme2_dry_ANOVA_sign_annotation$feature[s])[1,]
  }
}

write_tsv(RPN4_scheme2_dry_ANOVA_sign_annotation, "RPN4_scheme2_dry_ANOVA_sign_annotation.txt")






#
##
###### annotation levels, considering the full list of features:


MFsummary_topMONA_empty <- MFsummary[0,]
colnames(MFsummary_topMONA_empty) <- paste0(colnames(MFsummary), "_topMONA")
MFsummary_topMONA_empty[1:length(RPN_feat_info$group),] <- NA

MFsummary_topScore_empty <- MFsummary[0,]
colnames(MFsummary_topScore_empty) <- paste0(colnames(MFsummary), "_topScore")
MFsummary_topScore_empty[1:length(RPN_feat_info$group),] <- NA


RPN_feat_info_annotation <- cbind(mutate(RPN_feat_info, AnnoLevel = as.character(NA), is_topMONA_topScore = as.logical(NA)), MFsummary_topMONA_empty, MFsummary_topScore_empty)


for (i in 1:length(RPN_feat_info$group)) {
  if (RPN_feat_info$group[i] %in% MFsummary$group) {
    MFsummary_fil <- filter(MFsummary, group == RPN_feat_info$group[i])
    topMONA_row <- which.max(MFsummary_fil$individualMoNAScore)
    topScore_row <- which.max(MFsummary_fil$score)
    
    RPN_feat_info_annotation[i,paste0(colnames(MFsummary), "_topMONA")] <- MFsummary_fil[topMONA_row,colnames(MFsummary_fil)]
    RPN_feat_info_annotation[i,paste0(colnames(MFsummary), "_topScore")] <- MFsummary_fil[topScore_row,colnames(MFsummary_fil)]
    
    RPN_feat_info_annotation[i, "is_topMONA_topScore"] <- RPN_feat_info_annotation[i, "identifier_topMONA"] == RPN_feat_info_annotation[i, "identifier_topScore"]
    
    RPN_feat_info_annotation[i, "AnnoLevel"] <- ifelse(RPN_feat_info_annotation[i,"individualMoNAScore_topMONA"]>0.9, "2a",
                                                       ifelse(RPN_feat_info_annotation[i,"individualMoNAScore_topMONA"]<0.9 & RPN_feat_info_annotation[i,"individualMoNAScore_topMONA"]>0.7, "3a",
                                                              ifelse(RPN_feat_info_annotation[i,"individualMoNAScore_topMONA"]<0.7 & RPN_feat_info_annotation[i,"individualMoNAScore_topMONA"]>0.4, "3b",
                                                                     ifelse(mean(MFsummary_fil$neutral_formula == MFsummary_fil$neutral_formula[1]) == 1, "4a", "5"))))
    
  } else {
    RPN_feat_info_annotation[i, "AnnoLevel"] <- "5"
  }
}

write_tsv(RPN_feat_info_annotation, "RPLCNEG_patRoon_IPO_PCL_feat_info_annotation.txt")


## annotation levels, after filtering QCs:

RPN_feat_info_annotation_QCfiltered <- filter(RPN_feat_info_annotation, group %in% RPN_unique_feat_vector) 

write_tsv(RPN_feat_info_annotation_QCfiltered, "RPLCNEG_patRoon_IPO_PCL_feat_info_annotation_QCfiltered.txt")



## annotation levels, after filtering QCs, per sample scheme:


RPN_feat_info_annotation_QCfiltered_scheme1_dry <- filter(RPN_feat_info_annotation, group %in% RPN_feat_dry_ACN) 

write_tsv(RPN_feat_info_annotation_QCfiltered_scheme1_dry, "RPLCNEG_patRoon_IPO_PCL_feat_info_annotation_QCfiltered_scheme1_dry.txt")



RPN_feat_info_annotation_QCfiltered_scheme2_dry <- filter(RPN_feat_info_annotation, group %in% RPN_feat_dry_2U) 

write_tsv(RPN_feat_info_annotation_QCfiltered_scheme2_dry, "RPLCNEG_patRoon_IPO_PCL_feat_info_annotation_QCfiltered_scheme2_dry.txt")





## annotation levels, only for statistically significant:


RPN_feat_info_annotation_empty_dried_scheme1 <- RPN_feat_info_annotation[0,]
RPN_feat_info_annotation_empty_dried_scheme1[1:length(RPN4_scheme1_dry_ANOVA_sign$feature),] <- NA


RPN4_scheme1_dry_ANOVA_sign_annotat_withLevels <- cbind(RPN4_scheme1_dry_ANOVA_sign, RPN_feat_info_annotation_empty_dried_scheme1)

for (i in 1:length(RPN4_scheme1_dry_ANOVA_sign$feature)) {
  indices <- which(RPN_feat_info_annotation$group == RPN4_scheme1_dry_ANOVA_sign$feature[i])
  
  RPN4_scheme1_dry_ANOVA_sign_annotat_withLevels[i,colnames(RPN_feat_info_annotation_empty_dried_scheme1)] <- RPN_feat_info_annotation[indices,]
  
}

write_tsv(RPN4_scheme1_dry_ANOVA_sign_annotat_withLevels, "RPLCNEG_patRoon_IPO_PCL_scheme1_dry_ANOVA_sign_annot_Levels.txt")




RPN_feat_info_annotation_empty_dried_scheme2upper <- RPN_feat_info_annotation[0,]
RPN_feat_info_annotation_empty_dried_scheme2upper[1:length(RPN4_scheme2_dry_ANOVA_sign$feature),] <- NA


RPN4_scheme2_dry_ANOVA_sign_annotat_withLevels <- cbind(RPN4_scheme2_dry_ANOVA_sign, RPN_feat_info_annotation_empty_dried_scheme2upper)

for (i in 1:length(RPN4_scheme2_dry_ANOVA_sign$feature)) {
  indices <- which(RPN_feat_info_annotation$group == RPN4_scheme2_dry_ANOVA_sign$feature[i])
  
  RPN4_scheme2_dry_ANOVA_sign_annotat_withLevels[i,colnames(RPN_feat_info_annotation_empty_dried_scheme2upper)] <- RPN_feat_info_annotation[indices,]
  
}

write_tsv(RPN4_scheme2_dry_ANOVA_sign_annotat_withLevels, "RPLCNEG_patRoon_IPO_PCL_scheme2_dry_ANOVA_sign_annotat_withLevels.txt")






###
# combining significant features with top candidate in annotation

MFsummary <- read_tsv("MFsummary_WJ_no_formula_bis.txt")


MFsummary_colnames <- setNames(rep("", length(colnames(MFsummary))), colnames(MFsummary))
MFsummary_empty <- bind_rows(MFsummary_colnames)[0, ]




MFsummary_empty_for_RPN4_scheme1_dry_ANOVA_sign_annotation <- MFsummary_empty
MFsummary_empty_for_RPN4_scheme1_dry_ANOVA_sign_annotation[1:length(RPN4_scheme1_dry_ANOVA_sign$feature),] <- NA

RPN4_scheme1_dry_ANOVA_sign_annotation <- cbind(RPN4_scheme1_dry_ANOVA_sign, MFsummary_empty_for_RPN4_scheme1_dry_ANOVA_sign_annotation)

for (s in 1:length(RPN4_scheme1_dry_ANOVA_sign_annotation$feature)) {
  if (RPN4_scheme1_dry_ANOVA_sign_annotation$feature[s] %in% MFsummary$group) {
    RPN4_scheme1_dry_ANOVA_sign_annotation[s,7:39] <- filter(MFsummary, group == RPN4_scheme1_dry_ANOVA_sign_annotation$feature[s])[1,]
  }
}

write_tsv(RPN4_scheme1_dry_ANOVA_sign_annotation, "RPN4_scheme1_dry_ANOVA_sign_annotation.txt")




MFsummary_empty_for_RPN4_scheme2_dry_ANOVA_sign_annotation <- MFsummary_empty
MFsummary_empty_for_RPN4_scheme2_dry_ANOVA_sign_annotation[1:length(RPN4_scheme2_dry_ANOVA_sign$feature),] <- NA

RPN4_scheme2_dry_ANOVA_sign_annotation <- cbind(RPN4_scheme2_dry_ANOVA_sign, MFsummary_empty_for_RPN4_scheme2_dry_ANOVA_sign_annotation)

for (s in 1:length(RPN4_scheme2_dry_ANOVA_sign_annotation$feature)) {
  if (RPN4_scheme2_dry_ANOVA_sign_annotation$feature[s] %in% MFsummary$group) {
    RPN4_scheme2_dry_ANOVA_sign_annotation[s,7:39] <- filter(MFsummary, group == RPN4_scheme2_dry_ANOVA_sign_annotation$feature[s])[1,]
  }
}

write_tsv(RPN4_scheme2_dry_ANOVA_sign_annotation, "RPN4_scheme2_dry_ANOVA_sign_annotation.txt")






#
##
###### annotation levels, considering the full list of features:


MFsummary_topMONA_empty <- MFsummary[0,]
colnames(MFsummary_topMONA_empty) <- paste0(colnames(MFsummary), "_topMONA")
MFsummary_topMONA_empty[1:length(RPN_feat_info$group),] <- NA

MFsummary_topScore_empty <- MFsummary[0,]
colnames(MFsummary_topScore_empty) <- paste0(colnames(MFsummary), "_topScore")
MFsummary_topScore_empty[1:length(RPN_feat_info$group),] <- NA


RPN_feat_info_annotation <- cbind(mutate(RPN_feat_info, AnnoLevel = as.character(NA), is_topMONA_topScore = as.logical(NA)), MFsummary_topMONA_empty, MFsummary_topScore_empty)


for (i in 1:length(RPN_feat_info$group)) {
  if (RPN_feat_info$group[i] %in% MFsummary$group) {
    MFsummary_fil <- filter(MFsummary, group == RPN_feat_info$group[i])
    topMONA_row <- which.max(MFsummary_fil$individualMoNAScore)
    topScore_row <- which.max(MFsummary_fil$score)
    
    RPN_feat_info_annotation[i,paste0(colnames(MFsummary), "_topMONA")] <- MFsummary_fil[topMONA_row,colnames(MFsummary_fil)]
    RPN_feat_info_annotation[i,paste0(colnames(MFsummary), "_topScore")] <- MFsummary_fil[topScore_row,colnames(MFsummary_fil)]
    
    RPN_feat_info_annotation[i, "is_topMONA_topScore"] <- RPN_feat_info_annotation[i, "identifier_topMONA"] == RPN_feat_info_annotation[i, "identifier_topScore"]
    
    RPN_feat_info_annotation[i, "AnnoLevel"] <- ifelse(RPN_feat_info_annotation[i,"individualMoNAScore_topMONA"]>0.9, "2a",
                                                       ifelse(RPN_feat_info_annotation[i,"individualMoNAScore_topMONA"]<0.9 & RPN_feat_info_annotation[i,"individualMoNAScore_topMONA"]>0.7, "3a",
                                                              ifelse(RPN_feat_info_annotation[i,"individualMoNAScore_topMONA"]<0.7 & RPN_feat_info_annotation[i,"individualMoNAScore_topMONA"]>0.4, "3b",
                                                                     ifelse(mean(MFsummary_fil$neutral_formula == MFsummary_fil$neutral_formula[1]) == 1, "4a", "5"))))
    
  } else {
    RPN_feat_info_annotation[i, "AnnoLevel"] <- "5"
  }
}

write_tsv(RPN_feat_info_annotation, "RPLCNEG_patRoon_IPO_WJ_feat_info_annotation.txt")


## annotation levels, after filtering QCs:

RPN_feat_info_annotation_QCfiltered <- filter(RPN_feat_info_annotation, group %in% RPN_unique_feat_vector) 

write_tsv(RPN_feat_info_annotation_QCfiltered, "RPLCNEG_patRoon_IPO_WJ_feat_info_annotation_QCfiltered.txt")





## annotation levels, after filtering QCs, per sample scheme:


RPN_feat_info_annotation_QCfiltered_scheme1_dry <- filter(RPN_feat_info_annotation, group %in% RPN_feat_dry_ACN) 

write_tsv(RPN_feat_info_annotation_QCfiltered_scheme1_dry, "RPLCNEG_patRoon_IPO_WJ_feat_info_annotation_QCfiltered_scheme1_dry.txt")



RPN_feat_info_annotation_QCfiltered_scheme2_dry <- filter(RPN_feat_info_annotation, group %in% RPN_feat_dry_2U) 

write_tsv(RPN_feat_info_annotation_QCfiltered_scheme2_dry, "RPLCNEG_patRoon_IPO_WJ_feat_info_annotation_QCfiltered_scheme2_dry.txt")





## annotation levels, only for statistically significant:


RPN_feat_info_annotation_empty_dried_scheme1 <- RPN_feat_info_annotation[0,]
RPN_feat_info_annotation_empty_dried_scheme1[1:length(RPN4_scheme1_dry_ANOVA_sign$feature),] <- NA


RPN4_scheme1_dry_ANOVA_sign_annotat_withLevels <- cbind(RPN4_scheme1_dry_ANOVA_sign, RPN_feat_info_annotation_empty_dried_scheme1)

for (i in 1:length(RPN4_scheme1_dry_ANOVA_sign$feature)) {
  indices <- which(RPN_feat_info_annotation$group == RPN4_scheme1_dry_ANOVA_sign$feature[i])
  
  RPN4_scheme1_dry_ANOVA_sign_annotat_withLevels[i,colnames(RPN_feat_info_annotation_empty_dried_scheme1)] <- RPN_feat_info_annotation[indices,]
  
}

write_tsv(RPN4_scheme1_dry_ANOVA_sign_annotat_withLevels, "RPLCNEG_patRoon_IPO_WJ_scheme1_dry_ANOVA_sign_annotat_withLevels.txt")




RPN_feat_info_annotation_empty_dried_scheme2upper <- RPN_feat_info_annotation[0,]
RPN_feat_info_annotation_empty_dried_scheme2upper[1:length(RPN4_scheme2_dry_ANOVA_sign$feature),] <- NA


RPN4_scheme2_dry_ANOVA_sign_annotat_withLevels <- cbind(RPN4_scheme2_dry_ANOVA_sign, RPN_feat_info_annotation_empty_dried_scheme2upper)

for (i in 1:length(RPN4_scheme2_dry_ANOVA_sign$feature)) {
  indices <- which(RPN_feat_info_annotation$group == RPN4_scheme2_dry_ANOVA_sign$feature[i])
  
  RPN4_scheme2_dry_ANOVA_sign_annotat_withLevels[i,colnames(RPN_feat_info_annotation_empty_dried_scheme2upper)] <- RPN_feat_info_annotation[indices,]
  
}

write_tsv(RPN4_scheme2_dry_ANOVA_sign_annotat_withLevels, "RPLCNEG_patRoon_IPO_WJ_scheme2_dry_ANOVA_sign_annotat_withLevels.txt")



