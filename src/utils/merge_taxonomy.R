setwd("~/stage_lied/Data/")

JGI = read.csv("JGI/taxonomy2_JGI.csv", header = FALSE)
NCBI = read.csv("NCBI/taxonomy2_NCBI.csv", header = FALSE)

library(stringr)
library(fuzzyjoin)
#colnames(JGI)

JGI$V7 = as.character(JGI$V7)
NCBI$V7 = as.character(NCBI$V7)

merged = JGI %>% fuzzy_inner_join(NCBI, by = c("V7" = "V7"), match_fun = str_detect)

new_JGI = JGI[-which(JGI[,7]%in%merged[,7]),]

#new_JGI = new_JGI[,-c(8)]

#colnames(new_JGI) = c("V1","V2","V3","V4","V5","V6","V7","V8")
#join_df = rbind(new_JGI[,-c(8)], NCBI)
write.csv(new_JGI, "JGI_purged.csv")
