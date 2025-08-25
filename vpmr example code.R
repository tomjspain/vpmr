# Running vpmr code on an example dataframe
library(Rcpp)

setwd("M:/R")
sourceCpp("vpmr_cpp_v1.0.cpp")

con <- url("https://raw.githubusercontent.com/tomjspain/vpmr/main/example df.rds", "rb")
df <- readRDS(con)
close(con)

df_kendall <- rec_CIndex(df$ID, df$ObsCount, df$PredCount, type = "kendall")

df_kendall_jk <- jk_rec_CIndex(df$ID, df$ObsCount, df$PredCount, type = "kendall")

