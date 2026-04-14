# Running vpmr code on an example dataframe
library(Rcpp)
setwd("M:/R")
sourceCpp("vpmr_cpp_v2.1.cpp")
con <- url("https://raw.githubusercontent.com/tomjspain/vpmr/main/example%20df.rds", "rb")
df <- readRDS(con)
close(con)

# Point estimate
df_kendall <- rec_CIndex(df$ID, df$ObsCount, df$PredCount, type = "kendall")

# Leave-one-out jackknife
df_kendall_jk <- jk_rec_CIndex(df$ID, df$ObsCount, df$PredCount, type = "kendall")

# Grouped jackknife - partition (non-overlapping, group size = 10)
df_kendall_gjk_partition <- gjk_rec_CIndex_partition(df$ID, df$ObsCount, df$PredCount,
                                               type = "kendall", d = 10)

# Grouped jackknife - random overlapping (drop 10 per replicate, 200 replicates)
df_kendall_gjk_overlap <- gjk_rec_CIndex_overlap(df$ID, df$ObsCount, df$PredCount,
                                           type = "kendall", d = 10, n_reps = 200)

# Grouped jackknife - overlapping with pre-drawn drop sets
set.seed(42)
drop_sets <- do.call(rbind, replicate(200, sample(nrow(df), 10), simplify = FALSE))
df_kendall_gjk_overlap_sets <- gjk_rec_CIndex_overlap_sets(df$ID, df$ObsCount, df$PredCount,
                                                     type = "kendall", drop_sets = drop_sets)
