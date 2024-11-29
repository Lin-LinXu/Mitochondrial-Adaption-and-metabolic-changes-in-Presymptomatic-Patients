#!/usr/bin/env Rscript

# args: [1] string tissue

# Description: transform penalty scores to positive scores
#
# author Sascha Schäuble
# date of creation: Mon Jul  1 17:31:23 2024
# license: MIT

args = commandArgs(trailingOnly=TRUE)

TISSUE_STR <- args[1]

print(args)

library("tidyverse")
library("magrittr")
library("here")
library("data.table")
library("foreach")
library("doParallel")

PROJ_PATH <- here::here() %>% str_remove("src$")
DAT_PATH_RES <- paste0("../res/scores/", TISSUE_STR)
RES_PATH <- DAT_PATH_RES

DATE_STR <- format(Sys.time(), "%y%m%d")

FN <- "reactions.tsv" # Compass result file, needs to be generated and is not in github due to size
FN.res <- "reactions_transformed.tsv"
# ================================================================ #

#### Data wrangling ################################################
print("Read in data")
dat <- fread(file = paste0(DAT_PATH_RES, FN), showProgress = T)
dat.names <- dat %>% dplyr::select(V1)
dat <- dat %>% select(-V1)
print("Transform...")

log_transformed_df <- -log1p(dat)

min_vec <- apply(log_transformed_df, 1, min, simplify = TRUE)

print("Subtract min...")

min_subtracted_df <- log_transformed_df-min_vec

results <- cbind(dat.names, min_subtracted_df)
print(paste0("write output to ", RES_PATH, FN.res))
data.table::fwrite(results, file = paste0(RES_PATH, FN.res), sep = "\t")

# ================================================================ #
