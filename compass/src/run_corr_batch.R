#!/usr/bin/env Rscript

# args: [1] string tissue

# Description: run correlation for each cell type per tissue
#
# author Sascha Schäuble
# date of creation: Tue Jun 25 12:50:55 2024
# license: MIT

args = commandArgs(trailingOnly=TRUE)

library("tidyverse")
library("magrittr")
library("here")
library("janitor")
library("rstatix")
library("data.table")

TISSUE_STR <- args[1]

PROJ_PATH <- here::here() %>% str_remove("src$")
DAT_PATH0 <- paste0(PROJ_PATH, "dat/modelDat/")
DAT_PATH1 <- paste0(PROJ_PATH, "dat/snRNA/")
DAT_PATH_RES <- paste0(PROJ_PATH, "res/scores/", TISSUE_STR, "/")

FN0 <- "rxns.tsv"
FNRES <- "reactions_transformed.tsv"

DATE_STR <- format(Sys.time(), "%y%m%d")

FN.poi <- "poi.txt"
FN.dat <- "reactions_transformed.tsv" # due to size needs to be generated with Compass
FN.meta <- paste0(TISSUE_STR, "_meta.tsv") # need to unzip meta.zip first with passphrase
FN.Recon2.meta <- "rxns_meta.tsv"

# ================================================================ #


#### Data wrangling ################################################
poi <- read_tsv(paste0(DAT_PATH0, FN.poi), col_names = F) %>% pull()

recon2.meta <- read_tsv(paste0(DAT_PATH0, FN.Recon2.meta)) %>% dplyr::select(rxn_id, subsystem) %>% rename("rxn_ID" = "rxn_id")
recon2.meta %<>% filter(subsystem %in% poi)
recon2.meta$subsystem %>% table()

compass.dat  <- fread(paste0(DAT_PATH_RES, FN.dat), sep = "\t") %>% rename("rxn_ID" = "V1") %>%
  mutate(rxn_ID_base = rxn_ID %>% str_remove("_pos$") %>% str_remove("_neg$")) %>%
  left_join(recon2.meta, by = c("rxn_ID_base" = "rxn_ID")) %>% filter(subsystem %in% poi) %>%
  pivot_longer(-c("rxn_ID", "rxn_ID_base", "subsystem"), names_to = "Cell_ID")
compass.dat %>% glimpse()

meta <- read_tsv(paste0(DAT_PATH1, FN.meta))

compass.dat %<>% left_join(meta, by = "Cell_ID")

compass.dat %<>% select(-rxn_ID_base, -orig.ident) %>%
  mutate(
    rxn_ID = rxn_ID %>% as_factor(),
    Cell_ID = Cell_ID %>% as_factor(),
    Phenotype = Phenotype %>% as_factor(),
    Cell_type = Cell_type %>% as_factor()
  )

compass.dat %<>% 
  dplyr::select(-subsystem) %>%
  pivot_wider(names_from = rxn_ID, values_from = value)


print("compass dim")
compass.dat %>% dim()



# ================================================================ #

#### correlation per CT ############################################

rxnIDs <- colnames(compass.dat %>% dplyr::select(-c("Cell_ID", "Phenotype", "Cell_type")))


print(TISSUE_STR)
estNumCorrs <- (dim(compass.dat)[2]-3)^2 * length(compass.dat$Cell_type %>% levels()) * 2
print(paste("Estimated number of correlations that will be computed:", estNumCorrs))
# compass.dat %>% glimpse()
print("Starting correlation analysis...")
start_time <- Sys.time()

corr <- list()

for (ph in (compass.dat$Phenotype %>% levels())) {
  # compute separately for sham and CLP
  for (ct_i in (compass.dat$Cell_type %>% levels())) {
    # compute separately within and across CT
    print(paste(ph, ct_i, sep = "_"))
    corr[[paste(ph, ct_i, sep = "_")]] <- cor(
      x = compass.dat %>% filter(Phenotype == ph &
                                   Cell_type == ct_i) %>% dplyr::select(rxnIDs),
      method = "sp"
    ) %>% as_tibble(rownames = "rxn_IDi") %>% pivot_longer(-rxn_IDi, names_to = "rxn_IDj") %>% filter(rxn_IDi != rxn_IDj)
  }
}
corr <- rbindlist(corr, idcol = "condition_CT") %>% drop_na()
corr %<>% mutate(
  condition = condition_CT %>% str_extract("Sham|CLP") %>% as_factor(),
  Cell_type = condition_CT %>% str_remove("Sham_|CLP_") %>% as_factor()
  ) %>% 
    relocate(condition, .before = "condition_CT" ) %>%
    relocate(Cell_type, .before = "condition" ) %>% dplyr::select(-condition_CT)

end_time <- Sys.time()

end_time - start_time

corr %>% write_tsv(paste0(DAT_PATH_RES, TISSUE_STR, "_", "corr_", DATE_STR, ".tsv"))
# ================================================================ #
