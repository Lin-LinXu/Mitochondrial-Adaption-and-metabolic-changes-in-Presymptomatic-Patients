#!/usr/bin/env Rscript

# args: [1] string tissue
args = commandArgs(trailingOnly=TRUE)

print(args)

# Description: run statistical analysis for compass output
#
# author Sascha Schäuble
# date of creation: Tue July 08 14:12:26 2024
# license: MIT

library("tidyverse")
library("magrittr")
library("here")
library("janitor")
library("ggpubr")
library("data.table")
library("rstatix")
library("lsr")
library("foreach")
library("doParallel")



TISSUE_STR <- args[1]

PROJ_PATH <- here::here() %>% str_remove("src$")
DAT_PATH0 <- paste0(PROJ_PATH, "dat/modelDat/")
DAT_PATH1 <- paste0(PROJ_PATH, "dat/snRNA/")
RES_PATH <- paste0(PROJ_PATH, "res/scores/", TISSUE_STR, "/")

FN0 <- "rxns.tsv"
FN1 <- "poi.txt"
FN2 <- paste0(TISSUE_STR, "_meta.tsv")
FNRES <- "reactions_transformed.tsv"
TMPOUT <- "rxns_tmp.tsv"

DATE_STR <- format(Sys.time(), "%y%m%d")

THREADS <- detectCores()-2


print("Using the following path variables:")
print(RES_PATH)
print(paste0("Number of threads allocated: ", THREADS))
# ================================================================ #

#'
#' @description: run per cell type
#'
#' @param DF
#'
#' @return list with results and error
runStats <- function(PATH = DAT_PATH_RES, FN = FNRES, CELL_TYPE, PARALLEL = F, NUMCPU = THREADS) {
  
  log.rxn.err <- c()
  results <- list()
  
  if (PARALLEL) {
    registerDoParallel(NUMCPU)
    started.at = proc.time()
    dat.tmp <- fread(
      file = paste0(PATH, FN),
      select = c("V1", (
        dat.meta %>%
          filter(Cell_type == CELL_TYPE) %>%
          pull(Cell_ID)
      )),
      nThread = NUMCPU
    ) 
    
    dat.tmp.long <-
      dat.tmp %>% pivot_longer(-V1, names_to = "Cell_ID") %>%
      left_join(dat.meta, by = "Cell_ID")
    rm(dat.tmp)
    timetaken(started.at)
    started.at = proc.time()
    results <-
      foreach (
        r = 1:length(dat.m.rxns.filt$rxn_id),
        .errorhandling = 'pass',
        .verbose = T
      ) %dopar% {
        cbind(
          (
          dat.tmp.long %>% group_by(V1) %>%
          filter(V1 == dat.m.rxns.filt$rxn_id[r]) %>%
          rstatix::wilcox_test(formula = value ~ Phenotype, detailed = T)
          ), dat.tmp.long %>% group_by(V1) %>%
            filter(V1 == dat.m.rxns.filt$rxn_id[r]) %>%
            rstatix::cohens_d(formula = value ~ Phenotype) %>% dplyr::select(effsize,magnitude)
        )
      }
    timetaken(started.at)
    log.rxn.err <-
      dat.m.rxns.filt$rxn_id[lengths(results) == 5]
    results <-
      rbindlist(results[lengths(results) == 15],
                idcol = "rxnID")
  } else {
    # single core processing
    
    dat.tmp <- fread(file = paste0(PATH, FN),
                     select = c("V1", (
                       dat.meta %>% filter(Cell_type == CELL_TYPE) %>% pull(Cell_ID)
                     )))
    timetaken(started.at)
    dat.tmp.long <-
      dat.tmp %>% pivot_longer(-V1, names_to = "Cell_ID") %>% left_join(dat.meta, by = "Cell_ID")
    
    pb = txtProgressBar(
      min = 0,
      max = length(dat.m.rxns.filt$rxn_id),
      initial = 0,
      style = 3
    )
    started.at = proc.time()
    for (r in 1:length(dat.m.rxns.filt$rxn_id)) {
      skip_to_next <- F
      tryCatch(
        results[[r]] <- dat.tmp.long %>% group_by(V1) %>%
          filter(V1 == dat.m.rxns.filt$rxn_id[r]) %>%
          rstatix::wilcox_test(formula = value ~ Phenotype, detailed = T),
        error = function(e) {
          skip_to_next <<- T
        }
      )
      setTxtProgressBar(pb, r)
      
      if (skip_to_next) {
        log.rxn.err <- c(log.rxn.err, dat.m.rxns.filt$rxn_id[r])
        next
      }
    }
    results <- rbindlist(results, idcol = "rxnID")
    print(results %>% colnames())
    log.rxn.err <- log.rxn.err %>% as.character()
    close(pb)
  }
  
  return(list("results" = results, "errorRxns" = log.rxn.err))
  
}
# ================================================================ #

#### Data wrangling ################################################
print(paste("Registering", THREADS, "threads..."))
print("===")
print("Step 0: Reading in data...")

# tissue sample - cell type association
dat.meta <- read_tsv(paste0(DAT_PATH1, FN2))
dat.meta$Cell_type %<>% as_factor() 
dat.meta$Cell_type %>% table()

cellTypes <- dat.meta$Cell_type %>% as_factor() %>% levels()

# get reactions per pathway
dat.m.rxns <- read_tsv(paste0(DAT_PATH0, FN0), col_names = F)
dat.m.subs <- read_tsv(paste0(DAT_PATH0, FN1), col_names = F)
dat.m.rxns.filt <- dat.m.rxns %>% filter(X2 %in% dat.m.subs$X1)

colnames(dat.m.rxns.filt) <- c("rxn_id", "subsystem")
rm(dat.m.rxns, dat.m.subs)

# create result lists
res <- list()
res.err <- list()
# ================================================================ #


#### get stats #####################################################
print(as.logical(args[4]))
if (as.logical(args[4])) {
  print("===")
  print("Step 1: Computing stats...")
  # cellTypes
  for (ct in cellTypes) {
    print(ct)
    print(dat.meta$Cell_type %>% table())
    if ((ct %in% c("foo"))) {
      print("skipped")
      next
    }
#     print("test")
    print(paste0("Getting stats for cell type: ", ct ))
    started.at=proc.time()
    results <- runStats(CELL_TYPE = ct, PARALLEL = T)
    print(timetaken(started.at))
    res[[ct]] <- results$results
    res.err[[ct]] <- results$error
    print(head(res[[ct]]))
  }
  save.image(paste0(RES_PATH, TISSUE_STR, "_workspace_tmp_step1_", DATE_STR, ".RData"))
}

# ================================================================ #



#### get mean score difference per rxn #############################
if (as.logical(args[4])) {
  print("===")
  print("Step 2: Computing fold changes...")

  dat.Compass.means <- tibble()
  # neutrophils
  for (t in (dat.meta$Cell_type %>% levels())) {
    started.at = proc.time()
    dat.tmp <- fread(file = paste0(DAT_PATH_RES, FNRES),
                    select = c("V1", (
                      dat.meta %>%
                        filter(Cell_type == t) %>%
                        pull(Cell_ID)
                    )))
    timetaken(started.at)
    dat.tmp.long <- dat.tmp %>%
      pivot_longer(-V1, names_to = "Cell_ID") %>%
      left_join(dat.meta, by = "Cell_ID")
    dat.Compass.means <- rbind(
      dat.Compass.means,
      dat.tmp.long %>% group_by(Cell_type, Phenotype, V1) %>%
        summarise(
          meanCompass = mean(value), 
          sdCompass = sd(value))
    )
  }
  rm(dat.tmp, dat.tmp.long)
  dat.Compass.means %<>% rename("rxn_id" = "V1")
  dat.Compass.means %<>% pivot_wider(names_from = c(Phenotype), values_from = c(meanCompass, sdCompass)) %>% 
    mutate(log2FC = log2((meanCompass_CLP+1)/(meanCompass_Sham+1)))
  dat.Compass.means %<>% left_join(dat.m.rxns.filt, by = "rxn_id")

  dat.Compass.means.summary <- dat.Compass.means %>% 
    group_by(Cell_type, subsystem) %>% 
    summarise(total = n(),
              log2FCgr1 = (abs(log2FC) > 1) %>% sum(),
              log2FCgr05 = (abs(log2FC) > 0.5) %>% sum(),
              log2FCgr025 = (abs(log2FC) > 0.25) %>% sum(),
              log2FCgr01 = (abs(log2FC) > 0.1) %>% sum()
              )
  save.image(paste0(RES_PATH, TISSUE_STR, "_workspace_tmp_step2_", DATE_STR, ".RData"))
}
# ================================================================ #

#### combine results ###############################################
if (as.logical(args[4])) {
  print("===")
  print("Step 3: Combine Results...")

  res.stats <- res %>% rbindlist(idcol = "Cell_type")
  res.stats %>% glimpse()
  colnames(res.stats)[3] <- "rxn_id"
  res.stats %<>% left_join(dat.m.rxns.filt, by = "rxn_id") 
  res.stats %>% glimpse()
  res.stats %<>% mutate(Cell_type = Cell_type %>% as_factor(),
                            rxn_id = rxn_id %>% as_factor(),
                            group1 = group1 %>% as_factor(),
                            group2 = group2 %>% as_factor(),
                            subsystem = subsystem %>% as_factor())
  res.stats %<>% rstatix::adjust_pvalue(p.col = "p", method = "BH")
  res.stats %<>% left_join(dat.Compass.means, join_by(Cell_type,rxn_id, subsystem))
  dat.meta %>% glimpse()
  res.stats$Cell_type %>% table()
  dat.Compass.means$Cell_type %>% table()
  
  res.stats.signif <- res.stats %>% filter(p.adj <= 0.05)

  save.image(paste0(RES_PATH, TISSUE_STR, "_workspace_tmp_step3_", DATE_STR, ".RData"))
}
res.stats.signif <- res.stats %>% filter(p.adj <= 0.05)
# ================================================================ #

#### output tables #################################################
res.stats %>% write_tsv(paste0(RES_PATH, TISSUE_STR, "_", "stats_", DATE_STR, ".tsv"))

save.image(paste0(RES_PATH, TISSUE_STR, "_workspace_", DATE_STR, ".RData"))
# ================================================================ #

