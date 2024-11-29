# Description: generate chord diagrams per tissue and pathways of interest
#
# author Sascha Schäuble
# date of creation: Fri Jul  5 16:26:40 2024
# license: MIT

library("tidyverse")
library("magrittr")
library("here")
library("janitor")
library("ggpubr")
library("ggfortify")
library("RColorBrewer")
library("data.table")
library("circlize")
library("gridExtra")
library("ComplexHeatmap") # for legends
library("ggsci")
library("scales")

PROJ_PATH <- here::here() %>% str_remove("src$")
DAT_PATH0 <- paste0(PROJ_PATH, "dat/modelDat/")
DAT_PATH1 <- paste0(PROJ_PATH, "dat/corr/")
RES_PATH  <- paste0(PROJ_PATH, "res/")

DATE_STR <- format(Sys.time(), "%y%m%d")

FN.ws <- "tissue_workspace.RData" # change according to tissue
FN.corr <- "corr.dat.zip"
# FN.corr.passphrase <- "" # paste in passphrase token here to enable unzip

FN.poi <- "poi.txt"
FN.rxns.meta <- "rxns.tsv"

# ================================================================ #


#### functions #####################################################


#'
#' @description: draw a chord diagram
#'
#' @param LGD - TRUE/FALSE draw a legend?
#' @param CATLABELS - if LGD = TRUE then required (named) label vector for categories 
#' @param CATCOLMAP - if LGD = TRUE then required (named) color mapping vector for categories, 
#' CATLABELS and CATCOLMAP should have identical names of categories present in the data set
#' @param PWY - which pathway - "all" for all data
#'
#' @return chord_fun
drawChord <- function(DF,
                      LGD = F,
                      CATLABELS = NULL,
                      CATCOLMAP = NULL,
                      PWY,
                      ORDER = NULL,
                      PRINT = T,
                      PATH = RES_PATH,
                      FN,
                      TITLE = "",
                      SCALE = F,
                      ADD_LABEL = F,
                      ADD_CT_TRACK = F) {
  if (LGD) {
    stopifnot(!is.null(CATLABELS))
    stopifnot(!is.null(CATCOLMAP))
  }
  if (is.null(ORDER)) {
    # standard order
    ORDER <- c(
      "Glycerophospholipid metabolism",
      "Triacylglycerol synthesis",
      "Glycine, serine, alanine and threonine metabolism",
      "Lysine metabolism"
    )
  }
  NUMPWCOL <- length(ORDER)
  fun.drawChord <- function() {
    lgd_cat <- Legend(
      at = names(CATLABELS),
      labels = CATLABELS,
      type = "points",
      legend_gp = gpar(col = CATCOLMAP),
      title_position = "topleft",
      title = "Correlation CLP / Sham",
      size = unit(5, "mm")
    )
    GRIDCOL <- ggsci::pal_d3()(NUMPWCOL)
    names(GRIDCOL) <- ORDER
    my_chord <- chordDiagram(
      x = (DF %>% dplyr::select("from", "to", "value")),
      order = ORDER,
      # group = DF %>% select(from, Cell_type) %>% deframe(),
      grid.col = GRIDCOL,
      col = DF$cat_col,
      self.link = 1,
      annotationTrack = c("grid"),
      annotationTrackHeight = mm_h(c(12, 16)),
      scale = SCALE,
      # preAllocateTracks = list(track.height = 0.1)
      
    )
    
    if (ADD_CT_TRACK) {
      if (TITLE == "brain") {
        idx <- 1:8
      } else if (TITLE == "kidney") {
        idx <- 8:17
      } else if (TITLE == "liver") {
        idx <- c(8, 11, 17:22)
      } else if (TITLE == "wat") {
        idx <- c(8, 15, 23:26)
      } else {
        error()
      }
      
      my_cols <- hue_pal()(26)
      # names(my_cols) <- DF$Cell_type %>% levels()
      names(my_cols) <- c(
        "Microglial cells",
        "Glutamatergic neurons",
        "GABAergic neurons",
        "Oligodendrocytes",
        "Dopaminergic neurons",
        #5
        "Oligodendrocyte precursor cells",
        "Astrocytes",
        "Endothelial",
        "Distal convoluted tubule",
        "Proximal tubule",
        #10
        "T lymph / NK",
        "Ascending loop of Henle",
        "Collecting duct principal cell",
        "Podocytes",
        "Macrophages",
        #15
        "Collecting duct intercalated cell",
        "B lymph",
        "Hepatocytes",
        "Kupffer cells",
        "Hepatic stellate cells",
        #20
        "Cholangiocytes",
        "Neutrophils",
        "Adipocyte",
        "Pericytes",
        "Adipose stem and progenitor cells",
        #25
        "Epithelial"
      )
      my_cols <- my_cols[idx]
      
      ylim = get.cell.meta.data("ylim", sector.index = ORDER[1], track.index = 1)
      # y1 = ylim[1] + (ylim[2] - ylim[1]) * 0.2
      y1 = ylim[1] + (ylim[2] - ylim[1]) * 0.4
      y2 = ylim[2]
      for (i in seq_len(nrow(my_chord))) {
        if (my_chord$value1[i] > 0) {
          if (my_chord$col[i] == "#00000000") {
            # dummy_col <- "#00000000"
            dummy_col <- "white"
          } else {
            dummy_col <- my_cols[DF$Cell_type[i]]
          }
          circos.rect(
            xleft = my_chord[i, "x1"],
            ybottom = y1,
            xright = my_chord[i, "x1"] - abs(my_chord[i, "value1"]),
            ytop = y1 + (y2 - y1) * 0.9,
            col = dummy_col,
            border = "white",
            sector.index = my_chord$rn[i],
            track.index = 1
          )
          if (my_chord[i, "x1"] != my_chord[i, "x2"]) {
            circos.rect(
              my_chord[i, "x2"],
              y1,
              my_chord[i, "x2"] - abs(my_chord[i, "value1"]),
              y1 + (y2 - y1) * 0.9,
              col = dummy_col,
              border = "white",
              sector.index = my_chord$cn[i],
              track.index = 1
            )
          }
          
        }
      }
    }
    
    if (LGD) {
    draw(
      lgd_cat,
      x = unit(4, "mm"),
      y = unit(4, "mm"),
      just = c("left", "bottom")
    )
    }
    if (ADD_LABEL) {
      sectors <- unique(c(DF$from, DF$to))
      for (sector in sectors) {
        circos.text(
          x = get.cell.meta.data(sector.index = sector, name = "xcenter"),
          y = get.cell.meta.data(sector.index = sector, name = "ycenter"),
          labels = sector,
          sector.index = sector,
          facing = "downward",
          niceFacing = TRUE,
          adj = c(0, 0.5),
          cex = 1.3,
        )
      }
    }
    
    circos.clear()
  }

  if (PRINT) {
    plotCircos(
      PLOTFUN = fun.drawChord,
      PATH = PATH,
      FN = FN,
      TITLE = TITLE
    )
  }
  return(fun.drawChord())
}
# ================================================================ #

#'
#' @description: convenience function to save a circos plot
#'
#' @param variables
#'
#' @return none
plotCircos <- function(PLOTFUN, PATH = RES_PATH, FN, TITLE = "") {
  cairo_pdf(paste0(PATH, FN), width = 14, height = 10, family = "Arial")
  PLOTFUN()
  title(TITLE, cex = 1.5)
  dev.off()
}
# ================================================================ #

#'
#' @description: summarize information for correlations and the category of correlations
#'
#' @param variables
#'
#' @return variables
get_CorrelationCategories <- function(DF) {
  
  corr.dat.cat <- list()
  corr.lgd.labels <- list()
  for (ct in DF$Cell_type %>% as_factor() %>% levels()) {
    corr.dat.cat[[ct]] <- DF %>% filter(Cell_type == ct)
    corr.dat.cat[[ct]] %<>% mutate(category = if_else((CLP >=  0.3 & CLP <= 1   & Sham >= -0.3 & Sham <  0.3), "cat1", NA))       # +/0 "#0a9396"
    corr.dat.cat[[ct]] %<>% mutate(category = if_else((CLP >=  0.3 & CLP <= 1   & Sham >= -1   & Sham < -0.3), "cat2", category)) # +/-
    corr.dat.cat[[ct]] %<>% mutate(category = if_else((CLP >=  0.5 & CLP <= 1   & Sham >=  0.3 & Sham <  0.5), "cat3", category)) # ++/+
    corr.dat.cat[[ct]] %<>% mutate(category = if_else((CLP >= -0.3 & CLP <  0.3 & Sham >= -1   & Sham < -0.3), "cat4", category)) # 0/-
    corr.dat.cat[[ct]] %<>% mutate(category = if_else((CLP >=  0.3 & CLP <  0.5 & Sham >=  0.5 & Sham <= 1  ), "cat5", category)) # +/++
    corr.dat.cat[[ct]] %<>% mutate(category = if_else((CLP >= -0.3 & CLP <  0.3 & Sham >=  0.3 & Sham <= 1  ), "cat6", category)) # 0/+
    corr.dat.cat[[ct]] %<>% mutate(category = if_else((CLP >= -1   & CLP < -0.3 & Sham >=  0.3 & Sham <= 1  ), "cat7", category)) # -/+
    corr.dat.cat[[ct]] %<>% mutate(category = if_else((CLP >= -1   & CLP < -0.3 & Sham >= -0.3 & Sham <  0.3), "cat8", category)) # -/0
    # corr.dat.cat[[ct]] %<>% mutate(category = if_else((CLP >= -1   & CLP < -0.5 & Sham >= -0.5 & Sham <= -0.3  ), "cat9", category)) # --/-
    # corr.dat.cat[[ct]] %<>% mutate(category = if_else((CLP >= -0.5   & CLP < -0.3 & Sham >= -1 & Sham <  -0.5), "cat10", category)) # -/--
    corr.dat.cat[[ct]] %<>% mutate(category = if_else((
      (CLP >= -1 & CLP <  -0.5 & Sham >= -1 & Sham <  -0.5) | # --/--
        (CLP >= -0.5 & CLP <  -0.3 & Sham >= -0.5 & Sham <  -0.3) | # -/-
        (CLP >= -0.3 & CLP <   0.3 & Sham >= -0.3 & Sham <   0.3) | # 0/0
        (CLP >=  0.3 & CLP <   0.5 & Sham >=  0.3 & Sham <   0.5) | # +/+
        (CLP >=  0.5 & CLP <=  1   & Sham >=  0.5 & Sham <=  1) # ++/++
    ), "cat9", category)) # no change
    # ), "cat11", category)) # no change

  }

  # get summary of correlation types
  corr.dat.cat.df <- corr.dat.cat %>% data.table::rbindlist(idcol = "CT") %>% select(-Cell_type) %>% 
    rename("Cell_type" = "CT")

  corr.dat.cat.df.summary <- corr.dat.cat.df %>% group_by(Cell_type, subsystemi, subsystemj, category) %>% 
    summarise(count = n())
  # remove redundant A -> B, B -> A information as direction does not play a role
  corr.dat.cat.df.summary %<>%
    mutate(
      dummy_i = pmin(subsystemi, subsystemj),  # Smaller group name first
      dummy_j = pmax(subsystemi, subsystemj)   # Larger group name second
    ) %>% group_by(Cell_type, category) %>% 
    distinct(dummy_i, dummy_j, .keep_all = TRUE) %>%
    select(-subsystemi, -subsystemj) %>%
    rename(subsystemi = dummy_i, subsystemj = dummy_j) %>% ungroup()
  
  # add colors to categories
  cat_col_map <- c("#0a9396", "#005f73", "#94d2bd", "#cfe1b9", "#ee9b00", "#ca6702", "#ffbf69", "#dee1b9",
                   # "grey80"
                   # "#d29495", "#ee80e7",
                   "#00000000"
                   )
  names(cat_col_map) <- paste0("cat", 1:9)
  # names(cat_col_map) <- paste0("cat", 1:11)
  cat_labels <- c("+/0", "+/-", "++/+", "0/-", "+/++", "0/+", "-/+", "-/0", "no change") 
  # cat_labels <- c("+/0", "+/-", "++/+", "0/-", "+/++", "0/+", "-/+", "-/0", "--/-", "-/--", "no change") 
  names(cat_labels) <- names(cat_col_map)
  corr.dat.cat.df.summary %<>% 
    left_join(tibble(category = names(cat_col_map), cat_col = cat_col_map), by = "category")# %>% 
  corr.dat.cat.df.summary$Cell_type %<>% as_factor()
  corr.dat.cat.df.summary$category %<>% as_factor()
  corr.dat.cat.df.summary$subsystemi %<>% as_factor()
  corr.dat.cat.df.summary$subsystemj %<>% as_factor()
  
  corr.dat.cat.df.summary$Cell_type %>% levels()
  for (ct in corr.dat.cat.df.summary$Cell_type %>% levels()) {
    dummyNum <- corr.dat.cat.df.summary %>% filter(Cell_type == ct) %>% group_by(category) %>% summarise(count = sum(count)) %>% 
      deframe()
    dummyLabel <- c("+/0", "+/-", "++/+", "0/-", "+/++", "0/+", "-/+", "-/0", "no change")
    # dummyLabel <- c("+/0", "+/-", "++/+", "0/-", "+/++", "0/+", "-/+", "-/0", "--/-", "-/--", "no change")
    names(dummyLabel) <- paste0("cat", 1:9)
    # names(dummyLabel) <- paste0("cat", 1:11)
    idx <- match(names(dummyNum), names(dummyLabel))
    dummyLabel[idx] <- paste0(dummyLabel[idx],  if_else(dummyLabel[idx] == "no change","\t", "\t\t\t"), "n=",dummyNum)
      corr.lgd.labels[[ct]] <- c("+/0", "+/-", "++/+", "0/-", "+/++", "0/+", "-/+", "-/0", "no change")
    # corr.lgd.labels[[ct]] <- c("+/0", "+/-", "++/+", "0/-", "+/++", "0/+", "-/+", "-/0", "--/-", "-/--", "no change") 
    corr.lgd.labels[[ct]] <- dummyLabel
  }
  
  return(list(corr.links = corr.dat.cat.df.summary, corr.lgd = corr.lgd.labels))
}
# ================================================================ #


#### Data wrangling ################################################
poi <- read_tsv(paste0(DAT_PATH0, FN.poi), col_names = F) %>% pull()
rxns.meta <- read_tsv(paste0(DAT_PATH0, FN.rxns.meta), col_names = F) %>% rename("rxn_ID" = "X1", "subsystem" = "X2")

# unzip corr data
stopifnot(FN.corr.passphrase != "")
system(paste0("unzip -P ", FN.corr.passphrase, " ", DAT_PATH1, FN.corr, " -d ", DAT_PATH1 ))
FN.corr <- list.files(path = paste0(DAT_PATH1), pattern = "*.tsv")

corr.dat <- list()
for (i in FN.corr) {
  corr.dat[[i %>% str_extract("^[^_]+")]] <- read_tsv(paste0(DAT_PATH1, i))
  corr.dat[[i %>% str_extract("^[^_]+")]] %<>% left_join(rxns.meta, by = c("rxn_IDi" = "rxn_ID")) 
  corr.dat[[i %>% str_extract("^[^_]+")]] %<>% left_join(rxns.meta, by = c("rxn_IDj" = "rxn_ID"), suffix = c("i", "j")) 
  corr.dat[[i %>% str_extract("^[^_]+")]] %<>% filter(subsystemi %in% poi, subsystemj %in% poi) %>% 
    pivot_wider(names_from = condition, values_from = value) %>% drop_na()
}

corr.dat.summaries <- list()
for (i in names(corr.dat) ) {
  corr.dat.summaries[[i]] <- get_CorrelationCategories(DF = corr.dat[[i]])
}
cat_col_map <- c("#0a9396", "#005f73", "#94d2bd", "#cfe1b9", "#ee9b00", "#ca6702", "#ffbf69", "#dee1b9", 
                 # "#d29495", "#ee80e7",
                 "#00000000"
)
names(cat_col_map) <- paste0("cat", 1:9)
# names(cat_col_map) <- paste0("cat", 1:11)
cat_labels <- c("+/0", "+/-", "++/+", "0/-", "+/++", "0/+", "-/+", "-/0", "no change")
# cat_labels <- c("+/0", "+/-", "++/+", "0/-", "+/++", "0/+", "-/+", "-/0", "--/-", "-/--", "no change") 
names(cat_labels) <- names(cat_col_map)

# ================================================================ #

#### circos plots with correlations ################################
pb <- txtProgressBar(style = 3, max = length(names(corr.dat.summaries)))
it <- 0

# per tissue and cell type
my_chords <- list()
for (tis in names(corr.dat.summaries)) {
  for (ct in ( corr.dat.summaries[[tis]]$corr.links$Cell_type %>% levels() ) ) {
    tmp <- corr.dat.summaries[[tis]]$corr.links %>% filter(Cell_type == ct) %>%
      rename("from" = "subsystemi",
             "to" = "subsystemj",
             "value" = "count") %>%
      ungroup()

    my_chords[[paste0(tis, "_", ct)]] <- drawChord(
      DF = tmp,
      CATLABELS = corr.dat.summaries[[tis]]$corr.lgd[[ct]],
      CATCOLMAP = cat_col_map,
      PRINT = T,
      FN = paste0(tis, "_", janitor::make_clean_names(ct), "_", DATE_STR, ".pdf"),
      TITLE = ct,
      SCALE = F, 
      ADD_LABEL = T,
      LGD = T
    )
    my_chords[[paste0(tis, "_", ct)]]
  }
  it <- it +1
  setTxtProgressBar(pb, it)
}
close(pb)

# ignoring Cell type
my_chords.2 <- list()
for (tis in names(corr.dat.summaries)) {
    tmp <- corr.dat.summaries[[tis]]$corr.links %>% 
      rename("from" = "subsystemi",
             "to" = "subsystemj",
             "value" = "count") %>%
      group_by(from, to, category, cat_col) %>% 
      summarize(value = sum(value)) %>% 
      ungroup()

    ## labels
    dummyNum <- tmp %>% group_by(category) %>% summarise(count = sum(value)) %>% 
      deframe()
    dummyLabel <- rep("", 9)
    # dummyLabel <- rep("", 11)
    names(dummyLabel) <- paste0("cat", 1:9)
    # names(dummyLabel) <- paste0("cat", 1:11)
    idx <- match(names(dummyNum), names(dummyLabel))
    dummyLabel[idx] <- paste0(dummyLabel[idx],  dummyNum)
    
    my_chords.2[[tis]] <- drawChord(
      DF = tmp,
      CATLABELS = dummyLabel,
      CATCOLMAP = cat_col_map,
      PRINT = T,
      FN = paste0(tis, "_allCT", "_", DATE_STR, ".pdf"),
      TITLE = tis,
      SCALE = F, 
      ADD_LABEL = T,
      LGD = T
    )
    my_chords.2[[tis]]
}

# showing Cell type as additional layer
my_chords.3 <- list()
for (tis in names(corr.dat.summaries)) {
  tmp <- corr.dat.summaries[[tis]]$corr.links %>% 
    rename("from" = "subsystemi",
           "to" = "subsystemj",
           "value" = "count") %>%
    group_by(Cell_type, from, to, category, cat_col) %>% 
    summarize(value = sum(value)) %>% 
    ungroup()

  ## labels
  dummyNum <- tmp %>% group_by(category) %>% summarise(count = sum(value)) %>% 
    deframe()
  dummyLabel <- rep("", 9)
  names(dummyLabel) <- paste0("cat", 1:9)
  idx <- match(names(dummyNum), names(dummyLabel))
  dummyLabel[idx] <- paste0(dummyLabel[idx],  dummyNum)
  
  my_chords.3[[tis]] <- drawChord(
    DF = tmp, 
    CATLABELS = dummyLabel, 
    CATCOLMAP = cat_col_map,
    PRINT = T,
    FN = paste0(tis, "_allCTring", "_", DATE_STR, ".pdf"),
    TITLE = tis,
    SCALE = F, 
    ADD_LABEL = T,
    LGD = T, 
    ADD_CT_TRACK = T
  )
  my_chords.3[[tis]]
}

# ================================================================ #
