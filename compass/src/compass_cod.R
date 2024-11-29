# Description: visualization of Cohen's D results
#
# author Sascha Schäuble
# date of creation: Tue May 28 14:46:21 2024
# license: MIT

library("here")
library("tidyverse")
library("magrittr")
library("ggpubr")
library("ggridges")

# constants 
PROJ_PATH <- here() %>% str_remove("src$")
DAT_PATH0 <- paste0(PROJ_PATH, "dat/modelDat/")
DAT_PATH <- paste0(PROJ_PATH, "dat/")
RES_PATH <- paste0(PROJ_PATH, "res/")

DATE_STR <- format(Sys.time(), "%y%m%d")

FN.meta  <- "rxns_meta.tsv"
FN.poi <- "poi.txt"
FN.stats <- "stats.Rdata.zip" # needs to be decompressed first, see passphrase and code below
FN.passphrase <- "" # paste in passphrase token here to enable unzip | available upon request

SAVE_OUT <- F
# ================================================================ #

#### functions #####################################################

#'
#' @description: get ridge plots over cell types
#'
#' @param LOG2FC boolean | T: log2fC, F: cohensd 
#'
#' @return variables
get_ridgesCt <- function(DF, TISSUE, SUBS = NULL, CAPTION = NULL, LOG2FC = T) {
  df <- DF %>% filter(tissue == TISSUE)
  if (!is.null(SUBS)) {
    df %<>% filter(subsystem %in% SUBS)
  }
  
  if (LOG2FC) {

  highlight_ranges <- data.frame(
    xmin = c(min(df$log2FC), 0.1),
    xmax = c(-0.1, max(df$log2FC)),
    ymin = -Inf,
    ymax = Inf
  )
  
  p <- df %>% 
    ggplot(aes(x = log2FC, y = subsystem %>% fct_relevel(rev(sort(levels(.)))), 
               group = interaction(subsystem, Cell_type), fill = Cell_type)) +
    geom_rect(data = highlight_ranges, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    fill = "gray30", alpha = 0.3, inherit.aes = F) +
    geom_density_ridges(alpha = 0.7, panel_scaling = F, quantile_lines=TRUE) +
    theme_ridges() +
    theme(legend.position = "right") + 
    geom_vline(xintercept=0.1, linetype="dashed", color  ="gray10") + 
    geom_vline(xintercept=-0.1, linetype="dashed", color  ="gray10") + 
    theme_ridges(grid=T, line_size = 0.2) +
    labs(title = "COMPASS distribution for CLP vs Sham",
         subtitle = TISSUE,
         x = expression(log[2](FC(COMPASS~score))), y = "",
         caption = CAPTION
    )
  
  } else { # cohen's D
    
    highlight_neg <- data.frame(
      xmin = c(-0.2),
      xmax = c(0.2),
      ymin = -Inf,
      ymax = Inf
    )
    highlight_small <- data.frame(
      xmin = c(-0.2, 0.2),
      xmax = c(-0.5, 0.5),
      ymin = -Inf,
      ymax = Inf
    )
    highlight_moderate <- data.frame(
      xmin = c(-0.5, 0.5),
      xmax = c(-0.8, 0.8),
      ymin = -Inf,
      ymax = Inf
    )
    highlight_large <- data.frame(
      xmin = c(-0.8, 0.8),
      xmax = c(-Inf, Inf),
      ymin = -Inf,
      ymax = Inf
    )
    
    
    p <- df %>% 
      ggplot(aes(x = effsize, y = subsystem %>% fct_relevel(rev(sort(levels(.)))), 
                 group = interaction(subsystem, Cell_type), fill = Cell_type)) +
      geom_rect(data = highlight_small, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
                fill = "gray90", alpha = 1, inherit.aes = F)+
      geom_rect(data = highlight_small, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
                fill = "gray80", alpha = 1, inherit.aes = F) +
      geom_rect(data = highlight_moderate, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
                fill = "gray70", alpha = 1, inherit.aes = F) +
      geom_rect(data = highlight_large, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
                fill = "gray60", alpha = 1, inherit.aes = F) +
      geom_density_ridges(alpha = 0.7, panel_scaling = F, quantile_lines=TRUE, quantile_fun=function(effsize,...)mean(effsize)) +
      theme_ridges() +
      theme(legend.position = "right") + 
      geom_vline(xintercept=0.2, linetype="dashed", color  ="gray10") + 
      geom_vline(xintercept=-0.2, linetype="dashed", color  ="gray10") + 
      theme_ridges(grid=T, line_size = 0.2) +
      labs(title = "COMPASS distribution for CLP vs Sham",
           subtitle = TISSUE,
           x = "Cohen's d", y = "",
           caption = CAPTION
      )
  }
  
  return(p)
}
# ================================================================ #

# ================================================================ #




#### Data wrangling ################################################

subs.poi <- read_tsv(paste0(DAT_PATH0, FN.poi), col_names = F) %>% pull()

meta <- read_tsv(paste0(DAT_PATH0, FN.meta))

# unzip stats data
stopifnot(FN.passphrase != "")
system(paste0("unzip -P ", FN.passphrase, " ", DAT_PATH, FN.stats, " -d ", DAT_PATH ))
# load stats data
load(paste0(DAT_PATH, FN.stats %>% str_remove(".zip")))

# ================================================================ #

#### distribution | ridges #########################################

## tissues on one line per pathway 
highlight_ranges <- data.frame(
  xmin = c(-0.98, 0.1),
  xmax = c(-0.1, 0.55),
  ymin = -Inf,
  ymax = Inf

)
highlight_neg <- data.frame(
  xmin = c(-0.2),
  xmax = c(0.2),
  ymin = -Inf,
  ymax = Inf
)
highlight_small <- data.frame(
  xmin = c(-0.2, 0.2),
  xmax = c(-0.5, 0.5),
  ymin = -Inf,
  ymax = Inf
)
highlight_moderate <- data.frame(
  xmin = c(-0.5, 0.5),
  xmax = c(-0.8, 0.8),
  ymin = -Inf,
  ymax = Inf
)
highlight_large <- data.frame(
  xmin = c(-0.8, 0.8),
  xmax = c(-Inf, Inf),
  ymin = -Inf,
  ymax = Inf
)

p.subs.signif.ridge <- stats.df.signif %>% filter(subsystem != "Glycosphingolipid metabolism") %>% 
  ggplot(aes(x = effsize, y = subsystem %>% fct_relevel(rev(sort(levels(.)))), 
                            group = interaction(subsystem, tissue), fill = tissue)) +
  geom_rect(data = highlight_neg, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = "gray90", alpha = 1, inherit.aes = F) +
  geom_rect(data = highlight_small, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = "gray80", alpha = 1, inherit.aes = F) +
  geom_rect(data = highlight_moderate, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = "gray70", alpha = 1, inherit.aes = F) +
  geom_rect(data = highlight_large, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = "gray60", alpha = 1, inherit.aes = F) +
  geom_density_ridges(alpha = 0.7, panel_scaling = F) +

  ggsci::scale_fill_nejm(name = "Tissue") +

  theme_ridges() +
  theme(legend.position = "right") + 
  theme_ridges(grid=T, line_size = 0.2) +
  labs(

       x = "Cohen's d", y = "",

       )


cairo_pdf(file = paste0(RES_PATH,"ridges_signif_", DATE_STR, ".pdf"), family = "Arial", width = 10, height = 6)
p.subs.signif.ridge
dev.off()

# ================================================================ #


#### top rxns per pwy and tissue ###################################
TOPRANK <- 5
res.top5rxns <- stats.df %>% select(tissue, subsystem, Cell_type, rxn_id, effsize) %>%
  mutate(rxn_id_base = rxn_id %>% str_remove("_pos|_neg")) %>% 
  left_join(meta %>% select(-subsystem), by = c("rxn_id_base" = "rxn_id")) %>% 
  filter(!is.na(EC)) %>% 

  group_by(subsystem) %>%
  arrange(effsize) %>% 
  mutate(rank = row_number()) %>% 
  filter(rank <= TOPRANK | rank > (n() - TOPRANK)) %>%
  ungroup() %>%
  arrange(tissue, Cell_type, rank)
res.top5rxns$subsystem %<>% fct_drop()

res.top5rxns %>% 
  write_tsv(paste0(RES_PATH, "top5_rxnsPerTissueCT_", DATE_STR, ".tsv"))
stats.df %>% 
  mutate(rxn_id_base = rxn_id %>% str_remove("_pos|_neg")) %>%
  left_join(meta %>% select(rxn_id, EC), by = c("rxn_id_base" = "rxn_id")) %>% 
  filter(rxn_id %in% res.top5rxns$rxn_id) %>% 
  select(tissue, Cell_type, subsystem, rxn_id, effsize, magnitude, EC) %>% 
  write_tsv(paste0(RES_PATH, "top5_rxnsOverAllPwys_", DATE_STR, ".tsv"))
# ================================================================ #


#### volcano plot ##################################################

# add labels to be plotted
stats.df.volc <- stats.df
stats.df.volc %<>% left_join(res.top5rxns %>% 
                               select(tissue, Cell_type, subsystem, rxn_id, EC),
                             by = join_by(tissue, Cell_type, subsystem, rxn_id)) %>%
  rename("Top5label" = "EC")

xLIM <- c(stats.df.volc$effsize %>% min, stats.df.volc$effsize %>% max)

p.compass.volc <- list()
for (s in subs.poi) {
    p.compass.volc[[s]] <-
      stats.df.volc %>%
      filter(subsystem == s) %>%
      mutate(negLog10Wilc = -log10(p.adj)) %>%
      ggscatter(
        x = "effsize",
        y = "negLog10Wilc",
        color = "Cell_type",
        shape = "tissue",
        xlab = "Cohen's d",
        size = 3,
        title = s, 
        label = "Top5label",
        label.rectangle = T,
        repel = T
      ) %>% ggpar(legend = "right", xlim = xLIM, ylim = c(0, 300)) +
      geom_vline(
        xintercept = 0.2,
        linetype = "dashed",
        color  = "gray10"
      ) +
      geom_vline(
        xintercept = -0.2,
        linetype = "dashed",
        color  = "gray10"
      ) +
      geom_hline(
        yintercept = -log10(0.05),
        linetype = "dashed",
        color  = "gray10"
      ) +
      scale_shape_manual(values = c(5:8))
}

cairo_pdf(family = "Arial", filename = paste0(RES_PATH, "volc_poi_", DATE_STR, ".pdf"), width = 10, height = 5,
          onefile = T)
print(p.compass.volc$`Glycine, serine, alanine and threonine metabolism`)
print(p.compass.volc$`Lysine metabolism`)
print(p.compass.volc$`Glycerophospholipid metabolism`)
print(p.compass.volc$`Triacylglycerol synthesis`)
dev.off()
# ================================================================ #
