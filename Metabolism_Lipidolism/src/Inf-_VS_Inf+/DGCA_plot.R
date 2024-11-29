# Title: Differential correlation network plot - human
# Project: Mitochondrial Adaption at the Center of Early Metabolic Changes in Presymptomatic Sepsis Patients
# Maintainer: Lin-Lin Xu
# Last update: Nov 27, 2024
# License: MIT

## load library ----------------------------------------------------------------
library(tidyverse)
## Load data -------------------------------------------------------------------
node_file <- read_tsv("DGCA_human_nodes.tsv")
## Results ---------------------------------------------------------------------
re_list <- readRDS("Human_ddcor_megena_res_residual.rds")
names <- re_list$summary$modules$c1_19

Node_new <- node_file %>%
  mutate(hub = if_else(sign %in% "sign",sign,hub)) %>%
  mutate(label = if_else(hub %in% c("hub"),Name,"")) %>%
  mutate(label = str_replace_all(label, c("lysoPC" = "LPC", "lysoPE" = "LPE"))) %>%
  mutate(group1 = str_replace_all(group1, c("lysoPC" = "LPC", "lysoPE" = "LPE"))) %>%
  mutate(module = if_else(Name %in% names,"Yes","No")) %>%
  mutate(label = str_remove_all(label,"_X.*|_Y.*"))

clean_empty_strings <- function(x) {
  x[x != ""]
}

module_hubs <- re_list$summary$module.table %>%
  dplyr::select(module.id,module.hub) %>%
  rowwise() %>%
  mutate(module.hub = str_split(module.hub,"\\(\\d+\\),|\\(\\d+\\)$")) 
module_hubs$module.hub <- lapply(module_hubs$module.hub, clean_empty_strings)

Tab_S2 <- re_list$modules
# Initialize the Hub column with FALSE
Tab_S2$Hub <- FALSE

# Check the conditions and update the Hub column
for (i in 1:nrow(Tab_S2)) {
  module_id <- Tab_S2$Modules[i]
  gene <- Tab_S2$Genes[i]
  
  # Find the corresponding B list in df_previous
  b_list <- module_hubs$module.hub[module_hubs$module.id == module_id]
  
  # If B list is not empty, check if the gene is in the list
  if (length(b_list) > 0) {
    if (gene %in% unlist(b_list)) {
      Tab_S2$Hub[i] <- TRUE
    }
  }
}

library(xlsx)
write.xlsx(Tab_S2, "Supplementary_Tab2.xlsx",
           sheetName="Table S2",row.names = T,append=TRUE)



