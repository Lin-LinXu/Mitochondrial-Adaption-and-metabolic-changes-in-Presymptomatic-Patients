# Title: Ordinal Regression Plot
# Project: Mitochondrial Adaption at the Center of Early Metabolic Changes in Presymptomatic Sepsis Patients
# Maintainer: Lin-Lin Xu
# Last update: Nov 24, 2024
# License: MIT

## Load library ----------------------------------------------------------------
library(tidyverse)
library(janitor)
library('RColorBrewer')
library(patchwork)
library(ggprism)
set.seed(123)
# Read data --------------------------------------------------------------------
# Tidy data --------------------------------------------------------------------
Excel_data_list <- readRDS("Compound_data.rds")

Data_lst <- readRDS("human_data.rds")
Met_file <- Data_lst$Met_file
Lip_file <- Data_lst$Lip_file
Met_ID_info <- Excel_data_list$Metabolite_data_df$compound_info %>%
  dplyr::select(Name,CompoundID) %>%
  dplyr::rename(metabolite = CompoundID) %>%
  mutate(Name = str_replace_all(Name, c("Kynuerine" = "Kynurenine")))
Lip_ID_info <- Excel_data_list$Lipid_data_df$compound_info %>%
  dplyr::select(Name,CompoundID) %>%
  dplyr::rename(metabolite = CompoundID)

Data_time_lst <- readRDS("human_data_time.rds")
Met_file_time <- Data_time_lst$Met_file
Lip_file_time <- Data_time_lst$Lip_file

jco_pal = c("#3C5488FF","#F39B7FFF","#DC0000FF")
## main part -------------------------------------------------------------------
Ordinal_reg_all_result_list_BH <- readRDS("Ordinal_reg_BH.rds")
sign_metabolite <- Ordinal_reg_all_result_list_BH$Met_Ord$re_df_025$name
sign_lipid <- Ordinal_reg_all_result_list_BH$Lip_Ord$re_df_025$name
Lipid_class_info <- read_tsv("Lipid_class_info.tsv")
inf_vac <- c("SIRS-",
             "Uncomplicated Infection",
             "Sepsis")
# Helper function for string wrapping. 
# Default 20 character target width.
swr = function(string, nwrap=5) {
  paste(strwrap(string, width=nwrap), collapse="\n")
}
swr = Vectorize(swr)
## Metabolite ------------------------------------------------------------------
Plotdf_met <- Met_file_time %>%
  mutate(Name = str_replace_all(Name, c("Kynuerine" = "Kynurenine"))) %>%
  filter(CompoundID %in% sign_metabolite) %>%
  filter(Class2 %in% inf_vac) %>%
  mutate(Class2 = factor(Class2, levels = c("SIRS-",
                                            "Uncomplicated Infection",
                                            "Sepsis"))) %>%
  filter(Timepoint %in% c("-3","-1")) %>%
  mutate(Timepoint = factor(Timepoint, levels = c("-3","-1")))
Plotdf_met$Name <- swr(Plotdf_met$Name)

# Function to calculate the lower and upper whiskers
calculate_whiskers <- function(x) {
  q1 <- quantile(x, 0.25, na.rm = T)
  q3 <- quantile(x, 0.75, na.rm = T)
  iqr <- q3 - q1
  lower_whisker <- max(min(x, na.rm = T), q1 - 1.5 * iqr)
  upper_whisker <- min(max(x, na.rm = T), q3 + 1.5 * iqr)
  return(c(lower_whisker, upper_whisker))
}

# Calculate whiskers for each group
whiskers <- Plotdf_met %>%
  group_by(Name) %>%
  summarize(
    lower_whisker = calculate_whiskers(abundance)[1],
    upper_whisker = calculate_whiskers(abundance)[2]
  )

# Merge whiskers back to the original dataframe
Plotdf_met_new <- Plotdf_met %>%
  left_join(whiskers, by = "Name") %>%
  filter(abundance <= upper_whisker & abundance >= lower_whisker)


Plot1 <- ggplot2::ggplot(Plotdf_met,aes(x = Timepoint, y = abundance,fill = Class2)) +
  geom_boxplot(outliers = FALSE) +
  ggplot2::geom_point(aes(x = Timepoint, y = abundance, fill = Class2), data = Plotdf_met_new, size = 0.6, shape = 21, alpha = 0.50,position = position_jitterdodge()) +
  facet_wrap(vars(Name),scales = "free", nrow = 3) +
  labs(x = "", y="normalised intensity", 
       fill ="") +
  scale_fill_manual(values = jco_pal) +
  geom_hline(yintercept = 0, linetype="dashed", color = "red") +
  theme_prism(base_size = 16) +
  theme(legend.position="none") +
  theme(
    plot.title = element_text(size=16, face= "bold", colour= "black" ),
    axis.title.x = element_text(size=16,face= "bold",  colour = "black"),    
    axis.title.y = element_text(size=16, face="bold", colour = "black"),    
    axis.text.x = element_text(size=16, face="bold", colour = "black"), 
    axis.text.y = element_text(size=16, face="bold", colour = "black"), # bold
    strip.text.x = element_text(size = 16, face="bold", colour = "black" ),
    strip.text.y = element_text(size = 16, face="bold", colour = "black"),
    axis.line.x = element_line(color="black", size = 1),
    axis.line.y = element_line(color="black", size = 1),
    legend.background=element_blank(),
    legend.text = element_text(size=16, face= "bold"),
    legend.title = element_text(size=16, face= "bold")
  )


## Lipid -----------------------------------------------------------------------
Plotdf_lip <- Lip_file_time %>%
  filter(CompoundID %in% sign_lipid) %>%
  filter(Class2 %in% inf_vac) %>%
  filter(!Timepoint %in% "-2") %>%
  mutate(Timepoint = factor(Timepoint, levels = c("Pre","-3","-1"))) %>%
  left_join(Lipid_class_info,by = c("CompoundID" = "Name")) %>%
  mutate(Prefix = str_replace_all(Prefix, c("lysoPC" = "LPC",
                                            "lysoPE" = "LPE"))) %>%
  mutate(Name = str_replace_all(Name, c("lysoPC" = "LPC",
                                        "lysoPE" = "LPE")))
  
##* plot 1 ---------------------------------------------------------------------
Plotdf_lip$Name <- swr(Plotdf_lip$Name)

lipid_direction <- Ordinal_reg_all_result_list_BH$Lip_Ord$clm_df %>%
  filter(metabolite %in% sign_lipid) %>%
  dplyr::select(metabolite,estimate)

Lipid_groups <- c("TG","LPC","LPE","PC","CerP","SM",
                  "plasmenyl-PE","plasmenyl-PC","DG","PI","DGDG",
                  "PA","PS","MG")

Plotdf_lip_tmp <- Plotdf_lip %>%
  left_join(lipid_direction,by = c("CompoundID" = "metabolite")) %>%
  mutate(color = if_else(estimate < 0,"blue","red")) %>%
  mutate(Prefix = factor(Prefix,levels = Lipid_groups))

Plotdf_lip_tmp$Name <- swr(Plotdf_lip_tmp$Name)
  
p <- Plotdf_lip_tmp %>%
  mutate(Class2 = factor(Class2, levels = c("SIRS-",
                                            "Uncomplicated Infection",
                                            "Sepsis"))) %>%
  dplyr::rename(group = Prefix) %>%
  group_by(CompoundID, Class2, group,Timepoint,color) %>%
  summarise(mean=mean(abundance,na.rm = TRUE), sd=sd(abundance,na.rm = TRUE)) %>% 
  filter(!is.na(mean)) %>%
  filter(Class2 %in% "Sepsis") %>%
  ggplot2::ggplot(aes(x = Timepoint, y = mean, colour = color)) +
  geom_point() +
  geom_line(aes(x = Timepoint, y = mean, group = CompoundID, colour = color)) +
  labs(x="", y="normalised intensity",colour = "") +
  geom_hline(yintercept = 0, linetype="dashed", color = "#DC0000FF") +
  facet_wrap(vars(group),scales = "free", ncol = 7) +
  scale_color_manual(values = c("#3C5488FF","red")) +
  theme_prism(base_size = 16) +
  theme(legend.position="none") +
  theme(
    # LABELS APPEARANCE
    legend.position="none",
    plot.title = element_text(size=16, face= "bold", colour= "black" ),
    axis.title.x = element_text(size=16,face= "bold",  colour = "black"),    
    axis.title.y = element_text(size=16, face="bold", colour = "black"),    
    axis.text.x = element_text(size=16, face="bold", colour = "black"), 
    axis.text.y = element_text(size=16, face="bold", colour = "black"), # bold
    strip.text.x = element_text(size = 16, face="bold", colour = "black" ),
    strip.text.y = element_text(size = 16, face="bold", colour = "black"),
    axis.line.x = element_line(color="black", size = 1),
    axis.line.y = element_line(color="black", size = 1),
    legend.background=element_blank(),
    legend.text = element_text(size=16, face= "bold"),
    legend.title = element_text(size=16, face= "bold")
  )


pdf("Ordinal_reg.pdf",
    width = 17,height = 12)
Plot1 / p  + plot_layout(heights = c(3, 2))
dev.off()
