# Title: PLS-DA Control/Infect/Sepsis
# Project: Mitochondrial Adaption at the Center of Early Metabolic Changes in Presymptomatic Sepsis Patients
# Maintainer: Lin-Lin Xu
# Last update: Nov 24, 2024
# License: MIT

## Load library ----------------------------------------------------------------
library(tidyverse)
library(mixOmics)
library(ggprism)
library(vegan)
library(ggrepel)
library(patchwork)
jco_pal = c("#3C5488FF","#F39B7FFF","#DC0000FF")
## Main ------------------------------------------------------------------------
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
  dplyr::rename(metabolite = CompoundID) %>%
  mutate(Name = str_replace_all(Name, c("lysoPC" = "LPC", "lysoPE" = "LPE")))
# Meta data --------------------------------------------------------------------
meta_df <- read_csv2("mata_data1.csv") %>% 
  dplyr::select(.,-`Day post surgery before sepsis diagnosis or equivalent day post for comparator and SIRS Controls`) %>% 
  distinct(Classifier,Age,Gender, .keep_all = TRUE)
Clinical_df <- read_tsv("meta_data2.tsv") %>%
  mutate_at(7:length(.), as.numeric) %>%
  dplyr::select(.,-Group) 
data_imput_list <- readRDS("Clinical_data_imput_list.rds")
data_imput_df <- data_imput_list$ximp %>%
  as_tibble() %>%
  cbind(Clinical_df[,1:6],.) %>%
  pivot_longer(7:length(.),names_to = "parameter",values_to = "value") %>%
  filter(`Sample ID` %in% Met_file$`Sample ID`)
Clinical_mat <- data_imput_df %>%
  ungroup() %>%
  dplyr::select(.,-`Surgery day difference`,-`Ethnic Origin`,-`Admission Criteria`)
# Tidy data --------------------------------------------------------------------
Met_df <- Met_file %>%
  pivot_wider(names_from = "metabolite",values_from = "abundance") %>%
  mutate(Class2 = str_replace_all(Class2, c("SIRS-" = "Control"))) %>%
  mutate(Class2 = str_replace_all(Class2, c("Uncomplicated Infection" = "Un.Inf"))) %>%
  dplyr::filter(Class2 %in% c("Control","Un.Inf","Sepsis"))
Lip_df <- Lip_file %>%
  pivot_wider(names_from = "metabolite",values_from = "abundance") %>%
  mutate(Class2 = str_replace_all(Class2, c("SIRS-" = "Control"))) %>%
  mutate(Class2 = str_replace_all(Class2, c("Uncomplicated Infection" = "Un.Inf"))) %>%
  dplyr::filter(Class2 %in% c("Control","Un.Inf","Sepsis"))
same_patient <- intersect(Met_df$`Sample ID`,Lip_df$`Sample ID`)
Clinical_met_mat <- Clinical_mat %>%
  filter(`Sample ID` %in% Met_df$`Sample ID`) %>%
  pivot_wider(names_from = parameter,values_from = value) %>%
  column_to_rownames("Sample ID")
Clinical_lip_mat <- Clinical_mat %>%
  filter(`Sample ID` %in% Lip_df$`Sample ID`) %>%
  filter(`Sample ID` %in% same_patient) %>%
  pivot_wider(names_from = parameter,values_from = value) %>%
  column_to_rownames("Sample ID")

Lip_mat <- Lip_df %>%
  dplyr::select(.,-Class2) %>%
  column_to_rownames(var = "Sample ID")

Met_mat <- Met_df %>%
  dplyr::select(.,-Class2) %>%
  column_to_rownames(var = "Sample ID")
Ordinal_reg_all_result_list_BH <- readRDS("Ordinal_reg_BH.rds")
sign_metabolite <- Ordinal_reg_all_result_list_BH$Met_Ord$re_df_025$name
sign_lipid <- Ordinal_reg_all_result_list_BH$Lip_Ord$re_df_025$name
## Function --------------------------------------------------------------------
func_envfit <- function(pca_tmp,meta_data,Top_10){
  meta_data <- meta_data %>%
    dplyr::select(Top_10$name)
  fit <- envfit(pca_tmp,meta_data,perm=999)
  fit.df<- as.data.frame(fit$vectors$arrows*sqrt(fit$vectors$r))
  fit.df$name<-rownames(fit.df)
  fit.df$Rsquare <- fit$vectors$r
  fit.df$pvals <- fit$vectors$pvals
  fit.df$pvals.BH <- p.adjust(fit.df$pvals,method = "BH")
  fit.df.1 <- fit.df %>% 
    filter(pvals <= 0.05,pvals.BH <= 0.05)
  fit.df.2 <- fit.df %>% 
    filter(pvals >= 0.05)
  arrow_factor <- ordiArrowMul(fit)
  return_list <- list(
    df = fit.df.1,
    tmp = fit.df,
    arrow_factor = arrow_factor
  )
}
#* SIRS+ VS Sepsis -------------------------------------------------------------
#*** Metabolite data -----------------------------------------------------------
Log2_aft_met_sub_Y <- Met_df$Class2 %>%
  factor(.,levels = c("Control","Un.Inf","Sepsis"))
Log2_aft_met_sub_X <- Met_mat
Log2_aft_met_sub_X<- Log2_aft_met_sub_X %>%
  rownames_to_column(var = "rowname") %>%
  as_tibble(.name_repair = "unique") %>%
  column_to_rownames(var = "rowname")
X <- Log2_aft_met_sub_X
Y <- Log2_aft_met_sub_Y
groups <- c("Control","Un.Inf","Sepsis")
# mixOmics PLSDA ------------------------------------------------------------
set.seed(123)
srbct.plsda <- plsda(X, Y, ncomp = 10)
tmp <- plotLoadings(srbct.plsda, comp = 1,method = 'mean', contrib = 'max', plot = FALSE)  %>%
  arrange(rev(abs(importance)))
tb <- tmp %>% 
  rownames_to_column(var = "name") %>%
  dplyr::select(name,GroupContrib,importance) %>%
  as_tibble() 
loading_plot <- tb %>%
  left_join(Met_ID_info,by = c("name" = "metabolite")) %>%
  dplyr::select(.,-name) %>%
  dplyr::rename(name = Name) %>%
  head(10) %>%
  mutate(name = factor(name,levels = rev(name))) %>%
  arrange(name) %>%
  mutate(color = factor(GroupContrib,levels = c("Control","Un.Inf","Sepsis"))) %>%
  ggplot(aes(x = name,y = importance)) +
  geom_bar(aes(fill = color),stat = "identity") +
  scale_fill_manual(values = jco_pal) +
  labs(title="PLS-DA metabolites",x="",y="Contribution on comp 1") +
  theme_prism(base_size = 12) +
  coord_flip()

Top_10 <- tb %>%
  left_join(Met_ID_info,by = c("name" = "metabolite")) %>%
  head(10)
background = background.predict(srbct.plsda, comp.predicted=2, dist = "centroids.dist")

a <- plotIndiv(srbct.plsda, comp = 1:2,
               group = Y, ind.names = FALSE, # colour points by class
               background = background, # include prediction background for each class
               legend = TRUE, title = "")
plsda_df <- a$graph$data %>%
  mutate(x = as.numeric(x)) %>%
  mutate(y = as.numeric(y)) %>%
  mutate(group = str_replace_all(group, c("Uncomplicated Infection" = "Un.Inf")))
plsda_tmp <- data.frame(PC1=plsda_df$x,PC2=plsda_df$y)
Clinical_fit <- func_envfit(plsda_tmp,X,Top_10)
arrow_df <- Clinical_fit$df %>%
  left_join(Met_ID_info,by = c("name" = "metabolite")) %>%
  dplyr::select(.,-name) %>%
  dplyr::rename(name = Name)
tmp_df <- background$Control %>% as_tibble()
tmp_df1 <- background$`Un.Inf` %>% as_tibble()
tmp_df2 <- background$Sepsis %>% as_tibble()
PLSDA_Plot_metabolite <- ggplot(plsda_df)+geom_point(aes(x=x,y=y,color = group)) + 
  geom_polygon(aes(x = Var1,y = Var2),tmp_df,fill = "#3C5488FF",alpha = 0.1) +
  geom_polygon(aes(x = Var1,y = Var2),tmp_df1,fill = "#F39B7FFF",alpha = 0.1) +
  geom_polygon(aes(x = Var1,y = Var2),tmp_df2,fill = "#DC0000FF",alpha = 0.1) +
  geom_segment(data=arrow_df,aes(x=0,xend=PC1* Clinical_fit$arrow_factor * 10,
                                        y=0,yend=PC2* Clinical_fit$arrow_factor * 10,color = col),
               arrow = arrow(length = unit(0.2, "cm")),colour="black") +
  geom_text_repel(data=arrow_df,aes(x=PC1* Clinical_fit$arrow_factor * 10,
                                           y=PC2* Clinical_fit$arrow_factor * 10,label=name),
                  size=5,colour="black",fontface = "bold",max.overlaps = 10)+
  theme_prism(base_size = 16) + 
  scale_y_continuous(limits = c(min(tmp_df1$Var2,tmp_df$Var2,tmp_df2$Var2),max(tmp_df1$Var2,tmp_df$Var2,tmp_df2$Var2)),
                     expand = c(0,0),guide = guide_prism_minor()) +
  scale_x_continuous(limits = c(min(tmp_df1$Var1,tmp_df$Var1,tmp_df2$Var1),max(tmp_df1$Var1,tmp_df$Var1,tmp_df2$Var1)),
                     expand = c(0,0),guide = guide_prism_minor()) + 
  scale_color_manual(breaks=c("Control","Un.Inf","Sepsis"),
                     values = jco_pal) +
  labs(title="PLS-DA metabolites",x=a$graph$labels$x,y=a$graph$labels$y) +
  theme(plot.title.position = "plot",legend.text=element_text(size=16, face = "bold"))
metabolites_auroc <- auroc(srbct.plsda,roc.comp = 2)
#*** Lipids data ---------------------------------------------------------------
Log2_aft_lip_sub_Y <- Lip_df$Class2 %>%
  factor(.,levels = c("Control","Un.Inf","Sepsis"))
Log2_aft_lip_sub_X <- Lip_mat
Log2_aft_lip_sub_X<- Log2_aft_lip_sub_X %>%
  rownames_to_column(var = "rowname") %>%
  as_tibble(.name_repair = "unique") %>%
  column_to_rownames(var = "rowname")
X <- Log2_aft_lip_sub_X
Y <- Log2_aft_lip_sub_Y
groups <- c("Control","Un.Inf","Sepsis")
# mixOmics PLSDA ------------------------------------------------------------
set.seed(123)
srbct.plsda <- plsda(X, Y, ncomp = 10)
tmp <- plotLoadings(srbct.plsda, comp = 1,method = 'mean', contrib = 'max', plot = FALSE)  %>%
  arrange(rev(abs(importance)))
tb <- tmp %>% 
  rownames_to_column(var = "name") %>%
  dplyr::select(name,GroupContrib,importance) %>%
  as_tibble() 
loading_plot_lipid <- tb %>%
  left_join(Lip_ID_info,by = c("name" = "metabolite")) %>%
  dplyr::select(.,-name) %>%
  dplyr::rename(name = Name) %>%
  head(10) %>%
  mutate(name = factor(name,levels = rev(name))) %>%
  arrange(name) %>%
  mutate(color = factor(GroupContrib,levels = c("Control","Un.Inf","Sepsis"))) %>%
  ggplot(aes(x = name,y = importance)) +
  geom_bar(aes(fill = color),stat = "identity") +
  scale_fill_manual(values = jco_pal) +
  labs(title="PLS-DA lipids",x="",y="Contribution on comp 1") +
  theme_prism(base_size = 12) +
  theme(legend.position = "none") +
  coord_flip()

Top_10 <- tb %>%
  left_join(Lip_ID_info,by = c("name" = "metabolite")) %>%
  head(10)
background = background.predict(srbct.plsda, comp.predicted=2, dist = "centroids.dist")
a <- plotIndiv(srbct.plsda, comp = 1:2,
               group = Y, ind.names = FALSE, # colour points by class
               background = background, # include prediction background for each class
               legend = TRUE, title = "")
plsda_df <- a$graph$data %>%
  mutate(x = as.numeric(x)) %>%
  mutate(y = as.numeric(y)) %>%
  mutate(group = str_replace_all(group, c("Uncomplicated Infection" = "Un.Inf")))
plsda_tmp <- data.frame(PC1=plsda_df$x,PC2=plsda_df$y)
Clinical_fit_lip <- func_envfit(plsda_tmp,X,Top_10)
arrow_df_lip <- Clinical_fit_lip$df %>%
  left_join(Lip_ID_info,by = c("name" = "metabolite")) %>%
  dplyr::select(.,-name) %>%
  dplyr::rename(name = Name)

tmp_df <- background$Control %>% as_tibble()
tmp_df1 <- background$Un.Inf %>% as_tibble()
tmp_df2 <- background$Sepsis %>% as_tibble()
PLSDA_Plot_lipid <- ggplot(plsda_df)+geom_point(aes(x=x,y=y,color = group)) + 
  geom_polygon(aes(x = Var1,y = Var2),tmp_df,fill = "#3C5488FF",alpha = 0.1) +
  geom_polygon(aes(x = Var1,y = Var2),tmp_df1,fill = "#F39B7FFF",alpha = 0.1) +
  geom_polygon(aes(x = Var1,y = Var2),tmp_df2,fill = "#DC0000FF",alpha = 0.1) +
  geom_segment(data=arrow_df_lip,aes(x=0,xend=PC1* Clinical_fit_lip$arrow_factor * 3,
                                        y=0,yend=PC2* Clinical_fit_lip$arrow_factor * 3,color = col),
               arrow = arrow(length = unit(0.2, "cm")),colour="black") +
  geom_text_repel(data=arrow_df_lip,aes(x=PC1* Clinical_fit_lip$arrow_factor * 3,
                                           y=PC2* Clinical_fit_lip$arrow_factor * 3,label=name),
                  size=5,colour="black",fontface = "bold",max.overlaps = 30)+
  theme_prism(base_size = 16) + 
  scale_y_continuous(limits = c(min(tmp_df1$Var2,tmp_df$Var2,tmp_df2$Var2),max(tmp_df1$Var2,tmp_df$Var2,tmp_df2$Var2)),
                     expand = c(0,0),guide = guide_prism_minor()) +
  scale_x_continuous(limits = c(min(tmp_df1$Var1,tmp_df$Var1,tmp_df2$Var1),max(tmp_df1$Var1,tmp_df$Var1,tmp_df2$Var1)),
                     expand = c(0,0),guide = guide_prism_minor()) + 
  scale_color_manual(breaks=c("Control","Un.Inf","Sepsis"),
                     values = jco_pal) +
  labs(title="PLS-DA lipids",x=a$graph$labels$x,y=a$graph$labels$y) +
  theme(plot.title.position = "plot",legend.text=element_text(size=16, face = "bold"))
lipid_auroc <- auroc(srbct.plsda,roc.comp = 2)
## Assembly plot  --------------------------------------------------------------
Loadings <- loading_plot + loading_plot_lipid + plot_layout(guides = "collect")
PLSDA <- PLSDA_Plot_metabolite + PLSDA_Plot_lipid + plot_layout(guides = "collect")

auroc <- metabolites_auroc$graph.Comp2 + lipid_auroc$graph.Comp2 + 
  plot_annotation(tag_levels = list(c('metabolite', 'lipid')))
auroc <- metabolites_auroc$graph.Comp2 + 
  theme_bw() + 
  theme(legend.position = "none") +
  lipid_auroc$graph.Comp2 + 
  theme_bw() + 
  theme(legend.position = "none")
pdf("FigS3_PLSDA_aurocs.pdf",
    width = 6,height = 3)
auroc
dev.off()
pdf("Fig3_PLSDA.pdf",width = 16,height = 8)
PLSDA
dev.off()
pdf("FigS3_PLSDA_loadings.pdf",
    width = 9,height = 3)
Loadings
dev.off()