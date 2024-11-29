# Title: Functions
# Project: Mitochondrial Adaption at the Center of Early Metabolic Changes in Presymptomatic Sepsis Patients
# Maintainer: Lin-Lin Xu
# Last update: Nov 24, 2024
# License: MIT

## Partial correlation  --------------------------------------------------------
pcor_test <- function(df,x,y,p,q){
  df <- df %>%
    filter(across(where(is.numeric), ~ !is.na(.)))
  a <- df[[x]]
  b <- df[[y]]
  c <- df[[p]]
  d <- df[[q]] %>% as.factor() %>% as.integer()
  r <- tryCatch(
    {
      as_tibble(pcor.test(a,b,c(c,d),method = "spearman"))
    },
    error = function(e){
      data.frame(estimate = NA,
                 p.value = NA,
                 statistic = NA,
                 n = NA,
                 gp = NA,
                 Method = NA,
                 q.value = NA)
    }
  )
}
## calculate correlation 1 (normalized mean) -----------------------------------
Cal_cor <- function(df,df_info){

  Cor_met_cli_df <- df %>%
    ungroup() %>%
    dplyr::select(metabolite,meta_parameter,abundance,meta_value,Age,Sex) %>%
    group_by(metabolite,meta_parameter) %>%
    group_modify(~ pcor_test(.x,"abundance","meta_value","Age","Sex")) %>% 
    bind_rows() %>%
    ungroup() %>%
    mutate(q.value = p.adjust(p.value, method = "BH")) %>%
    left_join(df_info)
}
Organize_cor_re <- function(data_lst){
  
  Cor_met_cli_df <- data_lst$Cor_met_cli_df
  Cor_lip_cli_df <- data_lst$Cor_lip_cli_df
  
  Combine_all <- rbind(Cor_met_cli_df %>% mutate(type = "metabolite"),
                       Cor_lip_cli_df %>% mutate(type = "lipid"))
  
  Double_name <- Combine_all %>%
    filter(meta_parameter %in% "APPT") %>%
    .$Name %>%
    table() %>%
    as.data.frame() %>%
    `colnames<-`(c("Name","Freq")) %>%
    filter(Freq > 1)
  
  Plot_df <- Combine_all %>%
    mutate(Name = if_else(Name %in% Double_name$Name,paste0(Name,"_",metabolite),Name)) %>%
    mutate(Name = str_replace_all(Name, c("lysoPC" = "LPC", "lysoPE" = "LPE")))
  
  r <- Plot_df %>%
    dplyr::select(Name,meta_parameter,estimate) %>%
    pivot_wider(names_from = meta_parameter,values_from = estimate) %>%
    column_to_rownames(var = "Name") %>%
    as.matrix()
  
  p <- Plot_df %>%
    dplyr::select(Name,meta_parameter,p.value) %>%
    pivot_wider(names_from = meta_parameter,values_from = p.value) %>%
    column_to_rownames(var = "Name") %>%
    as.matrix()
  
  p.BH <- Plot_df %>%
    dplyr::select(Name,meta_parameter,q.value) %>%
    pivot_wider(names_from = meta_parameter,values_from = q.value) %>%
    column_to_rownames(var = "Name") %>%
    as.matrix()
  
  r[is.na(r)] <- 0
  p[is.na(p)] <- 1
  p.BH[is.na(p.BH)] <- 1
  
  re_list <- list(
    Plot_df = Plot_df,
    r = r,
    p = p,
    p.BH = p.BH
  )
  return(re_list)
}
Draw_cor_heatmap <- function(data_lst,meta_direction,Sign_signatures,class_level){
  r = data_lst$r
  p = data_lst$p
  p.BH = data_lst$p.BH
  Plot_df = data_lst$Plot_df
  
  column_df <- data.frame(name = colnames(r)) %>%
    left_join(meta_direction %>%
                dplyr::rename(name = parameter)) %>%
    dplyr::select(name,Sign) %>%
    column_to_rownames("name") %>%
    rename(group = Sign) %>%
    mutate(group = factor(group,levels = c("Sepsis","NonSepsis","neutral"))) %>%
    as.vector()
  columnPalette <- c("NonSepsis" = "#3C5488FF","Sepsis"= "#DC0000FF","neutral"= "grey")
  bar <- columnAnnotation(df = column_df,col = list(group =columnPalette),
                          show_annotation_name = F,
                          annotation_legend_param = list(title = "",nrow = 1))
  
  row_df <- data.frame(name = rownames(r)) %>%
    left_join(Plot_df %>% 
                dplyr::select(Name,metabolite) %>%
                dplyr::rename(name = Name)) %>%
    left_join(Sign_signatures %>% dplyr::select(.,-estimate)) %>%
    distinct() %>%
    dplyr::select(name,high_abundance,Signature) %>%
    column_to_rownames("name") %>%
    rename(group = high_abundance) %>%
    mutate(group = factor(group,levels = c("Sepsis","NonSepsis"))) %>%
    as.vector()
  rowPalette_group <- c("NonSepsis" = "#3C5488FF","Sepsis"= "#DC0000FF")
  rowPalette_test <- c("Ord" = "black")
  
  foo <- rowAnnotation(df = row_df,col = list(group = rowPalette_group,
                                              Signature = rowPalette_test),
                       show_legend = F,show_annotation_name = F,
                       annotation_legend_param = list(title = ""))
  row_df1 <- data.frame(name = rownames(r)) %>%
    left_join(Plot_df %>% 
                dplyr::select(Name,metabolite) %>%
                dplyr::rename(name = Name)) %>%
    left_join(Sign_signatures %>% dplyr::select(.,-estimate)) %>%
    distinct() %>%
    column_to_rownames("name") %>%
    dplyr::select(Compound_class) %>%
    as.vector()
  
  col_fun = circlize::colorRamp2(c(-1, 0, 1),  c("#01befe", "white", "#ff7d00"))
  Lip_p <- Heatmap(r, name = "correlation", 
                   col = col_fun, 
                   cell_fun = function(j, i, x, y, w, h, fill) {
                     if(p[i,j] <= 0.05) {
                       grid.points(x , y , pch = 1, size = unit(2, "mm"), gp = gpar(col = "black",lwd=2))
                     }
                     if(p.BH[i,j] <= 0.05) {
                       grid.points(x , y , pch = 2, size = unit(4, "mm"), gp = gpar(col = "black",lwd=2))
                     }
                   },
                   row_names_max_width = max_text_width(rownames(r),gp = gpar(fontsize = 12)),
                   column_names_max_height = max_text_width(colnames(r),gp = gpar(fontsize = 12)),
                   row_gap = unit(5, "mm"),column_gap = unit(5, "mm"),
                   row_split = factor(row_df1$Compound_class,levels = class_level),
                   row_title_rot = 0,
                   top_annotation = bar, right_annotation = foo,
                   cluster_rows = FALSE, cluster_columns = FALSE,
                   column_names_side = "top",column_names_rot = 45, 
                   show_row_names = T, show_column_names = T,
                   heatmap_legend_param = list(direction = "horizontal")
  )
}
Organize_cor_re_class <- function(data_lst){
  
  Cor_met_cli_df <- data_lst$Cor_met_cli_df
  Cor_lip_cli_df <- data_lst$Cor_lip_cli_df
  
  Combine_all <- rbind(Cor_met_cli_df %>% mutate(type = "metabolite"),
                       Cor_lip_cli_df %>% mutate(type = "lipid"))
  
  Plot_df <- Combine_all 
  
  r <- Plot_df %>%
    dplyr::select(metabolite,meta_parameter,estimate) %>%
    pivot_wider(names_from = meta_parameter,values_from = estimate) %>%
    column_to_rownames(var = "metabolite") %>%
    as.matrix()
  
  p <- Plot_df %>%
    dplyr::select(metabolite,meta_parameter,p.value) %>%
    pivot_wider(names_from = meta_parameter,values_from = p.value) %>%
    column_to_rownames(var = "metabolite") %>%
    as.matrix()
  
  p.BH <- Plot_df %>%
    dplyr::select(metabolite,meta_parameter,q.value) %>%
    pivot_wider(names_from = meta_parameter,values_from = q.value) %>%
    column_to_rownames(var = "metabolite") %>%
    as.matrix()
  
  r[is.na(r)] <- 0
  p[is.na(p)] <- 1
  p.BH[is.na(p.BH)] <- 1
  
  re_list <- list(
    Plot_df = Plot_df,
    r = r,
    p = p,
    p.BH = p.BH
  )
  return(re_list)
}
Draw_cor_heatmap_class <- function(data_lst,meta_direction,Sign_signatures,class_level){
  r = data_lst$r
  p = data_lst$p
  p.BH = data_lst$p.BH
  Plot_df = data_lst$Plot_df
  
  most_frequent <- function(x) {
    uniq_x <- unique(x)
    uniq_x[which.max(tabulate(match(x, uniq_x)))]
  }
  class_direction <- Sign_signatures %>%
    group_by(Compound_class) %>%
    mutate(high_abundance = most_frequent(high_abundance)) %>%
    ungroup() %>%
    dplyr::select(Compound_class,high_abundance) %>%
    distinct()
  
  
  column_df <- data.frame(name = colnames(r)) %>%
    left_join(meta_direction %>%
                dplyr::rename(name = parameter)) %>%
    dplyr::select(name,Sign) %>%
    column_to_rownames("name") %>%
    rename(group = Sign) %>%
    mutate(group = factor(group,levels = c("Sepsis","NonSepsis","neutral"))) %>%
    as.vector()
  columnPalette <- c("NonSepsis" = "#3C5488FF","Sepsis"= "#DC0000FF","neutral"= "grey")
  bar <- columnAnnotation(df = column_df,col = list(group =columnPalette),
                          show_annotation_name = F,
                          annotation_legend_param = list(title = "",nrow = 1))
  
  row_df <- data.frame(name = rownames(r)) %>%
    left_join(class_direction %>% 
                dplyr::rename(name = Compound_class)) %>%
    distinct() %>%
    column_to_rownames("name") %>%
    rename(group = high_abundance) %>%
    mutate(group = factor(group,levels = c("Sepsis","NonSepsis"))) %>%
    as.vector()
  rowPalette_group <- c("NonSepsis" = "#3C5488FF","Sepsis"= "#DC0000FF")

  foo <- rowAnnotation(df = row_df,col = list(group = rowPalette_group),
                       show_legend = F,show_annotation_name = F,
                       annotation_legend_param = list(title = ""))
  
  row_df1 <- data.frame(name = rownames(r)) %>%
    left_join(Plot_df %>% 
                dplyr::select(metabolite,type) %>%
                dplyr::rename(name = metabolite)) %>%
    distinct() %>%
    column_to_rownames("name") %>%
    dplyr::select(type) %>%
    as.vector()
  
  col_fun = circlize::colorRamp2(c(-1, 0, 1),  c("#01befe", "white", "#ff7d00"))
  Lip_p <- Heatmap(r, name = "correlation", 
                   col = col_fun, 
                   cell_fun = function(j, i, x, y, w, h, fill) {
                     if(p[i,j] <= 0.05) {
                       grid.points(x , y , pch = 1, size = unit(2, "mm"), gp = gpar(col = "black",lwd=2))
                     }
                     if(p.BH[i,j] <= 0.05) {
                       grid.points(x , y , pch = 2, size = unit(4, "mm"), gp = gpar(col = "black",lwd=2))
                     }
                   },
                   row_names_max_width = max_text_width(rownames(r),gp = gpar(fontsize = 12)),
                   column_names_max_height = max_text_width(colnames(r),gp = gpar(fontsize = 12)),
                   row_gap = unit(5, "mm"),column_gap = unit(5, "mm"),
                   row_split = factor(row_df1$type,levels = c("metabolite","lipid")),
                   row_title_rot = 0,
                   top_annotation = bar, right_annotation = foo,
                   cluster_rows = FALSE, cluster_columns = FALSE,
                   column_names_side = "top",column_names_rot = 45, 
                   show_row_names = T, show_column_names = T,
                   heatmap_legend_param = list(direction = "horizontal")
  )
}
Draw_cor_heatmap_reduced <- function(data_lst,meta_direction,Sign_signatures,class_level){
  r = data_lst$r %>% t()
  p = data_lst$p %>% t()
  p.BH = data_lst$p.BH %>% t()
  Plot_df = data_lst$Plot_df
  
  column_df <- data.frame(name = colnames(r)) %>%
    left_join(Plot_df %>% 
                dplyr::select(Name,metabolite) %>%
                dplyr::rename(name = Name)) %>%
    left_join(Sign_signatures %>% dplyr::select(.,-estimate)) %>%
    distinct() %>%
    dplyr::select(name,high_abundance) %>%
    column_to_rownames("name") %>%
    rename(group = high_abundance) %>%
    mutate(group = str_replace_all(group, c("^Sepsis$" = "Up Sepsis", "^NonSepsis$" = "Down Sepsis"))) %>%
    mutate(group = factor(group,levels = c("Up Sepsis","Down Sepsis"))) %>%
    as.vector()
  columnPalette_group <- c("Down Sepsis" = "#3C5488FF","Up Sepsis"= "#DC0000FF")
  bar <- columnAnnotation(df = column_df,col = list(group =columnPalette_group),
                          show_annotation_name = F,
                          annotation_legend_param = list(title = "",nrow = 1,
                                                         labels_gp = gpar(fontsize = 16,fontface = "bold"),
                                                         title_gp = gpar(fontsize = 16,fontface = "bold")
                                                         )
                          )
  
  column_df1 <- data.frame(name = colnames(r)) %>%
    left_join(Plot_df %>% 
                dplyr::select(Name,metabolite) %>%
                dplyr::rename(name = Name)) %>%
    left_join(Sign_signatures %>% dplyr::select(.,-estimate)) %>%
    distinct() %>%
    column_to_rownames("name") %>%
    dplyr::select(Compound_class) %>%
    mutate(Compound_class = str_replace_all(Compound_class, c("Amino acids and derivatives" = "Amino acids\nand derivatives"))) %>%
    as.vector()
  
  row_df <- data.frame(name = rownames(r)) %>%
    left_join(meta_direction %>%
                dplyr::rename(name = parameter)) %>%
    dplyr::select(name,Sign) %>%
    column_to_rownames("name") %>%
    rename(group = Sign) %>%
    mutate(group = factor(group,levels = c("Sepsis","NonSepsis","neutral"))) %>%
    as.vector()
  rowPalette <- c("NonSepsis" = "#3C5488FF","Sepsis"= "#DC0000FF","neutral"= "grey")
  
  foo <- rowAnnotation(df = row_df,col = list(group = rowPalette),
                       show_legend = F,show_annotation_name = F,
                       annotation_legend_param = list(title = "",
                                                      labels_gp = gpar(fontsize = 16,fontface = "bold"),
                                                      title_gp = gpar(fontsize = 16,fontface = "bold")))
  
  col_fun = circlize::colorRamp2(c(-1, 0, 1),  c("#01befe", "white", "#ff7d00"))
  Lip_p <- Heatmap(r, name = "correlation", 
                   col = col_fun, 
                   cell_fun = function(j, i, x, y, w, h, fill) {
                     if(p[i,j] <= 0.05) {
                       grid.points(x , y , pch = 1, size = unit(2, "mm"), gp = gpar(col = "black",lwd=2))
                     }
                     if(p.BH[i,j] <= 0.05) {
                       grid.points(x , y , pch = 2, size = unit(4, "mm"), gp = gpar(col = "black",lwd=2))
                     }
                   },
                   column_names_gp = grid::gpar(fontsize = 18,fontface = "bold"),
                   row_names_gp = grid::gpar(fontsize = 18,fontface = "bold"),
                   row_names_max_width = max_text_width(rownames(r),gp = gpar(fontsize = 18,fontface = "bold")),
                   column_names_max_height = max_text_width(colnames(r),gp = gpar(fontsize = 18,fontface = "bold")),
                   row_gap = unit(5, "mm"),column_gap = unit(5, "mm"),
                   column_split = factor(column_df1$Compound_class,levels = class_level),
                   column_title_gp = gpar(fontsize = 20, fontface = "bold"),
                   column_title_rot = 0,
                   bottom_annotation = bar, left_annotation = foo,
                   cluster_rows = FALSE, cluster_columns = FALSE,
                   row_names_side = "left",
                   column_names_side = "bottom",column_names_rot = 45, 
                   show_row_names = T, show_column_names = T,
                   heatmap_legend_param = list(direction = "horizontal",
                                               labels_gp = gpar(fontsize = 16,fontface = "bold"),
                                               title_gp = gpar(fontsize = 16,fontface = "bold"))
  )
}
