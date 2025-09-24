# library(ggplot2)

#' Plotting function WIP
#' @param quaint_results Data frame with summarized quaint results.
#' @param tree Species tree.
#' @param outgroup Outgroup tip names. 
#' @param plot_val Value you want to plot from the summarized results. Default is "mean_d".
#' @importFrom dplyr left_join %>% rename bind_rows
#' @importFrom ape compute.brlen
#' @importFrom circlize colorRamp2
#' @importFrom ComplexHeatmap Heatmap
#' @noRd
quaint_heatmap <- function(quaint_results,tree,outgroup,plot_val="mean_d") {
  tip_labels<-setdiff(tree$tip.label,outgroup)
  combos<-expand.grid(tip_labels,tip_labels)
  
  quaint_results_reversed<-quaint_results%>%
    rename(Var1=taxon2,Var2=taxon1)
  
  quaint_results<-quaint_results%>%
    rename(Var1=taxon1,Var2=taxon2)
  
  all_quaint_results<-bind_rows(quaint_results,quaint_results_reversed)
  
  combos_with_position<-combos%>%
    left_join(all_quaint_results,by=c("Var1","Var2"))
  
  # get order of tip labels
  edge_order <- tree$edge[,2]
  n_tips_vector <- 1:length(tree$tip.label)
  vector_order <- edge_order[edge_order %in% n_tips_vector]
  tip_order <- rev(tree$tip.label[vector_order])
  
  # set order of axes
  combos_with_position$tax1<-factor(combos_with_position$Var1,levels=rev(tip_order))
  combos_with_position$tax2<-factor(combos_with_position$Var2,levels=tip_order)
  
  combos_with_position$stat <- combos_with_position[[plot_val]]
  combos_with_position$stat_rounded<-round(combos_with_position$stat,2)
  
  ultra <- compute.brlen(tree, method = "Grafen", power = 0.5)
  stat_mat <- with(combos_with_position, tapply(stat, list(tax2, tax1), FUN = identity))
  hc <- as.hclust(ultra,use.edge.length = TRUE)
  dend <- as.dendrogram(hc) 

  col_fun <- colorRamp2(
    c(0,min(stat_mat[stat_mat>0],na.rm = TRUE), max(stat_mat,na.rm = TRUE)),
    c("white","#FFE6E6","red")
  )
  lab_order <- hc$labels
  stat_mat <- stat_mat[lab_order, lab_order]
  
  h <- Heatmap(
    stat_mat,
    col             = col_fun,
    name                 = plot_val,
    cluster_rows         = hc,
    cluster_columns      = hc,
    rect_gp              = grid::gpar(col = "black", lwd = 0.5),
    border               = TRUE,
    
    ## keep labels next to their dendrograms
    row_names_side       = "left",
    row_dend_side        = "left",
    column_names_side    = "top",
    column_dend_side     = "top",
    column_dend_height = unit(15, "mm"),
    column_names_max_height = unit(8, "cm"),
    row_dend_width = unit(15, "mm"),
    row_names_max_width = unit(8, "cm"))
  h
}
