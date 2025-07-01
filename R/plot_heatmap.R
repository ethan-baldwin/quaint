# library(ggplot2)

#' Plotting function WIP
#' @param quaint_results Data frame with summarized quaint results.
#' @param tree Species tree.
#' @param outgroup Outgroup tip names.
#' @param plot_val Value you want to plot.
#' @param cell_text Boolean. Whether you want text in the heatmap cells or not.
#' @noRd
plot_quaint_heatmap <- function(quaint_results,tree,outgroup,plot_val,cell_text=FALSE) {
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
  
  # plot heatmap
  heatmap<-ggplot(data=combos_with_position,aes(x=tax1,y=tax2,fill=stat))+
    geom_tile(color="black")+
    scale_fill_gradient2(low = "white", high = "red",na.value = "grey50",
                         name=plot_val) +
    theme_minimal()+
    theme(axis.text.x=element_text(angle=90,hjust = 0),
          axis.text = element_text(size = 20,face = "bold",color = "black"),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          axis.ticks = element_line(color = "black"),
          legend.key.height = unit(1, 'inch'))+
    scale_x_discrete(position="top") 
  # add text to the middle of cells if wanted
  if(cell_text==TRUE) {heatmap <- heatmap + geom_text(aes(x=tax1,y=tax2,label=stat_rounded),color="black",size=7)}
  
  
  heatmap
}

#' Plotting function WIP
#' @param quaint_results Data frame with summarized quaint results.
#' @param tree Species tree.
#' @param outgroup Outgroup tip names.
#' @param plot_val Value you want to plot.
#' @param cell_text Boolean. Whether you want text in the heatmap cells or not.
#' @noRd
get_plot_df <- function(quaint_results,tree,outgroup,plot_val,cell_text=FALSE) {
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
  
  return(combos_with_position)
}
