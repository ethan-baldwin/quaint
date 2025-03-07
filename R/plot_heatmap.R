library(ggplot2)

plot_quaint_heatmap <- function(quaint_results,tree,outgroup,plot_val) {
  tip_labels<-setdiff(tree$tip.label,outgroup)
  combos<-expand.grid(tip_labels,tip_labels)
  
  quaint_results_reversed<-quaint_results%>%
    rename(Var1=taxon2,Var2=taxon1)
  
  quaint_results<-quaint_results%>%
    rename(Var1=taxon1,Var2=taxon2)
  
  all_quaint_results<-bind_rows(quaint_results,quaint_results_reversed)
  
  combos_with_position<-combos%>%
    left_join(all_quaint_results,by=c("Var1","Var2"))
  
  # order of tip labels
  order <- c("S_purpurea_mont",
             "S_rosea_4",
             "S_minor_2",
             "S_minor_oke",
             "S_flava_2",
             "S_psittacina",
             "S_psittacina_oke",
             "S_oreophila",
             "S_leucophylla",
             "S_rubra",
             "S_jonesii",
             "S_rubra_gulfensis",
             "S_alata",
             "S_alabamensis",
             "S_rubra_wherryi")
  
  # vector of axis tip colors
  color_order <-   ifelse(grepl("jonesii|rubra|alabamensis|leucophylla|alata|oreophila",order, ignore.case = TRUE,fixed = FALSE), "#FF0000",
                          ifelse(grepl("purpurea|rosea|flava|minor|psittacina",order, ignore.case = TRUE,fixed = FALSE), "#4472C4",
                                 "Outgroup"))
  
  # shorten names
  combos_with_position$Var1 <- gsub("okefenokeensis","oke",combos_with_position$Var1,ignore.case = T)
  combos_with_position$Var1 <- gsub("_RVCAGA001","",combos_with_position$Var1,ignore.case = T)
  combos_with_position$Var1 <- gsub("_m003","",combos_with_position$Var1,ignore.case = T)
  combos_with_position$Var1 <- gsub("_m002","",combos_with_position$Var1,ignore.case = T)
  combos_with_position$Var1 <- gsub("_P008","",combos_with_position$Var1,ignore.case = T)
  combos_with_position$Var1 <- gsub("montana","mont",combos_with_position$Var1,ignore.case = T)
  
  combos_with_position$Var2 <- gsub("okefenokeensis","oke",combos_with_position$Var2,ignore.case = T)
  combos_with_position$Var2 <- gsub("_RVCAGA001","",combos_with_position$Var2,ignore.case = T)
  combos_with_position$Var2 <- gsub("_m003","",combos_with_position$Var2,ignore.case = T)
  combos_with_position$Var2 <- gsub("_m002","",combos_with_position$Var2,ignore.case = T)
  combos_with_position$Var2 <- gsub("_P008","",combos_with_position$Var2,ignore.case = T)
  combos_with_position$Var2 <- gsub("montana","mont",combos_with_position$Var2,ignore.case = T)
  
  # round proportion to two decimal places
  # combos_with_position$mean_d2<-round(combos_with_position$mean_d,2)
  
  # set order of axes
  # combos_with_position$tax1<-factor(combos_with_position$Var1,levels=order)
  combos_with_position$tax1<-factor(combos_with_position$Var1,levels=rev(order))
  combos_with_position$tax2<-factor(combos_with_position$Var2,levels=order)
  # combos_with_position$tax2<-factor(combos_with_position$Var2,levels=rev(order))
  
  combos_with_position$stat <- combos_with_position[[plot_val]]
  combos_with_position$stat_rounded<-round(combos_with_position$stat,2)
  # plot!
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
    scale_x_discrete(position="top") +
    geom_text(aes(x=tax1,y=tax2,label=stat_rounded),color="black",size=7)
  
  heatmap
  
  # add tip colors
  heatmap + theme(axis.text.x = element_text(color = rev(color_order)),
                  axis.text.y = element_text(color = color_order))
}

# quaint_results<-read.csv("quaint_summary.csv",header=TRUE)
# sp_tree <- read.tree("astral3_rooted.tre")
# outgroup <- c("heliamphora_pulchella_SRR25244091","heliamphora_ciliata_SRR24877724") 
# 
# plot_quaint_heatmap(quaint_results,sp_tree,outgroup)
# 
# 
# ggsave("heatmap_mean_d.png",width = 14,height = 12)

