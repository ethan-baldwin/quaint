find_column <- function(taxa_pair,taxa_quartet) {
  if (all(taxa_pair %in% taxa_quartet[1:2])|all(taxa_pair %in% taxa_quartet[3:4])) {
    return("12|34") 
  } else if (all(taxa_pair %in% taxa_quartet[2:3])|all(taxa_pair %in% c(taxa_quartet[1],taxa_quartet[4]))) {
    return("14|23")
  } else {return("13|24")}
  
}

parse_table <- function(position,combos,quartet_table,st) {
  current_combo <- combos[position,]
  ######### get order of taxa for quartet table ##########
  o <- st$tip.label[combos[position,1]]
  p1 <- st$tip.label[combos[position,2]]
  p2 <- st$tip.label[combos[position,3]]
  p3 <- st$tip.label[combos[position,4]]
  
  o <- st$tip.label[current_combo[[1]]]
  p1 <- st$tip.label[current_combo[[2]]]
  p2 <- st$tip.label[current_combo[[3]]]
  p3 <- st$tip.label[current_combo[[4]]]
  
  conc_taxa <- c(p2,p3)
  test_taxa <- c(p1,p2)
  
  qdf <- as.data.frame(quartet_table)
  sub_qt <- qdf[qdf[,p1]==1&qdf[,p2]==1&qdf[,p3]==1&qdf[,o]==1,]
  
  qt_taxa <- names(sub_qt)[which(sub_qt == 1, arr.ind=T)[, "col"]]
  v<-c("12|34","13|24","14|23") # vector of column names to exclude when determining new taxon's position
  qt_taxa <- setdiff(qt_taxa,v)
  
  # set column definitions
  conc_column <- find_column(conc_taxa,qt_taxa)
  abba_column <- find_column(test_taxa,qt_taxa)
  baba_column <- setdiff(c("12|34","13|24","14|23"),c(abba_column,conc_column))
  
  # get values from columns
  abba <- sub_qt[,abba_column]
  baba <- sub_qt[,baba_column]
  conc <- sub_qt[,conc_column]
  total <- abba+baba+conc
  d_stat <- (abba-baba)/(abba+baba)
  
  # add data to return df
  add2df <- c(o,p1,p2,p3,abba,baba,conc,total,d_stat)
  add2df
}

dstat_table <- function(test_taxa,species_tree,quartet_table,alpha=0.05) {
  mrca <- getMRCA(species_tree, test_taxa) 
  node_groups <- nodeGroups(species_tree,mrca) # create node groups to determine where taxon is
  t1 <- test_taxa[1]
  t2 <- test_taxa[2]
  # get tip number from tip label
  t1_num <- which(species_tree$tip.label==t1)
  t2_num <- which(species_tree$tip.label==t2)
  # get node group
  t1_group <- which(sapply(node_groups, is.element, el = t1_num))
  t2_group <- which(sapply(node_groups, is.element, el = t2_num))
  
  # get tip numbers for all outgroups and sister taxa
  sister_tips_t1 <- setdiff(node_groups[[t1_group]],c(t1_num)) 
  sister_tips_t2 <- setdiff(node_groups[[t2_group]],c(t2_num))
  og_tips <- node_groups[[setdiff(c(1,2,3),c(t1_group,t2_group))]]
  
  # check if taxa are sister
  if (length(sister_tips_t1) == 0 & length(sister_tips_t2) == 0) {
    stop("Test taxa are sister taxa and cannot be tested for introgression.")
  }
  
  # get all test combinations
  combos_1 <- expand.grid(og_tips,t1_num,t2_num,sister_tips_t2)
  combos_2 <- expand.grid(og_tips,t2_num,t1_num,sister_tips_t1)
  all_combos <- rbind(combos_1,combos_2)
  
  # parse quartet frequency df and prepare output data frame 
  num_rows <- 1:nrow(all_combos) # get number of rows
  data_list <- lapply(num_rows,parse_table,combos=all_combos,quartet_table=quartet_table,st=species_tree) # call function that parses quartet freq df
  return_df <- as.data.frame(do.call(rbind, data_list)) # convert output from lapply (list) into data frame
  cnames<-c(  "outgroup","taxon1","taxon2","taxon3",
              "abba","baba","concordant","total","d") # appropriate column names
  colnames(return_df) <- cnames # add column names to df
  return_df <- return_df %>% mutate_at(c("abba","baba","concordant","total","d"), as.numeric) # convert columns to numeric as approprtiate
  
  #perform chi squared tests
  # return_df <- rbind(return_df,data.frame(outgroup="Total",taxon1=t1,taxon2=t2,taxon3=NA,d=NA,t(colSums(return_df[,5:8]))))
  return_df <- return_df %>%
    rowwise() %>% 
    mutate(
      test_stat = suppressWarnings(chisq.test(c(abba,baba))$statistic),
      p_val = suppressWarnings(chisq.test(c(abba,baba))$p.value)) # perform chi squared test on all rows, including total row
  
  n_tests <- nrow(return_df)
  n_positive <- sum(return_df$d>0,na.rm = TRUE)
  n_positive_significant <- nrow(return_df[(return_df$d>0&return_df$p_val<alpha),])
  total_df <- data.frame(outgroup="Total",taxon1=t1,taxon2=t2,taxon3=NA,d=n_tests,t(colSums(return_df[,5:8])),test_stat=n_positive,p_val=n_positive_significant)
  return_df <- rbind(return_df,total_df) # add a row with totals for all topology columns
  # summary_df <- data.frame(n_total_tests = n_tests,n_positive_d = n_positive,n_positive_significant_tests = n_positive_significant)
  return_df
}

dstat_table_all <- function(test_taxa_list,species_tree,quartet_table,alpha=0.05) {
  pos <- 1:ncol(test_taxa_list)
  # test_taxa_list[,position]
  taxon_pair_df <- lapply(pos,function(x) dstat_table(allpairs[,x],rooted_st,quartet_table,alpha=alpha))
  taxon_pair_df <- as.data.frame(do.call(rbind, taxon_pair_df))
  # taxon_pair_df <- taxon_pair_df[taxon_pair_df$outgroup=="Total",]
  # taxon_pair_df <- subset(taxon_pair_df,select = -c(outgroup,taxon3))
  taxon_pair_df
}

get_valid_test_pairs <- function(species_tree,outgroup) {
  # remove outgroups
  taxa <- setdiff(species_tree$tip.label,outgroup)
  # get all pairwise combos
  pairwise_combos <- combn(taxa,2,simlify=FALSE)
  n_combos <- ncol(pairwise_combos)
  pairwise_combos[,sapply(1:n_combos,is_taxon_pair_valid,pairwise_combos,species_tree)]
}

is_taxon_pair_valid <- function(position,p_combos,tree) {
  taxa_pair <- p_combos[,position]
  mrca <- getMRCA(tree, taxa_pair) 
  node_groups <- nodeGroups(tree,mrca) # create node groups to determine where taxon is
  t1 <- taxa_pair[1]
  t2 <- taxa_pair[2]
  
  # get tip number from tip label
  t1_num <- which(tree$tip.label==t1)
  t2_num <- which(tree$tip.label==t2)
  
  # get node group
  t1_group <- which(sapply(node_groups, is.element, el = t1_num))
  t2_group <- which(sapply(node_groups, is.element, el = t2_num))
  
  # get tip numbers for all outgroups and sister taxa
  sister_tips_t1 <- setdiff(node_groups[[t1_group]],c(t1_num)) 
  sister_tips_t2 <- setdiff(node_groups[[t2_group]],c(t2_num))
  og_tips <- node_groups[[setdiff(c(1,2,3),c(t1_group,t2_group))]]
  
  !(length(sister_tips_t1)==0 & (length(og_tips)==0)) & !(length(sister_tips_t1)==0 & length(sister_tips_t2)==0) & !(length(sister_tips_t2)==0 & (length(og_tips)==0))
}
