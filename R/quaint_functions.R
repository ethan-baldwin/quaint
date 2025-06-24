library("MSCquartets")
library("combinat")
library("dplyr")
library("stringr")
library("tidyr")

find_column <- function(taxa_pair,taxa_quartet) {
  if (all(taxa_pair %in% taxa_quartet[1:2])|all(taxa_pair %in% taxa_quartet[3:4])) {
    return("12|34")
  } else if (all(taxa_pair %in% taxa_quartet[2:3])|all(taxa_pair %in% c(taxa_quartet[1],taxa_quartet[4]))) {
    return("14|23")
  } else {return("13|24")}

}

# parse_table <- function(position,combos,quartet_df,st) {
parse_table <- function(position,combos,quartet_table,st) {
  current_combo <- combos[position,]
  ######### get order of taxa for quartet table ##########
  # get names of each taxon
  o <- st$tip.label[current_combo[[1]]]
  p1 <- st$tip.label[current_combo[[2]]]
  p2 <- st$tip.label[current_combo[[3]]]
  p3 <- st$tip.label[current_combo[[4]]]

  # define taxa relationships
  conc_taxa <- c(p2,p3) # taxa that are concordant with the species tree when together in the quartet
  test_taxa <- c(p1,p2) # taxa we are testing for introgression

  # subset the quartet table to only the row corresponding to the current quartet of taxa
  # qdf <- as.data.frame(quartet_table)
  # sub_qt <- quartet_df[quartet_df[,p1]==1&quartet_df[,p2]==1&quartet_df[,p3]==1&quartet_df[,o]==1,]
  sub_qt <- quartet_table[quartet_table[,p1]==1&quartet_table[,p2]==1&quartet_table[,p3]==1&quartet_table[,o]==1,]

  # create list of taxa names in the order in which they appear in the quartet table
  # qt_taxa <- colnames(sub_qt)[which(sub_qt == 1, arr.ind=T)[, "col"]]
  qt_taxa <- names(sub_qt)[sub_qt==1]
  v<-c("12|34","13|24","14|23") # vector of column names to exclude when determining new taxon's position (incase the value of one of these columns is 1)
  qt_taxa <- setdiff(qt_taxa,v)

  # set column definitions
  conc_column <- find_column(conc_taxa,qt_taxa)
  abba_column <- find_column(test_taxa,qt_taxa)
  baba_column <- setdiff(c("12|34","13|24","14|23"),c(abba_column,conc_column))

  # get values from columns
  abba <- sub_qt[abba_column]
  baba <- sub_qt[baba_column]
  conc <- sub_qt[conc_column]
  total <- abba+baba+conc
  d_stat <- (abba-baba)/(abba+baba)
  f_stat <- (abba-baba)/(conc-baba)

  # create vector of data to return
  add2df <- c(o,p1,p2,p3,abba,baba,conc,total,d_stat,f_stat)
  add2df
}

quaint <- function(test_taxa,species_tree,quartet_table,qt_vector,outgroup=NULL) {

  mrca <- getMRCA(species_tree, test_taxa)
  node_groups <- nodeGroups(species_tree,mrca) # create node groups to determine where taxa are
  t1 <- test_taxa[1]
  t2 <- test_taxa[2]
  # get tip number from tip label
  t1_num <- which(species_tree$tip.label==t1)
  t2_num <- which(species_tree$tip.label==t2)
  # get node group
  t1_group <- which(sapply(node_groups, is.element, el = t1_num))
  t2_group <- which(sapply(node_groups, is.element, el = t2_num))

  # get tip numbers for all sister taxa
  sister_tips_t1 <- setdiff(node_groups[[t1_group]],c(t1_num))
  sister_tips_t2 <- setdiff(node_groups[[t2_group]],c(t2_num))

  # check if taxa are sister
  if (length(sister_tips_t1) == 0 & length(sister_tips_t2) == 0) {
    stop("Test taxa are sister taxa and cannot be tested for introgression.")
  }

  # get outgroup tips from tree
  og_tips_from_tree <- node_groups[[setdiff(c(1,2,3),c(t1_group,t2_group))]] #outgroup tips from comparing sister tips

  # if outgroup is defined, use that as the outgroup tips. otherwise, use outgroup tips from tree
  if(hasArg(outgroup)) {
    og_tips <- c()
    for(i in 1:length(outgroup)){
      og_tips[i] <- which(species_tree$tip.label==outgroup[i])
    }
    # make sure given outgroup is compatible with the test taxa
    if (length(intersect(og_tips, sister_tips_t1)) > 0 || length(intersect(og_tips, sister_tips_t2)) > 0) {
      stop(paste(c("The outgroup taxa are incompatible with testing for introgression with ",t1," and ",t2,". Make sure the outgroup is defined correctly.")))
    }
  } else {
    og_tips <- og_tips_from_tree
  }

  # get all quartets that can be used to test for introgression
  combos_1 <- expand.grid(og_tips,t1_num,t2_num,sister_tips_t2)
  combos_2 <- expand.grid(og_tips,t2_num,t1_num,sister_tips_t1)
  all_combos <- rbind(combos_1,combos_2)

  ##### subset quartet table into only those with quartets in all_combos
  # make all_combos into a vector of tip names
  ac <- all_combos
  ac[] = sapply(ac, function(x){ x[x] <- species_tree$tip.label[x]})
  ac_sorted_vector <- apply(ac, 1, function(x) paste(sort(x),collapse = ","))

  # if ran via quaint_all function, use the pre-subsetted quartet table (to save time)
  if(!hasArg(qt_vector)) {qt_vector <- get_qt_vector(quartet_table)}

  # further subset the quartet table to only the quartets needed for this test (to save time)
  subset_qt <- quartet_table[which(qt_vector %in% ac_sorted_vector),]

  # parse quartet frequency df and prepare output data frame
  num_rows <- 1:nrow(all_combos) # get number of rows
  data_list <- lapply(num_rows,parse_table,combos=all_combos,quartet_table=subset_qt,st=species_tree)
  return_df <- as.data.frame(do.call(rbind, data_list)) # convert output from lapply (list) into data frame
  cnames<-c(  "outgroup","taxon1","taxon2","taxon3",
              "abba","baba","concordant","total","d","f") # appropriate column names
  colnames(return_df) <- cnames # add column names to df
  return_df <- return_df %>% mutate_at(c("abba","baba","concordant","total","d","f"), as.numeric) # convert columns to numeric as approprtiate

  #perform chi squared tests
  return_df <- return_df %>%
    rowwise() %>%
    mutate(
      test_stat = ifelse(abba>1 & baba>1,suppressWarnings(chisq.test(c(abba,baba))$statistic),NA),
      p_val = ifelse(abba>1 & baba>1,suppressWarnings(chisq.test(c(abba,baba))$p.value),NA) # perform chi squared test on all rows, including total row
    )


  return_df
}

quaint_all <- function(species_tree,quartet_table,outgroup) {

  # make list of pairs of species that are can be tested w/ quartet test
  test_taxa_list <- get_valid_test_pairs(species_tree,outgroup)

  # create quartet vector for faster lookups
  qt_vector <- get_qt_vector(quartet_table)

  # run quaint on all taxon pairs
  pos <- 1:ncol(test_taxa_list)
  taxon_pair_df <- lapply(pos,function(x) quaint(test_taxa_list[,x],species_tree,quartet_table,qt_vector,outgroup = outgroup))
  taxon_pair_df <- as.data.frame(do.call(rbind, taxon_pair_df))

  taxon_pair_df$adj_p_val <- p.adjust(taxon_pair_df$p_val, method = "fdr")

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

get_qt_vector <- function(quartet_table) {
  # convert quartet table into vector of taxon names for faster processing
  drop<-c("12|34","13|24","14|23","1234") # drop these columns
  qt_dropped_cnames <- quartet_table[,!colnames(quartet_table) %in% drop] # remove columns in "drop"
  sorted_colnames <- sort(colnames(qt_dropped_cnames)) # get list of sorted column names (for more efficient comparisons later)
  qt_dropped_cnames <- qt_dropped_cnames[, sorted_colnames] # sort column names in matrix
  indices <- which(qt_dropped_cnames == 1, arr.ind = TRUE)
  qt_vector <- split(colnames(qt_dropped_cnames)[indices[, 2]], indices[, 1])
  qt_vector <- lapply(qt_vector,paste,collapse=",")
  qt_vector
}

summarize_quaint_table <- function(quaint_table,alpha = 0.05,use_adjusted_p=TRUE) {
  p_col <- if (use_adjusted_p) "adj_p_val" else "p_val"

  # add column with test pair as a factor
  quaint_table$pair <- factor(
    apply(
      quaint_table[ , c("taxon1", "taxon2")],          # take the two columns
      1,                                               # row by row
      function(x) paste(sort(x), collapse = "|")       # sort & paste
    )
  )
  
  # pull out the p‑values column
  pvals <- quaint_table[[p_col]]
  
  # add the significance columns
  quaint_table$d_sig <- ifelse(pvals < alpha, quaint_table$d, 0)
  quaint_table$f_sig <- ifelse(pvals < alpha, quaint_table$f, 0)
  
  # summarize the table by test pair
  summary_table <- quaint_table %>%
    group_by(pair) %>%
    summarise(
      total_abba = sum(abba),
      total_baba = sum(baba),
      total_concordant = sum(concordant),
      n_tests = n(),
      n_positive = sum(d > 0),
      n_positive_significant = sum(d > 0 & get(p_col) < alpha & !is.na(get(p_col)), na.rm = TRUE),
      # mean_d = ifelse(sum(d_sig >= 0) > 0, mean(d_sig[d_sig >= 0]), 0),
      mean_d = ifelse(sum(d_sig,na.rm = TRUE) > 0, mean(d_sig,na.rm=TRUE), 0),
      mean_f = mean(f_sig[d_sig>=0],na.rm=TRUE)
    ) %>%
    mutate(
      proportion = n_positive_significant/n_tests,
      d = (total_abba-total_baba)/(total_abba+total_baba)
    ) %>%
    separate(pair, into = c("taxon1", "taxon2"), sep = "\\|")

  summary_table
}
