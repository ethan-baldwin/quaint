library("MSCquartets")
library("combinat")
library("dplyr")
library("stringr")
# source("/Users/ethan/GitHub/sarracenia_genome_assembly/phylo/quartets/quartets_functions.R")
source("/Users/ethan/OneDrive - University of Georgia/Documents/Github/sarracenia_genome_assembly/phylo/quartets/quartets_functions.R")


setwd("~/Research/genome_assembly/introgression/quartets/")

# read data and get quartet frequency table
gt <- read.tree("merged.tre") # read gene tree
st <- read.tree("astral3.tre") # read species tree
og <- c("heliamphora_pulchella_SRR25244091","heliamphora_ciliata_SRR24877724","darlingtonia_californica_SRR24877818") # define outgroup
rooted_st <- root(st,outgroup = og) # root species tree
rooted_st <- drop.tip(rooted_st, "darlingtonia_californica_SRR24877818") # drop darlingtonia
all_t_names<-rooted_st$tip.label # get list of taxa names to calculate quartet frequencies
# all_t_names<-setdiff(st$tip.label, c("darlingtonia_californica_SRR24877818")) # remove
test_taxon_names <- c("S_jonesii","S_purpurea_montana") # define taxa pair you want to test introgression for
allpairs <- get_valid_test_pairs(rooted_st,og)
# rooted_st_2 <- drop.tip(rooted_st,"darlingtonia_californica_SRR24877818")

qt <- quartetTable(gt,taxonnames = all_t_names) # calculate quartet frequencies
# qt2 <- quartetTable(gt,taxonnames = all_t_names_2) # calculate quartet frequencies


test_taxon_names <- c("S_minor_2","S_purpurea_montana") # define taxa pair you want to test introgression for
d_table <- dstat_table(test_taxa = test_taxon_names,species_tree = rooted_st, quartet_table = qt)

d_all <- dstat_table_all(allpairs,rooted_st,qt,alpha = 0.01)
d_sum <- d_all[d_all$outgroup=="Total",]
d_sum$proportion <- d_sum$p_val/d_sum$d

# write.csv(d_all, file = "all_quartet_test_summary_row.csv")
write.csv(d_sum, file = "quartet_test_summary_01.csv")


write.csv(d_table, file = "montana_jonesii_d_table2.csv")

###### test on simulated tree ########

simtree1 <- read.tree("simtrees.1.newick")
simtree2 <- read.tree("simtrees.2.newick")
simtree3 <- read.tree("simtrees.3.newick")
simtree4 <- read.tree("simtrees.4.newick")
simtree3000 <- read.tree("simtrees.3000.newick")


remove_suffix <- function(tree) {
  tree$tip.label <- str_remove(tree$tip.label,'[_]1')
  return(tree)
}

simtree1 <- lapply(simtree1, remove_suffix)
class(simtree1)<-"multiPhylo"
simtree2 <- lapply(simtree2, remove_suffix)
class(simtree2)<-"multiPhylo"
simtree3 <- lapply(simtree3, remove_suffix)
class(simtree3)<-"multiPhylo"
simtree4 <- lapply(simtree4, remove_suffix)
class(simtree4)<-"multiPhylo"
simtree3000 <- lapply(simtree3000, remove_suffix)
class(simtree3000)<-"multiPhylo"

# write.tree(gt_renamed, file = "simtrees_renamed.tre")

qt1 <- quartetTable(simtree1,taxonnames = all_t_names) # calculate quartet frequencies
qt2 <- quartetTable(simtree2,taxonnames = all_t_names) # calculate quartet frequencies
qt3 <- quartetTable(simtree3,taxonnames = all_t_names) # calculate quartet frequencies
qt4 <- quartetTable(simtree4,taxonnames = all_t_names) # calculate quartet frequencies
qt3000 <- quartetTable(simtree3000,taxonnames = all_t_names) # calculate quartet frequencies


d_all_1 <- dstat_table_all(allpairs,rooted_st,qt1)
d_all_2 <- dstat_table_all(allpairs,rooted_st,qt2)
d_all_3 <- dstat_table_all(allpairs,rooted_st,qt3)
d_all_4 <- dstat_table_all(allpairs,rooted_st,qt4)
d_all_3000 <- dstat_table_all(allpairs,rooted_st,qt3000)

write.csv(d_all_1, file = "sim1_quartet_test.csv")
write.csv(d_all_2, file = "sim2_quartet_test.csv")
write.csv(d_all_3, file = "sim3_quartet_test.csv")
write.csv(d_all_4, file = "sim4_quartet_test.csv")
write.csv(d_all_3000, file = "sim3000_quartet_test.csv")

# run table parser to get d-stat table
test_taxon_names <- c("S_alabamensis","S_psittacina_P008") # define taxa pair you want to test introgression for
test_taxon_names <- c("S_alabamensis","S_rubra_RVCAGA001") # define taxa pair you want to test introgression for
test_taxon_names <- c("S_alabamensis","S_leucophylla") # define taxa pair you want to test introgression for
test_taxon_names <- c("S_alabamensis","S_rubra_gulfensis") # define taxa pair you want to test introgression for


d_table1 <- dstat_table(test_taxon_names,rooted_st,qt1)
d_table2 <- dstat_table(test_taxon_names,rooted_st,qt2)
d_table3 <- dstat_table(test_taxon_names,rooted_st,qt3)
d_table4 <- dstat_table(test_taxon_names,rooted_st,qt4)
