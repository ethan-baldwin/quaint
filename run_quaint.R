library("MSCquartets")
library("combinat")
library("dplyr")
library("stringr")

# This is pulling in the functions from "quaint_functions.R". Replace with path to this  
source("path/to/quaint/quaint_functions.R")
source("~/GitHub/sarracenia_genome_assembly/phylo/quartets/quartets_functions.R")

setwd("~/Research/genome_assembly/introgression/quartets/")

############### Read tree files and set outgroups. #################
# Read gene trees. This should be one file with all of the gene trees in newick format (one tree per line).
gt <- read.tree("merged.tre") 

# Read species tree. This should be the ROOTED species tree in newick format.
st <- read.tree("astral3.tre")

# Vector of outgroups names in your trees.
og <- c("heliamphora_pulchella_SRR25244091","heliamphora_ciliata_SRR24877724","darlingtonia_californica_SRR24877818") 

# In case your species tree isn't rooted.
# st <- root(st,og)

################### Calculate quartet table ##################

qt <- quartetTable(gt,taxonnames = st$tip.label) # calculate quartet frequencies

############### Testing for introgression between two taxa ################  

# vector containing names of taxa pair you want to test introgression for
test_taxon_names <- c("S_jonesii","S_purpurea_montana") 

# command to perform quartet frequency test for introgression
d_table <- dstat_table(test_taxa = test_taxon_names,species_tree = st, quartet_table = qt) 

# write results to csv
write.csv(d_table, file = "taxon1_taxon2_table.csv")

############### Testing for introgression between all taxon pairs ################  

d_all <- dstat_table_all(species_tree=st,outgroup=og,quartet_table=qt,alpha = 0.05)

# grab only the summary rows, add proportion of tests passed.
d_sum <- d_all[d_all$outgroup=="Total",]
d_sum$proportion <- d_sum$p_val/d_sum$d

# write results to csv
write.csv(d_sum, file = "quartet_test_summary.csv")




