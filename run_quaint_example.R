############### Read tree files and set outgroups. #################
# Read gene trees. This should be one file with all of the gene trees in newick format (one tree per line).
gene_trees <- read.tree("merged.tre") 

# Read species tree. This should be the ROOTED species tree in newick format.
sp_tree <- read.tree("astral3_rooted.tre")

# Vector of outgroup names in your trees.
outgroup <- c("heliamphora_pulchella_SRR25244091","heliamphora_ciliata_SRR24877724") 

# In case your species tree isn't rooted.
# st <- root(st,og)

################### Calculate quartet table ##################

quartet_table <- quartetTable(gene_trees,taxonnames = sp_tree$tip.label) # calculate quartet frequencies

############### Testing for introgression between two taxa ################  

# vector containing names of taxa pair you want to test introgression for
test_taxon_names <- c("S_jonesii","S_purpurea_montana") 

# command to perform quartet frequency test for introgression
quaint_table <- quaint(test_taxon_names,sp_tree,quartet_table) 

# write results to csv
write.csv(quaint_table, file = "S_jonesii_S_purpurea_montana_table.csv")

############### Testing for introgression between all taxon pairs ################  

qt_all_pairs <- quaint_all(sp_tree,quartet_table,outgroup)

# write full table to csv
write.csv(qt_all_pairs, file = "quaint_all_pairs.csv")

# summarize results by pair, using a p value cutoff of 0.05 for the chi sq tests
sum_table <- summarize_quaint_table(qt_all_pairs,alpha = 0.05)

write.csv(sum_table,file = "quaint_summary.csv")


