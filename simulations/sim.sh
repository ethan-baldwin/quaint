
ml ms
ml Seq-Gen

ms 15 500 -T -I 7 3 2 3 2 3 1 1 -m 1 2 2.5 -m 2 1 2.5 -m 2 3 2.5 -m 3 2 2.5 -m 4 5 2.5 -m 5 4 2.5 -m 5 6 2.5 -m 6 5 2.5 -em 2.0 3 4 2.5 -em 2.0 4 3 2.5 -m 7 6 0.1 | grep -v // > gene_trees_m.tre
seq-gen -l 1000 -s 0.1 -m GTR gene_trees_m.tre > seqs_m.phy
