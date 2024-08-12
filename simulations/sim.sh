
ml ms
ml Seq-Gen
ml ASTER/1.16-GCC-11.3.0

ms 15 500 -T -I 6 3 2 3 2 3 2 \
-ej 1 2 1 -ej 1 4 3 -ej 2 5 1 -ej 3 3 1 -ej 4 6 1 \
-em 0.5 1 5 0.5 | \
tail -n +4 | grep -v // |  grep -v ^$ > gene_trees_m.tre

# while read LINE; do echo $LINE > trees/${COUNTER}.tre; COUNTER=$((COUNTER + 1)); done < gene_trees_m_no_blanks.tre

COUNTER=1
while read LINE; do
  echo $LINE > trees/${COUNTER}.tre
  seq-gen -l 1000 -s 0.1 -m GTR trees/${COUNTER}.tre > seqs/seqs_${COUNTER}.phy
  # cd est_trees
  # iqtree2 -s ../seqs/seqs_${COUNTER}.phy -m MFP -B 1000 --prefix $COUNTER
  # cd ..
  COUNTER=$((COUNTER + 1))
done < gene_trees_m.tre

cat est_trees/*treefile > merged_est.tre
cat trees/*tre > merged_sim.tre

astral -i merged_sim.tre -o astral.sim.tre --root 15
astral -i merged_est.tre -o astral.est.tre --root 15


COUNTER=1
while read LINE; do
  echo $LINE > trees/${COUNTER}.tre
  seq-gen -l 1000 -s 0.1 -m GTR trees/${COUNTER}.tre > seqs/seqs_${COUNTER}.phy
  cd est_trees
  iqtree2 -s ../seqs/seqs_${COUNTER}.phy -m MFP -B 1000 --prefix $COUNTER
  cd ..
  COUNTER=$((COUNTER + 1))
done < gene_trees_m.tre

# seq-gen -l 1000 -s 0.1 -m GTR gene_trees_m.tre > seqs_m.phy

# -m 1 2 2.5 -m 2 1 2.5 -m 2 3 2.5 \
# -m 3 2 2.5 -m 4 5 2.5 -m 5 4 2.5 -m 5 6 2.5 -m 6 5 2.5 \
# -em 2.0 3 4 2.5 -em 2.0 4 3 2.5 | \

# -ma x 0 0 0 1 0 \
# x 0 0 0 0 0 \
# x 0 0 0 0 0 \
# x 0 0 0 0 0 \
# x 1 0 0 0 0 \
# x 0 0 0 0 0 | \
