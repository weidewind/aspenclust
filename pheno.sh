#!/usr/bin/bash

perl test_distances.pl --protein 9drugs --state nsyn --bigtag disttest --input mtub >disttest_mtub.log 2>disttest_mtub.errlog

for ph in 1 2 4 6 8 9 12 13 14
do
perl fdr_distances.pl --protein 9drugs --state nsyn --input mtub --bigtag disttest --phenotype $ph \
--pairsfile /export/home/popova/workspace/aspenclust/data/mtub/9drugs.site_pairs >fdrdist_mtub_${ph}.out 2>fdrdist_mtub_${ph}.err
done

 perl grep_pheno_distances_fdr_stepwise.pl --input mtub_test --protein 9drugs --state nsyn --bigtag disttest_stepwise \
 --sites --nopdb --pairsfile /export/home/popova/workspace/aspenclust/data/mtub/9drugs.site_pairs --numfakes 500 --numchunks 10