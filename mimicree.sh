#!/usr/bin/bash

tag=mimicree_start_00x20 #mimicree_non_recurrent_pos_epistasis
protein=mixed
datafolder=/export/home/popova/workspace/aspenclust/data/${tag}
outfolder=/export/home/popova/workspace/aspenclust/output/${tag}/disttest/nsyn
mkdir -p $outfolder

perl test_distances.pl --protein ${protein} --state nsyn --bigtag disttest --input ${tag} >${outfolder}/disttest_mtub.log 2>${outfolder}/disttest_mtub.errlog

perl fdr_distances.pl --protein ${protein} --state nsyn --input ${tag} --bigtag disttest  \
--pairsfile ${datafolder}/${protein}.site_pairs >${outfolder}/fdrdist.out 2>${outfolder}/fdrdist.err

perl grep_distances_fdr.pl --input ${tag} --protein ${protein} --state nsyn --bigtag disttest \
 --sites --nopdb --pairsfile ${datafolder}/${protein}.site_pairs