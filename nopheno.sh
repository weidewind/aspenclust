#!/usr/bin/bash

tag=mtub_no_domains_210423_epistat_nopheno
datafolder=/export/home/popova/workspace/aspenclust/data/${tag}
outfolder=/export/home/popova/workspace/aspenclust/output/${tag}/disttest/nsyn
mkdir -p $outfolder

perl test_distances.pl --protein 9drugs --state nsyn --bigtag disttest --input ${tag} >${outfolder}/disttest_mtub.log 2>${outfolder}/disttest_mtub.errlog

perl fdr_distances.pl --protein 9drugs --state nsyn --input ${tag} --bigtag disttest  \
--pairsfile ${datafolder}/9drugs.site_pairs >${outfolder}/fdrdist_mtub.out 2>${outfolder}/fdrdist_mtub.err

perl grep_distances_fdr.pl --input ${tag} --protein 9drugs --state nsyn --bigtag disttest \
 --sites --nopdb --pairsfile ${datafolder}/9drugs.site_pairs