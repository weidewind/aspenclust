#!/usr/bin/bash

infolder=mtub_no_domains_210423_epistat
perl grep_pheno_distances_fdr_stepwise.pl --input $infolder --protein 9drugs --state nsyn --bigtag disttest \
 --sites --nopdb --pairsfile /export/home/popova/workspace/aspenclust/data/${infolder}/9drugs.site_pairs --numfakes 500 --numchunks 10