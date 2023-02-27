#!/usr/bin/bash
folder_iqtree=/export/home/popova/workspace/aspenclust/output/covid/disttest/nsyn
folder_usher=/export/home/popova/workspace/aspenclust/output/covid/disttest_usher/nsyn
cfile=concordant_article_suppl
dfile=discordant_article_suppl
output=/export/home/popova/workspace/aspenclust/output/covid/disttest_iqtree_and_usher_article_suppl_merged_withBH.csv

head -1 ${folder_iqtree}/S_dist.fdr_results >${folder_iqtree}/${cfile}_with_disttest.csv
head -1 ${folder_iqtree}/S_dist.fdr_results >${folder_iqtree}/${dfile}_with_disttest.csv
cat ${folder_iqtree}/${cfile}.csv | xargs -I{} grep ^{}" " ${folder_iqtree}/S_dist.fdr_results >>${folder_iqtree}/${cfile}_with_disttest.csv
cat ${folder_iqtree}/${dfile}.csv | xargs -I{} grep ^{}" " ${folder_iqtree}/S_dist.fdr_results >>${folder_iqtree}/${dfile}_with_disttest.csv

head -1 ${folder_usher}/S_usher_dist.fdr_results >${folder_usher}/${cfile}_with_disttest.csv
head -1 ${folder_usher}/S_usher_dist.fdr_results >${folder_usher}/${dfile}_with_disttest.csv
cat ${folder_usher}/${cfile}.csv | xargs -I{} grep ^{}" " ${folder_usher}/S_usher_dist.fdr_results >>${folder_usher}/${cfile}_with_disttest.csv
cat ${folder_usher}/${dfile}.csv | xargs -I{} grep ^{}" " ${folder_usher}/S_usher_dist.fdr_results >>${folder_usher}/${dfile}_with_disttest.csv


#Rscript merge.R --zero 0.0025 --folder_usher ${folder_usher} --folder_iqtree ${folder_iqtree} --cfile ${cfile}_with_disttest.csv --dfile ${dfile}_with_disttest.csv --output $output

Rscript BH.R --zero 0.0025 --folder ${folder_usher} --file ${cfile}_with_disttest.csv --output ${folder_usher}/usher_${cfile}_withBH.csv
Rscript BH.R --zero 0.0025 --folder ${folder_usher} --file ${dfile}_with_disttest.csv --output ${folder_usher}/usher_${dfile}_withBH.csv
Rscript BH.R --zero 0.0025 --folder ${folder_iqtree} --file ${cfile}_with_disttest.csv --output ${folder_iqtree}/iqtree_${cfile}_withBH.csv
Rscript BH.R --zero 0.0025 --folder ${folder_iqtree} --file ${dfile}_with_disttest.csv --output ${folder_iqtree}/iqtree_${dfile}_withBH.csv