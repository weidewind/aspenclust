#!/usr/bin/bash

folder="/export/home/popova/workspace/aspenclust/data/little_ksu"
for protein in n1pandemic
do
	infile=$folder/passage_history/${protein}_passages_quaternary
	fafile=${infile}.fa
	cp $infile $fafile
	perl -pi -e 's/STRAIN/>STRAIN/g' $fafile
	perl -pi -e 's/\t/\n/g' $fafile
	perl -pi -e 's/^1/c/g' $fafile
	perl -pi -e 's/^0/a/g' $fafile
	perl -pi -e 's/^2/g/g' $fafile
	perl -pi -e 's/^3/t/g' $fafile

	xvfb-run megacc -a $folder/passage_history/ancestral_states/ancestral_seqs_ML_nucleotide_states.mao \
-d $fafile -f Fasta \
-t $folder/${protein}.l.r.newick -o $folder/passage_history/${protein}.ancstates
done
