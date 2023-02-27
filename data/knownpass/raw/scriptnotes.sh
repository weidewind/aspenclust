for p in h1 n1 h1pdm n1pdm h3 n2
# concatenates parts downloaded from gisaid
do
ls ${p}_* | xargs cat >${p}.longstrains.fa
perl -pi -e 's/\n(^[^>].*=.*)/$1/ge' ${p}.longstrains.fa ## fix all entries where names are split into several lines
done


## some of the sequences are shorter than the minlength used at gisaid for filtering. remove them
parallel --link python ../../../../flutils/select_noshorter_than.py --infile {1}.longstrains.fa --outfile {1}.minlength.validlongstrains.fa --minlength {2} ::: h1 n1 h1pdm n1pdm h3 n2 ::: 1695 1410 1698 1407 1698 1407

## generate short names. Make IDs from xENTIFIER or from EPI_ISL_x ? Now use xENTIFIER
parallel python ../../../../flutils/parse_passages_from_gisaid.py --idtype entifier --infile {1}.minlength.validlongstrains.fa --outfile {1}.passages --unparsedfile {1}.unparsed_pass ::: h1 n1 h1pdm n1pdm h3 n2
parallel python ../../../../flutils/cut_fasta_names_from_gisaid.py --idtype entifier --infile {1}.minlength.validlongstrains.fa --outfile {1}.strains.fa --unparsedfile {1}.unparsed_id ::: h1 n1 h1pdm n1pdm h3 n2
# Done: to get entifier, you must directly target entifier)) because there could be different number of elements divided by = sign


## find all sequences with the same name and check if they are identical
## with entifier as id there is no duplicates
parallel perl ../../../../flutils/find_and_resolve_duplicates.pl --infile {1}.strains.fa --print --outfile {1} ::: h1 n1 h1pdm n1pdm h3 n2
parallel grep -c "different" {1}.duplicates ::: h1 n1 h1pdm n1pdm h3 n2


## to many nonperfect sequences, won't remove them all
#parallel perl ../../../../flutils/remove_nonperfect_seqs.pl --infile {1}.strains.fa --outfile {1}.clean.strains.fa ::: h1 n1 h1pdm n1pdm h3 n2

# remove sequences with more than 2 ambiguous chars
parallel perl ../../../../flutils/remove_nonperfect_seqs.pl --infile {1}.strains.fa --outfile {1}.max2ambichar.strains.fa --max 2 --countfile {1}.count.ambiguous ::: h1 n1 h1pdm n1pdm h3 n2

# TODO remove seasonal h1n1 collected after 2010 (?) I saw one!

# TODO maybe find all strains that are not in previous alignments and add them to the said alignments? (at least add old seeds for new alignments) (at least check how many new strings do we get.)
# old data - noeggs. in h1n1pdm, epi_isl was used as identifier. in h3n2 - god knows what, i cannot find such substrings in long-named fasta
# TODO also we definitely do not have any new h1n1 seasonal strings, so just take old data and  count strings with known passaging, give that number and do not talk about h1n1 anymore. 

# or --thread -1 if you do not use parallel, to auto detect number of procs 
parallel "mafft --thread 15 --retree 1 --reorder {1}.max2ambichar.strains.fa >{1}.max2ambichar.strains.reordered.aligned" ::: h1 n1 h1pdm n1pdm h3 n2
# or
java -jar ~/local/macse_v2.05.jar -prog alignSequences -seq h1.strains.fa

#starts and ends were found by hand, using ksu sequences as reference, and hard-coded into trimmer.py
# todo maybe just add littleksu seq here and trim automatically)
#trim ends and delete gaps and realign
python trimmer.py
parallel perl -pi -e 's/-//g' {1}.max2ambichar.strains.reordered.aligned.trimmed ::: h1 n1 h1pdm n1pdm h3 n2
parallel "head -2 ../../little_ksu/{1}.ancestor.fa >{1}.max2ambichar.trimmed.withksu" ::: h1 n1 h1pdm n1pdm h3 n2
parallel "cat {1}.max2ambichar.strains.reordered.aligned.trimmed >>{1}.max2ambichar.trimmed.withksu" ::: h1 n1 h1pdm n1pdm h3 n2
# without --reorder to leave ksu strain first
parallel "mafft --thread 15 --retree 1 {1}.max2ambichar.trimmed.withksu >{1}.max2ambichar.trimmed.realigned.withksu" ::: h1 n1 h1pdm n1pdm h3 n2
# # Select seqs with the same pattern of gaps as in ksu. If the only thing that differs is missing ATG start, just add it
parallel "python select_equally_gapped.py --input {1}.max2ambichar.trimmed.realigned.withksu --output {1}.max2ambichar.trimmed.asksu2 --leftovers {1}.max2ambichar.trimmed.asksu.leftovers2" ::: h1 n1 h1pdm n1pdm h3 n2
parallel perl -pi -e 's/-//g' {1}.max2ambichar.trimmed.asksu2 ::: h1 n1 h1pdm n1pdm h3 n2

#select unique sequences (sequence+passaging)
parallel cd-hit -i {1}.max2ambichar.trimmed.asksu2 -o {1}.max2ambichar.asksu2.clustered -c 1 ::: h1 n1 h1pdm n1pdm h3 n2
parallel python select_nonredundant.py --input {1}.max2ambichar.trimmed.asksu2 --passages {1}.passages --clusters {1}.max2ambichar.asksu2.clustered.clstr --output {1}.uniq ::: h1 n1 h1pdm n1pdm h3 n2

#add one seq from littleksu alignment and realign
#then choose those aligned as reference littleksu (same gap indices)
python select_equally_gapped.py --input n2.max2ambichar.trimmed.realigned --output n2.max2ambichar.trimmed.asksu --leftovers n2.max2ambichar.trimmed.asksu.leftovers

# build trees
# automatically prints alignment of unique sequences in human-readable form!
parallel iqtree -s {1}.uniq -m GTR+I+G -nt AUTO -st CODON ::: h1 n1 h1pdm n1pdm h3 n2
# just trying to guess the best model:
iqtree -s ../h1.uniq -m MFP -nt AUTO -st CODON
# for seqs with ambigious characters, 
screen -U -S n1dis python fix_ambiguities.py --tree n1.uniq.treefile --phylip n1.uniq.uniqueseq.phy --output n1.disambig.fast
#some bullshit happens because at the end iqtree adds sequnces in h1, but not in the other proteins (because in all other cases this is not the final tree! but good enough for our purposes here), thus:
python fix_ambiguities.py --tree h1.uniq.treefile --fasta h1.uniq --output h1.disambig.fast
# since iqtree keeps maximum two strains with identical sequences, now we need to add the rest of the strains: read {1}.uniq.log, select missing strains and their non-missing twins, copy their seqs to the alignment
python addMissingTwins.py --fasta h1pdm.disambig.fast --iqlog h1pdm.uniq.log --output h1pdm.disambig.full
# get rid of strains with stop codons
parallel python deleteStopped.py --fasta {1}.disambig.full --output {1}.disambig.full.nostop ::: h1 n1 h3 n2 h1pdm n1pdm
# Finally! select strains with unambigious passaging info
parallel python selectKnown.py --passages {1}.passages --fasta {1}.disambig.full.nostop --output {1}.known.fa ::: h1 n1 h3 n2 h1pdm n1pdm
# add tags, so we can cluster sequences again and retain sequences with the same seq but different passaging history 
parallel python addTag.py --input {1}.known.fa --output {1}.known.withtag --passages {1}.passages ::: h3 n2 n1 h1 h1pdm n1pdm
#cluster:
cd-hit -i n2.known.withtag -o n2.known.clustered.withtag -c 0.996 # 3222
cd-hit -i h3.known.withtag -o h3.known.clustered.withtag -c 0.9955 # 2711
cd-hit -i n1pdm.known.withtag -o n1pdm.known.clustered.withtag -c 0.9965 # 2787
cd-hit -i h1pdm.known.withtag -o h1pdm.known.clustered.withtag -c 0.9965 # 2993

# delete tags
parallel cp {1}.known.fa clustered/{1}.readystrains.fa ::: h1 n1
parallel python cutTag.py --input {1}.known.clustered.withtag --output clustered/{1}.readystrains.fa ::: h3 n2 h1pdm n1pdm

# auto detected best codon model for h1 - MGK+F1X4+R5: 
#MGK - Nonsynonymous/synonymous (dn/ds) rate ratio with additional transition/transversion (ts/tv) rate ratio
#+F1X4 Unequal nucleotide frequencies but equal nt frequencies over three codon positions.
#+R5 - 5 categories of sites with different rates, no distribution is assumed:
# FreeRate model (Yang, 1995; Soubrier et al., 2012) that generalizes
# the +G model by relaxing the assumption of Gamma-distributed rates.
#+I Probably also should be added
# -asr - ancestral reconstruction -te -user-defined tree for asr
screen -U -S h1pdmgtr iqtree -s h1pdm.known.fa -m GTR+G+I -nt AUTO -asr

# find oldest strain and reroot tree, using it as an outgroup
parallel "python findOldest.py --datefile {1}.longstrains.fa --fasta clustered/{1}.readystrains.fa | xargs -I STRAIN python ../../../../ksuAncestorPipeline/reroot_tree.py --tree clustered/{1}.readystrains.fa.treefile --root STRAIN --outfile clustered/{1}.treefile.rerooted" ::: n1 h1 n1pdm h1pdm n2 h3
# convert char states for internal nodes to fasta
parallel python stateToFasta.py --input clustered/{1}.readystrains.fa.state --output clustered/{1}.internal.fa ::: h1 n1 h1pdm n1pdm n2 h3

# passages for our strains in fasta-format
parallel python passagesToFasta.py --fasta clustered/{1}.readystrains.fa --passages {1}.passages --output clustered/{1}.passages.fa ::: h1 n1 h1pdm n1pdm n2 h3

#recontsruct ancestral states by mega
./reconstrust_passage_states.sh

#parse mega output
parallel perl ../../../../flutils/parse_mega_ancestors_output.pl --protein {1} --statefile {1}.passages.asletters.ancestral_states.txt --nodemapfile {1}.passages.asletters.nodeMap.txt --tree_file {1}.treefile.rerooted --output {1}.passages.ancstates2 --folder clustered/ ::: h1 n1 h1pdm n1pdm h3 n2

# todo ... rename and cp splits..
parallel cp clustered/{1}.passages.ancstates2 ../megasplits/{1}.splits ::: h1 n1 h1pdm n1pdm n2 h3

parallel cp clustered/{1}.treefile.rerooted ../{1}.l.r.newick ::: h1 n1 h1pdm n1pdm n2 h3
parallel perl -pi -e 's/A/a/g' clustered/{1}.internal.fa ::: h1 n1 h1pdm n1pdm n2 h3
parallel perl -pi -e 's/C/c/g' clustered/{1}.internal.fa ::: h1 n1 h1pdm n1pdm n2 h3
parallel perl -pi -e 's/T/t/g' clustered/{1}.internal.fa ::: h1 n1 h1pdm n1pdm n2 h3
parallel perl -pi -e 's/G/g/g' clustered/{1}.internal.fa ::: h1 n1 h1pdm n1pdm n2 h3

parallel "cp clustered/{1}.readystrains.fa ../{1}.nointernals.fa" ::: h1 n1 h1pdm n1pdm n2 h3
parallel "cp clustered/{1}.internal.fa ../{1}.ancestor.fa" ::: h1 n1 h1pdm n1pdm n2 h3


#parallel  iqtree -s {1}.passages.fa -st MORPH -te {1}.treefile.rerooted -asr  ::: h1 n1 h1pdm n1pdm n2 h3
#parallel iqtree -s {1}.max2ambichar.trimmed.realigned -m MFP -nt 16 ::: h1 n1 h1pdm n1pdm h3 n2
#to find which sequences are clustered with the strange ones visible in iqtree alignment
#parallel cd-hit -i h3.max2ambichar.strains.fa -o {1}.max2ambichar.strains.clustered1.fa -c 1 ::: h1 n1 h1pdm n1pdm h3 n2
# strange sequences are found by hand in iqtree-produced .phy files, and checked for identical twins using cd-hit output. Aaand also hard-coded into trimmer.py

#todo root the tree!!


