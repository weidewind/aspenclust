#!/usr/bin/bash

for protein in h1 n1 h3 n2 h1pandemic n1pandemic
do
perl test_alleles.pl --input noeggs --stattype mean --protein $protein --simnumber 200 --norm countnorm --bigtag test200 &
done
wait
