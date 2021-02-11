#!/bin/bash
while getopts t:n:i: option
do
                case "${option}" in
                        t) tag=${OPTARG};;
                        n) number=${OPTARG};;
                        i) input=${OPTARG};;
                esac
done
echo $number
echo $input
echo $tag

for protein in h1 h3 n1 n2 n1pand h1pand
	do
	perl test_single_sites.pl --state nsyn --input $input --stattype median --protein $protein --simnumber $number --norm weightnorm --bigtag $tag --likelihood  --verbose --fromcsv &
	perl test_single_sites.pl --state syn --input $input --stattype median --protein $protein --simnumber $number --norm weightnorm --bigtag $tag --likelihood  --verbose --fromcsv &
        perl test_single_sites.pl --state nsyn --input $input --stattype mean --protein $protein --simnumber $number --norm weightnorm --bigtag $tag --likelihood  --verbose --fromcsv &
        perl test_single_sites.pl --state syn --input $input --stattype mean --protein $protein --simnumber $number --norm weightnorm --bigtag $tag --likelihood  --verbose --fromcsv &
	done
wait
