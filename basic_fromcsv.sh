#!/bin/bash
while getopts p:t:s:n:i: option
do
                case "${option}" in
                        p) protein=${OPTARG};;
                        t) tag=${OPTARG};;
			s) state=${OPTARG};;
			n) number=${OPTARG};;
			i) input=${OPTARG};;
                esac
done
norm=countnorm
echo $protein
echo $number
echo $input
echo $norm
if test -z "$protein"
then
	for protein in h1 h3 n1 n2 n1pand h1pand
	do
	perl test_single_sites.pl --state $state --input $input --stattype mean --protein $protein --simnumber $number --norm $norm --bigtag $tag --likelihood  --fromcsv &
	perl test_global.pl --state $state --input $input --stattype mean --protein $protein --simnumber $number --norm $norm --bigtag $tag --likelihood   --fromcsv &
	done
	wait
else
	perl test_single_sites.pl --state $state --input $input --stattype mean --protein $protein --simnumber $number --norm $norm --bigtag $tag --likelihood --fromcsv &
	perl test_global.pl --state $state --input $input --stattype mean --protein $protein --simnumber $number --norm $norm --bigtag $tag --likelihood --fromcsv &
	perl test_single_sites.pl --state $state --input $input --stattype median --protein $protein --simnumber $number --norm $norm --bigtag $tag --likelihood --fromcsv &
	perl test_global.pl --state $state --input $input --stattype median --protein $protein --simnumber $number --norm $norm --bigtag $tag --likelihood --fromcsv &
	wait
fi
