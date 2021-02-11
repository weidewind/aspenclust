#!/bin/bash
prot=h1pandemic
number=10000
bigt=parenttest
perl test_single_sites.pl --stattype median --protein $prot --simnumber $number --norm weightnorm --bigtag $bigt --likelihood --fromcsv &
perl test_single_sites.pl --stattype mean --protein $prot --simnumber $number --norm weightnorm --bigtag $bigt --likelihood --fromcsv &
perl test_single_sites.pl --stattype median --protein $prot --simnumber $number --norm countnorm --bigtag $bigt --likelihood --fromcsv &
perl test_single_sites.pl --stattype mean --protein $prot --simnumber $number --norm countnorm --bigtag $bigt --likelihood --fromcsv & 
perl test_global.pl --stattype median --protein $prot --simnumber $number --norm weightnorm --bigtag $bigt --likelihood --verbose --fromcsv &
perl test_global.pl --stattype mean --protein $prot --simnumber $number --norm weightnorm --bigtag $bigt --likelihood  --verbose --fromcsv &
perl test_global.pl --stattype median --protein $prot --simnumber $number --norm countnorm --bigtag $bigt --likelihood --verbose --fromcsv &
perl test_global.pl --stattype mean --protein $prot --simnumber $number --norm countnorm --bigtag $bigt --likelihood --verbose --fromcsv &
perl test_groups_batch.pl --stattype median --protein $prot --simnumber $number --norm weightnorm --bigtag $bigt --likelihood --verbose --fromcsv &
perl test_groups_batch.pl --stattype mean --protein $prot --simnumber $number --norm weightnorm --bigtag $bigt --likelihood  --verbose --fromcsv &
perl test_groups_batch.pl --stattype median --protein $prot --simnumber $number --norm countnorm --bigtag $bigt --likelihood --verbose --fromcsv &
perl test_groups_batch.pl --stattype mean --protein $prot --simnumber $number --norm countnorm --bigtag $bigt --likelihood --verbose --fromcsv &
wait
