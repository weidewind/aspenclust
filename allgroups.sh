#!/bin/bash
prot=h3
number=10000
bigt=parenttest
perl test_groups_batch.pl --stattype median --protein $prot --simnumber $number --norm weightnorm --bigtag $bigt --likelihood --verbose &
perl test_groups_batch.pl --stattype mean --protein $prot --simnumber $number --norm weightnorm --bigtag $bigt --likelihood  --verbose &
perl test_groups_batch.pl --stattype median --protein $prot --simnumber $number --norm countnorm --bigtag $bigt --likelihood --verbose &
perl test_groups_batch.pl --stattype mean --protein $prot --simnumber $number --norm countnorm --bigtag $bigt --likelihood --verbose &
wait
