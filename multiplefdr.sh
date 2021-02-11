prot=h1pandemic
script=fdr_single_sites.pl
screen -U -dm -S fdraddsyn_sitesnode3 perl $scrpit --input little_ksu --state syn --stattype mean --protein $prot --simnumber 1000  --fakenumber 20 --firstfake 1 --norm countnorm --bigtag fdrtest10_1000
screen -U -dm -S fdraddsyn_sitesnode3 perl $scrpit --input little_ksu --state syn --stattype mean --protein $prot --simnumber 1000  --fakenumber 20 --firstfake 21 --norm countnorm --bigtag fdrtest10_1000
screen -U -dm -S fdraddsyn_sitesnode3 perl $scrpit --input little_ksu --state syn --stattype mean --protein $prot --simnumber 1000  --fakenumber 20 --firstfake 41 --norm countnorm --bigtag fdrtest10_1000
screen -U -dm -S fdraddsyn_sitesnode3 perl $scrpit --input little_ksu --state syn --stattype mean --protein $prot --simnumber 1000  --fakenumber 20 --firstfake 61 --norm countnorm --bigtag fdrtest10_1000
screen -U -dm -S fdraddsyn_sitesnode3 perl $scrpit --input little_ksu --state syn --stattype mean --protein $prot --simnumber 1000  --fakenumber 20 --firstfake 81 --norm countnorm --bigtag fdrtest10_1000