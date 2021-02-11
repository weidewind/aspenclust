#!/usr/bin/bash

perl test_splits.pl --protein h3 --analysis single_sites --killist 137:AAC:D --site 137 --state nsyn --stattype mean --splits 4split_mega --input little_ksu --simnumber 100 --bigtag test_killistD &
perl test_splits.pl --protein h3 --analysis single_sites --killist 137:AAC:K --site 137 --state nsyn --stattype mean --splits 4split_mega --input little_ksu --simnumber 100 --bigtag test_killistK &
perl test_splits.pl --protein h3 --analysis single_sites --killist 137:AAC:T --site 137 --state nsyn --stattype mean --splits 4split_mega --input little_ksu --simnumber 100 --bigtag test_killistT &
perl test_splits.pl --protein h3 --analysis single_sites --killist 137:AAC:H --site 137 --state nsyn --stattype mean --splits 4split_mega --input little_ksu --simnumber 100 --bigtag test_killistH &
perl test_splits.pl --protein h3 --analysis single_sites --killist 137:AAC:S --site 137 --state nsyn --stattype mean --splits 4split_mega --input little_ksu --simnumber 100 --bigtag test_killistS &
wait
