#!/bin/bash
perl protein_clust.pl -p h3 --state nsyn --simnumber 1000 --output test &
perl protein_clust.pl -p h1 --state nsyn --simnumber 1000 --output test &
perl protein_clust.pl -p n1 --state nsyn --simnumber 1000 --output test &
perl protein_clust.pl -p n2 --state nsyn --simnumber 1000 --output test &
perl protein_clust.pl -p h1pandemic --state nsyn --simnumber 1000 --output test &
perl protein_clust.pl -p n1pandemic --state nsyn --simnumber 1000 --output test &
perl protein_clust.pl -p h3 --state syn --simnumber 1000 --output test &
perl protein_clust.pl -p h1 --state syn --simnumber 1000 --output test &
perl protein_clust.pl -p n1 --state syn --simnumber 1000 --output test &
perl protein_clust.pl -p n2 --state syn --simnumber 1000 --output test &
perl protein_clust.pl -p h1pandemic --state syn --simnumber 1000 --output test &
perl protein_clust.pl -p n1pandemic --state syn --simnumber 1000 --output test &
wait
