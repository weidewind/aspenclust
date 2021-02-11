#!/usr/bin/perl

use Parsers;

my %ns = parse_likelihoods_as_fasta("data/rerooted/n1.ancestor.likelihoods.csv", "csv");
my $out = "data/rerooted/n1.ancestor.fa";
open OUT, ">$out";
for my $nname (keys %ns){
	print OUT ">".$nname."\n";
	print OUT $ns{$nname}."\n";
}
close OUT;
