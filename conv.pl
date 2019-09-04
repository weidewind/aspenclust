#!/usr/bin/perl

use Parsers;

my %ns = parse_likelihoods_as_fasta("C:/Users/weidewind/workspace/aspenclust/data/toy.ancestor.likelihoods");
my $out = "C:/Users/weidewind/workspace/aspenclust/data/toy.ancestor.fa";
open OUT, ">$out";
for my $nname (keys %ns){
	print OUT ">".$nname."\n";
	print OUT $ns{$nname}."\n";
}
close OUT;
