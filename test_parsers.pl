#!/usr/bin/perl

use Parsers qw(parse_likelihoods parse_fasta_as_likelihoods parse_tree parse_fasta parse_likelihoods_csv);
use Data::Dumper;
use List::Util qw(sum shuffle);

my ($nodeseqs, $string) = parse_likelihoods_csv("../ksuAncestorPipeline/n1pandemic/ancestors.csv");
foreach my $key (keys %{$nodeseqs}){
	print Dumper ($nodeseqs->{$key}->[2]);
	print "\n";
	print sum(@{$nodeseqs->{$key}->[2]})."\n"; # must be 1
	foreach my $site(@{$nodeseqs->{$key}}){
		if ($nodeseqs->{$key}->[$site]->[2] > 0.000000001 & $nodeseqs->{$key}->[$site]->[2] < 1){
			print ($nodeseqs->{$key}->[$site]->[2])
		}
	}
}