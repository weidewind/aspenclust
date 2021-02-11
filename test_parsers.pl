#!/usr/bin/perl

use Parsers qw(parse_likelihoods parse_fasta_as_likelihoods parse_tree parse_fasta parse_likelihoods_csv);
use Data::Dumper;
use List::Util qw(sum shuffle);

my ($nodeseqs, $string) = parse_likelihoods_csv("../ksuAncestorPipeline/n1pandemic/ancestors.csv");
foreach my $key (keys %{$nodeseqs}){
	print "\n";
	print sum(@{$nodeseqs->{$key}->[2]})."\n"; # must be 1
	my $length = scalar @{$nodeseqs->{$key}}-1;
	foreach my $site(0..$length){
		my $probs = $nodeseqs->{$key}->[$site];
		if ($probs->[2] > 0.001 & $probs->[2] < 1){
			print Dumper ($probs);
		}
	}
}