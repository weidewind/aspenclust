#!/usr/bin/perl
use strict;
use Sampler qw(print_probs_for_pair_index sample_pheno_for_pair_index);
use Parsers qw(parse_site2pheno);
use Data::Dumper;
use DistanceMap;
use Getopt::ArgvFile;
use Getopt::Long;

my $protein = "toy";
my $state = "nsyn";
my $bigtag = "test";
my $input;
my $verbose;
my $fdr;
my $fakenum;
my $phenotype;


GetOptions (	
		'protein=s' => \$protein,
		'state=s' => \$state,
		'bigtag=s' => \$bigtag,
		'input=s' => \$input,
		'verbose' => \$verbose,
		'fdr' => \$fdr,
		'fakenum=i' => \$fakenum,
		'phenotype=s' => \$phenotype,
);

$| = 1;

my $args = {protein => $protein, state => $state, bigtag => $bigtag, input => $input, fdr => $fdr, fakenum => $fakenum, phenotype => $phenotype};
my $pheno_file = DistanceMap::get_site2pheno_filepath($args);
my ($s2p, $phenotypes) = parse_site2pheno($pheno_file);
print Dumper $s2p->[0..5];

my @pairs = ([0,5],[1,2], [2,3],[3,4],[1,5]);
print Dumper \@pairs;
my $sampler = Sampler->new(\@pairs, $phenotypes, $s2p);
print "probs for pair (1,2)\n";
$sampler->print_probs_for_pair_index(1);
for (my $i = 0; $i < 20; $i++){
	my $ph = $sampler->sample_pheno_for_pair_index(1);
	print $ph." "
}


