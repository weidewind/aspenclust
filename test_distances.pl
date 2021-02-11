#!/usr/bin/perl

use strict;
use DistanceMap;
use File::Spec;
use Cwd qw(abs_path cwd getcwd);
use lib getcwd(); # adds working directory to @INC
use Getopt::ArgvFile;
use Getopt::Long;
use Data::Dumper;
use Parsers;

my $protein = "toy";
my $state = "nsyn";
my $bigtag = "test";
my $input;
my $verbose;
my $fdr;
my $fakenum;
my $pairsfile; #contains pairs of sites, one tab-delimited pair per line (e.g. 123	144\n)

GetOptions (	
		'protein=s' => \$protein,
		'state=s' => \$state,
		'bigtag=s' => \$bigtag,
		'input=s' => \$input,
		'verbose' => \$verbose,
		'fdr' => \$fdr,
		'fakenum=i' => \$fakenum,
		'pairsfile=s' => \$pairsfile,
);

$| = 1;

my $args = {protein => $protein, state => $state, bigtag => $bigtag, input => $input, fdr => $fdr, fakenum => $fakenum};
print Dumper $args;
my $mutmap = DistanceMap->new($args);

print STDOUT "DistanceMap created\n";

if ($pairsfile){
	my @pairs;
	open P, "<$pairsfile" or die "cannot open file with pairs $pairsfile $!";
	my $header = <P>;
	while (<P>){
		my @sp = split(/\s+/);
		push @pairs, [$sp[0], $sp[1]];
	}
	close P;
	print Dumper (\@pairs);
	$mutmap->all_nonsequential_distances(\@pairs);
}
else{$mutmap->all_nonsequential_distances();}




