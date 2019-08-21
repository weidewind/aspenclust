#!/usr/bin/perl

use strict;
use Aspens;
use ProbsMutmap;
use Groups;
use File::Spec;
use Cwd qw(abs_path cwd getcwd);
use lib getcwd(); # adds working directory to @INC
use Getopt::ArgvFile;
use Getopt::Long;
use Data::Dumper;
use Parsers;

my $protein = "toy";
my $state = "nsyn";
my $simnumber = 100;
my $norm = "weightnorm"; #weightnorm or countnorm
my $bigtag = "test";
my $stattype = "median";
my $verbose;
my $sites;

GetOptions (	
		'protein=s' => \$protein,
		'state=s' => \$state,
		'simnumber=i' => \$simnumber,
		'norm=s' => \$norm,
		'bigtag=s' => \$bigtag,
		'sites=s' => \$sites,
		'verbose' => \$verbose,
		'stattype=s' => \$stattype,

);

$| = 1;

my $args = {protein => $protein, state => $state, bigtag => $bigtag};
my $mutmap = ProbsMutmap->new($args);
my @sites = split(',', $sites);
print STDOUT "Probsmutmap created\n";

unless (@sites) {@sites = (1..$mutmap->{static_length});}
Aspens::single_sites_stats({mutmap =>$mutmap, simnumber => $simnumber, sites => \@sites, verbose => $verbose, stattype => $stattype, norm =>$norm}); # Works!
#Aspens::print_scheme($mutmap,5);





