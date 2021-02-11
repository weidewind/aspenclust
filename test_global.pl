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
use Time::HiRes qw( clock );

use Parsers;

my $protein = "toy";
my $state = "nsyn";
my $simnumber = 100;
my $norm = "weightnorm"; #weightnorm or countnorm
my $bigtag = "test";
my $stattype = "median";
my $verbose;
my $likelihood;
my $input;
my $fromcsv;
my $date;

GetOptions (	
		'protein=s' => \$protein,
		'state=s' => \$state,
		'simnumber=i' => \$simnumber,
		'norm=s' => \$norm,
		'bigtag=s' => \$bigtag,
		'stattype=s' => \$stattype,
		'verbose' => \$verbose,
		'likelihood' => \$likelihood,
		'input=s' => \$input,
		'fromcsv' => \$fromcsv,
		'date' => \$date,
);
$| = 1;
my $args = {protein => $protein, state => $state, bigtag => $bigtag, likelihood => $likelihood, input => $input, fromcsv => $fromcsv, distance_in_date => $date};
my $mutmap = ProbsMutmap->new($args);
print STDOUT "Probsmutmap created\n";

my $time0 = clock();
Aspens::global_stats({mutmap => $mutmap, simnumber => $simnumber, stattype => $stattype, norm => $norm, verbose => $verbose});
my $timedone = clock();
my $exectime = $timedone-$time0;
print "Done in $exectime";



