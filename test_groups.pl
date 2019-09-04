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
my $stattype = "";
my $group = "1,2,3,4";
my $verbose;
my $likelihood;

GetOptions (	
		'protein=s' => \$protein,
		'state=s' => \$state,
		'simnumber=i' => \$simnumber,
		'norm=s' => \$norm,
		'bigtag=s' => \$bigtag,
		'stattype=s' => \$stattype,
		'verbose' => \$verbose,
		'group=s' => \$group,
		'likelihood' => \$likelihood,
);
$| = 1;
my $args = {protein => $protein, state => $state, bigtag => $bigtag, likelihood => $likelihood};
my $mutmap = ProbsMutmap->new($args);
my @sites = split(",", $group);
print STDOUT "Probsmutmap created\n";
## todo: can be rewritten for efficient computation of both types of stats
my $time0 = clock();
if ($stattype){
	Aspens::group_stats({mutmap => $mutmap, shuffletype => "labelshuffler", simnumber => $simnumber, stattype => $stattype, norm => $norm, verbose => $verbose, group => \@sites});
	Aspens::group_stats({mutmap => $mutmap, shuffletype => "siteshuffler",simnumber => $simnumber, stattype => $stattype, norm => $norm, verbose => $verbose, group => \@sites});
}
else{
	Aspens::group_stats({mutmap => $mutmap, shuffletype => "labelshuffler",simnumber => $simnumber, stattype => "mean", norm => $norm, verbose => $verbose, group => \@sites});
	Aspens::group_stats({mutmap => $mutmap, shuffletype => "labelshuffler",simnumber => $simnumber, stattype => "median", norm => $norm, verbose => $verbose, group => \@sites});
	Aspens::group_stats({mutmap => $mutmap, shuffletype => "siteshuffler",simnumber => $simnumber, stattype => "mean", norm => $norm, verbose => $verbose, group => \@sites});
	Aspens::group_stats({mutmap => $mutmap, shuffletype => "siteshuffler",simnumber => $simnumber, stattype => "median", norm => $norm, verbose => $verbose, group => \@sites});
}

my $timedone1 = clock();
my $exectime = $timedone1-$time0;
print "Done in $exectime";

