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
my $verbose;
my $likelihood;
my $fromcsv;

GetOptions (	
		'protein=s' => \$protein,
		'state=s' => \$state,
		'simnumber=i' => \$simnumber,
		'norm=s' => \$norm,
		'bigtag=s' => \$bigtag,
		'stattype=s' => \$stattype,
		'verbose' => \$verbose,
		'likelihood' => \$likelihood,
		'fromcsv' => \$fromcsv,
);
$| = 1;
my ($groups, $names) = Groups::only_groups($protein);
my %groupnames;
for (my $i = 0; $i < scalar @{$groups}; $i++){
	$groupnames{$names->[$i]} = $groups->[$i];
}
print Dumper(\%groupnames);
my $args = {protein => $protein, state => $state, bigtag => $bigtag, likelihood => $likelihood, fromcsv => $fromcsv};
my $mutmap = ProbsMutmap->new($args);
print STDOUT "Probsmutmap created\n";
## todo: can be rewritten for efficient computation of both types of stats
my $time0 = clock();
if ($stattype){
	my $out1 = Aspens::group_stats_batch({mutmap => $mutmap, shuffletype => "labelshuffler", simnumber => $simnumber, stattype => $stattype, norm => $norm, verbose => $verbose, groupnames => \%groupnames});
	my $out2 = Aspens::group_stats_batch({mutmap => $mutmap, shuffletype => "siteshuffler",simnumber => $simnumber, stattype => $stattype, norm => $norm, verbose => $verbose, groupnames => \%groupnames});
	print Dumper ($out1);
	}
else{
	Aspens::group_stats_batch({mutmap => $mutmap, shuffletype => "labelshuffler",simnumber => $simnumber, stattype => "mean", norm => $norm, verbose => $verbose, groupnames => \%groupnames});
	Aspens::group_stats_batch({mutmap => $mutmap, shuffletype => "labelshuffler",simnumber => $simnumber, stattype => "median", norm => $norm, verbose => $verbose, groupnames => \%groupnames});
	Aspens::group_stats_batch({mutmap => $mutmap, shuffletype => "siteshuffler",simnumber => $simnumber, stattype => "mean", norm => $norm, verbose => $verbose, groupnames => \%groupnames});
	Aspens::group_stats_batch({mutmap => $mutmap, shuffletype => "siteshuffler",simnumber => $simnumber, stattype => "median", norm => $norm, verbose => $verbose, groupnames => \%groupnames});
}

my $timedone1 = clock();
my $exectime = $timedone1-$time0;
print "Done in $exectime";

