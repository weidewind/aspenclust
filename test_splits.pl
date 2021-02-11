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
my $stattype = "mean";
my $sites;
my $verbose;
my $likelihood;
my $input;
my $fromcsv;
my $date;
my $splits;
my $analysis="all";
my $killist; # 137:ATT:M,C;144:TGG:F,Y,A


GetOptions (	
		'protein=s' => \$protein,
		'state=s' => \$state,
		'simnumber=i' => \$simnumber,
		'norm=s' => \$norm,
		'bigtag=s' => \$bigtag,
		'stattype=s' => \$stattype,
		'sites=s' => \$sites,
		'verbose' => \$verbose,
		'likelihood' => \$likelihood,
		'input=s' => \$input,
		'fromcsv' => \$fromcsv,
		'date' => \$date,
		'splits=s' => \$splits,
		'analysis=s' => \$analysis,
		'killist=s' =>\$killist,
);
$| = 1;
my $args = {protein => $protein, state => $state, bigtag => $bigtag, likelihood => $likelihood, input => $input, fromcsv => $fromcsv, distance_in_date => $date, splits => $splits};
my $mutmap = ProbsMutmap->new($args);
print STDOUT "Probsmutmap created\n";

my @sites = split(',', $sites);
unless (@sites) {@sites = (1..$mutmap->{static_length});}

my ($groups, $names) = Groups::only_groups($protein);
my %groupnames;
for (my $i = 0; $i < scalar @{$groups}; $i++){
        $groupnames{$names->[$i]} = $groups->[$i];
}

my $temp = $killist;
my %killist;
my @t = split(';', $temp);
foreach my $s (@t){
	my @s = split(":", $s);
	my @ders = split(',',$s[2]);
	$killist{$s[0]}{$s[1]} = \@ders;
}
print Dumper \%killist;

my $time0 = clock();

if ($analysis eq "all" || $analysis eq "global"){
	Aspens::global_stats({mutmap => $mutmap, simnumber => $simnumber, stattype => $stattype, norm => $norm, verbose => $verbose});
}
if ($analysis eq "all" || $analysis eq "single_sites"){
	Aspens::single_sites_stats({mutmap =>$mutmap, simnumber => $simnumber, sites => \@sites, verbose => $verbose, stattype => $stattype, norm =>$norm, killist => \%killist});
}
if ($analysis eq "all" || $analysis eq "groups"){
	Aspens::group_stats_batch({mutmap => $mutmap, shuffletype => "labelshuffler", simnumber => $simnumber, stattype => $stattype, norm => $norm, verbose => $verbose, groupnames => \%groupnames});
	Aspens::group_stats_batch({mutmap => $mutmap, shuffletype => "siteshuffler",simnumber => $simnumber, stattype => $stattype, norm => $norm, verbose => $verbose, groupnames => \%groupnames});
}
if ($analysis eq "all" || $analysis eq "alleles"){
	Aspens::alleles_stats({mutmap =>$mutmap, simnumber => $simnumber, sites => \@sites, verbose => $verbose, stattype => $stattype, norm =>$norm}); 
}

my $timedone = clock();
my $exectime = $timedone-$time0;
print "Done in $exectime";



