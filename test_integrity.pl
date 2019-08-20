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
use Carp::Assert;
use Parsers;

my $protein = "toy";
my $state = "nsyn";
my $simnumber = 100;
my $norm = "weightnorm"; #weightnorm or countnorm
my $bigtag = "test";
my $stattype = "median";
my $verbose = 1;
my $group = "1,2,3,4";

GetOptions (	
#		'protein=s' => \$protein,
#		'state=s' => \$state,
#		'simnumber=i' => \$simnumber,
		'norm=s' => \$norm,
		'bigtag=s' => \$bigtag,
		'stattype=s' => \$stattype,
#		'verbose' => \$verbose,
#		'group=s' => \$group,
);

my $args = {protein => $protein, state => $state, bigtag => $bigtag};
my $mutmap = ProbsMutmap->new($args);
my @sites = split(",", $group);

my $time0 = clock();
my $globalfile = Aspens::global_stats({mutmap => $mutmap, simnumber => $simnumber, stattype => $stattype, norm => $norm, verbose => $verbose});
my $timedone = clock();
my $exectime = $timedone-$time0;
print "Global test done in $exectime\n";
my $time0 = clock();
my $labelshufflerfile = Aspens::group_stats({mutmap => $mutmap, shuffletype => "labelshuffler", simnumber => $simnumber, stattype => $stattype, norm => $norm, verbose => $verbose, group => \@sites});
my $siteshufflerfile = Aspens::group_stats({mutmap => $mutmap, shuffletype => "siteshuffler",simnumber => $simnumber, stattype => $stattype, norm => $norm, verbose => $verbose, group => \@sites});
my $timedone = clock();
my $exectime = $timedone-$time0;
print "Group tests done in $exectime\n";

open G, "<$globalfile" or die "CAnnot open $globalfile $!";
my ($glsame, $gldiff);
while (<G>){
 if ($_ =~ "^Sum"){
	($glsame, $gldiff) = (split /\s+/, $_)[2,5]; 
	last;
 }
}
close G;

open G, "<$labelshufflerfile" or die "CAnnot open $labelshufflerfile $!";
my @array = <G>;
my @strs = grep {$_ =~ "^Sum"} (@array);
my ($grsame_labelshuffler, $grdiff_labelshuffler) = ( split /\s+/, $strs[0] )[2,5];
my ($csame_labelshuffler, $cdiff_labelshuffler) = (split /\s+/, $strs[1] )[2,5];
close G;

open G, "<$siteshufflerfile" or die "CAnnot open $siteshufflerfile $!";
@array = <G>;
@strs = grep {$_ =~ "^Sum"} (@array);
my ($grsame_siteshuffler, $grdiff_siteshuffler) = (split /\s+/, $strs[0])[2,5];
my ($csame_siteshuffler, $cdiff_siteshuffler) = (split /\s+/, $strs[1])[2,5];
close G;

print join(",", $glsame, $gldiff, $grsame_labelshuffler,  $grdiff_labelshuffler, $csame_labelshuffler, $cdiff_labelshuffler);
# weighted sum of observations (both for parallel and divergent subs) is the same for whole protein and for sum of group and its complement
assert(abs($glsame  - ($grsame_labelshuffler +$csame_labelshuffler)) < 0.0001, "For labelshuffler, the total number of observations for parallel substitutions (OPS) equals to the sum of OPS for group and its complement");
assert(abs($gldiff  - ($grdiff_labelshuffler +$cdiff_labelshuffler)) < 0.0001, "For labelshuffler, the total number of observations for divergent substitutions (ODS) equals to the sum of ODS for group and its complement");
assert(abs($glsame  - ($grsame_siteshuffler +$csame_siteshuffler)) < 0.0001, "For siteshuffler, the total number of observations for parallel substitutions (OPS) equals to the sum of OPS for group and its complement");
assert(abs($gldiff  - ($grdiff_siteshuffler +$cdiff_siteshuffler)) < 0.0001, "For siteshuffler, the total number of observations for divergent substitutions (ODS) equals to the sum of ODS for group and its complement");
if ($norm == "weightnorm"){
	assert(abs($glsame - 65) < 0.0001, "Total number of observations for parallel substitutions is close to 65");
	assert(abs($gldiff - 65) < 0.0001, "Total number of observations for divergent substitutions is close to 65");
}


