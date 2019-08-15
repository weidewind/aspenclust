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


my $args = {protein => "toy", state => "nsyn", bigtag => "test_debinned"};
my $mutmap = ProbsMutmap->new($args);
my $tag = "refactor";

## todo: can be rewritten for efficient computation of both types of stats
my $time0 = clock();
Aspens::groups_stats_labelshuffler($mutmap, 100, $tag, [1,2,3], "mean");
Aspens::groups_stats_labelshuffler($mutmap, 100, $tag, [1,2,3], "median");
Aspens::groups_stats_siteshuffler($mutmap, 100, $tag, [1,2,3], "median");
Aspens::groups_stats_siteshuffler($mutmap, 100, $tag, [1,2,3], "mean");
my $timedone = clock();
my $exectime = $timedone-$time0;
print "Done in $exectime";

