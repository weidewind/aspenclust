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
Aspens::groups_stats_labelshuffler($mutmap, 500, $tag, [1,2,3], "mean");
Aspens::groups_stats_labelshuffler($mutmap, 500, $tag, [1,2,3], "median");
my $timedone1 = clock();
Aspens::groups_stats_siteshuffler($mutmap, 500, $tag, [1,2,3], "median");
Aspens::groups_stats_siteshuffler($mutmap, 500, $tag, [1,2,3], "mean");
my $timedone2 = clock();
my $exectime1 = $timedone1-$time0;
my $exectime2 = $timedone2-$timedone1;
print "Done in $exectime1 $timedone2";

