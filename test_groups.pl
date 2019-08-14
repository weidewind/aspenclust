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


my $args = {protein => "toy", state => "nsyn", bigtag => "test_debinned"};
my $mutmap = ProbsMutmap->new($args);

## todo: can be rewritten for efficient computation of both types of stats
Aspens::groups_stats_labelshuffler($mutmap, 100, "testgroup", [1,2,3], "mean");
Aspens::groups_stats_labelshuffler($mutmap, 100, "testgroup", [1,2,3], "median");
Aspens::groups_stats_siteshuffler($mutmap, 100, "testgroup", [1,2,3], "median");
Aspens::groups_stats_siteshuffler($mutmap, 100, "testgroup", [1,2,3], "mean");


