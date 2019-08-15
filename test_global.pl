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

my $time0 = clock();
Aspens::global_stats($mutmap, 0.00005, 100);
my $timedone = clock();
my $exectime = $timedone-$time0;
print "Done in $exectime";



