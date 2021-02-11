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


my $args = {protein => "toy", state => "nsyn", bigtag => "test"};
my $mutmap = ProbsMutmap->new($args);


my @sites_subset = (2);
Aspens::single_sites_stats($mutmap, 0.00005, 100, \@sites_subset); # Works!





