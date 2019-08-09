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

my $protein = "";
my $state = 'nsyn';
my $output = '';        # option variable with default value
my $sites = "235,189,278";

GetOptions (    'protein=s' => \$protein,
                'state=s' => \$state,
                'output=s' => \$output,
                'sites=s' => \$sites,
        );


my $args = {protein => $protein, state => $state, bigtag => $output};
my $mutmap = ProbsMutmap->new($args);


my @sites = split(",", $sites);
foreach my $site (@sites){
	print $site."\n";
	Aspens::print_scheme($mutmap,$site);
}

