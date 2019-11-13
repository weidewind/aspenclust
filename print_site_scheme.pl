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
use IPC::System::Simple qw(capture);

my $protein = "";
my $state = 'nsyn';
my $output = '';        # option variable with default value
my $sites = "235,189,278";
my $type = "probtree";
my $fromcsv;
my $draw;
my $input;

GetOptions (    'protein=s' => \$protein,
                'state=s' => \$state,
                'output=s' => \$output,
                'sites=s' => \$sites,
		'type=s' => \$type,
		'input=s' => \$input,
		'fromcsv' => \$fromcsv,
		'draw' => \$draw,
        );



my $args = {protein => $protein, input => $input, state => $state, bigtag => $output, likelihood => 1, fromcsv => $fromcsv};
my $mutmap = ProbsMutmap->new($args);


my @sites = split(",", $sites);
foreach my $site (@sites){
	print $site."\n";
	if ($type eq "graph"){
		Aspens::print_scheme($mutmap,$site);
		if ($draw){die "cannot draw graph";}
	}
	elsif ($type eq "probtree"){
		Aspens::print_probtree_scheme($mutmap,$site);
		if ($draw){
			my $outpref = "output/$input/$output/$state/likelihoods/trees/$protein"."_".$site;
			my $command = "xvfb-run python drawProbAspenTree.py --treefile data/$input/$protein.ancestor.likelihoods.withinternals.newick --eventfile $outpref.probtree --output $outpref --width 870 --circle_size 20 --large";
	                my $logs = capture($command);
        	        print $logs."\n";
		}
	}
	else {
		die "Unknown type $type : exprected graph or probtree";
	}
}

