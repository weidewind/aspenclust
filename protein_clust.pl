#!/usr/bin/perl

# screen -U -S sites parallel perl protein_clust.pl --sites --simnumber 1000 --output include_sisters_1000 --protein {1} --state {2} ::: h1 h3 n1 n2 h1pandemic n1pandemic ::: nsyn syn

use Aspens;
use Mutmap;
use Groups;
use File::Spec;
use Cwd qw(abs_path cwd getcwd);
use lib getcwd(); # adds working directory to @INC
use Groups;
use Data::Dumper;
use Getopt::ArgvFile;
use Getopt::Long;

my $protein;
my $state = 'nsyn';
my $output = '';	# option variable with default value
my $simnumber = 10;
my $all;
my $sites;
my $subset;
my $groups;

GetOptions (	
		'protein=s' => \$protein,
		'state=s' => \$state,
		'output=s' => \$output,
		'simnumber=i' => \$simnumber,
		'subset=s' => \$subset,
		'sites' => \$sites,
		'groups' => \$groups,
		'all' => \ $all,
		);

$| = 1;

my $args = {protein => $protein, state => $state, bigtag => $output};
my $mutmap = Mutmap->new($args);
if ($all) {Aspens::logic_global_median_statistics($mutmap, 0.00005, $simnumber);}
if ($sites){ 
	if ($subset){
		my @sites_subset = split(/,/,$subset);
		Aspens::logic_median_statistics($mutmap, 0.00005, $simnumber,\@sites_subset);
	}
	else { Aspens::logic_median_statistics($mutmap, 0.00005, $simnumber);}
}
if ($groups){
 	my @groups_and_names = Groups::get_predefined_groups_and_names_for_protein($protein, $mutmap->aa_length());
	my @groups = @{$groups_and_names[0]};
	my @names = @{$groups_and_names[1]};
	my $item = 0;
	while($item < scalar @groups){
		my @group = @{$groups[$item]};
		my $name = $names[$item];
		Aspens::logic_medstat_groups_labelshuffler($mutmap, 0.00005, $simnumber, $name, \@group);
		Aspens::logic_medstat_groups($mutmap, 0.00005, $simnumber, $name, \@group);
		$item++;
	}
}
