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
use Parallel::ForkManager;

my $protein = "toy";
my $state = "nsyn";
my $fakenumber = 100;
my $simnumber = 1000;
my $norm = "weightnorm"; #weightnorm or countnorm
my $bigtag = "test";
my $stattype = "median";
my $verbose;
my $sites;
my $likelihood;
my $input;
my $fromcsv;
my $date;
my $greponly;
my $splits;
my $firstfake=1;

GetOptions (	
		'protein=s' => \$protein,
		'state=s' => \$state,
		'fakenumber=i' => \$fakenumber,
		'simnumber=i' => \$simnumber,
		'norm=s' => \$norm,
		'bigtag=s' => \$bigtag,
		'sites=s' => \$sites,
		'verbose' => \$verbose,
		'stattype=s' => \$stattype,
		'likelihood' => \$likelihood,
		'input=s' => \$input,
		'fromcsv' => \$fromcsv,
		'date' => \$date,
		'greponly' => \$greponly,
		'splits=s' => \$splits,
		'firstfake=i' => \$firstfake,
);

$| = 1;

my $args = {protein => $protein, state => $state, bigtag => $bigtag, likelihood => $likelihood, input => $input, fromcsv => $fromcsv, distance_in_date =>$date, splits => $splits};

unless ($greponly){
		my $mutmap = ProbsMutmap->new($args);
		my @sites = split(',', $sites);
		unless (@sites) {@sites = (1..$mutmap->{static_length});}
		print STDOUT "Probsmutmap created\n";

		my $pm = Parallel::ForkManager->new($fakenumber);
		 
		FAKES:
		for my $fakenum ($firstfake..($fakenumber+$firstfake-1)) {
			$pm->start and next FAKES;
			print "Starting fake $fakenum\n";
			srand($fakenum);
			Aspens::alleles_stats({mutmap =>$mutmap, simnumber => $simnumber, sites => \@sites, verbose => $verbose, stattype => $stattype, norm =>$norm, fake => $fakenum});
			print "Fake $fakenum done\n";
			$pm->finish;
		};

		$pm->wait_all_children;
}

my $results_folder = File::Spec->catdir(ProbsMutmap::get_output_base($args),$norm);
opendir(DIR, $results_folder);
my @files = readdir(DIR); 
unless (scalar @files > 0){die "No simulation files found in folder $results_folder\n";}
closedir(DIR);

my $fname = $protein."_".$state."_".$stattype."_alleles_FDR";
my $fmask = $protein."_".$state."_alleles_".$stattype."_statistics_FDR_";
my $joined = File::Spec->catfile($results_folder, $fname);
open JOIN, ">>$joined" or die "CAnnot create $joined\n";

for my $f (@files){
	my $ff = File::Spec->catfile($results_folder,$f);
	next unless (-f $ff && $f =~ /^$fmask/);
	open F, "<$ff" or die "CAnnot open $ff\n";
	while (<F>){
		if ($_ =~ /^[0-9]/) {print JOIN $_; }
	}
	close F;
}
close JOIN;



