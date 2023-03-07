#!/usr/bin/perl

use DistanceMap;
use File::Spec;
use Cwd qw(abs_path cwd getcwd);
use lib getcwd(); # adds working directory to @INC
use Getopt::ArgvFile;
use Getopt::Long;
use Data::Dumper;
use Parsers;
use Parallel::ForkManager;
use File::Path qw(make_path remove_tree);

my $protein = "toy";
my $state = "nsyn";
my $bigtag = "test";
my $phenotype;
my $input;
my $verbose;
my $switch;
my $pairsfile;

GetOptions (	
		'protein=s' => \$protein,
		'state=s' => \$state,
		'bigtag=s' => \$bigtag,
		'input=s' => \$input,
		'verbose' => \$verbose,
		'switch' => \$switch,
		'phenotype=s' => \$phenotype,
		'pairsfile=s' => \$pairsfile,
);

$| = 1;

my $args = {protein => $protein, state => $state, bigtag => $bigtag, input => $input, fdr => 1, phenotype => $phenotype};
my $folder = DistanceMap::get_input_base($args);
make_path(DistanceMap::get_output_base($args));
my $log = File::Spec->catfile(DistanceMap::get_output_base($args), "$protein.log");

open LOG, ">$log" or die "Cannot create log file $log $! \n";
print "Printing some data to logfile $log\n";

my @files = dirfiles($folder);
my @statfiles = grep {/^[0-9]+\.stat\.ident/} @files;
print LOG Dumper (\@statfiles );
my @fakenums;
for my $filename(@statfiles){
		push @fakenums, $filename =~ /([0-9]+)/;
}
print LOG Dumper (\@fakenums);
##
## Preparing commands for Forkmanager
my @commands;



foreach my $fakenum(@fakenums){	
	#check if this fake is already processed
	my $outfile = File::Spec->catfile(DistanceMap::get_output_base($args), $protein."_".$fakenum."_all_nonsequential");
	if (-e $outfile){
		print LOG "$outfile already exists. Skipping\n";
		next;
	}
	my $command = mycomm($fakenum, $pairsfile);
	push @commands, $command;
	print LOG $command."\n";	
}

##
## Forkmanager setup
my $manager = new Parallel::ForkManager(30);

$manager->run_on_start( 
	  sub {
		my $pid = shift;
		print LOG "Starting child processes under process id $pid\n";
	  }
	);
$manager->run_on_finish( 
	  sub {
		 my ( $pid, $exit_code, $signal, $core ) = @_;
		 if ( $core ) {
			print LOG "Process (pid: $pid) core dumped.\n";
		 } else {
			print LOG "Process (pid: $pid) exited with code $exit_code and signal $signal.\n";
		 }
	  }
   );
$manager->run_on_wait( 
	  sub {
		 print LOG "Waiting for all children to terminate... \n";
	  },
	  180 # time interval between checks
   ); 
##
## Launching a series of new iteration_gulps if needed 

foreach my $command (@commands) {
	  $manager->start and next;
	  system( $command );
	  $manager->finish;
   }
$manager->wait_all_children;

close LOG;

##
sub mycomm {
	my $fakenum = shift;
	my $pairsfile = shift;
	my $perlocation = "perl";
	my $exports = "";
	if ($switch) {
		$perlocation = "~/perl5/perlbrew/perls/perl-5.22.1/bin/perl";
	 	$exports = "export PERL5LIB=/export/home/popova/perl5/lib/perl5/x86_64-linux:/export/home/popova/perl5/lib/perl5:/export/home/popova/.perl/lib/perl5/5.22.1/x86_64-linux:/export/home/popova/.perl/lib/perl5/5.22.1:/export/home/popova/.perl/lib/perl5/x86_64-linux:/export/home/popova/.perl/lib/perl5:/export/home/popova/perl5/lib/perl5/x86_64-linux:/export/home/popova/perl5/lib/perl5:/export/home/popova/perl5/lib/perl5/x86_64-linux:/export/home/popova/perl5/lib/perl5:/export/home/popova/workspace/evopoisson:$PERL5LIB; ";
	}
	my $command = $exports.$perlocation." test_distances.pl --input $input --protein $protein --state $state --bigtag $bigtag --fdr --fakenum $fakenum";
	if ($pairsfile){
		$command .= " --pairsfile $pairsfile";
	}
	if ($phenotype){
		$command .= " --phenotype $phenotype";
	}
	return $command;
	
}

sub dirfiles{
        my $dirname = shift;
        opendir(DH, $dirname);
        my @files = readdir(DH);
        closedir(DH);
        if (scalar @files == 0){
                print LOG "No files found in $dirname!\n";
				print "No files found in $dirname!\n";
        }
        return @files
}
