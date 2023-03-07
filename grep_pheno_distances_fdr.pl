#!/usr/bin/perl
use strict;
use DistanceMap;
use File::Spec;
use Cwd qw(abs_path cwd getcwd);
use lib getcwd(); # adds working directory to @INC
use Getopt::ArgvFile;
use Getopt::Long;
use Data::Dumper;
use Statistics::Basic qw(:all);
use List::Util qw(sum);
$Statistics::Basic::IPRES = 6;
use Storable 'dclone';
use File::ReadBackwards;
use Parsers qw(parse_pairs parse_site2pheno);
use Sampler;




my $protein = "toy";
my $state = "nsyn";
my $bigtag = "test";
my $input;
my $sites;
my $groups;
my $nopdb;
my $pairsfile;

GetOptions (	
		'protein=s' => \$protein,
		'state=s' => \$state,
		'bigtag=s' => \$bigtag,
		'input=s' => \$input,
		'sites' => \$sites,
		'groups' => \$groups,
		'nopdb' => \$nopdb,
		'pairsfile=s' => \$pairsfile,
);

$| = 1;


my $args = {protein => $protein, state => $state, bigtag => $bigtag, input => $input};

my $log = File::Spec->catfile(DistanceMap::get_output_base($args), "$protein.grep_pheno_distances_fdr.log");
open LOG, ">$log" or die "Cannot create log file $log $! \n";
LOG->autoflush(1);
print "Printing some data to logfile $log\n";
my $errlog = File::Spec->catfile(DistanceMap::get_output_base($args), "$protein.grep_pheno_distances_fdr.errlog");
open STDERR, ">$errlog" or die "Cannot create log file $errlog $! \n";
print "Printing errors to $errlog\n";

if (!$sites || $groups || !$nopdb) {
	die "Sorry, for now only '--sites' option is available, and '--nopdb' is required. '--groups' option is not supported (yet) \n";
}

my $fdrfolder = DistanceMap::get_output_base({protein => $protein, state => $state, bigtag => $bigtag, input => $input, fdr => 1});
my $outfolder = DistanceMap::get_output_base($args) ;
my $datafile = File::Spec->catfile($outfolder, $protein."_all_nonsequential");
my %datahash;
open DAT, "<$datafile" or die "cannot open data file $datafile\n";
while (<DAT>){
	my ($site1, $site2, $totaldist, $totalpairs, $seqdist, $seqpairs, $nonseqdist, $nonseqpairs, $mean_nonseqdist) = split(/\s+/);
	$datahash{$site1}{$site2} = [$nonseqdist,$nonseqpairs];
#	$pairshash{$site1}{$site2} = $nonseqpairs;
}
close DAT;
print LOG "Finished reading the data\n";

my $pdbfile = File::Spec->catfile(DistanceMap::get_input_base($args), "pdb", $protein.".pdb_dists") unless $nopdb;
my $spfile = File::Spec->catfile(DistanceMap::get_input_base($args), $protein.".site_pairs");

my %fdrhash;
#my %fdrdetailed;
my @phenotypes = (0);
@phenotypes = dirfolders($fdrfolder);
if (scalar @phenotypes == 0){
	die "No phenotype folders found in $fdrfolder!\n";
}
foreach my $ph (@phenotypes){
	print LOG "Reading fakes for phenotype ".$ph.".. \n";
	#my($fdrhash_ph, $fdrdetailed_ph) = hash_folder(File::Spec->catfile($fdrfolder, $ph), $protein);
	my $fdrhash_ph = hash_folder(File::Spec->catfile($fdrfolder, $ph), $protein);
	$fdrhash{$ph} = $fdrhash_ph;
	#$fdrdetailed{$ph} = $fdrdetailed_ph;
}
my $pheno_file = DistanceMap::get_site2pheno_filepath($args);
my ($s2p, $phenotypes) = parse_site2pheno($pheno_file);
my @pairs = parse_pairs($pairsfile);
my $sampler = Sampler->new(\@pairs, \@phenotypes, $s2p);


if ($sites){
	my $fdroutput = File::Spec->catfile($outfolder, $protein."_dist.fdr_results");
	my %pdbdists = collect_pdbdists($pdbfile, $spfile) unless $nopdb;
	open OUT, ">$fdroutput" or die "cannot open $fdroutput for writing\n";
	print OUT "site1 site2 nonseqpairs meandist fdrmean fdrstddev fakesnum zscore lpval upval pdb\n";
	
	
	#my @sites1 = sort (keys %{$fdrhash{$phenotypes[0]}});
	#foreach my $site1 (@sites1){
		#my @sites2 = sort (keys %{$fdrhash{$phenotypes[0]}{$site1}});
		#foreach my $site2 (@sites2){
		for (my $i = 0; $i < scalar @pairs; $i++){
			my $site1 = $pairs[$i]->[0];
			my $site2 = $pairs[$i]->[1];
			## expect that it is the same for all phenotypes
			my $fakesnum = 0;
			foreach my $ph (@phenotypes){
				my $t = scalar @{$fdrhash{$ph}{$site1}{$site2}};
				if ($t > $fakesnum){
					$fakesnum = $t;
				}
			}
			# my $fakesnum = scalar @{$fdrhash{$phenotypes[0]}{$site1}{$site2}};
			# print "Sampling ".$fakesnum." fakes for sites ".$site1." and ".$site2."\n";
			my @fakes;
			for (my $f = 0; $f < $fakesnum; $f++){
				my $ph = $sampler->sample_pheno_for_pair_index($i);
				#my $ph = 1;
				push @fakes, $fdrhash{$ph}{$site1}{$site2}[$f];
			}
			@fakes = grep {/[0-9\.]+/} @fakes; # dunno why
			my $fdrsum = sum(@fakes);
			my $fdrmean = mean(@fakes);
			my $stddev = stddev(@fakes);
			my $nonseqpairs = $datahash{$site1}{$site2}[1];
			my $fakesnum = scalar @fakes;
			my $dat;
			my $pdb = "";
			## expect --nopdb
			my $pdb = $pdbdists{$site1}{$site2} unless $nopdb;
			if ($datahash{$site1}{$site2}[1] == 0){
				$dat = "NA";
				print OUT "$site1 $site2 $nonseqpairs $dat $fdrmean $stddev $fakesnum NA NA NA $pdb $fdrsum\n";
				next;
			}
			else { $dat = $datahash{$site1}{$site2}[0]/$datahash{$site1}{$site2}[1];}
			if ($fakesnum == 0 || $stddev == 0){
				print OUT "$site1 $site2 $nonseqpairs $dat $fdrmean $stddev $fakesnum NA NA NA $pdb $fdrsum\n";
				next;
			}
			my $zscore = ($dat-$fdrmean)/$stddev;
			my $lpval; #distance is greater that expected
			my $upval;
			foreach my $fdist (@fakes){
				if ($fdist>=$dat){$lpval++;}
				if ($fdist<=$dat){$upval++;}
			}
			my $lpval = $lpval/(scalar @fakes);
			my $upval = $upval/(scalar @fakes);
			
			print OUT "$site1 $site2 $nonseqpairs $dat $fdrmean $stddev $fakesnum $zscore $lpval $upval $pdb $fdrsum\n";
		}
		#}
	#}
	close OUT;
}

close LOG;
close ERRLOG;

sub hash_folder{
	my $fldr = shift;
	my $prt = shift;
	print LOG "Reading fakes from folder ".$fldr."\n";
	my @files = dirfiles($fldr);
	if (scalar @files == 0){die "No files found in $fldr!\n";}
	my @fdrfiles = grep {/^${prt}_[0-9]+_all_nonsequential/} @files;
	#print Dumper (\%datahash);
	my %fdrhash;
	#my @fdrdetailed;
	my $fcount = 0;
	foreach my $fdrfile(@fdrfiles){
		my $path = File::Spec->catfile($fldr, $fdrfile);
		tie *FH, 'File::ReadBackwards',  $path or die "Can't open $path: $!\n";
		my $lastline = <FH>;
		unless ($lastline =~ /^#/ ){ 
			print "File $path not ready yet. Skipping\n";
			untie *FH;
			next;
		}
		
		open FDR, "<$path" or die "Cannot open fdr file $path";
		my $header = <FDR>;
		while (<FDR>){
			my @splitter = split(/\s+/);
			last if (scalar @splitter < 9);
			if (! exists $fdrhash{$splitter[0]}{$splitter[1]}){
				$fdrhash{$splitter[0]}{$splitter[1]} = ();
			}
			push @{$fdrhash{$splitter[0]}{$splitter[1]}}, $splitter[8];
			#$fdrdetailed[$fcount]{$splitter[0]}{$splitter[1]} = [$splitter[6],$splitter[7]];
		}
		$fcount++;
	}
	close FDR;
	print LOG "Finished reading fakes from folder ".$fldr."\n";
	#return (\%fdrhash, \@fdrdetailed);
	return (\%fdrhash);
}


sub dirfiles{
        my $dirname = shift;
        opendir(DH, $dirname);
        my @files = readdir(DH);
        closedir(DH);
        return @files;
}

sub dirfolders{
	my $dirname = shift;
	opendir(DH, $dirname);
	my @files = readdir(DH);
	my @dirs = grep {-d "$dirname/$_" && ! /^\.{1,2}$/} @files;
	closedir(DH);
	return @dirs;
}
