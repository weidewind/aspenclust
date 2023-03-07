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
my $numfakes;
my $numchunks;

GetOptions (	
		'protein=s' => \$protein,
		'state=s' => \$state,
		'bigtag=s' => \$bigtag,
		'input=s' => \$input,
		'sites' => \$sites,
		'groups' => \$groups,
		'nopdb' => \$nopdb,
		'pairsfile=s' => \$pairsfile,
		'numfakes=s' => \$numfakes,
		'numchunks=s' => \$numchunks,
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


#my %fdrdetailed;
my @phenotypes = (0);
@phenotypes = dirfolders($fdrfolder);
if (scalar @phenotypes == 0){
	die "No phenotype folders found in $fdrfolder!\n";
}

my $pheno_file = DistanceMap::get_site2pheno_filepath($args);
my ($s2p, $phenotypes) = parse_site2pheno($pheno_file);
my @pairs = parse_pairs($pairsfile);
my $sampler = Sampler->new(\@pairs, \@phenotypes, $s2p);


if ($sites){
	my $fdroutput = File::Spec->catfile($outfolder, $protein."_dist.fdr_results");
	my %pdbdists = collect_pdbdists($pdbfile, $spfile) unless $nopdb;
	open OUT, ">$fdroutput" or die "cannot open $fdroutput for writing\n";
	print OUT "site1 site2 nonseqpairs meandist fdrmean fdrstddev fakesnum zscore lpval upval pdb fdrsum\n";
	
	my %files;
	foreach my $ph (@phenotypes){
		my $fldr = File::Spec->catfile($fdrfolder, $ph);
		my @fdrfiles = grep {/^${protein}_[0-9]+_all_nonsequential/} dirfiles($fldr);
		my @filepaths;
		foreach my $ff (@fdrfiles){
			push @filepaths, File::Spec->catfile($fldr, $ff);
		}
		$files{$ph} = \@filepaths;
	}
	
	print LOG "Numfakes ".$numfakes." numchunks ".$numchunks." \n";
	my @indices = 0 .. $numfakes-1;
	print LOG Dumper \@indices;
	my @chunks = array_to_chunks(\@indices, $numchunks);
	print LOG "here are the chunks: \n";
	print LOG Dumper \@chunks;
	my %chunks_stat;
	foreach my $file_indices_array (@chunks){
		my %fdrhash;
		foreach my $ph (@phenotypes){
			print LOG "Reading fakes for phenotype ".$ph." and files \n";
			## collect filepaths for this chunk of this phenotype's fakes
			## weird solution, but won't rewrite
			my @ph_files = @{$files{$ph}}[@{$file_indices_array}];
			# foreach my $i (@{$file_indices_array}){
				# push @ph_files, $files{$ph}->[$i];
			# }
			print LOG Dumper \@ph_files;
			my $fdrhash_ph = hash_files(\@ph_files);
			print LOG "Finished reading\n";
			$fdrhash{$ph} = $fdrhash_ph;
		}
		# print LOG "fdrhash looks like this\n";
		# print LOG Dumper \%fdrhash;
		for (my $i = 0; $i < scalar @pairs; $i++){
			my $site1 = $pairs[$i]->[0];
			my $site2 = $pairs[$i]->[1];
			
			## expect that it is the same for all phenotypes

			
			my @fakes;
			for (my $f = 0; $f < scalar @{$file_indices_array}; $f++){
				my $ph = $sampler->sample_pheno_for_pair_index($i);
				#my $ph = 1;
				push @fakes, $fdrhash{$ph}{$site1}{$site2}[$f];
			}
			@fakes = grep {/[0-9\.]+/} @fakes; # because it can be NA, and we are not interested in NAs
			
			# print LOG join(" ", @fakes[0..5])."\n";
			my $fdrmean = mean(@fakes);
			my $fdrsum = sum(@fakes);
			my $stddev = stddev(@fakes);
			my $nfakes = scalar @fakes;
			my $lpval; #distance is greater that expected
			my $upval;
			## if no relevant muts in data, mean distance for data is not defined and neither are lpval and upval
			unless ($datahash{$site1}{$site2}[1] == 0){
				my $dat = $datahash{$site1}{$site2}[0]/$datahash{$site1}{$site2}[1];
				foreach my $fdist (@fakes){
					if ($fdist>=$dat){$lpval++;}
					if ($fdist<=$dat){$upval++;}
				}
			}
			unless ($nfakes == 0){
				my $nfakes_old = exists $chunks_stat{$site1}{$site2}{"nfakes"} ? $chunks_stat{$site1}{$site2}{"nfakes"} : 0;
				my $stddev_old = exists $chunks_stat{$site1}{$site2}{"stddev"} ? $chunks_stat{$site1}{$site2}{"stddev"} : 0;
				my $fdrmean_old = exists $chunks_stat{$site1}{$site2}{"mean"} ? $chunks_stat{$site1}{$site2}{"mean"} : 0;
				
				my $nfakes_combined = $chunks_stat{$site1}{$site2}{"nfakes"} + $nfakes;
				my $fdrsum_combined = $chunks_stat{$site1}{$site2}{"fdrsum"} + $fdrsum;
				my $fdrmean_combined = $fdrsum_combined/$nfakes_combined;
				#my $stddev_combined = sqrt(($nfakes_old*($stddev_old**2 + ($fdrmean_old - $fdrmean_combined)**2) + $nfakes*($stddev**2 + ($fdrmean - $fdrmean_combined)**2))/$nfakes_combined);
				my $stddev_combined;
				if ($nfakes_old == 0){$stddev_combined = $stddev;}
				else {$stddev_combined = sqrt((($nfakes_old)*$stddev_old**2 + ($nfakes)*$stddev**2)/($nfakes_old + $nfakes));}
				$chunks_stat{$site1}{$site2}{"nfakes"} = $nfakes_combined;
				$chunks_stat{$site1}{$site2}{"fdrsum"} = $fdrsum_combined;
				$chunks_stat{$site1}{$site2}{"stddev"} = $stddev_combined;
				$chunks_stat{$site1}{$site2}{"mean"} = $fdrmean_combined;
				unless ($datahash{$site1}{$site2}[1] == 0){
					$chunks_stat{$site1}{$site2}{"lpval"} += $lpval;
					$chunks_stat{$site1}{$site2}{"upval"} += $upval;
				}
			}

		}
		# print LOG "after adding new files, combined stat for sites {100073}{168499} is \n";
		# print LOG Dumper $chunks_stat{100073}{168499};
	}

	for (my $i = 0; $i < scalar @pairs; $i++){
		my $site1 = $pairs[$i]->[0];
		my $site2 = $pairs[$i]->[1];
		
		my $nonseqpairs = $datahash{$site1}{$site2}[1];
		my $dat;
		my $pdb = "";
		## expect --nopdb
		my $pdb = $pdbdists{$site1}{$site2} unless $nopdb;
		
		## todo rename variables!
		my $fdrmean = $chunks_stat{$site1}{$site2}{"mean"};
		my $stddev = $chunks_stat{$site1}{$site2}{"stddev"};
		my $fakesnum = $chunks_stat{$site1}{$site2}{"nfakes"};
		my $fdrsum = $chunks_stat{$site1}{$site2}{"fdrsum"};
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
		my $lpval = $chunks_stat{$site1}{$site2}{"lpval"}/$fakesnum;
		my $upval = $chunks_stat{$site1}{$site2}{"upval"}/$fakesnum;
		
		print OUT "$site1 $site2 $nonseqpairs $dat $fdrmean $stddev $fakesnum $zscore $lpval $upval $pdb $fdrsum\n";
	}
	

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


sub hash_files{
	my @fdrfiles = @{$_[0]};
	my %fdrhash;
	my $fcount = 0;
	foreach my $path(@fdrfiles){
		tie *FH, 'File::ReadBackwards',  $path or die "Can't open $path: $!\n";
		my $lastline = <FH>;
		unless ($lastline =~ /^#/ ){ 
			untie *FH;
			die "File $path not ready yet. Skipping\n";
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
	return (\%fdrhash);
}

sub array_to_chunks{
	my @array = @{$_[0]};
	my $n = $_[1];
	my $total = scalar @array;
	my $chunklen = $total/$n;
	my $res = $total%$n;
	my @chunks;
	while (@array){
		my $tchunklen = $chunklen;
		if ($res > 0){
			$tchunklen += 1;
			$res -= 1;
		}
		my @ch = splice @array, 0, $tchunklen;
		push @chunks, \@ch;
	}
	return @chunks;
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
