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
$Statistics::Basic::IPRES = 6;
use Storable 'dclone';
use File::ReadBackwards;




my $protein = "toy";
my $state = "nsyn";
my $bigtag = "test";
my $input;
my $sites;
my $groups;
my $nopdb;

GetOptions (	
		'protein=s' => \$protein,
		'state=s' => \$state,
		'bigtag=s' => \$bigtag,
		'input=s' => \$input,
		'sites' => \$sites,
		'groups' => \$groups,
		'nopdb' => \$nopdb,
);

$| = 1;

my $args = {protein => $protein, state => $state, bigtag => $bigtag, input => $input};
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
print "Finished reading the data\n";


my @files = dirfiles($fdrfolder);
my @fdrfiles = grep {/^${protein}_[0-9]+_all_nonsequential/} @files;
#print Dumper (\%datahash);
my %fdrhash;
my @fdrdetailed;
my $fcount = 0;
foreach my $fdrfile(@fdrfiles){
	my $path = File::Spec->catfile($fdrfolder, $fdrfile);
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
		$fdrdetailed[$fcount]{$splitter[0]}{$splitter[1]} = [$splitter[6],$splitter[7]];
	}
	$fcount++;
}
close FDR;
print "Finished reading fakes\n";

my $pdbfile = File::Spec->catfile(DistanceMap::get_input_base($args), "pdb", $protein.".pdb_dists") unless $nopdb;
my $spfile = File::Spec->catfile(DistanceMap::get_input_base($args), $protein.".site_pairs");

if ($sites){
	my $fdroutput = File::Spec->catfile($outfolder, $protein."_dist.fdr_results");
	my %pdbdists = collect_pdbdists($pdbfile, $spfile) unless $nopdb;
	open OUT, ">$fdroutput" or die "cannot open $fdroutput for writing\n";
	print OUT "site1 site2 nonseqpairs meandist fdrmean fdrstddev fakesnum zscore lpval upval pdb\n";
	foreach my $site1 (sort (keys %fdrhash)){
		foreach my $site2 (sort (keys %{$fdrhash{$site1}})){
			my @fakes = @{$fdrhash{$site1}{$site2}};
			@fakes = grep {/[0-9\.]+/} @fakes;
			my $fakesnum = scalar @fakes;
			my $fdrmean = mean(@fakes);
			my $stddev = stddev(@fakes);
			my $nonseqpairs = $datahash{$site1}{$site2}[1];
			my $dat;
			my $pdb = "";
			my $pdb = $pdbdists{$site1}{$site2} unless $nopdb;
			if ($datahash{$site1}{$site2}[1] == 0){
				$dat = "NA";
				print OUT "$site1 $site2 $nonseqpairs $dat $fdrmean $stddev $fakesnum NA NA NA $pdb\n";
				next;
			}
			else { $dat = $datahash{$site1}{$site2}[0]/$datahash{$site1}{$site2}[1];}
			if ($fakesnum == 0 || $stddev == 0){
				print OUT "$site1 $site2 $nonseqpairs $dat $fdrmean $stddev $fakesnum NA NA NA $pdb\n";
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
			
			print OUT "$site1 $site2 $nonseqpairs $dat $fdrmean $stddev $fakesnum $zscore $lpval $upval $pdb\n";
		}
	}
	close OUT;
}

if ($groups){
		#fdr for neg and pos epistatically connected pairs
		my $posfile = File::Spec->catfile(DistanceMap::get_input_base($args), "epistat", $protein.".supertree.pos_cor2pcor.l90.n0.pos.tab");
		my $negfile = File::Spec->catfile(DistanceMap::get_input_base($args), "epistat", $protein.".supertree.all.neg.l0500.tab");
		my $fdroutput_groups = File::Spec->catfile($outfolder, $protein."_dist_groups.fdr_results");
		my $log = File::Spec->catfile($outfolder, $protein."_logs");
		open LOG, ">$log" or die "cannot open $log";
		open GOUT, ">$fdroutput_groups" or die "cannot open $fdroutput_groups for writing\n";
		print GOUT "group nonseqpairs meandist fdrmean fdrstddev zscore lpval upval\n";

		print LOG "\n$protein positive ";
		my @pospairs = collect_pairs($posfile);
		my @p = get_stats_for_group(\@pospairs,\@fdrdetailed, \%datahash);
		print GOUT "pos ".join(" ", @p)."\n";
		print LOG "\n$protein negative ";
		my @negpairs = collect_pairs($negfile);
		my @n = get_stats_for_group(\@negpairs,\@fdrdetailed, \%datahash);
		print GOUT "neg ".join(" ", @n)."\n";
		print LOG "\n$protein nonepistatic ";
		my @complpairs = collect_compl_pairs([\@pospairs, \@negpairs], \%datahash);
		my @empty = get_stats_for_group(\@complpairs,\@fdrdetailed, \%datahash);
		print GOUT "empty ".join(" ", @empty)."\n";

		my @closepairs = collect_close_pairs($pdbfile, $spfile);
		#print Dumper \@closepairs;
		my @c = get_stats_for_group(\@closepairs,\@fdrdetailed, \%datahash);
		print GOUT "pdb_close ".join(" ", @c)."\n";
		my @notclose = collect_compl_pairs([\@closepairs], \%datahash);
		my @nc = get_stats_for_group(\@notclose,\@fdrdetailed, \%datahash);
		print GOUT "not_pdb_close ".join(" ", @nc)."\n";
		close GOUT;
		close LOG;
}


sub get_stats_for_group{
	my @pospairs = @{$_[0]};
	my @fdrdetailed = @{$_[1]};
	my %datahash = %{$_[2]};
	my ($posdat, $posnonseqpairs) = collect_group_data(\@pospairs, \%datahash);
	my @posfakes = collect_group_fakes(\@pospairs, \@fdrdetailed);
	print LOG "in this group we have ".scalar @pospairs." pairs of sites\n";
	print LOG "fake stats for this group is as follows \n";
	print LOG Dumper(\@posfakes);
	my $fdrmean = mean(@posfakes);
	my $stddev = stddev(@posfakes);

	my $zscore = ($posdat-$fdrmean)/$stddev;
	my $lpval; #distance is greater that expected
	my $upval;
	foreach my $fdist (@posfakes){
		if ($fdist>=$posdat){$lpval++;}
		if ($fdist<=$posdat){$upval++;}
	}
	my $lpval = $lpval/(scalar @posfakes);
	my $upval = $upval/(scalar @posfakes);

	return ($posnonseqpairs, $posdat, $fdrmean, $stddev, $zscore, $lpval, $upval);
}



sub collect_group_fakes{
	my @pospairs = @{$_[0]};
	my @fdrdetailed = @{$_[1]};
	my @posfakes;
	foreach my $f (@fdrdetailed){
		my $posdist;
		my $poscount;
		foreach my $pos (@pospairs){
			my $inf = $f->{$pos->[0]}{$pos->[1]};
			$posdist += $inf->[0];
			$poscount += $inf->[1];
		}
		push @posfakes, $posdist/$poscount;
	}
	return @posfakes;
}

sub collect_group_data{
	my @pospairs = @{$_[0]};
	my %datahash = %{$_[1]};
	my $datdist;
	my $datcount;
	foreach my $pos (@pospairs){
		$datdist += $datahash{$pos->[0]}{$pos->[1]}[0];
		$datcount += $datahash{$pos->[0]}{$pos->[1]}[1];
	}
	return ($datdist/$datcount, $datcount);
}

# sub collect_empty_data{
	# my @groups = @{$_[0]};
	# my %datahash = %{$_[1]};
	# my $datdist;
	# my $datcount;
	# foreach my $site1 (keys %datahash){
		# foreach my $site2 (keys %{$datahash{$site1}}){
			# $datdist += $datahash{$site1 ->[0]}{$site2->[1]}[0];
			# $datcount += $datahash{$site1 ->[0]}{$site2->[1]}[1];
		# }
	# }
	# my $grdist;
	# my $grcount;
	# foreach my $group (@groups){
		# foreach my $pos (@{$group}){
			# $grdist += $datahash{$pos->[0]}{$pos->[1]}[0];
			# $grcount += $datahash{$pos->[0]}{$pos->[1]}[1];
		# }
	# }
	
	# return (($datdist-$grdist)/($datcount-$grcount), ($datcount-$grcount));
# }


sub collect_pairs {
	my $filename = shift;
	my @pairs;
	open F, "<$filename" or die "cannot open file $filename for collecting ";
	while (<F>){
		my ($site1, $site2) = split(/\s+/);
		push @pairs, [$site1, $site2];
	}
	close F;
	return @pairs;
}

sub collect_compl_pairs {
	my @sets_of_pairs = @{$_[0]};
	my %datahash = %{$_[1]};
	my @complpairs;
	my $dhash = dclone \%datahash;
	foreach my $set (@sets_of_pairs){
		foreach my $p (@{$set}){
			delete($dhash->{$p->[0]}{$p->[1]});
		}
	}
	foreach my $site1 (keys %{$dhash}){
		foreach my $site2 (keys %{$dhash->{$site1}}){
			push @complpairs, [$site1, $site2];
		}
	}
	return @complpairs;
}

sub collect_close_pairs {
	my $pdbfilename = shift;
	my $spfilename = shift;
	my @pairs;
	open PDB, "<$pdbfilename" or die "cannot open pdb file $pdbfilename for collecting ";
	open SP, "<$spfilename" or die "cannot open sp file $spfilename for collecting ";
	my $header = <SP>;
	while (<PDB>){
		my ($issame, $pdbdist) = split(/\s+/);
		my ($site1, $site2) = split(/\s+/, <SP>);
		if ($issame eq "same" & $pdbdist <5){
			push @pairs, [$site1, $site2];
		}
	}
	close PDB;
	close SP;
	return @pairs;
}

sub collect_pdbdists {
	my $pdbfilename = shift;
	my $spfilename = shift;
	my %dists;
	open PDB, "<$pdbfilename" or die "cannot open pdb file $pdbfilename for collecting ";
	open SP, "<$spfilename" or die "cannot open sp file $spfilename for collecting ";
	my $header = <SP>;
	while (<PDB>){
		my ($issame, $pdbdist) = split(/\s+/);
		my ($site1, $site2) = split(/\s+/, <SP>);
		$dists{$site1}{$site2} = $pdbdist;
	}
	close PDB;
	close SP;
	return %dists;
}

sub dirfiles{
        my $dirname = shift;
        opendir(DH, $dirname);
        my @files = readdir(DH);
        closedir(DH);
        if (scalar @files == 0){
                print "No files found in $dirname!\n";
        }
        return @files;
}
