#!/usr/bin/perl
use strict;
use File::Spec;

my $dirname = "output/parenttest/nsyn/likelihoods/weightnorm";

opendir(DH, $dirname);
my @files = readdir(DH);
closedir(DH);
if (scalar @files == 0){
	die "No files found in $dirname!\n";
}

foreach my $filename(sort @files){
	my $filepath =  File::Spec->catfile($dirname,$filename);
	next if (-d $filepath);
	next unless ($filename =~ /([a-zA-Z0-9]*)_n?syn_groups_(.*)shuffler_(.*)_(me[a-z]+)/);
	print $filepath."\n";
	my $prot = $1;
	my $shufflertype = $2;
	my $groupname = $3;
	my $stattype = $4;
	my $resfile = File::Spec->catfile($dirname, $prot."_groups_results_".$shufflertype);
	
	
	#print header
	unless (-e $resfile){ 
		open GROUPS, ">>$resfile" or die "Cannot open $resfile: $!\n";
		print GROUPS "protein,group,stat_type,group_count,group_same,group_diff,compl_count,compl_same,compl_diff,diffdiff,";
		if ($shufflertype eq "label"){
			print GROUPS "group_pvalue,enrich_pvalue,depl_pvalue\n";
		}
		elsif ($shufflertype eq "site"){
			print GROUPS "group_attraction_pvalue,gattr_and_crep_pvalue,group_repulsion_pvalue,grep_and_cattr_pvalue,group_attraction_enrichment,group_repulsion_enrichment\n";
		}
		close GROUPS;
	}
	open GROUPS, ">>$resfile" or die "Cannot open $resfile: $!\n";
	print GROUPS $prot.",".$groupname.",".$stattype.",";
	open IN, "<$filepath" or die "Cannot open $filepath $!";
	while (<IN>){
		if ($_ =~ /^group\s+[0-9]+/ || $_ =~ /^complement/){
			my @spl = split(/\s+/);
			print GROUPS $spl[1].",".$spl[2].",".$spl[3].",";
		}
		if ($_ =~ /^diffdiff/){
			my @spl = split(/\s+/);
			print GROUPS $spl[1].",";
		}
		if ($_ =~ /.*pvalue.*/){
			my @spl = split(/\s+/);
			print $_."\n";
			print GROUPS join(",", grep {/[0-9\.]+/} @spl);
			print GROUPS ",";
		}
	}
	close IN;
	print GROUPS "\n";
	close GROUPS;
}	

