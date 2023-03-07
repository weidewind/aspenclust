#!/usr/bin/perl
package Sampler;

use strict;
use File::Spec;
use Cwd qw(abs_path cwd getcwd);
use File::Path qw(make_path remove_tree);
use lib getcwd(); # adds working directory to @INC
use Getopt::ArgvFile;
use Getopt::Long;
use Data::Dumper;
use Parsers qw(parse_site2pheno);
use DistanceMap;
#use compare qw(mynew);


$| = 1;


sub new {
	my $class = shift;
	my @pairs = @{$_[0]};
	my @phenotypes = @{$_[1]};
	my @site2pheno = @{$_[2]};
	
	my $p0 = 0.01;
	my $npairs = scalar @pairs;
	my $npheno = scalar @phenotypes;
	
	my @pairs2phen_probs;
	for(my $i=0;$i<$npairs;$i++){
		my @phen_probs;
		my $norm=0;
		for(my $j=0;$j<$npheno;$j++){
			# print "pheno ".$phenotypes[$j].", pval for site1 ".$site2pheno[$pairs[$i]->[0]]->{$j}."pval for site1 ".$site2pheno[$pairs[$i]->[1]]->{$j};
			my $val=$site2pheno[$pairs[$i]->[0]]->{$phenotypes[$j]}*$site2pheno[$pairs[$i]->[1]]->{$phenotypes[$j]};
			$norm+=$val;
			push @phen_probs,$val;
		}
		# print "probs for pair ".$i.": ".$pairs[$i]->[0]." and ".$pairs[$i]->[1]."\n"; 
		# print Dumper \@phen_probs;
		my $n=0;
		for(my $j=0;$j<$npheno;$j++){
			$phen_probs[$j]/=$norm if $norm>0;
			$phen_probs[$j]=$p0 unless $phen_probs[$j]>0;
			$phen_probs[$j]=1/$phen_probs[$j];
			$n+=$phen_probs[$j];
		}
		for(my $j=0;$j<$npheno;$j++){
			$phen_probs[$j]/=$n;
		}
		for(my $j=1;$j<$npheno;$j++){
			$phen_probs[$j]+=$phen_probs[$j-1];
		}
		push @pairs2phen_probs,[@phen_probs];
	}
	my $self;
	$self = {
			pairs2phen_probs => \@pairs2phen_probs,
			pairs => \@pairs,
			phenotypes => \@phenotypes
	};
	bless $self, $class;
	return $self;
}


sub sample_pheno_for_pair_index{
	my $self = shift;
	my $pair_index = shift;
	my @probs=@{$self->{pairs2phen_probs}->[$pair_index]};
	die "\nThe function expects a vector of cumulative probabilities as an argument!" unless(sprintf("%.6f",$probs[-1])==1);
	my $k=0;
	if(scalar @probs > 1){
		my $smpl=rand;
		for($k;$k<@probs;$k++){
			last if $smpl<$probs[$k];
		}
	}
	return $self->{phenotypes}->[$k];
}

sub print_probs_for_pair_index{
	my $self = shift;
	my $pair_index = shift;
	my @probs=@{$self->{pairs2phen_probs}->[$pair_index]};
	my @pairs = @{$self->{pairs}};
	print Dumper \@pairs;
	print "pair index ".$pair_index." site1 ".$pairs[$pair_index]->[0]." site2 ".$pairs[$pair_index]->[1]."\n";
	print "phenotype prob\n";
	for (my $i = 0; $i<scalar @probs; $i++){
		print $self->{phenotypes}->[$i]." ".$probs[$i]."\n";
	}
}
