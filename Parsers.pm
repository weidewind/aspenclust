#!/usr/bin/perl

package Parsers;

use strict;
use vars qw(@ISA @EXPORT $VERSION);
use Exporter;
$VERSION = 1.00; # Or higher
@ISA = qw(Exporter);
@EXPORT = qw(parse_fasta parse_likelihoods parse_tree parse_likelihoods_as_fasta parse_fasta_as_likelihoods) ; # Symbols to autoexport (:DEFAULT tag)

use Bio::Phylo::IO;
use Bio::SeqIO;
use Data::Dumper;

# read .fasta into a hash
# returns 
# $nodeseqs{$nodename} = $string
sub parse_fasta {
	my $nodeseqs_file = shift;
	my %nodeseqs;
	my $seqio = Bio::SeqIO->new(-file => $nodeseqs_file, -format => "fasta");
	my $length;
	while ( my $seqobj = $seqio->next_seq ) {
		my $trimmed_id = (split(/\//, $seqobj->display_id))[0];
    	$nodeseqs{ $trimmed_id } = $seqobj->seq;
    	if (!$length){
    		$length = $seqobj->length();
    	}
	}
	return (\%nodeseqs, $length);
}

# returns 
# $nodeseqs{$nodename}[$ind] = (0,0,0.8,0.2)
# $keystring = "ACGT"

sub parse_likelihoods {
	my $file = shift;
	my %nodeseqs;
	open FILE, "<$file" or die "Cannot open $file\n"; 
	my $s = <FILE>;
	my @str = split(/\s+/,$s);
	my $nodename = substr($str[0],1)."_".$str[1];
	my $keystring = $str[2];
	
	while($s){
		while ($s && substr($s,0,1) ne ">"){
			my @str = split(/\s+/,$s);
			my $ind = shift @str;
			$nodeseqs{$nodename}[$ind-1] = \@str;
			$s = <FILE>;
		}
		# print ($nodename."\n");
		# print (Dumper($nodeseqs{$nodename}[1])."\n");
		last unless $s;
		my @str = split(/\s+/,$s);
		$nodename = substr($str[0],1)."_".$str[1];
		$s = <FILE>;

	}
	close FILE;
	return (\%nodeseqs, $keystring);
}

sub parse_likelihoods_as_fasta {
	my $file = shift;
	my ($nodeseqsprob, $keystring) = parse_likelihoods($file);
	my %nodeseqs;
	foreach my $node (keys %{$nodeseqsprob}){
		my $string = "";
		for (my $ind = 0; $ind < scalar @{$nodeseqsprob->{$node}}; $ind++){
				my $letter = substr($keystring, index_of_max($nodeseqsprob->{$node}[$ind]), 1);
				$string = $string.$letter;
		}
		$nodeseqs{$node} = $string;
	}
	return %nodeseqs;
}


# returns 
# $nodeseqs{$nodename}[$i] = (0,0,0.8,0.2)
# $keystring = "ACGT"
sub parse_fasta_as_likelihoods {
	my $nodeseqs_file = shift;
	my %nodeseqs;
	my $keystring = "ACGT";
	my $seqio = Bio::SeqIO->new(-file => $nodeseqs_file, -format => "fasta");
	my $length;
	while ( my $seqobj = $seqio->next_seq ) {
		my $trimmed_id = (split(/\//, $seqobj->display_id))[0];
		for my $i (0..$seqobj->length()-1){
			my $letter = substr($seqobj->seq, $i,1);
			my $probs = [0,0,0,0];
			my $ind = index($keystring, $letter);
			$probs->[$ind] = 1;
			$nodeseqs{ $trimmed_id } [$i] = $probs;
		}
	}
	return (\%nodeseqs, $keystring);
}

sub index_of_max {
	my @array = @{$_[0]};
	my $index = 0;
	my $maxval = $array[$index];
	for my $i ( 0 .. $#array )
	{
		$maxval < $array[$i] and
		$index = $i
			if  $maxval < $array[$i];
	}
	return $index;
}

sub parse_tree {
		my $tree_file = $_[0];
		open TREE, "< $tree_file" or die "Cannot open file ".$tree_file."\n";
		# get a newick string from some source
		my $tree_string;
		my $t;
		while($t = <TREE>){
			$t =~ s/\n$//;
			$tree_string = $tree_string.$t;
		}
		 close TREE;

		# Call class method parse from Bio::Phylo::IO
		# note: newick parser returns 'Bio::Phylo::Forest'
		# Call ->first to retrieve the first tree of the forest.
		my $tree = Bio::Phylo::IO->parse(
		  -string => $tree_string,
		  -format => 'newick'
		)->first;

		return $tree;	
} 


1;
