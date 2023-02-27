#!/usr/bin/perl
package Parsers;

use strict;
use vars qw(@ISA @EXPORT $VERSION);
use Exporter;
use List::Util qw(min max);
use compare;
$VERSION = 1.00; # Or higher
@ISA = qw(Exporter);
@EXPORT = qw(parse_fasta read_xpar parse_site_pair_distances parse_likelihoods parse_likelihoods_csv parse_tree parse_likelihoods_as_fasta parse_fasta_as_likelihoods parse_splits) ; # Symbols to autoexport (:DEFAULT tag)

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

sub read_xpar{
	my $xpar_file = $_[0];
	my $hash_of_nodes = $_[1];
	my $state = $_[2];
	open XPAR, "<$xpar_file" or die "Cannot open xpar file ".$xpar_file."\n";
	my %subs_on_node;
	my %nodes_with_sub;
	my $header = <XPAR>;
	while(<XPAR>){
			my @splitter = split(/\s/);
			my $nodname = $splitter[0];
			my @inds;
			if ($state eq "nsyn"){@inds = split(';',$splitter[4]);}
			elsif ($state eq "syn"){@inds = split(';',$splitter[3]);}
			else {die "Unknown state $state";}
			my %substs;
			foreach my $ind(@inds){
				my $ind = $ind =~ m/([0-9]+)/;
				$ind = $1;
				my $p=Supstitution->new();
				$p->position($ind);
				$p->ancestral_allele("ATG");
				$p->derived_allele("ATG");
				$substs{$ind} = $p;
				if (! exists $nodes_with_sub{$ind}){
					$nodes_with_sub{$ind} = ();
				}
	#			print "\n nodndame ".$_."\n";
	#			print "\nREF 1 ".ref($static_hash_of_nodes{$nodname})."\n";
	#			print "\nREF 2 ".ref(${$static_hash_of_nodes{$nodname}})."\n";
				push (@{$nodes_with_sub{$ind}}, \${$hash_of_nodes->{$nodname}}); #вытащить из дерева по имени
			}
			$subs_on_node{$nodname} = \%substs;

	}
	close XPAR;
	return (\%subs_on_node, \%nodes_with_sub);
}

sub parse_site_pair_distances {
	my $stat_file = shift;
	my $site_pairs_file = shift;
	my %hash;
	open STAT, "<$stat_file" or die "cannot open statfile $stat_file\n";
	open PAIRS, "<$site_pairs_file" or die "cannot open pairsfile $site_pairs_file \n";
	my $header = <PAIRS>;
	while (<PAIRS>){
		my ($fsite, $ssite) = split(/\s+/);
		my ($dist, $count) = split(/\s+/, <STAT>);
		my $f = min($fsite, $ssite);
		my $s = max($fsite, $ssite);
		## debugging
		# $hash{$f}{$s} = ();
		# $hash{$f}{$s}[0] += $dist;
		# $hash{$f}{$s}[1] += $count;
		if (!exists $hash{$f}{$s}){
			$hash{$f}{$s} = ();
		}
		$hash{$f}{$s}[0] += $dist;
		$hash{$f}{$s}[1] += $count;

	}
	close PAIRS;
	close STAT;
	return %hash;
}


# returns 
# $nodeseqs{$nodename}[$ind] = (0,0,0.8,0.2)
# $keystring = "ACGT"

sub parse_likelihoods {
	my $file = shift;
	my $csv = shift;
	if ($csv){ 
		return parse_likelihoods_csv($file);
	}
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


# parse likelihoods from mega cmd-line version 
sub parse_likelihoods_csv {
	my $file = shift;
	my %nodeseqs;
	open FILE, "<$file" or die "Cannot open $file\n"; 
	my $s = <FILE>;
	my @str = split(/N/,$s);
	shift @str;
	my @nodenames = map {"N".$_} @str;
	my $nodec = scalar @nodenames;
	my $keystring = "ACGT";
	my $c = 0;
	my %keyhash = map {$_ => $c++} split(//, $keystring); 
	print Dumper \%keyhash;
	while(<FILE>){
		my @values = split(/\s+/,$_);
		my $head = shift @values;
		my $valc = scalar @values;
		if ($valc != $nodec){die "The number of probabilities ( $valc )  does not match the number of nodes ( $nodec ) in file $file !";}
		my ($ind, $letter) = split(",", $head);
		for (my $i = 0; $i < scalar @nodenames; $i++){
			$nodeseqs{$nodenames[$i]}[$ind-1][$keyhash{$letter}] = sprintf("%.3f", $values[$i]);
		}
	}
	close FILE;
	return (\%nodeseqs, $keystring);
}

sub parse_likelihoods_as_fasta {
	my $file = shift;
	my $csv = shift;
	my ($nodeseqsprob, $keystring) = parse_likelihoods($file, $csv);
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
		open TREE, "< $tree_file" or die "Cannot open tree file ".$tree_file."\n";
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

sub parse_splits {
	my $splitsfile = shift;
	my %node_partition;
	
	open SPLITS, "<$splitsfile" or die "Cannot open $splitsfile : $! \n";
	while (<SPLITS>){
		my ($nname, $type) = split(/\s+/);
		$node_partition{$nname} = $type;
	}
	return %node_partition;
}




1;
