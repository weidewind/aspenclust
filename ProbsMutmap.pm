#!/usr/bin/perl

package ProbsMutmap;

use strict;
use Bio::Phylo::IO;
use Class::Struct;
use compareProbs qw(new substitutions);
use Bio::SeqIO;
use Cwd qw(abs_path cwd getcwd);
use File::Path qw(make_path remove_tree);
use lib getcwd(); # adds working directory to @INC
use Data::Dumper;
use IPC::System::Simple qw(capture);
use Parsers qw(parse_likelihoods parse_fasta_as_likelihoods parse_tree parse_fasta);
#use DistanceFinder qw(get_mrcn calc_true_patristic_distance node_distance);

sub likelihood_tag {
	my $l = shift;
	return "likelihoods" if $l;
	return "mostlikely" if !$l;
}

sub new {
	my $class = shift;
	my $args = shift;	
	my $output_base = File::Spec->catdir(getcwd(), "output", $args->{bigtag}, $args->{state}, likelihood_tag($args->{likelihood}));
	make_path($output_base);
	my $input_base = File::Spec->catdir(getcwd(), "data");
	
	my $treefile = File::Spec->catfile($input_base, $args->{protein}.".ancestor.likelihoods.withinternals.newick");
	my $fastafile = File::Spec->catfile($input_base, $args->{protein}.".nointernals.fa");
	my $internalseqs;
	my $keystring;
	if ($args->{likelihood}){
		my $probsfile = File::Spec->catfile($input_base, $args->{protein}.".ancestor.likelihoods");
		($internalseqs, $keystring) = parse_likelihoods($probsfile); # $nodeseqs{$nodename}[$ind] = (0,0,0.8,0.2)
	}
	else {
		my $probsfile = File::Spec->catfile($input_base, $args->{protein}.".ancestor.fa");
		($internalseqs, $keystring) = parse_fasta_as_likelihoods($probsfile); # $nodeseqs{$nodename}[$ind] = (0,0,0.8,0.2)
	}
	my $tree = parse_tree($treefile);
	my ($fasta, $length) = parse_fasta($fastafile);
	my @mutmaps = probmutmap($tree, $fasta, $internalseqs, $args->{state}, $keystring);
	my %distance_hash;
	my $self;
	$self = {
		static_output_base => $output_base,
		static_input_base => $input_base,
		static_protein => $args->{protein},
		static_length => $length,
		static_tree => $tree,
		static_treefile => $treefile,
		static_state => $args->{state},
		static_subs_on_node => $mutmaps[0], 
		static_nodes_with_sub => $mutmaps[1], 
		static_distance_hash => \%distance_hash,
	};
	bless $self, $class;
	#$self->set_distance_matrix();
	return $self;
}

sub graphpath{
	my $self = shift;
	my $ind = shift;
	my $filename = $self->{static_protein}."_".$ind.".graph";
	my $path = File::Spec->catdir($self->{static_output_base}, "graphs");
	make_path($path);
	my $path = File::Spec->catfile($path, $filename);
}

sub probmutmap {
	my $tree = $_[0];
	my %fasta = %{$_[1]}; #$fasta{$nodename} = $string
	my %internalseqs = %{$_[2]}; # $internalseqs{$nodename}[$ind] = (0,0,0.8,0.2)
	my $state = $_[3];
	my $keystring = $_[4];
	my %subs_on_node;
	my %nodes_with_sub;
	my $comparator = compareProbs->new($keystring),
	my @nodes = $tree -> get_nodes;
	foreach my $node(@nodes){
		if ($node->is_root()) {next;}
		my $name = $node ->get_name();
		my $seq;
		if ($node ->is_terminal()){
			$seq = $fasta{$name};
		}
		else{
			$seq = $internalseqs{$name};
		}
		my %subs = $comparator->substitutions($internalseqs{$node->get_ancestors()->[0]->get_name()},$seq, $state);		
		$subs_on_node{$name}=\%subs; # $subs{$ind} = (Substitution,Substitution,..)
		for	my $site_index(keys %subs){
			if (! exists $nodes_with_sub{$site_index}){
				$nodes_with_sub{$site_index} = ();
			}
			push (@{$nodes_with_sub{$site_index}}, \$node);
		}	
	}
	return (\%subs_on_node, \%nodes_with_sub);
}

sub set_same_ancestor_subs{
	my $self = shift;
	return if (exists $self->{static_same_ancestor_subs}); 
	my %subs_on_node = %{$self->{static_subs_on_node}};
	my %hash;
	
	foreach my $nname (keys %subs_on_node){
		foreach my $ind (keys $subs_on_node{$nname}){
			foreach my $sub (@{$subs_on_node{$nname}->{$ind}}){
				my $ancestor = $sub->{"Substitution::ancestral_allele"};
				push @{$hash{$ind}{$ancestor}}, $sub;
			}	
		}
	}
	
	$self->{static_same_ancestor_subs} = \%hash;
}

sub set_all_distances_probs{
	my $self = shift;
	my $hash = shift;
	my $site_index = shift;
	$self->{static_all_distances_probs}{$site_index} = $hash;
}

sub pathFinder {
	my $self = shift;
	my $args = shift;
	unless ($args) {return $self ->{static_output_base};}
	my $path = File::Spec->catdir($self ->{static_output_base}, $args->{norm});
	make_path($path);
	return $path;
}

1;