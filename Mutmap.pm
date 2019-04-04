#!/usr/bin/perl

package Mutmap;

use strict;
use Bio::Phylo::IO;
#use DistanceFinder qw(get_mrcn calc_true_patristic_distance node_distance);
use Class::Struct;
use compare qw(nsyn_substitutions syn_substitutions nsyn_substitutions_codons is_neighbour_changing);
use Bio::SeqIO;
use Cwd qw(abs_path cwd getcwd);
use File::Path qw(make_path remove_tree);
use lib getcwd(); # adds working directory to @INC
use Data::Dumper;
use IPC::System::Simple qw(capture);
use Parsers qw(parse_fasta parse_tree);

sub new {
	my ($class, $args) = @_;	
	my $output_base = File::Spec->catdir(getcwd(), "output", $args->{bigtag}, $args->{state});
	make_path($output_base);
	my $input_base = File::Spec->catdir(getcwd(), "data");
	
	my $treefile = File::Spec->catfile($input_base, $args->{protein}.".l.r.newick");
	my $fastafile = File::Spec->catfile($input_base, $args->{protein}.".all.fa");
	
	my $tree = parse_tree($treefile );
	my @arr = parse_fasta($fastafile);
	my %fasta = %{$arr[0]};
	my $length = $arr[1];
	my @mutmaps;
	if($args->{state} eq "syn"){
		 @mutmaps = synmutmap($tree, \%fasta);
	} 
	elsif($args->{state} eq "nsyn"){
		 @mutmaps = codonmutmap($tree, \%fasta);
	} 
	else {
		die "only syn or nsyn can be used as the second argument; unknown $args->{state} was used instead";
	}
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
		stataic_distance_hash => \%distance_hash,
	};
	bless $self, $class;
	#$self->set_distance_matrix();
	return $self;
}

sub aa_length {
	my $self = shift;
	return ($self->{static_length})/3;
}


sub mutmap {

	my $tree = $_[0];
	my %nodeseqs = %{$_[1]};
	my %subs_on_node;
	my %nodes_with_sub;
	
	my @nodes = $tree -> get_nodes;
	foreach my $node(@nodes){
		if ($node->is_root()) {next;}
		my $name = $node ->get_name();
		#my @nsyn = nsyn_substitutions($nodeseqs{$node->get_ancestors()->[0]->get_name()},
		#							  $nodeseqs{$name});
		my %nsyn = compare::nsyn_substitutions($nodeseqs{$node->get_ancestors()->[0]->get_name()},
									  $nodeseqs{$name});					  
		$subs_on_node{$name}=\%nsyn;
		for	my $site_index(keys %nsyn){
			if (! exists $nodes_with_sub{$site_index}){
				$nodes_with_sub{$site_index} = ();
			}
			push (@{$nodes_with_sub{$site_index}}, \$node);
		}	
	}

	return (\%subs_on_node, \%nodes_with_sub);
};

sub codonmutmap {

	my $tree = $_[0];
	my %nodeseqs = %{$_[1]};
	my %subs_on_node;
	my %nodes_with_sub;
	
	my @nodes = $tree -> get_nodes;
	foreach my $node(@nodes){
		if ($node->is_root()) {next;}
		my $name = $node ->get_name();
		#my @nsyn = nsyn_substitutions($nodeseqs{$node->get_ancestors()->[0]->get_name()},
		#							  $nodeseqs{$name});
		my %nsyn = compare::nsyn_substitutions_codons($nodeseqs{$node->get_ancestors()->[0]->get_name()},
									  $nodeseqs{$name});					  
		$subs_on_node{$name}=\%nsyn;
		for	my $site_index(keys %nsyn){
			if (! exists $nodes_with_sub{$site_index}){
				$nodes_with_sub{$site_index} = ();
			}
			push (@{$nodes_with_sub{$site_index}}, \$node);
		}	
	}

	return (\%subs_on_node, \%nodes_with_sub);
};



sub synmutmap {

	my $tree = $_[0];
	my %nodeseqs = %{$_[1]};
	my %subs_on_node;
	my %nodes_with_sub;
	
	my @nodes = $tree -> get_nodes;
	foreach my $node(@nodes){
		if ($node->is_root()) {next;}
		my $name = $node ->get_name();
		my %nsyn = compare::syn_substitutions($nodeseqs{$node->get_ancestors()->[0]->get_name()},
									  $nodeseqs{$name});					  
		$subs_on_node{$name}=\%nsyn;
		for	my $site_index(keys %nsyn){
			if (! exists $nodes_with_sub{$site_index}){
				$nodes_with_sub{$site_index} = ();
			}
			push (@{$nodes_with_sub{$site_index}}, \$node);
		}	
	}

	return (\%subs_on_node, \%nodes_with_sub);
}

1; 