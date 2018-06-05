#!/usr/bin/perl

package Mutmap;

use strict;
use Bio::Phylo::IO;
use Class::Struct;
use compare qw(nsyn_substitutions syn_substitutions nsyn_substitutions_codons is_neighbour_changing);
use Bio::SeqIO;
use Cwd qw(abs_path cwd getcwd);
use File::Path qw(make_path remove_tree);
use lib getcwd(); # adds working directory to @INC
use Data::Dumper;
use IPC::System::Simple qw(capture);

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
# read .fasta into a hash
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

 	sub set_distance_matrix {
 		my $self = shift;
 		my $prot = $self->{static_protein};
 		my $file = File::Spec->catfile($self->{static_input_base}, $prot."_distance_matrix.csv");
 		if (! (-e $file)){
 				print "Preparing distance matrix..\n";
 				my $logs = capture ('Rscript Distances.R --treefile '. $self->{static_treefile}.' --output '.$file);
 				print $logs."\n";
 		}
	 	open CSV, "<$file" or die "Cannot open file $file\n";
	 	my $header = <CSV>;
	 	$header =~ s/[\s\n\r\t]+$//s;
	 	my @nodelables = split(',', $header);
	 	while(<CSV>){
				$_ =~ s/[\s\n\r\t]+$//s;
	 			my @dists = split(',', $_);
	 			my $node = $dists[0];
	 			for (my $i = 1; $i < scalar @dists; $i++){
	 				$self->{static_distance_hash}{$node}{$nodelables[$i]} = $dists[$i];
	 			}
	 	}
	 	close CSV;
 		
 	}
	
		# R matrix distance instead of calc_true_patristic_distance since 2018
  #  sub node_distance {
	#	my $mutmap = shift;
   # 	my $node = shift;
	#	my $other_node = shift;
	#	print (scalar keys %{$mutmap->{static_distance_hash}{$node->get_name}});
	#	print "\t";
	#	print (scalar keys %{$mutmap->{static_distance_hash}{$other_node->get_name}});
	#	print "\n";
    #	return $mutmap->{static_distance_hash}{$node->get_name}{$other_node->get_name};
   # }
	
	 sub node_distance {
	 	my $mutmap = shift;
    	my $node = shift;
		my $other_node = shift;
    	if  (exists $mutmap ->{static_distance_hash}{$node}{$other_node}){
    		return $mutmap ->{static_distance_hash}{$node}{$other_node};
    	}
    	else {
    		my $dist = calc_true_patristic_distance($node, $other_node);
    		$mutmap ->{static_distance_hash}{$node}{$other_node} = $dist;
    		return $dist;
    	}
    	
    }
	
	## the only difference from calc_patristic_distance is that it uses get_mrcn instead of get_mrca
## If you give it two sequential nodes, it returns the distance between them 

 sub calc_true_patristic_distance {
        my ( $node, $other_node ) = @_;
        my $patristic_distance = 0;
        my $mrca    = get_mrcn($node, $other_node);
        my $mrca_id = $mrca->get_id;
        while ( $node->get_id != $mrca_id ) {
            my $branch_length = $node->get_branch_length;
            if ( defined $branch_length ) {
                $patristic_distance += $branch_length;
            }
            $node = $node->get_parent;
        }
        while ( $other_node and $other_node->get_id != $mrca_id ) {
            my $branch_length = $other_node->get_branch_length;
            if ( defined $branch_length ) {
                $patristic_distance += $branch_length;
            }
            $other_node = $other_node->get_parent;
        }
        return $patristic_distance;
    }
	
	## unlike original bio::phylo get_mrca, returns the node n1 closest to the root, if you give it two sequential nodes n1 and n2
## (in that case get_mrca returns the youngest ancestor of n1)

sub get_mrcn {
        my ( $node, $other_node ) = @_;
        if ( $node->get_id == $other_node->get_id ) {
            return $node;
        }
        my $self_anc  = $node->get_ancestors;
		unshift @{$self_anc}, $node;
        my $other_anc = $other_node->get_ancestors;
		unshift @{$other_anc}, $other_node;
        for my $i ( 0 .. $#{$self_anc} ) {
            my $self_anc_id = $self_anc->[$i]->get_id;
            for my $j ( 0 .. $#{$other_anc} ) {
                if ( $self_anc_id == $other_anc->[$j]->get_id ) {
                    return $self_anc->[$i];
                }
            }
        }

        return $self_anc->[-1];
    }

1;
