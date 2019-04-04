#!/usr/bin/perl
package DistanceFinder;

use strict;
use Bio::Phylo::IO;
use Cwd qw(abs_path cwd getcwd);
use File::Path qw(make_path remove_tree);
use lib getcwd(); # adds working directory to @INC
use Data::Dumper;
use IPC::System::Simple qw(capture);
use vars qw(@ISA @EXPORT $VERSION);
use Exporter;
$VERSION = 1.00; # Or higher
@ISA = qw(Exporter);
@EXPORT = qw(get_mrcn calc_true_patristic_distance node_distance); # Symbols to autoexport (:DEFAULT tag)


	# sub set_mutmap_distance_matrix {
 		# my $mutmap = shift;
 		# my $prot = $mutmap->{static_protein};
 		# my $file = File::Spec->catfile($mutmap->{static_input_base}, $prot."_distance_matrix.csv");
 		# if (! (-e $file)){
 				# print "Preparing distance matrix..\n";
 				# my $logs = capture ('Rscript Distances.R --treefile '. $mutmap->{static_treefile}.' --output '.$file);
 				# print $logs."\n";
 		# }
	 	# open CSV, "<$file" or die "Cannot open file $file\n";
	 	# my $header = <CSV>;
	 	# $header =~ s/[\s\n\r\t]+$//s;
	 	# my @nodelables = split(',', $header);
	 	# while(<CSV>){
				# $_ =~ s/[\s\n\r\t]+$//s;
	 			# my @dists = split(',', $_);
	 			# my $node = $dists[0];
	 			# for (my $i = 1; $i < scalar @dists; $i++){
	 				# $mutmap->{static_distance_hash}{$node}{$nodelables[$i]} = $dists[$i];
	 			# }
	 	# }
	 	# close CSV;
 		
 	# }
	
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
