#!/usr/bin/perl
package DistanceMap;

use strict;
use File::Spec;
use Cwd qw(abs_path cwd getcwd);
use File::Path qw(make_path remove_tree);
use lib getcwd(); # adds working directory to @INC
use Getopt::ArgvFile;
use Getopt::Long;
use Data::Dumper;
use DistanceFinder qw(get_mrcn calc_true_patristic_distance node_distance mock_distance);
use Parsers qw(parse_tree parse_fasta read_xpar parse_site_pair_distances);
#use compare qw(mynew);


$| = 1;

sub get_output_base {
	my $args = shift;
	my $output_base = File::Spec->catdir(getcwd(), "output", $args->{input}, $args->{bigtag}, $args->{state});
	if ($args->{fdr}){$output_base = File::Spec->catdir($output_base, "fdr");}
	if ($args->{phenotype}){$output_base = File::Spec->catdir($output_base, $args->{phenotype});}
	return $output_base;
}

sub get_input_base {
	my $args = shift;
	my $input_base = File::Spec->catdir(getcwd(), "data", $args->{input});
	if ($args->{fdr}){$input_base = File::Spec->catdir($input_base, "fdr", $args->{protein});}
	if ($args->{phenotype}){$input_base = File::Spec->catdir($input_base, $args->{phenotype});}
	return $input_base;
}

sub get_site2pheno_filepath {
	my $args = shift;
	my $filepath = File::Spec->catdir(getcwd(), "data", $args->{input}, "min_RR_SR.site2pheno");
	return $filepath;
}


sub new {
		my ($class, $args) = @_;	
		my $output_base = get_output_base($args);
		make_path($output_base);
		my $input_base = get_input_base($args);
		my $treefile = File::Spec->catfile($input_base, $args->{protein}.".l.r.newick");
		my $static_tree = parse_tree($treefile) or die "No tree at $treefile";
		my $xparfile;
		if($args->{fdr}){$xparfile = File::Spec->catfile($input_base, $args->{fakenum}.".xpar");}
		else{$xparfile = File::Spec->catfile($input_base, $args->{protein}.".xpar");}

		my $self;
		$self = {
				static_output_base => $output_base,
				static_input_base => $input_base,
				static_protein => $args->{protein},
				static_tree => $static_tree,
				static_treefile => $treefile,
				static_state  => $args->{state},
				static_fdr => $args->{fdr},
				static_phenotype => $args->{phenotype},
				static_fakenum => $args->{fakenum},
		};
		
		my %static_hash_of_nodes;	
		my @nodes = $static_tree -> get_nodes;
		foreach my $node(@nodes){
			my $name = $node ->get_name();
			$self->{static_hash_of_nodes}{$name} = \$node;
		}

		my @mutmaps = read_xpar($xparfile, $self->{static_hash_of_nodes}, $args->{state});
		$self->{static_subs_on_node} = $mutmaps[0];
		$self->{static_nodes_with_sub} = $mutmaps[1];
		$self->{static_length} = scalar keys $self->{static_nodes_with_sub};
		bless $self, $class;
		return $self;
}


sub global_pair_distances {
	my $self = shift;
	my @sites = @{$_[0]};
	# print "and now sites are ";
	# print Dumper \@sites;
	my $pair_counter = 0;
	my $distance_counter = 0;
	#print Dumper ($self->{static_nodes_with_sub});
	#print Dumper ($self->{static_hash_of_nodes});
	for (my $find = 0; $find < scalar @sites; $find++){
		my $fsite = $sites[$find];
		if (exists($self->{static_nodes_with_sub}{$fsite})){
			my @fnodes = @{$self->{static_nodes_with_sub}{$fsite}};
			for(my $sind = $find+1; $sind < scalar @sites; $sind++){
				my $ssite = $sites[$sind];
				if (exists($self->{static_nodes_with_sub}{$ssite})){
					my @snodes = @{$self->{static_nodes_with_sub}{$ssite}};
					if (scalar @fnodes > 0 and scalar @snodes > 0){
						# print ("site1 $fsite site2 $ssite\n");
						# print (scalar @fnodes);
						# print ("\n");
						# print (scalar @snodes);
						foreach my $fnode (@fnodes){
							foreach my $snode (@snodes){
								$pair_counter++;
								$distance_counter += node_distance($self, ${$fnode}, ${$snode});
							}
						}
					}
				}
			}
		}
	}
	return ($distance_counter, $pair_counter);
}

sub site_pair_all_distances {
	my $self = shift;
	my $fsite = shift;
	my $ssite = shift;
	my $pair_counter = 0;
	my $distance_counter = 0;
	
	if (exists($self->{static_nodes_with_sub}{$fsite}) & exists($self->{static_nodes_with_sub}{$ssite})){
		my @fnodes = @{$self->{static_nodes_with_sub}{$fsite}};
		my @snodes = @{$self->{static_nodes_with_sub}{$ssite}};
		# print ("site1 $fsite site2 $ssite\n");
		# print (scalar @fnodes);
		# print ("\n");
		# print (scalar @snodes);
		$pair_counter = (scalar @fnodes)*(scalar @snodes);
		foreach my $fnode (@fnodes){
			foreach my $snode (@snodes){
				$distance_counter += node_distance($self, ${$fnode}, ${$snode});
			}
		}
	}
	else {print "warning! static_nodes_with_sub does not contain $fsite or $ssite\n";}

	return ($distance_counter, $pair_counter);
}

sub all_nonsequential_distances {
	my $self = shift;
	my $pairs = shift;
	my $protein = $self->{static_protein};
	my $input_base = $self->{static_input_base};
	my $spfile = File::Spec->catfile($input_base, $protein.".site_pairs");

	my $statfile;
	my $output;
	print ("test self\n");
	print $self->{static_fdr}." ".!($self->{static_fdr})." ".$self->{static_fakenum}."\n";
	if($self->{static_fdr} != 1){
		print ("not an fdr simulation\n");
		$statfile = File::Spec->catfile($input_base, $protein.".stat.ident");
		$output = File::Spec->catfile($self->{static_output_base}, $protein."_all_nonsequential");
	}
	else{
		print ("fdr simulation\n");
		$statfile = File::Spec->catfile($input_base, $self->{static_fakenum}.".stat.ident");
		$output = File::Spec->catfile($self->{static_output_base},  $protein."_".$self->{static_fakenum}."_all_nonsequential");
	}
	
	print "statfile $statfile spfile $spfile output $output\n";
	open my $out, '>', $output or die "cannot open $output for writing\n";
	print $out "site1 site2 totaldist totalpairs seqdist seqpairs nonseqdist nonseqpairs mean_nonseqdist\n";
	my %seqhash = parse_site_pair_distances($statfile, $spfile);
	if ($pairs){
		foreach my $pair(@{$pairs}){
			my @seq = @{$seqhash{$pair->[0]}{$pair->[1]}};
			$self->print_diff($pair->[0], $pair->[1], $out, \@seq);
		}
	}
	else {
		foreach my $fsite(sort (keys %seqhash)){
			foreach my $ssite(sort (keys %{$seqhash{$fsite}})){
				my @seq = @{$seqhash{$fsite}{$ssite}};
				$self->print_diff($fsite, $ssite, $out, \@seq);
			}
		}
	}
	print $out "## Finished\n";
	print $out "## output_base ".$self->{static_output_base};
	close $out;
}

sub print_diff {
			my ($self, $fsite, $ssite, $out, $seq) = @_;
			my ($distance_counter, $pair_counter) = $self->site_pair_all_distances($fsite, $ssite);
			my $nonseqdist = $distance_counter-$seq->[0];
			my $nonseqpairs = $pair_counter-$seq->[1];
			my $meannonseq;
			if ($nonseqpairs > 0) {$meannonseq = $nonseqdist/$nonseqpairs;}
			else {$meannonseq = "NA";}
			print $out "$fsite $ssite $distance_counter $pair_counter ".$seq->[0]." ".$seq->[1]." ".$nonseqdist." ".$nonseqpairs." ".$meannonseq."\n";
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