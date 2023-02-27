#!/usr/bin/perl
## This class provides a map of mutations (node_name -> array of sites that mutated on the branch),  
## given a tree and sequences of all nodes (including internal ones)

## 


 # magic - ctrlf for magic constants
package Mutmapnolim;

use strict;
use Bio::Phylo::IO;


use Bio::Tools::CodonTable;
use TreeUtils::Phylo::FigTree;
use Bio::SeqIO;
use Data::Dumper;
use Bit::Vector;
use Try::Tiny;
use List::Util qw(sum min max);
use Const::Fast;
use Switch;
use List::Util qw/shuffle/; 
use Statistics::Basic qw(:all);
use Statistics::TTest;
use Statistics::Descriptive;
use Storable qw(store retrieve lock_retrieve);
use IPC::System::Simple qw(capture);
use Class::Struct;
use IO::Handle;
use Cwd qw(abs_path cwd getcwd);
use File::Spec;
use Clone 'clone';
use Sub::Identify ':all';
use File::Path qw(make_path remove_tree);
use autodie;
use Math::Round;
use Storable qw(dclone);
use AgeingStat;
use MeanStat;
use MedianStat;
use BPStat;
use Weeds;
use Textbits qw(concat cleave iterationFiles);
#use DnaUtilities::observation_vector qw(make_observation_vector shuffle_obsv);
use observation_vector qw(make_observation_vector shuffle_obsv);
#use DnaUtilities::compare qw(nsyn_substitutions syn_substitutions nsyn_substitutions_codons is_neighbour_changing);
use compare qw(nsyn_substitutions syn_substitutions nsyn_substitutions_codons get_synmuts new);
#use TreeUtils::Phylo::PhyloUtils qw(remove_zero_branches);
use shuffle_muts_on_tree qw(shuffle_mutations_on_tree prepare_shuffler StripConstrains Shuffler);
use shuffle_muts_on_tree_exp qw(shuffle_mutations_on_tree Constrains);
use PhyloUtils qw(remove_zero_branches);
use Groups;
use Memusage;
use Codeversion;

$| = 1;

	
	
	
	sub set_alignment_length {
		my $self = shift;
		$self->{static_alignment_length} = $_[0]; 
	}
		
	
	sub maxpath_tag{
		my $subtract_maxpath = $_[0];
		my $tag;
		if (defined $subtract_maxpath){
			if ($subtract_maxpath eq "y" || $subtract_maxpath eq "yes" || $subtract_maxpath == 1 ){
				$tag = "maxpath_subtracted";
			}
			elsif ($subtract_maxpath eq "n" || $subtract_maxpath eq "no" || $subtract_maxpath == 0 ) {
				$tag = "maxpath_not_subtracted";
			}
			else {die "Invalid subtract_maxpath: $subtract_maxpath";}
		}
		else {$tag = '';}
		
		return $tag;
	}


	sub state_tag {
		my $state = $_[0];
		my $tag;
		if ($state eq "s" || $state eq "syn") { $tag = "syn";}
		elsif ($state eq "n" || $state eq "nsyn") { $tag = "nsyn";}
		else {die "Unknown state $state; expected syn or nsyn";}
		return $tag;
	}
	
	sub syn_tag {
		my $syn = $_[0];
		my $tag;
		if ($syn == 1){
			$tag = "syn";
		}
		else {
			$tag = "nsyn";
		}
		return $tag;
	}

	sub temp_tag {
			return "unreadable";
	}
	
	sub stoppers_tag {
		my $skip_stoppers = shift;
		my $tag;
		if (defined $skip_stoppers){
			if ($skip_stoppers eq "y" || $skip_stoppers eq "yes" || $skip_stoppers == 1 ){
				$tag = "skip_stoppers";
			}
			elsif ($skip_stoppers eq "n" || $skip_stoppers eq "no" || $skip_stoppers == 0 ) {
				$tag = "with_stoppers";
			}
			else {die "Invalid skip_stoppers: $skip_stoppers";}
		}
		else {$tag = '';}
		
		return $tag;
	}
	
	sub neighbour_tag {
		my $no_neighbour_changing = shift;
		my $tag;
		if (defined $no_neighbour_changing){
			if ($no_neighbour_changing eq "y" || $no_neighbour_changing eq "yes" || $no_neighbour_changing == 1 ){
				$tag = "no_neighbour_changing";
			}
			elsif ($no_neighbour_changing eq "n" || $no_neighbour_changing eq "no" || $no_neighbour_changing == 0 ) {
				$tag = "with_neighbour_changing";
			}
			else {die "Invalid no_neighbour_changing: $no_neighbour_changing";}
		}
		else {$tag = '';}
		
		return $tag;
	}
	
	sub leaves_tag {
		my $no_leaves = shift;
		my $tag;
		if (defined $no_leaves){
			if ($no_leaves eq "y" || $no_leaves eq "yes" || $no_leaves == 1 ){
				$tag = "no_leaves";
			}
			elsif ($no_leaves eq "n" || $no_leaves eq "no" || $no_leaves == 0 ) {
				$tag = "with_leaves";
			}
			else {die "Invalid no_leaves: $no_leaves";}
		}
		else {$tag = '';}
		
		return $tag;
	}
	
	sub syn_lengths_tag {
		my $syn_lengths = shift;
		if ($syn_lengths){
			return "syn_lengths";
		}
		else {
			return '';
		}
	}
	

	sub printFooter {
		my $self = shift;
		my $outputStream = shift;
		print $outputStream "## protein ".$self->{static_protein}."\n";
		print $outputStream "## subtract_tallest ".$self->{static_subtract_tallest}."\n";
		print $outputStream "## state ".$self->{static_state}."\n";
		print $outputStream "## mutnum_control ".$self->{static_mutnum_control}."\n";
		print $outputStream "## omit neighbour-changing mutations (only for 'reversals', ancestor n-ch muts are not skipped. Only valid for syn state)? 1 if true ".$self->{static_no_neighbour_changing}."\n";
		print $outputStream "## omit mutations on terminal branches? 1 if true ".$self->{static_no_leaves}."\n";
		print $outputStream "## include halves of branches after mutations? 1 if true ".$self->{static_include_tips}."\n";
		print $outputStream "## skip stoppers (neighbour-changing background mutations are ignored (i.e. not used for subtree pruning))? Only valid for nsyn state. 1 if true ".$self->{static_skip_stoppers}."\n";
		print $outputStream "## compute branch lengths in syn substs? 1 if true ".$self->{static_syn_lengths}."\n";
		print $outputStream "## output_base ".$self->{static_output_base}."\n";
		if ($self->{realdata}){
			print $outputStream "## realdata restriction ".get_realdata_restriction($self->{realdata})."\n";
		}
		print $outputStream "## code version hash ".Codeversion->get_version()."\n";
	}
	
	sub pathFinder {
		my $args = shift;	
		my $output_base = File::Spec->catdir(getcwd(), "output", $args->{bigdatatag}, $args->{bigtag}, state_tag($args->{state}), maxpath_tag($args->{subtract_tallest}), stoppers_tag($args->{skip_stoppers}), neighbour_tag($args->{no_neighbour_changing}), syn_lengths_tag($args->{syn_lengths}), leaves_tag($args->{no_leaves})); 
		return $output_base;

	}
	
	sub createCodeversionFile {
		my $self = shift;
		my $script_name = shift;
		my $version_file = File::Spec->catfile($self->{static_output_base}, $script_name."_codeversion");
		open FILE, ">$version_file" or die "Cannot open $version_file: $!\n";
		print FILE (Codeversion::get_version());
		close FILE;
	}
	
	sub realdata_exists {
		my $args = shift;
		my $output_base = pathFinder($args);
		my $realdatapath = File::Spec->catfile($output_base, $args->{protein}."_".state_tag($args->{state})."_realdata");
		print "Checking if $realdatapath exists..";
		if (-f $realdatapath) {print "yes!\n";}
		else {print "no.\n";}
		return (-f $realdatapath);
	}
	
	# not to be confused with get_realdata_restriction, which needs realdata as an argument
	sub check_realdata_restriction{
		my $args = shift;
		my $output_base = pathFinder($args);
		my $realdatapath = File::Spec->catfile($output_base, $args->{protein}."_".state_tag($args->{state})."_realdata");
		my $realdata = lock_retrieve ($realdatapath) or die "Cannot retrieve ".$realdatapath;
		return get_realdata_restriction($realdata);
	}
	
	sub dataFinder {
		my $args = shift;	
		my $input_base = File::Spec->catdir(getcwd(), "data", $args->{bigdatatag}, syn_lengths_tag($args->{syn_lengths}));
		return $input_base;
	}

	sub set_tag {
		my $self = shift;
		my $tag = shift;
		my $output_subfolder = File::Spec->catdir($self->{static_output_base}, $tag);
		$self->{static_output_subfolder} = $output_subfolder;
		make_path($output_subfolder);
		print($output_subfolder);
		make_path(File::Spec->catdir($output_subfolder, temp_tag()));
	}
	
	sub get_tree {
		my $self = shift;
		return $self->{static_tree};
	}
	
	sub get_treefile {
		my $self = shift;
		return $self->{static_treefile};
	}
	
	sub new {
		my ($class, $args) = @_;	
		my $output_base = pathFinder ($args);
		my $input_base = dataFinder ($args);
		my $treefile = File::Spec->catfile($input_base, $args->{protein}.".l.r.newick");
		if ($args->{syn_lengths} && !(-e $treefile)) {print_tree_with_syn_lengths_static($args);}
		my $static_tree;
		if ($args->{tree_object}){$static_tree = $args->{tree_object};}
		else {$static_tree = parse_tree($treefile)  or die "No tree at $treefile";}
		
		my $self;
		make_path($output_base);
		make_path(File::Spec->catdir($output_base, temp_tag()));
		
		if ($args->{fromfile}){
			my $realdatapath;
			my $realdataname = $args->{protein}."_".state_tag($args->{state})."_realdata";
			if ($args->{fake}){
				 $realdatapath = File::Spec->catfile($output_base, $args->{tag}, $realdataname);
			}
			else { $realdatapath = File::Spec->catfile($output_base, $realdataname); }
			print "Creating mutmap from realdata $realdatapath\n";
			my $realdata = lock_retrieve ($realdatapath) or die "Cannot retrieve ".$realdatapath;
			if ($args->{fake}){
				unlink $realdatapath or warn "Could not unlink $realdatapath: $!";
			}
			
			#foreach my $ind(1..300){
			#if ($realdata->{static_nodes_with_sub}{$ind}){
			#my $debugnum = scalar @{$realdata->{static_nodes_with_sub}{$ind}};
			#print "Early news from nnew: numnodes for $ind is $debugnum\n";
			#}
			#if ($realdata->{static_nodes_with_sub}{$ind} && scalar @{$realdata->{static_nodes_with_sub}{$ind}} > 0){# possibility to have at least 2 mutations after an ancestral one
			#foreach my $node(@{$realdata->{static_nodes_with_sub}{$ind}}){
		#		print "Early News from nnew: nnode_name ".$$node->get_name()."\n";
		#	}
		#	}
		#	}
			
			

				$self = { 
					static_output_base => $output_base,
					static_input_base => $input_base,
					static_protein => $args->{protein},
					static_subtract_tallest => $args->{subtract_tallest},
					static_tree => $static_tree,
					static_treefile => $treefile,
					static_state => $args->{state},
					static_no_neighbour_changing =>$realdata->{no_neighbour_changing}, 
					static_mutnum_control => $args->{mutnum_control}, 
					static_no_leaves =>$realdata->{no_leaves},
					static_include_tips =>$realdata->{include_tips},
					static_skip_stoppers =>$realdata->{skip_stoppers},
					static_syn_lengths  =>$realdata->{syn_lengths},
					static_alignment_length => $realdata->{alignment_length}, 
					static_hash_of_nodes => $realdata->{hash_of_nodes}, 
					#static_distance_hash => $realdata->{distance_hash},
					static_subs_on_node => $realdata->{static_subs_on_node}, # we never use these two when we produce new mutmappers from file (they are taken from observaton_vectors)
					obs_vectors => $realdata->{obs_vectors}, #added on 17.11.2016
					static_nodes_with_sub => $realdata->{static_nodes_with_sub}, #
					static_background_subs_on_node => $realdata->{bkg_subs_on_node},
					static_background_nodes_with_sub => $realdata->{bkg_nodes_with_sub},
					static_comparator => compare->new(),
					realdata => $realdata,
				};
			
			
		#	foreach my $ind(1..300){
		#	if ($self->{static_nodes_with_sub}{$ind}){
		#	my $debugnum = scalar @{$self->{static_nodes_with_sub}{$ind}};
		#	print "news from nnew: numnodes for $ind is $debugnum\n";
		#	}
		#	if ($self->{static_nodes_with_sub}{$ind} && scalar @{$self->{static_nodes_with_sub}{$ind}} > 0){# possibility to have at least 2 mutations after an ancestral one
		##	foreach my $node(@{$self->{static_nodes_with_sub}{$ind}}){
		#		print "News from nnew: nnode_name ".$$node->get_name()."\n";
		#	}
		#	}
		#	}
			
		}
		else {
			my $fastafile = File::Spec->catfile($input_base, $args->{protein}.".all.fa");
			if ($args->{syn_lengths} && !(-e $fastafile)){
				my $copyargs = {%{$args}};
				delete $copyargs->{syn_lengths} ;
				$fastafile = File::Spec->catfile(dataFinder($copyargs), $args->{protein}.".all.fa");
			}
			my @arr = parse_fasta($fastafile);
			my %fasta = %{$arr[0]};
			my $alignment_length = $arr[1];
			my $static_protein  = $args->{protein};
			my %static_fasta = %fasta;
			my %static_hash_of_nodes;	
			my @nodes = $static_tree -> get_nodes;
			
			my @mutmaps;
			my @bkg_mutmaps;
			if($args->{state} eq "syn"){
				@mutmaps = synmutmap($static_tree, \%fasta);
				@bkg_mutmaps = codonmutmap($static_tree, \%fasta);
			} 
			elsif($args->{state} eq "nsyn"){
				@mutmaps = codonmutmap($static_tree, \%fasta);
				@bkg_mutmaps = synmutmap($static_tree, \%fasta);
			} 
			else {
				die "only syn or nsyn can be used as the second argument; unknown ".$args->{state}." was used instead";
			}

			# static_hash_of_nodes is here now

			$self = {
				static_output_base => $output_base,
				static_input_base => $input_base,
				static_protein => $args->{protein},
				static_alignment_length => $alignment_length, 
				static_subtract_tallest => $args->{subtract_tallest},
				static_no_neighbour_changing =>  $args->{no_neighbour_changing},
				static_mutnum_control => $args->{mutnum_control}, 
				static_no_leaves =>$args->{no_leaves},
				static_include_tips =>$args->{include_tips},
				static_skip_stoppers =>$args->{skip_stoppers},
				static_syn_lengths  =>$args->{syn_lengths},
				static_tree => $static_tree,
				static_treefile => $treefile,
				static_fasta => { %static_fasta },
				static_state  => $args->{state},
				static_hash_of_nodes => { %static_hash_of_nodes },
				static_subs_on_node => $mutmaps[0],
				static_nodes_with_sub => $mutmaps[1],
				static_background_subs_on_node => $bkg_mutmaps[0],
				static_background_nodes_with_sub => $bkg_mutmaps[1],
				static_comparator => compare->new(),
			};

## static_hash_of_nodes been here
			foreach my $node(@nodes){
				#if ($node->is_root()) {next;}
				my $name = $node ->get_name();
				$self->{static_hash_of_nodes}{$name} = \$node;
			}
		}	
		
		bless $self, $class;
		$self->set_tag($args->{tag});
		return $self;
	}
	






## according to http://www.biomedcentral.com/1471-2148/10/253
		const my @n1_decreasing => ("AGG", "TCG", "GAT", "CGT", "ACC", "GCC", "CAG", "GGG", "GGC");

		const my @h1_decreasing => ("ACG", "TCA", "CTC", "GCG", "GCA", "CCG", "TGC", "GTG");

		const my @n2_decreasing => ("AAT", "CTC", "GAG", "TCT", "ACT", "TGT", "CCG", "GGC", "GAC", "AAA", "TCA");

		const my @h3_decreasing => ("CTG", "CGC", "CCT", "TGC",  "GAC", "AGG", "TAT", "AAG", "GGG", "CGG");

		const my @n1_increasing => ("AGA", "ACA", "GGA", "CAC", "TCA", "CTT", "CAA", "AGT");

		const my @h1_increasing  => ("ACA", "GCC", "CCT", "TGT", "TCC", "AGC");

		const my @n2_increasing => ("AAC", "TCC", "GAA", "GTT", "TGC", "GAT", "AAG", "GCC", "ACA");
		
		const my @h3_increasing => ("TTG", "AGA", "TGT", "GCC",  "CTA", "GAT", "TAC", "CCG", "GGA", "AAA", "CCC");

		const my @all_codons => ("TCA", "TCC", "TCG", "TCT", "TTC", "TTT", "TTA", "TTG", "TAC", "TAT", "TAA", "TAG", "TGC", "TGT", "TGA", "TGG", "CTA", "CTC", "CTG", "CTT", "CCA", "CAT", "CAA", "CAG", "CGA", "CGC", "CGG", "CGT", "ATA", "ATC", "ATT", "ATG", "ACA", "ACC", "ACG", "ACT", "AAC", "AAT", "AAA", "AAG", "AGC", "AGT", "AGA", "AGG", "CCC", "CCG", "CCT", "CAC", "GTA", "GTC", "GTG", "GTT", "GCA", "GCC", "GCG", "GCT", "GAC", "GAT", "GAA", "GAG", "GGA", "GGC", "GGG", "GGT");

		

## returns a hash: key - codon, value - -1, if it is decreasing over time,  
##                                       1, if it is increasing, 
##										 0 otherwise.
sub codon_evolution{
	my $protein = $_[0];
	my %hash;
	@hash{@all_codons} = 0;
	switch($protein){
		case "n1" {	@hash{@n1_decreasing} = -1; 
			        @hash{@n1_increasing} = 1; }
		case "n2" {	@hash{@n2_decreasing} = -1; 
					@hash{@n2_increasing} = 1; }
		case "h1" {	@hash{@h1_decreasing} = -1; 
					@hash{@h1_increasing} = 1; }
		case "h3" {	@hash{@h3_decreasing} = -1; 
					@hash{@h3_increasing} = 1; }
		}
	
	return %hash;
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
					open TREE, "<$tree_file" or die "Cannot open file ".$tree_file."\n";
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


sub codonmutmap {

	my $tree = $_[0];
	my %nodeseqs = %{$_[1]};
	my %subs_on_node;
	my %nodes_with_sub;
	
	my @nodes = $tree -> get_nodes;
	foreach my $node(@nodes){
		if ($node->is_root()) {next;}
		my $name = $node ->get_name();
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



sub print_tree_with_syn_lengths_static {
	my $args = shift;
	$args->{syn_lengths} = 0;
	my $mutmap = Mutmapnolim->new($args);
	$mutmap->print_tree_with_syn_lengths();
}

sub print_tree_with_syn_lengths {
	my $self = shift;
	if ($self->{static_syn_lengths}){
		print "Already printed\n";
		return;
	}
	my $tree = clone($self->{static_tree});
	my $map;
	if ($self->{static_state} eq "nsyn"){
		$map = $self->{static_background_subs_on_node};
	}
	else {
		$map = $self->{static_subs_on_node};
	}
	my $outputdir = File::Spec->catdir($self->{static_input_base}, syn_lengths_tag(1));
	make_path($outputdir);
	my $treefile = File::Spec->catfile($outputdir, $self->{static_protein}.".l.r.newick");
	print_syn_lenths_tree ($tree, $map, $treefile);
}


sub print_syn_lenths_tree{
	my $tree = shift;
	my $map = shift;
	my $output = shift;
	
	my @nodes =  $tree-> get_nodes;
	foreach my $node (@nodes){
		$node -> set_branch_length(scalar keys %{$map->{$node->get_name}});
	}
	my $string = $tree->get_root->to_newick('-nodelabels' => 1);
	open FILE, ">$output" or die $!;
	print FILE $string;
	close FILE;
}

sub print_tree_with_nsyn_lengths{
	my $self = shift;
	unless ($self->{static_state} eq "nsyn"){($self->{static_state} = "nsyn");}
	my $map = $self->{static_subs_on_node};
	my $tree = clone($self->{static_tree});
	my $treefile = File::Spec->catfile($self->{static_input_base}, $self->{static_protein}.".l.r.updated.newick");
	
	my @nodes =  $tree-> get_nodes;
	foreach my $node (@nodes){
		$node -> set_branch_length(scalar keys %{$map->{$node->get_name}});
	}
	my $string = $tree->get_root->to_newick('-nodelabels' => 1);
	open FILE, ">$treefile" or die $!;
	print FILE $string;
	close FILE;
}


sub print_tree_with_allsubst_lengths{
	my $self = shift;
	unless ($self->{static_state} eq "nsyn"){();}
	my $map = $self->{static_subs_on_node};
	my $backmap = $self->{static_background_subs_on_node};
	my $tree = clone($self->{static_tree});
	my $treefile = File::Spec->catfile($self->{static_input_base}, $self->{static_protein}.".l.r.updated.newick");
	
	my @nodes =  $tree-> get_nodes;
	foreach my $node (@nodes){
		my $nsyns = scalar keys %{$map->{$node->get_name}};
		my $syns = scalar keys %{$backmap->{$node->get_name}};
		print $nsyns."and".$syns."\n";
		$node -> set_branch_length($nsyns+$syns);
	}
	my $string = $tree->get_root->to_newick('-nodelabels' => 1);
	open FILE, ">$treefile" or die $!;
	print FILE $string;
	close FILE;
}

sub mylength {
	my $self = shift;	
	my $length = ($self->{static_alignment_length})/3;
	return $length;
}
# sets static_sorted_nodnames and static_sorted_sites, retruns incidence_hash
sub incidence_matrix {
	my $self = shift;	
	my %matrix;
	my $length = $self->mylength();
	my @sorted_sites;
	my @sorted_nodnames;
	
	# select nodes with at least one mutation of the corresponding type (syn or nsyn, depending on the mutmap state)
	my @nodes = $self -> {static_tree} -> get_nodes;
	foreach my $node(@nodes){
		my $name = $node ->get_name();
		if (scalar (keys %{$self->{static_subs_on_node}{$name}}) > 0){
			push @sorted_nodnames, $name;
		}
	}
	
	my %empty_nodes_hash = map { $_ => 0 } @sorted_nodnames;
	my %incidence_hash;
	
	# select sites with at least 3 mutations of the corresponding type
	# upd - not sure if such sites should be excluded, changed minimum to 1
	foreach my $ind(1..$length){
		if ($self->{static_nodes_with_sub}{$ind}){
			my $debugnum = scalar @{$self->{static_nodes_with_sub}{$ind}};
			#print "news from incidence: numnodes for $ind is $debugnum\n";
		}
		if ($self->{static_nodes_with_sub}{$ind} && scalar @{$self->{static_nodes_with_sub}{$ind}} > 0){# possibility to have at least 2 mutations after an ancestral one
			my %site_incidence = %empty_nodes_hash;
			push @sorted_sites, $ind;
			#print " added $ind to sorted sites\n";
			foreach my $node(@{$self->{static_nodes_with_sub}{$ind}}){
				#print "News from incidence_matrix: nnode_name ".$$node->get_name()."\n";
				$site_incidence{$$node->get_name()} = 1;
			}
			$incidence_hash{$ind} = \%site_incidence;
		}
	}
	
	$self->{static_sorted_nodnames} = \@sorted_nodnames;
	$self->{static_sorted_sites} = \@sorted_sites;
	
	return %incidence_hash;
};

sub matrixPath {
	my $self = shift;
	my $statetag = state_tag($self->{static_state});
	my $matrix_file = File::Spec -> catfile($self->{static_output_base}, $self->{static_protein}."_".$statetag."_incidence_matrix");
	return $matrix_file;
}

sub xparPath {
	my $self = shift;
	my $xpar_file = File::Spec -> catfile($self->{static_output_base}, $self->{static_protein}.".xpar");
	return $xpar_file;
}

sub print_incidence_matrix {
	my $self = $_[0];
	my %incidence_hash = %{$_[1]};
	#my $path = $_[2];
	my $statetag = state_tag($self->{static_state});
	my $matrix_file = $self->matrixPath;
	my $sorted_sites_file = File::Spec -> catfile($self->{static_output_base}, $self->{static_protein}."_".$statetag."_sorted_sites");
	my $sorted_nodnames_file = File::Spec -> catfile($self->{static_output_base}, $self->{static_protein}."_".$statetag."_sorted_nodnames");
	
	# todo check if static_sorted_nodnames exists, if not - throw error and die (where did you take incidence_hash from?)
	unless ($self->{static_sorted_nodnames} && $self->{static_sorted_sites}) {die "There is no static_sorted_nodnames (or static_sorted_sites) in this mutmap. Where did you take incidence_hash from?\n";}
	
	open MATRIX, ">$matrix_file" or die "Cannot open file ".$matrix_file."\n";
	foreach my $nodname (@{$self->{static_sorted_nodnames}}){
		foreach my $ind (@{$self->{static_sorted_sites}}){
			print MATRIX $incidence_hash{$ind}->{$nodname};
		}
		print MATRIX "\n";
	}
	close MATRIX;
	
	open SSITES, ">$sorted_sites_file" or die "Cannot open file ".$sorted_sites_file."\n";
	foreach my $ind(@{$self->{static_sorted_sites}}){
		print SSITES $ind."\n";
	}
	close SSITES;
	
	open SNODES, ">$sorted_nodnames_file" or die "Cannot open file ".$sorted_nodnames_file."\n";
	foreach my $name(@{$self->{static_sorted_nodnames}}){
		print SNODES $name."\n";
	}
	close SNODES;
	
}


sub read_incidence_matrix {
	my $self = $_[0];
	my $matrix_file = $_[1];
	open MATRIX, "<$matrix_file" or die "Cannot open file ".$matrix_file."\n";
	my %subs_on_node;
	my %nodes_with_sub;
	my $line_index = 0;
	while(<MATRIX>){
		if (/^$/) {last;}
			my $nodname = $self->{static_sorted_nodnames}[$line_index];
			my @sites = split(',');
			my %substs;
			foreach my $s(@sites){
				my $ind = $self ->{static_sorted_sites}[$s-1];
#				print " $s is $ind\n";
				my $p=Substitution->new();
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
				push (@{$nodes_with_sub{$ind}}, \${$self ->{static_hash_of_nodes}{$nodname}}); #вытащить из дерева по имени
			}
			$subs_on_node{$nodname} = \%substs;
			$line_index++;

	}
	close MATRIX;
	return (\%subs_on_node, \%nodes_with_sub);
			
}

sub read_xpar{
	my $self = $_[0];
	my $xpar_file = $_[1];
	open XPAR, "<$xpar_file" or die "Cannot open file ".$xpar_file."\n";
	my %subs_on_node;
	my %nodes_with_sub;
	my $line_index = 0;
	my $header = <XPAR>;
	while(<XPAR>){
			my @splitter = split(/\s/);
			my $nodname = $splitter[0];
			my @nsyn_inds = split(';',$splitter[4]);
			my %substs;
			print "debug $nodname ";
			foreach my $ind(@nsyn_inds){
				print $ind ;
				my $p=Substitution->new();
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
				push (@{$nodes_with_sub{$ind}}, \${$self ->{static_hash_of_nodes}{$nodname}}); #вытащить из дерева по имени
			}
			print "\n";
			$subs_on_node{$nodname} = \%substs;

	}
	close XPAR;
	return (\%subs_on_node, \%nodes_with_sub);
	}

sub shuffle_incidence_matrix { 
	my $self = shift;
	my $file = $self->matrixPath;
 	if (-e $file){
 				print "Shuffling incidence matrix..\n";
 				my $logs = capture ('Rscript MatrixShuffler.R --file '.$file);
 				print $logs."\n";
 	}
 	else {
 				print "No incidence matrix found at $file!\n";
 	}		
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
    
    
## like calc_patristic_distance, but returns 0 for two sequential nodes    
    sub calc_my_distance {
        my ( $node, $other_node ) = @_;
        my $patristic_distance = 0;
        my $mrca    = get_mrcn($node, $other_node);
        my $mrca_id = $mrca->get_id;
        if ( $node->get_id == $mrca_id || $other_node->get_id == $mrca_id){
        	return 0;
        }
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
    

sub reversals_hist {
	my $self = shift;
	my $step = shift;
	my $separately = shift;
	my $tag = shift;
	my $site_nodes = shift; #array of arrays (site, node)
	my $groupfile = shift;

	my $root = $self->{static_tree}-> get_root;
	my @array;
	my %closest_ancestors;
	$self -> set_timestamps();
	$root->set_generic("-closest_ancestors" => \%closest_ancestors);
	my @args = ($step, $root);
	$self->visitor_coat ($root, \@array,\&reversion_ratio_visitor,\&no_check,\@args,0);
	my %hist;
	if (!@{$site_nodes}){
			foreach my $ancname (keys %{$self->{static_subtree_info}}){
				foreach my $site_index (keys %{$self->{static_subtree_info}{$ancname}}){
					 if ($self->{static_subtree_info}{$ancname}{$site_index}{"reversions"}){
						print $ancname.",".$site_index."\n";
						push @{$site_nodes}, [$site_index,$ancname];
					}
					
				}
			}
	}
	if (!$separately){
		foreach my $site_node (@{$site_nodes}){
			print $site_node->[0]."--".$site_node->[1]."\n";
			my $sitenodehist = $self->{static_subtree_info}{$site_node->[1]}{$site_node->[0]}{"reversions"};
			foreach my $bin (keys %{$sitenodehist}){
					$hist{$bin}[0] += $sitenodehist->{$bin}[0]; #reversions
					$hist{$bin}[1] += $sitenodehist->{$bin}[1]; #not reversions
			}
		}
	}
	
	my $dir = File::Spec -> catdir($self -> {static_output_base}, "reversals", $tag);
	make_path($dir);
	my $filepath = File::Spec -> catfile($dir, $self->{static_protein}."_reversals_hist");
	open HIST, ">$filepath";
	my %bins;
	if (!$separately){ 
		print HIST "bin,rev,norev\n";
		foreach my $bin (sort { $a <=> $b } keys %hist){
			$bins{$bin} = 1;
			print HIST $bin.",".$hist{$bin}[0].",".$hist{$bin}[1]."\n";
		}
	}
	else {
		print HIST "site_node,bin,rev,norev\n";
		foreach my $site_node (@{$site_nodes}){
			my $sitenodehist = $self->{static_subtree_info}{$site_node->[1]}{$site_node->[0]}{"reversions"};
			foreach my $bin (sort { $a <=> $b } keys %{$sitenodehist}){
				$bins{$bin} = 1;
				my $rev = $sitenodehist->{$bin}[0];
				unless ($rev) {$rev = 0;}
				my $norev = $sitenodehist->{$bin}[1];
				unless ($norev) {$norev = 0;}
				print HIST $site_node->[0]."_".$site_node->[1].",".$bin.",".$rev.",".$norev."\n";
			}
		}
	}

	if ($groupfile){
		print HIST "## set of site_nodes from ".$groupfile."\n";
	}
	close HIST;
	
	$self->visitor_coat ($root, \@array,\&entrenchment_visitor,\&no_check,\@args,0,"overwrite");
	my %alivehist;
	
	foreach my $site_node (@{$site_nodes}){
			foreach my $bin (keys %bins){
				print "maxdepth ".$self->{static_subtree_info}{$site_node->[1]}{$site_node->[0]}{"maxdepth"}.", bin*step ".($bin*$step)."\n";
				if ($self->{static_subtree_info}{$site_node->[1]}{$site_node->[0]}{"maxdepth"} >= $bin*$step){$alivehist{$bin} += 1;}
		}
	}
	
	my $filepath = File::Spec -> catfile($dir, $self->{static_protein}."_alive_hist");
	open HIST, ">$filepath";
	print HIST "bin,alive\n";
	foreach my $bin (sort { $a <=> $b } keys %alivehist){
		print HIST $bin.",".$alivehist{$bin}."\n";
	}

	close HIST;

	return %hist;
}

sub print_colored_trees {
	my $self = shift;
	my @allsites = @{@_[0]};;
	my $tag = $_[1];
	my $root = $self->{static_tree}-> get_root;
	my @array;
	my $state = $self->{static_state};
	my $myCodonTable   = Bio::Tools::CodonTable->new();
	
	my %closest_ancestors;
	$root->set_generic("-closest_ancestors" => \%closest_ancestors);
	my @args = ($root);
	$self->visitor_coat ($root, \@array,\&lrt_visitor,\&no_check,\@args,0);
	my @output_files;
	my $dir = File::Spec -> catdir($self -> {static_output_base}, "trees", $tag);
	make_path($dir);
	print($root->get_name()."\n");
	print($self -> {static_fasta}{$root->get_name()}."\n");
	my $rootseq = $self -> {static_fasta}{$root->get_name()};
	
	for (my $i = 0; $i < scalar @allsites; $i++){
		my $ind = $allsites[$i];
		
		my $rootletter;
		my $rootcodon;
		if($rootseq){
			$rootcodon = substr($rootseq, 3*($ind-1),3);
			$rootletter= $myCodonTable->translate($rootcodon);
		}
		else{
			my $i = 0;
			while($root->get_child($i)){
				my $node = $root->get_child($i);
				if (!($node->is_terminal)){
					my $nodeseq = $self -> {static_fasta}{$node->get_name()};
					$rootcodon = substr($nodeseq, 3*($ind-1),3);
					$rootletter= $myCodonTable->translate($rootcodon);
					last;
				}
				$i++;
			}
		}
		my $eventsfile = $self -> {static_protein}."_treescheme_".$ind;
		my $filepath = File::Spec -> catfile($dir, $eventsfile);
		open FILE, ">$filepath";
		if ($state eq "nsyn"){
			print FILE "root:".$root->get_name()."|".$rootletter."|".$rootletter."(".$rootcodon.")\n";
		}
		else{
			print FILE "root:".$root->get_name()."|".$rootcodon."|".$rootcodon."(".$rootletter.")\n";
		}
		
		foreach my $n ( @{$self -> {static_nodes_with_sub}{$ind}} ){
			my $ancnodename = $$n->get_name();
			my %subs;
			my %color;
			my $sub = ${$self -> {static_subs_on_node}{$ancnodename}}{$ind};
			$subs{$ancnodename} = $sub->{"Substitution::derived_allele"};
			$color{$ancnodename} = 1;
			if (exists $self->{static_subtree_info}{$ancnodename}{$ind}){
				foreach my $node ( keys %{$self->{static_subtree_info}{$ancnodename}{$ind}{"lrt"}}){
					if($self->{static_subtree_info}{$ancnodename}{$ind}{"lrt"}{$node}[0]){
						my $sub = ${$self -> {static_subs_on_node}{$node}}{$ind};
						print ($sub."\n");
						$subs{$node} =  $sub->{"Substitution::derived_allele"};
						$color{$node} = 1;
					}
					else {
						$color{$node} = 1;
					}
			}
			}

			my $ancsub = $subs{$ancnodename};
			my $aaancsub = $myCodonTable->translate($ancsub);
			if ($state eq "nsyn"){
				print FILE "Ancnode:$ancnodename"."|".$aaancsub."|".$aaancsub."(".$ancsub.")\n";
			}
			else {
				print FILE "Ancnode:$ancnodename"."|".$ancsub."|".$ancsub."(".$aaancsub.")\n";				
			}
			print FILE "Events:";
			foreach my $node (keys %subs){
				my $nodesub = $subs{$node};
				if ($state eq "nsyn"){
					my $aanodesub = $myCodonTable->translate($nodesub);
					print FILE $node."|".$aanodesub.",";
				}
				else{
					print FILE $node."|".$nodesub.",";
				}
			}
			print FILE "\n";
			print FILE "Subtree:";
			foreach my $node (keys %color){
				print FILE $node.",";
			}
			print FILE "\n";
		}
		
		my %subs;
		foreach my $n ( @{$self -> {static_background_nodes_with_sub}{$ind}} ){
			my $ancnodename = $$n->get_name();
			my $sub = ${$self -> {static_background_subs_on_node}{$ancnodename}}{$ind};
			my $str;
			if ($state eq "nsyn"){
				$str = $sub->{"Substitution::ancestral_allele"}."->".$sub->{"Substitution::derived_allele"};
			}
			else {
				my $ancc = $sub->{"Substitution::ancestral_allele"};
				my $derc = $sub->{"Substitution::derived_allele"};
				my $ancaa = $myCodonTable->translate($ancc);
				my $deraa = $myCodonTable->translate($derc);
				$str = $ancc."(".$ancaa.")->".$derc."(".$deraa.")";
			}
			$subs{$ancnodename} = $str;
		}
		print FILE "Synonymous:";
		foreach my $node (keys %subs){
			my $str = $subs{$node};
			print FILE $node."|".$str.",";
		}
		close FILE;
		push @output_files, $filepath;
	}
	return @output_files;
}

# prints tree with all mutations in the subtree of specified mutation (site, node). 
# If there is no such mutation, warns and proceeds.

sub print_subtree_with_mutations {
	my $self = shift;
	my @muts = @{@_[0]};
	my $tag = $_[1];
	my $ete = $_[2];
	my @pvalue_tags = @{$_[3]};
	my $root = $self->{static_tree}-> get_root;
	my @array;
	my $myCodonTable   = Bio::Tools::CodonTable->new();
	
	my %closest_ancestors;
	$root->set_generic("-closest_ancestors" => \%closest_ancestors);
	my @args = ($root);
	$self->visitor_coat ($root, \@array,\&lrt_visitor,\&no_check,\@args,0);
	my @output_files;
	for (my $i = 0; $i < scalar @muts; $i++){
		my ($ind, $ancnodename) = Textbits::cleave($muts[$i]);
		my $pvalue_tag = "";
		if (@pvalue_tags){$pvalue_tag = $pvalue_tags[$i];}
		if (!exists $self->{static_subtree_info}{$ancnodename}{$ind}){
				warn "there is no mutation at $ind , $ancnodename";
		}
		
		my %sites;
		my %color;
		my $sub = ${$self -> {static_subs_on_node}{$ancnodename}}{$ind};
		$sites{$ancnodename} = $ancnodename."_".$ind."_".$sub->{"Substitution::ancestral_allele"}.
							  "(".$myCodonTable->translate($sub->{"Substitution::ancestral_allele"}).")->".
							  $sub->{"Substitution::derived_allele"}.
							  "(".$myCodonTable->translate($sub->{"Substitution::derived_allele"}).")";
		$color{$ancnodename} = "-16776961";
		foreach my $node ( keys %{$self->{static_subtree_info}{$ancnodename}{$ind}{"lrt"}}){
			if($self->{static_subtree_info}{$ancnodename}{$ind}{"lrt"}{$node}[0]){
				my $sub = ${$self -> {static_subs_on_node}{$node}}{$ind};
				print ($sub."\n");
				$sites{$node} = $node."_".$ind."_".$sub->{"Substitution::ancestral_allele"}.
							  "(".$myCodonTable->translate($sub->{"Substitution::ancestral_allele"}).")->".
							  $sub->{"Substitution::derived_allele"}.
							  "(".$myCodonTable->translate($sub->{"Substitution::derived_allele"}).")";
				$color{$node} = "-16776961";
			}
			else {
				$color{$node} = "-16776961";
			}
		}
		
#		foreach my $node (@{$self->{static_background_nodes_with_sub}{$ind}}){
#			$sites{$node} = $node."_".$ind."_".$sub->{"Substitution::ancestral_allele"}.
#							  "(".$myCodonTable->translate($sub->{"Substitution::ancestral_allele"}).")->".
#							  $sub->{"Substitution::derived_allele"}.
#							  "(".$myCodonTable->translate($sub->{"Substitution::derived_allele"}).")";
#			$color{$node} = "-3407872";
#		}

		my $file = $self -> {static_protein}."_sites_".$ind."_".$ancnodename."_".$pvalue_tag.".tre";
		my $dir = File::Spec -> catdir($self -> {static_output_base}, "trees", $tag);
		make_path($dir);
		my $filepath = File::Spec -> catfile($dir, $file);
		open TREE, ">$filepath";
		print TREE "#NEXUS\n\nbegin trees;\n";
		print TREE "\ttree $ind = [&R] ";
		my $tree_name=tree2str($self -> {static_tree},sites => \%sites, color=>\%color);
		print TREE $tree_name;
		print TREE "\nend;\n";
		my $figblock = File::Spec -> catfile(getcwd(), "figtree_block");
		open BLOCK, "<$figblock" or die "Cannot open figtree_block: $!\n";
		while (<BLOCK>){
			print TREE $_;
		}
		close BLOCK;
		close TREE;
		push @output_files, $filepath;
		
		if ($ete){
			my $eventsfile = $self -> {static_protein}."_treescheme_".$ind."_".$ancnodename."_".$pvalue_tag;
			my $filepath = File::Spec -> catfile($dir, $eventsfile);
			open FILE, ">$filepath";
			print FILE "Ancnode:$ancnodename\n";
			print FILE "Events:";
			foreach my $node (keys %sites){
				print FILE $node.",";
			}
			print FILE "\n";
			print FILE "Subtree:";
			foreach my $node (keys %color){
				print FILE $node.",";
			}
			print FILE "\n";
			close FILE;
			push @output_files, $filepath;
		}

	}
	return @output_files;
}

sub print_events {
	my $self = shift;
	my @muts = @{@_[0]};
	my $root = $self->{static_tree}-> get_root;
	my @array;
	
	my %closest_ancestors;
	$root->set_generic("-closest_ancestors" => \%closest_ancestors);
	my @args = ($root);
	$self->visitor_coat ($root, \@array,\&lrt_visitor,\&no_check,\@args,0);
	
	my $file = $self -> {static_protein}."_events";
	my $dir = File::Spec -> catdir($self -> {static_output_base}, "distances");
	make_path($dir);
	my $filepath = File::Spec -> catfile($dir, $file);
	open FILE, ">$filepath";
	for (my $i = 0; $i < scalar @muts; $i++){
		my ($ind, $ancnodename) = Textbits::cleave($muts[$i]);
		if (!exists $self->{static_subtree_info}{$ancnodename}{$ind}){
				warn "there is no mutation at $ind , $ancnodename";
		}
		print FILE $ind.",".$ancnodename.",";
		foreach my $node ( keys %{$self->{static_subtree_info}{$ancnodename}{$ind}{"lrt"}}){
			if($self->{static_subtree_info}{$ancnodename}{$ind}{"lrt"}{$node}[0]){
				my $sub = ${$self -> {static_subs_on_node}{$node}}{$ind};
				print FILE $node." ";
			}
		}
		print FILE "\n";
	}
}

sub ancestral_aa {
	my $self = shift;
	my $site = shift;
	my $nodname = shift;
	my $myCodonTable   = Bio::Tools::CodonTable->new();
	my $sub = ${$self -> {static_subs_on_node}{$nodname}}{$site};
	return $myCodonTable->translate($sub->{"Substitution::derived_allele"})."(".$sub->{"Substitution::derived_allele"}.")";
}

# prints nexus tree, on wich all mutations in the specified site are shown 

sub print_static_tree_with_mutations{
	my $self = shift;
	my $site = shift;
	my $myCodonTable   = Bio::Tools::CodonTable->new();

	my %sites;
	my %color;
	foreach my $n(@{$self -> {static_nodes_with_sub}{$site}}){

		my $sub = ${$self -> {static_subs_on_node}{$$n->get_name()}}{$site};
		$sites{$$n->get_name()} = $$n->get_name()."_".$site."_".$sub->{"Substitution::ancestral_allele"}.
							  "(".$myCodonTable->translate($sub->{"Substitution::ancestral_allele"}).")->".
							  $sub->{"Substitution::derived_allele"}.
							  "(".$myCodonTable->translate($sub->{"Substitution::derived_allele"}).")";
		$color{$$n->get_name()} = "-16776961";
	}
	
	my $file = $self -> {static_protein}."_sites_".$site.".tre";
	my $dir = File::Spec -> catdir($self -> {static_output_base}, "trees");
	make_path($dir);
	my $filepath = File::Spec -> catfile($dir, $file);
	open TREE, ">$filepath";
	print TREE "#NEXUS\n\nbegin trees;\n";
	print TREE "\ttree $site = [&R] ";
	my $tree_name=tree2str($self -> {static_tree},sites => \%sites, color=>\%color);
	print TREE $tree_name;
	print TREE "\nend;";
	close TREE;

	foreach my $trn(@{$self -> {static_nodes_with_sub}{$site}}){
		print ($$trn->get_name()."\t");
		print (${$self -> {static_subs_on_node}{$$trn->get_name()}}{$site}->{"Substitution::derived_allele"});
		print "\n";
		foreach my $trr(@{$self -> {static_nodes_with_sub}{$site}}){
		
			print "\t".calc_true_patristic_distance($$trr, $$trn)."_";
			print (${$self -> {static_subs_on_node}{$$trr->get_name()}}{$site}->{"Substitution::derived_allele"}."_");
			print $$trr->get_name();
			print "\n";
		}
		print "\n";
	}
}




sub mean_ignore_nulls{
	if (!defined $_[0]){
		return 0;
	}

	my @arr = @{$_[0]};
	my $count = 0;
	my $sum = 0;
	for my $num(@arr){
		if (defined $num){
			$count++;
			$sum += $num;
		}
	}
	if ($count == 0){
		return 0;
	}
	return $sum/$count;
	
}





sub myclone {
	my $self = shift;
	my $clone = {
			static_output_base => $self->{static_output_base},
			static_protein => $self->{static_protein},
			static_tree =>  $self->{static_tree},
			static_fasta => $self->{static_fasta},
			static_state  => $self->{static_state},
			static_alignment_length  => $self->{static_alignment_length},
			static_hash_of_nodes => $self->{static_hash_of_nodes},
			#static_distance_hash => $self->{realdata}{"distance_hash"},
			static_background_subs_on_node => $self->{static_background_subs_on_node },
			static_background_nodes_with_sub => $self->{static_background_nodes_with_sub},
			obs_vectors => clone($self->{realdata}{"obs_vectors"})  #the only structure which can (and will) be changed
	};
	
	bless $clone, ref $self; #ref $self returns class of object $self
	return $clone;
}

sub mydebugclone {
	my $self = shift;
	print " Debugging mode: no shuffling, mutmap is a copy of the initial one\n";
	my $clone = {
			static_output_base => $self->{static_output_base},
			static_protein => $self->{static_protein},
			static_tree =>  $self->{static_tree},
			static_fasta => $self->{static_fasta},
			static_state  => $self->{static_state},
			static_alignment_length  => $self->{static_alignment_length},
			static_hash_of_nodes => $self->{static_hash_of_nodes},
			#static_distance_hash => $self->{realdata}{"distance_hash"},
			static_subs_on_node => $self->{static_subs_on_node },
			static_nodes_with_sub => $self->{static_nodes_with_sub},
			static_background_subs_on_node => $self->{static_background_subs_on_node },
			static_background_nodes_with_sub => $self->{static_background_nodes_with_sub},
			obs_vectors => clone($self->{realdata}{"obs_vectors"})  #the only structure which can (and will) be changed
	};
	
	bless $clone, ref $self; #ref $self returns class of object $self
	return $clone;
}



# outputs hash of hashes used for construction of observed_vector
sub get_hashes {
	my $self = shift;
	my %res_hash;
	#print ("what must be a hashref is a ".ref($self->{static_nodes_with_sub})."\n");
	my @nodes = $self->{static_tree}->get_nodes;
	foreach my $ind(keys $self->{static_nodes_with_sub}){
		my %x_hash;
		my %y_hash;
		foreach my $node(@nodes){
			my $nname = $node->get_name();
			$x_hash{$nname} = $node->get_branch_length;
			if ($self->{static_subs_on_node}{$nname}{$ind}){
				$y_hash{$nname} = 1;			
			}
			else {
				$y_hash{$nname} = 0;
			}
		}
		$res_hash{$ind}{"x"} = \%x_hash;
		$res_hash{$ind}{"y"} = \%y_hash;
	}
	return %res_hash;	
}



# constructs a hash of observation_vectors from our object
sub get_observation_vectors {
	my $self = shift;
	if ($self->{obs_vectors}){return %{$self->{obs_vectors}};}
	else {
		my %res_hash;
		my %hash = $self ->get_hashes();
		foreach my $ind (keys %hash){
			my @arr = make_observation_vector($hash{$ind}{"x"}, $hash{$ind}{"y"});
			$res_hash{$ind} = \@arr;
		}
		return %res_hash;
	}
}

sub shuffle_observation_vectors {
	my $obs_vectors = $_[0];
	my %shuffled_obs_vectors;
	foreach my $ind (keys %{$obs_vectors}){
		#debugging version
		#my @arr = @{$obs_vectors->{$ind}};
		# true working version
		my @arr = shuffle_obsv(\@{$obs_vectors->{$ind}});
		$shuffled_obs_vectors{$ind} = \@arr;

	}
	#print Dumper (\%shuffled_obs_vectors);
	return %shuffled_obs_vectors;
}

# constructs Mutmap object given the observation_vectors (shuffled)
sub read_observation_vectors {
	my $self = shift;
	my $obs_vectors = shift;
	my %subs_on_node;
	my %nodes_with_sub;
	#my $counter = 0;	
	foreach my $ind(keys %{$obs_vectors}){
		foreach my $set(@{$obs_vectors->{$ind}}){
		if ($set->[2] == 1){
			my $nodname = $set->[0];
			my $p=Substitution->new();
			$p->position($ind);
			$p->ancestral_allele("ATG");
			$p->derived_allele("ATG");
			$subs_on_node{$nodname}{$ind} = $p;
			if (! exists $nodes_with_sub{$ind}){
					$nodes_with_sub{$ind} = ();
			}
		#$counter++;			
			push (@{$nodes_with_sub{$ind}}, \${$self->{static_hash_of_nodes}{$nodname}});
		#print "TEST1 ".${$static_hash_of_nodes{$nodname}}->get_name()."\n"; 
		#print "TEST2 ".$nodname."\n";
		}
		}
	}
	#print "There are $counter mutations in this clone\n";	
	return (\%subs_on_node, \%nodes_with_sub);
}


# prints a table: nodename,nodelength,mutnum 
# used for checking if the amount of mutations is really proprotional to the length of the branch

sub shuffle_sanity_check {
	my $self = shift;
	my $tag = shift;
	my $file = File::Spec->catfile($self->{static_output_base}, "shuffler_check".$tag);
	open FILE, ">$file";
	print FILE "nodename,length,mutnum,bkgrmutnum\n";
	foreach my $nodname (keys %{$self->{static_background_subs_on_node}}){
			my $mutnum = scalar keys %{$self->{static_subs_on_node}{$nodname}};
			my $bkgrmutnum = scalar keys %{$self->{static_background_subs_on_node}{$nodname}};
			my $node = ${$self->{static_hash_of_nodes}{$nodname}};
			my $length = $node->get_branch_length;
			print FILE $nodname.",".$length.",".$mutnum.",".$bkgrmutnum."\n";
	}
	close FILE;

}


sub shuffle_mutator {
	my $self = shift;
	my $fake_type = shift;
	my @mock_mutmaps; 
	
	#todo (not working!)
	if ($fake_type eq "sim"){
		my $restriction = shift;
		my @sites = @{$_[0]};
		my $poisson = $_[1];
		my $continue = $_[2];
		my $rh_tree_constrains = $self->get_tree_constraints($restriction, \@sites);
		my $rh_out_tree = shuffle_muts_on_tree_exp::shuffle_mutations_on_tree($self->{static_tree}, $rh_tree_constrains, $poisson, $continue); #$continue - same type mutations do not block lines (there can be any number of mutations on one line)  
		@mock_mutmaps = $self->read_shuffled_tree($rh_out_tree);
	}
	elsif ($fake_type eq "matrix"){
		unless ($self->{static_sorted_sites}){
			$self->incidence_matrix; #sets sorted nodes and nodenames, otherwise we won't be able to read the shuffled matrix
		}
		print "Fake generation method - matrix\n";
		$self->shuffle_incidence_matrix();
		@mock_mutmaps = $self->read_incidence_matrix($self->matrixPath); 
	}
	elsif ($fake_type eq "wilke"){
		print "Fake generation method - epistat (wilke, null model 4)\n";
		@mock_mutmaps = $self->read_xpar($self->xparPath);
	}
	else {
		my %obs_vectors = $self->get_observation_vectors(); # created only once, reused afterwards
		my %shuffled_obs_vectors = shuffle_observation_vectors(\%obs_vectors);
		$self->{obs_vectors} = \%shuffled_obs_vectors;
		@mock_mutmaps = $self->read_observation_vectors(\%shuffled_obs_vectors); 
	}
	$self->{static_subs_on_node} = $mock_mutmaps[0];
	$self->{static_nodes_with_sub} = $mock_mutmaps[1];
	return $self; 
}

sub read_shuffled_tree {
	my $self = shift;
	my $rh_out_tree = shift; # $rh_out_tree->{root}->{$site}= ������ ������ (���� �����))
	my %subs_on_node;
	my %nodes_with_sub;
	my $rootname = $self->{static_tree}->get_root->get_name;
	
	foreach my $ind(keys %{$rh_out_tree->{$rootname}}){
		foreach my $nodename(@{$rh_out_tree->{$rootname}{$ind}}){
			my $p=Substitution->new();
			$p->position($ind);
			$p->ancestral_allele("ATG");
			$p->derived_allele("ATG");
			$subs_on_node{$nodename}{$ind} = $p;
			if (! exists $nodes_with_sub{$ind}){
					$nodes_with_sub{$ind} = ();
			}			
			push (@{$nodes_with_sub{$ind}}, \${$self->{static_hash_of_nodes}{$nodename}});
		}
	}
	
	return (\%subs_on_node, \%nodes_with_sub);
}

sub get_tree_constraints {
	my $self = shift;
	my $restriction = $_[0];
	my @group = @{$_[1]};
	my %group_hash;
	foreach my $site(@group){ $group_hash{$site} = 1; }
	my $subtree_info = $self->{"static_subtree_info"};
	my $constraints;
	
	my $rootname = $self->{static_tree}->get_root->get_name;
	my @nodes = $self->{static_tree}->get_nodes;
	my $totlen;
	foreach my $node (@nodes){ $totlen += $node->get_branch_length; }
	foreach my $site (keys %{$self->{static_nodes_with_sub}}){
		my $mutnum = scalar @{$self->{static_nodes_with_sub}{$site}};
		my $hazard = $mutnum/$totlen;
		my $constr = Constrains->new(number_of_mutations => $mutnum, hazard => $hazard);
		$constraints->{$rootname}{$site} = $constr; 
		#print "constraints hazard for $site ".$constraints->{$rootname}{$site}->hazard()."\n"; 
	}
	return $constraints;
}	

sub FDR_all {
	my $self = shift;
	my $number_of_fakes = shift;
	my $mock_mutmap = $self->myclone();
	for (my $i = 1; $i <= $number_of_fakes; $i++){
		$mock_mutmap->shuffle_mutator();
	}
}



# interface :)
sub iterations_gulp{
	my $self = shift;
	my $iterations = shift;
	my $tag = shift;
	my $verbose = shift;
	my $memusage = shift;
	my $restriction = shift;
	my $lifetime_restr = shift;
	my $onestrip = shift;
	my $shuffler_type = shift;
	my $debugmode = shift;
	my $poisson = shift;
	my $skip_stoppers_in_simulation = shift;
	$self->iterations_gulp_subtree_shuffling($iterations,$tag,$verbose,$memusage, $restriction, $lifetime_restr, $onestrip, $shuffler_type, $debugmode, $poisson, $skip_stoppers_in_simulation); # 0 - no lifetime restriction, 1 - one strip, "exp" or "strip" shuffler, "debug" - mock shuffler
}

# new simulation method 11.01.2017 
# works with each subtree separately
# tries to assign to subtree exactly the same number of mutations as in real data (proportionally to lengths of branches)
# subtrees are totally independent
sub iterations_gulp_subtree_shuffling {
	
	my $self = shift;
	my $iterations = shift;
	my $tag = shift;
	my $verbose = shift;
	my $memusage = shift;
	my $restriction = shift;
	my $lifetime = shift; # do we need lifetime restriction for shuffler?
	my $onestrip = shift;
	my $shufflertype = shift;
	my $debugmode = shift;
	my $poisson = shift;
	my $skip_stoppers_in_simulation = shift;
		
	my $mutnum_control = $self->{static_mutnum_control};	
	if ($verbose){print "Extracting realdata..\n";}	
	my $realdata = $self->{realdata};
	my $step = $realdata->{"step"}; #bin size
	unless (defined $step) {die "Oh no, bin size in realdata is not defined. Won't proceed with simulations.\n";}
	my $ancestor_nodes = $realdata->{"ancestor_nodes"};
	my $outdir;
	if ($self->{static_output_subfolder}){
		$outdir = $self->{static_output_subfolder};
	}
	else {
		$outdir = $self->{static_output_base};
	}
	my @simulated_hists;
	
	my $length = $self->mylength();
	my @sites = (1..$length); #todo debug (1..$length)
	# subtree_info has to be ready at this point. But that's obvious, since we have realdata,
	# and that's exactly where we take it from - realdata. 
	
	my $rh_constrains;
	my $shuffler;
	if ($shufflertype eq "strip"){
		$rh_constrains = $self->get_strip_constraints($lifetime, $restriction, \@sites); 
		$shuffler = shuffle_muts_on_tree::prepare_shuffler($self->{static_tree}, $rh_constrains, $onestrip);
	} 
	elsif ($shufflertype eq "exp"){
		$rh_constrains = $self->get_constraints($restriction, \@sites, $skip_stoppers_in_simulation); 
		Dumper ($rh_constrains);
	}
	else {
		die "Unknown shuffler type\n";
	} 
	for (my $i = 1; $i <= $iterations; $i++){
		my $rh_out_subtree;
		if ($debugmode){
			#$constraints->{$node_name}{$site} = $constr;  
			foreach my $node_name (keys %{$rh_constrains}){
				foreach my $site (keys %{$rh_constrains->{$node_name}}){
					foreach my $n (@{$self->{static_nodes_with_sub}{$site}}){
						push @{$rh_out_subtree->{$node_name}{$site}}, $$n->get_name; # all nodes with sub (including those before the ancestor and after stooppers), so it will produce errors when it gets into entrenchment_for_subtrees sub 
					}
				}
			}
		}
		else {
			if ($shufflertype eq "strip"){
				$rh_out_subtree = shuffle_muts_on_tree::shuffle_mutations_on_tree($shuffler); 
			}
			elsif ($shufflertype eq "exp"){
				print "Going to shuffle..\n";
				print "mutnum control parameter is $mutnum_control\n";
				$rh_out_subtree = shuffle_muts_on_tree_exp::shuffle_mutations_on_tree($self->{static_tree}, $rh_constrains, $poisson, $mutnum_control); 
			}
		}
		
		# $rh_out_subtree->{$name}->{$site}= ������ ������ (���� �����)) #todo
		my %hash;
		print " finished shuffling\n";
		$self->entrenchment_for_subtrees($rh_out_subtree, $step, $tag, $verbose, $lifetime, $skip_stoppers_in_simulation) 
		# >new iteration string and all the corresponding data  are printed inside this sub:
		# we used to launch visitor_coat in this sub (depth_groups_entrenchment_optimized_selection_alldepths) and take subtree_info from self. 
		# But now we have no tree - only a set of overlapping subtrees which cannot be converted into a tree.
		# So.. do we have to launch entrenchment_visitor for each subtree? Nope. Better to write a new algorithm (just tree search with restrictions) that will create the same output structure
	#	my %prehash = $mock_mutmap->depth_groups_entrenchment_optimized_selection_alldepths($step,$restriction,$ancestor_nodes, "overwrite", $tag, $verbose); #step (bin size), restriction, ancestor_nodes, should I overwrite static hash?
	}
	if ($memusage){
		my $locker = Memusage->get_locker($self);
		$locker->print_memusage();
	}	
}


sub get_strip_constraints {
	my $self = shift;
	my $lifetime = shift; # do we need lifetime restriction for shuffler?
	my $restriction = $_[0];
	my @group = @{$_[1]};
	my %group_hash;
	foreach my $site(@group){
		$group_hash{$site} = 1;
	}
	my $realdata = $self->{realdata};
	my $obs_hash = $self->get_obshash($realdata, $restriction);
	my $subtree_info = $realdata->{"subtree_info"};
	my $constraints;
	
	# all bkgr mutations, that would be stoppers for any ancestor node
	my %all_stoppers;
	my @nodes = $self->{static_tree}->get_nodes;
	foreach my $site_index(@group){
		foreach my $node(@nodes){
			my $nname = $node->get_name;
			if (!($self->has_no_background_mutation($nname, $site_index))){
				push @{$all_stoppers{$site_index}}, $node;
			}
		}
	}

	foreach my $site_node(keys %{$obs_hash}){
			my ($site, $node_name) = Textbits::cleave($site_node);
			my $maxdepth = $subtree_info->{$node_name}->{$site}->{"maxdepth"};
				if ($maxdepth > $restriction && $group_hash{$site}){
					my $mutnum = $subtree_info->{$node_name}->{$site}->{"totmuts"};
					#my $stoppers = $subtree_info->{$node_name}->{$site}->{"stoppers"}; # array of nodes (false: only those stoppers that are visible in realdata go here)
					my $stoppers = $all_stoppers{$site};
					my $constr = StripConstrains->new(number_of_mutations => $mutnum, stoppers => $stoppers);
					if ($lifetime){
						$constr->lifetime($maxdepth);
					}
					$constraints->{$node_name}{$site} = $constr;  
				}
		}
	return $constraints;
}	



sub get_constraints {
	my $self = shift;
	my $restriction = $_[0];
	my @group = @{$_[1]};
	my $skip_stoppers_in_simulation = $_[2];
	my %group_hash;
	foreach my $site(@group){
		$group_hash{$site} = 1;
	}
	my $realdata = $self->{realdata};
	my $obs_hash = $self->get_obshash($realdata, $restriction);
	my $subtree_info = $realdata->{"subtree_info"};
	my $constraints;
	
	# all bkgr mutations, that would be stoppers for any ancestor node
	my %all_stoppers;
	unless ($skip_stoppers_in_simulation && $self->{static_state} eq "nsyn"){
		my @nodes = $self->{static_tree}->get_nodes;
		foreach my $site_index(@group){
			foreach my $node(@nodes){
				my $nname = $node->get_name;
				if (!($self->has_no_background_mutation($nname, $site_index))){
					push @{$all_stoppers{$site_index}}, $node;
				}
			}
		}
	}

	foreach my $site_node(keys %{$obs_hash}){
			my ($site, $node_name) = Textbits::cleave($site_node);
			my $maxdepth = $subtree_info->{$node_name}->{$site}->{"maxdepth"};
			print "maxdepth for $site $node_name is $maxdepth\n";
				if ($maxdepth > $restriction && $group_hash{$site}){
					my $mutnum = $subtree_info->{$node_name}->{$site}->{"totmuts"};
					my $totlen = $subtree_info->{$node_name}->{$site}->{"totlengths"};
					print "for $site $node_name mutnum is $mutnum and totlen is $totlen\n";
					my $hazard = $mutnum/$totlen;
					print "hazard for $site $node_name is $hazard\n";
					my $stoppers = $all_stoppers{$site};
					my $constr = Constrains->new(number_of_mutations => $mutnum, stoppers => $stoppers, hazard => $hazard);
					$constraints->{$node_name}{$site} = $constr;  
					print "constraints hazard for $site $node_name ".$constraints->{$node_name}{$site}->hazard()."\n";
				}
		}
	return $constraints;
}	



sub get_obshash {
	my $self = shift;
	my $realdata = shift;
	my $restriction = shift;
	my $rr = get_realdata_restriction($realdata);
	print "realdata restriction is $rr and we need $restriction\n"; 
	unless(defined $rr && $rr <= $restriction ){
			die "realdata restriction is undefined or is greater than get_obshash restriction: ".$rr." > $restriction \n";
	}
	my $obs_hash = \%{$realdata->{"obs_hash".$rr}};
	if ($self->{weeds}){
		foreach my $site_node(keys %{$self->{weeds}{weeds}} ){
			delete $obs_hash->{$site_node};
		}
	}
	return $obs_hash;
}

# not to be confused with check_realdata_restriction, which gets constructor arguments as an argument
sub get_realdata_restriction {
	my $realdata = shift;
	#my @obshash_restriction = map { /^obs_hash(.*)/ ? $1 : () } (keys %{$realdata});
	#return $obshash_restriction[0];
	unless (defined $realdata->{restriction}) {die "realdata restriction is not defined!\n";}
	return $realdata->{restriction};
}

# 28.12 compute norm for given restriction and/or group
sub compute_norm {
	my $self = shift;
	my $restriction = $_[0];
	my @group;
	if ($_[1]){
		@group = @{$_[1]};
	}
	else {
		@group = (1..$self->mylength());
	}
	my %group_hash;
	foreach my $ind(@group){
		$group_hash{$ind} = 1;
	}
	my $realdata = $self->{realdata};
	my $obs_hash = $self->get_obshash($realdata, $restriction);
	my $subtree_info = $realdata->{"subtree_info"};
	
	my $norm;
	foreach my $site_node(keys %{$obs_hash}){
		my ($site, $node_name) = Textbits::cleave($site_node);
		if ($group_hash{$site}){
			my $maxdepth = $subtree_info->{$node_name}->{$site}->{"maxdepth"};
			if ($maxdepth > $restriction){
				foreach my $bin(keys %{$obs_hash->{$site_node}}){
					$norm += $obs_hash->{$site_node}->{$bin}->[0];
				}
			}
		}
	}
	return $norm;
}



# 13.09.2016 compute_norm for one site_node (for single site poisson analysis)
sub compute_norm_single_site {
	my $self = shift;
	my $site_node = shift;

	my $realdata = $self->{realdata};
	my $obs_hash = $self->get_obshash($realdata, 1000); # 1000 is supposed to be bigger than any restriction in realdata, so this will just silently give you obshash from realdata
	my $subtree_info = $realdata->{"subtree_info"};
	
	my $norm;
	foreach my $bin(keys %{$obs_hash->{$site_node}}){
		$norm += $obs_hash->{$site_node}->{$bin}->[0];
	}
	return $norm;
}


#26.02 N groups

sub print_nodes_in_analysis {
	my $self = shift;
	my $prot = $self->{static_protein};
	my $restriction = $_[0];
	my @groups = @{$_[1]};
	my @names = @{$_[2]};
	my $subtract_tallest = $self->{static_subtract_tallest};
	#$self->set_distance_matrix();
	$self->set_timestamps();
	#my %matrix = incidence_matrix(); #!
	my $i = 0;
	foreach my $group(@groups){
		print $names[$i]."\n";
		my %obs_hash = $self->nodeselector(1,$restriction,$subtract_tallest, $group, $names[$i]); #bin, restriction, subtract-tallest

		#%static_ring_hash = ();
		#%static_subtree_info = ();
		$i++;
	}
}



#5.11 for entrenchment_bootstrap_full_selection_vector
# analyze real data, prepare bkg_mutmap, and obs_vector
#28.12 hash is pruned: we do not keep mutation info if its maxdepth<50
#1.08 hash is not necessarily pruned - you can set restriction to 0 to get complete data
sub prepare_real_data {
	my $self = shift;
	my $args = shift;
	my $restriction = $args->{restriction};
	my $step = $args->{step};
	my $fake = $args->{fake};
	unless(defined $step) { $step = 0.5; }
	unless(defined $restriction) { $restriction = 50; }
	my $prot = $self->{static_protein};
	#$self -> set_distance_matrix();
	$self -> set_timestamps();
	my %matrix = $self->incidence_matrix(); 
	$self -> print_incidence_matrix(\%matrix);
	my $debugnum = scalar keys %{$self ->{static_nodes_with_sub}};
	print "Very early News from prepare: static_nodes_with_sub contains $debugnum keys\n";
	# used depth_groups_entrenchment_optimized_selector_alldepths but changed it for depth_groups_entrenchment_optimized_selector_alldepths_2, because the latter
	# keeps in obs_hash info about site index as well as about node name
	my %full_obs_hash = $self -> depth_groups_entrenchment_optimized_selector_alldepths_2($step, $restriction); # bin size
	my %ancestor_nodes;
	foreach my $ancnode(keys %full_obs_hash){
		my ($site, $node_name) = Textbits::cleave($ancnode);
		$ancestor_nodes{$node_name}{$site} = 1;
		#my @splitter = split(/_/, $ancnode);
		#$ancestor_nodes{$splitter[-1]} = 1;
		#print "Ancestor: ".$splitter[-1]."\n";
	}	
	my $restricted_norm;
	my %restricted_obs_hash;
	#my %bins;
	my $debugnum = scalar keys %full_obs_hash;
	print "Early news from prepare: there are $debugnum keys in full_obs_hash\n";
	print "Will restrict to subtrees longer than $restriction \n";
	foreach my $site_node(keys %full_obs_hash){
		my ($site, $node_name) = Textbits::cleave($site_node);
		my $maxdepth = $self -> {static_subtree_info}{$node_name}{$site}{"maxdepth"};
		#print "site $site nodename $node_name maxdepth $maxdepth \n";
		if ($maxdepth > $restriction){
			foreach my $bin(keys %{$full_obs_hash{$site_node}}){
				$restricted_norm += $full_obs_hash{$site_node}{$bin}[0];
				$restricted_obs_hash{$site_node}{$bin}[0] = $full_obs_hash{$site_node}{$bin}[0];
				$restricted_obs_hash{$site_node}{$bin}[1] = $full_obs_hash{$site_node}{$bin}[1];
			}
		}
#		print "MAXBIN $maxbin\n";
	}
#	my @bins = sort  { $a <=> $b } keys %bins;
#	
	
	
#	print " NORM50 $norm50  NORM100 $norm100  NORM150 $norm150 \n";
	my %obs_vectors = $self ->get_observation_vectors();
	my $debugnum = scalar keys %restricted_obs_hash;
	print "News from prepare: there are $debugnum keys in restricted_obs_hash\n";
	my %realdata = (
		"norm".$restriction => $restricted_norm,
		step => $step,
		restriction => $restriction,
	#	bins => \@bins,
		ancestor_nodes => \%ancestor_nodes,
		obs_vectors => \%obs_vectors,
		bkg_subs_on_node => $self -> {static_background_subs_on_node},
		bkg_nodes_with_sub => $self -> {static_background_nodes_with_sub},
		#distance_hash => $self -> {static_distance_hash},
		hash_of_nodes => $self -> {static_hash_of_nodes},
		subtree_info => $self -> {static_subtree_info},
		alignment_length => $self -> {static_alignment_length},
		no_neighbour_changing => $self -> {static_no_neighbour_changing},
		mutnum_control => $self -> {static_mutnum_control},
		no_leaves => $self -> {static_no_leaves},
		include_tips => $self -> {static_include_tips},
		skip_stoppers => $self -> {static_skip_stoppers},
		syn_lengths => $self -> {static_syn_lengths},
		"obs_hash".$restriction => \%restricted_obs_hash,
	);
	
	## added at 17.11.2016 for fake mutmaps (mutmap, produced from realdata, now can be used for printing it (these hashes are necessary for depth_.._2))
		$realdata{static_subs_on_node} = $self -> {static_subs_on_node}; # if it is used for fake, it will 
		$realdata{static_nodes_with_sub} = $self -> {static_nodes_with_sub};
	##
	
	
	my $realdatapath;
	if ($fake){
		$realdatapath = $self->{static_output_subfolder};
	} 
	else {
		$realdatapath = $self->{static_output_base};
	}
	$realdatapath = File::Spec->catfile($realdatapath, $prot."_".state_tag($self->{static_state})."_realdata");
	print "Saving real_data to $realdatapath\n";

	store \%realdata, $realdatapath;
	#$self->{static_ring_hash} = ();
	#$self->{static_subtree_info} = ();
	
	#%static_subs_on_node = ();
	#%static_background_subs_on_node = ();
	#%static_nodes_with_sub = ();
	#%static_background_nodes_with_sub = ();
}

	
	
sub select_ancestor_nodes {
	my $self = shift;
	my $restriction = $_[0];
	my @group = @{$_[1]};
	my %group_hash;
	foreach my $site(@group){
		$group_hash{$site} = 1;
	}
	my $realdata = $self->{realdata};
	my $obs_hash = $self->get_obshash($realdata, $restriction);
	my $subtree_info = $realdata->{"subtree_info"};
	my %group_nodes;
	#my $debugnum = scalar keys %{$obs_hash};
	#print "there are $debugnum keys in obshash\n";
	foreach my $site_node(keys %{$obs_hash}){
			my ($site, $node_name) = Textbits::cleave($site_node);
			my $maxdepth = $subtree_info->{$node_name}->{$site}->{"maxdepth"};
			my $mutnum = $subtree_info->{$node_name}->{$site}->{"totmuts"};
			#my $stoppers = $subtree_info->{$node_name}->{$site}->{"stoppers"}; # array of nodes
			if ($maxdepth > $restriction && $group_hash{$site}){
					push @{$group_nodes{$node_name}}, $mutnum; #mutnum control: keeping mutnums for every site in the group, which mutated at this node
					#$group_nodes{$node_name} += 1; # mutnum control: changed from = 1 to += 1
					#print "group_node ".$node_name."\n";
			}
	}
	#	my $count = scalar keys %group_nodes;
	#print "Total $count\n";
	return %group_nodes;

	
}	
	
sub select_ancestor_nodes_and_sites {
	my $self = shift;
	my $restriction = $_[0];
	my @group = @{$_[1]};
	my %group_hash;
	foreach my $site(@group){
		$group_hash{$site} = 1;
	}
	my $realdata = $self->{realdata};
	my $obs_hash = $self->get_obshash($realdata, $restriction);
	my $subtree_info = $realdata->{"subtree_info"};
	my %group_nodes;
	foreach my $site_node(keys %{$obs_hash}){
			my ($site, $node_name) = Textbits::cleave( $site_node);
			my $maxdepth = $subtree_info->{$node_name}->{$site}->{"maxdepth"};
			my $mutnum = $subtree_info->{$node_name}->{$site}->{"totmuts"};
			if ($maxdepth > $restriction && $group_hash{$site}){
					$group_nodes{$node_name}{$site} = $mutnum; #mutnum control: keeping mutnums for every site in the group, which mutated at this node
			}
	}
	return %group_nodes;
}		

 

 
 sub set_weeds {
 	my $self = shift;
 	my $fails_threshold = shift;
 	my $file = $self->weeds_file();
 	my $weeds = Weeds->new(File::Spec->catdir($self->{static_output_base}, $self->{static_protein}));
 	my $wweeds = $weeds->worstWeeds({fails_threshold => $fails_threshold});
 	$self->{weeds} = $wweeds;
 	$wweeds->printWeeds($file);
 	return $wweeds;
 }
 
 sub describe_weeds {
 	my $self = shift;
 	return unless -e $self->weeds_file();
 	opendir(DH, $self->{static_output_base});
	my $prot = $self->{static_protein};
	my @files = grep { /^${prot}_gulpselector_vector_boot_median_test_[0-9\.]+_single_sites/ } readdir(DH);
	closedir(DH);
	my @restr;
	foreach my $file(@files){
		$file =~ /.*_([0-9\.]+)_single_sites/;
		push @restr, $1;
	}
	my @sorted = sort { $a <=> $b } @restr;
 	my $file = File::Spec->catfile($self->{static_output_base}, $prot."_gulpselector_vector_boot_median_test_".$sorted[0]."_single_sites");
 	my $weeds = Weeds->readWeeds($self->weeds_file());
 	open FILE, "<$file" or die "Cannot open single sites file $file for reading: $!\n";
 	while (<FILE>){
 	}

 }
 
 sub weeds_file {
 	my $self = shift;
 	return File::Spec->catfile($self->{static_output_base},$self->{static_protein}."_weeds");
 }

 
sub count_iterations {
	my $self = shift;
	my $prot = $self->{static_protein};
	my $dir = $self->{static_output_base};
	my $dirname = File::Spec->catdir($dir, $prot); 
	my @files = Textbits::iterationFiles($dirname);
	unless (scalar @files > 0){
		return 0;
	}
	my $counter = 0;
	foreach my $gulp_filename(@files){
		next if (-d $gulp_filename);
		my $fullpath = File::Spec -> catfile($dirname, $gulp_filename);
		open GULP, "<$fullpath" or die "Cannot open $fullpath";
		while(<GULP>){
			if ($_ =~ /^>.*/){$counter++;}
		}
	}
	return $counter;
}

sub iterations_maxtag {
	my $self = shift;
	my $prot = $self->{static_protein};
	my $dir = $self->{static_output_base};
	my $dirname = File::Spec->catdir($dir,  $prot); 
	make_path ($dirname);
	opendir(DH, $dirname);
	my @files = readdir(DH);
	closedir(DH);
	unless (scalar @files > 0){
		return 0;
	}
	my @tags = 0;
	foreach my $gulp_filename(@files){
		next if (-d $gulp_filename);
		if ($gulp_filename =~ /.*_([0-9]+)$/){
			push @tags, $1;
		}
		else {
			print " tag in iterations file $gulp_filename ! It might be an error, and it might not.\nWon't change my behavior because of some stupid tag, just wanted you to be aware of it.\n";
		}
	}
	my $maxtag = List::Util::max(@tags);
	return $maxtag;
}


#27.01 writes to files immediately, does not hold hash of iterations in memory	
	
sub concat_and_divide_simult {
	my $self = shift;
	my $prot = $self->{static_protein};
	my @maxdepths = @{$_[0]};
	my @groups = @{$_[1]};
	my @group_names = @{$_[2]};
	my $mutnum_control = $self->{static_mutnum_control};
	my $subtract_maxpath = $self->{static_subtract_tallest};
	my $dir = $self->{static_output_base};
	my $subdir = $self->{static_output_subfolder};
	
	if (! defined $subdir){
		$subdir = $dir;
	}
	my $nodecount_file = File::Spec->catfile($subdir, $prot."_nodecount");
#	open NODECOUNT, ">$nodecount_file" or die "Cannot create $nodecount_file";
	
	my $realdata = $self->{realdata};
	my %hash;
	my $iteration_number = 0;
	
	my %norms;
	foreach my $md(@maxdepths){
		foreach my $group_number(0.. scalar @groups - 1){
			$norms{$md}[$group_number] = $self->compute_norm($md, $groups[$group_number]);
		}
		
	}
	
	my %group_hashes;
	foreach my $md(@maxdepths){
		foreach my $group_number(0.. scalar @groups - 1){
			my %node_hash = $self->select_ancestor_nodes($md, \@{$groups[$group_number]});
			foreach my $node_name (keys %node_hash){
				$group_hashes{$md}[$group_number]{$node_name} = \@{$node_hash{$node_name}}; #array of subtree totmuts - for each of the group sites, which mutated at this node
			}
		}
		
	}
	
	my %counter_hashes = %{ dclone(\%group_hashes) }; # tracks subtrees (node and number of muts) which where already found
	my %label_hashes; # keeps labels for hashes that are currently being filled in (used instead of {$simsite} label in _single_sites sub)
	my %sums;
	my %hash;
	
	my $dirname = File::Spec->catdir($dir, $prot); 
	my @files = Textbits::iterationFiles($dirname);
	
	my %filehandles;

	foreach my $md(@maxdepths){
		foreach my $group_number(0.. scalar @groups - 1){
			local *FILE;
			my $csvfile =  File::Spec->catfile($subdir, temp_tag(),$prot."_gulpselector_vector_".$md."_".$group_names[$group_number].".csv");
			open FILE, ">$csvfile" or die "Cannot create $csvfile";
			
			FILE->autoflush(1);
			$filehandles{$md}{$group_number} = *FILE;
		}
		
	}

	local *printer = sub {		
				foreach my $md(@maxdepths){ 
					foreach my $group_number(0..scalar @groups-1){
							
						unless (! exists $label_hashes{$md}[$group_number]){
							my $new_labels = $label_hashes{$md}[$group_number]{"current"} - $label_hashes{$md}[$group_number]{"printed"};
							if($new_labels > 1){
								#print "found ".($new_labels-1)." new labels for group ".$group_names[$group_number]."\n";
								my $next_to_print = $label_hashes{$md}[$group_number]{"printed"}+1;
								foreach my $label($next_to_print..($label_hashes{$md}[$group_number]{"current"}-1)){ # last hash (with current label) is uncomplete
									my $max = max(keys %{$hash{$md}[$group_number]{$label}});
									print "maxbin is $max\n";
									my @bins = (1..$max);
									#print " looking at label $label \n";
									if ($sums{$md}[$group_number]{$label} == 0){ # now it's impossible :)
										foreach my $bin(@bins){
											$hash{$md}[$group_number]{$label}{$bin}[0] = "NA";
											$hash{$md}[$group_number]{$label}{$bin}[1] = "NA";
										}
									}
									else {
										#print "CONC "."norm ".$norms{$md}[$group_number]."\n";
										#print "CONC "."sum ".$sums{$md}[$group_number]{$label}."\n";
										#print "CONC "."in hash, bin 12: ".$hash{$md}[$group_number]{$label}{12}[0]."\n";
										
										foreach my $bin(@bins){
											$hash{$md}[$group_number]{$label}{$bin}[0] = $hash{$md}[$group_number]{$label}{$bin}[0]*$norms{$md}[$group_number]/$sums{$md}[$group_number]{$label};
											$hash{$md}[$group_number]{$label}{$bin}[1] = $hash{$md}[$group_number]{$label}{$bin}[1]*$norms{$md}[$group_number]/$sums{$md}[$group_number]{$label};
										}
										#print "CONC expect  ".$hash{$md}[$group_number]{$label}{2}[1]." at bin 2  \n";
									}
									
									my $filehandle = $filehandles{$md}{$group_number};
									#print "CONC "."going to print something\n";
									print "CONC now expect  ".$hash{$md}[$group_number]{$label}{2}[1]." at bin 2  \n";
									#foreach my $bin(@bins){
									#	print $filehandle $bin.",".$bin.",";
									#}
									#print $filehandle "\n";
									foreach my $bin(@bins){
										print $filehandle $hash{$md}[$group_number]{$label}{$bin}[0].",".$hash{$md}[$group_number]{$label}{$bin}[1].",";
									}
									print $filehandle "\n";
									#print "CONC and now expect  ".$hash{$md}[$group_number]{$label}{2}[1]." at bin 2 (4th column) \n";
									delete $hash{$md}[$group_number]{$label};
									$label_hashes{$md}[$group_number]{"printed"} = $label;
								}
							}
						}
				
					}
				}
			};
	
	foreach my $gulp_filename(@files){
	next if (-d $gulp_filename);
	my $fullpath = File::Spec -> catfile($dirname, $gulp_filename);
	open GULP, "<$fullpath" or die "Cannot open $fullpath";

		my $str = <GULP>;
		while ($str){
				if ($str=~ /^>/) {  #skipping  ">iteration" lines 
					$iteration_number++;
					if ($iteration_number == 1 || ($iteration_number>=50 && $iteration_number%50 == 0)){
						print "Iteration number ".$iteration_number."\n";
					}
					$str = <GULP>;
				}
				#print ">iteration number $iteration_number\n";
				my $max_depth;
				my $simsite;
				my $simnode;
				my $simmutnum;

				#print "str is now $str\n";

				if ($str =~ /^site/){ # no more than a sanity check
					my @site_array = split(/\s+/, $str);
					$simsite = $site_array[1]; # careful
					$simnode = $site_array[3];
					$max_depth = $site_array[5];
					$simmutnum = $site_array[7];	
					#print "CONC ".$simsite." site,".$simnode." node,".$max_depth." maxdepth, ".$simmutnum."\n";
					$str = <GULP>; #5.02
				}
				else {print "something failed! wtf in reading iteration gulp\n"};
				
				my @bin_data;
				my $test_obs_summ = 0;
				my $test_exp_summ = 0;	
				while ($str =~ /^[0-9]/){
					my @onebin_array = split(/,/, $str);
					$test_obs_summ += $onebin_array[1];
					$test_exp_summ += $onebin_array[2];
					push @bin_data, \@onebin_array;
					$str = <GULP>;			
				}
				unless ($test_obs_summ - $test_exp_summ < 0.0001 && -$test_obs_summ + $test_exp_summ < 0.0001){
						print "Error! one site summtest failed! $simsite $simnode obssum $test_obs_summ, expsum $test_exp_summ\n";
				}
				#print "read all bins..str is now $str\n";
				# now $str is site.. or new iteration
				
				#print " Maxdepth $max_depth itnum $iteration_number bin ".$str_array[0]." exp ".$str_array[2]." obs ".$str_array[1]." \n";
				foreach my $md(@maxdepths){
					if ($max_depth > $md){
						#print "max_depth $max_depth , md $md\n";
						foreach my $group_number(0..scalar @groups-1){
							my $label;
							if (exists $label_hashes{$md}[$group_number]{"current"}) {$label = $label_hashes{$md}[$group_number]{"current"};}
							else {								
								$label = 1;
								$label_hashes{$md}[$group_number]{"current"} = 1;
							}
							if ($counter_hashes{$md}[$group_number]{$simnode}){ 
								my $mutnums = $counter_hashes{$md}[$group_number]{$simnode};
								#print " need ".scalar @{$mutnums}." sites for ".$group_names[$group_number]."\n";
								for (my $i = 0; $i < scalar @{$mutnums}; $i++){
									my $diff = abs($simmutnum - $mutnums->[$i])/$mutnums->[$i];
									#print "simmutnum $simmutnum need ".$mutnums->[$i]." difference ".$diff."\n";
									if ($diff <= $mutnum_control){
										## ! we need bin loop here, not outside!
										## hash hash hash
										#print " sum for ".$group_names[$group_number]." was ".$sums{$md}[$group_number]{$label}."\n";
										#print "CONC "."group number $group_number md $md node name $simnode simmutnum $simmutnum realmutnum".$mutnums->[$i]."\n";
										foreach my $bindat (@bin_data){
											$sums{$md}[$group_number]{$label} += $bindat->[1];
											$hash{$md}[$group_number]{$label}{$bindat->[0]}[1] += $bindat->[2];
											$hash{$md}[$group_number]{$label}{$bindat->[0]}[0] += $bindat->[1];
										}
										#print " sum for ".$group_names[$group_number]." is now ".$sums{$md}[$group_number]{$label}."\n";
										#print "before splicing ".scalar @{$mutnums}.", ";
										splice (@{$mutnums}, $i, 1);
										#print "after splicing ".scalar @{$mutnums}."\n";
										last; #we cant use one subtree more than once (in the same group)
									}
								}
								if (scalar @{$mutnums} == 0){
									#print "CONC no more mutnums for group number $group_number md $md node name $simnode, deleting this node\n";
									my $older = scalar keys %{$counter_hashes{$md}[$group_number]};
									#print "had $older , ";
									delete $counter_hashes{$md}[$group_number]{$simnode};
									my $newer = scalar keys %{$counter_hashes{$md}[$group_number]};
									#print "now has $newer \n";
								}
								if ((scalar keys %{$counter_hashes{$md}[$group_number]}) == 0){
									#print "Win! sum for label $label $md ".$group_names[$group_number]." is ".$sums{$md}[$group_number]{$label}."\n";
									#print "CONC no more nodes for group number $group_number md $md, starting new collection\n";
									$label_hashes{$md}[$group_number]{"current"} += 1; # label corresponds to number of full collections (undef, if there is 0, 1, if we got one. But collection with this label is not full!)
									$counter_hashes{$md}[$group_number] = dclone($group_hashes{$md}[$group_number]);
								}

							}
						}
					}
				}

				
	# because we iterate through site earlier $str = <GULP>;
			
			# maxbins are different for every iteration. Find maximum and use it.
		#	$maxbin = max($maxbin, $bin_data[$#bin_data]->[0]);
		#	print " here maxbin is $maxbin\n";
			#print "CONC "."maxbin $maxbin\n";	
			
			if ($iteration_number>=50 && $iteration_number%50 == 0){
				printer(); # print to .csv files every 50 iterations, delete printed from memory
			}
			
		}
		close GULP;	
	}
	
	## here 
	printer();
	
#	close NODECOUNT;
	
	foreach my $md(@maxdepths){
		foreach my $group_number(0.. scalar @groups - 1){	
					print "For ".$group_names[$group_number]." ".$md." we needed ";
					foreach my $snode (keys %{$group_hashes{$md}[$group_number]}){
						print "$snode :";
						foreach my $mnum (@{$group_hashes{$md}[$group_number]{$snode}}){
							print "$mnum".",";
						}
						print "; ";
					}
					print "\n";
					
					my $filehandle = $filehandles{$md}{$group_number};
					close $filehandle;
		}
		
	}
	
	
	
	
#	# write to file: out of this cycle
#			open CSV50, ">/export/home/popova/workspace/perlCoevolution/TreeUtils/Phylo/MutMap/".$prot."_gulpselector_vector_50_".$tag.".csv";
#			open CSV100, ">/export/home/popova/workspace/perlCoevolution/TreeUtils/Phylo/MutMap/".$prot."_gulpselector_vector_100_".$tag.".csv";
#			open CSV150, ">/export/home/popova/workspace/perlCoevolution/TreeUtils/Phylo/MutMap/".$prot."_gulpselector_vector_150_".$tag.".csv";
#			# normalization: in cycle or also out of it, if sums are kept in separate hash
#			foreach my $i (1..$iteration_number-1){
#				foreach my $bin(1..$maxbin){
#					print CSV50 $hash{50}[$i]{$bin}[0].",".$hash{50}[$i]{$bin}[1].",";
#					print CSV100 $hash{100}[$i]{$bin}[0].",".$hash{100}[$i]{$bin}[1].",";
#					print CSV150 $hash{150}[$i]{$bin}[0].",".$hash{150}[$i]{$bin}[1].",";
#				}
#				print CSV50 "\n";
#				print CSV100 "\n";
#				print CSV150 "\n";
#			}	
#			close CSV50;
#			close CSV100;
#			close CSV150;
}	
	
	
# todo. how to manage failed subtrees? 	
sub concat_and_divide_simult_for_mutnum_controlled {
	my $self = shift;
	my $prot = $self->{static_protein};
	my ($args) = @_;
	my @maxdepths = @{$args->{restriction_levels}};
	my @groups = @{$args->{groups}};
	my @group_names = @{$args->{group_names}};
	my $mutnum_control = $self->{static_mutnum_control};
	my $subtract_maxpath = $self->{static_subtract_tallest};
	my $dir = $self->{static_output_base};
	my $subdir = $self->{static_output_subfolder};
	
	if (! defined $subdir){
		$subdir = $dir;
	}
	my $nodecount_file = File::Spec->catfile($subdir, $prot."_nodecount");
#	open NODECOUNT, ">$nodecount_file" or die "Cannot create $nodecount_file";
	
	my $realdata = $self->{realdata};
	my %hash;
	my $iteration_number = 0;
	
	my %norms;
	foreach my $md(@maxdepths){
		foreach my $group_number(0.. scalar @groups - 1){
			$norms{$md}[$group_number] = $self->compute_norm($md, $groups[$group_number]);
			print "norm for $md $group_names[$group_number] ".$norms{$md}[$group_number]."\n";
		}
		
	}
	
	my %group_hashes;
	foreach my $md(@maxdepths){
		foreach my $group_number(0.. scalar @groups - 1){
			my %node_hash = $self->select_ancestor_nodes_and_sites($md, \@{$groups[$group_number]}); #
			my $sas;
			foreach my $node_name (keys %node_hash){
				$group_hashes{$md}[$group_number]{$node_name} = $node_hash{$node_name}; # $hash{$site}= totmut, now $group_hashes{$md}[$group_number]{$node_name}{$site} = totmut
			}
		}
		
	}
	
	my %counter_hashes = %{ dclone(\%group_hashes) }; # tracks subtrees (node and site) which where already found
	my %label_hashes; # keeps labels for hashes that are currently being filled in (used instead of {$simsite} label in _single_sites sub)
	my %sums;
	my %hash;
	
	my $dirname = File::Spec->catdir($dir, $prot); 
	my @files = Textbits::iterationFiles($dirname);
	
	my %filehandles;

	foreach my $md(@maxdepths){
		foreach my $group_number(0.. scalar @groups - 1){
			local *FILE;
			my $csvfile =  File::Spec->catfile($subdir, temp_tag(),$prot."_gulpselector_vector_".$md."_".$group_names[$group_number].".csv");
			open FILE, ">$csvfile" or die "Cannot create $csvfile";
			
			FILE->autoflush(1);
			$filehandles{$md}{$group_number} = *FILE;
		}
		
	}

	local *printer = sub {	
				foreach my $md(@maxdepths){ 
					foreach my $group_number(0..scalar @groups-1){
							
						unless (! exists $label_hashes{$md}[$group_number]){
							my $new_labels = $label_hashes{$md}[$group_number]{"current"} - $label_hashes{$md}[$group_number]{"printed"};
							if($new_labels > 1){
								#print "found ".($new_labels-1)." new labels for group ".$group_names[$group_number]."\n";
								my $next_to_print = $label_hashes{$md}[$group_number]{"printed"}+1;
								foreach my $label($next_to_print..($label_hashes{$md}[$group_number]{"current"}-1)){ # last hash (with current label) is uncomplete
									my $max = max(keys %{$hash{$md}[$group_number]{$label}});
									print "maxbin is $max\n";
									my @bins = (1..$max);
									#print " looking at label $label \n";
									if ($sums{$md}[$group_number]{$label} == 0){ # now it's impossible :)
										foreach my $bin(@bins){
											$hash{$md}[$group_number]{$label}{$bin}[0] = "NA";
											$hash{$md}[$group_number]{$label}{$bin}[1] = "NA";
										}
									}
									else {
										print "For $group_names[$group_number]\n";
										print "CONC "."norm ".$norms{$md}[$group_number]."\n";
										print "CONC "."sum ".$sums{$md}[$group_number]{$label}."\n";
										print "CONC "."in hash, bin 12: ".$hash{$md}[$group_number]{$label}{12}[0]."\n";
										
										foreach my $bin(@bins){
											$hash{$md}[$group_number]{$label}{$bin}[0] = $hash{$md}[$group_number]{$label}{$bin}[0]*$norms{$md}[$group_number]/$sums{$md}[$group_number]{$label};
											$hash{$md}[$group_number]{$label}{$bin}[1] = $hash{$md}[$group_number]{$label}{$bin}[1]*$norms{$md}[$group_number]/$sums{$md}[$group_number]{$label};
										}
										#print "CONC expect  ".$hash{$md}[$group_number]{$label}{2}[1]." at bin 2 (4th column) \n";
									}
									
									my $filehandle = $filehandles{$md}{$group_number};
									#print "CONC "."going to print something\n";
									#print "CONC now expect  ".$hash{$md}[$group_number]{$label}{2}[1]." at bin 2 (4th column) \n";
									#foreach my $bin(@bins){
									#	print $filehandle $bin.",".$bin.",";
									#}
									#print $filehandle "\n";
									foreach my $bin(@bins){
										print $filehandle $hash{$md}[$group_number]{$label}{$bin}[0].",".$hash{$md}[$group_number]{$label}{$bin}[1].",";
									}
									#print "CONC and now expect  ".$hash{$md}[$group_number]{$label}{2}[1]." at bin 2 (4th column) \n";
									print $filehandle "\n";
									delete $hash{$md}[$group_number]{$label};
									$label_hashes{$md}[$group_number]{"printed"} = $label;
								}
							}
						}
				
					}
				}
			};
	
	foreach my $gulp_filename(@files){
	next if (-d $gulp_filename);
	my $fullpath = File::Spec -> catfile($dirname, $gulp_filename);
	open GULP, "<$fullpath" or die "Cannot open $fullpath";

		my $str = <GULP>;
		while ($str){
				if ($str=~ /^>/) {  #skipping  ">iteration" lines 
					$iteration_number++;
					if ($iteration_number == 1 || ($iteration_number>=50 && $iteration_number%50 == 0)){
						print "Iteration number ".$iteration_number."\n";
					}
					$str = <GULP>;
				}
				#print ">iteration number $iteration_number\n";
				my $max_depth;
				my $simsite;
				my $simnode;
				my $simmutnum;

				#print "str is now $str\n";

				if ($str =~ /^site/){ # no more than a sanity check
					my @site_array = split(/\s+/, $str);
					$simsite = $site_array[1]; # careful
					$simnode = $site_array[3];
					$max_depth = $site_array[5];
					$simmutnum = $site_array[7];	
					print "CONC ".$simsite." site,".$simnode." node,".$max_depth." maxdepth, ".$simmutnum."\n";
					$str = <GULP>; #5.02
				}
				else {print "something failed! wtf in reading iteration gulp\n"};
				
				my @bin_data;
				my $test_obs_summ = 0;
				my $test_exp_summ = 0;	
				while ($str =~ /^[0-9]/){
					my @onebin_array = split(/,/, $str);
					$test_obs_summ += $onebin_array[1];
					$test_exp_summ += $onebin_array[2];
					push @bin_data, \@onebin_array;
					$str = <GULP>;			
				}
				unless ($test_obs_summ - $test_exp_summ < 0.0001 && -$test_obs_summ + $test_exp_summ < 0.0001){
						print "Error! one site summtest failed! $simsite $simnode obssum $test_obs_summ, expsum $test_exp_summ\n";
				}
				#print "read all bins..str is now $str\n";
				# now $str is site.. or new iteration
				
				#print " Maxdepth $max_depth itnum $iteration_number bin ".$str_array[0]." exp ".$str_array[2]." obs ".$str_array[1]." \n";
				foreach my $md(@maxdepths){
				#	if ($max_depth > $md){
						foreach my $group_number(0..scalar @groups-1){
							my $label;
							if (exists $label_hashes{$md}[$group_number]{"current"}) {$label = $label_hashes{$md}[$group_number]{"current"};}
							else {
								$label = 1;
								$label_hashes{$md}[$group_number]{"current"} = 1;
							}
							if ($counter_hashes{$md}[$group_number] && $counter_hashes{$md}[$group_number]{$simnode} && $counter_hashes{$md}[$group_number]{$simnode}{$simsite}){ 
								if ($simmutnum != $counter_hashes{$md}[$group_number]{$simnode}{$simsite}){
										print "Error! simmutnum for $simnode $simsite $simmutnum , expected ".$counter_hashes{$md}[$group_number]{$simnode}{$simsite}."\n";
								}
								if ($simmutnum){
									#print "found data for $simsite,$simnode $md $group_names[$group_number]\n";
										foreach my $bindat (@bin_data){
											$sums{$md}[$group_number]{$label} += $bindat->[1];
											$hash{$md}[$group_number]{$label}{$bindat->[0]}[1] += $bindat->[2];
											$hash{$md}[$group_number]{$label}{$bindat->[0]}[0] += $bindat->[1];
										}
										delete $counter_hashes{$md}[$group_number]{$simnode}{$simsite};
										my $l = scalar keys %{$counter_hashes{$md}[$group_number]{$simnode}};
										my $m = scalar keys %{$counter_hashes{$md}[$group_number]};
										#print "still need ".$m." for $md $group_names[$group_number] and $l for $simnode\n";
								}
								if (scalar keys %{$counter_hashes{$md}[$group_number]{$simnode}} == 0){
									delete $counter_hashes{$md}[$group_number]{$simnode};
									my $m = scalar keys %{$counter_hashes{$md}[$group_number]};
								#	print "now need ".$m." for $md $group_names[$group_number]\n";
								}
								if ((scalar keys %{$counter_hashes{$md}[$group_number]}) == 0){
									print "Win! sum for label $label $md ".$group_names[$group_number]." is ".$sums{$md}[$group_number]{$label}."\n";
									#print "CONC no more nodes for group number $group_number md $md, starting new collection\n";
									$label_hashes{$md}[$group_number]{"current"} += 1; # label corresponds to number of full collections (undef, if there is 0, 1, if we got one. But collection with this label is not full!)
									$counter_hashes{$md}[$group_number] = dclone($group_hashes{$md}[$group_number]);
									my $m = scalar keys %{$counter_hashes{$md}[$group_number]};
								#	print "and now need ".$m." for $md $group_names[$group_number]\n";
								}

							}
							
						}
				#	}
				}

				
	# because we iterate through site earlier $str = <GULP>;
			
			# maxbins are different for every iteration. Find maximum and use it.
			#$maxbin = max($maxbin, $bin_data[$#bin_data]->[0]);
			#print " here maxbin is $maxbin\n";
			#print "CONC "."maxbin $maxbin\n";	

			if ($iteration_number>=50 && $iteration_number%50 == 0){
				printer(); # print to .csv files every 50 iterations, delete printed from memory
			}
			
		}
		close GULP;	
	}
	
	## here 
	printer();
	
#	close NODECOUNT;
	
	foreach my $md(@maxdepths){
		foreach my $group_number(0.. scalar @groups - 1){	
					print "For ".$group_names[$group_number]." ".$md." we needed ";
					foreach my $snode (keys %{$group_hashes{$md}[$group_number]}){
						print "$snode :";
						foreach my $site (keys %{$group_hashes{$md}[$group_number]{$snode}}){
							print " $site muts ".$group_hashes{$md}[$group_number]{$snode}{$site}.",";
						}
						print "; ";
					}
					print "\n";
					
					my $filehandle = $filehandles{$md}{$group_number};
					close $filehandle;
		}
		
	}
	
	
}	
	



#13.09.2016 prints one file for one ancestor site_node (instead of group of nodes) 	
	
sub concat_and_divide_simult_single_sites {
	my $self = shift;
	my $prot = $self->{static_protein};
	my @maxdepths = @{$_[0]};
	my $md = min(@maxdepths);
	my $mutnum_control = $self->{static_mutnum_control};
	#my @groups = @{$_[1]};
	#my @group_names = @{$_[2]};
	my $subtract_maxpath = $self->{static_subtract_tallest};
	my $dir = $self->{static_output_base};
	my $subdir = $self->{static_output_subfolder};
	if (! defined $subdir){
		$subdir = $dir;
	}
	my $realdata = $self->{realdata};
	my $obs_hash = $self->get_obshash($realdata, List::Util::min(@maxdepths));
	my $subtree_info = $realdata->{"subtree_info"};
	
	my %hash;
	my $iteration_number = 1;
	
	my %norms; #previous: $norms{$md}[$group_number] = norm
#	foreach my $md(@maxdepths){
	foreach my $site_node(keys %{$obs_hash}){
			my ($site, $node_name) = Textbits::cleave($site_node);
			my $maxdepth = $subtree_info->{$node_name}->{$site}->{"maxdepth"};
			if ($maxdepth > $md){
				$norms{$site_node} = $self->compute_norm_single_site($site_node);
			}
	}	
#	}
	
	
	
	my $dirname = File::Spec->catdir($dir, $prot); 
	my @files = Textbits::iterationFiles($dirname);
	
		my %filehandles;
	my $lim = qx{echo `ulimit -n`};
	my $filenum = scalar keys %{$obs_hash};
	my $toomanyfiles = (($filenum+10)/$lim >= 1);
	unless ($toomanyfiles){
	#	foreach my $md(@maxdepths){
		foreach my $site_node(keys %{$obs_hash}){
				local *FILE;
				my $csvfile =  File::Spec->catfile($subdir, temp_tag(),$prot."_gulpselector_vector_".$md."_".$site_node.".csv");
				open FILE, ">$csvfile" or die "Cannot create $csvfile";
				FILE->autoflush(1);
				$filehandles{$site_node} = *FILE;
		}
			
	#	}
	}	
	
	
	foreach my $gulp_filename(@files){
	next if (-d $gulp_filename);
	my $fullpath = File::Spec -> catfile($dirname, $gulp_filename);
	open GULP, "<$fullpath" or die "Cannot open $fullpath";
	print "gulp file $gulp_filename\n";
	my $simsite; # careful	
	my $simnode;
	my $sim_site_node;
	my $str = <GULP>;
	while (<GULP>){

			my $max_depth;
			my @str_array;
			$str = $_;
			#my $str = <GULP>;
			
			my %sums;
			my %hash;
			
			my $test_obs_summ;
			my $test_exp_summ;
			
			while ($str =~ /^[^>]/){ 

				if ($str =~ /^site/){
					if ($test_obs_summ - $test_exp_summ < 0.0001 && -$test_obs_summ + $test_exp_summ < 0.0001){
					#	print "summtest ok\n";
					}
					else {
						print "Error! summtest failed! $simsite $simnode obssum $test_obs_summ, expsum $test_exp_summ\n";
					}
					$test_obs_summ = 0;
					$test_exp_summ = 0;	
				
					my @str_array = split(/\s+/, $str);
					# $str_array[0] is bin number
					$simsite = $str_array[1]; # careful 
					$simnode = $str_array[3];
					$sim_site_node = Textbits::concat($simsite, $simnode);
					$max_depth = $str_array[5];
					
					
					#print $site." site,".$node_name." node,".$max_depth." maxdepth\n";
					#5.02 my $str = <GULP>; 
					$str = <GULP>; #5.02
				}
				@str_array = split(/,/, $str);
				
				$test_obs_summ += $str_array[1];
				$test_exp_summ += $str_array[2];
	
				
				#print " Maxdepth $max_depth itnum $iteration_number bin ".$str_array[0]." exp ".$str_array[2]." obs ".$str_array[1]." \n";
#				foreach my $md(@maxdepths){
				#	if ($max_depth > $md){
						#foreach my $site_node(keys %{$obs_hash}){ # commented out at 13.06.2017
							#my ($obssite, $obsnode) = cleave($site_node);
							#if ($obsnode eq $simnode){ # 28.09.2016 
								#if($obssite eq $simsite) { # 30.01.2017 
									
									if ($norms{$sim_site_node}){  # checks for maxdepth in realdata # added sim_ at 13.06.2017
									#print "group number $group_number md $md node name $node_name\n";
										
										$sums{$sim_site_node} += $str_array[1]; #obssum # added sim_ at 13.06.2017 deleted $simsite at 15.06
										$hash{$sim_site_node}{$str_array[0]}[1] += $str_array[2]; # added sim_ at 13.06.2017 deleted $simsite at 15.06
										$hash{$sim_site_node}{$str_array[0]}[0] += $str_array[1]; # added sim_ at 13.06.2017 deleted $simsite at 15.06
										 
									}
								#}
							#}
						#}
			#		}
#				}

				
				$str = <GULP>;
				
			}
	 
	
			
			## parsed all information about this iteration 

			

			
print "iteration number $iteration_number\n";	
#print "sum50 $sum50 sum100 $sum100 sum150 $sum150 norm 50 $norm50 norm 100 $norm100 norm 150 $norm150\n";		
			
#			foreach my $md(@maxdepths){ 
				foreach my $site_node(keys %{$obs_hash}){
#					foreach my $simsite(keys %{$sums{$site_node}}){ # deleted  $simsite 15.06
						#print "maxdepth $md group number $group_number \n";
						my $max = max(keys %{$hash{$site_node}});
						#my @sorted = sort  { $a <=> $b } keys %{$hash{$site_node}};
						#foreach my $k(@sorted){
						#	print $k." ";
						#}
						print "\n";
						print "maxbin is $max\n";
						my @bins = (1..$max); #deleted  $simsite from $hash{$site_node}{$simsite} 15.06
						my $diff = abs($sums{$site_node} - $norms{$site_node})/$norms{$site_node}; #deleted  $simsite from $hash{$site_node}{$simsite} 15.06
						print "$site_node obssum ".$norms{$site_node}." simsum ".$sums{$site_node}."diff $diff \n"; #deleted  $simsite from $hash{$site_node}{$simsite} 15.06
						if ($sums{$site_node} == 0 || (defined($mutnum_control) && $diff > $mutnum_control)){
							foreach my $bin(@bins){
								$hash{$site_node}{$bin}[0] = "NA"; #deleted  $simsite from $hash{$site_node}{$simsite} 15.06
								$hash{$site_node}{$bin}[1] = "NA"; #deleted  $simsite from $hash{$site_node}{$simsite} 15.06
							}
						}
						else {
							foreach my $bin(@bins){
								#print "in hash: ".$hash[$group_number]{$bin}[0]."\n";
								#print "norm ".$norms[$group_number]."\n";
								#print "sum ".$sums[$group_number]."\n";
								$hash{$site_node}{$bin}[0] = $hash{$site_node}{$bin}[0]*$norms{$site_node}/$sums{$site_node};  #deleted  $simsite from $hash{$site_node}{$simsite}  and sums 15.06
								$hash{$site_node}{$bin}[1] = $hash{$site_node}{$bin}[1]*$norms{$site_node}/$sums{$site_node};  #deleted  $simsite from $hash{$site_node}{$simsite}  and sums 15.06
							}
						}
						
						my $filehandle;
						if ($toomanyfiles){
							my $csvfile =  File::Spec->catfile($subdir, temp_tag(),$prot."_gulpselector_vector_".$md."_".$site_node.".csv");
							open $filehandle, ">>$csvfile" or die "Cannot create $csvfile";
							$filehandle->autoflush(1);	
						}
						else {		
							$filehandle = $filehandles{$site_node};						
						}
						foreach my $bin(@bins){
							print $filehandle $hash{$site_node}{$bin}[0].",".$hash{$site_node}{$bin}[1].",";
						}
						print $filehandle "\n";
						if ($toomanyfiles){
							close $filehandle; 
						}
					
#					}
				}
#			}
			
			
			
#print ">iteration number ".$iteration_number."\n";		
#foreach my $md (@maxdepths)	{
#		foreach my $site_node(keys %{$obs_hash}){
#			foreach my $simsite(keys %{$sums{$site_node}}){	
#				print "md $md site_node $site_node simsite $simsite\n";
#				print $norms{$site_node}." realdata norm (realdata integral?)\n";
#				my $iter_integral;
#				foreach my $bin(keys %{$hash{$site_node}{$simsite}}){
#					$iter_integral += $hash{$site_node}{$simsite}{$bin}[0];
#				}
#				print $iter_integral." iteration integral normalized\n";
#				print $sums{$site_node}{$simsite}." iteration integral\n";
#			}
#		}
#}
			
			
			
			
			$iteration_number++;
			if ($iteration_number%50 == 0){
				print "iteration number ".$iteration_number."\n";
			}
		}
		
		close GULP;
		
	}
	
	
#	foreach my $md(@maxdepths){
	unless ($toomanyfiles){
		foreach my $site_node(keys %{$obs_hash}){
						my $filehandle = $filehandles{$site_node};
						close $filehandle;
		}
	}
		
#	}
	

}	
	

# counter from count_pvalues
sub group_counter {
	my $self = $_[0];
	my $prot = $self -> {static_protein};
	#my $prot = $_[0];
	my @restriction_levels = @{$_[1]};
	my @groups = @{$_[2]};
	my @group_names = @{$_[3]};
	my $dir = $self -> {static_output_base};
	
		my $countfile = File::Spec->catfile($dir, $prot."_count");
	open COUNTER, ">$countfile" or die "Cannot create $countfile";
	COUNTER->autoflush(1);
	
	my $realdata =  $self -> {realdata};
	#my $realdata = get_real_data();
	
	my $step =  $realdata->{"step"};
	unless (defined $step) {die "Oh no, realdata bin size is not defined. Won't proceed with pvalues\n";}
	#print "before cycle\n";
	
	
	my $obs_hash = $self->get_obshash($realdata, List::Util::min(@restriction_levels)); # if min $restriction is less than restriction in realdata, it will die
	my $subtree_info = $realdata->{"subtree_info"};
	for my $restriction(@restriction_levels){
		print "level $restriction\n";
		
		# only for @all@
		my $group_number = scalar @groups - 1;
		print " groups ".scalar @groups - 1;
		my %group_hash;
		print " size ".scalar @{$groups[$group_number]}."\n";
		foreach my $site(@{$groups[$group_number]}){
			#print "test-1 $site\n";
			$group_hash{$site} = 1;
		}
		my %obs_hash_restricted;
		my $norm_restricted;
		## copypaste from prepare: create restricted hash	
		# Beware! obs_hash50 from real_data really contains restricted data (maxdepth>50)
		my $mutcounter;
		foreach my $site_node(keys %{$obs_hash}){
			my ($site, $node_name) = Textbits::cleave($site_node);
			my $maxdepth = $subtree_info->{$node_name}->{$site}->{"maxdepth"};
			if ($maxdepth > $restriction && $group_hash{$site}){ #15.02 moved here
			foreach my $bin(keys %{$obs_hash->{$site_node}}){
					$norm_restricted += $obs_hash->{$site_node}->{$bin}->[0];
					$obs_hash_restricted{$site_node}{$bin}[0] = $obs_hash->{$site_node}->{$bin}->[0];
					$obs_hash_restricted{$site_node}{$bin}[1] = $obs_hash->{$site_node}->{$bin}->[1];
				}
		
			}
		#	print "MAXBIN $maxbin\n";
		}
		my $count = scalar keys %obs_hash_restricted;
		print COUNTER "$restriction all group $count muts $norm_restricted\n";
		print COUNTER "$restriction all group $count muts $norm_restricted\n";
			for (my $group_number = 0; $group_number < scalar @groups - 1; $group_number++){ 
	
			## group
			my %group_hash;
			
			foreach my $site(@{$groups[$group_number]}){
				$group_hash{$site} = 1;
			}
			my %obs_hash_restricted;
			my $norm_restricted;
			## copyaste from prepare: create restricted hash	
			# Beware! obs_hash50 from real_data really contains restricted data (maxdepth>50)
			foreach my $site_node(keys %{$obs_hash}){
				my ($site, $node_name) = Textbits::cleave($site_node);
				my $maxdepth = $subtree_info->{$node_name}->{$site}->{"maxdepth"};
				foreach my $bin(keys %{$obs_hash->{$site_node}}){
					if ($maxdepth > $restriction && $group_hash{$site}){
						$norm_restricted += $obs_hash->{$site_node}->{$bin}->[0];
						$obs_hash_restricted{$site_node}{$bin}[0] = $obs_hash->{$site_node}->{$bin}->[0];
						$obs_hash_restricted{$site_node}{$bin}[1] = $obs_hash->{$site_node}->{$bin}->[1];
					}
				}
			#	print "MAXBIN $maxbin\n";
			}
			
			## end f copypaste	
			my $count = scalar keys %obs_hash_restricted;
			print  COUNTER "$restriction ".$group_names[$group_number]." group $count muts $norm_restricted\n";
			print  COUNTER "$restriction ".$group_names[$group_number]." group $count muts $norm_restricted\n";
			print COUNTER   "$restriction ".$group_names[$group_number]."\n";
			foreach my $sn(keys %obs_hash_restricted){
				print  COUNTER $sn."\n";
			}
			$group_number++;
			}
		
	}
	close COUNTER;
}	
	
	
# 19.12 after concat_and_divide	
sub count_pvalues{	
	my $self = shift;
	my ($args) = @_;
	my @restriction_levels = @{$args->{restriction_levels}};
	my @groups = @{$args->{groups}};
	my @group_names = @{$args->{group_names}};
	my $fake = $args->{fake};
	my @stattypes = ("mean", "median");
	if ($args->{stattypes}){
		print " have stattypes!\n";
		@stattypes = @{$args->{stattypes}}; # "bp" (W,  test statistic of Barlow-Proschan�s test), "mean", "median" 
	}
	my $zscore = $args->{zscore};
	my $prot = $self -> {static_protein};
	my $dir = $self -> {static_output_base};
	my $outdir = $self -> {static_output_subfolder};
	
	my $countfile = File::Spec->catfile($outdir, $prot."_count");
	if ($fake) {
		open COUNTER, ">>$countfile" or die "Cannot create $countfile";
	}
	else {
		open COUNTER, ">$countfile" or die "Cannot create $countfile";
	}
	COUNTER->autoflush(1);
	
	my $realdata =  $self -> {realdata};
	#my $realdata = get_real_data();
	my $step =  $realdata->{"step"};
	unless (defined $step) {die "Oh no, realdata bin size is not defined. Won't proceed with pvalues\n";}
	#print "before cycle\n";
	
	
	my $obs_hash = $self->get_obshash($realdata, List::Util::min(@restriction_levels)); # if min $restriction is less than restriction in realdata, it will die
	#my $obs_hash = $realdata->{"obs_hash50"}; 
	my $subtree_info = $realdata->{"subtree_info"};
	for my $restriction(@restriction_levels){
		print "level $restriction\n";
		
		# only for @all@
		my $group_number = scalar @groups - 1;
		print " groups ".scalar @groups - 1;
		my %group_hash;
		print " size ".scalar @{$groups[$group_number]}."\n";
		foreach my $site(@{$groups[$group_number]}){
			#print "test-1 $site\n";
			$group_hash{$site} = 1;
		}
		my %obs_hash_restricted;
		my $norm_restricted;
		## copypaste from prepare: create restricted hash	
		# Beware! obs_hash50 from real_data really contains restricted data (maxdepth>50)
		foreach my $site_node(keys %{$obs_hash}){
			my ($site, $node_name) = Textbits::cleave($site_node);
			my $maxdepth = $subtree_info->{$node_name}->{$site}->{"maxdepth"};
		#	if (compare::is_neighbour_changing($self->$static_subs_on_node{$node_name}{$site}, 1) == 1) # fisk
			if ($maxdepth > $restriction && $group_hash{$site}){ #15.02 moved here
				foreach my $bin(keys %{$obs_hash->{$site_node}}){
					$norm_restricted += $obs_hash->{$site_node}->{$bin}->[0];
					$obs_hash_restricted{$site_node}{$bin}[0] = $obs_hash->{$site_node}->{$bin}->[0];
					$obs_hash_restricted{$site_node}{$bin}[1] = $obs_hash->{$site_node}->{$bin}->[1];
				}
		
			}
		#	print "MAXBIN $maxbin\n";
		}
		my $count = scalar keys %obs_hash_restricted;
		print COUNTER "$restriction all - $count\n";
		
			
		## end of copypaste	
			
		my $file = File::Spec->catfile($outdir, $prot."_gulpselector_vector_boot_median_test_".$restriction."_".$group_names[$group_number]);
		my $bfile = File::Spec->catfile($outdir, $prot."_boot_data_".$restriction."_".$group_names[$group_number]);
		my $bootfile;
		open $bootfile, ">$bfile" or die "Cannot create $bfile";
		
		my $outputfile;
		if ($fake){
			open $outputfile, ">>$file" or die "Cannot create $file";
		}
		else {
			open $outputfile, ">$file" or die "Cannot create $file";
		}
				
		my %histhash;
		foreach my $site_node(keys %obs_hash_restricted){
			foreach my $bin (keys %{$obs_hash_restricted{$site_node}}){
				$histhash{$bin}[0] += $obs_hash_restricted{$site_node}{$bin}[0];
				$histhash{$bin}[1] += $obs_hash_restricted{$site_node}{$bin}[1];
			}
		}
		print $outputfile "bin\tobs\texp\n";
		my @sorted_bins = sort { $a <=> $b } keys %histhash;
		foreach my $bin (@sorted_bins){
			print $outputfile $bin."\t".$histhash{$bin}[0]."\t".$histhash{$bin}[1]."\n";
		}
			
			
		my %flat_obs_hash;
		my %flat_exp_hash;
		#print " going to flat hash\n";
		foreach my $node (keys %obs_hash_restricted){
			foreach my $bin(keys %{$obs_hash_restricted{$node}}){
				$flat_obs_hash{$bin} += $obs_hash_restricted{$node}{$bin}[0]; 
				$flat_exp_hash{$bin} += $obs_hash_restricted{$node}{$bin}[1]; 
			}
		}
		
		my $test_obs_summ = sum(values %flat_obs_hash);
		my $test_exp_summ = sum(values %flat_exp_hash);
		#print ("  data  obs sum $test_obs_summ exp summ $test_exp_summ \n");
		unless (abs($test_obs_summ-$test_exp_summ) <0.00001 ){
			print "Error! data hist sum test for all failed! \n";
		}
		
		my %statdat;
		foreach my $stype (@stattypes){print $bootfile $stype." ";}
		print $bootfile "\n";
		foreach my $stype (@stattypes){
			print "stype $stype\n";
				my $stat = AgeingStat->new($stype);
				$stat->computeStats({obshash=>\%flat_obs_hash, exphash=>\%flat_exp_hash, step=>$step, zscore =>$zscore });
				$stat->printStats($outputfile);
				$statdat{$stype} = $stat;
				print $bootfile $stat->{'value'}." ";
		}
		print $bootfile "\n--------\n";
		
		my %pvals;	
		my %boots;
		
		if ($statdat{$stattypes[0]}->{'obs'} ne "NaN"){
		
		my $csvfile = File::Spec->catfile($outdir, temp_tag(),$prot."_gulpselector_vector_".$restriction."_".$group_names[$group_number].".csv");
		open CSVFILE, "<$csvfile" or die "Cannot open $csvfile";
		my $iteration = 0;
		while(<CSVFILE>){
			my %boot_obs_hash;
			my %boot_exp_hash;
			my @splitter = split(/,/, $_);
			if ($splitter[0] eq "NA"){ next; } 			# if no appropriate nodes were produced in this iteration, it is skipped
			
			for (my $i = 0; $i < scalar @splitter; $i++){ 
				my $bin = ($i/2)+1;
				my $obs = $splitter[$i];
				$boot_obs_hash{$bin} = $obs;
				$i++;
				my $exp = $splitter[$i];
				$boot_exp_hash{$bin} = $exp;
			}
			
			my $test_obs_summ = sum(values %boot_obs_hash);
			my $test_exp_summ = sum(values %boot_exp_hash);
			#print (" obs sum $test_obs_summ exp summ $test_exp_summ \n");
			unless (abs($test_obs_summ-$test_exp_summ) <0.00001 ){
				print "Error! boot hist sum test for all failed! $test_obs_summ obs, $test_exp_summ exp\n";
			}
			my %statboot;
			foreach my $stype (@stattypes){
				my $stat = AgeingStat->new($stype);
				$stat->computeStats({obshash=>\%boot_obs_hash, exphash=>\%boot_exp_hash, step=>$step, zscore =>$zscore });
				$statboot{$stype} = $stat;
			}

			foreach my $stype (@stattypes){
				if (nearest(.00000001,$statboot{$stype}->{value}) >= nearest(.00000001,$statdat{$stype}->{value})){
					$pvals{$stype}{"env"} += 1;
				}
				if (nearest(.00000001,$statboot{$stype}->{value}) <= nearest(.00000001,$statdat{$stype}->{value})){
					$pvals{$stype}{"epi"} += 1;
				}	
			}
			foreach my $stype (@stattypes){
				print $bootfile $statboot{$stype}->{value}." ";
				push @{$boots{$stype}}, $statboot{$stype}->{value};
			}
			print $bootfile "\n";
			
			$iteration++;
		}
	
		close CSVFILE;
		print $outputfile "Number of iterations: $iteration\n";
		print $outputfile "- pvalue_epistasis  pvalue_environment\n";
		if (! $iteration){
			print "Error! no valid iterations produced.\n";
		}
		else {
			print $bootfile "------\n";
			foreach my $stype (@stattypes){
				print $outputfile $stype."_stat ".($pvals{$stype}{"epi"}/$iteration)." ".($pvals{$stype}{"env"}/$iteration)."\n";
				print $bootfile array_mean(@{$boots{$stype}})." ";
			}

		}
		$self->printFooter($outputfile);
		close $outputfile;
		close $bootfile;
		
		}
		else {
			print $outputfile "hist sum is 0";
			close $outputfile;
			close $bootfile;
		}
		
	
		my @group_stattypes = grep { /^me/ } @stattypes;
		last unless (@group_stattypes);
		
		## now for the groups (all but "all")
		for (my $group_number = 0; $group_number < scalar @groups - 1; $group_number++){ # only for even
	
			## group
			my %group_hash;
			
			foreach my $site(@{$groups[$group_number]}){
				$group_hash{$site} = 1;
			}
			my %obs_hash_restricted;
			my $norm_restricted;
			## copyaste from prepare: create restricted hash	
			# Beware! obs_hash50 from real_data really contains restricted data (maxdepth>50)
			foreach my $site_node(keys %{$obs_hash}){
				my ($site, $node_name) = Textbits::cleave($site_node);
				my $maxdepth = $subtree_info->{$node_name}->{$site}->{"maxdepth"};
				foreach my $bin(keys %{$obs_hash->{$site_node}}){
					if ($maxdepth > $restriction && $group_hash{$site}){
						$norm_restricted += $obs_hash->{$site_node}->{$bin}->[0];
						$obs_hash_restricted{$site_node}{$bin}[0] = $obs_hash->{$site_node}->{$bin}->[0];
						$obs_hash_restricted{$site_node}{$bin}[1] = $obs_hash->{$site_node}->{$bin}->[1];
					}
				}
			#	print "MAXBIN $maxbin\n";
			}
			
			## end f copypaste	
			my $count = scalar keys %obs_hash_restricted;
			print  COUNTER "$restriction ".$group_names[$group_number]." group $count ";
			foreach my $sn (keys %obs_hash_restricted){
				print COUNTER $sn."\t";
			} 
			print COUNTER "\n";
			
			if ($count == 0){
				$group_number++;
				next;
			}
			my $file = File::Spec->catfile($outdir, $prot."_gulpselector_vector_boot_median_test_".$restriction."_".$group_names[$group_number]);
			my $outputfile;
			my $bfile = File::Spec->catfile($outdir,$prot."_boot_data_".$restriction."_".$group_names[$group_number]);
			open $bootfile, ">$bfile" or die "Cannot create $bfile";
			if ($fake){
				open $outputfile, ">>$file" or die "Cannot create $file";
			}
			else {
				open $outputfile, ">$file" or die "Cannot create $file";
			}
			

			my %histhash;
			foreach my $site_node(keys %obs_hash_restricted){
				foreach my $bin (keys %{$obs_hash_restricted{$site_node}}){
					$histhash{$bin}[0] += $obs_hash_restricted{$site_node}{$bin}[0];
					$histhash{$bin}[1] += $obs_hash_restricted{$site_node}{$bin}[1];
				}
			}
			print $outputfile "bin\tobs\texp\n";
			my @sorted_bins = sort { $a <=> $b } keys %histhash;
			foreach my $bin (@sorted_bins){
				print $outputfile $bin."\t".$histhash{$bin}[0]."\t".$histhash{$bin}[1]."\n";
			}
			
			my %flat_obs_hash;
			my %flat_exp_hash;

			foreach my $site_node (keys %obs_hash_restricted){
				foreach my $bin(keys %{$obs_hash_restricted{$site_node}}){
					$flat_obs_hash{$bin} += $obs_hash_restricted{$site_node}{$bin}[0];
					$flat_exp_hash{$bin} += $obs_hash_restricted{$site_node}{$bin}[1];
				}
			}
			
			my $test_obs_summ = sum(values %flat_obs_hash);
			my $test_exp_summ = sum(values %flat_exp_hash);
			unless (abs($test_obs_summ-$test_exp_summ) <0.00001 ){
				print "Error! data hist sum test for ".$group_names[$group_number]." failed! \n";
			}
			
			my %statdat;
			foreach my $stype (@group_stattypes){
				print $bootfile $stype." ";
			}
			print $bootfile "\n";
			
			foreach my $stype (@group_stattypes){
				my $stat = AgeingStat->new($stype);
				$stat->computeStats({obshash=>\%flat_obs_hash, exphash=>\%flat_exp_hash, step=>$step, zscore =>$zscore });
				$stat->printStats($outputfile);
				$statdat{$stype} = $stat;
				print $bootfile $stat->{'value'}." ";
			}
			print $bootfile "\n---------------\n";
			
			if($statdat{$group_stattypes[0]}->{'obs'} eq "NaN" || $statdat{$group_stattypes[0]}->{'exp'} eq "NaN") {
				print $outputfile " hist sum is 0";
				close $outputfile;
				$group_number++; # to skip complement for this group
				next;
			}
			#print going to read input file\n";		
				
			my $csvfile = File::Spec->catfile($outdir,temp_tag(),$prot."_gulpselector_vector_".$restriction."_".$group_names[$group_number].".csv");
			open CSVFILE, "<$csvfile" or die "Cannot open $csvfile";
			my $iteration = 0; # counter for meaningful iterations
			my %statboot;
			my %boots;
			my $itnumber = 0; # tracks iteration number, so that group and its complement are taken from the same iteration
			while(<CSVFILE>){
				my %boot_obs_hash;
				my %boot_exp_hash;
				my @splitter = split(/,/, $_);
				if ($splitter[0] eq "NA"){ # if no appropriate nodes were produced in this iteration, it is skipped
					$itnumber++;
					next;
				}

				for (my $i = 0; $i < scalar @splitter; $i++){
					my $bin = ($i/2)+1;
					my $obs = $splitter[$i];
					$boot_obs_hash{$bin} = $obs;
					$i++;
					my $exp = $splitter[$i];
					$boot_exp_hash{$bin} = $exp;
				}

				
				my $test_obs_summ = sum(values %boot_obs_hash);
				my $test_exp_summ = sum(values %boot_exp_hash);
				unless (abs($test_obs_summ-$test_exp_summ) <0.00001 ){
					print "Error! boot hist sum test for group ".$group_names[$group_number]."failed! $test_obs_summ obs, $test_exp_summ exp \n";
				}
				
				foreach my $stype (@group_stattypes){
					my $stat = AgeingStat->new($stype);
					$stat->computeStats({obshash=>\%boot_obs_hash, exphash=>\%boot_exp_hash, step=>$step, zscore =>$zscore });
					$statboot{$stype}[$itnumber] = $stat;
					print $bootfile $stat->{'value'}." ";
					push @{$boots{$stype}}, $stat->{'value'};
				}
				print $bootfile "\n";

				$itnumber++;
				$iteration++;
			}
			close CSVFILE;
			$self->printFooter($outputfile);	
			close $outputfile;
			print $bootfile "------------\n";
			foreach my $stype (@group_stattypes){
					print $bootfile array_mean(@{$boots{$stype}})." ";
			}
			close $bootfile;
			
				print "Number of meaningful iterations for group ".$group_names[$group_number]." is $iteration (haven't looked at the complement yet)\n";	
				if ($iteration == 0) {
					$group_number++;
					next;
				}
			##complement 
			
			my %complement_hash;
			$group_number++;
			foreach my $site(@{$groups[$group_number]}){ #next
				$complement_hash{$site} = 1;
			}
			my %complement_obs_hash_restricted;
			my $complement_norm_restricted;
			## copypaste from prepare: create restricted hash	
			# Beware! obs_hash50 from real_data really contains restricted data (maxdepth>50)
			foreach my $site_node(keys %{$obs_hash}){
				my ($site, $node_name) = Textbits::cleave($site_node);
				my $maxdepth = $subtree_info->{$node_name}->{$site}->{"maxdepth"};
				foreach my $bin(keys %{$obs_hash->{$site_node}}){
					if ($maxdepth > $restriction && $complement_hash{$site}){
						$complement_norm_restricted += $obs_hash->{$site_node}->{$bin}->[0];
						$complement_obs_hash_restricted{$site_node}{$bin}[0] = $obs_hash->{$site_node}->{$bin}->[0];
						$complement_obs_hash_restricted{$site_node}{$bin}[1] = $obs_hash->{$site_node}->{$bin}->[1];
					}
				}
			#	print "MAXBIN $maxbin\n";
			}
				
			## end of copypaste	
				
			my $count = scalar keys %complement_obs_hash_restricted;
			print  COUNTER " $restriction ".$group_names[$group_number]." complement $count "; 
			foreach my $sn (keys %complement_obs_hash_restricted){
				print COUNTER $sn."\t";
			} 
			print COUNTER "\n";
			my $file = File::Spec->catfile($outdir,$prot."_gulpselector_vector_boot_median_test_".$restriction."_".$group_names[$group_number]);
			my $bfile = File::Spec->catfile($outdir,$prot."_boot_data_".$restriction."_".$group_names[$group_number]);
			open $bootfile, ">$bfile" or die "Cannot create $bfile";
			my $outputfile;
			if ($fake){
				open $outputfile, ">>$file" or die "Cannot create $file";
			}
			else {
				open $outputfile, ">$file" or die "Cannot create $file";
			}
			my %complement_flat_obs_hash;
			my %complement_flat_exp_hash;
			#print " going to flat hash\n";
			foreach my $node (keys %complement_obs_hash_restricted){
				foreach my $bin (keys %{$complement_obs_hash_restricted{$node}}){
					$complement_flat_obs_hash{$bin} += $complement_obs_hash_restricted{$node}{$bin}[0];
					$complement_flat_exp_hash{$bin} += $complement_obs_hash_restricted{$node}{$bin}[1];
				}
			}
			
			my $test_obs_summ = sum(values %complement_flat_obs_hash);
			my $test_exp_summ = sum(values %complement_flat_exp_hash);
			#print (" data obs sum $test_obs_summ exp summ $test_exp_summ \n");
			unless (abs($test_obs_summ-$test_exp_summ) <0.00001 ){
				print "Error! data hist sum test for complement of ".$group_names[$group_number]." failed! \n";
			}
			
			#print " compung hist emdian \n"	;
			my %complement_statdat;
			
			foreach my $stype (@group_stattypes){
				print $bootfile $stype." ";
			}
			print $bootfile "\n";
			
			foreach my $stype (@group_stattypes){
				my $stat = AgeingStat->new($stype);
				$stat->computeStats({obshash=>\%complement_flat_obs_hash, exphash=>\%complement_flat_exp_hash, step=>$step, zscore =>$zscore });
				$stat->printStats($outputfile);
				$complement_statdat{$stype} = $stat;
				print $bootfile $stat->{'value'}." ";
			}
			print $bootfile "\n---------------\n";
			
			
	#		print $outputfile "\n observed median: $complement_obs_median expected median: $complement_exp_median observed mean: $complement_obs_mean expected mean: $complement_exp_mean\n";
			if($complement_statdat{$group_stattypes[0]}->{'obs'} eq "NaN" || $complement_statdat{$group_stattypes[0]}->{'exp'} eq "NaN") {
				print $outputfile " hist sum is 0";
				close $outputfile;
				next;
			}
			
			my $csvfile = File::Spec->catfile($outdir,temp_tag(),$prot."_gulpselector_vector_".$restriction."_".$group_names[$group_number].".csv");
			open CSVFILE, "<$csvfile" or die "Cannot open $csvfile";
			my $iteration = 0;
			my %complement_statboot;
			my %boots;
			my $itnumber = 0;
			while(<CSVFILE>){
				my %boot_obs_hash;
				my %boot_exp_hash;
				my @splitter = split(/,/, $_);
				if ($splitter[0] eq "NA"){ # if no appropriate nodes were produced in this iteration, it is skipped
					$itnumber++;
					next;
				}

				for (my $i = 0; $i < scalar @splitter; $i++){
					my $bin = ($i/2)+1;
					my $obs = $splitter[$i];
					$boot_obs_hash{$bin} = $splitter[$i];
					$i++;
					my $exp = $splitter[$i];
					$boot_exp_hash{$bin} = $splitter[$i];
				}
				
				my $test_obs_summ = sum(values %boot_obs_hash);
				my $test_exp_summ = sum(values %boot_exp_hash);
				unless (abs($test_obs_summ-$test_exp_summ) <0.00001 ){
					print "Error! boot hist sum test for complement ".$group_names[$group_number]." failed! $test_obs_summ obs, $test_exp_summ exp \n";
				}
				
				
				foreach my $stype (@group_stattypes){
					my $stat = AgeingStat->new($stype);
					$stat->computeStats({obshash=>\%boot_obs_hash, exphash=>\%boot_exp_hash, step=>$step, zscore =>$zscore });
					$complement_statboot{$stype}[$itnumber] = $stat;
					print $bootfile $stat->{'value'}." ";
					push @{$boots{$stype}}, $stat->{'value'};
				}
				print $bootfile "\n";
				
				
				$itnumber++;
				$iteration++;
			}
			close CSVFILE;
			print $bootfile "------------\n";
			foreach my $stype (@group_stattypes){
					print $bootfile array_mean(@{$boots{$stype}})." ";
			}
			close $bootfile;
		
			## 20.09.2016
			# delete iteration data from group, if there were no nodes in complement in this iteration, and vice versa		
			for (my $it = 0; $it < $itnumber; $it++){
				foreach my $stype (@group_stattypes){
					if (! defined $complement_statboot{$stype}[$it] || ! defined $statboot{$stype}[$it]){
						$complement_statboot{$stype}[$it] = undef;
						$statboot{$stype}[$it] = undef;
					} 
				}
			}
			# now group and complement contain undefs at the same indices
			foreach my $stype (@group_stattypes){
				$complement_statboot{$stype} = [ grep defined, @{$complement_statboot{$stype}} ] ;
				$statboot{$stype} = [ grep defined, @{$statboot{$stype}} ];
			}

			my $gcount = scalar @{$statboot{$group_stattypes[0]}};
			my $ccount = scalar @{$complement_statboot{$group_stattypes[0]}};
			
			print "Number of meaningful iterations (used for pvalue estimation) for stats ".$group_stattypes[0]." for group ".$group_names[$group_number]." is ".$ccount." or ".$gcount." first one is used for division, second - for iterating through simulations\n";
			
			my $updated_iteration_number = $ccount;
			if ($ccount != $gcount){
				print "Error! complement_boot ".$group_stattypes[0]." size is not equal to group_boot  size: ".$ccount." != ".$gcount."\n";
				print "It's ok, if you mix iterations and don't print NAs in .csv. I'll just compute the difference between all full arrays\n";
				$updated_iteration_number = min($gcount, $ccount);
			}
			
			my %pvals;
			for (my $i = 0; $i < $updated_iteration_number; $i++){
				foreach my $stype (@group_stattypes){
					if (nearest(0.00000001, $statboot{$stype}[$i]{'value'} - $complement_statboot{$stype}[$i]{'value'})
						>= nearest(0.00000001, $statdat{$stype}{'value'} - $complement_statdat{$stype}{'value'})){
						$pvals{$stype}{"env"}{"enrichment"} += 1;
					}
					if (nearest(0.00000001, ($statboot{$stype}->[$i]->{'value'} - $complement_statboot{$stype}->[$i]->{'value'}))
						<= nearest(0.00000001, ($statdat{$stype}->{'value'} - $complement_statdat{$stype}->{'value'}))){
						$pvals{$stype}{"env"}{"depletion"} += 1;
					}
					if (nearest(0.00000001, -($statboot{$stype}->[$i]->{'value'}) + $complement_statboot{$stype}->[$i]->{'value'})
						>= nearest(0.00000001, -($statdat{$stype}->{'value'}) + $complement_statdat{$stype}->{'value'})){
						$pvals{$stype}{"epi"}{"enrichment"} += 1;
					}
					if (nearest(0.00000001, -($statboot{$stype}->[$i]->{'value'}) + $complement_statboot{$stype}->[$i]->{'value'})
						<= nearest(0.00000001, -($statdat{$stype}->{'value'}) + $complement_statdat{$stype}->{'value'})){
						$pvals{$stype}{"epi"}{"depletion"} += 1;
					}
					if (nearest(0.00000001, $statboot{$stype}->[$i]->{'value'}) >= nearest(0.00000001,$statdat{$stype}->{'value'})){
						$pvals{$stype}{"env"}{"group"} += 1;
					}
					if (nearest(0.00000001, $statboot{$stype}->[$i]->{'value'}) <= nearest(0.00000001,$statdat{$stype}->{'value'})){
						$pvals{$stype}{"epi"}{"group"} += 1;
					}
				}
			}
			print $outputfile "Number of iterations: ".$updated_iteration_number."\n"; # and not $iteration. changed at 26.09.2016 
			if (! $updated_iteration_number){
				print "Error! no valid iterations produced.\n";
			}
			else {
				print $outputfile "- pvalue_epistasis_enrichment pvalue_environment_enrichment pvalue_epistasis pvalue_environment\n";
				foreach my $stype (@group_stattypes){
					print $outputfile $stype."_stat ".($pvals{$stype}{"epi"}{"enrichment"}/$updated_iteration_number)." ".($pvals{$stype}{"env"}{"enrichment"}/$updated_iteration_number)." ".($pvals{$stype}{"epi"}{"group"}/$updated_iteration_number)." ".($pvals{$stype}{"env"}{"group"}/$updated_iteration_number)."\n";
				}
			}
			$self->printFooter($outputfile);
			close $outputfile;	
			
			
		}
		
	}
	close COUNTER;
}

# 13.09.2016 pvalues for every ancestor node in the tree	
sub count_single_site_pvalues{	
	my $self = shift;
	my ($args) = @_;
	my @restriction_levels = @{$args->{restriction_levels}};
	my @stattypes = ("mean", "median");
	if ($args->{stattypes}){
		@stattypes = @{$args->{stattypes}}; # "bp" (W,  test statistic of Barlow-Proschan�s test), "mean", "median" 
	}
	my $zscore = $args->{zscore};
	
	my $prot = $self -> {static_protein};
	my $restriction = min(@restriction_levels);
	my $dir = $self -> {static_output_base};
	my $outdir = $self -> {static_output_subfolder};
	my $realdata =  $self -> {realdata};
	my $step = $realdata->{"step"}; #bin size
	unless (defined $step) {die "Oh no, bin size in realdata is not defined. Won't proceed with counting pvalues.\n";}
	#print "before cycle\n";
	
	

	
	
	my $obs_hash = $self->get_obshash($realdata, List::Util::min(@restriction_levels)); # if min $restriction is less than restriction in realdata, it will die
	my $subtree_info = $realdata->{"subtree_info"};

		print "level $restriction\n";
		my %obs_hash_restricted;
		my %norm_restricted;
		
		## create restricted hash	
		foreach my $site_node(keys %{$obs_hash}){
			my ($site, $node_name) = Textbits::cleave($site_node);
			my $maxdepth = $subtree_info->{$node_name}->{$site}->{"maxdepth"};
			if ($maxdepth > $restriction){ 
			foreach my $bin(keys %{$obs_hash->{$site_node}}){
					$norm_restricted{$site_node} += $obs_hash->{$site_node}->{$bin}->[0];
					$obs_hash_restricted{$site_node}{$bin}[0] = $obs_hash->{$site_node}->{$bin}->[0];
					$obs_hash_restricted{$site_node}{$bin}[1] = $obs_hash->{$site_node}->{$bin}->[1];
				}
		
			}
		}
		##
		my $count = scalar keys %obs_hash_restricted;
	    if (!$restriction) {$restriction = 0;}
		my $file = File::Spec->catfile($outdir, $prot."_gulpselector_vector_boot_median_test_".$restriction."_single_sites");
		open my $outputfile, ">$file" or die "Cannot create $file";
		{ my $ofh = select $outputfile;
	 		 $| = 1;
	  		select $ofh;
		}

		foreach my $site_node (sort (keys %obs_hash_restricted)){
			my %flat_obs_hash;
			my %flat_exp_hash;
			foreach my $bin(keys %{$obs_hash_restricted{$site_node}}){
				$flat_obs_hash{$bin} += $obs_hash_restricted{$site_node}{$bin}[0]; 
				$flat_exp_hash{$bin} += $obs_hash_restricted{$site_node}{$bin}[1]; 
			}
			print $outputfile "site_node\tbin\tobs\texp\n";
			my @sorted_bins = sort { $a <=> $b } keys $obs_hash_restricted{$site_node};
			foreach my $bin (@sorted_bins){
				if (defined $obs_hash_restricted{$site_node}{$bin}[0] && defined $obs_hash_restricted{$site_node}{$bin}[1]){
					print $outputfile $site_node."\t".$bin."\t".$obs_hash_restricted{$site_node}{$bin}[0]."\t".$obs_hash_restricted{$site_node}{$bin}[1]."\n";
				}
			}
			print $outputfile "\n site_node: $site_node\n";
			my %statdat;
			foreach my $stype (@stattypes){
				my $stat = AgeingStat->new($stype);
				$stat->computeStats({obshash=>\%flat_obs_hash, exphash=>\%flat_exp_hash, step=>$step, zscore =>$zscore });
				$stat->printStats($outputfile);
				$statdat{$stype} = $stat;
			}
				
			my %pvals;	
			
			if ($statdat{$stattypes[0]}->{'obs'} ne "NaN"){
			
			my $csvfile = File::Spec->catfile($outdir, temp_tag(),$prot."_gulpselector_vector_".$restriction."_".$site_node.".csv");
			open CSVFILE, "<$csvfile" or die "Cannot open $csvfile";
			my $iteration = 0;
			my @array_obs_minus_exp;
			while(<CSVFILE>){
				my %boot_obs_hash;
				my %boot_exp_hash;
				my @splitter = split(/,/, $_);
				if ($splitter[0] eq "NA"){ # if no appropriate nodes were produced in this iteration, it is skipped
					next;
				}
	
				for (my $i = 0; $i < scalar @splitter; $i++){
					my $bin = ($i/2)+1;
					my $obs = $splitter[$i];
					$boot_obs_hash{$bin} = $splitter[$i];
					$i++;
					my $exp = $splitter[$i];
					$boot_exp_hash{$bin} = $exp;
				}
				
				print ("iteration ".$iteration."\n");
				my $test_obs_summ = sum(values %boot_obs_hash);
				my $test_exp_summ = sum(values %boot_exp_hash);
				
				unless (abs($test_obs_summ-$test_exp_summ) <0.00001 ){
					print "Error! boot hist sum test for $site_node failed! $test_obs_summ obs, $test_exp_summ exp\n";
					my $binnum = (scalar @splitter)/2;
					print "there are $binnum bins here\n";
				}
				my %statboot;
				foreach my $stype (@stattypes){
					my $stat = AgeingStat->new($stype);
					$stat->computeStats({obshash=>\%boot_obs_hash, exphash=>\%boot_exp_hash, step=>$step, zscore =>$zscore });
					$statboot{$stype} = $stat;
				}
	
				foreach my $stype (@stattypes){
					if (nearest(.00000001,$statboot{$stype}->{value}) >= nearest(.00000001,$statdat{$stype}->{value})){
						$pvals{$stype}{"env"} += 1;
					}
					if (nearest(.00000001,$statboot{$stype}->{value}) <= nearest(.00000001,$statdat{$stype}->{value})){
						$pvals{$stype}{"epi"} += 1;
					}	
				} 
				$iteration++;
			}

			close CSVFILE;
			
			my ($psite, $pnode_name) = Textbits::cleave($site_node);
			my $pmaxdepth = $subtree_info->{$pnode_name}->{$psite}->{"maxdepth"};
			my $pmutcount = sum(values %flat_obs_hash);
			print $outputfile "Number of iterations: $iteration\n";
				if ($iteration > 0){
					print $outputfile "#\tsite_node\tmutations\tmaxlength\t";
					foreach my $stype (@stattypes){
						print $outputfile "pvalue_epistasis(".$stype.")\t";
					}
					foreach my $stype (@stattypes){
						print $outputfile "pvalue_environment(".$stype.")\t";
					}
					print $outputfile "iterations\n";
					print $outputfile ">\t".$site_node."\t".$pmutcount."\t".$pmaxdepth."\t";
					foreach my $stype (@stattypes){
						print $outputfile ($pvals{$stype}{"epi"}/$iteration)."\t";
					}
					foreach my $stype (@stattypes){
						print $outputfile ($pvals{$stype}{"env"}/$iteration)."\t";
					}
					print $outputfile $iteration."\n";
				}
				else {
					print $outputfile "No iterations found for site_node $site_node !\n";
				}
			}
		else {
			print $outputfile "hist sum is 0";	
		}
		}
			$self -> printFooter($outputfile);
		close $outputfile;	


}

# 13.09.2016 pvalues for every ancestor node in the tree	
sub count_sites_as_groups_pvalues{	
	my $self = shift;
	my ($args) = @_;
	my @restriction_levels = @{$args->{restriction_levels}};
	my @stattypes = ("mean", "median");
	if ($args->{stattypes}){
		print " have stattypes!\n";
		@stattypes = @{$args->{stattypes}}; # "bp" (W,  test statistic of Barlow-Proschanӳ test), "mean", "median" 
	}
	my $zscore = $args->{zscore};
	my $prot = $self -> {static_protein};
	my $restriction = min(@restriction_levels);
	my $dir = $self -> {static_output_base};
	my $outdir = $self -> {static_output_subfolder};
	my $realdata =  $self -> {realdata};
	my $step = $realdata->{"step"}; #bin size
	unless (defined $step) {die "Oh no, bin size in realdata is not defined. Won't proceed with counting pvalues.\n";}
	#print "before cycle\n";

	my $obs_hash = $self->get_obshash($realdata, List::Util::min(@restriction_levels)); # if min $restriction is less than restriction in realdata, it will die
	my $subtree_info = $realdata->{"subtree_info"};

		print "level $restriction\n";
		my %obs_hash_restricted;
		my %norm_restricted;
		
		## create restricted hash	
		foreach my $site_node(keys %{$obs_hash}){
			my ($site, $node_name) = Textbits::cleave($site_node);
			my $maxdepth = $subtree_info->{$node_name}->{$site}->{"maxdepth"};
			if ($maxdepth > $restriction){ 
			foreach my $bin(keys %{$obs_hash->{$site_node}}){
					$norm_restricted{$site} += $obs_hash->{$site_node}->{$bin}->[0];
					$obs_hash_restricted{$site}{$bin}[0] += $obs_hash->{$site_node}->{$bin}->[0];
					$obs_hash_restricted{$site}{$bin}[1] += $obs_hash->{$site_node}->{$bin}->[1];
				}
		
			}
		}
		##
		my $count = scalar keys %obs_hash_restricted;
	    if (!$restriction) {$restriction = 0;}
		my $file = File::Spec->catfile($outdir, $prot."_gulpselector_vector_boot_median_test_".$restriction."_sites_as_groups");
		open my $outputfile, ">$file" or die "Cannot create $file";
		{ my $ofh = select $outputfile;
	 		 $| = 1;
	  		select $ofh;
		}

		foreach my $site (sort (keys %obs_hash_restricted)){
			my %flat_obs_hash;
			my %flat_exp_hash;
			foreach my $bin(keys %{$obs_hash_restricted{$site}}){
				$flat_obs_hash{$bin} += $obs_hash_restricted{$site}{$bin}[0]; 
				$flat_exp_hash{$bin} += $obs_hash_restricted{$site}{$bin}[1]; 
			}
			print $outputfile "site\tbin\tobs\texp\n";
			my @sorted_bins = sort { $a <=> $b } keys $obs_hash_restricted{$site};
			foreach my $bin (@sorted_bins){
				if (defined $obs_hash_restricted{$site}{$bin}[0] && defined $obs_hash_restricted{$site}{$bin}[1]){
					print $outputfile $site."\t".$bin."\t".$obs_hash_restricted{$site}{$bin}[0]."\t".$obs_hash_restricted{$site}{$bin}[1]."\n";
				}
			}
			print $outputfile "\n site: $site\n";
			my %statdat;
			foreach my $stype (@stattypes){
				my $stat = AgeingStat->new($stype);
				$stat->computeStats({obshash=>\%flat_obs_hash, exphash=>\%flat_exp_hash, step=>$step, zscore =>$zscore });
				$stat->printStats($outputfile);
				$statdat{$stype} = $stat;
			}
				
			my %pvals;	
			
			if ($statdat{$stattypes[0]}->{'obs'} ne "NaN"){
			
			my $csvfile = File::Spec->catfile($outdir, temp_tag(),$prot."_gulpselector_vector_".$restriction."_".$site.".csv");
			open CSVFILE, "<$csvfile" or die "Cannot open $csvfile";
			my $iteration = 0;
			my @array_obs_minus_exp;
			while(<CSVFILE>){
				my %boot_obs_hash;
				my %boot_exp_hash;
				my @splitter = split(/,/, $_);
				if ($splitter[0] eq "NA"){ # if no appropriate nodes were produced in this iteration, it is skipped
					next;
				}
	
				for (my $i = 0; $i < scalar @splitter; $i++){
					my $bin = ($i/2)+1;
					my $obs = $splitter[$i];
					$boot_obs_hash{$bin} = $splitter[$i];
					$i++;
					my $exp = $splitter[$i];
					$boot_exp_hash{$bin} = $exp;
				}
				
				print ("iteration ".$iteration."\n");
				my $test_obs_summ = sum(values %boot_obs_hash);
				my $test_exp_summ = sum(values %boot_exp_hash);
				
				unless (abs($test_obs_summ-$test_exp_summ) <0.00001 ){
					print "Error! boot hist sum test for $site failed! $test_obs_summ obs, $test_exp_summ exp\n";
					my $binnum = (scalar @splitter)/2;
					print "there are $binnum bins here\n";
				}
				my %statboot;
				foreach my $stype (@stattypes){
					my $stat = AgeingStat->new($stype);
					$stat->computeStats({obshash=>\%boot_obs_hash, exphash=>\%boot_exp_hash, step=>$step, zscore =>$zscore });
					$statboot{$stype} = $stat;
				}
	
				foreach my $stype (@stattypes){
					if (nearest(.00000001,$statboot{$stype}->{value}) >= nearest(.00000001,$statdat{$stype}->{value})){
						$pvals{$stype}{"env"} += 1;
					}
					if (nearest(.00000001,$statboot{$stype}->{value}) <= nearest(.00000001,$statdat{$stype}->{value})){
						$pvals{$stype}{"epi"} += 1;
					}	
				} 
				$iteration++;
			}

			close CSVFILE;
			
		#	my $pmaxdepth = $subtree_info->{$pnode_name}->{$psite}->{"maxdepth"};
			my $pmutcount = sum(values %flat_obs_hash);
			print $outputfile "Number of iterations: $iteration\n";
				if ($iteration > 0){
					print $outputfile "#\tsite_node\tmutations\tmaxlength\t";
					foreach my $stype (@stattypes){
						print $outputfile "pvalue_epistasis(".$stype.")\t";
					}
					foreach my $stype (@stattypes){
						print $outputfile "pvalue_environment(".$stype.")\t";
					}
					print $outputfile "iterations\n";
					print $outputfile ">\t".$site."\t".$pmutcount."\tNA\t"; # $pmaxdepth is NA
					foreach my $stype (@stattypes){
						print $outputfile ($pvals{$stype}{"epi"}/$iteration)."\t";
					}
					foreach my $stype (@stattypes){
						print $outputfile ($pvals{$stype}{"env"}/$iteration)."\t";
					}
					print $outputfile $iteration."\n";
				}
				else {
					print $outputfile "No iterations found for site $site !\n";
				}
			}
		else {
			print $outputfile "hist sum is 0";	
		}
		}
			$self -> printFooter($outputfile);
		close $outputfile;	


}


sub no_check{
	return 1;
}


# prints protein_for_LRT files
sub print_data_for_LRT {
	my $self = shift;
	# my $dir = File::Spec->catdir(getcwd(),"likelihood", $self->{static_state}); # before august 2016 refactoring 
	my $dir = File::Spec -> catfile($self->{static_output_base}, "likelihood");
	make_path($dir);
	my $filename = File::Spec->catfile($dir, ($self->{static_protein})."_for_LRT.csv");
	open my $file, ">$filename" or die "Cannot create $filename";
	my $root = $self->{static_tree}-> get_root;
	my @array;
	my @group = (1..$self->mylength());
	
	my %closest_ancestors;
	$root->set_generic("-closest_ancestors" => \%closest_ancestors);
	my @args = ($root);
	$self->visitor_coat ($root, \@array,\&lrt_visitor,\&no_check,\@args,0);
	print $file "site,ancestor_node,t_branch_start,t_branch_end,event_indicator\n";
	foreach my $ind (@group){
		foreach my $ancnode(@{$self->{static_nodes_with_sub}{$ind}}){
			if(ref($ancnode) eq "REF"){
				$ancnode = ${$ancnode};
			}
			my $ancnodename = $ancnode->get_name();
			foreach my $node ( keys %{$self->{static_subtree_info}{$ancnodename}{$ind}{"lrt"}}){
				print $file $ind.",".$ancnodename.",".$self->{static_subtree_info}{$ancnodename}{$ind}{"lrt"}{$node}[1].",".
				$self->{static_subtree_info}{$ancnodename}{$ind}{"lrt"}{$node}[2].",";
				my $event = 0;
				if($self->{static_subtree_info}{$ancnodename}{$ind}{"lrt"}{$node}[0]) {$event = 1};
				print $file "$event\n";
			}
		}
	}
	$self->printFooter($file);	
	close $file;
}



sub entrenchment_for_subtrees{
	my $self = shift;
	my $rh_out_subtree = shift;
	my $step = shift;
	my $tag = shift;
	my $verbose = shift;
	my $lifetime = shift;
	my$skip_stoppers_in_simulation = shift;
	my $restriction = 0;
	# $rh_out_subtree->{$name}->{$site}= array of exits
	my $dir = File::Spec->catdir($self->{static_output_base}, $self->{static_protein});
	make_path($dir);
	my $filename = File::Spec->catfile($dir, $self->{static_protein}."_for_enrichment_".$tag);
	if ($verbose){print "Going to print in file $filename \n";}
	open my $file, ">>$filename" or die "Cannot create $filename\n";
	my %hist;
	print $file ">new iteration\n";
	print ">new iteration\n";
	foreach my $nodename (keys %{$rh_out_subtree}){
		my $node = ${$self->{static_hash_of_nodes}{$nodename}};
		my %alive =  map {$_ => 1} keys %{$rh_out_subtree->{$nodename}};
		$node->set_generic("-alive" => \%alive);
		my %subtree_info;
		my @array;
		my $mutations;
		foreach my $site (keys %{$rh_out_subtree->{$nodename}}){
			my %site_muts = map {$_ => 1} @{$rh_out_subtree->{$nodename}{$site}}; #array of nodenames
			$mutations->{$site} = \%site_muts;
			#print "ancnode $nodename site $site , exits: ";
			#foreach my $nname(keys %site_muts){
			#	print $nname."\t";
			#}
			#print "\n";
		}
		
		my @args = (undef, $step, $node, \%subtree_info, $mutations, $lifetime, $skip_stoppers_in_simulation); # undef - site_index
		#subtree info is fresh and clean
		$self->my_visit_depth_first ($node, \@array,\&subtree_info,\&no_check,\@args,0);
		foreach my $site (keys %{$rh_out_subtree->{$nodename}}){
			my $site_node = Textbits::concat($site, $nodename);
			my $exits = scalar @{$rh_out_subtree->{$nodename}{$site}};
			if ($exits != $subtree_info{$site}{"totmuts"}){
				print "Error! for $site $nodename entrenchment found ".$subtree_info{$site}{"totmuts"}.", and shuffler planted $exits (not an error in debugmode, because in that case exits comprise all mutations in site (including those before ancestor and after stoppers)\n";
				foreach my $ex (@{$rh_out_subtree->{$nodename}{$site}}){
					my $exnode = ${$self->{static_hash_of_nodes}{$ex}};
					print " ex ".$ex." depth ".get_sequential_distance($node, $exnode)."\n";
				}
			}
			#print " for $site $nodename data maxdepth is ".$self->{realdata}{subtree_info}{$nodename}{$site}{"maxdepth"}." and sim maxdepth is ".$subtree_info{$site}{"maxdepth"}."\n";
			if($self->{realdata}{subtree_info}{$nodename}{$site}{"maxdepth"} < $subtree_info{$site}{"maxdepth"}){
				print "Possible error (not an error for nolim, when you do not care if simulated mutations fall out of real site_node maxdepth limits) $site $nodename went out of bounds: data maxdepth is ".$self->{realdata}{subtree_info}{$nodename}{$site}{"maxdepth"}." and sim maxdepth is ".$subtree_info{$site}{"maxdepth"}."\n";
			}
			if($self->{realdata}{subtree_info}{$nodename}{$site}{"totmuts"} != $subtree_info{$site}{"totmuts"}){
				print "Debugmode-only totmut error: $site $nodename realdata totmuts ".$self->{realdata}{subtree_info}{$nodename}{$site}{"totmuts"}." and subtree_info totmuts is ".$subtree_info{$site}{"totmuts"}."\n";
			}
			my $realdata_totmuts;
			my $realdata_totlen = $self->{realdata}{subtree_info}{$nodename}{$site}{"totlengths"};
			foreach my $bin (keys %{$self->{realdata}{subtree_info}{$nodename}{$site}{"hash"}}){
				$realdata_totmuts += $self->{realdata}{subtree_info}{$nodename}{$site}{"hash"}{$bin}[0];
				#$realdata_totlen += $self->{realdata}{subtree_info}{$nodename}{$site}{"hash"}{$bin}[1];
			}
			#print  $site_node." ".$subtree_info{$site}{"maxdepth"}." total square $square ".$subtree_info{$site}{"totmuts"}." must be $exits \n";
			my $total_muts;
			my $total_length = $subtree_info{$site}{"totlengths"};
			foreach my $bin (sort {$a <=> $b} keys %{$subtree_info{$site}{"hash"}}){
				$total_muts += $subtree_info{$site}{"hash"}{$bin}[0];
				#$total_length += $subtree_info{$site}{"hash"}{$bin}[1];
			}	
			if ($total_muts != $subtree_info{$site}{"totmuts"}){
				print "Error! total_muts $total_muts is not equal to subtree_info total_muts ".$subtree_info{$site}{"totmuts"}."\n";
			}
			if ($total_length != $realdata_totlen){
				print "Debugmode-only length error: $site $nodename total_length $total_length is not equal to realdata_totlen $realdata_totlen \n";
			}
		
			
			#print "depth ".$static_depth_hash{$site}{$node->get_name()}."\n";
			#print "maxdepth ".$static_subtree_info{$node->get_name()}{$site}{"maxdepth"}."\n";
			if ($subtree_info{$site}{"maxdepth"} > $restriction){ #10.10 restriction returns. 

				#print $node->get_name()." ".$site." TOTALS: $total_muts, $total_length, maxdepth ".$static_subtree_info{$node->get_name()}{$site}{"maxdepth"}."\n";
				if ($total_length > 0){
					#print "total muts $total_muts \n";
					print $file "site $site node $nodename maxdepth ".$subtree_info{$site}{"maxdepth"}." muts ".$total_muts." total_length ".$total_length."\n";
					my $check_local_lengths_sum;
					my $check_total_obs;
					my $check_total_exp;
					foreach my $bin (sort {$a <=> $b} (keys %{$subtree_info{$site}{"hash"}})){
							my $local_length = $subtree_info{$site}{"hash"}{$bin}[1];
							$check_local_lengths_sum += $local_length;
							$hist{$site_node}{$bin}[0] += $subtree_info{$site}{"hash"}{$bin}[0]; #observed
							$hist{$site_node}{$bin}[1] += $total_muts*$local_length/$total_length; #expected	
							#print "adding to obs bin $bin ".$static_subtree_info{$node->get_name()}{$site}{"hash"}{$bin}[0]."\n";
						    if (!$hist{$site_node}{$bin}[0]){
						    	$hist{$site_node}{$bin}[0] += 0;
						    }
						    if (!$hist{$site_node}{$bin}[1]){
						    	$hist{$site_node}{$bin}[1] += 0;
						    }

	 						print $file "$bin,".$hist{$site_node}{$bin}[0].",".$hist{$site_node}{$bin}[1]."\n";
	 						$check_total_obs += $hist{$site_node}{$bin}[0];
	 						$check_total_exp += $hist{$site_node}{$bin}[1];
					}

					unless ($total_length == $check_local_lengths_sum){
						print "Error! local length sumtest failed! total $total_length, local_sum $check_local_lengths_sum\n";
					}
					
					unless ($check_total_obs-$check_total_exp < 0.001 && -$check_total_obs+$check_total_exp < 0.001 ){
						print "Error! obsexp sumtest failed! total obs $check_total_obs, total exp $check_total_exp total_muts $total_muts site_node $site_node \n";
					}
				}
			}
		}
	}
			
		# todo delete all "-alive" hashes in this tree! So, we have to traverse the tree twice anyway. 
		# Nope! we completely reset these hashes with every traversing, so we do not need to delete them
		
		close $file;
		
}

# used by iterations_gulp_subtree_shuffling (entrenchment_for_subtrees)
# works with individual subtrees defined by an ancestor node (simultaneously with all sites that changed there)
 	sub subtree_info {
 		my $self = $_[0];
 		my $node = $_[1];
		my $ind = $_[2]->[0]; # no ind here: 
		my $step = $_[2]->[1];
		my $starting_node = $_[2]->[2];
		my $subtree_info = $_[2]->[3];
		my $mutations = $_[2]->[4];
		my $lifetime = $_[2]->[5];
		my $skip_stoppers_in_simulation = $_[2]->[6];
		if ($node eq $starting_node){
			#print " \n equality: ".$starting_node ->get_name."\t".$node ->get_name."\n";
			return;
		}
		if ($starting_node -> is_terminal){
			#print " \n terminal: ".$starting_node ->get_name."\n";
			return;
		}
		my %alive = %{$node->get_parent->get_generic("-alive")};  # a copy of existing hash! since it's one-dim, it is perfectly ok that it's only a shallow copy 
		my $nname = $node->get_name();
		my $nlength = $node->get_branch_length;
		my $startnname = $starting_node->get_name();
		my $depth = $_[3] - $starting_node->get_branch_length ;
		if (abs($depth - get_sequential_distance($starting_node, $node)) > 0.00001 ){
			print "Error! depth is strange : $depth ".get_sequential_distance($starting_node, $node)."\n";
		}
			
 		#print "depth $depth step $step bin ".(bin($depth,$step))."\n";
 		foreach my $site_index(keys %alive){
 			#print "$site_index is alive at ".$node->get_name()."!\n";
 				if (!$skip_stoppers_in_simulation && !($self->has_no_background_mutation($nname, $site_index))){
 				#	print "deleting $site_index ".$starting_node->get_name()."from alive: background mutation found at ".$node->get_name()."\n";
 					my $found;
 					foreach my $stopper (@{$self->{realdata}{subtree_info}{$startnname}{$site_index}{"stoppers"}}){
 						if ($nname eq $stopper->get_name()){
 							#print "Found stopper ".$stopper->get_name()." at $site_index ".$startnname."\n";
 							$found = 1;
 						}
 					}
 					unless ($found){
 						print "not_an_error_Iguess: ".$nname." is not in realdata stoppers list for ".$startnname." $site_index (in simulation we may pass a foreground mutation and see some previously hidden stoppers) \n";
 					}
 					delete $alive{$site_index};
 					next;
 				}
 				if ($lifetime && $self->{realdata}{subtree_info}{$startnname}{$site_index}{"maxdepth"} < $depth){
 				#	print "deleting $site_index from alive: allowed depth is ".$self->{realdata}{subtree_info}{$starting_node->get_name()}{$site_index}{"maxdepth"}.", node ".$node->get_name()."depth is $depth\n";
 					delete $alive{$site_index};
 					next;
 				}
 				$subtree_info->{$site_index}{"hash"}{bin($depth-$nlength/2,$step)}[1] += $nlength;
 				$subtree_info->{$site_index}{"totlengths"} += $nlength; 
		 		if ($mutations->{$site_index}{$nname}){ #take muts from hash! 
		 			$subtree_info->{$site_index}{"hash"}{bin($depth-$nlength/2,$step)}[0] += 1; # changed at 27.02.2017 (-halfbranch)
		 		#	print "Found a mutation at $site_index ".$node->get_name()."\n";
		 			$subtree_info->{$site_index}{"totmuts"} += 1;
		 			$subtree_info->{$site_index}{"totlengths"} -= $nlength/2;
		 	#		print "subtree_info totmuts for starting_node ".$starting_node->get_name()." node ".$node->get_name()." site $site_index is ".$subtree_info->{$site_index}{"totmuts"}."\n";
		 			delete $alive{$site_index};
		 		}	 	
		 		#	print "subtree maxdepth for ".$starting_node->get_name()." $site_index is ".$subtree_info->{$site_index}{"maxdepth"}."\n";
		 		if (!($subtree_info->{$site_index}{"maxdepth"}) || $subtree_info->{$site_index}{"maxdepth"} < $depth){
		 			$subtree_info->{$site_index}{"maxdepth"} = $depth;
		 		}
		 		#print "addded to 1\n";
 		}
 		$node->set_generic("-alive" => \%alive);
 	}

# 21.12 Differs from depth_groups_entrenchment_optimized_selector_alldepths in that it keeps both site and node in obs_hash
# (yes, that's one line) 

sub depth_groups_entrenchment_optimized_selector_alldepths_2 {
	my $self = shift;
	my $step = shift;
	my $restriction = shift;
	unless (defined $restriction) {$restriction = 0};

	my $root = $self ->{static_tree}-> get_root;
	my @array;
	my %hist;
	print "real data\n";
	
	my @group;
	if ($_[3]){
		@group = @{$_[3]};
	}
	else {
		@group = (1..$self->mylength());
	}
	

	my %closest_ancestors;
	$root->set_generic("-closest_ancestors" => \%closest_ancestors);
	my @args = ( $step, $root);
	$self->visitor_coat ($root, \@array,\&entrenchment_visitor,\&no_check,\@args,0, "overwrite"); # 27.02.2017 added "overwrite" flag
	my $debugnum = scalar keys %{$self ->{static_nodes_with_sub}};
	#print "News from depth..2: static_nodes_with_sub contains $debugnum keys\n";
	my $debugnum = scalar keys %{$self ->{static_subtree_info}};
	#print 	"News from depth..2: static_subtree_info contains $debugnum keys (nodes)\n";
	foreach my $ind (@group){
		foreach my $nod(@{$self ->{static_nodes_with_sub}{$ind}}){
			#print "nod is ".$nod." ref(nod) is ".ref($nod)."\n";
			my $node;
			if(ref($nod) eq "REF"){
				$node = ${$nod};
			}
			else {$node = $nod;}
			my $site_node = Textbits::concat($ind,$node->get_name());
			#print "site_node $site_node \n";
			my $total_muts = $self ->{static_subtree_info}{$node->get_name()}{$ind}{"totmuts"};
			my $total_length= $self ->{static_subtree_info}{$node->get_name()}{$ind}{"totlengths"};
			#print "depth ".$static_depth_hash{$ind}{$node->get_name()}."\n";
			#print " maxdepth is ".$self ->{static_subtree_info}{$node->get_name()}{$ind}{"maxdepth"}. " and restriction is $restriction\n";
			if ($self ->{static_subtree_info}{$node->get_name()}{$ind}{"maxdepth"} > $restriction){
				my %subtract_hash;
				
				if ($self ->{static_subtract_tallest}){
					my $tallest_tip = ${$self ->{static_hash_of_nodes}{$self ->{static_subtree_info}{$node->get_name()}{$ind}{"maxdepth_node"}}};	
					my @path_to_tallest_tip = @{$tallest_tip->get_ancestors()};	
					my $index_of_ancestor_node;
					my $path_length = $tallest_tip->get_branch_length; # todo: check if ancestors contain the node itself
				
					for(my $n = 0; $n < scalar @path_to_tallest_tip; $n++){		
						if ($path_to_tallest_tip[$n]->get_name eq $node->get_name){
							$index_of_ancestor_node = $n;
							last;
						}
						my $depth = get_sequential_distance($node, $path_to_tallest_tip[$n]);
						$subtract_hash{bin($depth,$step)} += $path_to_tallest_tip[$n]->get_branch_length;
						$path_length += $path_to_tallest_tip[$n]->get_branch_length;
					}
				}
				
				foreach my $bin (keys %{$self ->{static_subtree_info}{$node->get_name()}{$ind}{"hash"}}){
				#	$total_muts += $self ->{static_subtree_info}{$node->get_name()}{$ind}{"hash"}{$bin}[0];
				#	$total_length += $self ->{static_subtree_info}{$node->get_name()}{$ind}{"hash"}{$bin}[1];
					if ($self ->{static_subtract_tallest} && $subtract_hash{$bin}){
						$total_length -= $subtract_hash{$bin};
					}
				}	
			#	print $node->get_name()." ".$ind." TOTALS: $total_muts, $total_length, maxdepth ".$self ->{static_subtree_info}{$node->get_name()}{$ind}{"maxdepth"}."\n";
				#if ($total_length > 0 && $total_muts/$total_length < 0.005){
				if ($total_length > 0 && $total_muts > 0){ # 21.10 added total_muts > 0
			#	print "total muts $total_muts \n";
			#	print "site $ind node ".$node->get_name()." maxdepth ".$self ->{static_subtree_info}{$node->get_name()}{$ind}{"maxdepth"}."\n"; # commented out 16.09
					foreach my $bin (keys %{$self ->{static_subtree_info}{$node->get_name()}{$ind}{"hash"}}){
						#if ($total_length > 0 && $static_subtree_info{$node->get_name()}{$ind}{"hash"}{$bin}[1] > 0){ #there are some internal nodes with 0-length terminal daughter branches
						 if ($total_length > 0){ 
							my $local_length = $self ->{static_subtree_info}{$node->get_name()}{$ind}{"hash"}{$bin}[1];
							if ($self->{static_subtract_tallest} && $subtract_hash{$bin}){
								$local_length -= $subtract_hash{$bin};
							}
							$hist{$site_node}{$bin}[0] += $self ->{static_subtree_info}{$node->get_name()}{$ind}{"hash"}{$bin}[0]; #observed
							$hist{$site_node}{$bin}[1] += $total_muts*$local_length/$total_length; #expected						
						}
						if (!$hist{$site_node}{$bin}[0]){
							$hist{$site_node}{$bin}[0] += 0;
						}
						if (!$hist{$site_node}{$bin}[1]){
							$hist{$site_node}{$bin}[1] += 0;
						}
	# print "$bin,".$hist{$site_node}{$bin}[0].",".$hist{$site_node}{$bin}[1]."\n";  # commented out 16.09
				}
				}
			}
			
		}
	}	
	
	#foreach my $bin (sort {$a <=> $b} keys %hist){
	#	print " up to ".$bin*$step."\t".$hist{$bin}[0]."\t".$hist{$bin}[1]."\n";
	#}
	return %hist;	
}








#26.02 almost depth_groups_entrenchment_optimized_selector_alldepths_2, but prints only total nodecounts for each site_node
sub nodeselector {
	my $self = shift;
	my $step = $_[0];
	my $restriction = $_[1];

	my $root = $self->{static_tree}-> get_root;
	my @array;
	my %hist;
	
	my @group;
	if ($_[3]){
		@group = @{$_[3]};
	}
	else {
		@group = (1..$self->mylength());
	}
	 my $name = $_[4];

	my %closest_ancestors;
	$root->set_generic("-closest_ancestors" => \%closest_ancestors);
	my @args = ( $step, $root);
	$self->visitor_coat ($root, \@array,\&entrenchment_visitor,\&no_check,\@args,0);
	my %sitecounts;
	my %mutcounts;
	foreach my $ind (@group){
		foreach my $nod(@{$self->{static_nodes_with_sub}{$ind}}){
			my $node;
			if(ref($nod) eq "REF"){
				$node = ${$nod};
			}
			else {$node = $nod;}
			my $site_node = Textbits::concat($ind, $node->get_name());
			my $total_muts = $self->{static_subtree_info}{$node->get_name()}{$ind}{"totmuts"};
			my $total_length = $self->{static_subtree_info}{$node->get_name()}{$ind}{"totlengths"};
			
	#print "depth ".$static_depth_hash{$ind}{$node->get_name()}."\n";
			if ($self->{static_subtree_info}{$node->get_name()}{$ind}{"maxdepth"} > $restriction){
				
				my %subtract_hash;
				
				if ($self->{static_subtract_tallest}){
					#print "just checking ".$static_subtree_info{$node->get_name()}{$ind}{"maxdepth_node"}."\n";
					my $tallest_tip = ${$self ->{static_hash_of_nodes}{$self->{static_subtree_info}{$node->get_name()}{$ind}{"maxdepth_node"}}};	
					my @path_to_tallest_tip = @{$tallest_tip->get_ancestors()};	
					my $index_of_ancestor_node;
					my $path_length = $tallest_tip->get_branch_length; # todo: check if ancestors contain the node itself
				
					for(my $n = 0; $n < scalar @path_to_tallest_tip; $n++){		
						if ($path_to_tallest_tip[$n]->get_name eq $node->get_name){
							$index_of_ancestor_node = $n;
							last;
						}
						my $depth = get_sequential_distance($node, $path_to_tallest_tip[$n]);
						$subtract_hash{bin($depth,$step)} += $path_to_tallest_tip[$n]->get_branch_length;
						$path_length += $path_to_tallest_tip[$n]->get_branch_length;
					}
				}
				
				foreach my $bin (keys %{$self->{static_subtree_info}{$node->get_name()}{$ind}{"hash"}}){
				#	$total_muts += $self->{static_subtree_info}{$node->get_name()}{$ind}{"hash"}{$bin}[0];
				#	$total_length += $self->{static_subtree_info}{$node->get_name()}{$ind}{"hash"}{$bin}[1];
					if ($self->{static_subtract_tallest} && $subtract_hash{$bin}){
						$total_length -= $subtract_hash{$bin};
					}
				}	
				#print $node->get_name()." ".$ind." TOTALS: $total_muts, $total_length, maxdepth ".$static_subtree_info{$node->get_name()}{$ind}{"maxdepth"}."\n";
				#if ($total_length > 0 && $total_muts/$total_length < 0.005){
				if ($total_length > 0 && $total_muts > 0){ # 21.10 added total_muts > 0
				#print "total muts $total_muts \n";
				my $md = $self->{static_subtree_info}{$node->get_name()}{$ind}{"maxdepth"};
				if ($md > 50){
					$sitecounts{50}{$ind} = 1;
				}
				if ($md > 100){
					$sitecounts{100}{$ind} = 1;
				}
				if ($md > 150){
					$sitecounts{150}{$ind} = 1;
				}				
				print "site $ind node ".$node->get_name()." maxdepth ".$self->{static_subtree_info}{$node->get_name()}{$ind}{"maxdepth"};
				my %totalcounts; # 26.02 counting nodes in analysis	
					foreach my $bin (keys %{$self->{static_subtree_info}{$node->get_name()}{$ind}{"hash"}}){
						#if ($total_length > 0 && $static_subtree_info{$node->get_name()}{$ind}{"hash"}{$bin}[1] > 0){ #there are some internal nodes with 0-length terminal daughter branches
						 if ($total_length > 0){ 
							my $local_length = $self->{static_subtree_info}{$node->get_name()}{$ind}{"hash"}{$bin}[1];
							if ($self->{static_subtract_tallest} && $subtract_hash{$bin}){
								$local_length -= $subtract_hash{$bin};
							}
							$hist{$site_node}{$bin}[0] += $self->{static_subtree_info}{$node->get_name()}{$ind}{"hash"}{$bin}[0]; #observed
							$hist{$site_node}{$bin}[1] += $total_muts*$local_length/$total_length; #expected
							$totalcounts{$site_node} += $self->{static_subtree_info}{$node->get_name()}{$ind}{"hash"}{$bin}[0];	# 26.02 counting nodes in analysis
							
							if ($md > 50){
								$mutcounts{50} += $self->{static_subtree_info}{$node->get_name()}{$ind}{"hash"}{$bin}[0];	# 26.02 counting nodes in analysis
							}	
							if ($md > 100){
								$mutcounts{100} += $self->{static_subtree_info}{$node->get_name()}{$ind}{"hash"}{$bin}[0];	# 26.02 counting nodes in analysis
							}
							if ($md > 150){
								$mutcounts{150} += $self->{static_subtree_info}{$node->get_name()}{$ind}{"hash"}{$bin}[0];	# 26.02 counting nodes in analysis
							}							
						}
						if (!$hist{$site_node}{$bin}[0]){
							$hist{$site_node}{$bin}[0] += 0;
						}
						if (!$hist{$site_node}{$bin}[1]){
							$hist{$site_node}{$bin}[1] += 0;
						}

				}
					 print " ".$totalcounts{$site_node}."\n";
				}
			}
			
		}
	}	
	
		my $sites50 = scalar keys %{$sitecounts{50}};
		my $sites100 = scalar keys %{$sitecounts{100}};
		my $sites150 = scalar keys %{$sitecounts{150}};
		print $name." ".$sites50." ".$mutcounts{50}." ".$sites100." ".$mutcounts{100}." ".$sites150." ".$mutcounts{150}." "."\n";
	
	return %hist;	
}

sub my_median{
my @values = @{$_[0]};	
my $median;
my $mid = int @values/2;
my @sorted_values = sort { $a <=> $b } @values;
if (@values % 2) {
    $median = $sorted_values[ $mid ];
} else {
    $median = ($sorted_values[$mid-1] + $sorted_values[$mid])/2;
} 
return $median;
}






## for hist->interval->site_index
sub hist_median_group {
	my @pre_hist = @{$_[0]};
	my @group = @{$_[1]};
	
	my @hist;
	foreach my $ind(@group){
		for (my $interval = 0; $interval < scalar @pre_hist; $interval++){
			$hist[$interval] += $pre_hist[$interval]->[$ind];
		}
	}
	
	return hist_median(\@hist);
}

sub print_hist {
	my @hist = @{$_[0]};

	my $counter = 0;
	print "\n";
	for (my $interval = 0; $interval < scalar @hist; $interval++){
			print $hist[$interval]."\t";
			$counter+=$hist[$interval];
		}
	print "\n";

}




sub median_difference{
	my @a1 = @{$_[0]};
	my @a2 = @{$_[1]};
	my $stat = Statistics::Descriptive::Full->new();
	$stat->add_data(@a1);
	my $median1 = $stat->median();

	$stat = Statistics::Descriptive::Full->new();
	$stat->add_data(@a2);
	my $median2 = $stat->median();
	my $diff = $median1-$median2;
	print "median1 $median1 median2 $median2 diff $diff\n";
	return $diff;
}

# one hist for one ancestor aa in one site
# the chosen one (used by r scripts for drawing plots)
# circles, not rings! ()
# mutations at the ends of branches !

sub egor_smart_site_entrenchment {
	my $self = shift;
	my $verbose = shift;
	my $step = 1;
	my $root = $self->{static_tree}-> get_root;
	my $file = File::Spec->catfile($self->{static_output_base}, "egor_smart_".$self->{static_protein}.".csv");
	open my $plotcsv, ">$file" or die "Cannot create $file \n";
	my @array;
	print $plotcsv "radius,site,node,density,cum_muts,cum_length\n";
	my $hash_ready;
	if (exists $self->{static_ring_hash}){
		warn "Static_ring_hash is ready, egor_smart_site_entrenchment won't change it\n";
		$hash_ready = 1;
	}

	for (my $ind = 1; $ind < $self->mylength(); $ind++){
		if ($verbose) {print "$ind\n"};
		foreach my $nod(@{$self->{static_nodes_with_sub}{$ind}}){
			my $node = ${$nod};
			if (!$node->is_terminal){
			if ($verbose) {print $node->get_name()."\n"};
			my %hist;
			#my $ancestor = ${$static_subs_on_node{$node->get_name()}}{$ind}->{"Substitution::ancestral_allele"};
			my @args = ($ind, $step, $node);
			unless ($hash_ready) {$self->my_visit_depth_first ($node, \@array,\&update_ring,\&has_no_mutation,\@args,0);} #visitor_coat cannot be used inside a cycle; we still do not want to mess the hash up, so we check for its existance before this cycle
			
			my $cumulative_muts;
			my $cumulative_length;
			my @sorted_keys = sort {$a <=> $b} (keys %{$self->{static_ring_hash}{$ind}{$node->get_name()}});
			foreach my $bin (@sorted_keys){
				my $muts_in_bin = $self->{static_ring_hash}{$ind}{$node->get_name()}{$bin}[0];
				my $length_of_bin = $self->{static_ring_hash}{$ind}{$node->get_name()}{$bin}[1];
				$cumulative_muts += $muts_in_bin;
				$cumulative_length += $length_of_bin;
				if ($verbose) {print "bin $bin observed ".$self->{static_ring_hash}{$ind}{$node->get_name()}{$bin}[0]." totmut $cumulative_muts muts in bin $muts_in_bin totlen $cumulative_length bin len $length_of_bin\n"};
				if ($cumulative_length > 0){ #there are some internal nodes with 0-length terminal daughter branches
					#	$hist{$bin}{$ancestor} += $self->{static_ring_hash}{$ind}{$node->get_name()}{$bin}[0]/$self->{static_ring_hash}{$ind}{$node->get_name()}{$bin}[1]; #density

					$hist{$bin} = $cumulative_muts/$cumulative_length; #density
				}
				else {
					$hist{$bin} = 0;
				}
				if ($muts_in_bin > 0){	
					print $plotcsv "$bin,$ind,".$node->get_name().",".$hist{$bin}.",".$cumulative_muts.",".$cumulative_length."\n";
				}
				
			}
			

		}
		}
		
	}

	$self->printFooter($plotcsv);
	close $plotcsv;
	
}

# last egor plots sub. NOT the chosen one (egor_h1.csv are printed by some other method (discovered by comparison))
# honest rings, can be used for mann-kendall analysis
# ! mutations at the ends of branches !

sub egor_diff_rings_site_entrenchment {
	my $self = shift;
	my $step = shift;
	my $cumulative = shift;
	print ($self->{static_protein}."\n");
	my $root = $self->{static_tree}-> get_root;
	my $tag;
	if ($cumulative) { $tag = "cumulative"; }
	else { $tag = "diff_rings"; }
	my $file = File::Spec->catfile($self->{static_output_base}, "egor_".$tag."_".$self->{static_protein}.".csv");
	open my $plotcsv, ">$file" or die "Cannot create $file \n";
	my @array;
	print $plotcsv "radius,site,node,density,cum_muts,cum_length\n";
	
	my $hash_ready;
	if (exists $self->{static_ring_hash}){
		$hash_ready = 1;
		warn "Static_ring_hash is ready, egor_diff_rings_site_entrenchment won't change it\n";
	}

	for (my $ind = 1; $ind < $self->mylength(); $ind++){
		#print "$ind\n";
		foreach my $nod(@{$self->{static_nodes_with_sub}{$ind}}){
			my $node = ${$nod};
			if (!$node->is_terminal){
			#print $node->get_name()."\n";
			my %hist;
			#my $ancestor = ${$static_subs_on_node{$node->get_name()}}{$ind}->{"Substitution::ancestral_allele"};
			my @args = ($ind, $step, $node);
			unless ($hash_ready) {$self->my_visit_depth_first ($node, \@array,\&update_ring,\&has_no_mutation,\@args,0);} #visitor_coat cannot be used inside a cycle; we still do not want to mess the hash up, so we check for its existance
			my $cumulative_muts;
			my $cumulative_length;
			my @sorted_keys = sort {$a <=> $b} (keys %{$self->{static_ring_hash}{$ind}{$node->get_name()}});
			foreach my $bin (@sorted_keys){
				#print "bin $bin observed ".$self->{static_ring_hash}{$ind}{$node->get_name()}{$bin}[0]." totmut $cumulative_muts totlen $cumulative_length\n";
				my $muts_in_bin = $self->{static_ring_hash}{$ind}{$node->get_name()}{$bin}[0];
				my $length_of_bin = $self->{static_ring_hash}{$ind}{$node->get_name()}{$bin}[1];
				$cumulative_muts += $muts_in_bin;
				$cumulative_length += $length_of_bin;
				if ($cumulative_length > 0){ #there are some internal nodes with 0-length terminal daughter branches
					#	$hist{$bin}{$ancestor} += $static_ring_hash{$ind}{$node->get_name()}{$bin}[0]/$static_ring_hash{$ind}{$node->get_name()}{$bin}[1]; #density
					$hist{$bin} = $cumulative_muts/$cumulative_length; #density
				}
				else {
					$hist{$bin} = 0;
				}
				if ($muts_in_bin > 0 || $bin == $sorted_keys[0] || $bin == $sorted_keys[-1]){	
				#if ($muts_in_bin > 0){	
					print $plotcsv "$bin,$ind,".$node->get_name().",".$hist{$bin}.",".$cumulative_muts.",".$cumulative_length."\n";
					unless ($cumulative){
						$cumulative_muts = 0;
						$cumulative_length = 0;
					}
				}
				
			}
			

		}
		}
	}
	
	close $plotcsv;
	
	
}


# decorator for my_visit_depth_first, checks for existance of corresponding hash and prevents unintentional changes in it (or deletes it and overwrites)
# must not be used in loop context!
## !!! does not work as expected, should not be used for preventing changes in pre-existing data. 
sub visitor_coat {
		my $self = shift;
		my $node = $_[0];
		my @array = @{$_[1]};
		my $action_callback = $_[2];
		my $check_callback = $_[3];
		my $callback_args = $_[4];
		my $depth = $_[5];
		my $overwrite = $_[6];
		
		my $visitor_name = sub_name($action_callback);
		
		if ($visitor_name eq "update_ring"){
			if (exists $self->{static_ring_hash}){
				if ($overwrite){
					delete $self->{static_ring_hash};
				}
				else {return;}	
			}	
		}
		elsif ($visitor_name eq "entrenchment_visitor" || $visitor_name eq "lrt_visitor"){
			if (exists $self->{static_subtree_info}){
				if ($overwrite){
					delete $self->{static_subtree_info};
				}
				else {return;} # added at 27.02.2017 or we will add to already existing structure again.
			# commented out at 27.09.2016	
			#		foreach my $node(keys $self->{static_subtree_info}){
			#				foreach my $site(keys $self->{static_subtree_info}{$node}){
			#					if ($visitor_name eq "entrenchment_visitor"){
			#						if (exists $self->{static_subtree_info}{$node}{$site}{"hash"}){
			#							if ($overwrite){
		#									delete $self->{static_subtree_info}{$node}{$site}{"hash"};
	#										delete $self->{static_subtree_info}{$node}{$site}{"maxdepth"};
	#										delete $self->{static_subtree_info}{$node}{$site}{"maxdepth_node"};
	#									}
	#									else {return;}
	#								}
	#									
	#							}
	#							elsif ($visitor_name eq "lrt_visitor"){
	#								if (exists $self->{static_subtree_info}{$node}{$site}{"lrt"}){
	#									if ($overwrite){
	#										delete $self->{static_subtree_info}{$node}{$site}{"lrt"};
	#									}
	#									else {return;}
	#								}
	#							}
	#						}
	#				}
			}
		}
		else {
			print "Warning: visitor_coat cannot perform any check for $visitor_name subroutine. Launching my_visit_depth_first without checking for previous launches\n";
		}
		$self->my_visit_depth_first($node, \@array, \&$action_callback, \&$check_callback, $callback_args, $depth);
}

 sub my_visit_depth_first {
		my $self = shift;
		my $node = $_[0];
		my @array = @{$_[1]};
		my $action_callback = $_[2];
		my $check_callback = $_[3];
		my $callback_args = $_[4];
		my $depth = $_[5];
		push @array, $node;
		my $len = $node -> get_branch_length;
		$depth += $len;
		&$action_callback($self, $node, $callback_args, $depth);
		if (! $node->is_terminal && &$check_callback($self, $node, $callback_args)){
			my $i = 0;
			while($node->get_child($i)){
				@array = $self -> my_visit_depth_first($node->get_child($i), \@array, \&$action_callback, \&$check_callback, $callback_args, $depth);
				$i++;
			}
		}	
		$node = pop @array;
#print $node->get_name()."\t".$depth."\n";
		$depth -= $len;
		return @array;
    }
  
sub print_xparr {
	my $self = shift;
	my $tag = shift;
	my $root = $self ->{static_tree}-> get_root;
	my @array;

	my $dir = File::Spec -> catdir($self -> {static_output_base}, "xparr");
	make_path($dir);
	my $filepath = File::Spec -> catfile($dir, $self->{static_protein}.".xparr");
	open FILE, ">$filepath" or die "Cannot open $filepath\n";
	print FILE "child\tparent\tlength\tsyn\tnonsyn\tnonsyn\tphenotypes:Wilke\n";
	close FILE;
	$self->my_visit_breadth_first($root,\@array, \&xparr_callback, $filepath);
}  

sub xparr_callback {
	my $self = shift;
	my $node = shift;
	my $file = shift;
	open FILE, ">>$file" or die "Cannot open $file\n";
	my $nsynsubs = join(",", keys %{$self->{static_subs_on_node}{$node->get_name}});
	my $synsubs = join(",", keys %{$self->{static_background_subs_on_node}{$node->get_name}});
	my $parent_name;
	unless ($node->is_root) {$parent_name = $node->get_parent->get_name;}
	print join("\t", $node->get_name,$parent_name,$node -> get_branch_length,$synsubs,$nsynsubs)."\n";
	print FILE join("\t", $node->get_name,$parent_name,$node -> get_branch_length,$synsubs,$nsynsubs,$nsynsubs)."\n";
	close FILE;
	
}
   
 sub my_visit_breadth_first {
		my $self = shift;
		my $node = $_[0];
		my @array = @{$_[1]};
		my $action_callback = $_[2];
		my $callback_args = $_[3];

		&$action_callback($self, $node, $callback_args);
		if (!$node->is_terminal){
			my $i = 0;
			while($node->get_child($i)){
				push @array, $node->get_child($i);
				$i++;
			}
		}	
		$node = shift @array;
		if (defined $node){
			@array = $self -> my_visit_breadth_first($node, \@array, \&$action_callback, $callback_args);
	#print $node->get_name()."\t".$depth."\n";
			return @array;
		}
    }
    
    # old version, does not account for  mutations of the other type.
    # Used in update ring (for making plots) (that's ok, because plotting subroutines use newer has_no_mutation method to account for background mutations )
 	sub has_no_same_type_mutation { #has_no_same_type_mutation
 		my $self = $_[0];
 		my $node = $_[1];
 		my $site_index = $_[2]->[0];
 		my $starting_node = $_[2]->[2];
 		
 		if ($node eq $starting_node){
 			return 1;
 		}
 		if (${ $self->{static_subs_on_node}{$node->get_name()}}{$site_index}){
 			return 0;
 		}
 		else {
 			return 1;
 		}
 	}  
 	
 	
 	# also accounts for mutations of the other type (synonimous for nsyn and non-synonimous for syn)
 	sub has_no_mutation{
 		my $self = $_[0];
 		my $node = $_[1];
 		my $site_index = $_[2]->[0];
 		my $starting_node = $_[2]->[2];
 		my $comparator = $self->{static_comparator};

 		if ($node eq $starting_node){
 			return 1;
 		}
 		if (${$self->{static_subs_on_node}{$node->get_name()}}{$site_index}){
 			return 0;
 		}
 		
 		if ($self->{static_state} eq "nsyn"){
 			if ($comparator->is_neighbour_changing(${$self->{static_background_subs_on_node}{$node->get_name()}}{$site_index}, 1) == 1){
 				return 0;
 			}
 			else {
 				return 1;
 			}
 		}
 		else {
 			if (${$self->{static_background_subs_on_node}{$node->get_name()}}{$site_index}){
 				return 0;
 			}
 			else {
 				return 1;
 			}
 		}
 	}  
 
 	
 	sub has_no_background_mutation {
 		my $self = $_[0];
 		my $nodename = $_[1];
 		my $site_index = $_[2];
 		my $comparator = $self->{static_comparator};
 		if ($self->{static_state} eq "nsyn"){
 			if ($self->{static_skip_stoppers}){
 				return 1;
 			}
 			else {
 				if ($comparator->is_neighbour_changing(${$self->{static_background_subs_on_node}{$nodename}}{$site_index}, 1) == 1){
 					return 0;
 				}
 				else {
 					return 1;
 				}
 			}
 		}
 		else {
 			if (${$self->{static_background_subs_on_node}{$nodename}}{$site_index}){
 				return 0;
 			}
 			elsif ($self->{static_no_neighbour_changing} && $comparator->is_neighbour_changing(${$self->{static_subs_on_node}{$nodename}}{$site_index}, 1) == 1){ # added at 10.01.2017
 				return 0;
 			}
 			else {
 				return 1;
 			}
 		}
 	}
 	

 	

 	
 	#sub max ($$) { $_[$_[0] < $_[1]] }
 	
 	sub update_ring {
 		my $self = $_[0];
 		my $node = $_[1];
		my $site_index = $_[2]->[0];
		my $step = $_[2]->[1];
		my $starting_node = $_[2]->[2];
		my $depth = $_[3] - $starting_node->get_branch_length ;
		if ($node eq $starting_node){
			#print " \n equality: ".$starting_node ->get_name."\t".$node ->get_name."\n";
			return;
		}
		if ($starting_node -> is_terminal){
			#print " \n terminal: ".$starting_node ->get_name."\n";
			return;
		}
 		#print "depth $depth step $step bin ".(bin($depth,$step))."\n";
 		if (!($self->has_no_same_type_mutation($_[1], \@{$_[2]}))){ #has_no_same_type_mutation
 			$self->{static_ring_hash}{$site_index}{$starting_node->get_name()}{bin($depth,$step)}[0] += 1;
 			#$self->{static_ring_hash}{$site_index}{$starting_node->get_name()}{bin($depth,$step)}[1] -= ($node->get_branch_length)/2 # 19.09.2016 mutations happen in the middle of a branch; todo
 			#print "addded to 0\n";
 		}
 		#my $newdepth = $self->{static_distance_hash}{$anc_node->get_name()}{$node->get_name()} - ($node->get_branch_length)/2; # 19.09.2016 mutations happen in the middle of a branch; todo
 		$self->{static_ring_hash}{$site_index}{$starting_node->get_name()}{bin($depth,$step)}[1] += $node->get_branch_length;
 		#print "addded to 1\n";
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
 	
 # since 01.06.2017 - instead of distance matrix	
sub	set_timestamps {
 		my $self = shift;
		my $tree = $self->{static_tree};
		$tree->visit_breadth_first(
			-in   => sub{
				my $node=shift;
				if($node->is_root){
					$node->set_generic('time' => 0);
				}else{
					my $pnode=$node->get_parent;
					my $time=$pnode->get_generic('time');
					$time+=$node->get_branch_length;
					$node->set_generic('time' => $time);
				}
			}
		);	
	}
	
sub get_sequential_distance {
	my $ancnode = shift;
	my $node = shift;
	return $node->get_generic('time') - $ancnode->get_generic('time');
}

 	
 	
 	
 	
   # track_tallest is needed for finding longest path in the subtree and subtracting its length
   # Added at 08.10 for testing whether this will improve correspondence between simulation_observed and simulation_expected.
   # stopped using it at 15.09.2016    
   # 15.09.2016 version: halves of branches with foreground mutations are trimmed 
   # 19.09.2016: corrected
   	sub entrenchment_visitor {
 		my $self = shift;
 		my $node = $_[0];
 		my $step = $_[1]->[0];
 		my $subtract_tallest = $self->{static_subtract_tallest};
 		my $no_neighbour_changing = $self->{static_no_neighbour_changing};
 		my $no_leaves = $self->{static_no_leaves};
 		my $include_tips =  $self->{static_include_tips};
		if (!$node->is_root){
		my $nname = $node->get_name();
		my $nlength = $node->get_branch_length;
		my %closest_ancestors = %{$node->get_parent->get_generic("-closest_ancestors")}; # closest_ancestors: ancestor mutation node for this node, key is a site number
		my $comparator = $self->{static_comparator};
		if (%closest_ancestors){
			## pasted here 27.02.2017
			my @ancestors = keys %closest_ancestors;	
			foreach my $site_index(@ancestors){
				if (!($self->has_no_background_mutation($nname, $site_index))){
					## copy-pasted at 27.02.2017 (now we do not add branches with bkgr mutations to square (or maxdepth))
					my $anc_node = $closest_ancestors{$site_index};
					push @{$self->{static_subtree_info}{$anc_node->get_name()}{$site_index}{"stoppers"}}, $node;
					##
					delete $closest_ancestors{$site_index};
				}
			}	
			##
			foreach my $site_index(keys %closest_ancestors){ 
				my $anc_node = $closest_ancestors{$site_index};
				my $ancname = $anc_node->get_name();
				my $depth = get_sequential_distance($anc_node,$node);
				my $halfdepth = get_sequential_distance($anc_node,$node) - ($nlength)/2; 
				$self->{static_subtree_info}{$ancname}{$site_index}{"hash"}{bin($halfdepth,$step)}[1] += $nlength;
				$self->{static_subtree_info}{$ancname}{$site_index}{"totlengths"} += $nlength; # for lambda computation only!		
			#	print "anc ".$anc_node->get_name()." site ".$site_index." node ".$node->get_name()." depth $depth bin ".bin($depth,$step)." branchlength ".$node->get_branch_length."\n";
				# commented out at 27.02.2017 and copy-pasted higher
				my $current_maxdepth = $self->{static_subtree_info}{$ancname}{$site_index}{"maxdepth"};
				if ($current_maxdepth < $depth){
						$self->{static_subtree_info}{$ancname}{$site_index}{"maxdepth"} = $depth;
						if ($subtract_tallest){
							$self->{static_subtree_info}{$ancname}{$site_index}{"maxdepth_node"} = $nname;
						}
				}					
			}
			
			## deleted from here 27.02.2017
		}
		
		# for nodes with mutations
		foreach my $site_index(keys %{$self->{static_subs_on_node}{$nname}}){
			
			if ($closest_ancestors{$site_index}){
				my $anc_node = $closest_ancestors{$site_index};
				my $ancname = $anc_node->get_name();
				my $halfdepth = get_sequential_distance($anc_node,$node) - ($nlength)/2; #19.09.2016 
			#	my $fulldepth = get_sequential_distance($anc_node,$node); #19.09.2016 
			#	print " ancestor ".$anc_node->get_name(). " node ".$node->get_name()." depth $depth\n";
				if (!$no_neighbour_changing || ($no_neighbour_changing && ! $comparator->is_neighbour_changing($self->{static_subs_on_node}{$nname}{$site_index}, 1))){
					if (!$no_leaves || ($no_leaves && !($node->is_terminal()))){
						$self->{static_subtree_info}{$ancname}{$site_index}{"hash"}{bin($halfdepth,$step)}[0] += 1; #19.09.2016 
						$self->{static_subtree_info}{$ancname}{$site_index}{"totmuts"} += 1; #21.12.2016
					}
				}
			#	$self->{static_subtree_info}{$ancname}{$site_index}{"hash"}{bin($halfdepth,$step)}[1] += $nlength; # #29.08 halves are back again 19.09.2016  15.09.2016 version: halves of branches with foreground mutations are trimmed (the only thing I changed here) 
			#	$self->{static_subtree_info}{$ancname}{$site_index}{"hash"}{bin($fulldepth,$step)}[1] -= $nlength; #19.09.2016 we added this length before, but shouldn't have done it
				$self->{static_subtree_info}{$ancname}{$site_index}{"totlengths"} -= ($nlength)/2; # for lambda computation only!		
			#	print "mutation! anc ".$anc_node->get_name()." site ".$site_index." node ".$node->get_name()." depth $depth bin ".bin($depth,$step)." branchlength ".$node->get_branch_length."\n";
				
			}
			$closest_ancestors{$site_index} = $node;
			
		}
		
		$node->set_generic("-closest_ancestors" => \%closest_ancestors);
		}
#!
 	}

     	sub reversion_ratio_visitor {
 		my $self = shift;
 		my $node = $_[0];
 		my $step = $_[1]->[0];
 		my $no_neighbour_changing = $self->{static_no_neighbour_changing};
		my $myCodonTable   = Bio::Tools::CodonTable->new();
		if (!$node->is_root){
		my $nname = $node->get_name();
		my $nlength = $node->get_branch_length;
		my %closest_ancestors = %{$node->get_parent->get_generic("-closest_ancestors")}; # closest_ancestors: ancestor mutation node for this node, key is a site number
		my $comparator = $self->{static_comparator};
		if (%closest_ancestors){
			my @ancestors = keys %closest_ancestors;
			foreach my $site_index(@ancestors){
				my $anc_node = $closest_ancestors{$site_index};		
				if (!($self->has_no_background_mutation($nname, $site_index))){
					my $anc_node = $closest_ancestors{$site_index};
					push @{$self->{static_subtree_info}{$anc_node->get_name()}{$site_index}{"stoppers"}}, $node;
					delete $closest_ancestors{$site_index};
				}
			}	
		}
		# for nodes with mutations
		foreach my $site_index(keys %{$self->{static_subs_on_node}{$nname}}){
			if ($closest_ancestors{$site_index}){
				my $anc_node = $closest_ancestors{$site_index};
				my $ancname = $anc_node->get_name();
				my $halfdepth = get_sequential_distance($anc_node,$node) - ($nlength)/2; #19.09.2016 
				print get_sequential_distance($anc_node,$node)." ".(($nlength)/2)."\n";
				if (!$no_neighbour_changing || ($no_neighbour_changing && ! $comparator->is_neighbour_changing($self->{static_subs_on_node}{$nname}{$site_index}, 1))){
						my $allele0 = $self->{static_subs_on_node}{$ancname}{$site_index}{"Substitution::ancestral_allele"};
						my $allele2 = $self->{static_subs_on_node}{$nname}{$site_index}{"Substitution::derived_allele"};
						if ($self->{static_state} eq "nsyn"){
							$allele0 = $myCodonTable->translate($allele0);
							$allele2 = $myCodonTable->translate($allele2);
						}
						$self->{static_subtree_info}{$ancname}{$site_index}{"repertoire"}{$allele2}  =1;
						if ($allele0 eq $allele2){
							$self->{static_subtree_info}{$ancname}{$site_index}{"reversions"}{bin($halfdepth,$step)}[0] += 1; #reversion
						}
						else {
							$self->{static_subtree_info}{$ancname}{$site_index}{"reversions"}{bin($halfdepth,$step)}[1] += 1; #not a reversion
						}
				}				
			}
			$closest_ancestors{$site_index} = $node;
		}
		$node->set_generic("-closest_ancestors" => \%closest_ancestors);
		}
#!
 	}

   
   
   	sub lrt_visitor {
   		my $self = shift;
 		my $node = $_[0];

		if (!$node->is_root){
			my %closest_ancestors = %{$node->get_parent->get_generic("-closest_ancestors")};
			my $nname = $node->get_name();
			if (%closest_ancestors){
				foreach my $site_index(keys %closest_ancestors){ 
					my $anc_node = $closest_ancestors{$site_index};
					my $depth =  get_sequential_distance($anc_node,$node);
					$self->{static_subtree_info}{$anc_node->get_name()}{$site_index}{"lrt"}{$nname}[1] = $depth-($node->get_branch_length); #23.03 t_branch_start
					$self->{static_subtree_info}{$anc_node->get_name()}{$site_index}{"lrt"}{$nname}[2] = $depth; #23.03 t_branch_end					
				}
				
				my @ancestors = keys %closest_ancestors;	
				foreach my $site_index(@ancestors){
					if (!($self->has_no_background_mutation($nname, $site_index))){
						delete $closest_ancestors{$site_index};
					}
				}	
			}
			
			foreach my $site_index(keys %{$self->{static_subs_on_node}{$nname}}){
				if ($closest_ancestors{$site_index}){
					my $anc_node = $closest_ancestors{$site_index};
					$self->{static_subtree_info}{$anc_node->get_name()}{$site_index}{"lrt"}{$nname}[0] = 1;
				}
				$closest_ancestors{$site_index} = $node;
			}
			
			$node->set_generic("-closest_ancestors" => \%closest_ancestors);
		}

 	}
   
   sub predefined_groups_and_names {
 		my $self = shift;
 		my $prot = $self->{static_protein};
 		return Groups::get_predefined_groups_and_names_for_protein($prot, $self->mylength());
   }
   
   sub user_groups_and_names {
 		my $self = shift;
		my $filename = shift;
 		my $prot = $self->{static_protein};
		my $filepath = File::Spec->catdir(getcwd(), "data", "groups", $filename);
 		return Groups::get_user_groups_and_names_for_protein($prot, $self->mylength(),$filepath);
   }
   
   sub fake_predefined_groups_and_names {
   	 	my $self = shift;
 		my $prot = $self->{static_protein};
 		my $state = $self->{static_state};
 		return Groups::get_fake_predefined_groups_and_names_for_protein($prot, $self->mylength(), $state);
   }

	sub protein_no_group {
 		my $self = shift;
 		my $prot = $self->{static_protein};
 		return Groups::get_no_groups_for_protein($prot, $self->mylength());
   }
   
 # changed at 21.09.2016  
   sub bin {
   	my $depth = $_[0];
   	my $step = $_[1];
   	
   	my $bin = int($depth/$step);
   	if ((int($depth/$step) == $depth/$step && $depth != 0) || $step == 1 || $step == 0.5){ # 0 goes to 0 bin, if step is 0.5 or 1, and to 1 bin otherwise
   		$bin -= 1;
   	}
   	return $bin+1;
   }
   
   sub bin_new {
   	   	my $depth = $_[0];
   		my $step = $_[1];
   	#	print "Warning: not using step option for binning (if you did not explicitly set step to 1, it will cause havoc)\n";
   		return $depth;
   }
   
   sub array_mean {
	return sum(@_)/@_;
   }

   

    

1;
