#!/usr/bin/perl

package Aspens;

use strict;
use Bio::Phylo::IO;
use TreeUtils::Phylo::PhyloUtils qw(remove_zero_branches);
use DistanceFinder qw(get_mrcn calc_true_patristic_distance node_distance);
use Bio::Tools::CodonTable;
use TreeUtils::Phylo::FigTree;
use Bio::SeqIO;
use Data::Dumper;
use Bit::Vector;
use Try::Tiny;
use List::Util qw(sum shuffle);
use Const::Fast;
use Switch;
use Statistics::Basic qw(:all);
use Statistics::TTest;
use Statistics::Descriptive;
use Storable qw(store retrieve lock_retrieve);

use Class::Struct;
use IO::Handle;
use Cwd qw(abs_path cwd getcwd);
use File::Spec;

$| = 1;

sub print_scheme{
	my $self = shift;
	my $ind  = shift;
	my $filepath  = $self->graphpath($ind);
	$self->set_same_ancestor_subs();
	my %hash = %{$self->{static_same_ancestor_subs}};
	my $myCodonTable   = Bio::Tools::CodonTable->new();
	
	open OUT, ">$filepath" or die "Cannot open $filepath $!";
	foreach my $anc (keys %{$hash{$ind}}){
		print OUT "Ancestor:".$myCodonTable->translate($anc)."(".$anc.")\n";
		print OUT "Derived:";
		my @subs;
		foreach my $sub (@{$hash{$ind}->{$anc}}){
			push @subs, $sub->{"Substitution::derived_allele"}."|".$sub->{"Substitution::probability"};
		}
		print OUT join(",", @subs)."\n";
	}
	close OUT;
}


#global analysis

sub global_stats{
	my $mutmap = shift;
	my $step = shift;
	my $iterate = shift;
	
	my $tree = $mutmap->{static_tree};
	my %subs_on_node = %{$mutmap->{static_subs_on_node}};
	my %nodes_with_sub = %{$mutmap->{static_nodes_with_sub}};

	my $outfile = File::Spec->catfile($mutmap->{static_output_base}, $mutmap->{static_protein}."_".$mutmap->{static_state}."_global_median_statistics");	
	
	open OUT, ">$outfile" or die "cannot create output file $outfile: $!";
	my @bootstrap_median_diff;
	my @bins = ();
	my @array = (1..$mutmap->{static_length});

	foreach my $ind(@array){
			print OUT $ind."\n";
			my %distr = find_all_distances_probs($mutmap, $ind);
print (Dumper(\%distr));
			my @site_bins = distr_to_stathist_probs(\%distr, $step);
			
			if (defined $site_bins[0]->[1] && defined $site_bins[1]->[1]){
				for( my $interval = 0; $interval < scalar @{$site_bins[0]}; $interval++){
					$bins[0]->[$interval] += $site_bins[0]->[$interval];
#					print OUT $site_bins[0]->[$interval];
#					print OUT "\t";
#					print OUT $bins[0]->[$interval];
#					print OUT "\t";
				}
				print OUT "\n";
				for (my $interval = 0; $interval < scalar @{$site_bins[1]}; $interval++){
					$bins[1]->[$interval] += $site_bins[1]->[$interval];
#					print OUT $site_bins[1]->[$interval];
#					print OUT "\t";			
#					print OUT $bins[1]->[$interval];
#					print OUT "\t";
				}
#				print OUT "\n";
			}
	}
	
	print OUT "Same\n";
	for( my $interval = 0; $interval < scalar @{$bins[0]}; $interval++){
		print OUT $interval."\t".$bins[0]->[$interval]."\n";
	}
	print OUT "Different\n";
	for( my $interval = 0; $interval < scalar @{$bins[1]}; $interval++){
		print OUT $interval."\t".$bins[1]->[$interval]."\n";
	}
	my $same_median = hist_median(\@{$bins[0]});
	my $diff_median = hist_median(\@{$bins[1]});
	my $obs_difference = $diff_median-$same_median;

	
	my @bootstrap_median_diff;
	for (my $t = 0; $t < $iterate; $t++){
	    print "shuffling $t\n";
		my @shuffler_bins;
		for my $ind (@array){
			my %shuffled_distr = find_all_distances_probs($mutmap, $ind, "shuffle");
			my @shuffler_site_bins = distr_to_stathist_probs (\%shuffled_distr, $step);
				if (defined $shuffler_site_bins[0]->[1] && defined $shuffler_site_bins[1]->[1]){

				for(my $interval = 0; $interval < scalar @{$shuffler_site_bins[0]}; $interval++){
					$shuffler_bins[0]->[$interval] += $shuffler_site_bins[0]->[$interval];
				}
				for(my $interval = 0; $interval < scalar @{$shuffler_site_bins[1]}; $interval++){
					$shuffler_bins[1]->[$interval] += $shuffler_site_bins[1]->[$interval];
				}
			}
			
		}
		my $shd = hist_median(\@{$shuffler_bins[1]});
		my $shs = hist_median(\@{$shuffler_bins[0]});
		my $diff = $shd-$shs;
		print OUT "boot $shd $shs $diff \n";
		push @bootstrap_median_diff, $diff;
	}
	
	my @sorted_bootstrap = sort {$a <=> $b} @bootstrap_median_diff;
	my $pvalue = 0;
	for (my $i = 0; $i < $iterate; $i++){
		if($sorted_bootstrap[$i] >= $obs_difference){
			$pvalue = ($iterate - $i)/$iterate;
			last;
		}
	}
	print OUT "same_median\tdiff_median\tmedian_difference\tpvalue\n";
	print OUT ">\t";
	print OUT $same_median."\t"; #this have to be the median of "same" statistics
	print OUT $diff_median."\t"; #this have to be the median of "diff" statistics
	print OUT $obs_difference."\t";
	print OUT $pvalue."\n";
	close OUT;
}



# analyze sites
sub single_sites_stats {
	my $mutmap = shift;
	my $step = shift;
	my $iterate = shift;
	my $sites = shift;
	my $verbose = shift;
	
	$sites = [1..$mutmap->{static_length}] unless ($sites) ;
	
	my $tree = $mutmap->{static_tree};
	my %subs_on_node = %{$mutmap->{static_subs_on_node}};
	my %nodes_with_sub = %{$mutmap->{static_nodes_with_sub}};

	my $outfile = File::Spec->catfile($mutmap->{static_output_base}, $mutmap->{static_protein}."_".$mutmap->{static_state}."_sites_median_statistics");	
	
	open OUT, ">$outfile" or die "cannot create output file $outfile: $!";
	my @bootstrap_median_diff;

	for my $ind(@{$sites}){
		next unless ($nodes_with_sub{$ind});
		my @bootstrap_median_diff;
		print OUT ">debug site $ind \n";
		my %distr = find_all_distances_probs($mutmap, $ind);
		if ($verbose) {
			print OUT "output from find_all_distances_probs. Same [dist, prob?], Diff[dist, prob], [same count, diff count]\n";
			print OUT (Dumper(\%distr));
			print OUT "node,anc,der,prob\n";
			foreach my $node (@{$nodes_with_sub{$ind}}){
				my @sorted =  sort { $b->{"Substitution::probability"} <=> $a->{"Substitution::probability"} } @{$subs_on_node{${$node}->get_name()}->{$ind}};
				foreach my $sub(@sorted){
					print OUT ${$node}->get_name().",".$sub->{"Substitution::ancestral_allele"}.",".$sub->{"Substitution::derived_allele"}.",".$sub->{"Substitution::probability"}."\n";
				}
			}
		}
		print OUT "key\tsame_count\tdiff_count\n";
		foreach my $k (keys %distr){
			print OUT $k."\t".$distr{$k}[2][0]."\t".$distr{$k}[2][1]."\n"; # distr{$anc$der$num}[2] = (same, diff)
		}
		my @bins = distr_to_stathist_probs(\%distr, $step);
		if ($verbose) {
			my @hbins = compress_array_of_bins_to_hash(\@bins);
			print OUT "This is a compressed output from distr_to_stathist_probs (this data structure is only used for debugging purposes): {interval-> weighted same}, {interval-> weighted diff}\n";
			print OUT (Dumper(\@hbins));
		}
		my $same_median = hist_median(\@{$bins[0]});
		my $diff_median = hist_median(\@{$bins[1]});
		my $obs_difference = $diff_median-$same_median;
		print OUT ">site\tsame_median\tdiff_median\tmedian_difference\tpvalue\n";
		print OUT $ind."\t";
		print OUT $same_median."\t"; #this have to be the median of "same" statistics
		print OUT $diff_median."\t"; #this have to be the median of "diff" statistics
		print OUT $obs_difference."\t";
		
		for (my $t = 0; $t < $iterate; $t++){
			my %shuffled_distr = find_all_distances_probs($mutmap, $ind, "shuffle");
			my @shuffler_bins = distr_to_stathist_probs(\%shuffled_distr, $step);
			push @bootstrap_median_diff, hist_median(\@{$shuffler_bins[1]})-hist_median(\@{$shuffler_bins[0]});
		}
		
		my @sorted_bootstrap = sort {$a <=> $b} @bootstrap_median_diff;
		my $pvalue = 0;
		for (my $i = 0; $i < $iterate; $i++){
			if($sorted_bootstrap[$i] >= $obs_difference){
				$pvalue = ($iterate - $i)/$iterate;
				last;
			}
		}
		print OUT $pvalue;
		if ($pvalue < 0.01){
			print OUT "\tSignif";
		}
		print OUT "\n";
		}
	close OUT;
}

# for each interval prints two numbers (for each substitution, i.e. same ancestor and same derived aa or codon): observed number of convergent mutations 
# and expected from the of divergent mutations in this interval
# each subst (a, d) is analyzed separately
# counts each pair only once


## bootstrap: standard shuffler, shuffles labels on sites
sub logic_medstat_groups_labelshuffler {
	my $mutmap = shift;
	my $step = shift;
	my $iterate = shift;
	my $groupname = shift;
	my @group = @{$_[0]};
	my @array_of_sites = (1..$mutmap->{static_length});
	
	my $tree = $mutmap->{static_tree};
	my %subs_on_node = %{$mutmap->{static_subs_on_node}};
	my %nodes_with_sub = %{$mutmap->{static_nodes_with_sub}};

	my $outfile = File::Spec->catfile($mutmap->{static_output_base}, 
$mutmap->{static_protein}."_".$mutmap->{static_state}."_groups_labelshuffler_".$groupname);	
	open OUT, ">$outfile" or die "cannot create output file $outfile: $!";


	my @complement;
	
	my @bootstrap_median_diff;
	my @bins;
	my @meaningful_sites;
	
	
#my @arr = @h1_host_shift;
print OUT "site\tsame_median\tdiff_median\tmedian_difference\tpvalue\n";
	for my $ind (@array_of_sites){
		print OUT $ind."\n";
		my %distr = find_all_distances($mutmap, $ind);
		my @site_bins = distr_to_stathist(\%distr, $step);
		if (defined $site_bins[0]->[1] && defined $site_bins[1]->[1]){
			push @meaningful_sites, $ind;	
			for( my $interval = 0; $interval < scalar @{$site_bins[0]}; $interval++){
				$bins[0]->[$interval]->[$ind] = $site_bins[0]->[$interval];
			}
			for (my $interval = 0; $interval < scalar @{$site_bins[1]}; $interval++){
				$bins[1]->[$interval]->[$ind] = $site_bins[1]->[$interval];
			}
		}
	}
	
	my %group_hash;
	foreach my $gs(@group){$group_hash{$gs} = 1;}
	my @group;
	
	foreach my $ms(@meaningful_sites){
		if (exists $group_hash{$ms}){
			push @group, $ms;
		}
		else {
			push @complement, $ms;
		}
	}

	my $same_median_group = hist_median_group(\@{$bins[0]}, \@group);
	my $diff_median_group = hist_median_group(\@{$bins[1]}, \@group);
	my $obs_difference_group = $diff_median_group-$same_median_group;
	print OUT "\t size \t same median \t diff median \t difference\n";
	print OUT "group\t";
	print OUT scalar @group;
	print OUT "\t";
	print OUT $same_median_group."\t"; #this have to be the median of "same" statistics
	print OUT $diff_median_group."\t"; #this have to be the median of "diff" statistics
	print OUT $obs_difference_group."\n";
	
	my $same_median_complement = hist_median_group(\@{$bins[0]}, \@complement);
	my $diff_median_complement = hist_median_group(\@{$bins[1]}, \@complement);
	my $obs_difference_complement = $diff_median_complement-$same_median_complement;
	print OUT "complement\t";
	print OUT scalar @complement;
	print OUT "\t";
	print OUT $same_median_complement."\t"; #this have to be the median of "same" statistics
	print OUT $diff_median_complement."\t"; #this have to be the median of "diff" statistics
	print OUT $obs_difference_complement."\n";
	
	my $diffdiff = $obs_difference_group - $obs_difference_complement;
	print OUT $diffdiff."\n";
	
	my @bootstrap_median_diff;
	my @group_bootstrap;
	
	my $group_count;
	my $enrich_count;
	my $depl_count;
	
	for (my $t = 0; $t < $iterate; $t++){
		my @shuffler_bins;
		for my $ind (@array_of_sites){
			my %distr = find_all_distances($mutmap, $ind, "shuffle");
			my @site_bins = distr_to_stathist(\%distr, $step);
			if (defined $site_bins[0]->[1] && defined $site_bins[1]->[1]){
				push @meaningful_sites, $ind;	
				for( my $interval = 0; $interval < scalar @{$site_bins[0]}; $interval++){
					$shuffler_bins[0]->[$interval]->[$ind] = $site_bins[0]->[$interval];
				}
				for (my $interval = 0; $interval < scalar @{$site_bins[1]}; $interval++){
					$shuffler_bins[1]->[$interval]->[$ind] = $site_bins[1]->[$interval];
				}
			}
		}
		my $bootstrap_difference_group = hist_median_group(\@{$shuffler_bins[1]}, \@group)-hist_median_group(\@{$shuffler_bins[0]}, \@group);
		my $bootstrap_difference_complement = hist_median_group(\@{$shuffler_bins[1]}, \@complement)-hist_median_group(\@{$shuffler_bins[0]}, 
\@complement);	
		if ($bootstrap_difference_group >= $obs_difference_group){
			$group_count++;
		}
		if ($bootstrap_difference_group-$bootstrap_difference_complement >= $diffdiff){
			$enrich_count++;
		}
		if ($bootstrap_difference_group-$bootstrap_difference_complement <= $diffdiff){
			$depl_count++;
		}
	}
	
	print OUT "group pvalue ".$group_count/$iterate."\n";
	print OUT "enrichment pvalue ".$enrich_count/$iterate."\n";
	print OUT "depletion pvalue ".$depl_count/$iterate."\n";
	close OUT;
}

## bootstrap: randomly chooses group of sites
sub logic_medstat_groups {
	my $mutmap = shift;
	my $step = shift;
	my $iterate = shift;
	my $groupname = shift;
	my @group = @{$_[0]};
	my @array_of_sites = (1..$mutmap->{static_length});
	
	my $tree = $mutmap->{static_tree};
	my %subs_on_node = %{$mutmap->{static_subs_on_node}};
	my %nodes_with_sub = %{$mutmap->{static_nodes_with_sub}};

	my $outfile = File::Spec->catfile($mutmap->{static_output_base}, $mutmap->{static_protein}."_".$mutmap->{static_state}."_groups_".$groupname);	
	open OUT, ">$outfile" or die "cannot create output file $outfile: $!";
	my @complement;

	my @bootstrap_median_diff;
	my @bins;
	my %sites_hash;
	my @meaningful_sites;
	
	print OUT "site\tsame_median\tdiff_median\tmedian_difference\tpvalue\n";
	for (my $ind = 1; $ind <566; $ind++){
		print OUT $ind."\n";
		my %distr = find_all_distances($mutmap, $ind);
		my @site_bins = distr_to_stathist(\%distr, $step);
			if (defined $site_bins[0]->[1] && defined $site_bins[1]->[1]){
			push @meaningful_sites, $ind;	
			for( my $interval = 0; $interval < scalar @{$site_bins[0]}; $interval++){
				$bins[0]->[$interval]->[$ind] = $site_bins[0]->[$interval];
			}
			for (my $interval = 0; $interval < scalar @{$site_bins[1]}; $interval++){
				$bins[1]->[$interval]->[$ind] = $site_bins[1]->[$interval];
			}			
		}
	}
	
	my %group_hash;
	foreach my $gs(@group){
		$group_hash{$gs} = 1;
	}
	my @group;
	
	foreach my $ms(@meaningful_sites){
		if (exists $group_hash{$ms}){
			push @group, $ms;
		}
		else {
			push @complement, $ms;
		}
	}
	
	my $same_median_group = hist_median_group(\@{$bins[0]}, \@group);
	my $diff_median_group = hist_median_group(\@{$bins[1]}, \@group);
	my $obs_difference_group = $diff_median_group-$same_median_group;
	print OUT "\t same median \t diff median \t difference\n";
	print OUT "group\t";
	print OUT $same_median_group."\t"; #this have to be the median of "same" statistics
	print OUT $diff_median_group."\t"; #this have to be the median of "diff" statistics
	print OUT $obs_difference_group."\n";
	
	my $same_median_complement = hist_median_group(\@{$bins[0]}, \@complement);
	my $diff_median_complement = hist_median_group(\@{$bins[1]}, \@complement);
	my $obs_difference_complement = $diff_median_complement-$same_median_complement;
	print OUT "complement\t";
	print OUT $same_median_complement."\t"; #this have to be the median of "same" statistics
	print OUT $diff_median_complement."\t"; #this have to be the median of "diff" statistics
	print OUT $obs_difference_complement."\n";
	
	my $diffdiff = $obs_difference_group - $obs_difference_complement;
	print OUT $diffdiff."\n";
	
	#my @bootstrap_median_diff_group;
	#my @bootstrap_median_diff_complement;
	my $counter1 = 0;
	my $counter2 = 0;
	my $counter3 = 0;
	my $counter4 = 0;
	my $counter5;
	my $counter6;

	for (my $t = 0; $t < $iterate; $t++){
		my @bootstrap_group = shuffle @meaningful_sites;
		my @bootstrap_complement = splice (@bootstrap_group, scalar @group, scalar @meaningful_sites - scalar @group);
		
		my $same_median_group = hist_median_group(\@{$bins[0]}, \@bootstrap_group);
		my $diff_median_group = hist_median_group(\@{$bins[1]}, \@bootstrap_group);
		my $same_median_complement = hist_median_group(\@{$bins[0]}, \@bootstrap_complement);
		my $diff_median_complement = hist_median_group(\@{$bins[1]}, \@bootstrap_complement);

		if ($diff_median_group-$same_median_group - $diff_median_complement+$same_median_complement >= $diffdiff){
			$counter5++;
		}
		if ($diff_median_group-$same_median_group - $diff_median_complement+$same_median_complement <= $diffdiff){
			$counter6++;
		}
		if ($diff_median_group-$same_median_group >= $obs_difference_group){ 
			$counter1++;
			if ($diff_median_complement-$same_median_complement <= $obs_difference_complement){
				$counter2++;
			}
		}
		
		if ($diff_median_group-$same_median_group <= $obs_difference_group){ 
			$counter3++;
			if ($diff_median_complement-$same_median_complement >= $obs_difference_complement){
				$counter4++;
			}
		}
	}
	
	print OUT "pvalue e  ".$counter1/$iterate." pvalue enrichment  ".$counter2/$iterate."\n"; 
	print OUT "pvalue d ".$counter3/$iterate." pvalue depletion ".$counter4/$iterate."\n";
	print OUT "pvalue diffdiff enrichment ".$counter5/$iterate." pvalue diffdiff depletion ".$counter6/$iterate."\n";
	close OUT;
}


# merged find_all_distances_except_seq_and_sis_radius and shuffler

sub find_all_distances {
	my $mutmap = shift;
	my $site_index = shift;
	my $shuffle = shift;
	my $tree = $mutmap->{static_tree};
	my $state = $mutmap->{static_state};
	my %subs_on_node = %{$mutmap->{static_subs_on_node}};
	my %hash;
	return %hash unless ($mutmap->{static_nodes_with_sub}{$site_index});
	my @nodes = @{$mutmap->{static_nodes_with_sub}{$site_index}};
	
	my @distances_same;
	my @distances_diff;
	my $myCodonTable   = Bio::Tools::CodonTable->new();
	my %hash_of_nodes;
	#split the array of nodes by the ancestor codon 
	foreach my $node (@nodes){
		my $ancestor = ${$subs_on_node{${$node}->get_name()}}{$site_index}->{"Substitution::ancestral_allele"}; #todo to array of subst
		push @{$hash_of_nodes{$ancestor}}, $node;
	}
#print ("site $site_index\n");
#print Dumper(\%hash_of_nodes);
	foreach my $ancestor (keys %hash_of_nodes){
		my @nodes_subset = @{$hash_of_nodes{$ancestor}};
		my @shuffled;
		if ($shuffle) {@shuffled = shuffle @nodes_subset;}
		else {@shuffled = @nodes_subset;}
#print Dumper(\@nodes_subset);
		for (my $i = 0; $i < scalar @nodes_subset; $i++){
			my $sub1 = ${$subs_on_node{${$nodes_subset[$i]}->get_name()}}{$site_index};
			my $derived1;
			if ($state eq "nsyn"){
				$derived1 = $myCodonTable -> translate($sub1->{"Substitution::derived_allele"});
			}
			elsif ($state eq "syn"){
				$derived1 = ($sub1->{"Substitution::derived_allele"});
			}
			else {
				die "wrong argument in sub shuffler; only syn or nsyn accepted as the 6th arg";
			}
			my $count_same = 1; # to add the node1 itself
			my $count_diff;
#print "node 1:  $ancestor $derived1 \n";
			for (my $j = 0; $j < scalar @nodes_subset; $j++){
				if ($j == $i){ next; }
			#	if (value_is_in_array(${$shuffled[$j]}, \@{ ${$shuffled[$i]}->get_sisters })){ #!
			#		next; 
			#	}
				my $sub2 = ${$subs_on_node{${$nodes_subset[$j]}->get_name()}}{$site_index}; ##mistake found: nodes instead of nodes_subset
				my $derived2;
				if ($state eq "nsyn"){
					$derived2 = $myCodonTable -> translate($sub2->{"Substitution::derived_allele"});
				}
				elsif ($state eq "syn"){
					$derived2 = ($sub2->{"Substitution::derived_allele"});
				}
				else {
					die "wrong argument in sub shuffler; only syn or nsyn accepted as the 6th arg";
				}
#print "node 2:   $ancestor $derived2 \n";
			#	my $dist = calc_my_distance(${$shuffled[$i]}, ${$shuffled[$j]});
			    my $dist = node_distance($mutmap, ${$shuffled[$i]}, ${$shuffled[$j]}); #!
			   # my $dist =  calc_true_patristic_distance(${$shuffled[$i]}, ${$shuffled[$j]});
#print " dist $dist\n";
				if ($dist > 0){ # ie these nodes are not sequential; does not work since 02 06 2015
					if (!exists $hash{"$ancestor$derived1$i"} ){
						my @same = ();
						my @diff = ();
						$hash{"$ancestor$derived1$i"} = (\@same, \@diff);
					}
					if ($derived1 eq $derived2){
						push @{ ($hash{"$ancestor$derived1$i"})->[0] }, $dist;
						$count_same++;
					}
					else {
						push @{ ($hash{"$ancestor$derived1$i"})->[1] }, $dist;
						$count_diff++;
					}
				}
				else {
					#print "Panic! Distance between ".${$shuffled[$i]}->get_name." and ".${$shuffled[$j]}->get_name." is $dist \n";
				}

		}
		
#print " count same: $count_same, count diff: $count_diff \n";
		push @{ ($hash{"$ancestor$derived1$i"})->[2] }, $count_same;
		push @{ ($hash{"$ancestor$derived1$i"})->[2] }, $count_diff;

	}
	}

	return %hash;		
}

#$hash{"$ancestor$derived1$i"}->[0]  =  (($samedist,$pairweight),($samedist,$pairweight).. )
#$hash{"$ancestor$derived1$i"}->[1]  =  (($diffdist,$pairweight),($diffdist,$pairweight).. )
#$hash{"$ancestor$derived1$i"}->[2] = ($samecount, $diffcount)

sub find_all_distances_probs {
	my $mutmap = shift;
	my $site_index = shift;
	my $shuffle = shift;
	my $tree = $mutmap->{static_tree};
	my $state = $mutmap->{static_state};
	my %subs_on_node = %{$mutmap->{static_subs_on_node}};
	my %hash;
	return %hash unless ($mutmap->{static_nodes_with_sub}{$site_index});
	my @nodes = @{$mutmap->{static_nodes_with_sub}{$site_index}};
	
	my @distances_same;
	my @distances_diff;
	my $myCodonTable   = Bio::Tools::CodonTable->new();
	my %hash_of_nodes;
	# collect nodes with substitutions from a certain allele
#print Dumper(\@nodes);
	foreach my $node (@nodes){
		foreach my $sub(@{$subs_on_node{${$node}->get_name()}->{$site_index}}){
				my $ancestor = $sub->{"Substitution::ancestral_allele"}; 
				$hash_of_nodes{$ancestor}{${$node}->get_name()} = $node;
			}
	}
#print ("site $site_index\n");	
#print Dumper(\%hash_of_nodes);
	foreach my $ancestor (keys %hash_of_nodes){
#print $ancestor."\n";	
		my @nodes_subset =  values %{$hash_of_nodes{$ancestor}};
		my @shuffled;
		if ($shuffle) {@shuffled = shuffle @nodes_subset;}
		else {@shuffled = @nodes_subset;}
		for (my $i = 0; $i < scalar @nodes_subset; $i++){
#print ${$nodes_subset[$i]}->get_name()."\n" unless $shuffle;
#print Dumper($subs_on_node{${$nodes_subset[$i]}->get_name()}{$site_index}) unless $shuffle;
			foreach my $sub1 (@{$subs_on_node{${$nodes_subset[$i]}->get_name()}{$site_index}}){
					next if $sub1->{"Substitution::ancestral_allele"} ne $ancestor; 
					my $weight1 = $sub1->{"Substitution::probability"};
					my $derived1 = $sub1->{"Substitution::derived_allele"};
					my $count_same = $weight1; # to add the node1 itself
					my $count_diff;
					my $count_same_pairs;
					my $count_diff_pairs;
#print "node 1: ".${$nodes_subset[$i]}->get_name()." $ancestor $derived1 weight $weight1\n" unless $shuffle;
					for (my $j = 0; $j < scalar @nodes_subset; $j++){
						if ($j == $i){ next; }
					#	if (value_is_in_array(${$shuffled[$j]}, \@{ ${$shuffled[$i]}->get_sisters })){ #!
					#		next; 
					#	}
							foreach my $sub2 (@{$subs_on_node{${$nodes_subset[$j]}->get_name()}{$site_index}}){
								next if $sub2->{"Substitution::ancestral_allele"} ne $ancestor;
								my $weight2 = $sub2->{"Substitution::probability"};
								my $derived2 = ($sub2->{"Substitution::derived_allele"});

							#	my $dist = calc_my_distance(${$shuffled[$i]}, ${$shuffled[$j]});
								my $dist = node_distance($mutmap, ${$shuffled[$i]}, ${$shuffled[$j]}); #!
								
								my $pairweight = $weight1*$weight2;
								if (${$nodes_subset[$j]}->get_parent eq ${$nodes_subset[$i]}->get_parent){
									print(${$nodes_subset[$j]}->get_parent);
									$pairweight = $pairweight/$sub2->{"Substitution::ancestor_probability"};
								}
							   # my $dist =  calc_true_patristic_distance(${$shuffled[$i]}, ${$shuffled[$j]});
#print "node 2: ".${$nodes_subset[$j]}->get_name()." $ancestor $derived2 dist $dist weight2 $weight2 pairweight $pairweight \n" unless $shuffle;
								if ($pairweight > 0){ # ie these nodes are not sequential; does not work since 02 06 2015
									if (!exists $hash{"$ancestor$derived1$i"} ){
										my @same = ();
										my @diff = ();
										$hash{"$ancestor$derived1$i"} = (\@same, \@diff);
									}
									if ($derived1 eq $derived2){
										push @{ ($hash{"$ancestor$derived1$i"})->[0] }, [$dist,$pairweight];
										$count_same = $count_same+$weight2;
									}
									else {
										push @{ ($hash{"$ancestor$derived1$i"})->[1] }, [$dist,$pairweight];
										$count_diff = $count_diff+$weight2;
									}
								}
								else {
									print "Panic! Distance between ".${$shuffled[$i]}->get_name." and ".${$shuffled[$j]}->get_name." is $dist \n";
								}

							}


					}
#					 print " $ancestor$derived1$i count same: $count_same, count diff: $count_diff \n" unless $shuffle;
							push @{ ($hash{"$ancestor$derived1$i"})->[2] }, $count_same;
							push @{ ($hash{"$ancestor$derived1$i"})->[2] }, $count_diff;

		}

	}
	}

	return %hash;		
}



## careful! it does not ignore sequentials, if method is calc_true_patristic_distance instead og calc_my_distance
sub find_all_distances_except_seq_and_sis_radius_old {
	my $mutmap = shift;
	my $site_index = shift;
	my $tree = $mutmap->{static_tree};
	my $state = $mutmap->{static_state};
	my %subs_on_node = %{$mutmap->{static_subs_on_node}};
	my %hash;
	return %hash unless ($mutmap->{static_nodes_with_sub}{$site_index});
	my @nodes = @{$mutmap->{static_nodes_with_sub}{$site_index}};
	
#print "site $site_index\n";	
 
	
	my @distances_same;
	my @distances_diff;
	my $myCodonTable   = Bio::Tools::CodonTable->new();
	
	for (my $i = 0; $i < scalar @nodes; $i++){
		
		my $sub1 = ${$subs_on_node{${$nodes[$i]}->get_name()}}{$site_index};
		my $derived1;
		if ($state eq "nsyn"){
			$derived1 = $myCodonTable -> translate($sub1->{"Substitution::derived_allele"});
		}
		elsif ($state eq "syn"){
			$derived1 = ($sub1->{"Substitution::derived_allele"});
		}
		else {
			die "wrong argument in sub shuffler; only syn or nsyn accepted as the 6th arg";
		}
		my $ancestor1 = $sub1->{"Substitution::ancestral_allele"};
		
		my $count_same = 1; # to add the node1 itself
		my $count_diff;
#print "node 1:  $ancestor1 $derived1 \n";
		for (my $j = 0; $j < scalar @nodes; $j++){
			if ($j == $i){ next; }
			if (value_is_in_array(${$nodes[$j]}, \@{ ${$nodes[$i]}->get_sisters })){ 
			#	print " Hello sister!\n";
				next; 
			}
			my $sub2 = ${$subs_on_node{${$nodes[$j]}->get_name()}}{$site_index};
			my $derived2;
			if ($state eq "nsyn"){
				$derived2 = $myCodonTable -> translate($sub2->{"Substitution::derived_allele"});
			}
			elsif ($state eq "syn"){
				$derived2 = ($sub2->{"Substitution::derived_allele"});
			}
			else {
				die "wrong argument in sub shuffler; only syn or nsyn accepted as the 6th arg";
			}
			my $ancestor2 = $sub2->{"Substitution::ancestral_allele"};

			if ($ancestor1 ne $ancestor2 ){ next; }
#print "node 2:  $ancestor2 $derived2 \n";
			my $dist = calc_true_patristic_distance(${$nodes[$i]}, ${$nodes[$j]});
			if ($dist > 0){ # ie these nodes are not sequential
				if (!exists $hash{"$ancestor1$derived1$i"} ){
					my @same = ();
					my @diff = ();
					$hash{"$ancestor1$derived1$i"} = (\@same, \@diff);
				}
				if ($derived1 eq $derived2){
					push @{ ($hash{"$ancestor1$derived1$i"})->[0] }, $dist;
					$count_same++;
				}
				else {
					push @{ ($hash{"$ancestor1$derived1$i"})->[1] }, $dist;
					$count_diff++;
				}
			}
			else {
			#	print "WOA! sequential nodes here!\n";
			}

		}
#print "$ancestor1 $derived1 $i \n";
#print " count same: $count_same, count diff: $count_diff \n";
		push @{ ($hash{"$ancestor1$derived1$i"})->[2] }, $count_same;
		push @{ ($hash{"$ancestor1$derived1$i"})->[2] }, $count_diff;

	}

	return %hash;		
}




    
	sub value_is_in_array{
		my $value = $_[0];
		my @array = @{$_[1]};
		
		foreach my $v(@array){
			if ($value eq $v){
				return 1;
			}
		}
		return 0;
	}
	
	
# 2019: takes output from find_all_distances_probs, ie
#$hash{"$ancestor$derived1$i"}->[0]  =  (($samedist,$pairweight),($samedist,$pairweight).. )
#$hash{"$ancestor$derived1$i"}->[1]  =  (($diffdist,$pairweight),($diffdist,$pairweight).. )
#$hash{"$ancestor$derived1$i"}->[2] = ($samecount, $diffcount)

# as of 06/08/2019 For one site sum of $distr{$subst}->[0] = sum of $distr{$subst}->[1] = total number of non-lonely substitutions (does not depend on weights!)
sub distr_to_stathist_probs {
		my %distr = %{$_[0]};
		my $step = $_[1];
		my @bins;
		my %hash;
		my %pruned_distr;
		
		if (!%distr){return @bins;}
		# neva july normalization
		foreach my $subst(keys %distr){
#print "subst $subst \n";
			if (defined $distr{$subst}->[0] && defined $distr{$subst}->[1]) {
				my $same_size = $distr{$subst}->[2]->[0];
				my $diff_size = $distr{$subst}->[2]->[1];
				$pruned_distr{$subst} = $distr{$subst};
				my $ancestor_derived = $subst =~ s/[0-9]//gr; 
				$hash{$ancestor_derived} = 1;
			}
		}
		if (!%pruned_distr){return @bins;}
		my $mutgroups_count = scalar keys %hash;
#print "pruned distr\n";
#print (Dumper (\%pruned_distr));
		foreach my $subst(keys %pruned_distr){
#print "pruned subst $subst \n";
					my $same_size = $pruned_distr{$subst}->[2]->[0];
					my $diff_size = $pruned_distr{$subst}->[2]->[1];			
						
					my %stemp;
					my %dtemp;
					my $maxbin = 0;
					my $d_weights_sum;
					my $s_weights_sum;
					foreach my $s_distance_weight (@{$pruned_distr{$subst}->[0]}){
						my $b = bin($s_distance_weight->[0],$step);
						$stemp{$b} = $stemp{$b}+$s_distance_weight->[1];
#print "s_weights_sum was $s_weights_sum "; 							
						$s_weights_sum +=$s_distance_weight->[1];
#print "and now is	$s_weights_sum \n";				
						$maxbin = $b if ($b > $maxbin);
					}
					foreach my $d_distance_weight (@{$pruned_distr{$subst}->[1]}){
						my $b = bin($d_distance_weight->[0],$step);
						$dtemp{$b} = $dtemp{$b}+$d_distance_weight->[1];
#print "d_weights_sum was $d_weights_sum "; 								
						$d_weights_sum +=$d_distance_weight->[1];
#print "and now is	$d_weights_sum \n";	
						$maxbin = $b if ($b > $maxbin);
					}
#print "s_weights_sum $s_weights_sum d_weights_sum $d_weights_sum\n";
					foreach my $interval(0..$maxbin){
						# my $scount = $stemp{$interval}/$same_size; # normalize to get average count of mutations in i-radius of a mutation
						# my $dcount = $dtemp{$interval}/$same_size; 
						# if (!defined $scount){
							# $scount = 0; # can be 0 for mutation on one branch
						# }
						# if (!defined $dcount){
							# $dcount = 0;
						# }
					# #	print "interval $interval scount $scount same_size $same_size diff_size $diff_size mutgroups_count $mutgroups_count \n";
						# $bins[0]->[$interval] += $scount/(($same_size-1)*$mutgroups_count); # $mutgroups_count - for integral over all intervals to be 1
						# $bins[1]->[$interval] += $dcount/($diff_size*$mutgroups_count);
						
						$bins[0]->[$interval] += $stemp{$interval}/$s_weights_sum; # $mutgroups_count - for integral over all intervals to be 1
						$bins[1]->[$interval] += $dtemp{$interval}/$d_weights_sum;
					}

		}

for ( my $i  = 0; $i < scalar @{$bins[0]}; $i++){
print "interval $i s ".$bins[0]->[$i]."\n" if $bins[0]->[$i]>0;
}
for ( my $i  = 0; $i < scalar @{$bins[1]}; $i++){
print "interval $i d ".$bins[1]->[$i]."\n" if $bins[1]->[$i]>0;
}
		return @bins;
}


sub distr_to_stathist_probs_2 {
		my %distr = %{$_[0]};
		my $step = $_[1];
		my @bins;
		my %hash;
		my %pruned_distr;
		
		if (!%distr){return @bins;}
		# neva july normalization
		foreach my $subst(keys %distr){
#print "subst $subst \n";
			if (defined $distr{$subst}->[0] && defined $distr{$subst}->[1]) {
				my $same_size = $distr{$subst}->[2]->[0];
				my $diff_size = $distr{$subst}->[2]->[1];
				$pruned_distr{$subst} = $distr{$subst};
				my $ancestor_derived = $subst =~ s/[0-9]//gr; 
				$hash{$ancestor_derived} = 1;
			}
		}
		if (!%pruned_distr){return @bins;}
		my $mutgroups_count = scalar keys %hash;
#print "pruned distr\n";
#print (Dumper (\%pruned_distr));
		foreach my $subst(keys %pruned_distr){
#print "pruned subst $subst \n";
					my $same_size = $pruned_distr{$subst}->[2]->[0];
					my $diff_size = $pruned_distr{$subst}->[2]->[1];			
						
					my %stemp;
					my %dtemp;
					my $maxbin = 0;
					my $d_sum;
					my $s_sum;
					foreach my $s_distance_weight (@{$pruned_distr{$subst}->[0]}){
						my $b = bin($s_distance_weight->[0],$step);
						$stemp{$b} = $stemp{$b}+$s_distance_weight->[1];
#print "s_weights_sum was $s_weights_sum "; 							
						$s_sum +=1; # add 1 instead of weight
#print "and now is	$s_weights_sum \n";				
						$maxbin = $b if ($b > $maxbin);
					}
					foreach my $d_distance_weight (@{$pruned_distr{$subst}->[1]}){
						my $b = bin($d_distance_weight->[0],$step);
						$dtemp{$b} = $dtemp{$b}+$d_distance_weight->[1];
#print "d_weights_sum was $d_weights_sum "; 								
						$d_sum +=1; # add 1 instead of weight
#print "and now is	$d_weights_sum \n";	
						$maxbin = $b if ($b > $maxbin);
					}
#print "s_weights_sum $s_weights_sum d_weights_sum $d_weights_sum\n";
					foreach my $interval(0..$maxbin){						
						$bins[0]->[$interval] += $stemp{$interval}/$s_sum; # $mutgroups_count - for integral over all intervals to be 1
						$bins[1]->[$interval] += $dtemp{$interval}/$d_sum;
					}

		}

for ( my $i  = 0; $i < scalar @{$bins[0]}; $i++){
print "interval $i s ".$bins[0]->[$i]."\n" if $bins[0]->[$i]>0;
}
for ( my $i  = 0; $i < scalar @{$bins[1]}; $i++){
print "interval $i d ".$bins[1]->[$i]."\n" if $bins[1]->[$i]>0;
}
		return @bins;
}


sub compress_array_of_bins_to_hash {
		my @bins = @{$_[0]};
		my @shortbins;
		print scalar (@{$bins[0]})."\n";
		foreach my $interval(0..scalar (@{$bins[0]}) -1){
			if ($bins[0]->[$interval] || $bins[1]->[$interval]){
				$shortbins[0]->{$interval} = $bins[0]->[$interval];
				$shortbins[1]->{$interval} = $bins[1]->[$interval];
			}
		}
		return @shortbins;
}

sub distr_to_stathist {
		my %distr = %{$_[0]};
		my $step = $_[1];
		my @bins;
		my %hash;
		my %pruned_distr;

		# neva july normalization
		foreach my $subst(keys %distr){
		#	print "subst $subst \n";
			if (defined $distr{$subst}->[0] && defined $distr{$subst}->[1]) {
				my $same_size = $distr{$subst}->[2]->[0];
				my $diff_size = $distr{$subst}->[2]->[1];
				if ($same_size > 1 && $diff_size > 0){
					$pruned_distr{$subst} = $distr{$subst};
					my $ancestor_derived = $subst =~ s/[0-9]//gr; 
					$hash{$ancestor_derived} = 1;
				}
			}
		}
		
		my $mutgroups_count = scalar keys %hash;

		foreach my $subst(keys %pruned_distr){
#print "pruned subst $subst \n";
					my $same_size = $pruned_distr{$subst}->[2]->[0];
					my $diff_size = $pruned_distr{$subst}->[2]->[1];
				
						
						my %stemp;
						my %dtemp;
						my $maxbin = 0;
						foreach my $s_distance (@{$pruned_distr{$subst}->[0]}){
							my $b = bin($s_distance,$step);
							$stemp{$b} = $stemp{$b}+1;
							$maxbin = $b if ($b > $maxbin);
						}
						foreach my $d_distance (@{$pruned_distr{$subst}->[1]}){
							my $b = bin($d_distance,$step);
							$dtemp{$b} = $dtemp{$b}+1;
							$maxbin = $b if ($b > $maxbin);
						}
#print "maxbin $maxbin\n";
						foreach my $interval(0..$maxbin){
							my $scount = $stemp{$interval}/$same_size; # normalize to get average count of mutations in i-radius of a mutation
							my $dcount = $dtemp{$interval}/$same_size; 
							if (!defined $scount){
								$scount = 0; # can be 0 for mutation on one branch
							}
							if (!defined $dcount){
								$dcount = 0;
							}
						#	print "interval $interval scount $scount same_size $same_size diff_size $diff_size mutgroups_count $mutgroups_count \n";
							$bins[0]->[$interval] += $scount/(($same_size-1)*$mutgroups_count); # $mutgroups_count - for integral over all intervals to be 1
							$bins[1]->[$interval] += $dcount/($diff_size*$mutgroups_count);
						}

		}
		

		return @bins;
}

#takes an array of probabilities for 0,1,2...
sub hist_median{
	my @hist = @{$_[0]};
	my $summ = sum (@hist);
	return 0 if $summ == 0; # since 2018 :)
	my $head = 0;
	my $interval = 0;
	my $median = 0;
	
	while ($head < $summ/2){
		$head += $hist[$interval];
		$median = $interval;
		$interval++;
	}
	
	if ($head == $summ/2){
		$median += 0.5;
	}
#print_hist(\@hist);
	return $median;
}

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

   sub bin {
   	my $depth = $_[0];
   	my $step = $_[1];
   	
   	my $bin = int($depth/$step);
   	if ((int($depth/$step) == $depth/$step && $depth != 0) || $step == 1 || $step == 0.5){ # 0 goes to 0 bin, if step is 0.5 or 1, and to 1 bin otherwise
   		$bin -= 1;
   	}
   	return $bin+1;
   }

1;
