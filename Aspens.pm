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
use List::MoreUtils qw(uniq);
use Array::Utils qw(:all);

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
	my $average = shift;
	
	my $tree = $mutmap->{static_tree};
	my %subs_on_node = %{$mutmap->{static_subs_on_node}};
	my %nodes_with_sub = %{$mutmap->{static_nodes_with_sub}};

	my $outfile = File::Spec->catfile($mutmap->{static_output_base}, $mutmap->{static_protein}."_".$mutmap->{static_state}."_global_".$average."_statistics");	
	
	open OUT, ">$outfile" or die "cannot create output file $outfile: $!";
	my @bootstrap_median_diff;
	
	my @temp = collect_distances($mutmap, "merge");
	my @bins = @{$temp[0]};

	my $same_median = hist_average(\%{$bins[0]},$average);
	my $diff_median = hist_average(\%{$bins[1]}, $average);
	my $obs_difference = $diff_median-$same_median;
	
	my @bootstrap_median_diff;
	for (my $t = 0; $t < $iterate; $t++){
	    print "shuffling $t\n";
		my @temp = collect_distances($mutmap, "merge", "shuffle",);
		my @shuffler_bins = @{$temp[0]};
		my $shd = hist_average(\%{$shuffler_bins[1]}, $average);
		my $shs = hist_average(\%{$shuffler_bins[0]}, $average);
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
	print OUT "same_$average\tdiff_$average\t".$average."_difference\tpvalue\n";
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
	my $average = shift;
	
	$sites = [1..$mutmap->{static_length}] unless ($sites) ;
	
	my $tree = $mutmap->{static_tree};
	my %subs_on_node = %{$mutmap->{static_subs_on_node}};
	my %nodes_with_sub = %{$mutmap->{static_nodes_with_sub}};

	my $outfile = File::Spec->catfile($mutmap->{static_output_base}, $mutmap->{static_protein}."_".$mutmap->{static_state}."_sites_".$average."_statistics");	
	
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
		my @bins = distr_to_stathist_probs(\%distr);
		if ($verbose) {
			print OUT "Output from distr_to_stathist_probs : {interval-> weighted same}, {interval-> weighted diff}\n";
			print OUT (Dumper(\@bins));
		}
		my $same_median = hist_average(\%{$bins[0]}, $average);
		my $diff_median = hist_average(\%{$bins[1]}, $average);
		my $obs_difference = $diff_median-$same_median;
		print OUT ">same_$average\tdiff_$average\t".$average."_difference\tpvalue\n";
		print OUT $ind."\t";
		print OUT $same_median."\t"; #this have to be the median of "same" statistics
		print OUT $diff_median."\t"; #this have to be the median of "diff" statistics
		print OUT $obs_difference."\t";
		
		for (my $t = 0; $t < $iterate; $t++){
			my %shuffled_distr = find_all_distances_probs($mutmap, $ind, "shuffle");
			my @shuffler_bins = distr_to_stathist_probs(\%shuffled_distr);
			push @bootstrap_median_diff, hist_average(\%{$shuffler_bins[1]}, $average)-hist_average(\%{$shuffler_bins[0]}, $average);
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
sub groups_stats_labelshuffler {
	my $mutmap = shift;
	my $iterate = shift;
	my $groupname = shift;
	my @group = @{$_[0]};
	my $average = $_[1]; #mean or median

	my $outfile = File::Spec->catfile($mutmap->{static_output_base}, $mutmap->{static_protein}."_".$mutmap->{static_state}."_groups_labelshuffler_".$groupname."_".$average);	
	my ($diffdiff, $obs_difference_group, $obs_difference_complement, $group, $complement, $bins) = group_stats($mutmap, \@group, $average, $outfile);

	my $group_count;
	my $enrich_count;
	my $depl_count;
	
	for (my $t = 0; $t < $iterate; $t++){
		my @temp = collect_distances($mutmap,"", "shuffle");
		my @shuffler_bins = @{$temp[0]};
		my $bootstrap_difference_group = hist_group(\%{$shuffler_bins[1]}, $group, $average)-hist_group(\%{$shuffler_bins[0]}, $group, $average);
		my $bootstrap_difference_complement = hist_group(\%{$shuffler_bins[1]}, $complement, $average)-hist_group(\%{$shuffler_bins[0]}, $complement, $average);	
		if ($bootstrap_difference_group >= $obs_difference_group){$group_count++;}
		if ($bootstrap_difference_group-$bootstrap_difference_complement >= $diffdiff){$enrich_count++;}
		if ($bootstrap_difference_group-$bootstrap_difference_complement <= $diffdiff){$depl_count++;}
	}
	
	open OUT, ">>$outfile" or die "cannot create output file $outfile: $!";
	print OUT "group pvalue\t".$group_count/$iterate."\n";
	print OUT "enrichment pvalue\t".$enrich_count/$iterate."\n";
	print OUT "depletion pvalue\t".$depl_count/$iterate."\n";
	close OUT;
}


## bootstrap: randomly chooses group of sites
sub groups_stats_siteshuffler {
	my $mutmap = shift;
	my $iterate = shift;
	my $groupname = shift;
	my @group = @{$_[0]};
	my $average = $_[1];
	
	my $outfile = File::Spec->catfile($mutmap->{static_output_base}, $mutmap->{static_protein}."_".$mutmap->{static_state}."_groups_siteshuffler_".$groupname."_".$average);	

	my ($diffdiff, $obs_difference_group, $obs_difference_complement, $group, $complement, $bins) = group_stats($mutmap, \@group, $average, $outfile);
	my @meaningful_sites = (@{$group}, @{$complement});
	my $counter1 = 0;
	my $counter2 = 0;
	my $counter3 = 0;
	my $counter4 = 0;
	my $counter5;
	my $counter6;

	for (my $t = 0; $t < $iterate; $t++){
		my @bootstrap_group = shuffle @meaningful_sites;
		my @bootstrap_complement = splice (@bootstrap_group, scalar @{$group}, scalar @meaningful_sites - scalar @{$group});
		
		my $same_median_group = hist_group(\%{$bins->[0]}, \@bootstrap_group, $average);
		my $diff_median_group = hist_group(\%{$bins->[1]}, \@bootstrap_group, $average);
		my $same_median_complement = hist_group(\%{$bins->[0]}, \@bootstrap_complement, $average);
		my $diff_median_complement = hist_group(\%{$bins->[1]}, \@bootstrap_complement, $average);

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
	
	open OUT, ">>$outfile" or die "cannot create output file $outfile: $!";	
	print OUT "pvalue e\t".$counter1/$iterate."\tpvalue enrichment\t".$counter2/$iterate."\n"; 
	print OUT "pvalue d\t".$counter3/$iterate."\tpvalue depletion\t".$counter4/$iterate."\n";
	print OUT "pvalue diffdiff enrichment\t".$counter5/$iterate."\tpvalue diffdiff depletion\t".$counter6/$iterate."\n";
	close OUT;
}

sub group_stats {
	my $mutmap = shift;
	my @group = @{$_[0]};
	my $average = $_[1];
	my $outfile = $_[2];

	my @array_of_sites = (1..$mutmap->{static_length});
	my @complement = array_diff(@array_of_sites, @group);
	my @temp = collect_distances($mutmap);
	my @bins = @{$temp[0]};
	my @meaningful_sites = @{$temp[1]};
	my @group = intersect(@group, @meaningful_sites);
	my @complement = intersect(@complement, @meaningful_sites);
	
	open OUT, ">$outfile" or die "Cannot open $outfile: $!";
	
	my $same_median_group = hist_group(\%{$bins[0]}, \@group, $average);
	my $diff_median_group = hist_group(\%{$bins[1]}, \@group, $average);
	my $obs_difference_group = $diff_median_group-$same_median_group;
	print OUT "\tsize\tsame $average\tdiff $average\tdifference\n";
	print OUT "group\t";
	print OUT scalar @group;
	print OUT "\t";
	print OUT $same_median_group."\t"; #this have to be the median of "same" statistics
	print OUT $diff_median_group."\t"; #this have to be the median of "diff" statistics
	print OUT $obs_difference_group."\n";
	
	my $same_median_complement = hist_group(\%{$bins[0]}, \@complement, $average);
	my $diff_median_complement = hist_group(\%{$bins[1]}, \@complement, $average);
	my $obs_difference_complement = $diff_median_complement-$same_median_complement;
	print OUT "complement\t";
	print OUT scalar @complement;
	print OUT "\t";
	print OUT $same_median_complement."\t"; #this have to be the median of "same" statistics
	print OUT $diff_median_complement."\t"; #this have to be the median of "diff" statistics
	print OUT $obs_difference_complement."\n";
	
	my $diffdiff = $obs_difference_group - $obs_difference_complement;
	print OUT "diffdiff\t".$diffdiff."\n";
	close OUT;
	return ($diffdiff, $obs_difference_group, $obs_difference_complement, \@group, \@complement, \@bins);
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
								my $pairweight = pairweight(${$nodes_subset[$j]},${$nodes_subset[$i]}, $sub1, $sub2);

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

# as of 06/08/2019 For one site sum of $distr{$subst}->[0] = sum of $distr{$subst}->[1] = total number of non-lonely substitutions (does not depend on 
#weights!)
sub distr_to_stathist_probs {
		my %distr = %{$_[0]};
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
		foreach my $subst(keys %pruned_distr){
					my $same_size = $pruned_distr{$subst}->[2]->[0];
					my $diff_size = $pruned_distr{$subst}->[2]->[1];			
						
					my %stemp;
					my %dtemp;
					my $d_weights_sum;
					my $s_weights_sum;
					foreach my $s_distance_weight (@{$pruned_distr{$subst}->[0]}){
						my $b = $s_distance_weight->[0];
						$stemp{$b} = $stemp{$b}+$s_distance_weight->[1];						
						$s_weights_sum +=$s_distance_weight->[1];		
					}
					foreach my $d_distance_weight (@{$pruned_distr{$subst}->[1]}){
						my $b = $d_distance_weight->[0];
						$dtemp{$b} = $dtemp{$b}+$d_distance_weight->[1];								
						$d_weights_sum +=$d_distance_weight->[1];
					}
					#print "Printing unique distances\n";
					my @alldists = uniq(keys %stemp, keys %dtemp); 
					#print Dumper (\@alldists);
					foreach my $interval(@alldists){
						$bins[0]->{$interval} += $stemp{$interval}/$s_weights_sum; # $mutgroups_count - for integral over all intervals to be 1
						$bins[1]->{$interval} += $dtemp{$interval}/$d_weights_sum;
					}

		}

# for ( my $i  = 0; $i < scalar @{$bins[0]}; $i++){
# print "interval $i s ".$bins[0]->{$i}."\n" if $bins[0]->{$i}>0;
# }
# for ( my $i  = 0; $i < scalar @{$bins[1]}; $i++){
# print "interval $i d ".$bins[1]->{$i}."\n" if $bins[1]->{$i}>0;
# }
		return @bins;
}


sub collect_distances{
	my $mutmap = shift;
	my $merge = shift;
	my $shuffle = shift;
	my @array = (1..$mutmap->{static_length});
	my @bins;
	my @meaningful_sites;
	for my $ind (@array){
		my %distr = find_all_distances_probs($mutmap, $ind, $shuffle);
		my @site_bins = distr_to_stathist_probs(\%distr); 
		if (defined $site_bins[0] && defined $site_bins[1]){
			push @meaningful_sites, $ind;
			foreach my $interval (keys %{$site_bins[0]}){
				if ($merge){
					$bins[0]->{$interval} += $site_bins[0]->{$interval};
				}
				else{
					$bins[0]->{$interval}->[$ind] = $site_bins[0]->{$interval};
				}
			}
			foreach my $interval (keys %{$site_bins[1]}){
				if ($merge){
					$bins[1]->{$interval} += $site_bins[1]->{$interval};
				}
				else{
					$bins[1]->{$interval}->[$ind] = $site_bins[1]->{$interval};
				}
			}
		}
	}
	return (\@bins, \@meaningful_sites);
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

sub hist_average {
	if ($_[1] eq "mean"){return (hist_mean($_[0]));}
	elsif ($_[1] eq "median"){return (hist_median($_[0]));}
}

#takes a hash of probabilities for distances
sub hist_median{
	my %hist = %{$_[0]};
	my $summ = sum (values %hist);
	return 0 if $summ == 0; # since 2018 :)
	my $head = 0;
	my $i = 0;
	my $median = 0;
	my @sorted_dists = sort(keys %hist);
	
	while ($head < $summ/2){
		$head += $hist{$sorted_dists[$i]};
		$median = $sorted_dists[$i];
		$i++;
	}
	
	if ($head == $summ/2){
		$median = ($sorted_dists[$i-1]+$sorted_dists[$i])/2;
	}
#print_hist(\@hist);
	return $median;
}

sub hist_mean {
	my %hist = %{$_[0]};
	my $summ = sum (values %hist);
	return 0 if $summ == 0; # since 2018 :)
	my $t;
	foreach my $dist (keys %hist){
		$t += $dist*$hist{$dist};
	}
	return $t/$summ;
}

sub hist_group {
	my %pre_hist = %{$_[0]};
	my @group = @{$_[1]};
	my $type = $_[2];
	my %hist;
	foreach my $ind(@group){
		foreach my $interval(keys %pre_hist){
			$hist{$interval} += $pre_hist{$interval}->[$ind];
		}
	}
	return (hist_average(\%hist, $type));
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
   
sub pairweight {
		my $n1 = shift;
		my $n2 = shift;
		my $sub1 = shift;
		my $sub2 = shift;
		my $w1 = $sub1->{"Substitution::probability"};
		my $w2 = $sub2->{"Substitution::probability"};
		my $pairweight = $w1*$w2;
		if ($n1->get_parent eq $n2->get_parent){
			$pairweight = $pairweight/$sub2->{"Substitution::ancestor_probability"};
		}
		return $pairweight;
}   



1;
