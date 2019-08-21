#!/usr/bin/perl

package Aspens;

use strict;
use Bio::Phylo::IO;
use TreeUtils::Phylo::PhyloUtils qw(remove_zero_branches);
use DistanceFinder qw(get_mrcn calc_true_patristic_distance node_distance mock_distance);
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
	my $args = shift;
	my @names = qw(mutmap simnumber stattype norm verbose);
	my ($mutmap, $simnumber, $stattype, $normalization, $verbose) = map {$args->{$_}} @names;
	
	my $tree = $mutmap->{static_tree};
	my %subs_on_node = %{$mutmap->{static_subs_on_node}};
	my %nodes_with_sub = %{$mutmap->{static_nodes_with_sub}};

	my $path = $mutmap->pathFinder({norm => $normalization});
	my $outfile = File::Spec->catfile($path, $mutmap->{static_protein}."_".$mutmap->{static_state}."_global_".$stattype."_statistics");	
	
	open OUT, ">$outfile" or die "cannot create output file $outfile: $!";
	my @bootstrap_median_diff;
	
	print "collecting distances\n";
	my @temp = collect_distances({mutmap => $mutmap, merge => "merge", shuffle=>"", norm=>$normalization});
	my @bins = @{$temp[0]};
	if ($verbose) {
		close OUT;
		my $header = "Histogram\n";
		print_histogram(\@bins, $outfile, $header);
		open OUT, ">>$outfile" or die "cannot create output file $outfile: $!";
	}

	my $same_median = hist_average(\%{$bins[0]},$stattype);
	my $diff_median = hist_average(\%{$bins[1]}, $stattype);
	my $obs_difference = $diff_median-$same_median;
	
	my @bootstrap_median_diff;
	for (my $t = 0; $t < $simnumber; $t++){
	    print "shuffling $t\n";
		my @temp = collect_distances({mutmap => $mutmap, merge => "merge", shuffle=>"shuffle", norm=>$normalization});
		my @shuffler_bins = @{$temp[0]};
		my $shd = hist_average(\%{$shuffler_bins[1]}, $stattype);
		my $shs = hist_average(\%{$shuffler_bins[0]}, $stattype);
		my $diff = $shd-$shs;
		print OUT "boot $shd $shs $diff \n";
		push @bootstrap_median_diff, $diff;
	}
	
	my @sorted_bootstrap = sort {$a <=> $b} @bootstrap_median_diff;
	my $pvalue = 0;
	for (my $i = 0; $i < $simnumber; $i++){
		if($sorted_bootstrap[$i] >= $obs_difference){
			$pvalue = ($simnumber - $i)/$simnumber;
			last;
		}
	}
	print OUT "same_$stattype\tdiff_$stattype\t".$stattype."_difference\tpvalue\n";
	print OUT ">\t";
	print OUT $same_median."\t"; #this have to be the median of "same" statistics
	print OUT $diff_median."\t"; #this have to be the median of "diff" statistics
	print OUT $obs_difference."\t";
	print OUT $pvalue."\n";
	close OUT;
	return $outfile;
}


# analyze sites
sub single_sites_stats {
	my $args = shift;
	my @names = qw(mutmap simnumber sites verbose stattype norm);
	my ($mutmap, $simnumber, $sites, $verbose, $stattype, $norm) = map {$args->{$_}} @names;
	
	$sites = [1..$mutmap->{static_length}] unless ($sites) ;
	my %subs_on_node = %{$mutmap->{static_subs_on_node}};
	my %nodes_with_sub = %{$mutmap->{static_nodes_with_sub}};

	my $path = $mutmap->pathFinder({norm => $norm});
	my $outfile = File::Spec->catfile($path, $mutmap->{static_protein}."_".$mutmap->{static_state}."_sites_".$stattype."_statistics");	
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
		
		my @bins = distr_to_stathist_probs(\%distr, $norm); 
		if ($verbose) {
			close OUT;
			my $header = "Histogram for site $ind\n";
			print_histogram(\@bins, $outfile, $header);
			open OUT, ">>$outfile" or die "cannot create output file $outfile: $!";
		}
		my $same_median = hist_average(\%{$bins[0]}, $stattype);
		my $diff_median = hist_average(\%{$bins[1]}, $stattype);
		my $obs_difference = $diff_median-$same_median;
		print OUT ">same_$stattype\tdiff_$stattype\t".$stattype."_difference\tpvalue\n";
		print OUT $ind."\t";
		print OUT $same_median."\t"; #this have to be the median of "same" statistics
		print OUT $diff_median."\t"; #this have to be the median of "diff" statistics
		print OUT $obs_difference."\t";
		
		for (my $t = 0; $t < $simnumber; $t++){
			my %shuffled_distr = find_all_distances_probs($mutmap, $ind, "shuffle");
			my @shuffler_bins = distr_to_stathist_probs(\%shuffled_distr, $norm);
			push @bootstrap_median_diff, hist_average(\%{$shuffler_bins[1]}, $stattype)-hist_average(\%{$shuffler_bins[0]}, $stattype);
		}
		
		my @sorted_bootstrap = sort {$a <=> $b} @bootstrap_median_diff;
		my $pvalue = 0;
		for (my $i = 0; $i < $simnumber; $i++){
			if($sorted_bootstrap[$i] >= $obs_difference){
				$pvalue = ($simnumber - $i)/$simnumber;
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
	return $outfile;
}

# for each interval prints two numbers (for each substitution, i.e. same ancestor and same derived aa or codon): observed number of convergent mutations 
# and expected from the of divergent mutations in this interval
# each subst (a, d) is analyzed separately
# counts each pair only once


## bootstrap: standard shuffler, shuffles labels on sites
sub group_stats {
	my $args = shift;
	my @names = qw(mutmap shuffletype simnumber groupname group stattype norm outfile verbose);
	my ($mutmap, $shuffletype, $simnumber, $groupname, $group,  $stattype, $norm, $outfile, $verbose) = map {$args->{$_}} @names;
	
	my $path = $mutmap->pathFinder({norm => $norm});
	my $outfile = File::Spec->catfile($path, $mutmap->{static_protein}."_".$mutmap->{static_state}."_groups_".$shuffletype."_".$groupname."_".$stattype);	
	my ($diffdiff, $obs_difference_group, $obs_difference_complement, $group, $complement, $bins) = group_realstats({mutmap=>$mutmap, group=>$group, stattype=>$stattype, norm => $norm, verbose=> $verbose,outfile=>$outfile});

	if ($shuffletype eq "labelshuffler"){
		my $group_count;
		my $enrich_count;
		my $depl_count;
		
		for (my $t = 0; $t < $simnumber; $t++){
			my @temp = collect_distances({mutmap => $mutmap, shuffle=>"shuffle", norm=>$norm});
			my @shuffler_bins = @{$temp[0]};
			my $bootstrap_difference_group = hist_group_average(\%{$shuffler_bins[1]}, $group, $stattype)-hist_group_average(\%{$shuffler_bins[0]}, $group, $stattype);
			my $bootstrap_difference_complement = hist_group_average(\%{$shuffler_bins[1]}, $complement, $stattype)-hist_group_average(\%{$shuffler_bins[0]}, $complement, $stattype);	
			if ($bootstrap_difference_group >= $obs_difference_group){$group_count++;}
			if ($bootstrap_difference_group-$bootstrap_difference_complement >= $diffdiff){$enrich_count++;}
			if ($bootstrap_difference_group-$bootstrap_difference_complement <= $diffdiff){$depl_count++;}
		}
		
		open OUT, ">>$outfile" or die "cannot create output file $outfile: $!";
		print OUT "group pvalue\t".$group_count/$simnumber."\n";
		print OUT "enrichment pvalue\t".$enrich_count/$simnumber."\n";
		print OUT "depletion pvalue\t".$depl_count/$simnumber."\n";
	}
	elsif ($shuffletype eq "siteshuffler"){
		my @meaningful_sites = (@{$group}, @{$complement});
		my $counter1;
		my $counter2;
		my $counter3;
		my $counter4;
		my $counter5;
		my $counter6;

		for (my $t = 0; $t < $simnumber; $t++){
			my @bootstrap_group = shuffle @meaningful_sites;
			my @bootstrap_complement = splice (@bootstrap_group, scalar @{$group}, scalar @meaningful_sites - scalar @{$group});
			
			my $same_median_group = hist_group_average(\%{$bins->[0]}, \@bootstrap_group, $stattype);
			my $diff_median_group = hist_group_average(\%{$bins->[1]}, \@bootstrap_group, $stattype);
			my $same_median_complement = hist_group_average(\%{$bins->[0]}, \@bootstrap_complement, $stattype);
			my $diff_median_complement = hist_group_average(\%{$bins->[1]}, \@bootstrap_complement, $stattype);

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
		print OUT "pvalue e\t".$counter1/$simnumber."\tpvalue enrichment\t".$counter2/$simnumber."\n"; 
		print OUT "pvalue d\t".$counter3/$simnumber."\tpvalue depletion\t".$counter4/$simnumber."\n";
		print OUT "pvalue diffdiff enrichment\t".$counter5/$simnumber."\tpvalue diffdiff depletion\t".$counter6/$simnumber."\n";
	}
	else {die "Unknown group shuffler type $shuffletype";}
	close OUT;
	return $outfile;
}



sub group_realstats {
	my $args = shift;
	my @names = qw(mutmap group verbose stattype norm outfile);
	my ($mutmap, $group, $verbose, $stattype, $norm, $outfile) = map {$args->{$_}} @names;

	my @array_of_sites = (1..$mutmap->{static_length});
	my @complement = array_diff(@array_of_sites, @{$group});
	my @temp = collect_distances({mutmap => $mutmap, norm => $norm});
	my @bins = @{$temp[0]};
	my @meaningful_sites = @{$temp[1]};
	my @group = intersect(@{$group}, @meaningful_sites);
	my @complement = intersect(@complement, @meaningful_sites);
	
	my @groupbins = splicebins(\@bins, \@group);
	my @complbins = splicebins(\@bins, \@complement);
	
	open OUT, ">$outfile" or die "Cannot open $outfile: $!";
	if ($verbose) {
		close OUT;
		my $header = "Histogram for group\n";
		print_histogram(\@groupbins, $outfile, $header);
		$header = "Histogram for complement\n";
		print_histogram(\@complbins, $outfile, $header);
		open OUT, ">>$outfile" or die "cannot create output file $outfile: $!";
	}
	#my $same_median_group = hist_group_average(\%{$bins[0]}, \@group, $stattype);
	#my $diff_median_group = hist_group_average(\%{$bins[1]}, \@group, $stattype);
	my $same_median_group = hist_average(\%{$groupbins[0]},  $stattype);
	my $diff_median_group = hist_average(\%{$groupbins[1]},  $stattype);
	my $obs_difference_group = $diff_median_group-$same_median_group;
	print OUT "\tsize\tsame $stattype\tdiff $stattype\tdifference\n";
	print OUT "group\t";
	print OUT scalar @group;
	print OUT "\t";
	print OUT $same_median_group."\t"; #this have to be the median of "same" statistics
	print OUT $diff_median_group."\t"; #this have to be the median of "diff" statistics
	print OUT $obs_difference_group."\n";
		# my $same_median_complement = hist_group_average(\%{$bins[0]}, \@complement, $stattype);
	# my $diff_median_complement = hist_group_average(\%{$bins[1]}, \@complement, $stattype);
	my $same_median_complement = hist_average(\%{$complbins[0]}, $stattype);
	my $diff_median_complement = hist_average(\%{$complbins[1]}, $stattype);
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
	if (exists $mutmap->{static_all_distances_probs}{$site_index}){
		return %{$mutmap->{static_all_distances_probs}{$site_index}} unless $shuffle;
	}
	
	my @nodes = @{$mutmap->{static_nodes_with_sub}{$site_index}};
	my @distances_same;
	my @distances_diff;
	my $myCodonTable   = Bio::Tools::CodonTable->new();
	my %hash_of_nodes;
	# collect nodes with substitutions from a certain allele
	#print Dumper(\@nodes);

	## cashed in ref2 branch with almost no effect on running time
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
							foreach my $sub2 (@{$subs_on_node{${$nodes_subset[$j]}->get_name()}{$site_index}}){
								next if $sub2->{"Substitution::ancestral_allele"} ne $ancestor;
								my $weight2 = $sub2->{"Substitution::probability"};
								my $derived2 = ($sub2->{"Substitution::derived_allele"});
								my $dist = node_distance($mutmap, ${$shuffled[$i]}, ${$shuffled[$j]}); #!
								my $pairweight = pairweight(${$nodes_subset[$j]},${$nodes_subset[$i]}, $sub1, $sub2);

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
							#print " $ancestor$derived1$i count same: $count_same, count diff: $count_diff \n" unless $shuffle;
							push @{ ($hash{"$ancestor$derived1$i"})->[2] }, $count_same;
							push @{ ($hash{"$ancestor$derived1$i"})->[2] }, $count_diff;

			}

		}
	}
	$mutmap->set_all_distances_probs(\%hash, $site_index) unless $shuffle;
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
		my $normalization = $_[1];
		unless ($normalization) {$normalization = "weightnorm";}
		my @bins;
		my %hash;
		my %pruned_distr;
		
		if (!%distr){return @bins;}
		# neva july normalization
		foreach my $subst(keys %distr){
			if (defined $distr{$subst}->[0] && defined $distr{$subst}->[1]) {
				$pruned_distr{$subst} = $distr{$subst};
				my $ancestor_derived = $subst =~ s/[0-9]//gr; 
				$hash{$ancestor_derived} = 1;
			}
		}
		if (!%pruned_distr){return @bins;}
		my $mutgroups_count = scalar keys %hash;
		foreach my $subst(keys %pruned_distr){
					my %stemp;
					my %dtemp;
					my $d_weights_sum;
					my $s_weights_sum;
					foreach my $s_dist_and_weight (@{$pruned_distr{$subst}->[0]}){
						my $d = $s_dist_and_weight->[0];
						$stemp{$d} = $stemp{$d}+$s_dist_and_weight->[1];
						switch ($normalization){
							case "weightnorm" {$s_weights_sum +=$s_dist_and_weight->[1]}
							case "countnorm" {$s_weights_sum +=1}
							else {die "Normalization option $normalization unknown"}
						}	
					}
					foreach my $d_dist_and_weight (@{$pruned_distr{$subst}->[1]}){
						my $d = $d_dist_and_weight->[0];
						$dtemp{$d} = $dtemp{$d}+$d_dist_and_weight->[1];								
						switch ($normalization){
							case "weightnorm" {$d_weights_sum +=$d_dist_and_weight->[1]}
							case "countnorm" {$d_weights_sum +=1}
							else {die "Normalization option $normalization unknown"}
						}	
					}
					#print "Printing unique distances\n";
					my @alldists = uniq(keys %stemp, keys %dtemp); 
					#print Dumper (\@alldists);
					foreach my $interval(@alldists){
						$bins[0]->{$interval} += $stemp{$interval}/$s_weights_sum; # $mutgroups_count - for integral over all intervals to be 1
						$bins[1]->{$interval} += $dtemp{$interval}/$d_weights_sum;
					}

		}

		return @bins;
}

sub collect_distances{
	my $args = shift;
	my @names = qw(mutmap merge shuffle norm);
	my ($mutmap, $merge, $shuffle, $normalization) = map {$args->{$_}} @names;
	my @array = (1..$mutmap->{static_length});
	my @bins;
	my @meaningful_sites;
	for my $ind (@array){
		my %distr = find_all_distances_probs($mutmap, $ind, $shuffle);
		my @site_bins = distr_to_stathist_probs(\%distr, $normalization); 
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



sub splicebins {
	my @bins = @{$_[0]};
	my @group = @{$_[1]};
	my @groupbins;
	
	foreach my $interval (keys %{$bins[0]}){
		foreach my $ind (@group){
			$groupbins[0]->{$interval} += $bins[0]->{$interval}->[$ind];
			$groupbins[1]->{$interval} += $bins[1]->{$interval}->[$ind];
		}
	}
	return @groupbins;
}

sub print_histogram {
			my @bins = @{$_[0]};
			my $outfile = $_[1];
			my $header = $_[2];
			open OUT, ">>$outfile" or die "CAnnot open $outfile $!";
			print OUT $header;
			print OUT "Distance,Same,Diff\n";
			my @sortedsame = sort (keys %{$bins[0]});
			foreach my $dist (@sortedsame){
				print OUT $dist.",".$bins[0]{$dist}.",".$bins[1]{$dist}."\n"
			}
			my $sumsame = sum (values %{$bins[0]});
			my $sumdiff = sum (values %{$bins[1]});
			print OUT "Sum same: $sumsame Sum diff: $sumdiff\n";
			close OUT;
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

sub hist_group_average {
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
