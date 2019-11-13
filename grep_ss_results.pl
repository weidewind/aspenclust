#!/usr/bin/perl
use strict;
use File::Spec;

my @prots = ("h3", "h1", "n1", "n2","n1pand","h1pand"); 
my $state = "syn";
my $dirname = "output/rerooted/full/".$state."/likelihoods/weightnorm";
my $stattype = "mean";
foreach my $i (@prots){
	my $outfile = File::Spec->catfile($dirname, $i."_sites_results_".$stattype);
	open OUT, ">$outfile" or die "CAnnot open $outfile $!";
	my $path  = File::Spec->catfile($dirname, $i."_".$state."_sites_".$stattype."_statistics");
	open F, "<$path" or warn "Cannot open $path $!";
	while (<F>){
		if ($_ =~ /^key/){
			my %hash;
			$_ = <F>;
			while (!($_ =~ /^>/ || $_ =~ /^Hist/) ){
				my @splitter = split('\s+');
				my $letters =$splitter[0];
				$letters =~ m/([a-zA-Z]+).*/; 
				$hash{$1} = $splitter[1];
				$_ = <F>;
			}		
			my @sorted_keys = sort (keys %hash);
			my $string = "";
			my $sep = ",";
			for (my $j = 0; $j < scalar keys %hash; $j++){
				if ($j+1 == scalar keys %hash){$sep = "";}
				elsif (substr($sorted_keys[$j],0,3) ne substr($sorted_keys[$j+1],0,3)){$sep = ";"}
				else {$sep = ",";}
				$string = $string.$sorted_keys[$j].$hash{$sorted_keys[$j]}.$sep;
			}
			while(!($_ =~ /^>/)){
				$_ = <F>
			}
			$_ = <F>;
			chomp($_);
			$_ =~ s/\s+Signif//g;
			print OUT $_."\t".$string."\n";


	
	}

	}	
	close F;
	close OUT;	
}
