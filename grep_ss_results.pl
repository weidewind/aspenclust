#!/usr/bin/perl
use strict;
use File::Spec;

my @prots = ("h3","h1","n1","n2"); 
my $dirname = "output/parenttest/nsyn/likelihoods/weightnorm";
foreach my $i (@prots){
	my $outfile = File::Spec->catfile($dirname, $i."_sites_results_mean");
	open OUT, ">$outfile" or die "CAnnot open $outfile $!";
	my $path  = File::Spec->catfile($dirname, $i."_nsyn_sites_mean_statistics");
	open F, "<$path" or die "Cannot open $path $!";
	while (<F>){
		if ($_ =~ /^key/){
			my %hash;
			$_ = <F>;
			while (!($_ =~ /^>/) ){
				my @splitter = split('\s+');
				$hash{substr($splitter[0],0,4)} = $splitter[1];
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
			$_ = <F>;
			chomp($_);
			$_ =~ s/\s+Signif//g;
			print OUT $_."\t".$string."\n";


	
	}

	}	
	close F;
	close OUT;	
}
