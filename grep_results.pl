#!/usr/bin/perl
use Getopt::Long;
use Getopt::ArgvFile;
use File::Spec;
use Cwd qw(abs_path cwd getcwd);

my $input;
my $all;
my $sites;
my $groups;
GetOptions (	
		'input=s' =>\$input, # name of a folder inside output directory, which contains "nsyn" and "syn" folders
		'all' => \$all,
		'sites' => \$sites,
		'groups' => \$groups,
	);

my @states = ("nsyn", "syn");
if ($all){
	my %allprot_pvalues;
	for my $state (@states){
		add_allprot_pvalues(\%allprot_pvalues, $input, $state);
	}
	my $dirname = dirname($input, $state);
	my $allprot_output =  File::Spec->catfile($dirname,"allprot.csv");
	open CSV, ">$allprot_output" or die "cannot create file $allprot_output - $!\n";
	print CSV "protein,nsyn,syn\n";
	foreach my $prot (sort keys %allprot_pvalues){
		print CSV $prot.",".$allprot_pvalues{$prot}{"nsyn"}.",".$allprot_pvalues{$prot}{"syn"}."\n";
	}
	close CSV;	
}
if ($sites){
	my %sites_pvalues;
	for my $state (@states){
		add_sites_pvalues(\%sites_pvalues, $input, $state);
	}
	my $dirname = dirname($input, $state);
	my $sites_output =  File::Spec->catfile($dirname,"sites.csv");
	open CSV, ">$sites_output" or die "cannot create file $sites_output - $!\n";
	print CSV "protein,state,site,pvalue\n";
	foreach my $prot (sort keys %sites_pvalues){
		foreach my $state (@states){
			foreach my $site (sort keys %{$sites_pvalues{$prot}{$state}}){
				my $pvalue = $sites_pvalues{$prot}{$state}{$site};
				#if ($pvalue < 0.05) { print CSV $prot.",".$state.",".$site.",".$pvalue."\n"; }
				print CSV $prot.",".$state.",".$site.",".$pvalue."\n"; 
			}
		}
	}
	close CSV;
}

if ($groups){
	my %groups_pvalues;
	for my $state (@states){
		add_groups_pvalues(\%groups_pvalues, $input, $state);
	}
	my $dirname = dirname($input, $state);
	my $groups_output =  File::Spec->catfile($dirname,"groups.csv");
	open CSV, ">$groups_output" or die "cannot create file $groups_output - $!\n";
	print CSV "protein,state,group,size,group_same_median,group_diff_median,compl_same_median,compl_diff_median,diffdiff,labelshuffler_group_pvalue,labelshuffler_enrichment_pvalue,labelshuffler_depletion_pvalue,siteshuffler_attraction_group_pvalue,siteshuffler_repulsion_group_pvalue,siteshuffler_enrichment_pvalue,siteshuffler_depletion_pvalue\n";
	foreach my $prot (sort keys %groups_pvalues){
		foreach my $state (@states){
			foreach my $group (sort keys %{$groups_pvalues{$prot}{$state}}){
				print "$prot $state group $group \n";
				my $lsh_pvalues_string = $groups_pvalues{$prot}{$state}{$group}{"labelshuffler"};
				my $ssh_pvalues_string = $groups_pvalues{$prot}{$state}{$group}{"siteshuffler"};
				print CSV $prot.",".$state.",".$group.",".$lsh_pvalues_string.$ssh_pvalues_string."\n"; 
			}
		}
	}
	close CSV;
}

sub add_allprot_pvalues{
	my $allprot_pvalues = shift;
	my $input = shift;
	my $state = shift;
	my $dirname = dirname($input, $state); 
	my @files = dirfiles($dirname);
	my @allprot_files = grep {/^[a-zA-Z0-9]+_${state}_global_median_statistics/} @files;
	for my $filename(@allprot_files){
		my $filepath =  File::Spec->catfile($dirname,$filename);
		open ALLPROT, "<$filepath" or die "cannot open $filepath $!";
		my $lastline;
		$lastline = $_ while <ALLPROT>;
		close ALLPROT;
		next unless $lastline =~ /^>/;
		my $prot = (split(/_/,$filename))[0];
		my $pvalue = (split(/\s++/,$lastline))[-1];
		$allprot_pvalues->{$prot}{$state} = $pvalue;
	}
}

sub add_sites_pvalues{
	my $sites_pvalues = shift;
	my $input = shift;
	my $state = shift;
	my $dirname = dirname($input, $state); 
	my @files = dirfiles($dirname);
	my @sites_files = grep {/^[a-zA-Z0-9]+_${state}_sites_median_statistics/} @files;
	for my $filename(@sites_files){
		my $filepath =  File::Spec->catfile($dirname,$filename);
		open SITES, "<$filepath" or die "cannot open $filepath $!";
		my $prot = (split(/_/,$filename))[0];
		while (<SITES>){
			my $str = <SITES> if $_ =~ /^>site/;
			my @splitter = split(/\s++/,$str);
			next if ($splitter[1].$splitter[2].$splitter[3] eq "000"); # all substs are the same or all different
			my $site = $splitter[0];
			my $pvalue = $splitter[-1];
			if ($pvalue eq "Signif"){$pvalue = $splitter[-2];}
			$sites_pvalues->{$prot}{$state}{$site} = $pvalue;
		}
		close SITES;
	}
}

sub add_groups_pvalues{
	my $groups_pvalues = shift;
	my $input = shift;
	my $state = shift;
	my $dirname = dirname($input, $state); 
	my @files = dirfiles($dirname);
	my @group_files = grep {/^[a-zA-Z0-9]+_${state}_groups_[a-zA-Z0-9_]+/} @files;
	my @labelshuffler_files = grep {/^[a-zA-Z0-9_]+_labelshuffler_[a-zA-Z0-9_]+/} @group_files;
	my @siteshuffler_files = grep {! /^[a-zA-Z0-9_]+_labelshuffler_[a-zA-Z0-9_]+/} @group_files;
	for my $filename(@labelshuffler_files){
	print $filename." labelshuffler \n";
		my $filepath =  File::Spec->catfile($dirname,$filename);
		open F, "<$filepath" or die "cannot open $filepath $!";
		my @fname = (split(/_/,$filename));
		next if ($fname[-1] eq "complement");
		my $prot = $fname[0];
		my $group = join("_",@fname[4..$#fname]);
		while (<F>){
			next unless $_ =~ /^group/;
			my @splitter = split(/\s++/,$_);
			my $pvalue_string = join(",",@splitter[1,2,3]); # group size, group same median, group diff median
			my $str = <F>;
			@splitter = split(/\s++/,$str);
			$pvalue_string = $pvalue_string.",".join(",",@splitter[2,3]); # compl same median, compl diff median
			my $i = 0;
			while ($i < 4){
				$str = <F>;
				chomp $str;
				my $val = (split(/\s++/,$str))[-1];
				$pvalue_string = $pvalue_string.",".$val; # diffdiff, 3 pvalues 
				$i++;
			}
			$groups_pvalues->{$prot}{$state}{$group}{"labelshuffler"} = $pvalue_string;
			print $pvalue_string." labelshuffler string \n";
			last;
		}
		close F;
	}
	for my $filename(@siteshuffler_files){
	print $filename." siteshuffler \n";
		my $filepath =  File::Spec->catfile($dirname,$filename);
		open F, "<$filepath" or die "cannot open $filepath $!";
		my @fname = (split(/_/,$filename));
		next if ($fname[-1] eq "complement");
		my $prot = $fname[0];
		my $group = join("_",@fname[3..$#fname]);
		while (<F>){
			next unless $_ =~ /^pvalue/;
			my $str = $_;
			my $pvalue_string;
			my $i = 0;
			while ($i < 2){
				chomp $str;
				my $val = (split(/\s++/,$str))[2];
				$pvalue_string = $pvalue_string.",".$val; # group pvalues - attraction and repulsion 
				$i++;
				$str = <F>;
			}
			chomp $str;
			@splitter = split(/\s++/,$str);
			$pvalue_string = $pvalue_string.",".join(",",@splitter[3,7]);
			$groups_pvalues->{$prot}{$state}{$group}{"siteshuffler"} = $pvalue_string;
			print $pvalue_string." siteshuffler string \n";
			last;
		}
		close F;
	}
}

sub dirfiles{
	my $dirname = shift; 
	opendir(DH, $dirname);
	my @files = readdir(DH);
	closedir(DH);
	if (scalar @files == 0){
		print "No files found in $dirname!\n";
	}
	return @files
}

sub dirname {
	my $input = shift;
	my $state = shift;
	return File::Spec->catdir(getcwd(), "output", $input, $state); 
}
