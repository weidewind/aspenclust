#!/usr/bin/perl

use FigTree;
use Parsers;
use File::Spec;
use Cwd qw(abs_path cwd getcwd);
use lib getcwd(); # adds working directory to @INC
use Getopt::ArgvFile;
use Getopt::Long;
use Data::Dumper;



my $treefile;
my $xparfile;
my $pair;
my $output;


GetOptions (	
		'treefile=s' => \$treefile,
		'xparfile=s' => \$xparfile,
		'pair=s' => \$pair,
		'output=s' => \$output,
);

$| = 1;

my ($site1, $site2) = split(/,/, $pair);
my $tree = parse_tree($treefile);
my @nodes1;
my @nodes2;
my %sites;
my %color;
open XPAR, "<$xparfile" or die "cannot open xpar $xparfile";
my $header = <XPAR>;
while (<XPAR>){
			my @splitter = split(/\s/);
			my $nodname = $splitter[0];
			if ($splitter[4] =~ /[A-Z]$site1[A-Z]/){
				$sites{$nodname} = $nodname."_".$site1;
				$color{$nodname} = "-16763905";
				if ($splitter[4] =~ /[A-Z]$site2[A-Z]/){
					$sites{$nodname} .= "_&_".$site2;
					$color{$nodname} = "-6750055";
				}
			}
			elsif ($splitter[4] =~ /[A-Z]$site2[A-Z]/){
				$sites{$nodname} = $nodname."_".$site2;
				$color{$nodname} = "-3407872";
			}
}
close XPAR;


open TREE, ">$output";
print TREE "#NEXUS\n\nbegin trees;\n";
print TREE "\ttree $ind = [&R] ";
my $tree_name=FigTree::tree2str($tree,sites => \%sites, color=>\%color);
print TREE $tree_name;
print TREE "\nend;\n";
my $figblock = File::Spec -> catfile(getcwd(), "figtree_block");
open BLOCK, "<$figblock" or die "Cannot open figtree_block: $!\n";
while (<BLOCK>){
	print TREE $_;
}
close BLOCK;
close TREE;