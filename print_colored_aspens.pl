#!/usr/bin/perl

use File::Spec;
use Cwd qw(abs_path cwd getcwd);
use lib getcwd(); #adds working directory to @INC
use Mutmapnolim (realdata_exists, check_realdata_restriction);
use Getopt::Long;
use Getopt::ArgvFile;
use File::Path qw(make_path remove_tree);
use IPC::System::Simple qw(capture); 



my $protein = "h3";
my $state = 'nsyn';
my $input = '';
my $output = '';	# option variable with default value
my $tag = '';
my $sites = "235,189,278";
my $skip_stoppers;
my $scale;
my $circle_size=20;
my $notext;
my $passagefile;

GetOptions (	'protein=s' => \$protein,
		'state=s' => \$state,
		'skip_stoppers' => \$skip_stoppers,
		'input=s' => \$input,
		'output=s' => \$output,
		'tag=s' => \$tag,
		'sites=s' => \$sites,
		'scale=i' => \$scale,
		'circle_size=i'=> \$circle_size,
		'notext' => \$notext, #do not add text labels (derived codon) to internal nodes (used for overlaying)
		'no_print' => \$no_print, # do not print trees, only produce events file
		'passagefile=s' => \$passagefile, # filepath
	);

my @sites = split(/,/, $sites);
my $mutmap = Mutmapnolim->new({bigdatatag => $input, bigtag => $output, protein => $protein, state => $state, skip_stoppers => $skip_stoppers});
my @tree_files = $mutmap->print_colored_trees(\@sites, $tag);

unless ($no_print){
   foreach my $file(@tree_files){
	next if $file =~ m/tre$/;
	my $treefile = $mutmap->get_treefile;
	my $outfile = $file;
	my $command = "xvfb-run python drawAspenTree.py --treefile $treefile --eventfile $file --output $outfile --width 870 --circle_size $circle_size --large";
	if ($scale){ $command = $command." --scale $scale";}
	if ($notext) { $command = $command." --notext";}
	if ($passagefile) { $command = $command." --passagefile $passagefile";}
	print $command."\n";
	my $logs = capture($command);
	print $logs."\n";
   }
}
