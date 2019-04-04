#!/usr/bin/perl

package ProbsMutmap;

use strict;
use Bio::Phylo::IO;
use Class::Struct;
use compare qw(nsyn_substitutions syn_substitutions nsyn_substitutions_codons is_neighbour_changing);
use Bio::SeqIO;
use Cwd qw(abs_path cwd getcwd);
use File::Path qw(make_path remove_tree);
use lib getcwd(); # adds working directory to @INC
use Data::Dumper;
use IPC::System::Simple qw(capture);



sub probsmutmap {
	
}