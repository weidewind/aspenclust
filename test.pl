#!/usr/bin/perl

use strict;
use Aspens;
use ProbsMutmap;
use Groups;
use File::Spec;
use Cwd qw(abs_path cwd getcwd);
use lib getcwd(); # adds working directory to @INC
use Getopt::ArgvFile;
use Getopt::Long;
use Data::Dumper;

use Parsers;


my $args = {protein => "toy", state => "nsyn", bigtag => "test"};

# my $input_base = File::Spec->catdir(getcwd(), "data");
# my $probsfile = File::Spec->catfile($input_base, $args->{protein}.".ancestor.likelihoods");
# my %internalseqs = parse_likelihoods_as_fasta($probsfile);
# foreach my $node(keys %internalseqs){
	# print(">".$node."\n".$internalseqs{$node}."\n");
# }

my $mutmap = ProbsMutmap->new($args);

#Aspens::logic_global_median_statistics($mutmap, 0.00005, 5);
# print(Dumper($mutmap->{static_subs_on_node}{"STRAIN709_1932"}));
my @sites_subset = (2);
Aspens::logic_median_statistics($mutmap, 0.00005, 100, \@sites_subset);




# my @anc_seq;
# $anc_seq[0] = [0,0,0.8,0.2];
# $anc_seq[1] = [1,0,0,0];
# $anc_seq[2] = [0.5,0,0.3,0.2];
# sub_probabilities(\@anc_seq, "", "ACGT");


# sub sub_probabilities{
	# my $anc_seq=shift; # $anc_seq[$ind] = (0,0,0.8,0.2), $ind = 0,1,2
	# my $seq=shift;
	# my $keystring = shift;
	# my $myCodonTable = Bio::Tools::CodonTable->new();
	
	# my %anc_aas;	# $anc_aas{$aa} = $probability
	# for (my $i = 0; $i < 4; $i++){
			# for (my $j = 0; $j < 4; $j++){
				# for (my $h = 0; $h < 4; $h++){
					# my $cod_prob = $anc_seq->[0]->[$i]*$anc_seq->[1]->[$j]*$anc_seq->[2]->[$h];
					# if ($cod_prob>0){
						# my $aa = $myCodonTable->translate(substr($keystring,$i,1).substr($keystring,$j,1).substr($keystring,$h,1));
						# print("$aa $cod_prob\n");
						# $anc_aas{$aa} = $anc_aas{$aa}+$cod_prob;
					# }
					
				# }
		# }
	# }
	# print (Dumper(\%anc_aas));
	
	
	
	# #my %subs; 
# }
