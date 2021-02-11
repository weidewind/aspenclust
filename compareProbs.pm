package compareProbs;

use strict;
use vars qw(@ISA @EXPORT $VERSION);
use Exporter;
$VERSION = 1.00; # Or higher
@ISA = qw(Exporter);
@EXPORT = qw(new count_substitutions substitutions syn_substitutions nsyn_substitutions_codons is_neighbour_changing); # Symbols to autoexport (:DEFAULT tag)

use Bio::Tools::CodonTable;
use Class::Struct;
use List::Util;
use Carp;
use Devel::StackTrace;
use Data::Dumper;

struct Substitution => {
	position => '$',
	ancestral_allele => '$',
	derived_allele => '$',
	probability => '$',
	ancestor_probability => '$'
};




sub new {
	my $class = shift;
	my $keystring = shift;
	my $codon_table = Bio::Tools::CodonTable->new();
	
	my %possible_subs = (
			A => ['T', 'C', 'G'],
			T => ['A', 'C', 'G'],
			G => ['T', 'C', 'A'],
			C => ['T', 'A', 'G']
		);
	my %neighbour_hash = ();
	my $self;
	$self = { 	keystring => $keystring,
				codon_table => $codon_table,
				possible_subs => { %possible_subs },
				neighbour_hash => { %neighbour_hash }
	};
	bless $self, $class;
	return $self;
}


sub get_codon_table{
	my $self = shift;
	return $self->{codon_table};
}

sub get_possible_subs {
	my $self = shift;
	return $self->{possible_subs};
}

sub is_unambiguous {
	my @cod = ${$_[0]};
	return all { $_ == 1 || $_ == 0 } @cod;
}


sub cod_probabilities {
	my $self = shift;
	my $seq=shift; # $seq[$ind] = (0,0,0.8,0.2), $ind = 0,1,2 OR!!! $seq = "TGA"
	my $keystring = $self->{keystring};
	
	my %cods;	# $aas{$aa} = $probability
	
	if (ref($seq) ne 'ARRAY'){
		$cods{$seq} = 1;
	}
	else{
		for (my $i = 0; $i < 4; $i++){
				for (my $j = 0; $j < 4; $j++){
					for (my $h = 0; $h < 4; $h++){
						my $cod_prob = $seq->[0]->[$i]*$seq->[1]->[$j]*$seq->[2]->[$h];
						if ($cod_prob>0.0){
							my $cod = substr($keystring,$i,1).substr($keystring,$j,1).substr($keystring,$h,1);
							
							$cods{$cod} = $cod_prob;
						}
						
					}
			}
		}
	}
	# my $sum;
	# if (scalar keys %cods > 1){
		# foreach my $cod (keys %cods){ 
			# print("$cod $cods{$cod}\t");
			# $sum += $cods{$cod}; 
			# }
		# print ($sum."\n");
	# }
	return %cods;
}

# $sub_probs{$anccod}{$dercod} = $prob
sub sub_probabilities{
	my $self = shift;
	my $anc_seq = shift; # $anc_seq[$ind] = (0,0,0.8,0.2), $ind = 0,1,2 
	my $der_seq = shift; # $der_seq[$ind] = (0,0,0.8,0.2), $ind = 0,1,2 OR!!!  = char, if derived string is terminal
	my $state = shift;
	my $myCodonTable = $self->{codon_table};
	
	my %anc_probs = $self->cod_probabilities($anc_seq); # $anc_probs{$cod} = $cod_prob
	my %der_probs = $self->cod_probabilities($der_seq);
	
	my %sub_probs;
	foreach my $anccod (keys %anc_probs){
		foreach my $dercod (keys %der_probs){
			if ($anccod ne $dercod){
				my $ancaa = $myCodonTable->translate($anccod);
				my $deraa = $myCodonTable->translate($dercod);
				if ($state eq "nsyn"){
					if ($ancaa ne $deraa){
						$sub_probs{$anccod}{$deraa} = $sub_probs{$anccod}{$deraa}+$anc_probs{$anccod}*$der_probs{$dercod};
					}
				}
				else{
					if ($ancaa eq $deraa){
						$sub_probs{$anccod}{$dercod} = $anc_probs{$anccod}*$der_probs{$dercod};
					}
				}
			}
		}
	}

	return %sub_probs; 
}

sub populate_with_subs{
	my $ra_nsyn = $_[0]; # $ra_nsyn{$ind} = (Substitution,Substitution,..)
	my %sub_probs = %{$_[1]};
	my %anc_probs = %{$_[2]};
	my $ind = $_[3];
	if (scalar keys %sub_probs > 0){
		$ra_nsyn->{$ind} = [];
	}
	
	for my $anc (keys %sub_probs){
		for my $der (keys %{$sub_probs{$anc}}){
			my $p=Substitution->new();
			$p->probability($sub_probs{$anc}{$der});	
			$p->ancestor_probability($anc_probs{$anc});
			$p->position($ind);
			$p->ancestral_allele($anc);
			$p->derived_allele($der);
			push $ra_nsyn->{$ind}, $p;
		}
	}
	return $ra_nsyn;
}

## returns $ra_nsyn{$ind} = (Substitution,Substitution,..)
sub substitutions{
	my $self=shift;
	my $anc_seq=shift; # $anc_seq[$ind] = (0,0,0.8,0.2)
	my $seq=shift;
	my $state=shift;
	my $keystring=$self->{keystring}; # ACGT
	
	return -1 unless defined($anc_seq)&&defined($seq);
	
	my %ra_nsyn; # $ra_nsyn{$ind} = (Substitution,Substitution,..)
	my @anc_seq = @{$anc_seq};
	my $alen= scalar @anc_seq;
	my $len;
	if (ref($seq) ne 'ARRAY'){
		$len = length $seq ;
	}
	else {
		my @seq = @{$seq};
		$len = scalar @seq;
	}

	die "\nError in substitutions(): sequences have different lengths:\n".$alen."\n".$len."\n" if $alen!=$len;
	for(my $i=0;$i<=$len-3;$i+=3){
		my $ind = ($i/3)+1;
		my @acod=@anc_seq[$i,$i+1,$i+2];
		my $cod;
		if (ref($seq) ne 'ARRAY'){
			$cod = substr($seq,$i,3) ;
		}
		else{
			my @seq = @{$seq};
			my @cod=@seq[$i,$i+1,$i+2];
			$cod = \@cod;
		}
		my %sub_probs = $self->sub_probabilities(\@acod,$cod,$state);
		my %anc_probs = $self->cod_probabilities(\@acod); # $anc_probs{$cod} = $cod_prob
		populate_with_subs(\%ra_nsyn, \%sub_probs, \%anc_probs, $ind);
	};
	return  %ra_nsyn;
}


# tells if synonimous substitution changes the range of one-symbol neighbours
sub is_neighbour_changing {
	my $self = shift;
	my $subst = $_[0];
	my $full = $_[1];
	if (!$full){
		$full = 0;
	}
	if (exists $self->{neighbour_hash}{$subst->{"Substitution::ancestral_allele"}}{$subst->{"Substitution::derived_allele"}}{$full}){
		return $self->{neighbour_hash}{$subst->{"Substitution::ancestral_allele"}}{$subst->{"Substitution::derived_allele"}}{$full};
	}
	else {
	
	my $codonTable = $self->{codon_table};
	
	my %anc_neighbours = $self->get_neighbours($subst->{"Substitution::ancestral_allele"});
	my %der_neighbours = $self->get_neighbours($subst->{"Substitution::derived_allele"});
	
	my $answer = 0;
	if (length(keys %anc_neighbours) != length(keys %der_neighbours)){
		$answer = 1;
	}
	else {
		if ($full == 1){ # no_neighbour_changing option, prohibits quantitive change (not only qualitative)
			foreach my $k (keys %anc_neighbours){
				if (!$der_neighbours{$k} || $der_neighbours{$k} ne $anc_neighbours{$k}){
					$answer = 1;
				}
			}
		}
		else {
			foreach my $k (keys %anc_neighbours){
				if (!$der_neighbours{$k}){
					$answer = 1;
				}
			}
		}
	}
	
	$self->{neighbour_hash}{$subst->{"Substitution::ancestral_allele"}}{$subst->{"Substitution::derived_allele"}}{$full} = $answer;
	return $answer;
	}

}

sub test_is_neighbour_changing{
			my $p=Substitution->new();
			$p->position(5);
			$p->ancestral_allele("CCT");
			$p->derived_allele("CCC");
			print is_neighbour_changing($p, 1);
			my $p=Substitution->new();
			$p->position(5);
			$p->ancestral_allele("CCA");
			$p->derived_allele("CCC");
			print is_neighbour_changing($p, 1);
			$p->position(5);
			$p->ancestral_allele("CCT");
			$p->derived_allele("CCC");
			print is_neighbour_changing($p, 1);
}

# returns a hash: key - amino acid, which can be reached in one step, value - number of paths leading to that amino acid
sub get_neighbours {
	my $self = shift;
	my $codon = $_[0];
	my %neighbours;
	my $codonTable = $self->{codon_table};
	my $possible_subs = $self->{possible_subs};
	for (my $i = 0; $i < 3; $i++){
		my $str = $codon;
		foreach my $letter (@{$possible_subs->{substr($codon, $i, 1)}}){
			substr($str, $i, 1) = $letter;
			$neighbours{$codonTable->translate($str)}++;
		}
	}
	return %neighbours;
}

sub get_synmuts {
	my $self = shift;
	my $codon = $_[0];
	my $synmuts;
	my $codonTable = $self->{codon_table};
	my $possible_subs = $self->{possible_subs};
	for (my $i = 0; $i < 3; $i++){
		my $str = $codon;
		foreach my $letter (@{$possible_subs->{substr($codon, $i, 1)}}){
			substr($str, $i, 1) = $letter;
			if ($codonTable->translate($str) eq $codonTable->translate($codon)){
				my $type = muttype($letter, substr($codon, $i, 1));
			#	print "codon $codon str $str type $type\n";
				$synmuts->{$type}{$str} = 0;
			}
		}
	}
	return $synmuts;
}

sub muttype {
	my $anc = shift;
	my $der = shift;
#	print "anc $anc der $der\n";
	if ($anc eq 'A' || $anc eq 'G'){
		if ($der eq 'G' || $der eq 'A') {return "ts";}
		else {return "tv";}
	}
	elsif ($anc eq 'T' || $anc eq 'C'){
		if ($der eq 'C' || $der eq 'T') {return "ts";}
		else {return "tv";}
	}
	else { die "Undefined letters: $anc $der\n";}
}


sub test_get_synmuts {
	use Data::Dumper;
	my $synmuts = get_synmuts("CGA");
	print Dumper $synmuts;
}
#test_get_synmuts();

sub test_get_neighbours {
	my %neigh = get_neighbours("ATG");
	foreach my $n(keys %neigh){
		print $n."\t".$neigh{$n}."\n";
	}
}

#test_is_neighbour_changing();

1;