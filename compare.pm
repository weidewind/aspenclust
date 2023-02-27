package compare;

use strict;
use vars qw(@ISA @EXPORT $VERSION);
use Exporter;
$VERSION = 1.00; # Or higher
@ISA = qw(Exporter);
@EXPORT = qw(comparator get_synmuts count_substitutions nsyn_substitutions syn_substitutions nsyn_substitutions_codons is_neighbour_changing); # Symbols to autoexport (:DEFAULT tag)

use Bio::Tools::CodonTable;
use Class::Struct;
use Carp;
use Devel::StackTrace;


struct Supstitution => {
	position => '$',
	ancestral_allele => '$',
	derived_allele => '$'
};




sub comparator {
	my $class = shift;
	my $codon_table = Bio::Tools::CodonTable->new();
	my %possible_subs = (
			A => ['T', 'C', 'G'],
			T => ['A', 'C', 'G'],
			G => ['T', 'C', 'A'],
			C => ['T', 'A', 'G']
		);
	my %neighbour_hash = ();
	my $self;
	$self = { 
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


#This function counts synonimous and non synonymous substitutions between two sequences
sub count_substitutions{
	my $anc_seq=shift;
	my $seq=shift;
	my $rCDS=shift;
	my $ra_syn=shift;
	my $ra_nsyn=shift;
	@{$ra_syn}=();
	@{$ra_nsyn}=();
	return -1 unless defined($anc_seq)&&defined($seq);
	my $len=length $anc_seq;
	if(defined($rCDS)){
		$rCDS->[1]=$len-($len%3)-1 if $rCDS->[1]<0; 
	}else{
		$rCDS=[(0,$len-($len%3)-1)];
	};
	die "\nError in count_substitutions(): sequences have different lengths:\n$anc_seq\n$seq" if $len!=length $seq;
	my $myCodonTable   = Bio::Tools::CodonTable->new();
	for(my $i=$rCDS->[0];$i<=$rCDS->[1]-2;$i+=3){
		my $acod=substr $anc_seq,$i,3;
		my $cod=substr $seq,$i,3;
		next if $acod=~m/-/ || $cod=~m/-/;
		next if $acod eq $cod;
		next if $myCodonTable->is_ter_codon($acod) || $myCodonTable->is_ter_codon($cod);
		my $aaa=$myCodonTable->translate($acod);
		my $aa=$myCodonTable->translate($cod);
		$cod=substr $seq,1,$i;
		my $n=0;
		while($cod=~m/-/g){$n++};
		$n=($i-$rCDS->[0]-$n)/3+1;
		if($aaa eq $aa){
			push @{$ra_syn},$n;
		}else{
			my $p=Supstitution->new();
			$p->position($n);
			$p->ancestral_allele($aaa);
			$p->derived_allele($aa);
			push @{$ra_nsyn},$p;
		};
	};
	return @{$ra_syn}+@{$ra_nsyn}; # sum of lengths
};

sub nsyn_substitutions{
	my $anc_seq=shift;
	my $seq=shift;
	my $rCDS=shift;
	#my $ra_nsyn=shift;
	my %ra_nsyn;
	#@{$ra_nsyn}=();
	#%{$ra_nsyn} = {};
	return -1 unless defined($anc_seq)&&defined($seq);
	my $len=length $anc_seq;
	if(defined($rCDS)){
		$rCDS->[1]=$len-($len%3)-1 if $rCDS->[1]<0; 
	}else{
		$rCDS=[(0,$len-($len%3)-1)];
	};
	die "\nError in count_substitutions(): sequences have different lengths:\n$anc_seq\n$seq" if $len!=length $seq;
	my $myCodonTable   = Bio::Tools::CodonTable->new();
	for(my $i=$rCDS->[0];$i<=$rCDS->[1]-2;$i+=3){
		my $acod=substr $anc_seq,$i,3;
		my $cod=substr $seq,$i,3;
		next if $acod=~m/-/ || $cod=~m/-/;
		next if $acod eq $cod;
		next if $myCodonTable->is_ter_codon($acod) || $myCodonTable->is_ter_codon($cod);
		my $aaa=$myCodonTable->translate($acod);
		my $aa=$myCodonTable->translate($cod);
		$cod=substr $seq,1,$i;
		my $n=0;
		while($cod=~m/-/g){$n++};
		$n=($i-$rCDS->[0]-$n)/3+1;
		if($aaa ne $aa){
			my $p=Supstitution->new();
			$p->position($n);
			$p->ancestral_allele($aaa);
			$p->derived_allele($aa);
			#push @{$ra_nsyn},$p;
			print ("mutmap: ".$n."\t".$aa."\t".$aaa."\n");
			$ra_nsyn{$n} = $p;
		};
	};
	##return @{$ra_nsyn}; 
	return  %ra_nsyn;
}

sub nsyn_substitutions_codons{
	my $anc_seq=shift;
	my $seq=shift;
	my $rCDS=shift;
	#my $ra_nsyn=shift;
	my %ra_nsyn;
	#@{$ra_nsyn}=();
	#%{$ra_nsyn} = {};
	return -1 unless defined($anc_seq)&&defined($seq);
	my $len=length $anc_seq;
	if(defined($rCDS)){
		$rCDS->[1]=$len-($len%3)-1 if $rCDS->[1]<0; 
	}else{
		$rCDS=[(0,$len-($len%3)-1)];
	};
	die "\nError in count_substitutions(): sequences have different lengths:\n$anc_seq\n$seq" if $len!=length $seq;
	my $myCodonTable   = Bio::Tools::CodonTable->new();
	for(my $i=$rCDS->[0];$i<=$rCDS->[1]-2;$i+=3){
		my $acod=substr $anc_seq,$i,3;
		my $cod=substr $seq,$i,3;
		next if $acod=~m/-/ || $cod=~m/-/;
		next if (lc $acod eq lc $cod);
		next if $myCodonTable->is_ter_codon($acod) || $myCodonTable->is_ter_codon($cod);
		my $aaa=$myCodonTable->translate($acod);
		my $aa=$myCodonTable->translate($cod);
		my $codprev=substr $seq,1,$i;
		my $n=0;
		while($codprev=~m/-/g){$n++};
		$n=($i-$rCDS->[0]-$n)/3+1;
		if($aaa ne $aa){
			my $p=Supstitution->new();
			$p->position($n);
			$p->ancestral_allele($acod);
			$p->derived_allele($cod);
			#push @{$ra_nsyn},$p;
			#print ("mutmap: ".$n."\t".$aa."\t".$aaa."\n");
			$ra_nsyn{$n} = $p;
		};
	};
	##return @{$ra_nsyn}; 
	return  %ra_nsyn;
}

sub syn_substitutions{
	my $anc_seq=shift;
	my $seq=shift;
	my $rCDS=shift;
	#my $ra_nsyn=shift;
	my %ra_nsyn;
	#@{$ra_nsyn}=();
	#%{$ra_nsyn} = {};
	return -1 unless defined($anc_seq)&&defined($seq);
	my $len=length $anc_seq;
	if(defined($rCDS)){
		$rCDS->[1]=$len-($len%3)-1 if $rCDS->[1]<0; 
	}else{
		$rCDS=[(0,$len-($len%3)-1)];
	};
	die "\nError in count_substitutions(): sequences have different lengths:\n$anc_seq\n$seq" if $len!=length $seq;
	my $myCodonTable   = Bio::Tools::CodonTable->new();
	for(my $i=$rCDS->[0];$i<=$rCDS->[1]-2;$i+=3){
		my $acod=substr $anc_seq,$i,3;
		my $cod=substr $seq,$i,3;
		next if $acod=~m/-/ || $cod=~m/-/;
		next if (lc $acod eq lc $cod);
		next if $myCodonTable->is_ter_codon($acod) || $myCodonTable->is_ter_codon($cod);
		my $aaa=$myCodonTable->translate($acod);
		my $aa=$myCodonTable->translate($cod);
		my $codprev=substr $seq,1,$i;
		my $n=0;
		while($codprev=~m/-/g){$n++};
		$n=($i-$rCDS->[0]-$n)/3+1;
		if($aaa eq $aa){
			my $p=Supstitution->new();
			$p->position($n);
			$p->ancestral_allele($acod);
			$p->derived_allele($cod);
			#push @{$ra_nsyn},$p;
			#print ("mutmap: ".$n."\t".$acod."\t".$cod."\n");
			$ra_nsyn{$n} = $p;
		};
	};
	##return @{$ra_nsyn}; 
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
	if (exists $self->{neighbour_hash}{$subst->{"Supstitution::ancestral_allele"}}{$subst->{"Supstitution::derived_allele"}}{$full}){
		return $self->{neighbour_hash}{$subst->{"Supstitution::ancestral_allele"}}{$subst->{"Supstitution::derived_allele"}}{$full};
	}
	else {
	
	my $codonTable = $self->{codon_table};
	
	my %anc_neighbours = $self->get_neighbours($subst->{"Supstitution::ancestral_allele"});
	my %der_neighbours = $self->get_neighbours($subst->{"Supstitution::derived_allele"});
	
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
	
	$self->{neighbour_hash}{$subst->{"Supstitution::ancestral_allele"}}{$subst->{"Supstitution::derived_allele"}}{$full} = $answer;
	return $answer;
	}

}

sub test_is_neighbour_changing{
			my $p=Supstitution->new();
			$p->position(5);
			$p->ancestral_allele("CCT");
			$p->derived_allele("CCC");
			print is_neighbour_changing($p, 1);
			my $p=Supstitution->new();
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