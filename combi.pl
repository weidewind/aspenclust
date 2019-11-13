#!/usr/bin/perl

package Aspens;

use strict;
use List::Util qw(sum reduce);


sub number_of_variants {
	my %subs = %{$_[0]};
	
	my $places = sum(values %subs);
	my $res = 1;
	for my $aa (keys %subs){
		my $count = $subs{$aa};
		$places -= $count;
		last if $places == 0;
		my $binomial = upfact($count,$places)/fact($places);
		print $binomial." ";
		$res *= $binomial;
	}
	print "\n";
	return $res;
}

# factorial
sub fact {
	my $n = shift;
	return reduce { $a * $b } 1 .. $n;
	}
	
# (n+1)*..*(n+m)	
sub upfact {
	my $n = shift;
	my $m = shift;
	return reduce { $a * $b } $n+1 .. $n+$m;
}	

print fact(3);
print "\n";
print upfact(3,2);
print "\n";
print upfact(3,4);
print "\n";
print number_of_variants({'a' => 2, 'b' => 3});
print "\n";
print number_of_variants({'a' => 2, 'b' => 3, 'c' =>2});
print "\n";
print number_of_variants({'a' => 2});