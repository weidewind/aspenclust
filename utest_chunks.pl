#!/usr/bin/perl
use Data::Dumper;

sub array_to_chunks{
	my @array = @{$_[0]};
	my $n = $_[1];
	my $total = scalar @array;
	my $chunklen = $total/$n;
	my @chunks;
	while (@array){
		my @ch = splice @array, 0, $chunklen;
		push @chunks, \@ch;
	}
	return @chunks;
}

sub new_array_to_chunks{
	my @array = @{$_[0]};
	my $n = $_[1];
	my $total = scalar @array;
	my $chunklen = $total/$n;
	my $res = $total%$n;
	my @chunks;
	while (@array){
		my $tchunklen = $chunklen;
		if ($res > 0){
			$tchunklen += 1;
			$res -= 1;
		}
		my @ch = splice @array, 0, $tchunklen;
		push @chunks, \@ch;
	}
	return @chunks;
}

#my @pairs = (0,1,2,3,4,5,6,7,8,9,10);
my $n = 10;
my @pairs = 0 .. $n-1;

my @pairchunks = new_array_to_chunks(\@pairs, 4);
print Dumper \@pairchunks;

print 2**3;

# print join(" ", @pairs[0..5]);
# my %hash;
# for (my $i = 0; $i <10; $i++){
	# $hash{$i}{$i+1} += 1;
# }
# print Dumper \%hash;