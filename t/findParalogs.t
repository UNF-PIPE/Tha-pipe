#!/usr/bin/perl
use strict;
use warnings;

use Test::More;

use findParalogs qw( findParalogs );

my $test_path = "t/test_data/paralog_data/test1.fasta";
open my $IN, '<', $test_path or die "Can't find test file: $?, $!";

my $alignment;
while(<$IN>){
	$alignment .= $_;
}

my @paralogs = findParalogs($alignment);
my $par = 0;
print "Found the following paralogs:\n";
for(@paralogs){
	print "$_\n";
	if($_ eq "Camel"){
		$par=1;
	}
}

ok($par, "Test data contains the paralog Camel");

done_testing();
