#!/usr/bin/perl
use strict;
use warnings;

use Test::More;

use findParalogs qw( findParalogs );

my $test_path = "/home/simon/UNF/fylogeniPipeline/fyl_ex/aa_Landpl_aligned_noTerminator.fasta";
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

ok($par, "Paralog array contains Camel");

done_testing();
