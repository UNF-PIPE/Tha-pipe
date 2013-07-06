#!/usr/bin/perl
use strict;
use warnings;

use Test::More;

use findParalogs qw( findParalogs );

#Read alignment file
my $test_path = "t/test_data/paralog_data/paralogTest.fasta";
open my $IN, '<', $test_path or die "Can't find test file: $?, $!";
my $alignment = join("", <$IN>);

#Find the paralogs
my @paralogs = findParalogs($alignment, 1);

#Perform test
my $CamelPar = 0;
my $GinkgoOrth = 1;
my $AntirPar = 0;
print "Found the following paralogs:\n";
for(@paralogs){
	print "$_\n";
	if($_ eq "Camel"){
		$CamelPar=1;
	}
        if($_ eq "Ginkgo"){
                $GinkgoOrth = 0;
        }
        if($_ eq "Antir"){
                $AntirPar = 1;
        }
}

ok($CamelPar, "The two Camel genes are paralogs");
ok($GinkgoOrth, "The two Ginkgo genes are orthologs");
ok($AntirPar, "There are paralogs among the many Antir genes");

done_testing();
