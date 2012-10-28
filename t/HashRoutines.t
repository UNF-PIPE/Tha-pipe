#!/usr/bin/perl

use strict;
use warnings;

use Test::More;

use HashRoutines qw(make_gbk_hash mkHash);

my $test_gbk_dir = "t/test_data/gbk_data/";
my $test_genome_path = "/home/data/NCBI-genomes/Bartonella_bacilliformis.fasta";

my %gbk_hash = make_gbk_hash($test_gbk_dir);

#my $tmp = $gbk_hash{'YP_190472.1'}[5];
#print "$tmp\n";
#print "gbk_hash{'YP_190472.1'}";

ok( exists $gbk_hash{'YP_988349.1'}, "Test data contains group YP_988349.1");
ok( ref($gbk_hash{'YP_988349.1'}) eq 'ARRAY', 'Group YP_988349.1 is an array ref' ); 
ok( $gbk_hash{'YP_988349.1'}[0]  == 10700,   'Group YP_988349.1 has a first entry of 10700');


ok(exists $gbk_hash{'YP_190472.1'}, "Test data contains group YP_190472.1");
ok( ref($gbk_hash{'YP_190472.1'}) eq 'ARRAY', 'Group YP_190472.1 is an array ref' );
ok($gbk_hash{'YP_190472.1'}[2] == 1, 'Group YP_190472.1 has a third entry of 1 (i.e. it is complement)');

ok(exists $gbk_hash{'YP_190472.1'}[5], "there is an element in the sixth position");
ok( $gbk_hash{'YP_190472.1'}[5] eq 'NC_006677.1', 'genome ID found and stored in hash, and can be retrieved') ;


done_testing();
