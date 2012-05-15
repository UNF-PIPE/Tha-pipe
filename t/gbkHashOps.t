#!/usr/bin/perl

use strict;
use warnings;

use Test::More;


use gbkHashOps qw(make_gbk_hash get_gene);

my $test_gbk_dir = "t/test_data/gbk_data/";
my $test_genome_path = "/home/data/NCBI-genomes/Bartonella_bacilliformis.fasta";

my %gbk_hash = make_gbk_hash( $test_gbk_dir  );

my %gene = get_gene( $test_genome_path, \%gbk_hash, 'YP_988349.1',0,0);

#print %gene;


ok( exists $gbk_hash{'YP_988349.1'}, "Test data contains group YP_988349.1");
ok( ref($gbk_hash{'YP_988349.1'}) eq 'ARRAY', 'Group YP_988349.1 is an array ref' ); 
ok( $gbk_hash{'YP_988349.1'}[0]  == 10700,   'Group YP_988349.1 has a first entry of 10700');




ok(exists $gbk_hash{'YP_190472.1'}, "Test data contains group YP_190472.1");
ok( ref($gbk_hash{'YP_190472.1'}) eq 'ARRAY', 'Group YP_190472.1 is an array ref' );
ok($gbk_hash{'YP_190472.1'}[2] == 1, 'Group YP_190472.1 has a third entry of 1 (i.e. it is complement)');




done_testing();
