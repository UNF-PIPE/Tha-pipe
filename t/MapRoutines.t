#!/usr/bin/perl
use strict;
use warnings;
use Test::More;

use MapRoutines qw( get_gene );
use HashRoutines  qw( mkHash make_gbk_hash );

my %genome = mkHash("t/test_data/findAltStart_data/genomes/");
my %annotation = make_gbk_hash("t/test_data/gbk_data/"); 

my $YP_988343 = "ATGCGTGAAATTATTTTTGATACTGAAACAACAGGTTTAGATAAAGATAACGACCGTATTATTGAAATCGGTTGTGTGGAAATGGTTGATCGTTATCTTACAGGACGCCAATTCCATGTCTATTTGAATCCGCAAGGGGTTATTATTCCTGACGAGGTTGTCGCAATCCACGGATTGACCAATGAGCGTTTAAAAGGTGAAAAAAAATTTGATGATATTGCTGATGAGCTTCTAGAATTTATAGATGGCGCAATGATGATTGCTCATAATGCAAATTTTGATATTAGCTTTCTTAATGCAGAATTAAAACGGGTAAATAAGCCACTGATCAGCATTGATAATGTCATCGACACATTAGCTATGGCGCGGCGAAAATTTCCTATGGGGCCTAATTCTCTTGATGTTTTATGTAAACGTTTTGGGATTGATAATAGTCACCGTATTCTTCATGGTGCTCTTCTTGATGCAGAAATTCTTGCTGATGTCTATATTGAATTAATTGGGGGTAAACAGGGAACATTAGGATTTTATAACAGTAATGGTCGTCATTTAAACACTGAAAATGGCAGCAATATTTCTTATGTCTTTAAAATTCGCCCACAGGCCCTGCCTCCAAGGTTAAGTGAACAAGAAAAAAGCATGCATGCTGACCTGGTCAACCAAATAGGCAAAAAAGCTTTATGGAATCAGTTCAAGATCCCTTAA";
my $ID = "YP_988343.1";
my $gene = get_gene(\%genome, \%annotation, $ID, 0, 0);
my $gene_extended = get_gene(\%genome, \%annotation, $ID, 10, 3);

ok( $YP_988343 eq $gene, "The gene YP_988343.1 are correctly retrieved from the genome");
ok( index ($gene_extended, $YP_988343) != -1, "The gene seq for YP_988343.1 is a substring of the extended seq");
#print $gene_extended . "\n";
#print $YP_988343 . "\n";
done_testing();
