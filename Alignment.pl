#!/usr/bin/perl
use strict;
#use diagnostics;
use warnings;
use Bio::SeqIO;
use Bio::Perl;
#use Bio::DB::GenBank;

my %hash;
open my $IN, "<", $ARGV[0] or die "Can't open the file";
while (<$IN>) {
	my @array = split(/:\s/,$_);
	$array[0] =~ s/my_prefix//;
	my @array1 = split(/\s/,$array[1]);
	$hash {$array[0]} = [@array1];	
}

my $file = ("/home/bjorn/pipeprojekt/NCBI-proteoms/Bartonella_bacilliformis_NC_008783.faa");
my $proteome = Bio::SeqIO->new(-file => $file, -format => 'fasta');
while (my $accession = $proteome->next_seq) {
	my $acc = $accession->accession_number();
	my $seq = $accession->seq;
	print "$acc \n";
}
