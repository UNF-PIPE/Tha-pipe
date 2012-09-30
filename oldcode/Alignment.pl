# or die (""aj")!/usr/bin/perl
use strict;
#use diagnostics;
use warnings;
use Bio::SeqIO;
use Bio::Perl;
use File::Slurp;
#use Bio::DB::GenBank;

my %hash;
open my $IN, "<", $ARGV[0] or die "Can't open the file";
while (<$IN>) {
	my @array = split(/:\s/,$_);
	$array[0] =~ s/my_prefix//;
	my @array1 = split(/\s/,$array[1]);
	$hash {$array[0]} = [@array1];	
}

my %genomes;
my $protDir = $ARGV[1]; 
my @files = read_dir($protDir);
for my $file2 ( @files ) {
	my $file = "$protDir$file2";
	my $proteome = Bio::SeqIO->new(-file => "<$file", -format => 'fasta');

	while (my $accession1 = $proteome->next_seq) {
		my $header = $accession1->id;
		my $seq = $accession1->seq;
		my @split = split(/\|/,$header);
		$genomes {$split[3]} = $seq;
	}
}
print $genomes{"YP_192887.1"};

