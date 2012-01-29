#!/usr/bin/perl
use strict;
#use diagnostics;
use warnings;
use Bio::SeqIO;
use Bio::DB::GenBank;

my %hash;
open my $IN, "<", $ARGV[0] or die "Can't open the file";
while (<$IN>) {
	my @array = split(/:\s/,$_);
	$array[0] =~ s/my_prefix//;
	my @array1 = split(/\s/,$array[1]);
	$hash {$array[0]} = [@array1];	
}

 
my $db_obj = Bio::DB::GenBank->new;

my $seq_obj = $db_obj->get_Seq_by_acc("YP_988971.1"); 
$seq_obj = Bio::SeqIO->new(-file => '>sequence.fasta', -format => 'fasta' ); 
system('ls -l | grep seq');
