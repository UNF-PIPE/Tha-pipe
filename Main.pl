#!/usr/local/bin/perl

use strict;
use warnings;
use Bio::SeqIO;
use Getopt::Long;

my %gbkhash;

#Test. The "wrap around" code goes here
&make_gbk_hash("/home/data/NCBI-annotation/NC_008783.gbk");
print &get_gene("/home/data/NCBI-genomes/Bartonella_bacilliformis.fasta", \%gbkhash, $ARGV[0]);

sub get_gene{ 
	#Arguments: 1 = path to genome, 2 = %gbkhash, 3 = protein_id
	my $genome = Bio::SeqIO->new(-file => $_[0], -format => "fasta"); #Makes a SeqIO-object from the genome
	my $gene;

	while (my $seq = $genome->next_seq){
		#Start and stop positions from %gbkhash
		my $start = $_[1]{$_[2]}[0];
		my $stop = $_[1]{$_[2]}[1];

		#If gene is on complementary strand, take the reverse complement of the seq
		if($_[1]{$_[2]}[2] == 0){
			$gene = $seq->subseq($start, $stop);
		}
		else{
			my $gene_obj = $seq->trunc($start, $stop);
			my $reversed = $gene_obj->revcom;
			$gene = $reversed->seq;
		}
	}
	return $gene;
}

#Subroutine for making a hash of genbank annotation file	
sub make_gbk_hash{	
	my $genbankfile = $_[0]; #downloadable from ftp://ftp.ncbi.nih.gov/genomes/Bacteria/Bartonella_bacilliformis_KC583_uid58533/
	my $all_lines;

	my $ACID; #Accesion ID
	my $start; #start location for ACID;
	my $stop; #stop location for ACID;
	my $complement; #values true or false
	my $GI; # GI number
	my $geneID; #GeneID

	open GBFILE, $genbankfile or die $!;		
	my @array = <GBFILE>;
	chomp(@array);
	$all_lines = join("", @array);
	#print $all_lines;
	my @array2= split(/\/translation/,$all_lines);
	foreach(@array2){
		if ($_ =~ m/CDS\s+(\d+)\.\.(\d+).\s+/){
			#print "$1\t";
			$start = $1;
			$stop = $2;
			$complement = 0;	
		}
		elsif ($_ =~ m/CDS\s+complement\((\d+)\.\.(\d+)\)\s+/) {
			#print "$1 $2\t";
			$start = $1;
			$stop = $2;
			$complement = 1;
		}
		if ($_ =~ m/\/protein_id=\"(.*)\"\s+\/db_xref=\"GI:(.*)\"\s+\/db_xref=\"GeneID:(.*)\"/) {
			#print "$1\n";
			$ACID = $1;			
			$GI = $2;
			$geneID = $3;
		}
		#Store in hash with AccessionID as key
		$gbkhash{$ACID}=[$start, $stop, $complement, $GI, $geneID];
	}
} #endbracket of subroutine make_gbk_hash
