#!/usr/local/bin/perl

#Script takes a feature file for a genome, and stores it in a hash;
#This script also tests to retrive info from the stored values in the hash;

use strict;
use warnings;
use Getopt::Long;

my $genbankfile= "NC_008783.gbk"; #downloadable from ftp://ftp.ncbi.nih.gov/genomes/Bacteria/Bartonella_bacilliformis_KC583_uid58533/
my %gbkhash;
my $all_lines;

my $ACID; #Accesion ID
my $start; #start location for ACID;
my $stop; #stop location for ACID;
my $complement; #values true or false
my $GI; # GI number
my $geneID; #GeneID

&make_gbk_hash;

#Test to get values from a specific Accession ID
my @retrieved_array = @{$gbkhash{"YP_989461.1"}};
print "YP_989461.1\t";
print join("\t",@retrieved_array);
print "\n";



#Subroutine for making a hash of genbank annotation file	
sub make_gbk_hash() {	
	open GBFILE, "$genbankfile" or die $!;		
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
