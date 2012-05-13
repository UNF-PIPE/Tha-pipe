#!/usr/local/bin/perl

use strict;
use warnings;
use File::Slurp;

use Exporter;
use base qw( Exporter );
our @EXPORT_OK = qw(make_gbk_hash get_gene);

#Subroutine for making a hash of genbank annotation file	
sub make_gbk_hash{	
	my $genbankDir = $_[0]; #downloadable from ftp://ftp.ncbi.nih.gov/genomes/Bacteria/Bartonella_bacilliformis_KC583_uid58533/
	my @files = read_dir($genbankDir);
	my $all_lines;
	my %gbkhash;

	my $ACID; #Accesion ID
	my $start; #start location for ACID;
	my $stop; #stop location for ACID;
	my $complement; #values true or false
	my $GI; # GI number
	my $geneID; #GeneID

	foreach my $file (@files){
	print $file."\n";
	open GBFILE, $genbankDir.$file or die $!;		
	my @array = <GBFILE>;
	chomp(@array);
	$all_lines = join("", @array);
	#print $all_lines;
	my @array2= split(/\/translation/,$all_lines);
	foreach(@array2){
		if ($_ =~ m/CDS\s+(\d+)\.\.(\d+).\s+/){
			print "$1\t";
			$start = $1;
			$stop = $2;
			$complement = 0;	
		}
		elsif ($_ =~ m/CDS\s+complement\((\d+)\.\.(\d+)\)\s+/) {
			print "$1 $2\t";
			$start = $1;
			$stop = $2;
			$complement = 1;
		}
		if ($_ =~ m/\/protein_id=\"(.*)\"\s+\/db_xref=\"GI:(.*)\"\s+\/db_xref=\"GeneID:(.*)\"/) {
			print "$1\n";
			$ACID = $1;			
			$GI = $2;
			$geneID = $3;
		}
		#Store in hash with AccessionID as key
		$gbkhash{$ACID}=[$start, $stop, $complement, $GI, $geneID];
	}
	}
	return %gbkhash
} #endbracket of subroutine make_gbk_hash

sub get_gene{ 
    #The soubroutine will retrieve a subseq (e.g. a gene) from a genome
    #Arguments: 1 = path to genome, 2 = %gbkhash, 3 = protein_id, 4 = extension of startposition backwards (optional), 5 = extension of stop position (optional)
    my $genome = Bio::SeqIO->new(-file => $_[0], -format => "fasta"); #Makes a SeqIO-object from the genome
    my $gene;

    while (my $seq = $genome->next_seq){
	#Start and stop positions from %gbkhash
	my $start = $_[1]{$_[2]}[0] - $_[3];
	my $stop = $_[1]{$_[2]}[1] + $_[4];

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



