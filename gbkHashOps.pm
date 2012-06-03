#!/usr/local/bin/perl

use strict;
use warnings;
use File::Slurp;
use Bio::SeqIO;
use Getopt::Long;

use Exporter;
use base qw( Exporter );
our @EXPORT_OK = qw(make_gbk_hash get_gene);

#Subroutine for making a hash of genbank annotation files	
sub make_gbk_hash{	
	my $genbankDir = $_[0]; #One example gbk file is downloadable from ftp://ftp.ncbi.nih.gov/genomes/Bacteria/Bartonella_bacilliformis_KC583_uid58533/
	my @files = read_dir($genbankDir);
	my $all_lines;
	my %gbkhash;

	my $ACID; #Accesion ID
	my $start; #start location for ACID;
	my $stop; #stop location for ACID;
	my $complement; #values true or false
	my $GI; # GI number
	my $geneID; #GeneID
	my $ge;

	foreach my $file (@files){
		print $file."\n";
		open GBFILE, $genbankDir.$file or die $!;		
		my @array = <GBFILE>;
		chomp(@array);
		my $toprow = shift(@array);
		
		if($toprow =~ m/^LOCUS\s*(\S*)\s/) {
			$ge = "$1";
		}	
		else{$ge = "unknown";}

		my @toprow_split = split(/\t/,$toprow);
		my $genomeID = $toprow_split[2];
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
			$gbkhash{$ACID}=[$start, $stop, $complement, $GI, $geneID, $ge];
		}
	}
	return %gbkhash
}#endbracket of subroutine make_gbk_hash

sub get_gene{
    my ($genome_path,$orthoHash,$ID,$ext_start,$ext_stop) = @_;
    #The soubroutine will retrieve a subseq (e.g. a gene) from a genome
    #Arguments: 1 = path to genome, 2 = %gbkhash, 3 = protein_id, 4 = extension of startposition backwards (optional), 5 = extension of stop position (optional)
    my $genome = Bio::SeqIO->new(-file => $genome_path, -format => "fasta"); #Makes a SeqIO-object from the genome
    my $gene;

    while (my $seq = $genome->next_seq){
        #Start and stop positions from %gbkhash
        my $start = $orthoHash->{$ID}[0] - $ext_start;
        my $stop = $orthoHash->{$ID}[1] + $ext_stop;

        #If gene is on complementary strand, take the reverse complement of the seq
        if($orthoHash->{$ID}[2] == 0){
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
1;
