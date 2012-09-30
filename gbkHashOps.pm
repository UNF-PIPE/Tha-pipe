#!/usr/local/bin/perl

use strict;
use warnings;
use File::Slurp;
use Bio::SeqIO;
use Getopt::Long;

use Exporter;
use base qw( Exporter );
our @EXPORT_OK = qw(get_gene);

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
