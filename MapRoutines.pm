package MapRoutines;
use strict;
use warnings;
use Bio::SeqIO;
use IO::String;
use File::Slurp;

use Exporter;
use base qw( Exporter );
our @EXPORT_OK = qw( get_gene translateToProt );

sub get_gene{
    my ($genomeHash,$gbkHash,$ID,$ext_start,$ext_stop) = @_; 
    #The soubroutine will retrieve a subseq (e.g. a gene) from a genome
    #Arguments: 1 = path to genome, 2 = %gbkhash, 3 = protein_id, 4 = extension of startposition backwards (optional), 5 = extension of stop position (optional)
    my $specieID = $gbkHash->{$ID}[5];
	my $genome_string = ">$specieID\n" . $genomeHash->{$specieID};
    my $io = IO::String->new($genome_string);
	my $genome = Bio::SeqIO->new(-fh => $io, -format => 'fasta');
    my $gene;

    while (my $seq = $genome->next_seq){
        #Start and stop positions from %gbkhash
        my ($start, $stop);

        #If gene is on complementary strand, take the reverse complement of the seq
        if($gbkHash->{$ID}[2] == 0){ 
            $start = $gbkHash->{$ID}[0] - $ext_start;
            $stop = $gbkHash->{$ID}[1] + $ext_stop;
            
            $gene = $seq->subseq($start, $stop);
        }   
        else{
            $start = $gbkHash->{$ID}[0] - $ext_stop;
            $stop = $gbkHash->{$ID}[1] + $ext_start;
            
            my $gene_obj = $seq->trunc($start, $stop);
            my $reversed = $gene_obj->revcom;
            $gene = $reversed->seq;
        }   
    }   
    return $gene;
}

sub translateToProt {
	my $sequence = $_[0];
	my $DNA = Bio::PrimarySeq->new ( -seq => $sequence ,
					 #-id  => 'YP_989117.1',
					 #-accession_number => 'X78121',
					 -alphabet => 'dna',
					 -is_circular => 0 );
	my $protein = $DNA->translate;
  #print substr $protein->seq,0,-1 . "\n";
	return substr $protein->seq,0,-1;
}
1;
