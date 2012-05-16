package findAltStart;
use strict;
use warnings;
use Bio::SeqIO;
use IO::String;

use Exporter;
use base qw( Exporter );
our @EXPORT_OK = qw( findAltStart findGaps );

sub findAltStart {
	my ($gapSeqs, $noOfGaps, $gbk) = @_;
	#my @gapSeqs = @{$_[0]};
	#my $noOfGaps = $_[1];
	#my @ext_seqs;
	my $taxID;
	my %extSeqs;

	for my $ProtID (@{$gapSeqs}){
		#$_ =~ m/^>gi\|(.*)\|(.*)/;
		#$ProtID = $1;
		#$taxID = $2;
		#print $ProtID;

		#Temporary if
		#if($_ eq 'YP_988340.1'){
		get_gene("/home/data/NCBI-genomes/Bartonella_bacilliformis.fasta", $gbk,$ProtID,4*$noOfGaps,0) =~ m/(ATG\D+)/;
		#print $1;
		$extSeqs{$ProtID} = $1
		#$ext_seqs[$count] = $1;  
		#}
	}
	return %extSeqs;

}

sub findGaps {
	my $limit = $_[1];
	my $io = IO::String->new($_[0]);
	my $alignment = Bio::SeqIO->new(-fh => $io, -format => 'fasta');
	my @gap_ids;

	while(my $sequence = $alignment->next_seq()){
		my $seq_string = $sequence->seq;
		if ($seq_string =~ m/^-{$limit}/) {
			#print ">" . $sequence->id . "\n";
			#print $seq_string . "\n\n";
			$sequence->id =~ m/^gi\|(.*)\|/;
			my $id = $1;
			push(@gap_ids, $id);
		}
	}
	return @gap_ids;
}

sub get_gene{
    my ($genome_path,$gbk,$ID,$ext_start,$ext_stop) = @_; 
    #The soubroutine will retrieve a subseq (e.g. a gene) from a genome
    #Arguments: 1 = path to genome, 2 = %gbkhash, 3 = protein_id, 4 = extension of startposition backwards (optional), 5 = extension of stop position (optional)
    my $genome = Bio::SeqIO->new(-file => $genome_path, -format => "fasta"); #Makes a SeqIO-object from the genome
    my $gene;

    while (my $seq = $genome->next_seq){
        #Start and stop positions from %gbkhash
        my $start = $gbk->{$ID}[0] - $ext_start;
        my $stop = $gbk->{$ID}[1] + $ext_stop;

        #If gene is on complementary strand, take the reverse complement of the seq
        if($gbk->{$ID}[2] == 0){ 
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
