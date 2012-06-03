package findAltStart;
use strict;
use warnings;
use Bio::SeqIO;
use IO::String;

use Exporter;
use base qw( Exporter );
our @EXPORT_OK = qw( findAltStart findGaps );

#Extends the sequences privided to the nearest ATG
sub findAltStart {
	my ($gapSeqs, $noOfGaps, $gbk) = @_;
	my $taxID;
	my %extSeqs;
	#testvars
	my $len1;
	my $len2;
	my $originalSeq;
	my $lengthDiffStart;
	my $lengthDiffStop;
	my $lengthDiffStartStop;
	my %logSeqs;
	my $workingSeq;
	my $check = 1;
	#end testvars
	for my $ProtID (@{$gapSeqs}){
		$workingSeq = get_gene("/home/data/NCBI-genomes/Bartonella_bacilliformis.fasta", $gbk,$ProtID,4*$noOfGaps,0);
		$originalSeq = get_gene("/home/data/NCBI-genomes/Bartonella_bacilliformis.fasta", $gbk,$ProtID,0,0);
		while ($check == 1){
			$workingSeq =~ m/(ATG\D+)/;
			$workingSeq = $1;
			$lengthDiffStart = length($workingSeq) - length($originalSeq);
			if ($lengthDiffStart == 0) { #The extended sequence is the same as the original sequence
				$extSeqs{$ProtID} = $workingSeq;
				$check = 0;
			}
			if ($lengthDiffStart % 3 == 0) { #The new ATG is in frame with the old one
				$1 =~ m/(T[AG][AG]\D+)/;
				$len2 = length($1);
				if ($len2 != 0) { #The stop codon is not the last one in the sequence
					$lengthDiffStartStop = length($workingSeq) - $len2; #Length difference between new start and in-between stop
					if ($lengthDiffStartStop % 3 == 0) { #The stop codon is in frame with the new start
						#Keep looking for a closer, better ATG
						$workingSeq = substr $workingSeq, -$len2; #Begin at the position of this stop codon
						$check = 1;
					}
				}
				else { #If it is the last stop codon, return the extended sequence
					$extSeqs{$ProtID} = $workingSeq;
					$check = 0;
				}
			}
			else {
				#log the sequence
				$logSeqs{$ProtID} = $workingSeq;
				#Keep looking for a closer, better ATG
				$workingSeq = substr $workingSeq, -$len2;
				$check = 1; 
			}
		}
	}
	return %extSeqs;
}

#Returns the ids of sequences with more than the specified number of gaps
sub findGaps {
	my $limit = $_[1];
	my $io = IO::String->new($_[0]);
	my $alignment = Bio::SeqIO->new(-fh => $io, -format => 'fasta');
	my @gap_ids;

	while(my $sequence = $alignment->next_seq()){
		my $seq_string = $sequence->seq;
		if ($seq_string =~ m/^-{$limit}/) {
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
