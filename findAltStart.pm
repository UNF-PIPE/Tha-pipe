package findAltStart;
use strict;
use warnings;
use Bio::SeqIO;
use IO::String;
use File::Slurp;

use Exporter;
use base qw( Exporter );
our @EXPORT_OK = qw( findAltStart findGaps);

#Extends the sequences provided to the nearest ATG
sub findAltStart {

	my ($gapSeqs, $noOfGaps, $gbkHash, $genomeHash) = @_;
	my $taxID;
	my %extSeqs;
	my %logSeqs;
	
	my $len1;
	my $len2;
	my $originalSeq;
	my $lengthDiffStart;
	my $lengthDiffStop;
	my $lengthDiffStartStop;
	my $workingSeq;
	my $check = 1;
	my $counter = 0;
	my $origLength;
	for my $ProtID (@{$gapSeqs}){

		$counter = 0;
		$workingSeq = get_gene($genomeHash, $gbkHash,$ProtID,4*$noOfGaps,0); #The gene 
			#extended with the specified number of gaps
		$originalSeq = get_gene($genomeHash, $gbkHash,$ProtID,0,0); #The non-extended gene
    
		unless ((my $sub = substr $workingSeq, 4*$noOfGaps) eq $originalSeq) { #unless get_gene 
				#returns two completely different sequences (parts without extension)
			print "\n\n get_gene mismatch \n\n";
#			print "o: $originalSeq \n";
#			print "w_sub: $sub \n";
			next;
		}
		$origLength = length($originalSeq);
		while ($check == 1){
			$counter++;
			$workingSeq =~ m/([AGT]TG\D+)/;
			$workingSeq = $1;
			$lengthDiffStart = length($workingSeq) - $origLength;
			if ($lengthDiffStart == 0) { #The extended sequence is the same as the
					# original sequence
				$extSeqs{$ProtID} = translateToProt($originalSeq);
        #print $extSeqs{$ProtID} . "\n";

				$check = 0;
			}
			if ($lengthDiffStart < 0) {
				print "\n Too short workingSeq: lengthDiffStart = $lengthDiffStart \n originalSeq: \n $originalSeq \n\n workingSeq: \n $workingSeq \n";
			}
			if ($lengthDiffStart % 3 == 0) { #The new ATG is in 
					#frame with the old one
				
				$1 =~ m/(T[AG][AG]\D+)$/g;
				$len2 = length($1);
				if ($len2 != 0) { #The stop codon is not the last one in the sequence
					
					$lengthDiffStartStop = length($workingSeq) - $len2; #Length difference 
						#between new start and in-between stop
					if (($lengthDiffStartStop % 3 == 0) && ($lengthDiffStartStop > 0)) { #The 
							#stop codon is in frame with the new start (and downstream of it)
						
						#Keep looking for a closer, better start codon	
						$workingSeq = substr $workingSeq, (length($workingSeq)-$len2);

						$check = 1;
					}
					else { #If no downstream stop codon if found that is in frame with
							# the new start
#						print "\n About to translate... lengthDiffStart = $lengthDiffStart \n";
						#Return the translated extended sequence
						$extSeqs{$ProtID} = translateToProt($workingSeq);

						$check = 0;	
					}
				}
				else { #If it is the last stop codon
#					print "\n About to translate... lengthDiffStart = $lengthDiffStart \n";
					#Return the translated extended sequence
					$extSeqs{$ProtID} = translateToProt($workingSeq);

					$check = 0;
				}
			}
			else { #If the found start codon is not in fram with the original one

				#log the sequence 
#				$logSeqs{$ProtID} = $workingSeq;

				#Keep looking for a closer, better start codon

				$workingSeq = substr $workingSeq, 1;

				$check = 1; 
			}
		}
	}
	return %extSeqs;
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
        my $start = $gbkHash->{$ID}[0] - $ext_start;
        my $stop = $gbkHash->{$ID}[1] + $ext_stop;

        #If gene is on complementary strand, take the reverse complement of the seq
        if($gbkHash->{$ID}[2] == 0){ 
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
