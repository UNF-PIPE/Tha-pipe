package findAltStart;
use strict;
use warnings;
use Bio::SeqIO;
use IO::String;
use File::Slurp;

use Exporter;
use base qw( Exporter );
our @EXPORT_OK = qw( findAltStart findGaps mkHash);

#Extends the sequences privided to the nearest ATG
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
	my @all;
	for my $ProtID (@{$gapSeqs}){
		$counter = 0;
		$workingSeq = get_gene($genomeHash, $gbkHash,$ProtID,4*$noOfGaps,0); #The gene extended with the number of gaps
##		print "first workingSeq \n\n";
##		print $workingSeq;
##		print "\n\n";
		$originalSeq = get_gene($genomeHash, $gbkHash,$ProtID,0,0); #The non-extended gene
		while ($check == 1){
			$counter++;
			print $counter . "\n\n";
			$workingSeq =~ m/(ATG\D+)/;
			$workingSeq = $1;
			$lengthDiffStart = length($workingSeq) - length($originalSeq);
			if ($lengthDiffStart == 0) { #The extended sequence is the same as the original sequence
				$extSeqs{$ProtID} = translateToProt($workingSeq);
				$check = 0;
			}
			if ($lengthDiffStart % 3 == 0) { #The new ATG is in frame with the old one
				@all = $1 =~ m/(T[AG][AG]\D+)$/g;
				$len2 = length($1);
				if ($len2 != 0) { #The stop codon is not the last one in the sequence
					$lengthDiffStartStop = length($workingSeq) - $len2; #Length difference between new start and in-between stop
					if ($lengthDiffStartStop % 3 == 0) { #The stop codon is in frame with the new start
						#Keep looking for a closer, better ATG
						$workingSeq = substr $workingSeq, -$len2; #Begin at the position of this stop codon
						#$workingSeq =
						$check = 1;
					}
					else {
						print $ProtID . "\n\n";
						for ( @all) {
							print $_ . "\n";
						} 
						print "slut på detta varv \n\n";	
					}
				}
				else { #If it is the last stop codon, return the extended sequence
					$extSeqs{$ProtID} = translateToProt($workingSeq);
					$check = 0;
				}
			}
			else {
				#log the sequence if the found ATG is not in frame with the original one 
				$logSeqs{$ProtID} = $workingSeq;
				#Keep looking for a closer, better ATG
				$workingSeq =~ m/(T[AG][AG]\D+)/;
				$len2 = length($1);
				$workingSeq = substr $workingSeq, -$len2;
				$check = 1; 
			}
		}
	}
##	print "original seq \n\n";
##	print $originalSeq;
##	print "\n\n";
##	print "preliminary extended seq \n\n";
##	print $workingSeq;
##	print "\n\n";
##	print "logSeqs \n\n";
##	print %logSeqs;
##	print "\n\n";
	print "extSeqs \n\n";
	print %extSeqs;
	print "\n\n";
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
	return $protein->seq;

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

sub mkHash {
      my %hash;                         #Predefine hash 
      my $dir = $_[0];                  #Define directory
      my @files = read_dir($dir);       #Contains all files in the directory 
      for my $file ( @files ) {         #Loop trough all files 
            my $filePath = "$dir$file"; #Current path to file
            my $hash = Bio::SeqIO->new(-file => "<$filePath", -format => 'fasta'); #Insert file into SeqIO object
            while (my $accession1 = $hash->next_seq) { # Go through the SeqIO objects
                 my $header = $accession1->id;
                 my $seq = $accession1->seq;
                 my @split = split(/\|/,$header);
                 $hash {$split[3]} = $seq; 
            }   
      }   
      return %hash;
}

1;
