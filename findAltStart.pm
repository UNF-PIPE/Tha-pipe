package findAltStart;
use strict;
use warnings;
use Bio::SeqIO;
use IO::String;
use File::Slurp;
use MapRoutines qw(get_gene translateToProt);
use Exporter;
use base qw( Exporter );
our @EXPORT_OK = qw( findAltStart findGaps);

#Extends the sequences provided to the nearest ATG
sub findAltStart {
	my ($gapSeqs, $noOfGaps, $gbkHash, $genomeHash) = @_;
	my $taxID;
	my %extSeqs;
	my %logSeqs;
	
	my $originalSeq;
	my $lengthDiffStart;
	my $lengthDiffStop;
	my $lengthDiffStartStop;
	my $workingSeq;
	my $check = 1;
	my $origLength;
	my $candidateSeq;
	my $checkForStop;
	my @arr;
	for my $ProtID (@{$gapSeqs}){
		$checkForStop = 0;
		$workingSeq = get_gene($genomeHash, $gbkHash,$ProtID,4*$noOfGaps,0); #The gene 
			#extended with the specified number of gaps
		$originalSeq = get_gene($genomeHash, $gbkHash,$ProtID,0,0); #The non-extended gene
   		$candidateSeq = $originalSeq; 
		unless ((my $sub = substr $workingSeq, 4*$noOfGaps) eq $originalSeq) { #unless get_gene 
				#returns two completely different sequences (parts without extension)
			print "\n\n get_gene mismatch \n\n";
			next;
		}
		$origLength = length($originalSeq);
		while ($check == 1){
			$workingSeq =~ m/([AGT]TG\D+)/;
			$workingSeq = $1;
			$lengthDiffStart = length($workingSeq) - $origLength;
			if ($lengthDiffStart == 0) { #The extended sequence is the same as the original sequence
				$extSeqs{$ProtID} = translateToProt($candidateSeq);				
				last;
			}
			if ($lengthDiffStart < 0) {
				print "\n Too short workingSeq: lengthDiffStart = $lengthDiffStart \n originalSeq: \n $originalSeq \n\n workingSeq: \n $workingSeq \n";
			}
			if ($lengthDiffStart % 3 == 0) { #The new ATG is in 
					#frame with the old one
				@arr = findInternalStop($workingSeq, $candidateSeq);
				#If an in frame stop codon was found
				if ($arr[2] == 1) {
					$workingSeq = $arr[0];
				}
				#If no in frame stop codon was found
				elsif ($arr[2] == 2) {
					$workingSeq = $arr[0];
					$candidateSeq = $arr[1];
					last;
				}
				#If the only stop corresponds to the last stop in the gene
				elsif ($arr[2] == 3) {
					$workingSeq = $arr[0];
					$extSeqs{$ProtID} = translateToProt($arr[0]);
					last
				}
				else {
					print "\n err \n";
				}
			}
			else { #If the found start codon is not in fram with the original one
				#Keep looking for a closer, better start codon
				$workingSeq = substr $workingSeq, 1;
			}
		}
	}
	return %extSeqs;
}

#Help routine to findAltStart
sub findInternalStop { #Used to find a stop codon internal to the extended part of a sequence in the findAltStart sub.
	my $wSeq = $_[0];
	my $cSeq = $_[1];
	my @output = ($wSeq, $cSeq, 0);
	my $len;
	my $lengthDiffStartStop;

	$wSeq =~ m/(T[AG][AG]\D+)$/g;
	$len = length($wSeq);
	if ($len != 0) { #The stop codon is not the last one in the sequence				
		$lengthDiffStartStop = length($wSeq) - $len; #Length difference 
						#between new start and in-between stop
		if (($lengthDiffStartStop % 3 == 0) && ($lengthDiffStartStop > 0)) { #The 
			#stop codon is in frame with the new start (and downstream of it)
			#Keep looking for a closer, better start codon	
			$wSeq = substr $wSeq, (length($wSeq)-$len);
			@output = ($wSeq, $cSeq, 1);
			return @output;
		}
		else { #If no downstream stop codon if found that is in frame with
				# the new start
				#Return the translated extended sequence $extSeqs{$ProtID} = translateToProt($workingSeq);
				$cSeq = $wSeq;						
				$wSeq = substr $wSeq, 1;						
				@output = ($wSeq, $cSeq, 2);
				return @output;
		}		
	}
	else { #If it is the last stop codon
			#Return the translated extended sequence
			#$extSeqs{$ProtID} = translateToProt($wSeq);
			@output = ($wSeq, $cSeq, 3);
			return @output;
	}
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

1;
