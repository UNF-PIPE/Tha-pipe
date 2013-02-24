package findAltStart;
use strict;
use warnings;
use Bio::SeqIO;
use IO::String;
use File::Slurp;
use MapRoutines qw(get_gene translateToProt);
use Exporter;
use base qw(Exporter);
use Bio::Tree::Tree;
our @EXPORT_OK = qw(findAltStart findGaps);

sub findAltStart {
	my %out_sequences;
	my ($gapSeqs, $noOfGaps, $gbkHash, $genomeHash) = @_;
	my ($original_gene, $extended_gene, $extension);
	my @codon_arrays;
	my $suitable_start_pos;
		
	for my $ProtID (@{$gapSeqs}) {
		#Fetch gene and extension
		$original_gene = get_gene($genomeHash,$gbkHash,$ProtID,0,0);
		my $translated_orig = translateToProt($original_gene);
		if ($translated_orig =~ m/\*\D+/) {
			die "Error: Original gene contains internal stop codon, igonered";
		} 
		$extended_gene = get_gene($genomeHash,$gbkHash,$ProtID,3*$noOfGaps,0);
		$extension = substr($extended_gene,0,3*$noOfGaps);
		
		#Build start and stop codon arrays
		@codon_arrays = buildCodonArrays($extension);
		
		if (@codon_arrays) {	
			#Compare codon arrays to find a suitable start codon
			$suitable_start_pos = compareCodonArrays(\@codon_arrays);
		} else {
			$suitable_start_pos = -1;
		}

		#Modify extension if necessary. Translate and add to output hash
		if ($suitable_start_pos == -1) {
			$out_sequences{$ProtID} = $translated_orig;
		} else {
			my $translated_ext = translateToProt(substr($extended_gene, $suitable_start_pos));
			if ($translated_ext =~ m/\*\D+/) {
				die "Error: Extended protein contains internal stop codon, igonered";
			} else {
				$out_sequences{$ProtID} = $translated_ext;
			}
		}
	}
	return %out_sequences;
}

sub buildCodonArrays {
	my $extension = $_[0];
	my $position = 0;
	my $end = length($extension) - 1;
	my @start_codons;
	my @stop_codons;
	
	while ($position <= $end-3) {
		if (inFrame($position, $end+1)) {
			if ((substr($extension, $position, 3)) =~ m/([ATG]TG)/) {
				push(@start_codons, $position);}
			elsif ((substr($extension, $position, 3)) =~ m/(T[AG][AG])/) {
				push(@stop_codons, $position);}
		}
		$position += 1;
	}
	return (\@start_codons,\@stop_codons);
}	

#Checks if a codon in an extension is in frame with the original gene
sub inFrame {
	my $codon_start_position = $_[0];
	my $beginning_of_gene = $_[1]; #Needs to be end of extension +1

	return (($beginning_of_gene - $codon_start_position) % 3  == 0);
}

#Finds a suitable start codon from arrays of discovered start and stop codons. Returns -1 if none was found.
sub compareCodonArrays {
	my ($codon_array_ref) = shift;
	my @start_codons;
	my @stop_codons;
	my @closer_starts = (-1);
	if (defined $codon_array_ref->[0]) {
		@start_codons = @{$codon_array_ref->[0]};
	}
	if (defined $codon_array_ref->[1]) {
		@stop_codons = @{$codon_array_ref->[1]};
	}

	if ($#start_codons + 1 != 0) {
		if (@stop_codons == 0) {
			#Pick first start codon from 5' end of extension
			return $start_codons[0];
		}
		else { 
			#Pick first start codon closer to gene then the closest stop codon
			push(@closer_starts, grep {$_ > $stop_codons[-1]} @start_codons);
                        return $closer_starts[1];
		}
	}
	else {return -1};
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
