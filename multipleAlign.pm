package multipleAlign;
use strict;
use warnings;
use Bio::SeqIO;
use IO::String;

use Exporter;
use base qw( Exporter );
our @EXPORT_OK = qw(multipleAlign);

sub multipleAlign {
	my ($orthologsRef, $proteomes, $new_sequences) = @_;
	my $ortseq = "";			#Reset string
	my $alignStr = "";			#Reset string
	foreach my $ortrow  (@{$orthologsRef}) { #Loop array containing orthologs  
		$ortseq = $ortseq . "\>gi\|@{$ortrow}[1]\|@{$ortrow}[0]\n"; 	#Append fasta header containing species name and gene-ID 
		if (ref $new_sequences && exists $new_sequences->{$ortrow->[1]} ) {
			$ortseq .= $new_sequences->{$ortrow->[1]} . "\n";
		}
		else {
			$ortseq .= $proteomes->{$ortrow->[1]} . "\n"; 
		}
	}
	my @ortalign = qx(echo  '$ortseq' \| muscle -quiet);	#Pipe ortholog fasta into clustalO
	foreach my $rows (@ortalign) {
		$alignStr =  $alignStr . $rows;	#Append all rows into a single string
	}
	return $alignStr;
}

