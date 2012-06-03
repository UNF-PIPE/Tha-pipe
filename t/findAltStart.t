#!/usr/bin/perl
use strict;
use warnings;

use Bio::SeqIO;
use IO::String;
use Test::More;

use findAltStart qw( findAltStart findGaps mkHash);

my $test_data_path = "t/test_data/findAltStart_data/";
my %gbkTestHash = (
	"YP_003037636.1" => ["3590035", "3590829", "0" ], 
	"YP_989117.1" => ["848590", "849420", "0", "", "","CP000524.1"], 
	"YP_191040.1" => ["647552", "648352", "1"],
	"YP_003225146.1" => ["8889", "9719", "0"],
);
my $YP_9891171seq = "MVDYIEYNETPLKFNGKIRIFDDYAFAEMRKVGQIAAECLDALTDIIKPGITTQEIDDFIFIFGAERGALPADLNYRGYSHSCCTSINHVVCHGIPNKKSLQEGDIVNVDVTFILNGWHGDSSRMYPVGKVKRAAERLLEITHECLMRGIEAVKPGATTGDIGAAIQRYAESERCSVVRDFCGHGIGQLFHDAPNILHYGNPGEGEELKQGMIFTIEPMINLGKPQVKILSDGWTAVTRDRSLSAQYEHTIGVTDQGCEIFTQSPKNIFYIPNSCA";

#Opens the alignment file and concatenates it into one string.
open my $aln_file, '<', $test_data_path . "test2_aln.fasta" or die "Can't find test file: $?, $!";
my $aln = join("", <$aln_file>);

#Run the subroutines
my %genomeHash = mkHash($test_data_path . "genomes/");
my @gapSeqs = findGaps($aln, 13);
my %extSeqs = findAltStart(\@gapSeqs, 13, \%gbkTestHash, \%genomeHash);

#print join(" ", @gapSeqs);
#print $extSeqs{"YP_989117.1"} . "\n";
#print $aln;

#Translate the retrieved DNA seq back to protein
my $extendedSeq = Bio::PrimarySeq->new ( -seq => $extSeqs{"YP_989117.1"} ,
                                         -id  => 'YP_989117.1',
                                         #-accession_number => 'X78121',
                                         -alphabet => 'dna',
                                         -is_circular => 0 );
my $extendedSeq_tr = $extendedSeq->translate;


ok(join("", @gapSeqs) =~ "YP_989117.1", "The seq YP_989117.1 contains more than 13 gaps");
ok(!(join("", @gapSeqs) =~ "YP_003225146.1"), "The seq YP_003225146.1 contains less than 13 gaps");
ok($extendedSeq_tr->seq =~ $YP_9891171seq, "The seq retrieved from the genome for YP_989117.1 is correct");

done_testing();
