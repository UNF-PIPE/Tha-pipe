package mkHash;
use strict;
use warnings;
use Bio::SeqIO;
use IO::String;
use File::Slurp;

use Exporter;
use base qw( Exporter );
our @EXPORT_OK = qw(mkHash);

sub mkHash {
	my %hash; 				#Predefine hash	
	my $dir = $_[0]; 			#Define directory
	my @files = read_dir($dir); 	#Contains all files in the directory 
	for my $file ( @files ) { 		#Loop trough all files 
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
