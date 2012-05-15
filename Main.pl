#!/usr/local/bin/perl

use strict;
use warnings;
use Bio::SeqIO;
use Getopt::Long;
use File::Slurp;
use Bio::TreeIO;
use IO::String;
use Bio::AlignIO;
use Bio::Align::ProteinStatistics;
use Bio::Tree::DistanceFactory;
use Bio::Tree::TreeI;
use Parallel::ForkManager;


use OrthoMCLParser qw( parse_orthomcl_file );
use findParalogs qw( findParalogs );
use gbkHashOps qw( make_gbk_hash get_gene );

#Get the parameters for ParserOrthoMCLgroups
my $speciesPerLine;
my $proteinsPerSpecies;
my $orthomcl_groups_file;
GetOptions ("min|minSpeciesPerLine=s" => \$speciesPerLine, 'max|maxProteinsPerSpecies=s' => \$proteinsPerSpecies, "g|groupsfile=s" => \$orthomcl_groups_file);

#Test. The "wrap around" code goes here
my %orthoHash = parse_orthomcl_file(
    $orthomcl_groups_file, $speciesPerLine, $proteinsPerSpecies
);

my %prot = proteoms("/home/data/NCBI-proteoms/");

my $cores = `cat /proc/cpuinfo | grep 'processor'| wc -l`;
my $pm = new Parallel::ForkManager($cores++);
while (my ($key,$value) = each %orthoHash) {
	my $pid = $pm->start and next;

	my $alignSequence = multipleAlign(\@{$value}, \%prot); 
	#my $alignSequence = multipleAlign(\@{$value}, \%prot, \%INSERT NEW SEQUENCE HASH HERE <-----); 
	if(my @gap_ids = findGaps($alignSequence,11)){
		print $gap_ids[0] . "\n";
		print $alignSequence . "\n";
	}

	$pm->finish; # Terminates the child process
	#makeTree($alignSequence);
}
$pm->wait_all_children;
#&multAlign("/home/data/NCBI-proteoms/", \%orthoHash);
#print $orthoHash{1244}[1][0] . "\n";

#my %gbkhash = &make_gbk_hash("/home/data/NCBI-annotation/NC_008783.gbk");
#print &get_gene("/home/data/NCBI-genomes/Bartonella_bacilliformis.fasta", \%gbkhash, $ARGV[0]);


sub proteoms {
	my %proteoms; 				#Predefine proteoms hash	
	my $protDir = $_[0]; 			#Define proteome directory
	my @files = read_dir($protDir); 	#Contains all files in prot directory 
	for my $file ( @files ) { 		#Loop trough all proteome files 
        	my $filePath = "$protDir$file"; #Current path to proteome file
        	my $proteome = Bio::SeqIO->new(-file => "<$filePath", -format => 'fasta'); #Insert file into SeqIO object
        	while (my $accession1 = $proteome->next_seq) { # Go through the SeqIO objects
               		my $header = $accession1->id;
               		my $seq = $accession1->seq;
			my @split = split(/\|/,$header);
               		$proteoms {$split[3]} = $seq; #Save "my_prefixXXXX as key, sequence as value
       		}
	}
	return %proteoms;
}

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
