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
use mkHash qw(mkHash);
use multipleAlign qw(multipleAlign);

#Get the parameters for ParserOrthoMCLgroups
my $speciesPerLine;
my $proteinsPerSpecies;
my $orthomcl_groups_file;
GetOptions ("min|minSpeciesPerLine=s" => \$speciesPerLine, 'max|maxProteinsPerSpecies=s' => \$proteinsPerSpecies, "g|groupsfile=s" => \$orthomcl_groups_file);

#Test. The "wrap around" code goes here
my %orthoHash = parse_orthomcl_file(
    $orthomcl_groups_file, $speciesPerLine, $proteinsPerSpecies
);

my %prot = mkHash("/home/data/NCBI-proteoms/");
my %genome = mkHash("/home/data/NCBI-genomes/");

my $cores = `cat /proc/cpuinfo | grep 'processor'| wc -l`;
my $pm = new Parallel::ForkManager($cores++);
while (my ($key,$value) = each %orthoHash) {
	my $pid = $pm->start and next;

	my $alignSequence = multipleAlign(\@{$value}, \%prot); 
	#my $alignSequence = multipleAlign(\@{$value}, \%prot, \%INSERT NEW SEQUENCE HASH HERE <-----); 
	if(my @gap_ids = findGaps($alignSequence,11)){
		#print $gap_ids[0] . "\n";
		#print $alignSequence . "\n";
	}

	$pm->finish; # Terminates the child process
	#makeTree($alignSequence);
}
$pm->wait_all_children;

#&multAlign("/home/data/NCBI-proteoms/", \%orthoHash);
#print $orthoHash{1244}[1][0] . "\n";

#my %gbkhash = &make_gbk_hash("/home/data/NCBI-annotation/NC_008783.gbk");
#print &get_gene("/home/data/NCBI-genomes/Bartonella_bacilliformis.fasta", \%gbkhash, $ARGV[0]);

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
