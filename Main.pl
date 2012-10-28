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

use findParalogs qw( findParalogs makeTree );
use HashRoutines qw( mkHash make_gbk_hash parse_orthomcl_file );
use MapRoutines qw( get_gene );
use multipleAlign qw(multipleAlign);
use findAltStart qw( findAltStart findGaps );

#Get the parameters for ParserOrthoMCLgroups
my $speciesPerLine;
my $proteinsPerSpecies;
my $orthomcl_groups_file;
GetOptions ("min|minSpeciesPerLine=s" => \$speciesPerLine, 'max|maxProteinsPerSpecies=s' => \$proteinsPerSpecies, "g|groupsfile=s" => \$orthomcl_groups_file);

#Test. The "wrap around" code goes here
my %orthoHash = parse_orthomcl_file($orthomcl_groups_file, $speciesPerLine, $proteinsPerSpecies);
my %prot = mkHash("/home/data/NCBI-proteoms/");
my %genome = mkHash("/home/data/NCBI-genomes/");
my %annotation = make_gbk_hash("/home/data/NCBI-annotation/"); 

#my $cores = `cat /proc/cpuinfo | grep 'processor'| wc -l`;
#my $pm = new Parallel::ForkManager($cores++);
while (my ($key,$value) = each %orthoHash) {
	#my $pid = $pm->start and next;
	my $alignedSequences = multipleAlign(\@{$value}, \%prot); 
	my @gapSeqs = findGaps($alignedSequences, 13);
    	my %extSeqs = findAltStart(\@gapSeqs, 13, \%annotation, \%genome);
	$alignedSequences = multipleAlign(\@{$value}, \%prot, \%extSeqs); 

	#$pm->finish; # Terminates the child process
	my $tree = makeTree($alignedSequences);
    	my $treeOut = $tree->as_text('newick');
    	print $treeOut . "\n";
}
#$pm->wait_all_children;
