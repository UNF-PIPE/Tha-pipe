#!/usr/bin/perl

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

use findParalogs qw( findParalogs );
use HashRoutines qw( mkHash make_gbk_hash parse_orthomcl_file );
use MapRoutines qw( get_gene );
use multipleAlign qw(multipleAlign);
use findAltStart qw( findAltStart findGaps );

#Declare parameters
my ($speciesPerLine, $proteinsPerSpecies, $orthomcl_groups_file, $proteoms, $genomes, $annotations, $out, $speciePos_genome, $speciePos_proteome, $gapLength, $cores);

#Set default values
$proteoms = "t/test_data/proteome_data/NCBI-proteoms/";
$genomes = "t/test_data/genome_data/NCBI-genomes/";
$annotation = "t/test_data/genome_data/NCBI-annotation/";
$speciePos_genome = 3; #The position in the genome fasta headers that contain the specie identifier
$speciePos_proteome = 3; #The position in the proteome fasta headers that contain the specie identifier
$gapLength = 13; #Sequences startin g with more than this number of gaps will be sent to findAltStart
$cores = `cat /proc/cpuinfo | grep 'processor'| wc -l`; #on Linux
#$cores = `sysctl -n hw.ncpu`; #on Mac

GetOptions ("min|minSpeciesPerLine=s" => \$speciesPerLine, 'max|maxProteinsPerSpecies=s' => \$proteinsPerSpecies, "g|groupsfile=s" => \$orthomcl_groups_file, "p|proteoms=s" => \$proteoms, "gs|genomes=s" => \$genomes, "a|annotations=s" => \$annotations, "o|outname=s" => \$out, "hg|speciePos_genome=s" => \$speciePos_genome, "hp|speciePos_proteome=s" => \$speciePos_proteome,"gap|gapLength=s" => \$gapLength, "c|cores=s" => \$cores);

#Parse the orthomcl-file and create a hash with all ortholog groups that meet the criteria. 
my %orthoHash = parse_orthomcl_file($orthomcl_groups_file, $speciesPerLine, $proteinsPerSpecies);

#Make hashes for proteoms, genomes and annotation files. Keys are gene IDs, species and gene IDs respectively
my %prot = mkHash($proteoms,$speciePos_proteome);
my %genome = mkHash($genomes,$speciePos_genome);
my %annotation = make_gbk_hash($annotations);

#Get the number of processor cores and create a forkmanager object for multithreading usage.
my $pm = new Parallel::ForkManager($cores++);

while (my ($key,$value) = each %orthoHash) {
	my $pid = $pm->start and next;
	my $alignedSequences = multipleAlign(\@{$value}, \%prot);

        #Find sequences starting with more than n gaps. These will be passed on to finAltStart, which will try to find alternative start codons.
	my @gapSeqs = findGaps($alignedSequences, $gapLength);
    	my %extSeqs = findAltStart(\@gapSeqs, $gapLength, \%annotation, \%genome);
        if(@gapSeqs){
                print "GAP!\t$key\n";
        }
        
        #Re-align the extended sequences
	$alignedSequences = multipleAlign(\@{$value}, \%prot, \%extSeqs); 
        
        # Check for paralogs
        my @Paralogs = findParalogs($alignedSequences, 2);
        my $outName;
        if(@Paralogs){
                #$outName = "/home/simon/Research/phylogeny_pipeline/results/130325/$key" . "_par.fasta";
                $outName = $out . $key . "_par.fasta";
                print "PARALOG!\t$outName\n";
                print join("\n", @Paralogs) . "\n\n";
        }
        else{
                $outName = $out . $key . ".fasta";
                #$outName = "/home/simon/Research/phylogeny_pipeline/results/130325/$key" . ".fasta"; 
        }
        
        #Write to file
        open(my $outfile, ">", $outName);
        print $outfile $alignedSequences;

        #my $tree = makeTree($alignedSequences);
        #my $treeOut = $tree->as_text('newick');
        #print $treeOut . "\n";
	$pm->finish; # Terminates the child process
}
$pm->wait_all_children;
