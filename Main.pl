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

#Get the parameters for ParserOrthoMCLgroups
my ($speciesPerLine, $proteinsPerSpecies, $orthomcl_groups_file, $proteoms, $genomes, $annotations, $out, $speciePos);
#Set default values
$speciePos = 2;
GetOptions ("min|minSpeciesPerLine=s" => \$speciesPerLine, 'max|maxProteinsPerSpecies=s' => \$proteinsPerSpecies, "g|groupsfile=s" => \$orthomcl_groups_file, "p|proteoms=s" => \$proteoms, "gs|genomes=s" => \$genomes, "a|annotations=s" => \$annotations, "o|outname=s" => \$out, "h|specieHeader_pos=s" => \$speciePos);

#Parse the orthomcl-file and create a hash with all ortholog groups that meet the criteria. 
my %orthoHash = parse_orthomcl_file($orthomcl_groups_file, $speciesPerLine, $proteinsPerSpecies);

#Make hashes for proteoms, genomes and annotation files. Keys are gene IDs, species and gene IDs respectively
my %prot = mkHash($proteoms);
my %genome = mkHash($genomes);
my %annotation = make_gbk_hash($annotations);
#my %prot = mkHash("t/test_data/proteome_data/NCBI-proteoms/");
#my %genome = mkHash("t/test_data/genome_data/NCBI-genomes/");
#my %annotation = make_gbk_hash("t/test_data/genome_data/NCBI-annotation/"); 

#Get the number of processor cores and create a forkmanager object for multithreading usage.
my $cores = `cat /proc/cpuinfo | grep 'processor'| wc -l`; #on Linux
#my $cores = `sysctl -n hw.ncpu`; #on Mac
my $pm = new Parallel::ForkManager($cores++);

while (my ($key,$value) = each %orthoHash) {
	my $pid = $pm->start and next;
	my $alignedSequences = multipleAlign(\@{$value}, \%prot);

        #Find sequences starting with more than n gaps. These will be passed on to finAltStart, which will try to find alternative start codons.
#	my @gapSeqs = findGaps($alignedSequences, 13);
#    	my %extSeqs = findAltStart(\@gapSeqs, 13, \%annotation, \%genome);
        
        #Re-align the extended sequences
#	$alignedSequences = multipleAlign(\@{$value}, \%prot, \%extSeqs); 
        
        # Check for paralogs
        my @Paralogs = findParalogs($alignedSequences, $speciePos);
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
