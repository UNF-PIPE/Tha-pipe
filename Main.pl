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
use POSIX qw(strftime);

use findParalogs qw( findParalogs );
use HashRoutines qw( mkHash make_gbk_hash parse_orthomcl_file );
use MapRoutines qw( get_gene );
use multipleAlign qw(multipleAlign);
use findAltStart qw( findAltStart findGaps );

#Declare parameters
my ($speciesPerLine, $proteinsPerSpecies, $orthomcl_groups_file, $proteoms, $genomes, $annotations, $out, $speciePos_genome, $speciePos_proteome, $gapLength, $cores, $logfile, $altStart);

#Set default values
$proteoms = "t/test_data/proteome_data/NCBI-proteoms/";
$genomes = "t/test_data/genome_data/NCBI-genomes/";
$annotations = "t/test_data/genome_data/NCBI-annotation/";
$speciePos_genome = 3; #The position in the genome fasta headers that contain the specie identifier
$speciePos_proteome = 3; #The position in the proteome fasta headers that contain the specie identifier
$gapLength = 13; #Sequences startin g with more than this number of gaps will be sent to findAltStart
$cores = `cat /proc/cpuinfo | grep 'processor'| wc -l`; #on Linux
#$cores = `sysctl -n hw.ncpu`; #on Mac
my $date = strftime "%F", localtime;
$logfile = $date . ".log";
$altStart = 1;

GetOptions ("min|minSpeciesPerLine=s" => \$speciesPerLine, 'max|maxProteinsPerSpecies=s' => \$proteinsPerSpecies, "g|groupsfile=s" => \$orthomcl_groups_file, "p|proteoms=s" => \$proteoms, "gs|genomes=s" => \$genomes, "a|annotations=s" => \$annotations, "o|outname=s" => \$out, "hg|speciePos_genome=s" => \$speciePos_genome, "hp|speciePos_proteome=s" => \$speciePos_proteome,"gap|gapLength=s" => \$gapLength, "c|cores=s" => \$cores, "log|logfile=s" => \$logfile,"as|altStart=s" => \$altStart);

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
        my @groupLog; #Array to store log info per ortholog group
        
        #Align the sequences
	my $alignedSequences = multipleAlign(\@{$value}, \%prot);

        #Find sequences starting with more than n gaps. These will be passed on to finAltStart, which will try to find alternative start codons.
	my @gapSeqs = findGaps($alignedSequences, $gapLength);
        if(@gapSeqs){
                my $gapInfo = "Sequences starting with more than " . $gapLength . " nr of gaps: " . join(" ", @gapSeqs);
                push(@groupLog, $gapInfo);
        }
        
        #Look for alternative start codons. &findAltStart returns the aligned sequences (of which some may have been extended) as a hash. 
        #The id:s of those that actually were extended are returned in an array.
        my ($foundSeqs, $extSeqs) = findAltStart(\@gapSeqs, $gapLength, \%annotation, \%genome, $altStart);
        if(@{$extSeqs}){
               my $extInfo = "Sequences that were extended: " . join(" ", @{$extSeqs});
               push(@groupLog, $extInfo);
        }
        
        #Re-align the extended sequences
	$alignedSequences = multipleAlign(\@{$value}, \%prot, $foundSeqs); 
        
        # Check for paralogs
        my @Paralogs = findParalogs($alignedSequences, 2);
        my $outName;
        if(@Paralogs){
                my $parInfo = "Species having paralogous proteins in this group: " . join(" ", @Paralogs);
                push(@groupLog, $parInfo);
        }
        
        #Write to file
        $outName = $out . $key . ".fasta";
        open(my $outfile, ">", $outName);
        print $outfile $alignedSequences;

        #Store the logged info in the hash %log for later printing
        if(@groupLog){
                my $tmpFile = $out . $key . "_tmp.log";
                open(my $tmpLog, ">", $tmpFile);
                print $tmpLog $key . "\n";
                print $tmpLog join("\n", @groupLog) . "\n\n";
                #$log{$key} = \@groupLog;
        }

        $pm->finish; # Terminates the child process
}
$pm->wait_all_children;

#Write log info to file
$logfile = $out . $logfile; 
open (logFile, ">", $logfile);
my @outFiles = read_dir $out;
for(@outFiles){
        if ($_ =~ m/_tmp\.log$/){
                my $tempfile = $out . $_;
                open (TEMPFILE, "<", $tempfile);
                while(<TEMPFILE>){
                        print logFile $_; 
                }
                close(TEMPFILE);
                unlink $tempfile;
        }
}
close(logFile);
