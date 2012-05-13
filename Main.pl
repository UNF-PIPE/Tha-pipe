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
		
	}

	$pm->finish; # Terminates the child process
	#makeTree($alignSequence);
}
$pm->wait_all_children;
#&multAlign("/home/data/NCBI-proteoms/", \%orthoHash);
#print $orthoHash{1244}[1][0] . "\n";

#my %gbkhash = &make_gbk_hash("/home/data/NCBI-annotation/NC_008783.gbk");
#print &get_gene("/home/data/NCBI-genomes/Bartonella_bacilliformis.fasta", \%gbkhash, $ARGV[0]);

sub get_gene{
    my ($genome_path,$orthoHash,$ID,$ext_start,$ext_stop) = @_;
    #The soubroutine will retrieve a subseq (e.g. a gene) from a genome
    #Arguments: 1 = path to genome, 2 = %gbkhash, 3 = protein_id, 4 = extension of startposition backwards (optional), 5 = extension of stop position (optional)
    my $genome = Bio::SeqIO->new(-file => $genome_path, -format => "fasta"); #Makes a SeqIO-object from the genome
    my $gene;

    while (my $seq = $genome->next_seq){
        #Start and stop positions from %gbkhash
        my $start = $orthoHash->{$ID}[0] - $ext_start;
        my $stop = $orthoHash->{$ID}[1] + $ext_stop;

        #If gene is on complementary strand, take the reverse complement of the seq
        if($orthoHash->{$ID}[2] == 0){
            $gene = $seq->subseq($start, $stop);
        }
        else{
            my $gene_obj = $seq->trunc($start, $stop);
            my $reversed = $gene_obj->revcom;
            $gene = $reversed->seq;
        }
    }
    return $gene; 
}


#Subroutine for making a hash of genbank annotation file	
sub make_gbk_hash{	
	my $genbankfile = $_[0]; #downloadable from ftp://ftp.ncbi.nih.gov/genomes/Bacteria/Bartonella_bacilliformis_KC583_uid58533/
	my $all_lines;
	my %gbkhash;

	my $ACID; #Accesion ID
	my $start; #start location for ACID;
	my $stop; #stop location for ACID;
	my $complement; #values true or false
	my $GI; # GI number
	my $geneID; #GeneID

	open GBFILE, $genbankfile or die $!;		
	my @array = <GBFILE>;
	chomp(@array);
	$all_lines = join("", @array);
	#print $all_lines;
	my @array2= split(/\/translation/,$all_lines);
	foreach(@array2){
		if ($_ =~ m/CDS\s+(\d+)\.\.(\d+).\s+/){
			#print "$1\t";
			$start = $1;
			$stop = $2;
			$complement = 0;	
		}
		elsif ($_ =~ m/CDS\s+complement\((\d+)\.\.(\d+)\)\s+/) {
			#print "$1 $2\t";
			$start = $1;
			$stop = $2;
			$complement = 1;
		}
		if ($_ =~ m/\/protein_id=\"(.*)\"\s+\/db_xref=\"GI:(.*)\"\s+\/db_xref=\"GeneID:(.*)\"/) {
			#print "$1\n";
			$ACID = $1;			
			$GI = $2;
			$geneID = $3;
		}
		#Store in hash with AccessionID as key
		$gbkhash{$ACID}=[$start, $stop, $complement, $GI, $geneID];
	}
	return %gbkhash
} #endbracket of subroutine make_gbk_hash


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
