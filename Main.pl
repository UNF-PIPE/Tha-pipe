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

#Get the parameters for ParserOrthoMCLgroups
my $speciesPerLine;
my $proteinsPerSpecies;
my $orthomcl_groups_file;
GetOptions ("min|minSpeciesPerLine=s" => \$speciesPerLine, 'max|maxProteinsPerSpecies=s' => \$proteinsPerSpecies, "g|groupsfile=s" => \$orthomcl_groups_file);

#Test. The "wrap around" code goes here
open(my $IN, $orthomcl_groups_file); #Open filehandle to the groups file
chomp(my @groups = <$IN>); #Read the file into an array
my %orthoHash = ParserOrthoMCLgroups(\@groups,$speciesPerLine,$proteinsPerSpecies);
my %prot = proteoms("/home/data/NCBI-proteoms/");

while (my ($key,$value) = each %orthoHash) {
	my $alignSequence = multipleAlign(\@{$value}, \%prot); 
	makeTree($alignSequence);
}

#&multAlign("/home/data/NCBI-proteoms/", \%orthoHash);
#print $orthoHash{1244}[1][0] . "\n";

#my %gbkhash = &make_gbk_hash("/home/data/NCBI-annotation/NC_008783.gbk");
#print &get_gene("/home/data/NCBI-genomes/Bartonella_bacilliformis.fasta", \%gbkhash, $ARGV[0]);

sub get_gene{ 
    #Arguments: 1 = path to genome, 2 = %gbkhash, 3 = protein_id
    my $genome = Bio::SeqIO->new(-file => $_[0], -format => "fasta"); #Makes a SeqIO-object from the genome
    my $gene;

    while (my $seq = $genome->next_seq){
	#Start and stop positions from %gbkhash
	my $start = $_[1]{$_[2]}[0];
	my $stop = $_[1]{$_[2]}[1];

	#If gene is on complementary strand, take the reverse complement of the seq
	if($_[1]{$_[2]}[2] == 0){
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

sub ParserOrthoMCLgroups {
	#takes as input the lines of a OrthoMCL groups output file in an array with one group/line per element, each as space-separated strings
	#returns those groups passing filters
	my ($groupsRef,$speciesPerLine,$proteinsPerSpecies) = @_;
	my @groups = @{$groupsRef};

	my @passedGroups = ();

	#loop over groups, for each one check if it passes filters
	GROUP: foreach my $group (@groups) {
		my ($prefix,@members) = split(/\s/,$group);
		my %speciesCount = ();	#collect the species we've seen in this group: value is number of proteins we've seen

		#loop over group members
		MEMBER: foreach my $member (@members) {
			my ($species,$protein) = split(/\|/,$member);
			if (defined($speciesCount{$species})) {
				if ($speciesCount{$species} >= $proteinsPerSpecies) {	#if we've already reached the max number of proteins for this species; fail group
					next GROUP;
				} else {
					$speciesCount{$species}++;
				}
			} else {
				$speciesCount{$species} = 1;
			}
		}
		
		unless(keys(%speciesCount) >= $speciesPerLine) { next GROUP; }	#if we haven't found enough species; fail group

		#if we reach this point, we have passed the filters - add this group to the passed group array
		push(@passedGroups,$group);
	}

	#return @passedGroups;
	my %orthoHash;
	foreach my $group (@passedGroups) {
		#Stores positives in a 2d hasharray [X][Y] where Y=0-1 for a spec. X gives the organism-protein pair
		$group =~ /my_prefix(\d+)\:/;
	 	my $orthoID = $1;
	 	my @group = split(/\||\s/,$group);
		my $lineSize = @group;
		my @clusterArray;
	    
	 	for(my $i=1;$i <=($lineSize-2); $i += 2){
				push(@clusterArray, ["$group[$i]" , "$group[$i+1]"]);	
	    }
	    $orthoHash{$orthoID} = \@clusterArray; #makes the hash an array ref.
	}
	return %orthoHash
}

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
	my $ortseq = "";			#Reset string
	my $alignStr = "";			#Reset string
	foreach my $ortrow  (@{$_[0]}) { #Loop array containing orthologs  
		$ortseq = $ortseq . "\>gi\|@{$ortrow}[1]\|@{$ortrow}[0]\n"; 	#Append fasta header containing species name and gene-ID 
		$ortseq = $ortseq .  " $_[1]{@{$ortrow}[1]} \n"; 		#Append sequence
	}
	#print $ortseq . "\n";
	my @ortalign = qx(echo  '$ortseq' \| muscle -quiet);	#Pipe ortholog fasta into clustalO
	foreach my $rows (@ortalign) {
		$alignStr =  $alignStr . $rows;	#Append all rows into a single string
	}
	return $alignStr;
}

sub makeTree {
	my $io = IO::String->new($_[0]);	#Convert string into io-object
	my $alnio = Bio::AlignIO->new(-fh => $io, -format=>'fasta'); #Make NJ tree
	my $dfactory = Bio::Tree::DistanceFactory->new(-method => 'NJ');
	my $stats = Bio::Align::ProteinStatistics->new;
	my $treeout = Bio::TreeIO->new(-format => 'newick');
	while( my $aln = $alnio->next_aln ) {
		my $mat = $stats->distance(-method => 'Kimura', -align  => $aln);
		my $tree = $dfactory->make_tree($mat);
		my $treeOut;
		Bio::TreeIO->new( -format => 'newick', -fh => IO::String->new(\$treeOut)  )->write_tree( $tree ); #Save tree to $treeOut
		print $treeOut . "\n";
	}
}

sub findGaps {
	my $limit = $_[1];
	my $io = IO::String->new($_[0]);
	my $alignment = Bio::SeqIO->new(-fh => $io, -format => 'fasta');
	
	while(my $sequence = $alignment->next_seq()){
		my $seq_string = $sequence->seq;
		if (($seq_string =~ tr/-/-/) <= $limit) {
			print ">" . $sequence->id . "\n";
			print $seq_string . "\n\n";
		}
	}
}
