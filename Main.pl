#!/usr/local/bin/perl

use strict;
use warnings;
use Bio::SeqIO;
use Getopt::Long;

#Get the parameters for ParserOrthoMCLgroups
my $speciesPerLine;
my $proteinsPerSpecies;
my $orthomcl_groups_file;
GetOptions ("minSpeciesPerLine=s" => \$speciesPerLine, 'maxProteinsPerSpecies=s' => \$proteinsPerSpecies, "groupsfile=s" => \$orthomcl_groups_file);

#Test. The "wrap around" code goes here
open(my $IN, $orthomcl_groups_file); #Open filehandle to the groups file
chomp(my @groups = <$IN>); #Read the file into an array
my %orthoHash = ParserOrthoMCLgroups(\@groups,$speciesPerLine,$proteinsPerSpecies);
print $orthoHash{1244}[1][1] . "\n";

#my %gbkhash = &make_gbk_hash("/home/data/NCBI-annotation/NC_008783.gbk");
#print &get_gene("/home/data/NCBI-genomes/Bartonella_bacilliformis.fasta", \%gbkhash, $ARGV[0]);

#To here

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
