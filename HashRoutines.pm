#!/usr/local/bin/perl
package HashRoutines;
use strict;
use warnings;
use File::Slurp;
use Bio::SeqIO;
use Getopt::Long;

use Exporter;
use base qw( Exporter );
our @EXPORT_OK = qw(make_gbk_hash mkHash parse_orthomcl_file parse_groups);

#Subroutine for making a hash of genbank annotation files	
sub make_gbk_hash{	
	my $genbankDir = $_[0]; #One example gbk file is downloadable from ftp://ftp.ncbi.nih.gov/genomes/Bacteria/Bartonella_bacilliformis_KC583_uid58533/
	my @files = read_dir($genbankDir);
	my $all_lines;
	my %gbkhash;

	my $ACID; #Accesion ID
	my $start; #start location for ACID;
	my $stop; #stop location for ACID;
	my $complement; #values true or false
	my $GI; # GI number
	my $geneID; #GeneID
	my $ge;

	foreach my $file (@files){
		#print $file."\n";
		open GBFILE, $genbankDir.$file or die $!;		
		my @array = <GBFILE>;
		chomp(@array);
		#my $toprow = shift(@array);
	  my @toprows = @array[0..10];	
#		if($toprow =~ m/^LOCUS\s*(\S*)\s/) {
#			$ge = "$1";
#		}
          foreach (@toprows) {
                  if($_ =~ m/^VERSION\s*(\S*)\s/) {
                          $ge = "$1";
                  }	
          }

		#my @toprow_split = split(/\t/,$toprow);
		#my $genomeID = $toprow_split[2];
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
			$gbkhash{$ACID}=[$start, $stop, $complement, $GI, $geneID, $ge];
		}
	}
	return %gbkhash
}#endbracket of subroutine make_gbk_hash


sub mkHash {
	my %hash; 				#Predefine hash	
	my ($dir,$speciesPos) = @_; 			#Get arguments
	my @files = read_dir($dir); 	#Contains all files in the directory 
	for my $file ( @files ) { 		#Loop trough all files 
        	my $filePath = "$dir$file"; #Current path to file
        	my $hash = Bio::SeqIO->new(-file => "<$filePath", -format => 'fasta'); #Insert file into SeqIO object
        	while (my $accession1 = $hash->next_seq) { # Go through the SeqIO objects
               		my $header = $accession1->id;
               		my $seq = $accession1->seq;
			            my @split = split(/\|/,$header);
               		$hash {$split[$speciesPos]} = $seq; 
       		}
	}
	return %hash;
}

sub parse_orthomcl_file {
    # File name is first argument, remove that from @_
    my $file = shift;

    # Open filehandle to the groups file
    open my $IN, '<', $file or die "Can't find orthomcl file: $?, $!";

    # Read the file into an array
    chomp(my @groups = <$IN>);

    # takes as input the lines of a OrthoMCL groups output file in an array
    # with one group/line per element, each as space-separated strings returns
    # those groups passing filters. @_ contains the rest of the parameters.
    return parse_groups( \@groups, @_ );
}

sub parse_groups {
    my ($groups, $speciesPerLine, $proteinsPerSpecies) = @_;
    my %orthoHash;

    # loop over groups, for each one check if it passes filters
    GROUP: foreach my $group (@$groups) {
        my ($prefix, @members) = split /\s+/, $group;
        # collect the species we've seen in this group: value is number of
        # proteins we've seen
        my %speciesCount = ();

        # loop over group members
        MEMBER: foreach my $member (@members) {
            my ($species, $protein) = split /\|/,$member;
            $speciesCount{$species}++;
            # if we've already reached the max number of proteins for this
            # species; fail group
            if ( $speciesCount{$species} > $proteinsPerSpecies ) {
                next GROUP;
            }
        }

        # if we haven't found enough species; fail group
        unless ( keys %speciesCount >= $speciesPerLine ) {
            next GROUP;
        } 

        # if we reach this point, we have passed the filters - add this group
        # to the orthoHash
        my ($orthoID) = $prefix =~ /(\d+):/;

        my @clusterArray;
        for my $protein ( @members ) {
            push @clusterArray, [ split /\|/, $protein ];
        }
        $orthoHash{$orthoID} = \@clusterArray; # makes the hash an array ref.
    }

    return %orthoHash
}
1;


