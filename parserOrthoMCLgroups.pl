#!/usr/bin/perl -w

#reads OrthoMCL cluster/group file on STDIN, passes it to ParserOrthoMCLgroups

use warnings;
use strict;
use Getopt::Long;

my $speciesPerLine;
my $proteinsPerSpecies;
my %orthoHash; #change name, stupid
GetOptions ("minSpeciesPerLine=s" => \$speciesPerLine, 'maxProteinsPerSpecies=s' => \$proteinsPerSpecies);

chomp(my @input = <STDIN>);
#my ($speciesPerLine,$proteinsPerSpecies) = @ARGV;	#the minimum number of species required per line, and the maximum number of proteins allowed per species

my @passedGroups = ParserOrthoMCLgroups(\@input,$speciesPerLine,$proteinsPerSpecies);

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

exit;

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

	return @passedGroups;
}
			
		
