#!/usr/bin/perl -w

#reads OrthoMCL cluster/group file on STDIN, passes it to ParserOrthoMCLgroups

use strict;
use Getopt::Long;

my $speciesPerLine;
my $proteinsPerSpecies;
GetOptions ("minSpeciesPerLine=s" => \$speciesPerLine, 'maxProteinsPerSpecies=s' => \$proteinsPerSpecies);

chomp(my @input = <STDIN>);
#my ($speciesPerLine,$proteinsPerSpecies) = @ARGV;	#the minimum number of species required per line, and the maximum number of proteins allowed per species

my @passedGroups = ParserOrthoMCLgroups(\@input,$speciesPerLine,$proteinsPerSpecies);

foreach my $group (@passedGroups) {
	print "$group\n";
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
			
		
