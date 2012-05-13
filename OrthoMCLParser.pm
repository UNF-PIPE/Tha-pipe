package OrthoMCLParser;
use strict;
use warnings;

use Exporter;
use base qw( Exporter );
our @EXPORT_OK = qw( parse_orthomcl_file );


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

    my @passedGroups = ();

    # loop over groups, for each one check if it passes filters
    GROUP: foreach my $group (@$groups) {
        my ($prefix, @members) = split /\s/, $group;
        # collect the species we've seen in this group: value is number of
        # proteins we've seen
        my %speciesCount = ();

        # loop over group members
        MEMBER: foreach my $member (@members) {
            my ($species, $protein) = split /\|/,$member;
            $speciesCount{$species}++;
            # if we've already reached the max number of proteins for this
            # species; fail group
            if ( $speciesCount{$species} >= $proteinsPerSpecies ) {
                next GROUP;
            }
        }

        # if we haven't found enough species; fail group
        unless ( keys %speciesCount >= $speciesPerLine ) {
            next GROUP;
        } 

        # if we reach this point, we have passed the filters - add this group
        # to the passed group array
        push @passedGroups, $group;
    }

    #return @passedGroups;
    my %orthoHash;
    foreach my $group (@passedGroups) {
        # Stores positives in a 2d hasharray [X][Y] where Y=0-1 for a spec. X
        # gives the organism-protein pair
        $group =~ /my_prefix(\d+)\:/;
        my $orthoID = $1;
        my @group = split /\||\s/, $group;
        my $lineSize = @group;
        my @clusterArray;
        
        for(my $i=1;$i <=($lineSize-2); $i += 2){
                push @clusterArray, ["$group[$i]" , "$group[$i+1]"];
        }
        $orthoHash{$orthoID} = \@clusterArray; #makes the hash an array ref.
    }
    return %orthoHash
}

1; # Return true at the end
