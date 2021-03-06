package OrthoMCLParser;
use strict;
use warnings;

use Exporter;
use base qw( Exporter );
our @EXPORT_OK = qw( parse_orthomcl_file parse_groups );


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

1; # Return true at the end
