package findParalogs;
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

use Exporter;
use base qw( Exporter );
our @EXPORT_OK = qw( findParalogs );

sub makeTree {
        #Input: A string with a multiple alignment
        my $io = IO::String->new($_[0]);        #Convert string into io-object
        my $alnio = Bio::AlignIO->new(-fh => $io, -format=>'fasta'); #Make AlignIO object
        my $dfactory = Bio::Tree::DistanceFactory->new(-method => 'NJ');
        my $stats = Bio::Align::ProteinStatistics->new;
        my $treeout = Bio::TreeIO->new(-format => 'newick');
        while( my $aln = $alnio->next_aln ) {
                my $mat = $stats->distance(-method => 'Kimura', -align  => $aln);
                my $tree = $dfactory->make_tree($mat);
                #my $treeOut = $tree->as_text('newick');
                #print $treeOut . "\n";
                return $tree;
        }
}

sub findParalogs {
	my $tree = makeTree($_[0]);
	my @leaves = $tree->get_leaf_nodes;
	my %species;
	my @paralogs;
	for(@leaves){
		my $id = $_->id;
		my $specie = (split(/\|/, $id))[1];
		$species{$specie}++;
	}
	for(keys %species){
		if($species{$_} > 1){
			if( checkIfParalog($tree,$_) ){
				push(@paralogs,$_);
			}
		}	
	}
	return @paralogs; 	
}

sub checkIfParalog{
	my ($tree, $species) = @_;
	my @leaves = $tree->get_leaf_nodes;
	my @parCandidates;
	for(@leaves){
                my $id = $_->id;
                my $parCandidate = (split(/\|/, $id))[1];
		if($parCandidate eq $species){
			push(@parCandidates, $_);
		}
        }
        my $LCA = $tree->get_lca(-nodes => \@parCandidates);

        for my $child ( $LCA->get_all_Descendents ) {
                if($child->is_Leaf){
                        my $id = $child->id;
                        $id = (split(/\|/, $id))[1];
                        if($id ne $species){
                                return 1;
                        }
                }
        }
        return 0;
}
1;
