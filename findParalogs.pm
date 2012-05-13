package findParalogs;
use strict;
use warnings;

use Exporter;
use base qw( Exporter );
our @EXPORT_OK = qw( makeTree findParalogs );

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
	
}
