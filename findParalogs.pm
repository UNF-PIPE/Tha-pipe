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
	my $tree = makeTree($_[0]);
	my @leaves = $tree->get_leaf_nodes;
	my %species;
	my @paralogs;
	for(@leaves){
		my $id = $_->id;
		my $specie = split(/\|/, $id)[1];
		$species{$specie}++;
	}
	for(keys %species){a
		if($species{$_}==2){
			if( checkIfParalog($tree,$_) ){
				push(@paralogs,$_);
			}
		}	
		elsif($species{$_}>2){
			if( checkIfParalogsComplicated($tree, $_) ){
				push(@paralogs, $_);
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
                my $parCandidate = split(/\|/, $id)[1];
		if($parCandidate eq $species){
			push(@parCandidates, $_);
		}
        }

	if($parCandidates[0]->ancestor eq $parCandidates[1]->ancestor){
		return 0;
	}
	return 1; 


}



sub checkIfParalogsComplicated{
	my ($tree, $species) = @_;
	my @leaves = $tree->get_leaf_nodes;
	my @parCandidates;
	for(@leaves){
                my $id = $_->id;
                my $parCandidate = split(/\|/, $id)[1];
		if($parCandidate eq $species){
			push(@parCandidates, $_);
		}
        }
	
	for my $i (0 .. (length($parCandidates)-1){
		for my $j (($i+1) .. length($parCandidates)){
			print $i.$j."\n";
		}
	}	

}






for my $i (0 .. (5-1){
	for my $j (($i+1) .. 5){
		print $i.$j."\n";
	}
}	

