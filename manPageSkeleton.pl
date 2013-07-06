#!/usr/bin/perl
use strict;
use Getopt::Long;
use Pod::Usage;

if (!defined($ARGV[0])) {
        exec( 'perldoc', $0 ); 
}

=head1 NAME

Tha pipah

=head1 SYNOPSIS

Main.pl [-min int] [-max int] [-g file] [-h int] [-p directory] [-gs directory] [-a directory] [-o directory] 

=head1 DESCRIPTION

-

=head1 OPTIONS

=head2 -min, -minSpeciesPerLine

Minimum species per ortholog group in the ortholog file.

=head2 -max, -maxProteinsPerSpecies

Maximum proteins per species.

=head2 -g, -groupsfile

File that contains all the orthologs from OrthoMCL.

=head2 -p, -proteoms

Proteom fasta file directory. One file for each species.

=head2 -gs, -genomes

Genome fasta file directory. One file for each species.

=head2 -a, -annotations

Annotation file directory (.gbk), One file for each species.

=head2 -o, -outname

Existing output folder for multiple alignments.

=head2 -h, -specieHeader

Position of species name in the fasta header.

h|specieHeader_pos

=head1 FILES

=head1 REQUIREMENTS 

-

=head1 EXAMPLES

perl Main.pl -min 4 -max 1 -g groups.txt -p ./t/test_data/proteome_data/NCBI-proteoms/ -gs t/test_data/genome_data/NCBI-genomes/ -a t/test_data/genome_data/NCBI-annotation/ -o /home/simon/projects/pipan/results/2013-07-06_paralog_bugFix/test_res_max1/

=head1 AUTHOR

Simon Forsberg, Joakim Karlsson, Bjorn Viklund and Matilda Aslin

=cut
