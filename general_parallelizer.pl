#!/usr/bin/perl
# Author: Martin Dahlo
#
# Usage: perl scriptname.pl <data-to-be-processed file> <number of parts (nodes)> <output directory>
#
# Ex.
# perl scriptname.pl projectMap/sampleId/rawData/reads.fq 64 projectMap/sampleId/alignment/
# or
# perl scriptname.pl readsToSnpcal 11 snpResults

use warnings;
use strict;
use POSIX qw/ceil/;
use File::Path;
use Cwd;





# should these be set by arguments at execution?
######################################################################
##########      SETTINGS: MAKE SURE THESE ARE CORRECT       ##########
######################################################################
my $projId = "b2010074"; # UPPMAX project id
my $hours = "24"; # number of hours to book the nodes for





=pod
This scipt will take any file, divide it into a specified number of equaly sized parts
and submit one job to the queue system for each part.

Simply write the functional content of the sbatch file to specify what will be done with the divided file.

Pseudo code:

* count the number of lines in a file
* calculate how many lines there should be in each part of the file
* split the file in as many parts as specified
* create and submit a sbatch file for each part of the file


=cut

# Usage message
my $usage = "Usage: perl scriptname.pl <data-to-be-processed file> <number of parts (nodes)> <output directory>\n";


# get the file name
my $infile = shift or die $usage;  
my $parts = shift or die $usage;
my $outdir = shift or die $usage;  

# adjust paths to be general
$infile = absPath($infile); # adjust to absoulte path
$outfile = remSlash(absPath($outdir)); # adjust to absolute path, then remove trailing slashes


# create the output dir if needed
mkpath("$outdir/tempfiles");

# count the number of lines
open(WC, "wc -l $infile |") or die "Failed: $!\n";
my ($lines) = <WC> =~ /^(\d+)\s/;


# calculate the correct number of Lines Per Part
my $lpp = ceil($lines/$parts) + (4 - (ceil($lines/$parts)%4)); # the number of lines per part, rounded upward towards the the first number evenly divided by 4 (to avoid breaking fastq files, where 1 read = 4 lines)


# split the file, placing the part in $outdir/tempfiles with the prefix 'tempfile'
system("split --lines=$lpp $infile $outdir/tempfiles/tempfile");


# loop and process each part
foreach my $file (<$outdir/tempfiles/*>){
  # get the filename only (without path)
  my $fileName =  pop @{[ split /\//, $file ]}; # get the file name of the current file
  
  # write a submission file
  open SB, ">", "$outdir/tempfiles/sbatch" or die $!;
  
  print SB <<EOF;
#! /bin/bash -l
#SBATCH -A $projId
#SBATCH -p node -n 1
#SBATCH -t $hours:00:00
#SBATCH -J $fileName
#SBATCH -o $outdir/$fileName.out
#SBATCH -e $outdir/$fileName.err
# more options if wanted

# load the modules you need
# Ex.
# module load bioinfo-tools
# module load bowtie


# the commands you want performed on each file
# Ex.
# cd /home/myDir/glob/projectX
# bowtie -C --sam -p 8 myIndexPrefix -f $file $outdir/$fileName.aligned.sam

EOF

  # submit the file
  close(SB);
  system("sbatch $outdir/tempfiles/sbatch");

}

# print to screen
print "\nAll jobs submitted. Results for each part will be stored in $outdir\n";





####   Subroutines   ####
 
# Removes the / sign, if any, at the end of directory names
sub remSlash{
  my $str = shift;
  if($str =~ m/\/$/){
    chop($str);
  }
  return $str;
}

# Adjust relative path to absolute paths
sub absPath{
  my $str = shift;
  if($str !~ m/^\//){ # if the first char is not a / ie. a relative path
    $str = &Cwd::cwd()."/$str"; # attach the current work directory infront of the relative path
  }
  return $str;  
}
