#!/usr/bin/perl -w
use strict;

my @array = ("a","b","c","d","e");
print "Length: " . ($#array+1) . "\n";

for my $i (0 .. $#array-1) {
	for my $j ($i+1 .. $#array) {
		print "$array[$i] vs $array[$j]" . "\n";	
	}
}
