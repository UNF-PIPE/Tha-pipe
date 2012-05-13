#!/usr/bin/perl

use strict;
use warnings;

use Test::More;

use create_gbkhash qw(make_gbk_hash);

my $test_data = "test_data/gbk_data/";


my %res = make_gbk_hash( $test_data  );


ok( exists $res{'YP_988349.1'}, "Test data contains group YP_988349.1"};


ok( ref($res{'YP_988349.1'}) eq 'ARRAY', 'Group YP_988349.1 is an array ref' ); 

ok( $res{'YP_988349.1'}[0]  == 10700,   'Group YP_988349.1 has a first entry of 10700');


done_testing();
