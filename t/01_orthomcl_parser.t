#!/usr/bin/perl
use strict;
use warnings;

use Test::More;
use Data::Printer;

use OrthoMCLParser qw( parse_groups );

my $test_data = [
    'my_prefix0001: Spec1|PR000 Spec2|PR100 Spec3|Pr200',
];

my %res = parse_groups( $test_data, 3, 3 );
ok( exists $res{'0001'}, 'Test data contains group 0001' );
ok( ref($res{'0001'}) eq 'ARRAY', 'Group 0001 is an array ref' );
ok( $res{'0001'}[0][0] eq 'Spec1', 'Index 0 is species' );
ok( $res{'0001'}[0][1] eq 'PR000', 'Index 1 is protein name' );

done_testing();

