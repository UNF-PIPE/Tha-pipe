#!/usr/bin/perl
use strict;
use warnings;

use Test::More;

use OrthoMCLParser qw( parse_groups );

my $test_data = [
    'my_prefix0001: Spec1|PR000 Spec2|PR100 Spec3|Pr200',
    'no_standard0002: Spec1|PR001 Spec2|PR101 Spec3|Pr201',
];

my %res = parse_groups( $test_data, 3, 1 );

ok( exists $res{'0001'}, 'Test data contains group 0001' );
ok( ref($res{'0001'}) eq 'ARRAY', 'Group 0001 is an array ref' );
is( $res{'0001'}[0][0], 'Spec1', 'Index 0 is species' );
is( $res{'0001'}[0][1], 'PR000', 'Index 1 is protein name' );

ok( exists $res{'0002'}, 'Test data contains group 0002 (non standard prefix)' );

done_testing();

