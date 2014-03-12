#!perl -T

use Test::More tests => 1;

BEGIN {
    use_ok( 'deltaMass' ) || print "Bail out!\n";
}

diag( "Testing deltaMass $deltaMass::VERSION, Perl $], $^X" );
