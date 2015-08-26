#!/usr/local/bin/perl
=pod

=head1 NAME

modbaseformat.mt.pl - reformats MODTIE display_complexes() output
to modtie database format

=head1 VERSION
fpd060316_1607

=cut

use strict;
use warnings;
use modtie qw/getstdinput format_modbaseoutput/ ;

main() ;

sub main {

   my $input = getstdinput({ARGV => \@ARGV, mode => 'format_modbaseoutput'}) ;
   format_modbaseoutput($input) ;

}
