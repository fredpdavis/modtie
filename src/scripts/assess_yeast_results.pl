#!/usr/local/bin/perl
=pod

=head1 NAME

assess_yeast_predictions.pl - compares MODTIE yeast predictions to experimental
interactions stored in BIND interactions

=head1 VERSION

fpd060131_1914

=cut

use strict;
use warnings;
use modtie qw/getstdinput assess_yeast_results/ ;

main() ;

sub main {

   my $input = getstdinput({ARGV => \@ARGV, mode => 'assess_yeast_results'}) ;
   assess_yeast_results($input) ;

}
