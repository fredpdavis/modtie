#!/usr/local/bin/perl
=pod

=head1 NAME

build_potential_postcalc.pl - post-process counted statistics for modtie statistical potential

=head1 VERSION

fpd060117_1156

=cut

use strict;
use warnings;
use modtie qw/getstdinput buildpotential_postcalc/ ;

main() ;

sub main {

   my $input = getstdinput({mode => 'buildpotential_postcalc', ARGV => \@ARGV}) ;
   buildpotential_postcalc($input) ;

}
