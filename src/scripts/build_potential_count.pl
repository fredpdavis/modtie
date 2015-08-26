#!/usr/local/bin/perl
=pod

=head1 NAME

build_potential_counts.pl - countup statistics for modtie statistical potential

=head1 VERSION
fpd060117_0156

=cut

use strict;
use warnings;
use modtie qw/getstdinput buildpotential_count/ ;

main() ;

sub main {

   my $input = getstdinput({mode => 'buildpotential_count'}) ;
   buildpotential_count($input) ;

}
