#!/usr/local/bin/perl
=pod

=head1 NAME

calc_model_pdbstart.pl - greps first ATOM record residue number

=head1 VERSION

fpd100925_2017 

=cut

use strict;
use warnings;
use modtie qw/getstdinput model_2_pdbstart/ ;

main() ;

sub main {

   my $input = getstdinput({mode => 'model_2_pdbstart'}) ;
   model_2_pdbstart($input) ;

}
