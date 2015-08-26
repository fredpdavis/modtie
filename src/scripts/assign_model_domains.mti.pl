#!/usr/local/bin/perl
=pod

=head1 NAME

assign_model_domains.pl - assigns SCOP domains to modbase models

=head1 VERSION

fpd051207_0021
orig:fpd051121_2056

=cut

use strict;
use warnings;
use modtie qw/getstdinput model_2_domains/ ;

main() ;

sub main {

   my $input = getstdinput({mode => 'model_2_domains'}) ;
   model_2_domains($input) ;

}
