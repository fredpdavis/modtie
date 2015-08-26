#!/usr/local/bin/perl
=pod

=head1 NAME

calc_p_b_prior.pl - script to calculate P(B) prior from the Yeast-GFP data

what is the probability that two members of interacting family pairs actually interact?

our estimate: what is the probability that two members of interacting family pairs co-localize?

=head1 VERSION

fpd060121_1740

=cut

use strict;
use warnings;

use modtie qw/getstdinput calc_p_b_prior/;

main() ;

sub main {

   my $input = getstdinput({mode => 'calc_p_b_prior', ARGV => \@ARGV}) ;

   calc_p_b_prior($input) ;

}
