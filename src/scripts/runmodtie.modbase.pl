#!/usr/local/bin/perl
=pod

=head1 NAME

modtie_full.pl - script to make simplified modtie prediction calls

=head1 VERSION

fpd051207_0021
orig:fpd051121_2056

=cut

use strict;
use warnings;

use modtie qw/getstdinput runmodtie_modbase/;

main() ;

sub main {

   my $input = getstdinput({mode => 'runmodtie_modbase', ARGV => \@ARGV}) ;
   runmodtie_modbase($input);

}
