#!/usr/local/bin/perl
=pod

=head1 NAME

modtie_roc.pl - script to roc modtie potentials

=head1 VERSION

fpd060121_1740

=cut

use strict;
use warnings;

use modtie qw/getstdinput runmodtie_roc/;

main() ;

sub main {

   my $input = getstdinput({mode => 'runmodtie_roc', ARGV => \@ARGV}) ;
   runmodtie_roc($input);

}
