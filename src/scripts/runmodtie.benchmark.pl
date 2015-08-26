#!/usr/local/bin/perl
=pod

=head1 NAME

modtie_benchmark.pl - script to benchmark modtie potentials

=head1 VERSION

fpd060121_1740

=cut

use strict;
use warnings;

use modtie qw/getstdinput runmodtie_benchmark/;

main() ;

sub main {

   my $input = getstdinput({mode => 'runmodtie_benchmark', ARGV => \@ARGV}) ;
   runmodtie_benchmark($input);

}
