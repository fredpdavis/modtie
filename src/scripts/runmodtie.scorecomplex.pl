#!/usr/local/bin/perl
=pod

=head1 NAME

runmodtie.scorecomplex.pl - script to score a single complex given domain definitions

=head1 VERSION

fpd060107_1441

=cut

use strict;
use warnings;

use modtie qw/getstdinput runmodtie_scorecomplex/;

main() ;

sub main {

   my $input = getstdinput({mode => 'runmodtie_scorecomplex', ARGV => \@ARGV}) ;
   runmodtie_scorecomplex($input);

}
