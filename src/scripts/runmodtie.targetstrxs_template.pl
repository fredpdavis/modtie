#!/usr/local/bin/perl
=pod

=head1 NAME

runmodtie.targetstrxs_template.pl - script to score a single complex given domain definitions

=head1 VERSION

fpd060107_1441

=cut

use strict;
use warnings;

use modtie qw/getstdinput runmodtie_targetstrxs_template/;

main() ;

sub main {

   my $input = getstdinput({mode => 'runmodtie_targetstrxs_template', ARGV => \@ARGV}) ;
   runmodtie_targetstrxs_template($input);

}
