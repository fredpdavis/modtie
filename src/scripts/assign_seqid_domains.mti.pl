#!/usr/local/bin/perl
=pod

=head1 NAME

assign_seq_domains.pl - assigns SCOP domains to modbase sequences

=head1 VERSION

fpd051208_0257

=cut

use strict;
use warnings;
use modtie qw/getstdinput seqid_2_domains/ ;

main() ;

sub main {

   my $input = getstdinput({mode => 'seqid_2_domains'}) ;
   seqid_2_domains($input) ;

}
