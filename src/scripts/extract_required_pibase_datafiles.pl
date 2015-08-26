#!/usr/local/bin/perl
=pod

=head1 NAME

extract_required_pibase_datafiles.pl - extract required data files PIBASE

=head1 VERSION

fpd100921_0935 

=cut

use strict;
use warnings;

use modtie qw/extract_required_pibase_datafiles/ ;

main() ;

sub main {

   modtie::extract_required_pibase_datafiles() ;

}
