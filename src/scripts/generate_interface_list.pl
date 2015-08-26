#!/usr/local/bin/perl
=pod

=head1 NAME

generate_interface_list.pl - script to make template interface list from PIBASE

=head1 VERSION

fpd100921_0935 

=cut

use strict;
use warnings;

use modtie qw/generate_interface_list/ ;

main() ;

sub main {

   modtie::generate_interface_list() ;

}
