#!/usr/local/bin/perl
=pod

=head1 NAME

cut_target_domains.pl - carves domains out of model pdb files

=head1 VERSION

fpd051219_0222

=cut

use strict;
use warnings;
use modtie qw/cut_domains/ ;

main() ;

sub main {

   cut_domains() ;

}
