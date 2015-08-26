#!/usr/local/bin/perl
=pod

=head1 NAME

databaseid_2_seqid.pl -  given a list of database ids (eg uniprot, GI codes), returns modbase seq_id

=head1 VERSION

fpd060520_1730

=cut

use strict;
use warnings;
use DBI ;

use modtie qw/databaseid_2_seqid/;

main() ;

sub main {

#   my $input = getstdinput({mode => 'runmodtie_xset', ARGV => \@ARGV}) ;
#   runmodtie_xset($input);

   databaseid_2_seqid() ;

}
