=head1 NAME

pibase.pm - PIBASE functions necessary for MODTIE

=head1 DESCRIPTION

The pibase.pm perl library is a subset of routines from the original PIBASE
perl library that are used for MODTIE interaction with PIBASE files

=head1 AUTHOR

Fred P. Davis, HHMI-JFRC (davisf@janelia.hhmi.org)

=head1 LICENCE AND COPYRIGHT

Copyright 2005,2010 Fred P. Davis (davisf@janelia.hhmi.org).
See the file COPYING for copying permission.

This file is part of MODTIE.

MODTIE is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

MODTIE is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MODTIE.  If not, see <http://www.gnu.org/licenses/>.

=head1 SUBROUTINES

=cut


package modtie::pibase ;

use strict ;
use warnings ;
use modtie ;
use vars qw($VERSION @ISA @EXPORT @EXPORT_OK) ;

require Exporter ;
use Carp qw/croak/ ;
use File::Temp qw/tempfile tempdir/;
use File::Path qw/mkpath/ ;

@ISA = qw/Exporter/ ;
@EXPORT_OK = qw/connect_pibase/ ;
push @EXPORT_OK, qw/safe_copy safe_move sid_2_domdir replace_undefs/;
push @EXPORT_OK, qw/subset_extract get_salign parse_modeller_ali/ ;
push @EXPORT_OK, qw/todload_bdp_ids/ ;


=head2 connect_pibase($dbspecs)

   Title:       connect_pibase()
   Function:    Connects to the pibase database.
   Args:        $_->{db}	database name
                $_->{user}	user name
                $_->{pass}	password
   Returns:     DBI database handle to pibaes

=cut

sub connect_pibase {

   require DBI ;

#   my $modtie_specs = shift ;
#   modtie::complete_modtie_specs($modtie_specs) ;
   my $modtie_specs = modtie::set_modtie_specs() ;
   my $specs = $modtie_specs->{pibase_specs} ;

   my $curhost = Sys::Hostname::hostname()  ;
   my $dbname ;

   if ($specs->{host} =~ /^$curhost/) {
      $dbname = "DBI:mysql:database=$specs->{db}";
   } else {
      $dbname = "DBI:mysql:database=$specs->{db}".';'."host=$specs->{host}"; } 

   if (exists $specs->{mysql_socket}) {
      $dbname.=';mysql_socket='.$specs->{mysql_socket};}

   my $dbh = DBI->connect( $dbname, $specs->{user}, $specs->{pass},
                           {RaiseError => 1, AutoCommit => 1} ) ;

   return ($dbh, $specs->{db}, $specs->{root}) ;

}


=head2 mysql_fetchcols(dbh, query)

   Function:    Processes an n column query and returns a list of array references, where each array reference holds all the valeus for a given column.
   Args:        $_[0] = dbh - DBI database handle
                $_[1] = SQL select query

   Returns:     @_ - list of arrayref
                $a[i]->[j] ith column, jth row

=cut

sub mysql_fetchcols {

   my $dbh = shift ;
   my $query = shift ;
   my $vals = shift ;

   my $sth ;
   if ($query->can("execute")) {
      $sth = $query ;
      $sth->execute(@{$vals}) ;
   } else {
      $query = qq{$query} ;
      $sth = $dbh->prepare($query) ;
      $sth->execute() ;
   }

   my @results ;
   while (my @currow = $sth->fetchrow()) {
      foreach my $j (0 .. $#currow) {
         push @{$results[$j]}, $currow[$j] ; } }

   return @results ;

}



=head2 mysql_singleval(dbh, query)

   Function:    Processes a 1 column, 1 row query and returns a scalar
                  containing the value.
   Args:        $_[0] = dbh - DBI database handle
                $_[1] = query - SQL format
   Return:      $a - scalar
                $a = "result"

=cut

sub mysql_singleval {

   my ($dbh, $query) = @_ ;

   $query = qq{$query} ;
   my $sth = $dbh->prepare($query) ;
   $sth->execute() ;

   my $ans = $sth->fetchrow() ;

   return $ans ;

}



=head2 sid_2_domdir()

   Title:       sid_2_domdir()
   Function:    returns the directory name where the PDB file of the
      specified domain resides.
   Args:        $_ = subset_id
   Returns:     directory name

=cut

sub sid_2_domdir {

   my $in = shift ;
   my $sid = $in->{sid} ;

   my $modtie_specs ;
   if (exists $in->{specs}) {
      $modtie_specs = $in->{specs} ;
   } else {
      $modtie_specs = modtie::set_modtie_specs();
   }
   my $specs = $modtie_specs->{pibase_specs} ;


   my ($bdp) = ($sid =~ /BDP([0-9]+)/) ;
   my $dirnum = POSIX::floor($bdp / 100) ;
   my $dir = $specs->{subsets_dir}."/$dirnum/$bdp" ;

   return $dir ;

}



=head2 get_weirdo_resser()

   Title:       get_weirdo_resser()
   Function:    residue_info() altered to read in non-standard amino acids
                  (HETATOM OR ATOM) in as regular amino acids according to
                  MODELLER's restyp.lib mapping
   Args:        $_->{pdb_fn} - pdb file name
                $_->{outfile} - output file name [optional]
                $_->{identifier} - bdp_id identifier [optional]

   Returns: (if outfile is not specified)
      $_->{chain_no}->[i] - serial chain number
      $_->{chain_id}->[i] - chain identifier
      $_->{resno_serial}->[i] - serial residue number
      $_->{resno}->[i] - residue number
      $_->{resno_int}->[i] - integer portion of residue number
      $_->{resna}->[i] - residue name
      $_->{chain_type}->[i] - chain_type ('p' or 'n')
      $_->{weird_res_ser}->{resno_serial."\n".chain_no} = 1 if resna = MSE or MEX
      $_->{weird_res_raw}->{resno."\n".chain_id} = 1 if resna = MSE or MEX

   Output file: Residue listing ($_->{outfile}
      1. bdp identifier
      2. serial chain number
      3. chain identifier
      4. serial residue number
      5. residue number
      6. integer portion of residue number
      7. residue name
      8. chain type ('p' or 'n')

=cut

sub get_weirdo_resser {
#IMPORTED FROM: pibase::PDB::residues::get_weirdo_resser()

   my $params = shift ;

   my $specs;
   if (!exists $params->{specs}) {
      $specs = modtie::set_modtie_specs() ;
   } else {
      $specs = $params->{specs} ;
   }

   if (!exists ($params->{pdb_fn})) {
      return (0, 'PDB file not specified]') ; }

   my $pdb_file = $params->{pdb_fn} ;

   my $outfile ;
   if (exists $params->{outfile}) {
      $outfile = $params->{outfile} ; }

   my $identifier ;
   if (exists $params->{identifier}) {
      $identifier = $params->{identifier} ;
   } else {
      $identifier = 'bdp_id' ;
   }

# Build hash with single letter abbreviations for standard residue names from PDB ATOM records.

   my %rescode = (
      'ALA' => 'A' ,
      'ARG' => 'R' ,
      'ASN' => 'N' ,
      'ASP' => 'D' ,
      'CYS' => 'C' ,
      'GLN' => 'Q' ,
      'GLU' => 'E' ,
      'GLY' => 'G' ,
      'HIS' => 'H' ,
      'HSD' => 'H' ,
      'ILE' => 'I' ,
      'LEU' => 'L' ,
      'LYS' => 'K' ,
      'MET' => 'M' ,
      'PHE' => 'F' ,
      'PRO' => 'P' ,
      'SER' => 'S' ,
      'THR' => 'T' ,
      'TRP' => 'W' ,
      'TYR' => 'Y' ,
      'VAL' => 'V',
      'UNK' => 'X',
      '  C' => 'c',
      '  G' => 'g',
      '  A' => 'a',
      '  T' => 't',
      '  U' => 'u',
      '  I' => 'i',
      ' +C' => 'c',
      ' +G' => 'g',
      ' +A' => 'a',
      ' +T' => 't',
      ' +U' => 'u',
      ' +I' => 'i',
      'MSE' => 'M',
      'MEX' => 'C',
   ) ;

# got last 2 mappings from src/commands/read_model.F90
# special mappings from modlib/restyp.lib
#      'CSH' => 'C', 
#      'PR0' => 'P', 
#      'PRZ' => 'P',  #heterocyclic aromatic compound - not really a pro, but never used anyways - 
#      'ASX' => 'B', 
#      'GLX' => 'Z', 
##      'CSS' => 'C',  - not used see 1h32 models
#      'CYX' => 'C', 
#      'MSE' => 'M', 
#      'MEX' => 'C', 

   my @chain_id ;
   my @start_res ;
   my @end_res ;
   my @num_atm ;
   my @num_res ;
   my @num_het ;
   my @chain_seq ;
   my @chain_type ;

   my @chain_resno ; # array of residue list
   my @chain_resno_int ;
   my @chain_resna ;

   my $lastchain = '' ;
   my $lastresno = '' ;

# Read a line from the pdb file.

   if ($pdb_file =~ /\.gz$/) {
      open(PDBFILE, $specs->{binaries}->{zcat}." ".$pdb_file." | ") ;
   } else {
      open (PDBFILE, $pdb_file) ;
   }
   while (my $line = <PDBFILE>) {

      chomp $line;

      my $rectype ;

# If the line is an ATOM or HETATM record, set rectype to the record type.
# - unless an MSE or MEX entry which is set as an ATOM record.

      if (($line =~ /^ATOM/) ||
          ( ($line =~ /^HETATM/) &&
	    length($line) > 20 &&
	    ((substr($line, 17, 3) eq 'MSE')||
	     (substr($line, 17, 3) eq 'MEX')) )) {
         $rectype = 'ATOM';
      } elsif ($line =~ /^HETATM/) {
         $rectype = 'HETATM'; }

# If record type has been set,

      if (defined $rectype) {

# Extract:
# * chain id (22)
# * residue number (23 - 27) - note: includes insertion code
# * residue name (18 - 20)

         my $curchain = substr($line, 21, 1) ;
         my $curresno = substr($line, 22, 5) ; $curresno =~ s/ //g ;
         my $curresna = substr($line, 17, 3) ;

         my $rescode ;
         if (exists ($rescode{$curresna})) {
   	    $rescode = $rescode{$curresna} ; }

# Initialize the new chain flag.
# If the current chain id is different from the last chain id, set the new chain flag.

         my $newchain = 0;
         if ($curchain ne $lastchain) {
   	    $newchain = 1; }

# Initialize the new residue flag.
# If the current residue number is different from the last residue number OR the new chain flag has been set, set the new residue flag.

         my $newres = 0;
         if (($curresno ne $lastresno) ||
	     $newchain) {
            $newres = 1; }

# If the new chain flag has been set,


         if ($newchain) {

# Create an entry in the chain arrays for the current chain.
# Store the current chain id.

   	    push @chain_id, $curchain ;

# Initialize the start residue, end residue, chain seqeuncee, number of residues, atoms, and HETATMS (defaults 0), and chaintype (default 'p')

   	    push @start_res, $curresno ;
   	    push @end_res, $curresno ;
   	    push @chain_seq, "" ;
   	    push @num_res, 0 ;
   	    push @num_atm, 0 ;
   	    push @num_het, 0;
   	    push @chain_type, 'p' ;

            push @chain_resno, [ ] ;
            push @chain_resno_int, [ ] ;
            push @chain_resna, [ ] ;

	 }


# If the new reisdue flag has been set and the current line is an ATOM record,
         if ($newres && ($rectype eq 'ATOM')) {

# Set the end residue number of the current chain to the current residue number.
   	    $end_res[$#end_res] = $curresno ;

# If the residue name is known (single letter abbreviation known),
   	    if (defined $rescode) {

# If the residue code contains any lower case letters (nucleotides), set the chain type to 'n'.
               if ($rescode =~ /[a-z]/) {
   	          $chain_type[$#chain_type] = 'n' ; }

# Concatenate the current residue letter abbreviation to the stored sequence for this chain.

	       push @{$chain_resno[$#chain_type]}, $curresno ;
	       my ($curresno_int, undef) = residue_int($curresno) ;
	       push @{$chain_resno_int[$#chain_type]}, $curresno_int ;
	       push @{$chain_resna[$#chain_type]}, $curresna ;

   	       $chain_seq[$#chain_id] .= $rescode{$curresna} ;

# Increment the number of residues in this chain.

   	       $num_res[$#num_res]++ ;

	    }

	 }

# If the record type is 'ATOM' and the single letter abbreviation is known, increment the number of atoms.


         if (($rectype eq 'ATOM') &&
	     (defined $rescode)) {
   	    $num_atm[$#num_atm]++; }

# Elsif the record type is 'HETATM', increment the number of HETATMs.

         elsif ($rectype eq 'HETATM') {
   	    $num_het[$#num_het]++; }

# Store the current residue number as the last residue number.
# Store the current chain id as the last chain id.

         $lastresno = $curresno ;
         $lastchain = $curchain ;

      }

# Elsif the line is an ENDMDL entry, stop reading from the pdb file.

      elsif ($line =~ /^ENDMDL/) {
         last ; }

   }
   close(PDBFILE) ;

# Initialize the chain and residue count.
   my $ch_count = 1;

# Open the output file.
   my $outfile_fh ;
   if (defined $outfile) {
      if ($outfile !~ /^\>/) {
         $outfile = '>'.$outfile ; }
      open($outfile_fh, $outfile) ; }

# Iterate through the chains.
   my $results ;


   my $weird_res_ser ;
   my $weird_res_raw ;
   foreach my $j (0 .. $#chain_id) {

# If the number of atoms is greater than 0,

      if (($num_atm[$j] > 0 ) && (defined $chain_resno_int[$j]->[0])) {

	 my $res_count = 1 ;

# If the chain sequence is '', goto the next chain.

         if ($chain_seq[$j] eq '') {
            next; }

# Designate the fields to be displayed based on the current display format.

         foreach my $k ( 0 .. $#{$chain_resno[$j]} ) {
	    if (defined $outfile_fh) {
	       my @outvals = ( $identifier, $ch_count, $chain_id[$j],$res_count,
	                    $chain_resno[$j]->[$k], $chain_resno_int[$j]->[$k],
			    $chain_resna[$j]->[$k], $chain_type[$j] ) ;
               print $outfile_fh join("\t", @outvals)."\n" ;
	    } else {
	       push @{$results->{chain_no}}, $ch_count ;
	       push @{$results->{chain_id}}, $chain_id[$j] ;
	       push @{$results->{resno_serial}}, $res_count;
	       push @{$results->{resno}}, $chain_resno[$j]->[$k];
	       push @{$results->{resno_int}}, $chain_resno_int[$j]->[$k];
	       push @{$results->{resna}}, $chain_resna[$j]->[$k];
	       push @{$results->{chain_type}}, $chain_type[$j]  ;

	       if ($chain_resna[$j]->[$k] eq 'MSE' ||
	           $chain_resna[$j]->[$k] eq 'MEX') {
	          $weird_res_ser->{$res_count."\n".$ch_count}++ ;
	          $weird_res_raw->{$chain_resno[$j]->[$k]."\n".$chain_id[$j]}++ ;
	       }
	    }
	    $res_count++ ;
	 }

# Increment the chain count.

         $ch_count++;

      }
   }
   if (defined $outfile_fh) { close($outfile_fh) ; }

   $results->{weird_res_ser} = $weird_res_ser ;
   $results->{weird_res_raw} = $weird_res_raw ;
   return ($results) ;

}


=head2 SUB subset_extract()

   Function: extracts specified residues/chains from PDB file
   Args:        $_[0] = PDB file name
                $_[1] = output PDB file name
                $_[2] = chain identifier
                $_[3] = start residue number
                $_[4] = end residue number

   Returns:     $_->[i] arrayref of errors

   Files IN:    PDB file ($_[0])
   Files OUT:   subset PDB file ($_[1])

=cut

sub subset_extract {
#IMPORTED FROM: pibase::PDB::subsets::subset_extract()

   my $in = shift ;
   my $pdb_fn = $in->{in_fn} ;
   my $cutpdb_fn= $in->{out_fn} ;
   my $chain_id = $in->{chain} ;
   my $resno_start = $in->{start} ;
   my $resno_end = $in->{end} ;

   my $binaries = modtie::locate_binaries() ;

   my ($temp_fh, $localbdp) = tempfile(SUFFIX => ".pdb") ; close($temp_fh) ;
   if ($pdb_fn =~ /\.gz$/) {
      my $tcom = $binaries->{zcat}." $pdb_fn > ".$localbdp ;
      system($tcom) ;
   } else {
      safe_copy($pdb_fn, $localbdp) ;
   }

   if (!-e $localbdp) {
      return {error_fl=> "ERROR: pdb file access error. couldnt copy locally"};}


   my $errors ;

   if ($binaries->{'subset_extractor'} eq 'ERROR') {
      croak("ERROR subset_extract(): subset_extractor binary not found") ; }

   if ($pdb_fn eq $cutpdb_fn) {
      croak("ERROR subset_extract(): ".
            "the original and cut pdb filenames are the same") ; }

   my $compress_out_fl = 0 ;
   if ($cutpdb_fn =~ /\.gz$/) {
      $compress_out_fl = 1 ;
      $cutpdb_fn =~ s/\.gz$// ;
   }

   my ($subsetdef_fh, $subsetdef_fn) = tempfile("subsetdef.XXXXXX") ;
   my $errfile_fn = $subsetdef_fn.".err" ;

   foreach my $j ( 0 .. $#{$chain_id} ) {
      my @outvals = ($chain_id->[$j], $resno_start->[$j], $resno_end->[$j]) ;
      print $subsetdef_fh join("\t", @outvals)."\n" ;
   }
   close ($subsetdef_fh) ;

   my $tcom = "$binaries->{subset_extractor} $localbdp ".
              "< $subsetdef_fn 2> $errfile_fn >$cutpdb_fn" ;
   system($tcom) ;
   if ($compress_out_fl) {
      system("gzip ".$cutpdb_fn) ;
      $cutpdb_fn .= '.gz' ;
   }
   unlink $subsetdef_fn ;

   if (-s $errfile_fn) {
      push @{$errors}, "subset_extractor error, see: $errfile_fn" ;
   } else {
      unlink $errfile_fn ;
   }

   {
      my ($newdev, $newino) = stat($localbdp) ;
      my ($olddev, $oldino) = stat($pdb_fn) ;
      if ($newdev != $olddev || $newino != $oldino) {
         unlink $localbdp ; }
   }

   return $errors ;
}


=head2 parse_modeller_ali()

   Title:       parse_modeller_ali()
   Function:    reads in a modeller PIR format alignment and returns residue
                number equivalence hashes
   Args:        $_->{ali_fn} alignment file
                $_->{modpipe_newstyle_orderswitch}
                  - 1 (Default) if new style MODPIPE run
                  - reordered sequences in the alignment file
   Results:     ->{seq} = $seq ;
                ->{resno_start} = $resno_start ;
                ->{resno_end} = $resno_end ;
                ->{chain_start} = $chain_start ;
                ->{chain_end} = $chain_end ;
                ->{alipos_2_serial} = $alipos_2_serresno ;
                ->{alipos_2_chainno} = $alipos_2_chainno ;
                ->{alipos_2_resna} = $alipos_2_resna ;
                ->{maxlength} = $maxlength;

=cut

sub parse_modeller_ali {

   my $params = shift ;
   my $ali_fn = $params->{ali_fn} ;

   if ((!defined $ali_fn) || (!-e $ali_fn)) {
      return {error => ['PDB file not specified or does not exist']} ; }

   open (ALIF, $ali_fn) ;
   my $headers ;
   my $seq ;
   my $cur_seq = -1 ;

   while ( my $line = <ALIF> ) {

      chomp $line;

      if (($line =~ /^structure/) || ($line =~ /^sequence/)) {

         if (exists $params->{modpipe_newstyle_orderswitch} &&
            $params->{modpipe_newstyle_orderswitch} == 1) {

            if ($line =~ /^structure/) { $cur_seq = 0 ; }
            elsif ($line =~ /^sequence/) { $cur_seq = 1 ; }

         } else {
	    $cur_seq++ ;
         }

	 $seq->[$cur_seq] = '' ;
         $headers->[$cur_seq] = $line ;

      } elsif ( ($line !~ /^$/) && ($line !~ /^\>/) &&
                ($line !~ /^C;/) && ($line !~ /^R;/) ) {

         $seq->[$cur_seq] .= $line ;

      }
   }
   close(ALIF) ;

   foreach my $j ( 0 .. $#{$seq}) {
      $seq->[$j] =~ s/\*$// ; }

   my ( $resno_start, $resno_end, $chain_start, $chain_end ) ;
   foreach my $j ( 0 .. $#{$headers} ) {
      my @t = split(':', $headers->[$j]) ;

      $t[2] =~ s/ //g ; $t[4] =~ s/ //g ;
      $t[2] =~ s/\.//g ; $t[4] =~ s/\.//g ;
      $t[3] =~ s/\.//g ; $t[5] =~ s/\.//g ;

      $resno_start->[$j] = $t[2] ;
      $chain_start->[$j] = $t[3] ;

#      print STDERR " from $headers->[$j]:\n" ;
#      print STDERR " modeller.pm: entry $j: STARTER IS $t[2] on $t[3]\n" ;
#      print STDERR " modeller.pm: entry $j: END     IS $t[4] on $t[5]\n" ;

      $resno_end->[$j] = $t[4] ;
      $chain_end->[$j] = $t[5] ;
   }

   my $alipos_2_serresno ;
   my $alipos_2_chainno ;
   my $alipos_2_resna ;
   my $maxlength = length($seq->[0]) ;

   foreach my $j ( 0 .. $#{$seq} ) {
      my $curres = 1 ;
      my $curchain = 1 ;
      my $curlength = length($seq->[$j]) ;
      if ($curlength > $maxlength) {$maxlength = $curlength} ;
      foreach my $k ( 0 .. ($curlength - 1) ) {
         my $curchar = substr($seq->[$j], $k, 1) ;
	 if ( $curchar eq "\\") {
	    $curres = 1 ;
	    $curchain++ ;
	 } elsif ( $curchar ne '-' ) {
            $alipos_2_serresno->[$j]->[$k] = $curres;
	    $alipos_2_chainno->[$j]->[$k] = $curchain ;
	    $alipos_2_resna->[$j]->[$k] = $curchar;
            $curres++ ;
         }

      }
   }

   my $results ;
   $results->{seq} = $seq ;
   $results->{resno_start} = $resno_start ;
   $results->{resno_end} = $resno_end ;
   $results->{chain_start} = $chain_start ;
   $results->{chain_end} = $chain_end ;
   $results->{alipos_2_serial} = $alipos_2_serresno ;
   $results->{alipos_2_chainno} = $alipos_2_chainno ;
   $results->{alipos_2_resna} = $alipos_2_resna ;
   $results->{maxlength} = $maxlength;
   return $results ;

}


=head2 get_salign (modeller_bin, bdp_file)

   Title:       get_salign()
   Function:    Calls MODELLER.SALIGN to structurally align two pdb files
   Args:        $_->{pdb_fn_1} - name of pdb file 1
                $_->{pdb_fn_2} - name of pdb file 2
                $_->{modeller_bin} - name of MODELLER binary file

=cut

sub get_salign {

# dont takes picks directly - no facility in MODELLER to pick discontiguous
# segments in both structures

   my $params = shift ;

   my $modeller_bin = $params->{modeller_bin} ;
   my $fn ;
   $fn->{pdb}->[0] = $params->{pdb_fn_1} ;
   $fn->{pdb}->[1] = $params->{pdb_fn_2} ;

# Specify the temporary alignment TOP file, and the output alignment file.

   my ($ali_top_fn, $ali_top_fh) ;
   my $temp_fh ;

   my ($fh) ;
   ($fh->{top}, $fn->{top}) =
      tempfile( "resequiv.align.XXXXXXX", SUFFIX => ".top") ;

# Generate the actual TOP file.


   my @afd ;
   foreach my $j ( 0 .. 1) {
      $fn->{pdbn}->[$j] = $fn->{pdb}->[$j] ; $fn->{pdbn}->[$j] =~ s/^.*\///g ;
      if ($fn->{pdb}->[$j] =~ /\//) {
         $fn->{pdbb}->[$j] = $fn->{pdb}->[$j] ;
         $fn->{pdbb}->[$j] =~ s/\/[^\/]+$//g ;
         push @afd, $fn->{pdbb}->[$j] ;
      }
   }

   my $t = timestamp() ;
   print {$fh->{top}} "# modtie::pibase::salign()\n#".$t."\n\n";
   print {$fh->{top}} "SET OUTPUT_CONTROL = 1 1 1 1 1\n" ;
   if ($#afd >= 0) {
      my $atom_file_dir  = ''; $atom_file_dir = join(':', @afd);
      print {$fh->{top}} "SET ATOM_FILES_DIRECTORY = \'$atom_file_dir\'\n\n" ; }

   print {$fh->{top}} "SET RMS_CUTOFFS = 3.5 3.5 60 60 15 60 60 60 60 60 60\n";
   print {$fh->{top}} "SET GAP_PENALTIES_3D = 0 3\n";
   print {$fh->{top}} "SET GAP_GAP_SCORE = 0, GAP_RESIDUE_SCORE = 0\n";
   print {$fh->{top}} "SET RR_FILE = \'\$(LIB)/as1.sim.mat\'\n" ;


   print {$fh->{top}} "READ_MODEL FILE = \'$fn->{pdbn}->[0]\'\n";
   print {$fh->{top}} "SEQUENCE_TO_ALI ALIGN_CODES = \'1_$fn->{pdbn}->[0]\',";
   print {$fh->{top}} " ATOM_FILES = \'$fn->{pdbn}->[0]\'\n" ;

   print {$fh->{top}} "READ_MODEL FILE = \'$fn->{pdbn}->[1]\'\n";
   print {$fh->{top}} "SEQUENCE_TO_ALI ADD_SEQUENCE = on, ";
   print {$fh->{top}} "ALIGN_CODES = ALIGN_CODES \'2_$fn->{pdbn}->[1]\', ";
   print {$fh->{top}} "ATOM_FILES = ATOM_FILES \'$fn->{pdbn}->[1]\'\n" ;

   print {$fh->{top}} "SET FEATURE_WEIGHTS = 1. 0. 0. 0. 1. 0., ";
   print {$fh->{top}} "GAP_PENALTIES_1D = -450 -50\n";
   print {$fh->{top}} "SALIGN OUTPUT = \'ALIGNMENT QUALITY\', IMPROVE_ALIGNMENT = on, FIT = on\n" ;


   print {$fh->{top}} "SET FEATURE_WEIGHTS = 1. .5 1. 1. 1. 0., ";
   print {$fh->{top}} "GAP_PENALTIES_1D = -450 -50\n";
   print {$fh->{top}} "SALIGN OUTPUT = \'ALIGNMENT QUALITY\', IMPROVE_ALIGNMENT = on, FIT = on\n" ;


   print {$fh->{top}} "SET WRITE_FIT = off, WRITE_WHOLE_PDB = off\n";
   print {$fh->{top}} "SET FEATURE_WEIGHTS = 1. 1. 1. 1. 1. 0., GAP_PENALTIES_1D = -450. -50.\n";
   print {$fh->{top}} "SALIGN OUTPUT = \'ALIGNMENT QUALITY\', IMPROVE_ALIGNMENT = on, FIT = on\n";
   print {$fh->{top}} "COMPARE OUTPUT = \'SHORT\'\n" ;

   close ($fh->{top}) ;

# Specify the location of the MODELLER LOG file.

   $fn->{'log'} = $fn->{top} ; $fn->{'log'} =~ s/top$/log/ ;

# Run the TOP file through MODELLER.

   system("$modeller_bin $fn->{top} >/dev/null 2>&1") ;

   if (-s $fn->{'log'}) {
      return {ali_log => $fn->{'log'}, ali_top => $fn->{'top'}} ;
   } else {
      return {error => 'error'} ;
   }

}



=head2 todload_bdp_ids(@results_type)

   Name: todload_bdp_ids()
   Function:    Returns bdp_id and depending on results_type
                specified, its relation to bdp_path and pdb_id
                in a variety of forms.
      Analogous to load_bdp_ids() with tables-on-disk instead of DBI
   Return:      query results
   Args:        $_[0] = DBI dbh handle
                $_[1] = results type

=over

=item * path_2_bdp_id (hash) [default]

=item * bdp_id_2_path (hash)

=item * bdp_id_2_pdb_id (hash)

=item * bdp_id_2_raw_pdb (hash)

=item * pdb_id_2_bdp_id (hash)

=item * bdp_id (array)

=back

=cut

sub todload_bdp_ids {

   my @result_types = @_ ;

   my $tablespecs = table_spec("bdp_files") ;
   my ($bdp_id, $file_path, $file_base, $pdb_id, $raw_pdb)=
      rawselect_tod("SELECT bdp_id, file_path, file_base, pdb_id, raw_pdb ".
                    "FROM bdp_files") ;

   if ($#result_types < 0) {
      push @result_types, 'path_2_bdp_id' ; }

   my @results ;

   foreach my $j ( 0 .. $#result_types) {

      my $ans ;
      if ($result_types[$j] eq 'path_2_bdp_id') {
         foreach my $k ( 0 .. $#{$bdp_id}) {
            $ans->{$file_path->[$k]} = $bdp_id->[$k] ; }
      } elsif ($result_types[$j] eq 'bdp_id_2_path') {
         foreach my $k ( 0 .. $#{$bdp_id}) {
            $ans->{$bdp_id->[$k]} = $file_path->[$k] ; }
      } elsif ($result_types[$j] eq 'bdp_id_2_pdb_id') {
         foreach my $k ( 0 .. $#{$bdp_id}) {
            $ans->{$bdp_id->[$k]} = $pdb_id->[$k] ; }
      } elsif ($result_types[$j] eq 'bdp_id_2_raw_pdb') {
         foreach my $k ( 0 .. $#{$bdp_id}) {
            $ans->{$bdp_id->[$k]} = $raw_pdb->[$k] ; }
      } elsif ($result_types[$j] eq 'pdb_id_2_bdp_id') {
         foreach my $k ( 0 .. $#{$bdp_id}) {
            if ($raw_pdb->[$k] == 1) {
               $ans->{$pdb_id->[$k]} = $bdp_id->[$k] ; }}
      } elsif ($result_types[$j] eq 'bdp_id') {
         push @{$ans}, @{$bdp_id} ;
      }

      push @results, $ans ;
   }

   return @results ;

}



=head2 full_table_specs()

   Title:       full_table_specs()
   Args:        none
   Returns:     $->{table_name}->{prikey} = primary key field
                $->{table_name}->{field_name}->[i] = name of ith field
                $->{table_name}->{field_spec}->[i] = type of ith field

=cut

sub full_table_specs {

   my $tables ;
   $tables->{scop_cla}->{prikey} = "(scop_id,chain_id,start_resno,end_resno)" ;
   $tables->{scop_cla}->{field_name}->[0] = "scop_id" ;
   $tables->{scop_cla}->{field_spec}->[0] = "char(7) not null" ;
   $tables->{scop_cla}->{field_name}->[1] = "pdb_id" ;
   $tables->{scop_cla}->{field_spec}->[1] = "char(4) not null" ;
   $tables->{scop_cla}->{field_name}->[2] = "chain_id" ;
   $tables->{scop_cla}->{field_spec}->[2] = "char(1) not null" ;
   $tables->{scop_cla}->{field_name}->[3] = "start_resno" ;
   $tables->{scop_cla}->{field_spec}->[3] = "char(10) not null" ;
   $tables->{scop_cla}->{field_name}->[4] = "end_resno" ;
   $tables->{scop_cla}->{field_spec}->[4] = "char(10) not null" ;
   $tables->{scop_cla}->{field_name}->[5] = "class_id" ;
   $tables->{scop_cla}->{field_spec}->[5] = "char(1) not null" ;
   $tables->{scop_cla}->{field_name}->[6] = "fold_id" ;
   $tables->{scop_cla}->{field_spec}->[6] = "smallint unsigned not null" ;
   $tables->{scop_cla}->{field_name}->[7] = "superfam_id" ;
   $tables->{scop_cla}->{field_spec}->[7] = "smallint unsigned not null" ;
   $tables->{scop_cla}->{field_name}->[8] = "fam_id" ;
   $tables->{scop_cla}->{field_spec}->[8] = "smallint unsigned not null" ;
   $tables->{scop_cla}->{field_name}->[9] = "cl_id" ;
   $tables->{scop_cla}->{field_spec}->[9] = "mediumint unsigned not null" ;
   $tables->{scop_cla}->{field_name}->[10] = "cf_id" ;
   $tables->{scop_cla}->{field_spec}->[10] = "mediumint unsigned not null" ;
   $tables->{scop_cla}->{field_name}->[11] = "sf_id" ;
   $tables->{scop_cla}->{field_spec}->[11] = "mediumint unsigned not null" ;
   $tables->{scop_cla}->{field_name}->[12] = "fa_id" ;
   $tables->{scop_cla}->{field_spec}->[12] = "mediumint unsigned not null" ;
   $tables->{scop_cla}->{field_name}->[13] = "dm_id" ;
   $tables->{scop_cla}->{field_spec}->[13] = "mediumint unsigned not null" ;
   $tables->{scop_cla}->{field_name}->[14] = "sp_id" ;
   $tables->{scop_cla}->{field_spec}->[14] = "mediumint unsigned not null" ;
   $tables->{scop_cla}->{field_name}->[15] = "px_id" ;
   $tables->{scop_cla}->{field_spec}->[15] = "mediumint unsigned not null" ;
   $tables->{interface_contacts_prototype}->{prikey} = "(bdp_id,subset_id_1,subset_id_2,chain_no_1,resno_1,chain_no_2,resno_2)" ;
   $tables->{interface_contacts_prototype}->{field_name}->[0] = "bdp_id" ;
   $tables->{interface_contacts_prototype}->{field_spec}->[0] = "integer unsigned not null" ;
   $tables->{interface_contacts_prototype}->{field_name}->[1] = "subset_id_1" ;
   $tables->{interface_contacts_prototype}->{field_spec}->[1] = "char(50) not null" ;
   $tables->{interface_contacts_prototype}->{field_name}->[2] = "subset_id_2" ;
   $tables->{interface_contacts_prototype}->{field_spec}->[2] = "char(50) not null" ;
   $tables->{interface_contacts_prototype}->{field_name}->[3] = "chain_no_1" ;
   $tables->{interface_contacts_prototype}->{field_spec}->[3] = "integer unsigned not null" ;
   $tables->{interface_contacts_prototype}->{field_name}->[4] = "chain_id_1" ;
   $tables->{interface_contacts_prototype}->{field_spec}->[4] = "char(1) not null" ;
   $tables->{interface_contacts_prototype}->{field_name}->[5] = "resno_1" ;
   $tables->{interface_contacts_prototype}->{field_spec}->[5] = "char(10) not null" ;
   $tables->{interface_contacts_prototype}->{field_name}->[6] = "resna_1" ;
   $tables->{interface_contacts_prototype}->{field_spec}->[6] = "char(3) not null" ;
   $tables->{interface_contacts_prototype}->{field_name}->[7] = "chain_no_2" ;
   $tables->{interface_contacts_prototype}->{field_spec}->[7] = "integer unsigned not null" ;
   $tables->{interface_contacts_prototype}->{field_name}->[8] = "chain_id_2" ;
   $tables->{interface_contacts_prototype}->{field_spec}->[8] = "char(1) not null" ;
   $tables->{interface_contacts_prototype}->{field_name}->[9] = "resno_2" ;
   $tables->{interface_contacts_prototype}->{field_spec}->[9] = "char(10) not null" ;
   $tables->{interface_contacts_prototype}->{field_name}->[10] = "resna_2" ;
   $tables->{interface_contacts_prototype}->{field_spec}->[10] = "char(3) not null" ;
   $tables->{interface_contacts_prototype}->{field_name}->[11] = "min_dist" ;
   $tables->{interface_contacts_prototype}->{field_spec}->[11] = "float not null" ;
   $tables->{interface_contacts_prototype}->{field_name}->[12] = "num_contacts" ;
   $tables->{interface_contacts_prototype}->{field_spec}->[12] = "integer unsigned not null" ;
   $tables->{interface_contacts_prototype}->{field_name}->[13] = "num_contacts_4" ;
   $tables->{interface_contacts_prototype}->{field_spec}->[13] = "integer unsigned not null" ;
   $tables->{interface_contacts_prototype}->{field_name}->[14] = "num_contacts_4p5" ;
   $tables->{interface_contacts_prototype}->{field_spec}->[14] = "integer unsigned not null" ;
   $tables->{interface_contacts_prototype}->{field_name}->[15] = "num_contacts_5" ;
   $tables->{interface_contacts_prototype}->{field_spec}->[15] = "integer unsigned not null" ;
   $tables->{interface_contacts_prototype}->{field_name}->[16] = "num_contacts_5p5" ;
   $tables->{interface_contacts_prototype}->{field_spec}->[16] = "integer unsigned not null" ;
   $tables->{pqs_ranking}->{prikey} = "(pqs_id)" ;
   $tables->{pqs_ranking}->{field_name}->[0] = "pqs_id" ;
   $tables->{pqs_ranking}->{field_spec}->[0] = "char(10) not null" ;
   $tables->{pqs_ranking}->{field_name}->[1] = "pdb_id" ;
   $tables->{pqs_ranking}->{field_spec}->[1] = "char(4) not null" ;
   $tables->{pqs_ranking}->{field_name}->[2] = "split_no" ;
   $tables->{pqs_ranking}->{field_spec}->[2] = "smallint unsigned" ;
   $tables->{pqs_ranking}->{field_name}->[3] = "ramach_most_fav" ;
   $tables->{pqs_ranking}->{field_spec}->[3] = "float" ;
   $tables->{pqs_ranking}->{field_name}->[4] = "morris_g" ;
   $tables->{pqs_ranking}->{field_spec}->[4] = "float" ;
   $tables->{pqs_ranking}->{field_name}->[5] = "b_main" ;
   $tables->{pqs_ranking}->{field_spec}->[5] = "float" ;
   $tables->{pqs_ranking}->{field_name}->[6] = "b_side" ;
   $tables->{pqs_ranking}->{field_spec}->[6] = "float" ;
   $tables->{pqs_ranking}->{field_name}->[7] = "no_atoms" ;
   $tables->{pqs_ranking}->{field_spec}->[7] = "mediumint unsigned" ;
   $tables->{pqs_ranking}->{field_name}->[8] = "no_hoh" ;
   $tables->{pqs_ranking}->{field_spec}->[8] = "smallint unsigned" ;
   $tables->{pqs_ranking}->{field_name}->[9] = "no_hetatm" ;
   $tables->{pqs_ranking}->{field_spec}->[9] = "smallint unsigned" ;
   $tables->{pdb_entry_type}->{prikey} = "(pdb_id)" ;
   $tables->{pdb_entry_type}->{field_name}->[0] = "pdb_id" ;
   $tables->{pdb_entry_type}->{field_spec}->[0] = "char(4) not null" ;
   $tables->{pdb_entry_type}->{field_name}->[1] = "entry_type" ;
   $tables->{pdb_entry_type}->{field_spec}->[1] = "char(20) not null" ;
   $tables->{pdb_entry_type}->{field_name}->[2] = "experiment_type" ;
   $tables->{pdb_entry_type}->{field_spec}->[2] = "char(11) not null" ;
   $tables->{cath_domain_list}->{prikey} = "(domain_name)" ;
   $tables->{cath_domain_list}->{field_name}->[0] = "domain_name" ;
   $tables->{cath_domain_list}->{field_spec}->[0] = "char(7) binary not null" ;
   $tables->{cath_domain_list}->{field_name}->[1] = "pdb_id" ;
   $tables->{cath_domain_list}->{field_spec}->[1] = "char(4) not null" ;
   $tables->{cath_domain_list}->{field_name}->[2] = "chain_id" ;
   $tables->{cath_domain_list}->{field_spec}->[2] = "char(1) binary not null" ;
   $tables->{cath_domain_list}->{field_name}->[3] = "domain_no" ;
   $tables->{cath_domain_list}->{field_spec}->[3] = "tinyint unsigned not null" ;
   $tables->{cath_domain_list}->{field_name}->[4] = "class" ;
   $tables->{cath_domain_list}->{field_spec}->[4] = "smallint unsigned not null" ;
   $tables->{cath_domain_list}->{field_name}->[5] = "arch" ;
   $tables->{cath_domain_list}->{field_spec}->[5] = "smallint unsigned not null" ;
   $tables->{cath_domain_list}->{field_name}->[6] = "topol" ;
   $tables->{cath_domain_list}->{field_spec}->[6] = "smallint unsigned not null" ;
   $tables->{cath_domain_list}->{field_name}->[7] = "homol" ;
   $tables->{cath_domain_list}->{field_spec}->[7] = "smallint unsigned not null" ;
   $tables->{cath_domain_list}->{field_name}->[8] = "s35no" ;
   $tables->{cath_domain_list}->{field_spec}->[8] = "smallint unsigned not null" ;
   $tables->{cath_domain_list}->{field_name}->[9] = "s95no" ;
   $tables->{cath_domain_list}->{field_spec}->[9] = "smallint unsigned not null" ;
   $tables->{cath_domain_list}->{field_name}->[10] = "s100no" ;
   $tables->{cath_domain_list}->{field_spec}->[10] = "smallint unsigned not null" ;
   $tables->{cath_domain_list}->{field_name}->[11] = "domain_length" ;
   $tables->{cath_domain_list}->{field_spec}->[11] = "smallint unsigned not null" ;
   $tables->{interface_secstrx_basic_contacts_tables}->{prikey} = "(bdp_id)" ;
   $tables->{interface_secstrx_basic_contacts_tables}->{field_name}->[0] = "bdp_id" ;
   $tables->{interface_secstrx_basic_contacts_tables}->{field_spec}->[0] = "integer unsigned not null" ;
   $tables->{interface_secstrx_basic_contacts_tables}->{field_name}->[1] = "table_name" ;
   $tables->{interface_secstrx_basic_contacts_tables}->{field_spec}->[1] = "char(40)" ;
   $tables->{interface_secstrx_basic_contacts_tables}->{field_name}->[2] = "source_file" ;
   $tables->{interface_secstrx_basic_contacts_tables}->{field_spec}->[2] = "char(255)" ;
   $tables->{pdb_release}->{prikey} = "(pdb_id)" ;
   $tables->{pdb_release}->{field_name}->[0] = "pdb_id" ;
   $tables->{pdb_release}->{field_spec}->[0] = "char(4) not null" ;
   $tables->{pdb_release}->{field_name}->[1] = "release_date" ;
   $tables->{pdb_release}->{field_spec}->[1] = "date" ;
   $tables->{external_data_sources}->{prikey} = "(data_source,table_name)" ;
   $tables->{external_data_sources}->{field_name}->[0] = "data_source" ;
   $tables->{external_data_sources}->{field_spec}->[0] = "char(50) not null" ;
   $tables->{external_data_sources}->{field_name}->[1] = "table_name" ;
   $tables->{external_data_sources}->{field_spec}->[1] = "char(100) not null" ;
   $tables->{external_data_sources}->{field_name}->[2] = "upload_date" ;
   $tables->{external_data_sources}->{field_spec}->[2] = "date not null" ;
   $tables->{external_data_sources}->{field_name}->[3] = "data_date" ;
   $tables->{external_data_sources}->{field_spec}->[3] = "date not null" ;
   $tables->{external_data_sources}->{field_name}->[4] = "version" ;
   $tables->{external_data_sources}->{field_spec}->[4] = "char(50)" ;
   $tables->{external_data_sources}->{field_name}->[5] = "url" ;
   $tables->{external_data_sources}->{field_spec}->[5] = "char(100)" ;
   $tables->{external_data_sources}->{field_name}->[6] = "comments" ;
   $tables->{external_data_sources}->{field_spec}->[6] = "char(100)" ;
   $tables->{bdp_numdom}->{prikey} = "(bdp_id,subset_source_id)" ;
   $tables->{bdp_numdom}->{field_name}->[0] = "bdp_id" ;
   $tables->{bdp_numdom}->{field_spec}->[0] = "integer unsigned not null" ;
   $tables->{bdp_numdom}->{field_name}->[1] = "subset_source_id" ;
   $tables->{bdp_numdom}->{field_spec}->[1] = "integer not null" ;
   $tables->{bdp_numdom}->{field_name}->[2] = "num_dom" ;
   $tables->{bdp_numdom}->{field_spec}->[2] = "integer not null" ;
   $tables->{bdp_files}->{prikey} = "(bdp_id)" ;
   $tables->{bdp_files}->{field_name}->[0] = "bdp_id" ;
   $tables->{bdp_files}->{field_spec}->[0] = "integer unsigned auto_increment" ;
   $tables->{bdp_files}->{field_name}->[1] = "file_path" ;
   $tables->{bdp_files}->{field_spec}->[1] = "char(255) not null" ;
   $tables->{bdp_files}->{field_name}->[2] = "file_base" ;
   $tables->{bdp_files}->{field_spec}->[2] = "char(30) not null" ;
   $tables->{bdp_files}->{field_name}->[3] = "pdb_id" ;
   $tables->{bdp_files}->{field_spec}->[3] = "char(4) not null" ;
   $tables->{bdp_files}->{field_name}->[4] = "raw_pdb" ;
   $tables->{bdp_files}->{field_spec}->[4] = "bool" ;
   $tables->{subsets}->{prikey} = "(subset_id)" ;
   $tables->{subsets}->{field_name}->[0] = "subset_id" ;
   $tables->{subsets}->{field_spec}->[0] = "char(50) binary not null" ;
   $tables->{subsets}->{field_name}->[1] = "bdp_id" ;
   $tables->{subsets}->{field_spec}->[1] = "integer" ;
   $tables->{subsets}->{field_name}->[2] = "pdb_id" ;
   $tables->{subsets}->{field_spec}->[2] = "char(4)" ;
   $tables->{subsets}->{field_name}->[3] = "description" ;
   $tables->{subsets}->{field_spec}->[3] = "char(250)" ;
   $tables->{subsets}->{field_name}->[4] = "subset_source_id" ;
   $tables->{subsets}->{field_spec}->[4] = "integer not null" ;
   $tables->{subsets}->{field_name}->[5] = "class" ;
   $tables->{subsets}->{field_spec}->[5] = "char(70) not null" ;
   $tables->{pdb_obsolete}->{prikey} = "(old_pdb_id)" ;
   $tables->{pdb_obsolete}->{field_name}->[0] = "old_pdb_id" ;
   $tables->{pdb_obsolete}->{field_spec}->[0] = "char(4) not null" ;
   $tables->{pdb_obsolete}->{field_name}->[1] = "new_pdb_id" ;
   $tables->{pdb_obsolete}->{field_spec}->[1] = "char(4)" ;
   $tables->{pdb_obsolete}->{field_name}->[2] = "obsolete_date" ;
   $tables->{pdb_obsolete}->{field_spec}->[2] = "date" ;
   $tables->{interface_secstrx_contacts_tables}->{prikey} = "(bdp_id)" ;
   $tables->{interface_secstrx_contacts_tables}->{field_name}->[0] = "bdp_id" ;
   $tables->{interface_secstrx_contacts_tables}->{field_spec}->[0] = "integer unsigned not null" ;
   $tables->{interface_secstrx_contacts_tables}->{field_name}->[1] = "table_name" ;
   $tables->{interface_secstrx_contacts_tables}->{field_spec}->[1] = "char(40)" ;
   $tables->{interface_secstrx_contacts_tables}->{field_name}->[2] = "source_file" ;
   $tables->{interface_secstrx_contacts_tables}->{field_spec}->[2] = "char(255)" ;
   $tables->{bdp_secstrx_tables}->{prikey} = "(bdp_id)" ;
   $tables->{bdp_secstrx_tables}->{field_name}->[0] = "bdp_id" ;
   $tables->{bdp_secstrx_tables}->{field_spec}->[0] = "integer unsigned not null" ;
   $tables->{bdp_secstrx_tables}->{field_name}->[1] = "table_name" ;
   $tables->{bdp_secstrx_tables}->{field_spec}->[1] = "char(40)" ;
   $tables->{bdp_secstrx_tables}->{field_name}->[2] = "source_file" ;
   $tables->{bdp_secstrx_tables}->{field_spec}->[2] = "char(255)" ;
   $tables->{subsets_files}->{prikey} = "(subset_id)" ;
   $tables->{subsets_files}->{field_name}->[0] = "subset_id" ;
   $tables->{subsets_files}->{field_spec}->[0] = "char(50) binary" ;
   $tables->{subsets_files}->{field_name}->[1] = "file_path" ;
   $tables->{subsets_files}->{field_spec}->[1] = "char(200) not null" ;
   $tables->{bdp_xpack}->{prikey} = "(bdp_id)" ;
   $tables->{bdp_xpack}->{field_name}->[0] = "bdp_id" ;
   $tables->{bdp_xpack}->{field_spec}->[0] = "integer unsigned not null" ;
   $tables->{bdp_xpack}->{field_name}->[1] = "xpack" ;
   $tables->{bdp_xpack}->{field_spec}->[1] = "tinyint(1) not null" ;
   $tables->{bindingsite_contacts_tables}->{prikey} = "(bdp_id)" ;
   $tables->{bindingsite_contacts_tables}->{field_name}->[0] = "bdp_id" ;
   $tables->{bindingsite_contacts_tables}->{field_spec}->[0] = "integer unsigned not null" ;
   $tables->{bindingsite_contacts_tables}->{field_name}->[1] = "cutoff" ;
   $tables->{bindingsite_contacts_tables}->{field_spec}->[1] = "float not null" ;
   $tables->{bindingsite_contacts_tables}->{field_name}->[2] = "table_name" ;
   $tables->{bindingsite_contacts_tables}->{field_spec}->[2] = "char(40)" ;
   $tables->{bindingsite_contacts_tables}->{field_name}->[3] = "source_file" ;
   $tables->{bindingsite_contacts_tables}->{field_spec}->[3] = "char(255)" ;
   $tables->{bdp_residues_tables}->{prikey} = "(bdp_id)" ;
   $tables->{bdp_residues_tables}->{field_name}->[0] = "bdp_id" ;
   $tables->{bdp_residues_tables}->{field_spec}->[0] = "integer unsigned not null" ;
   $tables->{bdp_residues_tables}->{field_name}->[1] = "table_name" ;
   $tables->{bdp_residues_tables}->{field_spec}->[1] = "char(40)" ;
   $tables->{bdp_residues_tables}->{field_name}->[2] = "source_file" ;
   $tables->{bdp_residues_tables}->{field_spec}->[2] = "char(255)" ;
   $tables->{bindingsite_secstrx_basic_contacts_prototype}->{prikey} = "(bdp_id,subset_id_1,subset_id_2,subset_id)" ;
   $tables->{bindingsite_secstrx_basic_contacts_prototype}->{field_name}->[0] = "bdp_id" ;
   $tables->{bindingsite_secstrx_basic_contacts_prototype}->{field_spec}->[0] = "integer unsigned not null" ;
   $tables->{bindingsite_secstrx_basic_contacts_prototype}->{field_name}->[1] = "subset_id_1" ;
   $tables->{bindingsite_secstrx_basic_contacts_prototype}->{field_spec}->[1] = "char(50) not null" ;
   $tables->{bindingsite_secstrx_basic_contacts_prototype}->{field_name}->[2] = "subset_id_2" ;
   $tables->{bindingsite_secstrx_basic_contacts_prototype}->{field_spec}->[2] = "char(50) not null" ;
   $tables->{bindingsite_secstrx_basic_contacts_prototype}->{field_name}->[3] = "subset_id" ;
   $tables->{bindingsite_secstrx_basic_contacts_prototype}->{field_spec}->[3] = "char(50) not null" ;
   $tables->{bindingsite_secstrx_basic_contacts_prototype}->{field_name}->[4] = "sse_id_1" ;
   $tables->{bindingsite_secstrx_basic_contacts_prototype}->{field_spec}->[4] = "integer unsigned not null" ;
   $tables->{bindingsite_secstrx_basic_contacts_prototype}->{field_name}->[5] = "sse_id_2" ;
   $tables->{bindingsite_secstrx_basic_contacts_prototype}->{field_spec}->[5] = "integer unsigned not null" ;
   $tables->{bindingsite_secstrx_basic_contacts_prototype}->{field_name}->[6] = "num_res_1" ;
   $tables->{bindingsite_secstrx_basic_contacts_prototype}->{field_spec}->[6] = "integer unsigned not null" ;
   $tables->{bindingsite_secstrx_basic_contacts_prototype}->{field_name}->[7] = "num_res_2" ;
   $tables->{bindingsite_secstrx_basic_contacts_prototype}->{field_spec}->[7] = "integer unsigned not null" ;
   $tables->{bindingsite_secstrx_basic_contacts_prototype}->{field_name}->[8] = "num_res_12" ;
   $tables->{bindingsite_secstrx_basic_contacts_prototype}->{field_spec}->[8] = "integer unsigned not null" ;
   $tables->{bindingsite_secstrx_basic_contacts_prototype}->{field_name}->[9] = "sse_1" ;
   $tables->{bindingsite_secstrx_basic_contacts_prototype}->{field_spec}->[9] = "enum('H', 'B', 'T', ' ') not null" ;
   $tables->{bindingsite_secstrx_basic_contacts_prototype}->{field_name}->[10] = "sse_2" ;
   $tables->{bindingsite_secstrx_basic_contacts_prototype}->{field_spec}->[10] = "enum('H', 'B', 'T', ' ') not null" ;
   $tables->{bindingsite_secstrx_basic_contacts_tables}->{prikey} = "(bdp_id)" ;
   $tables->{bindingsite_secstrx_basic_contacts_tables}->{field_name}->[0] = "bdp_id" ;
   $tables->{bindingsite_secstrx_basic_contacts_tables}->{field_spec}->[0] = "integer unsigned not null" ;
   $tables->{bindingsite_secstrx_basic_contacts_tables}->{field_name}->[1] = "table_name" ;
   $tables->{bindingsite_secstrx_basic_contacts_tables}->{field_spec}->[1] = "char(40)" ;
   $tables->{bindingsite_secstrx_basic_contacts_tables}->{field_name}->[2] = "source_file" ;
   $tables->{bindingsite_secstrx_basic_contacts_tables}->{field_spec}->[2] = "char(255)" ;
   $tables->{interface_sasa}->{prikey} = "(subset_id_1,subset_id_2)" ;
   $tables->{interface_sasa}->{field_name}->[0] = "bdp_id" ;
   $tables->{interface_sasa}->{field_spec}->[0] = "integer unsigned not null" ;
   $tables->{interface_sasa}->{field_name}->[1] = "subset_id_1" ;
   $tables->{interface_sasa}->{field_spec}->[1] = "char(50) not null" ;
   $tables->{interface_sasa}->{field_name}->[2] = "subset_id_2" ;
   $tables->{interface_sasa}->{field_spec}->[2] = "char(50) not null" ;
   $tables->{interface_sasa}->{field_name}->[3] = "dsasa_all" ;
   $tables->{interface_sasa}->{field_spec}->[3] = "float not null" ;
   $tables->{interface_sasa}->{field_name}->[4] = "dsasa_sc" ;
   $tables->{interface_sasa}->{field_spec}->[4] = "float not null" ;
   $tables->{interface_sasa}->{field_name}->[5] = "dsasa_mc" ;
   $tables->{interface_sasa}->{field_spec}->[5] = "float not null" ;
   $tables->{interface_sasa}->{field_name}->[6] = "dsasa_polar" ;
   $tables->{interface_sasa}->{field_spec}->[6] = "float not null" ;
   $tables->{interface_sasa}->{field_name}->[7] = "dsasa_nonpolar" ;
   $tables->{interface_sasa}->{field_spec}->[7] = "float not null" ;
   $tables->{pqs_biolist}->{prikey} = "(pdb_id)" ;
   $tables->{pqs_biolist}->{field_name}->[0] = "pdb_id" ;
   $tables->{pqs_biolist}->{field_spec}->[0] = "char(4) not null" ;
   $tables->{pqs_biolist}->{field_name}->[1] = "action" ;
   $tables->{pqs_biolist}->{field_spec}->[1] = "char(30) not null" ;
   $tables->{pqs_biolist}->{field_name}->[2] = "reason" ;
   $tables->{pqs_biolist}->{field_spec}->[2] = "char(40)" ;
   $tables->{pqs_biolist}->{field_name}->[3] = "num_split" ;
   $tables->{pqs_biolist}->{field_spec}->[3] = "smallint unsigned" ;
   $tables->{pqs_biolist}->{field_name}->[4] = "quaternary_state" ;
   $tables->{pqs_biolist}->{field_spec}->[4] = "char(50)" ;
   $tables->{pqs_biolist}->{field_name}->[5] = "memo" ;
   $tables->{pqs_biolist}->{field_spec}->[5] = "char(50)" ;
   $tables->{patch_residues_tables}->{prikey} = "(bdp_id)" ;
   $tables->{patch_residues_tables}->{field_name}->[0] = "bdp_id" ;
   $tables->{patch_residues_tables}->{field_spec}->[0] = "integer unsigned not null" ;
   $tables->{patch_residues_tables}->{field_name}->[1] = "cutoff" ;
   $tables->{patch_residues_tables}->{field_spec}->[1] = "float not null" ;
   $tables->{patch_residues_tables}->{field_name}->[2] = "table_name" ;
   $tables->{patch_residues_tables}->{field_spec}->[2] = "char(40)" ;
   $tables->{patch_residues_tables}->{field_name}->[3] = "source_file" ;
   $tables->{patch_residues_tables}->{field_spec}->[3] = "char(255)" ;
   $tables->{interface_size}->{prikey} = "(subset_id_1,subset_id_2)" ;
   $tables->{interface_size}->{field_name}->[0] = "bdp_id" ;
   $tables->{interface_size}->{field_spec}->[0] = "integer unsigned not null" ;
   $tables->{interface_size}->{field_name}->[1] = "subset_id_1" ;
   $tables->{interface_size}->{field_spec}->[1] = "char(50) not null" ;
   $tables->{interface_size}->{field_name}->[2] = "subset_id_2" ;
   $tables->{interface_size}->{field_spec}->[2] = "char(50) not null" ;
   $tables->{interface_size}->{field_name}->[3] = "num_res_1" ;
   $tables->{interface_size}->{field_spec}->[3] = "integer unsigned not null" ;
   $tables->{interface_size}->{field_name}->[4] = "num_res_2" ;
   $tables->{interface_size}->{field_spec}->[4] = "integer unsigned not null" ;
   $tables->{interface_size}->{field_name}->[5] = "subset_size_1" ;
   $tables->{interface_size}->{field_spec}->[5] = "integer unsigned not null" ;
   $tables->{interface_size}->{field_name}->[6] = "subset_size_2" ;
   $tables->{interface_size}->{field_spec}->[6] = "integer unsigned not null" ;
   $tables->{interface_secstrx_profile}->{prikey} = "(subset_id_1,subset_id_2)" ;
   $tables->{interface_secstrx_profile}->{field_name}->[0] = "bdp_id" ;
   $tables->{interface_secstrx_profile}->{field_spec}->[0] = "integer unsigned not null" ;
   $tables->{interface_secstrx_profile}->{field_name}->[1] = "subset_id_1" ;
   $tables->{interface_secstrx_profile}->{field_spec}->[1] = "char(50) not null" ;
   $tables->{interface_secstrx_profile}->{field_name}->[2] = "subset_id_2" ;
   $tables->{interface_secstrx_profile}->{field_spec}->[2] = "char(50) not null" ;
   $tables->{interface_secstrx_profile}->{field_name}->[3] = "num_sse_1" ;
   $tables->{interface_secstrx_profile}->{field_spec}->[3] = "integer unsigned not null" ;
   $tables->{interface_secstrx_profile}->{field_name}->[4] = "num_sse_2" ;
   $tables->{interface_secstrx_profile}->{field_spec}->[4] = "integer unsigned not null" ;
   $tables->{interface_secstrx_profile}->{field_name}->[5] = "num_sse_12" ;
   $tables->{interface_secstrx_profile}->{field_spec}->[5] = "integer unsigned not null" ;
   $tables->{interface_secstrx_profile}->{field_name}->[6] = "sse_contact_vector" ;
   $tables->{interface_secstrx_profile}->{field_spec}->[6] = "char(36) not null" ;
   $tables->{interface_secstrx_profile}->{field_name}->[7] = "sse_vector_1" ;
   $tables->{interface_secstrx_profile}->{field_spec}->[7] = "char(8) not null" ;
   $tables->{interface_secstrx_profile}->{field_name}->[8] = "sse_vector_2" ;
   $tables->{interface_secstrx_profile}->{field_spec}->[8] = "char(8) not null" ;
   $tables->{interface_secstrx_profile}->{field_name}->[9] = "sse_contact_vector_norm" ;
   $tables->{interface_secstrx_profile}->{field_spec}->[9] = "integer unsigned not null" ;
   $tables->{interface_secstrx_profile}->{field_name}->[10] = "sse_vector_1_norm" ;
   $tables->{interface_secstrx_profile}->{field_spec}->[10] = "integer unsigned not null" ;
   $tables->{interface_secstrx_profile}->{field_name}->[11] = "sse_vector_2_norm" ;
   $tables->{interface_secstrx_profile}->{field_spec}->[11] = "integer unsigned not null" ;
   $tables->{interface_secstrx_profile}->{field_name}->[12] = "num_sse_1_H" ;
   $tables->{interface_secstrx_profile}->{field_spec}->[12] = "integer unsigned not null" ;
   $tables->{interface_secstrx_profile}->{field_name}->[13] = "num_sse_1_B" ;
   $tables->{interface_secstrx_profile}->{field_spec}->[13] = "integer unsigned not null" ;
   $tables->{interface_secstrx_profile}->{field_name}->[14] = "num_sse_1_E" ;
   $tables->{interface_secstrx_profile}->{field_spec}->[14] = "integer unsigned not null" ;
   $tables->{interface_secstrx_profile}->{field_name}->[15] = "num_sse_1_G" ;
   $tables->{interface_secstrx_profile}->{field_spec}->[15] = "integer unsigned not null" ;
   $tables->{interface_secstrx_profile}->{field_name}->[16] = "num_sse_1_I" ;
   $tables->{interface_secstrx_profile}->{field_spec}->[16] = "integer unsigned not null" ;
   $tables->{interface_secstrx_profile}->{field_name}->[17] = "num_sse_1_T" ;
   $tables->{interface_secstrx_profile}->{field_spec}->[17] = "integer unsigned not null" ;
   $tables->{interface_secstrx_profile}->{field_name}->[18] = "num_sse_1_S" ;
   $tables->{interface_secstrx_profile}->{field_spec}->[18] = "integer unsigned not null" ;
   $tables->{interface_secstrx_profile}->{field_name}->[19] = "num_sse_1_u" ;
   $tables->{interface_secstrx_profile}->{field_spec}->[19] = "integer unsigned not null" ;
   $tables->{interface_secstrx_profile}->{field_name}->[20] = "num_sse_2_H" ;
   $tables->{interface_secstrx_profile}->{field_spec}->[20] = "integer unsigned not null" ;
   $tables->{interface_secstrx_profile}->{field_name}->[21] = "num_sse_2_B" ;
   $tables->{interface_secstrx_profile}->{field_spec}->[21] = "integer unsigned not null" ;
   $tables->{interface_secstrx_profile}->{field_name}->[22] = "num_sse_2_E" ;
   $tables->{interface_secstrx_profile}->{field_spec}->[22] = "integer unsigned not null" ;
   $tables->{interface_secstrx_profile}->{field_name}->[23] = "num_sse_2_G" ;
   $tables->{interface_secstrx_profile}->{field_spec}->[23] = "integer unsigned not null" ;
   $tables->{interface_secstrx_profile}->{field_name}->[24] = "num_sse_2_I" ;
   $tables->{interface_secstrx_profile}->{field_spec}->[24] = "integer unsigned not null" ;
   $tables->{interface_secstrx_profile}->{field_name}->[25] = "num_sse_2_T" ;
   $tables->{interface_secstrx_profile}->{field_spec}->[25] = "integer unsigned not null" ;
   $tables->{interface_secstrx_profile}->{field_name}->[26] = "num_sse_2_S" ;
   $tables->{interface_secstrx_profile}->{field_spec}->[26] = "integer unsigned not null" ;
   $tables->{interface_secstrx_profile}->{field_name}->[27] = "num_sse_2_u" ;
   $tables->{interface_secstrx_profile}->{field_spec}->[27] = "integer unsigned not null" ;
   $tables->{interface_secstrx_profile}->{field_name}->[28] = "num_sse_12_HH" ;
   $tables->{interface_secstrx_profile}->{field_spec}->[28] = "integer unsigned not null" ;
   $tables->{interface_secstrx_profile}->{field_name}->[29] = "num_sse_12_HB" ;
   $tables->{interface_secstrx_profile}->{field_spec}->[29] = "integer unsigned not null" ;
   $tables->{interface_secstrx_profile}->{field_name}->[30] = "num_sse_12_HE" ;
   $tables->{interface_secstrx_profile}->{field_spec}->[30] = "integer unsigned not null" ;
   $tables->{interface_secstrx_profile}->{field_name}->[31] = "num_sse_12_HG" ;
   $tables->{interface_secstrx_profile}->{field_spec}->[31] = "integer unsigned not null" ;
   $tables->{interface_secstrx_profile}->{field_name}->[32] = "num_sse_12_HI" ;
   $tables->{interface_secstrx_profile}->{field_spec}->[32] = "integer unsigned not null" ;
   $tables->{interface_secstrx_profile}->{field_name}->[33] = "num_sse_12_HT" ;
   $tables->{interface_secstrx_profile}->{field_spec}->[33] = "integer unsigned not null" ;
   $tables->{interface_secstrx_profile}->{field_name}->[34] = "num_sse_12_HS" ;
   $tables->{interface_secstrx_profile}->{field_spec}->[34] = "integer unsigned not null" ;
   $tables->{interface_secstrx_profile}->{field_name}->[35] = "num_sse_12_Hu" ;
   $tables->{interface_secstrx_profile}->{field_spec}->[35] = "integer unsigned not null" ;
   $tables->{interface_secstrx_profile}->{field_name}->[36] = "num_sse_12_BB" ;
   $tables->{interface_secstrx_profile}->{field_spec}->[36] = "integer unsigned not null" ;
   $tables->{interface_secstrx_profile}->{field_name}->[37] = "num_sse_12_BE" ;
   $tables->{interface_secstrx_profile}->{field_spec}->[37] = "integer unsigned not null" ;
   $tables->{interface_secstrx_profile}->{field_name}->[38] = "num_sse_12_BG" ;
   $tables->{interface_secstrx_profile}->{field_spec}->[38] = "integer unsigned not null" ;
   $tables->{interface_secstrx_profile}->{field_name}->[39] = "num_sse_12_BI" ;
   $tables->{interface_secstrx_profile}->{field_spec}->[39] = "integer unsigned not null" ;
   $tables->{interface_secstrx_profile}->{field_name}->[40] = "num_sse_12_BT" ;
   $tables->{interface_secstrx_profile}->{field_spec}->[40] = "integer unsigned not null" ;
   $tables->{interface_secstrx_profile}->{field_name}->[41] = "num_sse_12_BS" ;
   $tables->{interface_secstrx_profile}->{field_spec}->[41] = "integer unsigned not null" ;
   $tables->{interface_secstrx_profile}->{field_name}->[42] = "num_sse_12_Bu" ;
   $tables->{interface_secstrx_profile}->{field_spec}->[42] = "integer unsigned not null" ;
   $tables->{interface_secstrx_profile}->{field_name}->[43] = "num_sse_12_EE" ;
   $tables->{interface_secstrx_profile}->{field_spec}->[43] = "integer unsigned not null" ;
   $tables->{interface_secstrx_profile}->{field_name}->[44] = "num_sse_12_EG" ;
   $tables->{interface_secstrx_profile}->{field_spec}->[44] = "integer unsigned not null" ;
   $tables->{interface_secstrx_profile}->{field_name}->[45] = "num_sse_12_EI" ;
   $tables->{interface_secstrx_profile}->{field_spec}->[45] = "integer unsigned not null" ;
   $tables->{interface_secstrx_profile}->{field_name}->[46] = "num_sse_12_ET" ;
   $tables->{interface_secstrx_profile}->{field_spec}->[46] = "integer unsigned not null" ;
   $tables->{interface_secstrx_profile}->{field_name}->[47] = "num_sse_12_ES" ;
   $tables->{interface_secstrx_profile}->{field_spec}->[47] = "integer unsigned not null" ;
   $tables->{interface_secstrx_profile}->{field_name}->[48] = "num_sse_12_Eu" ;
   $tables->{interface_secstrx_profile}->{field_spec}->[48] = "integer unsigned not null" ;
   $tables->{interface_secstrx_profile}->{field_name}->[49] = "num_sse_12_GG" ;
   $tables->{interface_secstrx_profile}->{field_spec}->[49] = "integer unsigned not null" ;
   $tables->{interface_secstrx_profile}->{field_name}->[50] = "num_sse_12_GI" ;
   $tables->{interface_secstrx_profile}->{field_spec}->[50] = "integer unsigned not null" ;
   $tables->{interface_secstrx_profile}->{field_name}->[51] = "num_sse_12_GT" ;
   $tables->{interface_secstrx_profile}->{field_spec}->[51] = "integer unsigned not null" ;
   $tables->{interface_secstrx_profile}->{field_name}->[52] = "num_sse_12_GS" ;
   $tables->{interface_secstrx_profile}->{field_spec}->[52] = "integer unsigned not null" ;
   $tables->{interface_secstrx_profile}->{field_name}->[53] = "num_sse_12_Gu" ;
   $tables->{interface_secstrx_profile}->{field_spec}->[53] = "integer unsigned not null" ;
   $tables->{interface_secstrx_profile}->{field_name}->[54] = "num_sse_12_II" ;
   $tables->{interface_secstrx_profile}->{field_spec}->[54] = "integer unsigned not null" ;
   $tables->{interface_secstrx_profile}->{field_name}->[55] = "num_sse_12_IT" ;
   $tables->{interface_secstrx_profile}->{field_spec}->[55] = "integer unsigned not null" ;
   $tables->{interface_secstrx_profile}->{field_name}->[56] = "num_sse_12_IS" ;
   $tables->{interface_secstrx_profile}->{field_spec}->[56] = "integer unsigned not null" ;
   $tables->{interface_secstrx_profile}->{field_name}->[57] = "num_sse_12_Iu" ;
   $tables->{interface_secstrx_profile}->{field_spec}->[57] = "integer unsigned not null" ;
   $tables->{interface_secstrx_profile}->{field_name}->[58] = "num_sse_12_TT" ;
   $tables->{interface_secstrx_profile}->{field_spec}->[58] = "integer unsigned not null" ;
   $tables->{interface_secstrx_profile}->{field_name}->[59] = "num_sse_12_TS" ;
   $tables->{interface_secstrx_profile}->{field_spec}->[59] = "integer unsigned not null" ;
   $tables->{interface_secstrx_profile}->{field_name}->[60] = "num_sse_12_Tu" ;
   $tables->{interface_secstrx_profile}->{field_spec}->[60] = "integer unsigned not null" ;
   $tables->{interface_secstrx_profile}->{field_name}->[61] = "num_sse_12_SS" ;
   $tables->{interface_secstrx_profile}->{field_spec}->[61] = "integer unsigned not null" ;
   $tables->{interface_secstrx_profile}->{field_name}->[62] = "num_sse_12_Su" ;
   $tables->{interface_secstrx_profile}->{field_spec}->[62] = "integer unsigned not null" ;
   $tables->{interface_secstrx_profile}->{field_name}->[63] = "num_sse_12_uu" ;
   $tables->{interface_secstrx_profile}->{field_spec}->[63] = "integer unsigned not null" ;
   $tables->{interface_secstrx_profile}->{field_name}->[64] = "numres_sse_1_H" ;
   $tables->{interface_secstrx_profile}->{field_spec}->[64] = "integer unsigned not null" ;
   $tables->{interface_secstrx_profile}->{field_name}->[65] = "numres_sse_1_B" ;
   $tables->{interface_secstrx_profile}->{field_spec}->[65] = "integer unsigned not null" ;
   $tables->{interface_secstrx_profile}->{field_name}->[66] = "numres_sse_1_E" ;
   $tables->{interface_secstrx_profile}->{field_spec}->[66] = "integer unsigned not null" ;
   $tables->{interface_secstrx_profile}->{field_name}->[67] = "numres_sse_1_G" ;
   $tables->{interface_secstrx_profile}->{field_spec}->[67] = "integer unsigned not null" ;
   $tables->{interface_secstrx_profile}->{field_name}->[68] = "numres_sse_1_I" ;
   $tables->{interface_secstrx_profile}->{field_spec}->[68] = "integer unsigned not null" ;
   $tables->{interface_secstrx_profile}->{field_name}->[69] = "numres_sse_1_T" ;
   $tables->{interface_secstrx_profile}->{field_spec}->[69] = "integer unsigned not null" ;
   $tables->{interface_secstrx_profile}->{field_name}->[70] = "numres_sse_1_S" ;
   $tables->{interface_secstrx_profile}->{field_spec}->[70] = "integer unsigned not null" ;
   $tables->{interface_secstrx_profile}->{field_name}->[71] = "numres_sse_1_u" ;
   $tables->{interface_secstrx_profile}->{field_spec}->[71] = "integer unsigned not null" ;
   $tables->{interface_secstrx_profile}->{field_name}->[72] = "numres_sse_2_H" ;
   $tables->{interface_secstrx_profile}->{field_spec}->[72] = "integer unsigned not null" ;
   $tables->{interface_secstrx_profile}->{field_name}->[73] = "numres_sse_2_B" ;
   $tables->{interface_secstrx_profile}->{field_spec}->[73] = "integer unsigned not null" ;
   $tables->{interface_secstrx_profile}->{field_name}->[74] = "numres_sse_2_E" ;
   $tables->{interface_secstrx_profile}->{field_spec}->[74] = "integer unsigned not null" ;
   $tables->{interface_secstrx_profile}->{field_name}->[75] = "numres_sse_2_G" ;
   $tables->{interface_secstrx_profile}->{field_spec}->[75] = "integer unsigned not null" ;
   $tables->{interface_secstrx_profile}->{field_name}->[76] = "numres_sse_2_I" ;
   $tables->{interface_secstrx_profile}->{field_spec}->[76] = "integer unsigned not null" ;
   $tables->{interface_secstrx_profile}->{field_name}->[77] = "numres_sse_2_T" ;
   $tables->{interface_secstrx_profile}->{field_spec}->[77] = "integer unsigned not null" ;
   $tables->{interface_secstrx_profile}->{field_name}->[78] = "numres_sse_2_S" ;
   $tables->{interface_secstrx_profile}->{field_spec}->[78] = "integer unsigned not null" ;
   $tables->{interface_secstrx_profile}->{field_name}->[79] = "numres_sse_2_u" ;
   $tables->{interface_secstrx_profile}->{field_spec}->[79] = "integer unsigned not null" ;
   $tables->{interface_secstrx_profile}->{field_name}->[80] = "numres_sse_12_HH" ;
   $tables->{interface_secstrx_profile}->{field_spec}->[80] = "integer unsigned not null" ;
   $tables->{interface_secstrx_profile}->{field_name}->[81] = "numres_sse_12_HB" ;
   $tables->{interface_secstrx_profile}->{field_spec}->[81] = "integer unsigned not null" ;
   $tables->{interface_secstrx_profile}->{field_name}->[82] = "numres_sse_12_HE" ;
   $tables->{interface_secstrx_profile}->{field_spec}->[82] = "integer unsigned not null" ;
   $tables->{interface_secstrx_profile}->{field_name}->[83] = "numres_sse_12_HG" ;
   $tables->{interface_secstrx_profile}->{field_spec}->[83] = "integer unsigned not null" ;
   $tables->{interface_secstrx_profile}->{field_name}->[84] = "numres_sse_12_HI" ;
   $tables->{interface_secstrx_profile}->{field_spec}->[84] = "integer unsigned not null" ;
   $tables->{interface_secstrx_profile}->{field_name}->[85] = "numres_sse_12_HT" ;
   $tables->{interface_secstrx_profile}->{field_spec}->[85] = "integer unsigned not null" ;
   $tables->{interface_secstrx_profile}->{field_name}->[86] = "numres_sse_12_HS" ;
   $tables->{interface_secstrx_profile}->{field_spec}->[86] = "integer unsigned not null" ;
   $tables->{interface_secstrx_profile}->{field_name}->[87] = "numres_sse_12_Hu" ;
   $tables->{interface_secstrx_profile}->{field_spec}->[87] = "integer unsigned not null" ;
   $tables->{interface_secstrx_profile}->{field_name}->[88] = "numres_sse_12_BB" ;
   $tables->{interface_secstrx_profile}->{field_spec}->[88] = "integer unsigned not null" ;
   $tables->{interface_secstrx_profile}->{field_name}->[89] = "numres_sse_12_BE" ;
   $tables->{interface_secstrx_profile}->{field_spec}->[89] = "integer unsigned not null" ;
   $tables->{interface_secstrx_profile}->{field_name}->[90] = "numres_sse_12_BG" ;
   $tables->{interface_secstrx_profile}->{field_spec}->[90] = "integer unsigned not null" ;
   $tables->{interface_secstrx_profile}->{field_name}->[91] = "numres_sse_12_BI" ;
   $tables->{interface_secstrx_profile}->{field_spec}->[91] = "integer unsigned not null" ;
   $tables->{interface_secstrx_profile}->{field_name}->[92] = "numres_sse_12_BT" ;
   $tables->{interface_secstrx_profile}->{field_spec}->[92] = "integer unsigned not null" ;
   $tables->{interface_secstrx_profile}->{field_name}->[93] = "numres_sse_12_BS" ;
   $tables->{interface_secstrx_profile}->{field_spec}->[93] = "integer unsigned not null" ;
   $tables->{interface_secstrx_profile}->{field_name}->[94] = "numres_sse_12_Bu" ;
   $tables->{interface_secstrx_profile}->{field_spec}->[94] = "integer unsigned not null" ;
   $tables->{interface_secstrx_profile}->{field_name}->[95] = "numres_sse_12_EE" ;
   $tables->{interface_secstrx_profile}->{field_spec}->[95] = "integer unsigned not null" ;
   $tables->{interface_secstrx_profile}->{field_name}->[96] = "numres_sse_12_EG" ;
   $tables->{interface_secstrx_profile}->{field_spec}->[96] = "integer unsigned not null" ;
   $tables->{interface_secstrx_profile}->{field_name}->[97] = "numres_sse_12_EI" ;
   $tables->{interface_secstrx_profile}->{field_spec}->[97] = "integer unsigned not null" ;
   $tables->{interface_secstrx_profile}->{field_name}->[98] = "numres_sse_12_ET" ;
   $tables->{interface_secstrx_profile}->{field_spec}->[98] = "integer unsigned not null" ;
   $tables->{interface_secstrx_profile}->{field_name}->[99] = "numres_sse_12_ES" ;
   $tables->{interface_secstrx_profile}->{field_spec}->[99] = "integer unsigned not null" ;
   $tables->{interface_secstrx_profile}->{field_name}->[100] = "numres_sse_12_Eu" ;
   $tables->{interface_secstrx_profile}->{field_spec}->[100] = "integer unsigned not null" ;
   $tables->{interface_secstrx_profile}->{field_name}->[101] = "numres_sse_12_GG" ;
   $tables->{interface_secstrx_profile}->{field_spec}->[101] = "integer unsigned not null" ;
   $tables->{interface_secstrx_profile}->{field_name}->[102] = "numres_sse_12_GI" ;
   $tables->{interface_secstrx_profile}->{field_spec}->[102] = "integer unsigned not null" ;
   $tables->{interface_secstrx_profile}->{field_name}->[103] = "numres_sse_12_GT" ;
   $tables->{interface_secstrx_profile}->{field_spec}->[103] = "integer unsigned not null" ;
   $tables->{interface_secstrx_profile}->{field_name}->[104] = "numres_sse_12_GS" ;
   $tables->{interface_secstrx_profile}->{field_spec}->[104] = "integer unsigned not null" ;
   $tables->{interface_secstrx_profile}->{field_name}->[105] = "numres_sse_12_Gu" ;
   $tables->{interface_secstrx_profile}->{field_spec}->[105] = "integer unsigned not null" ;
   $tables->{interface_secstrx_profile}->{field_name}->[106] = "numres_sse_12_II" ;
   $tables->{interface_secstrx_profile}->{field_spec}->[106] = "integer unsigned not null" ;
   $tables->{interface_secstrx_profile}->{field_name}->[107] = "numres_sse_12_IT" ;
   $tables->{interface_secstrx_profile}->{field_spec}->[107] = "integer unsigned not null" ;
   $tables->{interface_secstrx_profile}->{field_name}->[108] = "numres_sse_12_IS" ;
   $tables->{interface_secstrx_profile}->{field_spec}->[108] = "integer unsigned not null" ;
   $tables->{interface_secstrx_profile}->{field_name}->[109] = "numres_sse_12_Iu" ;
   $tables->{interface_secstrx_profile}->{field_spec}->[109] = "integer unsigned not null" ;
   $tables->{interface_secstrx_profile}->{field_name}->[110] = "numres_sse_12_TT" ;
   $tables->{interface_secstrx_profile}->{field_spec}->[110] = "integer unsigned not null" ;
   $tables->{interface_secstrx_profile}->{field_name}->[111] = "numres_sse_12_TS" ;
   $tables->{interface_secstrx_profile}->{field_spec}->[111] = "integer unsigned not null" ;
   $tables->{interface_secstrx_profile}->{field_name}->[112] = "numres_sse_12_Tu" ;
   $tables->{interface_secstrx_profile}->{field_spec}->[112] = "integer unsigned not null" ;
   $tables->{interface_secstrx_profile}->{field_name}->[113] = "numres_sse_12_SS" ;
   $tables->{interface_secstrx_profile}->{field_spec}->[113] = "integer unsigned not null" ;
   $tables->{interface_secstrx_profile}->{field_name}->[114] = "numres_sse_12_Su" ;
   $tables->{interface_secstrx_profile}->{field_spec}->[114] = "integer unsigned not null" ;
   $tables->{interface_secstrx_profile}->{field_name}->[115] = "numres_sse_12_uu" ;
   $tables->{interface_secstrx_profile}->{field_spec}->[115] = "integer unsigned not null" ;
   $tables->{scop_interface_clusters}->{prikey} = "(bdp_id,subset_id_1,subset_id_2)" ;
   $tables->{scop_interface_clusters}->{field_name}->[0] = "bdp_id" ;
   $tables->{scop_interface_clusters}->{field_spec}->[0] = "integer unsigned not null" ;
   $tables->{scop_interface_clusters}->{field_name}->[1] = "subset_id_1" ;
   $tables->{scop_interface_clusters}->{field_spec}->[1] = "char(50) not null" ;
   $tables->{scop_interface_clusters}->{field_name}->[2] = "subset_id_2" ;
   $tables->{scop_interface_clusters}->{field_spec}->[2] = "char(50) not null" ;
   $tables->{scop_interface_clusters}->{field_name}->[3] = "scopclass_pair" ;
   $tables->{scop_interface_clusters}->{field_spec}->[3] = "char(150) not null" ;
   $tables->{scop_interface_clusters}->{field_name}->[4] = "cluster_level" ;
   $tables->{scop_interface_clusters}->{field_spec}->[4] = "enum('fam','sf') not null" ;
   $tables->{scop_interface_clusters}->{field_name}->[5] = "interface_class" ;
   $tables->{scop_interface_clusters}->{field_spec}->[5] = "char(150) not null" ;
   $tables->{scop_interface_clusters}->{field_name}->[6] = "cluster_no" ;
   $tables->{scop_interface_clusters}->{field_spec}->[6] = "integer not null" ;
   $tables->{scop_interface_clusters}->{field_name}->[7] = "member_no" ;
   $tables->{scop_interface_clusters}->{field_spec}->[7] = "integer not null" ;
   $tables->{interface_contacts_tables}->{prikey} = "(bdp_id)" ;
   $tables->{interface_contacts_tables}->{field_name}->[0] = "bdp_id" ;
   $tables->{interface_contacts_tables}->{field_spec}->[0] = "integer unsigned not null" ;
   $tables->{interface_contacts_tables}->{field_name}->[1] = "cutoff" ;
   $tables->{interface_contacts_tables}->{field_spec}->[1] = "float not null" ;
   $tables->{interface_contacts_tables}->{field_name}->[2] = "table_name" ;
   $tables->{interface_contacts_tables}->{field_spec}->[2] = "char(40)" ;
   $tables->{interface_contacts_tables}->{field_name}->[3] = "source_file" ;
   $tables->{interface_contacts_tables}->{field_spec}->[3] = "char(255)" ;
   $tables->{bdp_chains}->{prikey} = "(bdp_id,real_chain_no)" ;
   $tables->{bdp_chains}->{field_name}->[0] = "bdp_id" ;
   $tables->{bdp_chains}->{field_spec}->[0] = "integer unsigned not null" ;
   $tables->{bdp_chains}->{field_name}->[1] = "real_chain_no" ;
   $tables->{bdp_chains}->{field_spec}->[1] = "integer unsigned not null" ;
   $tables->{bdp_chains}->{field_name}->[2] = "real_chain_id" ;
   $tables->{bdp_chains}->{field_spec}->[2] = "char(1) not null" ;
   $tables->{bdp_chains}->{field_name}->[3] = "chain_type" ;
   $tables->{bdp_chains}->{field_spec}->[3] = "enum('p', 'n') not null" ;
   $tables->{bdp_chains}->{field_name}->[4] = "pdb_chain_no" ;
   $tables->{bdp_chains}->{field_spec}->[4] = "integer unsigned" ;
   $tables->{bdp_chains}->{field_name}->[5] = "pdb_chain_id" ;
   $tables->{bdp_chains}->{field_spec}->[5] = "char(1)" ;
   $tables->{bdp_chains}->{field_name}->[6] = "start_resno" ;
   $tables->{bdp_chains}->{field_spec}->[6] = "char(10) not null" ;
   $tables->{bdp_chains}->{field_name}->[7] = "start_resno_int" ;
   $tables->{bdp_chains}->{field_spec}->[7] = "integer not null" ;
   $tables->{bdp_chains}->{field_name}->[8] = "end_resno" ;
   $tables->{bdp_chains}->{field_spec}->[8] = "char(10) not null" ;
   $tables->{bdp_chains}->{field_name}->[9] = "end_resno_int" ;
   $tables->{bdp_chains}->{field_spec}->[9] = "integer not null" ;
   $tables->{bdp_chains}->{field_name}->[10] = "num_res" ;
   $tables->{bdp_chains}->{field_spec}->[10] = "integer unsigned not null" ;
   $tables->{bdp_chains}->{field_name}->[11] = "num_atoms" ;
   $tables->{bdp_chains}->{field_spec}->[11] = "integer unsigned not null" ;
   $tables->{bdp_chains}->{field_name}->[12] = "num_hetatm" ;
   $tables->{bdp_chains}->{field_spec}->[12] = "integer unsigned not null" ;
   $tables->{bdp_chains}->{field_name}->[13] = "sequence" ;
   $tables->{bdp_chains}->{field_spec}->[13] = "text not null" ;
   $tables->{pqs_list}->{prikey} = "(pqs_file_base)" ;
   $tables->{pqs_list}->{field_name}->[0] = "pqs_file_base" ;
   $tables->{pqs_list}->{field_spec}->[0] = "char(50) not null" ;
   $tables->{pqs_list}->{field_name}->[1] = "pdb_id" ;
   $tables->{pqs_list}->{field_spec}->[1] = "char(4) not null" ;
   $tables->{subsets_sequence}->{prikey} = "(subset_id)" ;
   $tables->{subsets_sequence}->{field_name}->[0] = "bdp_id" ;
   $tables->{subsets_sequence}->{field_spec}->[0] = "integer not null" ;
   $tables->{subsets_sequence}->{field_name}->[1] = "subset_id" ;
   $tables->{subsets_sequence}->{field_spec}->[1] = "char(50) binary not null" ;
   $tables->{subsets_sequence}->{field_name}->[2] = "num_chains" ;
   $tables->{subsets_sequence}->{field_spec}->[2] = "integer not null" ;
   $tables->{subsets_sequence}->{field_name}->[3] = "num_res" ;
   $tables->{subsets_sequence}->{field_spec}->[3] = "integer not null" ;
   $tables->{subsets_sequence}->{field_name}->[4] = "sequence" ;
   $tables->{subsets_sequence}->{field_spec}->[4] = "text not null" ;
   $tables->{patch_residues_prototype}->{prikey} = "(bdp_id,subset_id,chain_no,resno)" ;
   $tables->{patch_residues_prototype}->{field_name}->[0] = "bdp_id" ;
   $tables->{patch_residues_prototype}->{field_spec}->[0] = "integer unsigned not null" ;
   $tables->{patch_residues_prototype}->{field_name}->[1] = "subset_id" ;
   $tables->{patch_residues_prototype}->{field_spec}->[1] = "char(50) not null" ;
   $tables->{patch_residues_prototype}->{field_name}->[2] = "chain_no" ;
   $tables->{patch_residues_prototype}->{field_spec}->[2] = "integer unsigned not null" ;
   $tables->{patch_residues_prototype}->{field_name}->[3] = "chain_id" ;
   $tables->{patch_residues_prototype}->{field_spec}->[3] = "char(1) not null" ;
   $tables->{patch_residues_prototype}->{field_name}->[4] = "resno" ;
   $tables->{patch_residues_prototype}->{field_spec}->[4] = "char(10) not null" ;
   $tables->{patch_residues_prototype}->{field_name}->[5] = "resna" ;
   $tables->{patch_residues_prototype}->{field_spec}->[5] = "char(3) not null" ;
   $tables->{patch_residues_prototype}->{field_name}->[6] = "num_contacts" ;
   $tables->{patch_residues_prototype}->{field_spec}->[6] = "integer unsigned" ;
   $tables->{bindingsite_sse_topology}->{prikey} = "(bdp_id,subset_id_1,subset_id_2,subset_id)" ;
   $tables->{bindingsite_sse_topology}->{field_name}->[0] = "bdp_id" ;
   $tables->{bindingsite_sse_topology}->{field_spec}->[0] = "integer unsigned not null" ;
   $tables->{bindingsite_sse_topology}->{field_name}->[1] = "subset_id_1" ;
   $tables->{bindingsite_sse_topology}->{field_spec}->[1] = "char(50) not null" ;
   $tables->{bindingsite_sse_topology}->{field_name}->[2] = "subset_id_2" ;
   $tables->{bindingsite_sse_topology}->{field_spec}->[2] = "char(50) not null" ;
   $tables->{bindingsite_sse_topology}->{field_name}->[3] = "subset_id" ;
   $tables->{bindingsite_sse_topology}->{field_spec}->[3] = "char(50) not null" ;
   $tables->{bindingsite_sse_topology}->{field_name}->[4] = "num_domains" ;
   $tables->{bindingsite_sse_topology}->{field_spec}->[4] = "integer unsigned not null" ;
   $tables->{bindingsite_sse_topology}->{field_name}->[5] = "num_edges" ;
   $tables->{bindingsite_sse_topology}->{field_spec}->[5] = "integer unsigned not null" ;
   $tables->{bindingsite_sse_topology}->{field_name}->[6] = "nodelist" ;
   $tables->{bindingsite_sse_topology}->{field_spec}->[6] = "text not null" ;
   $tables->{bindingsite_sse_topology}->{field_name}->[7] = "nodelist_sse" ;
   $tables->{bindingsite_sse_topology}->{field_spec}->[7] = "text not null" ;
   $tables->{bindingsite_sse_topology}->{field_name}->[8] = "edgelist" ;
   $tables->{bindingsite_sse_topology}->{field_spec}->[8] = "text not null" ;
   $tables->{bindingsite_sse_topology}->{field_name}->[9] = "edgelist_sse" ;
   $tables->{bindingsite_sse_topology}->{field_spec}->[9] = "text not null" ;
   $tables->{bdp_residues_prototype}->{prikey} = "(bdp_id,chain_no,resno_serial)" ;
   $tables->{bdp_residues_prototype}->{field_name}->[0] = "bdp_id" ;
   $tables->{bdp_residues_prototype}->{field_spec}->[0] = "integer unsigned not null" ;
   $tables->{bdp_residues_prototype}->{field_name}->[1] = "chain_no" ;
   $tables->{bdp_residues_prototype}->{field_spec}->[1] = "integer unsigned not null" ;
   $tables->{bdp_residues_prototype}->{field_name}->[2] = "chain_id" ;
   $tables->{bdp_residues_prototype}->{field_spec}->[2] = "char(1) not null" ;
   $tables->{bdp_residues_prototype}->{field_name}->[3] = "resno_serial" ;
   $tables->{bdp_residues_prototype}->{field_spec}->[3] = "integer unsigned not null" ;
   $tables->{bdp_residues_prototype}->{field_name}->[4] = "resno" ;
   $tables->{bdp_residues_prototype}->{field_spec}->[4] = "char(10) not null" ;
   $tables->{bdp_residues_prototype}->{field_name}->[5] = "resno_int" ;
   $tables->{bdp_residues_prototype}->{field_spec}->[5] = "integer not null" ;
   $tables->{bdp_residues_prototype}->{field_name}->[6] = "resna" ;
   $tables->{bdp_residues_prototype}->{field_spec}->[6] = "char(3) not null" ;
   $tables->{bdp_residues_prototype}->{field_name}->[7] = "chain_type" ;
   $tables->{bdp_residues_prototype}->{field_spec}->[7] = "enum('p', 'n') not null" ;
   $tables->{pdb_patchres_bychain}->{prikey} = "(pdb_id,chain_id,patch_no)" ;
   $tables->{pdb_patchres_bychain}->{field_name}->[0] = "pdb_id" ;
   $tables->{pdb_patchres_bychain}->{field_spec}->[0] = "char(4) not null" ;
   $tables->{pdb_patchres_bychain}->{field_name}->[1] = "chain_id" ;
   $tables->{pdb_patchres_bychain}->{field_spec}->[1] = "char(1) not null" ;
   $tables->{pdb_patchres_bychain}->{field_name}->[2] = "patch_no" ;
   $tables->{pdb_patchres_bychain}->{field_spec}->[2] = "integer not null" ;
   $tables->{pdb_patchres_bychain}->{field_name}->[3] = "resnos" ;
   $tables->{pdb_patchres_bychain}->{field_spec}->[3] = "text not null" ;
   $tables->{interface_continuity}->{prikey} = "(subset_id_1,subset_id_2)" ;
   $tables->{interface_continuity}->{field_name}->[0] = "bdp_id" ;
   $tables->{interface_continuity}->{field_spec}->[0] = "integer unsigned not null" ;
   $tables->{interface_continuity}->{field_name}->[1] = "subset_id_1" ;
   $tables->{interface_continuity}->{field_spec}->[1] = "char(50) not null" ;
   $tables->{interface_continuity}->{field_name}->[2] = "subset_id_2" ;
   $tables->{interface_continuity}->{field_spec}->[2] = "char(50) not null" ;
   $tables->{interface_continuity}->{field_name}->[3] = "seq_segments_1" ;
   $tables->{interface_continuity}->{field_spec}->[3] = "integer unsigned not null" ;
   $tables->{interface_continuity}->{field_name}->[4] = "strx_patches_1" ;
   $tables->{interface_continuity}->{field_spec}->[4] = "integer unsigned not null" ;
   $tables->{interface_continuity}->{field_name}->[5] = "seq_segments_2" ;
   $tables->{interface_continuity}->{field_spec}->[5] = "integer unsigned not null" ;
   $tables->{interface_continuity}->{field_name}->[6] = "strx_patches_2" ;
   $tables->{interface_continuity}->{field_spec}->[6] = "integer unsigned not null" ;
   $tables->{interface_continuity}->{field_name}->[7] = "cutoff" ;
   $tables->{interface_continuity}->{field_spec}->[7] = "float not null" ;
   $tables->{bdp_interaction_topology_graph}->{prikey} = "(bdp_id,subset_source_id)" ;
   $tables->{bdp_interaction_topology_graph}->{field_name}->[0] = "bdp_id" ;
   $tables->{bdp_interaction_topology_graph}->{field_spec}->[0] = "integer unsigned not null" ;
   $tables->{bdp_interaction_topology_graph}->{field_name}->[1] = "subset_source_id" ;
   $tables->{bdp_interaction_topology_graph}->{field_spec}->[1] = "integer unsigned not null" ;
   $tables->{bdp_interaction_topology_graph}->{field_name}->[2] = "file_path" ;
   $tables->{bdp_interaction_topology_graph}->{field_spec}->[2] = "char(200) not null" ;
   $tables->{interface_planarity}->{prikey} = "(subset_id_1,subset_id_2)" ;
   $tables->{interface_planarity}->{field_name}->[0] = "bdp_id" ;
   $tables->{interface_planarity}->{field_spec}->[0] = "integer unsigned not null" ;
   $tables->{interface_planarity}->{field_name}->[1] = "subset_id_1" ;
   $tables->{interface_planarity}->{field_spec}->[1] = "char(50) not null" ;
   $tables->{interface_planarity}->{field_name}->[2] = "subset_id_2" ;
   $tables->{interface_planarity}->{field_spec}->[2] = "char(50) not null" ;
   $tables->{interface_planarity}->{field_name}->[3] = "planarity" ;
   $tables->{interface_planarity}->{field_spec}->[3] = "float not null" ;
   $tables->{subsets_residues_tables}->{prikey} = "(bdp_id)" ;
   $tables->{subsets_residues_tables}->{field_name}->[0] = "bdp_id" ;
   $tables->{subsets_residues_tables}->{field_spec}->[0] = "integer unsigned not null" ;
   $tables->{subsets_residues_tables}->{field_name}->[1] = "table_name" ;
   $tables->{subsets_residues_tables}->{field_spec}->[1] = "char(40)" ;
   $tables->{subsets_residues_tables}->{field_name}->[2] = "source_file" ;
   $tables->{subsets_residues_tables}->{field_spec}->[2] = "char(255)" ;
   $tables->{interface_resvector}->{prikey} = "(subset_id_1,subset_id_2)" ;
   $tables->{interface_resvector}->{field_name}->[0] = "bdp_id" ;
   $tables->{interface_resvector}->{field_spec}->[0] = "integer unsigned not null" ;
   $tables->{interface_resvector}->{field_name}->[1] = "subset_id_1" ;
   $tables->{interface_resvector}->{field_spec}->[1] = "char(50) not null" ;
   $tables->{interface_resvector}->{field_name}->[2] = "subset_id_2" ;
   $tables->{interface_resvector}->{field_spec}->[2] = "char(50) not null" ;
   $tables->{interface_resvector}->{field_name}->[3] = "contact_vector" ;
   $tables->{interface_resvector}->{field_spec}->[3] = "char(210) not null" ;
   $tables->{interface_resvector}->{field_name}->[4] = "res_vector_1" ;
   $tables->{interface_resvector}->{field_spec}->[4] = "char(20) not null" ;
   $tables->{interface_resvector}->{field_name}->[5] = "res_vector_2" ;
   $tables->{interface_resvector}->{field_spec}->[5] = "char(20) not null" ;
   $tables->{interface_resvector}->{field_name}->[6] = "contact_vector_norm" ;
   $tables->{interface_resvector}->{field_spec}->[6] = "integer not null" ;
   $tables->{interface_resvector}->{field_name}->[7] = "res_vector_1_norm" ;
   $tables->{interface_resvector}->{field_spec}->[7] = "integer not null" ;
   $tables->{interface_resvector}->{field_name}->[8] = "res_vector_2_norm" ;
   $tables->{interface_resvector}->{field_spec}->[8] = "integer not null" ;
   $tables->{cath_domain_description}->{prikey} = "(domain_name,segment_id)" ;
   $tables->{cath_domain_description}->{field_name}->[0] = "domain_name" ;
   $tables->{cath_domain_description}->{field_spec}->[0] = "char(7) binary not null" ;
   $tables->{cath_domain_description}->{field_name}->[1] = "pdb_id" ;
   $tables->{cath_domain_description}->{field_spec}->[1] = "char(4) not null" ;
   $tables->{cath_domain_description}->{field_name}->[2] = "chain_id" ;
   $tables->{cath_domain_description}->{field_spec}->[2] = "char(1) binary not null" ;
   $tables->{cath_domain_description}->{field_name}->[3] = "domain_no" ;
   $tables->{cath_domain_description}->{field_spec}->[3] = "tinyint unsigned not null" ;
   $tables->{cath_domain_description}->{field_name}->[4] = "class_id" ;
   $tables->{cath_domain_description}->{field_spec}->[4] = "smallint unsigned not null" ;
   $tables->{cath_domain_description}->{field_name}->[5] = "arch_id" ;
   $tables->{cath_domain_description}->{field_spec}->[5] = "smallint unsigned not null" ;
   $tables->{cath_domain_description}->{field_name}->[6] = "topol_id" ;
   $tables->{cath_domain_description}->{field_spec}->[6] = "smallint unsigned not null" ;
   $tables->{cath_domain_description}->{field_name}->[7] = "homol_id" ;
   $tables->{cath_domain_description}->{field_spec}->[7] = "smallint unsigned not null" ;
   $tables->{cath_domain_description}->{field_name}->[8] = "class" ;
   $tables->{cath_domain_description}->{field_spec}->[8] = "char(30) not null" ;
   $tables->{cath_domain_description}->{field_name}->[9] = "arch" ;
   $tables->{cath_domain_description}->{field_spec}->[9] = "char(50) not null" ;
   $tables->{cath_domain_description}->{field_name}->[10] = "topol" ;
   $tables->{cath_domain_description}->{field_spec}->[10] = "char(150) not null" ;
   $tables->{cath_domain_description}->{field_name}->[11] = "homol" ;
   $tables->{cath_domain_description}->{field_spec}->[11] = "char(150) not null" ;
   $tables->{cath_domain_description}->{field_name}->[12] = "domain_length" ;
   $tables->{cath_domain_description}->{field_spec}->[12] = "smallint unsigned not null" ;
   $tables->{cath_domain_description}->{field_name}->[13] = "domain_sequence_header" ;
   $tables->{cath_domain_description}->{field_spec}->[13] = "char(25) not null" ;
   $tables->{cath_domain_description}->{field_name}->[14] = "domain_sequence" ;
   $tables->{cath_domain_description}->{field_spec}->[14] = "text not null" ;
   $tables->{cath_domain_description}->{field_name}->[15] = "segment_id" ;
   $tables->{cath_domain_description}->{field_spec}->[15] = "tinyint not null" ;
   $tables->{cath_domain_description}->{field_name}->[16] = "start_resno" ;
   $tables->{cath_domain_description}->{field_spec}->[16] = "char(10) not null" ;
   $tables->{cath_domain_description}->{field_name}->[17] = "end_resno" ;
   $tables->{cath_domain_description}->{field_spec}->[17] = "char(10) not null" ;
   $tables->{cath_domain_description}->{field_name}->[18] = "segment_length" ;
   $tables->{cath_domain_description}->{field_spec}->[18] = "smallint unsigned not null" ;
   $tables->{cath_domain_description}->{field_name}->[19] = "segment_sequence_header" ;
   $tables->{cath_domain_description}->{field_spec}->[19] = "char(25) not null" ;
   $tables->{cath_domain_description}->{field_name}->[20] = "segment_sequence" ;
   $tables->{cath_domain_description}->{field_spec}->[20] = "text not null" ;
   $tables->{interface_secstrx_contacts_prototype}->{prikey} = "(bdp_id,subset_id_1,subset_id_2)" ;
   $tables->{interface_secstrx_contacts_prototype}->{field_name}->[0] = "bdp_id" ;
   $tables->{interface_secstrx_contacts_prototype}->{field_spec}->[0] = "integer unsigned not null" ;
   $tables->{interface_secstrx_contacts_prototype}->{field_name}->[1] = "subset_id_1" ;
   $tables->{interface_secstrx_contacts_prototype}->{field_spec}->[1] = "char(50) not null" ;
   $tables->{interface_secstrx_contacts_prototype}->{field_name}->[2] = "subset_id_2" ;
   $tables->{interface_secstrx_contacts_prototype}->{field_spec}->[2] = "char(50) not null" ;
   $tables->{interface_secstrx_contacts_prototype}->{field_name}->[3] = "sse_id_1" ;
   $tables->{interface_secstrx_contacts_prototype}->{field_spec}->[3] = "integer unsigned not null" ;
   $tables->{interface_secstrx_contacts_prototype}->{field_name}->[4] = "sse_id_2" ;
   $tables->{interface_secstrx_contacts_prototype}->{field_spec}->[4] = "integer unsigned not null" ;
   $tables->{interface_secstrx_contacts_prototype}->{field_name}->[5] = "num_res_1" ;
   $tables->{interface_secstrx_contacts_prototype}->{field_spec}->[5] = "integer unsigned not null" ;
   $tables->{interface_secstrx_contacts_prototype}->{field_name}->[6] = "num_res_2" ;
   $tables->{interface_secstrx_contacts_prototype}->{field_spec}->[6] = "integer unsigned not null" ;
   $tables->{interface_secstrx_contacts_prototype}->{field_name}->[7] = "num_res_12" ;
   $tables->{interface_secstrx_contacts_prototype}->{field_spec}->[7] = "integer unsigned not null" ;
   $tables->{interface_secstrx_contacts_prototype}->{field_name}->[8] = "sse_1" ;
   $tables->{interface_secstrx_contacts_prototype}->{field_spec}->[8] = "enum('H', 'B', 'E', 'G', 'I', 'T', 'S', ' ') not null" ;
   $tables->{interface_secstrx_contacts_prototype}->{field_name}->[9] = "sse_2" ;
   $tables->{interface_secstrx_contacts_prototype}->{field_spec}->[9] = "enum('H', 'B', 'E', 'G', 'I', 'T', 'S', ' ') not null" ;
   $tables->{interface_sc}->{prikey} = "(subset_id_1,subset_id_2)" ;
   $tables->{interface_sc}->{field_name}->[0] = "bdp_id" ;
   $tables->{interface_sc}->{field_spec}->[0] = "integer unsigned not null" ;
   $tables->{interface_sc}->{field_name}->[1] = "subset_id_1" ;
   $tables->{interface_sc}->{field_spec}->[1] = "char(50) not null" ;
   $tables->{interface_sc}->{field_name}->[2] = "subset_id_2" ;
   $tables->{interface_sc}->{field_spec}->[2] = "char(50) not null" ;
   $tables->{interface_sc}->{field_name}->[3] = "sc" ;
   $tables->{interface_sc}->{field_spec}->[3] = "float" ;
   $tables->{interface_sc}->{field_name}->[4] = "median_dist" ;
   $tables->{interface_sc}->{field_spec}->[4] = "float" ;
   $tables->{interface_contacts_special_prototype}->{prikey} = "(bdp_id,subset_id_1,subset_id_2,chain_no_1,resno_1,chain_no_2,resno_2,type)" ;
   $tables->{interface_contacts_special_prototype}->{field_name}->[0] = "bdp_id" ;
   $tables->{interface_contacts_special_prototype}->{field_spec}->[0] = "integer unsigned not null" ;
   $tables->{interface_contacts_special_prototype}->{field_name}->[1] = "subset_id_1" ;
   $tables->{interface_contacts_special_prototype}->{field_spec}->[1] = "char(50) not null" ;
   $tables->{interface_contacts_special_prototype}->{field_name}->[2] = "subset_id_2" ;
   $tables->{interface_contacts_special_prototype}->{field_spec}->[2] = "char(50) not null" ;
   $tables->{interface_contacts_special_prototype}->{field_name}->[3] = "chain_no_1" ;
   $tables->{interface_contacts_special_prototype}->{field_spec}->[3] = "integer unsigned not null" ;
   $tables->{interface_contacts_special_prototype}->{field_name}->[4] = "chain_id_1" ;
   $tables->{interface_contacts_special_prototype}->{field_spec}->[4] = "char(1) not null" ;
   $tables->{interface_contacts_special_prototype}->{field_name}->[5] = "resno_1" ;
   $tables->{interface_contacts_special_prototype}->{field_spec}->[5] = "char(10) not null" ;
   $tables->{interface_contacts_special_prototype}->{field_name}->[6] = "resna_1" ;
   $tables->{interface_contacts_special_prototype}->{field_spec}->[6] = "char(3) not null" ;
   $tables->{interface_contacts_special_prototype}->{field_name}->[7] = "chain_no_2" ;
   $tables->{interface_contacts_special_prototype}->{field_spec}->[7] = "integer unsigned not null" ;
   $tables->{interface_contacts_special_prototype}->{field_name}->[8] = "chain_id_2" ;
   $tables->{interface_contacts_special_prototype}->{field_spec}->[8] = "char(1) not null" ;
   $tables->{interface_contacts_special_prototype}->{field_name}->[9] = "resno_2" ;
   $tables->{interface_contacts_special_prototype}->{field_spec}->[9] = "char(10) not null" ;
   $tables->{interface_contacts_special_prototype}->{field_name}->[10] = "resna_2" ;
   $tables->{interface_contacts_special_prototype}->{field_spec}->[10] = "char(3) not null" ;
   $tables->{interface_contacts_special_prototype}->{field_name}->[11] = "type" ;
   $tables->{interface_contacts_special_prototype}->{field_spec}->[11] = "enum('salt', 'hbond', 'ssbond')" ;
   $tables->{cath_names}->{prikey} = "(class,arch,topol,homol)" ;
   $tables->{cath_names}->{field_name}->[0] = "class" ;
   $tables->{cath_names}->{field_spec}->[0] = "smallint unsigned not null" ;
   $tables->{cath_names}->{field_name}->[1] = "arch" ;
   $tables->{cath_names}->{field_spec}->[1] = "smallint unsigned not null" ;
   $tables->{cath_names}->{field_name}->[2] = "topol" ;
   $tables->{cath_names}->{field_spec}->[2] = "smallint unsigned not null" ;
   $tables->{cath_names}->{field_name}->[3] = "homol" ;
   $tables->{cath_names}->{field_spec}->[3] = "smallint unsigned not null" ;
   $tables->{cath_names}->{field_name}->[4] = "representative_dom" ;
   $tables->{cath_names}->{field_spec}->[4] = "char(7) not null" ;
   $tables->{cath_names}->{field_name}->[5] = "description" ;
   $tables->{cath_names}->{field_spec}->[5] = "tinytext" ;
   $tables->{interface_contacts_special_tables}->{prikey} = "(bdp_id)" ;
   $tables->{interface_contacts_special_tables}->{field_name}->[0] = "bdp_id" ;
   $tables->{interface_contacts_special_tables}->{field_spec}->[0] = "integer unsigned not null" ;
   $tables->{interface_contacts_special_tables}->{field_name}->[1] = "cutoff" ;
   $tables->{interface_contacts_special_tables}->{field_spec}->[1] = "float not null" ;
   $tables->{interface_contacts_special_tables}->{field_name}->[2] = "table_name" ;
   $tables->{interface_contacts_special_tables}->{field_spec}->[2] = "char(40)" ;
   $tables->{interface_contacts_special_tables}->{field_name}->[3] = "source_file" ;
   $tables->{interface_contacts_special_tables}->{field_spec}->[3] = "char(255)" ;
   $tables->{interface_secstrx_prototype}->{prikey} = "(bdp_id,subset_id_1,subset_id_2,subset_id)" ;
   $tables->{interface_secstrx_prototype}->{field_name}->[0] = "bdp_id" ;
   $tables->{interface_secstrx_prototype}->{field_spec}->[0] = "integer unsigned not null" ;
   $tables->{interface_secstrx_prototype}->{field_name}->[1] = "subset_id_1" ;
   $tables->{interface_secstrx_prototype}->{field_spec}->[1] = "char(50) not null" ;
   $tables->{interface_secstrx_prototype}->{field_name}->[2] = "subset_id_2" ;
   $tables->{interface_secstrx_prototype}->{field_spec}->[2] = "char(50) not null" ;
   $tables->{interface_secstrx_prototype}->{field_name}->[3] = "subset_id" ;
   $tables->{interface_secstrx_prototype}->{field_spec}->[3] = "char(50) not null" ;
   $tables->{interface_secstrx_prototype}->{field_name}->[4] = "sse_id" ;
   $tables->{interface_secstrx_prototype}->{field_spec}->[4] = "integer unsigned not null" ;
   $tables->{interface_secstrx_prototype}->{field_name}->[5] = "num_res" ;
   $tables->{interface_secstrx_prototype}->{field_spec}->[5] = "integer unsigned not null" ;
   $tables->{interface_secstrx_prototype}->{field_name}->[6] = "sse" ;
   $tables->{interface_secstrx_prototype}->{field_spec}->[6] = "enum('H', 'B', 'E', 'G', 'I', 'T', 'S', ' ') not null" ;
   $tables->{scop_hie}->{prikey} = "(self_sun_id)" ;
   $tables->{scop_hie}->{field_name}->[0] = "self_sun_id" ;
   $tables->{scop_hie}->{field_spec}->[0] = "mediumint unsigned not null" ;
   $tables->{scop_hie}->{field_name}->[1] = "parent_sun_id" ;
   $tables->{scop_hie}->{field_spec}->[1] = "mediumint unsigned" ;
   $tables->{scop_hie}->{field_name}->[2] = "kids_sun_id" ;
   $tables->{scop_hie}->{field_spec}->[2] = "text" ;
   $tables->{subsets_sasa}->{prikey} = "(subset_id)" ;
   $tables->{subsets_sasa}->{field_name}->[0] = "bdp_id" ;
   $tables->{subsets_sasa}->{field_spec}->[0] = "integer unsigned not null" ;
   $tables->{subsets_sasa}->{field_name}->[1] = "subset_id" ;
   $tables->{subsets_sasa}->{field_spec}->[1] = "char(50) not null" ;
   $tables->{subsets_sasa}->{field_name}->[2] = "sasa_all" ;
   $tables->{subsets_sasa}->{field_spec}->[2] = "float not null" ;
   $tables->{subsets_sasa}->{field_name}->[3] = "sasa_sc" ;
   $tables->{subsets_sasa}->{field_spec}->[3] = "float not null" ;
   $tables->{subsets_sasa}->{field_name}->[4] = "sasa_mc" ;
   $tables->{subsets_sasa}->{field_spec}->[4] = "float not null" ;
   $tables->{subsets_sasa}->{field_name}->[5] = "sasa_polar" ;
   $tables->{subsets_sasa}->{field_spec}->[5] = "float not null" ;
   $tables->{subsets_sasa}->{field_name}->[6] = "sasa_nonpolar" ;
   $tables->{subsets_sasa}->{field_spec}->[6] = "float not null" ;
   $tables->{interface_vol}->{prikey} = "(subset_id_1,subset_id_2)" ;
   $tables->{interface_vol}->{field_name}->[0] = "bdp_id" ;
   $tables->{interface_vol}->{field_spec}->[0] = "integer unsigned not null" ;
   $tables->{interface_vol}->{field_name}->[1] = "subset_id_1" ;
   $tables->{interface_vol}->{field_spec}->[1] = "char(50) not null" ;
   $tables->{interface_vol}->{field_name}->[2] = "subset_id_2" ;
   $tables->{interface_vol}->{field_spec}->[2] = "char(50) not null" ;
   $tables->{interface_vol}->{field_name}->[3] = "vol_12" ;
   $tables->{interface_vol}->{field_spec}->[3] = "float not null" ;
   $tables->{interface_vol}->{field_name}->[4] = "delta_vol" ;
   $tables->{interface_vol}->{field_spec}->[4] = "float not null" ;
   $tables->{bdp_interaction_topology}->{prikey} = "(bdp_id,subset_source_id)" ;
   $tables->{bdp_interaction_topology}->{field_name}->[0] = "bdp_id" ;
   $tables->{bdp_interaction_topology}->{field_spec}->[0] = "integer unsigned not null" ;
   $tables->{bdp_interaction_topology}->{field_name}->[1] = "subset_source_id" ;
   $tables->{bdp_interaction_topology}->{field_spec}->[1] = "integer unsigned not null" ;
   $tables->{bdp_interaction_topology}->{field_name}->[2] = "num_domains" ;
   $tables->{bdp_interaction_topology}->{field_spec}->[2] = "integer unsigned not null" ;
   $tables->{bdp_interaction_topology}->{field_name}->[3] = "num_edges" ;
   $tables->{bdp_interaction_topology}->{field_spec}->[3] = "integer unsigned not null" ;
   $tables->{bdp_interaction_topology}->{field_name}->[4] = "num_domain_classes" ;
   $tables->{bdp_interaction_topology}->{field_spec}->[4] = "integer unsigned not null" ;
   $tables->{bdp_interaction_topology}->{field_name}->[5] = "nodestring" ;
   $tables->{bdp_interaction_topology}->{field_spec}->[5] = "text not null" ;
   $tables->{bdp_interaction_topology}->{field_name}->[6] = "edgestring" ;
   $tables->{bdp_interaction_topology}->{field_spec}->[6] = "text not null" ;
   $tables->{bindingsite_contacts_prototype}->{prikey} = "(bdp_id,subset_id_1,subset_id,resno_1,resno_2)" ;
   $tables->{bindingsite_contacts_prototype}->{field_name}->[0] = "bdp_id" ;
   $tables->{bindingsite_contacts_prototype}->{field_spec}->[0] = "integer unsigned not null" ;
   $tables->{bindingsite_contacts_prototype}->{field_name}->[1] = "subset_id_1" ;
   $tables->{bindingsite_contacts_prototype}->{field_spec}->[1] = "char(50) not null" ;
   $tables->{bindingsite_contacts_prototype}->{field_name}->[2] = "subset_id_2" ;
   $tables->{bindingsite_contacts_prototype}->{field_spec}->[2] = "char(50) not null" ;
   $tables->{bindingsite_contacts_prototype}->{field_name}->[3] = "subset_id" ;
   $tables->{bindingsite_contacts_prototype}->{field_spec}->[3] = "char(50) not null" ;
   $tables->{bindingsite_contacts_prototype}->{field_name}->[4] = "chain_id_1" ;
   $tables->{bindingsite_contacts_prototype}->{field_spec}->[4] = "char(1) not null" ;
   $tables->{bindingsite_contacts_prototype}->{field_name}->[5] = "resno_1" ;
   $tables->{bindingsite_contacts_prototype}->{field_spec}->[5] = "char(10) not null" ;
   $tables->{bindingsite_contacts_prototype}->{field_name}->[6] = "resna_1" ;
   $tables->{bindingsite_contacts_prototype}->{field_spec}->[6] = "char(3) not null" ;
   $tables->{bindingsite_contacts_prototype}->{field_name}->[7] = "chain_id_2" ;
   $tables->{bindingsite_contacts_prototype}->{field_spec}->[7] = "char(1) not null" ;
   $tables->{bindingsite_contacts_prototype}->{field_name}->[8] = "resno_2" ;
   $tables->{bindingsite_contacts_prototype}->{field_spec}->[8] = "char(10) not null" ;
   $tables->{bindingsite_contacts_prototype}->{field_name}->[9] = "resna_2" ;
   $tables->{bindingsite_contacts_prototype}->{field_spec}->[9] = "char(3) not null" ;
   $tables->{bindingsite_contacts_prototype}->{field_name}->[10] = "min_dist" ;
   $tables->{bindingsite_contacts_prototype}->{field_spec}->[10] = "float not null" ;
   $tables->{bindingsite_contacts_prototype}->{field_name}->[11] = "num_contacts" ;
   $tables->{bindingsite_contacts_prototype}->{field_spec}->[11] = "integer unsigned not null" ;
   $tables->{bindingsite_contacts_prototype}->{field_name}->[12] = "num_contacts_4" ;
   $tables->{bindingsite_contacts_prototype}->{field_spec}->[12] = "integer unsigned not null" ;
   $tables->{bindingsite_contacts_prototype}->{field_name}->[13] = "num_contacts_4p5" ;
   $tables->{bindingsite_contacts_prototype}->{field_spec}->[13] = "integer unsigned not null" ;
   $tables->{bindingsite_contacts_prototype}->{field_name}->[14] = "num_contacts_5" ;
   $tables->{bindingsite_contacts_prototype}->{field_spec}->[14] = "integer unsigned not null" ;
   $tables->{bindingsite_contacts_prototype}->{field_name}->[15] = "num_contacts_5p5" ;
   $tables->{bindingsite_contacts_prototype}->{field_spec}->[15] = "integer unsigned not null" ;
   $tables->{interface_sse_topology}->{prikey} = "(bdp_id,subset_id_1,subset_id_2)" ;
   $tables->{interface_sse_topology}->{field_name}->[0] = "bdp_id" ;
   $tables->{interface_sse_topology}->{field_spec}->[0] = "integer unsigned not null" ;
   $tables->{interface_sse_topology}->{field_name}->[1] = "subset_id_1" ;
   $tables->{interface_sse_topology}->{field_spec}->[1] = "char(50) not null" ;
   $tables->{interface_sse_topology}->{field_name}->[2] = "subset_id_2" ;
   $tables->{interface_sse_topology}->{field_spec}->[2] = "char(50) not null" ;
   $tables->{interface_sse_topology}->{field_name}->[3] = "num_domains" ;
   $tables->{interface_sse_topology}->{field_spec}->[3] = "integer unsigned not null" ;
   $tables->{interface_sse_topology}->{field_name}->[4] = "num_edges" ;
   $tables->{interface_sse_topology}->{field_spec}->[4] = "integer unsigned not null" ;
   $tables->{interface_sse_topology}->{field_name}->[5] = "nodelist" ;
   $tables->{interface_sse_topology}->{field_spec}->[5] = "text not null" ;
   $tables->{interface_sse_topology}->{field_name}->[6] = "nodelist_sse" ;
   $tables->{interface_sse_topology}->{field_spec}->[6] = "text not null" ;
   $tables->{interface_sse_topology}->{field_name}->[7] = "edgelist" ;
   $tables->{interface_sse_topology}->{field_spec}->[7] = "text not null" ;
   $tables->{interface_sse_topology}->{field_name}->[8] = "edgelist_sse" ;
   $tables->{interface_sse_topology}->{field_spec}->[8] = "text not null" ;
   $tables->{subsets_source}->{prikey} = "(subset_source_id,subset_source)" ;
   $tables->{subsets_source}->{field_name}->[0] = "subset_source_id" ;
   $tables->{subsets_source}->{field_spec}->[0] = "integer unsigned auto_increment not null" ;
   $tables->{subsets_source}->{field_name}->[1] = "subset_source" ;
   $tables->{subsets_source}->{field_spec}->[1] = "char(50) not null" ;
   $tables->{subsets_source}->{field_name}->[2] = "version" ;
   $tables->{subsets_source}->{field_spec}->[2] = "char(20) not null" ;
   $tables->{subsets_class}->{prikey} = "(subset_source_id,class)" ;
   $tables->{subsets_class}->{field_name}->[0] = "subset_source_id" ;
   $tables->{subsets_class}->{field_spec}->[0] = "integer unsigned not null" ;
   $tables->{subsets_class}->{field_name}->[1] = "class" ;
   $tables->{subsets_class}->{field_spec}->[1] = "char(70) not null" ;
   $tables->{subsets_class}->{field_name}->[2] = "description" ;
   $tables->{subsets_class}->{field_spec}->[2] = "char(250)" ;
   $tables->{bdp_numres}->{prikey} = "(bdp_id)" ;
   $tables->{bdp_numres}->{field_name}->[0] = "bdp_id" ;
   $tables->{bdp_numres}->{field_spec}->[0] = "integer unsigned not null" ;
   $tables->{bdp_numres}->{field_name}->[1] = "num_res" ;
   $tables->{bdp_numres}->{field_spec}->[1] = "integer not null" ;
   $tables->{scop_des}->{prikey} = "(sun_id)" ;
   $tables->{scop_des}->{field_name}->[0] = "sun_id" ;
   $tables->{scop_des}->{field_spec}->[0] = "mediumint unsigned not null" ;
   $tables->{scop_des}->{field_name}->[1] = "entry_type" ;
   $tables->{scop_des}->{field_spec}->[1] = "enum('cl', 'cf', 'sf', 'fa', 'dm' ,'sp', 'px') not null" ;
   $tables->{scop_des}->{field_name}->[2] = "class_id" ;
   $tables->{scop_des}->{field_spec}->[2] = "char(1) not null" ;
   $tables->{scop_des}->{field_name}->[3] = "fold_id" ;
   $tables->{scop_des}->{field_spec}->[3] = "smallint unsigned not null" ;
   $tables->{scop_des}->{field_name}->[4] = "superfam_id" ;
   $tables->{scop_des}->{field_spec}->[4] = "smallint unsigned not null" ;
   $tables->{scop_des}->{field_name}->[5] = "fam_id" ;
   $tables->{scop_des}->{field_spec}->[5] = "smallint unsigned not null" ;
   $tables->{scop_des}->{field_name}->[6] = "scop_id" ;
   $tables->{scop_des}->{field_spec}->[6] = "char(7)" ;
   $tables->{scop_des}->{field_name}->[7] = "description" ;
   $tables->{scop_des}->{field_spec}->[7] = "char(250)" ;
   $tables->{pqs_asalist}->{prikey} = "(pqs_id)" ;
   $tables->{pqs_asalist}->{field_name}->[0] = "pqs_id" ;
   $tables->{pqs_asalist}->{field_spec}->[0] = "char(10) not null" ;
   $tables->{pqs_asalist}->{field_name}->[1] = "pdb_id" ;
   $tables->{pqs_asalist}->{field_spec}->[1] = "char(4) not null" ;
   $tables->{pqs_asalist}->{field_name}->[2] = "split_no" ;
   $tables->{pqs_asalist}->{field_spec}->[2] = "smallint unsigned" ;
   $tables->{pqs_asalist}->{field_name}->[3] = "contact_type" ;
   $tables->{pqs_asalist}->{field_spec}->[3] = "enum('OK', 'XPACK') not null" ;
   $tables->{pqs_asalist}->{field_name}->[4] = "homo_hetero" ;
   $tables->{pqs_asalist}->{field_spec}->[4] = "enum('HETERO', 'HOMO') not null" ;
   $tables->{pqs_asalist}->{field_name}->[5] = "quaternary_state" ;
   $tables->{pqs_asalist}->{field_spec}->[5] = "char(30)" ;
   $tables->{pqs_asalist}->{field_name}->[6] = "num_chain" ;
   $tables->{pqs_asalist}->{field_spec}->[6] = "smallint unsigned" ;
   $tables->{pqs_asalist}->{field_name}->[7] = "num_resid" ;
   $tables->{pqs_asalist}->{field_spec}->[7] = "smallint unsigned" ;
   $tables->{pqs_asalist}->{field_name}->[8] = "num_hetatm" ;
   $tables->{pqs_asalist}->{field_spec}->[8] = "smallint unsigned" ;
   $tables->{pqs_asalist}->{field_name}->[9] = "delta_asa" ;
   $tables->{pqs_asalist}->{field_spec}->[9] = "float" ;
   $tables->{pqs_asalist}->{field_name}->[10] = "delta_solvat_e" ;
   $tables->{pqs_asalist}->{field_spec}->[10] = "float" ;
   $tables->{pqs_asalist}->{field_name}->[11] = "num_buried" ;
   $tables->{pqs_asalist}->{field_spec}->[11] = "smallint unsigned" ;
   $tables->{pqs_asalist}->{field_name}->[12] = "num_salt_bridges" ;
   $tables->{pqs_asalist}->{field_spec}->[12] = "smallint unsigned" ;
   $tables->{pqs_asalist}->{field_name}->[13] = "num_disulfide" ;
   $tables->{pqs_asalist}->{field_spec}->[13] = "smallint unsigned" ;
   $tables->{pqs_asalist}->{field_name}->[14] = "rmsd" ;
   $tables->{pqs_asalist}->{field_spec}->[14] = "char(10)" ;
   $tables->{pqs_asalist}->{field_name}->[15] = "assembly_formula" ;
   $tables->{pqs_asalist}->{field_spec}->[15] = "char(75)" ;
   $tables->{pqs_asalist}->{field_name}->[16] = "calpha_only" ;
   $tables->{pqs_asalist}->{field_spec}->[16] = "bool" ;
   $tables->{interface_secstrx_basic_contacts_prototype}->{prikey} = "(bdp_id,subset_id_1,subset_id_2)" ;
   $tables->{interface_secstrx_basic_contacts_prototype}->{field_name}->[0] = "bdp_id" ;
   $tables->{interface_secstrx_basic_contacts_prototype}->{field_spec}->[0] = "integer unsigned not null" ;
   $tables->{interface_secstrx_basic_contacts_prototype}->{field_name}->[1] = "subset_id_1" ;
   $tables->{interface_secstrx_basic_contacts_prototype}->{field_spec}->[1] = "char(50) not null" ;
   $tables->{interface_secstrx_basic_contacts_prototype}->{field_name}->[2] = "subset_id_2" ;
   $tables->{interface_secstrx_basic_contacts_prototype}->{field_spec}->[2] = "char(50) not null" ;
   $tables->{interface_secstrx_basic_contacts_prototype}->{field_name}->[3] = "sse_id_1" ;
   $tables->{interface_secstrx_basic_contacts_prototype}->{field_spec}->[3] = "integer unsigned not null" ;
   $tables->{interface_secstrx_basic_contacts_prototype}->{field_name}->[4] = "sse_id_2" ;
   $tables->{interface_secstrx_basic_contacts_prototype}->{field_spec}->[4] = "integer unsigned not null" ;
   $tables->{interface_secstrx_basic_contacts_prototype}->{field_name}->[5] = "num_res_1" ;
   $tables->{interface_secstrx_basic_contacts_prototype}->{field_spec}->[5] = "integer unsigned not null" ;
   $tables->{interface_secstrx_basic_contacts_prototype}->{field_name}->[6] = "num_res_2" ;
   $tables->{interface_secstrx_basic_contacts_prototype}->{field_spec}->[6] = "integer unsigned not null" ;
   $tables->{interface_secstrx_basic_contacts_prototype}->{field_name}->[7] = "num_res_12" ;
   $tables->{interface_secstrx_basic_contacts_prototype}->{field_spec}->[7] = "integer unsigned not null" ;
   $tables->{interface_secstrx_basic_contacts_prototype}->{field_name}->[8] = "sse_1" ;
   $tables->{interface_secstrx_basic_contacts_prototype}->{field_spec}->[8] = "enum('H', 'B', 'T', ' ') not null" ;
   $tables->{interface_secstrx_basic_contacts_prototype}->{field_name}->[9] = "sse_2" ;
   $tables->{interface_secstrx_basic_contacts_prototype}->{field_spec}->[9] = "enum('H', 'B', 'T', ' ') not null" ;
   $tables->{subsets_residues_prototype}->{prikey} = "(bdp_id,chain_no,resno_serial,subset_id)" ;
   $tables->{subsets_residues_prototype}->{field_name}->[0] = "bdp_id" ;
   $tables->{subsets_residues_prototype}->{field_spec}->[0] = "integer unsigned not null" ;
   $tables->{subsets_residues_prototype}->{field_name}->[1] = "chain_no" ;
   $tables->{subsets_residues_prototype}->{field_spec}->[1] = "integer unsigned not null" ;
   $tables->{subsets_residues_prototype}->{field_name}->[2] = "chain_id" ;
   $tables->{subsets_residues_prototype}->{field_spec}->[2] = "char(1) not null" ;
   $tables->{subsets_residues_prototype}->{field_name}->[3] = "resno_serial" ;
   $tables->{subsets_residues_prototype}->{field_spec}->[3] = "integer unsigned not null" ;
   $tables->{subsets_residues_prototype}->{field_name}->[4] = "resno" ;
   $tables->{subsets_residues_prototype}->{field_spec}->[4] = "char(10) not null" ;
   $tables->{subsets_residues_prototype}->{field_name}->[5] = "resno_int" ;
   $tables->{subsets_residues_prototype}->{field_spec}->[5] = "integer not null" ;
   $tables->{subsets_residues_prototype}->{field_name}->[6] = "subset_id" ;
   $tables->{subsets_residues_prototype}->{field_spec}->[6] = "char(50) not null" ;
   $tables->{intersubset_contacts}->{prikey} = "(subset_id_1,subset_id_2)" ;
   $tables->{intersubset_contacts}->{field_name}->[0] = "bdp_id" ;
   $tables->{intersubset_contacts}->{field_spec}->[0] = "integer unsigned not null" ;
   $tables->{intersubset_contacts}->{field_name}->[1] = "subset_id_1" ;
   $tables->{intersubset_contacts}->{field_spec}->[1] = "char(50) binary not null" ;
   $tables->{intersubset_contacts}->{field_name}->[2] = "class_1" ;
   $tables->{intersubset_contacts}->{field_spec}->[2] = "char(70) not null" ;
   $tables->{intersubset_contacts}->{field_name}->[3] = "subset_id_2" ;
   $tables->{intersubset_contacts}->{field_spec}->[3] = "char(50) binary not null" ;
   $tables->{intersubset_contacts}->{field_name}->[4] = "class_2" ;
   $tables->{intersubset_contacts}->{field_spec}->[4] = "char(70) not null" ;
   $tables->{intersubset_contacts}->{field_name}->[5] = "num_contacts" ;
   $tables->{intersubset_contacts}->{field_spec}->[5] = "integer unsigned not null" ;
   $tables->{intersubset_contacts}->{field_name}->[6] = "cutoff" ;
   $tables->{intersubset_contacts}->{field_spec}->[6] = "float not null" ;
   $tables->{intersubset_contacts}->{field_name}->[7] = "num_contacts_4" ;
   $tables->{intersubset_contacts}->{field_spec}->[7] = "integer unsigned not null" ;
   $tables->{intersubset_contacts}->{field_name}->[8] = "num_contacts_4p5" ;
   $tables->{intersubset_contacts}->{field_spec}->[8] = "integer unsigned not null" ;
   $tables->{intersubset_contacts}->{field_name}->[9] = "num_contacts_5" ;
   $tables->{intersubset_contacts}->{field_spec}->[9] = "integer unsigned not null" ;
   $tables->{intersubset_contacts}->{field_name}->[10] = "num_contacts_5p5" ;
   $tables->{intersubset_contacts}->{field_spec}->[10] = "integer unsigned not null" ;
   $tables->{intersubset_contacts}->{field_name}->[11] = "hbond" ;
   $tables->{intersubset_contacts}->{field_spec}->[11] = "integer unsigned not null" ;
   $tables->{intersubset_contacts}->{field_name}->[12] = "salt" ;
   $tables->{intersubset_contacts}->{field_spec}->[12] = "integer unsigned not null" ;
   $tables->{intersubset_contacts}->{field_name}->[13] = "ssbond" ;
   $tables->{intersubset_contacts}->{field_spec}->[13] = "integer unsigned not null" ;
   $tables->{intersubset_contacts}->{field_name}->[14] = "chains" ;
   $tables->{intersubset_contacts}->{field_spec}->[14] = "enum('same', 'diff', 'both') not null" ;
   $tables->{subsets_vol}->{prikey} = "(subset_id)" ;
   $tables->{subsets_vol}->{field_name}->[0] = "subset_id" ;
   $tables->{subsets_vol}->{field_spec}->[0] = "char(50) binary not null" ;
   $tables->{subsets_vol}->{field_name}->[1] = "volume" ;
   $tables->{subsets_vol}->{field_spec}->[1] = "float not null" ;
   $tables->{subsets_details}->{prikey} = "(subset_id,segment_id)" ;
   $tables->{subsets_details}->{field_name}->[0] = "subset_id" ;
   $tables->{subsets_details}->{field_spec}->[0] = "char(50) binary not null" ;
   $tables->{subsets_details}->{field_name}->[1] = "segment_id" ;
   $tables->{subsets_details}->{field_spec}->[1] = "char(50) not null" ;
   $tables->{subsets_details}->{field_name}->[2] = "chain_no" ;
   $tables->{subsets_details}->{field_spec}->[2] = "integer unsigned" ;
   $tables->{subsets_details}->{field_name}->[3] = "chain_id" ;
   $tables->{subsets_details}->{field_spec}->[3] = "char(1) not null" ;
   $tables->{subsets_details}->{field_name}->[4] = "start_resno" ;
   $tables->{subsets_details}->{field_spec}->[4] = "char(10) not null" ;
   $tables->{subsets_details}->{field_name}->[5] = "start_resno_int" ;
   $tables->{subsets_details}->{field_spec}->[5] = "integer" ;
   $tables->{subsets_details}->{field_name}->[6] = "end_resno" ;
   $tables->{subsets_details}->{field_spec}->[6] = "char(10) not null" ;
   $tables->{subsets_details}->{field_name}->[7] = "end_resno_int" ;
   $tables->{subsets_details}->{field_spec}->[7] = "integer" ;
   $tables->{subsets_env}->{prikey} = "(subset_id)" ;
   $tables->{subsets_env}->{field_name}->[0] = "bdp_id" ;
   $tables->{subsets_env}->{field_spec}->[0] = "integer not null" ;
   $tables->{subsets_env}->{field_name}->[1] = "subset_id" ;
   $tables->{subsets_env}->{field_spec}->[1] = "char(50) binary not null" ;
   $tables->{subsets_env}->{field_name}->[2] = "num_domains_bdp" ;
   $tables->{subsets_env}->{field_spec}->[2] = "integer not null" ;
   $tables->{subsets_env}->{field_name}->[3] = "num_interactions" ;
   $tables->{subsets_env}->{field_spec}->[3] = "integer not null" ;
   $tables->{subsets_env}->{field_name}->[4] = "num_domains_chain" ;
   $tables->{subsets_env}->{field_spec}->[4] = "integer not null" ;
   $tables->{bdp_secstrx_prototype}->{prikey} = "(bdp_id)" ;
   $tables->{bdp_secstrx_prototype}->{field_name}->[0] = "bdp_id" ;
   $tables->{bdp_secstrx_prototype}->{field_spec}->[0] = "integer unsigned not null" ;
   $tables->{bdp_secstrx_prototype}->{field_name}->[1] = "chain_id" ;
   $tables->{bdp_secstrx_prototype}->{field_spec}->[1] = "char(1) not null" ;
   $tables->{bdp_secstrx_prototype}->{field_name}->[2] = "resno" ;
   $tables->{bdp_secstrx_prototype}->{field_spec}->[2] = "char(10) not null" ;
   $tables->{bdp_secstrx_prototype}->{field_name}->[3] = "sse" ;
   $tables->{bdp_secstrx_prototype}->{field_spec}->[3] = "enum('H', 'B', 'E', 'G', 'I', 'T', 'S', ' ') not null" ;
   $tables->{bdp_secstrx_prototype}->{field_name}->[4] = "sse_basic" ;
   $tables->{bdp_secstrx_prototype}->{field_spec}->[4] = "enum('H', 'B', 'T', ' ') not null" ;
   $tables->{bdp_secstrx_prototype}->{field_name}->[5] = "sse_id" ;
   $tables->{bdp_secstrx_prototype}->{field_spec}->[5] = "integer unsigned not null" ;
   $tables->{bdp_secstrx_prototype}->{field_name}->[6] = "sse_basic_id" ;
   $tables->{bdp_secstrx_prototype}->{field_spec}->[6] = "integer unsigned not null" ;
   $tables->{pdb_entries}->{prikey} = "(pdb_id)" ;
   $tables->{pdb_entries}->{field_name}->[0] = "pdb_id" ;
   $tables->{pdb_entries}->{field_spec}->[0] = "char(4) not null" ;
   $tables->{pdb_entries}->{field_name}->[1] = "header" ;
   $tables->{pdb_entries}->{field_spec}->[1] = "char(100)" ;
   $tables->{pdb_entries}->{field_name}->[2] = "accession_date" ;
   $tables->{pdb_entries}->{field_spec}->[2] = "date" ;
   $tables->{pdb_entries}->{field_name}->[3] = "compound" ;
   $tables->{pdb_entries}->{field_spec}->[3] = "text" ;
   $tables->{pdb_entries}->{field_name}->[4] = "source" ;
   $tables->{pdb_entries}->{field_spec}->[4] = "char(150)" ;
   $tables->{pdb_entries}->{field_name}->[5] = "author" ;
   $tables->{pdb_entries}->{field_spec}->[5] = "text" ;
   $tables->{pdb_entries}->{field_name}->[6] = "resolution" ;
   $tables->{pdb_entries}->{field_spec}->[6] = "float" ;
   $tables->{pdb_entries}->{field_name}->[7] = "experiment_type" ;
   $tables->{pdb_entries}->{field_spec}->[7] = "char(90)" ;
   $tables->{interface_secstrx_tables}->{prikey} = "(bdp_id)" ;
   $tables->{interface_secstrx_tables}->{field_name}->[0] = "bdp_id" ;
   $tables->{interface_secstrx_tables}->{field_spec}->[0] = "integer unsigned not null" ;
   $tables->{interface_secstrx_tables}->{field_name}->[1] = "table_name" ;
   $tables->{interface_secstrx_tables}->{field_spec}->[1] = "char(40)" ;
   $tables->{interface_secstrx_tables}->{field_name}->[2] = "source_file" ;
   $tables->{interface_secstrx_tables}->{field_spec}->[2] = "char(255)" ;

   return $tables ;
}


=head2 table_spec(@tablelist)

   Title:       table_spec()
   Function:    Return mysql DDL format table specs.
   Returns:     $_ hashref pointing from tablename to specs
                  $_->{i} = specs for ith table

=cut

sub table_spec {

   my @tables = @_ ;
   my $all_tables = full_table_specs() ;

   my $ans ;

   foreach my $j ( 0 ..$#tables) {
      if (exists $all_tables->{$tables[$j]}) {
         push @{$ans}, {name => $tables[$j], specs => $all_tables->{$tables[$j]}} ;
      }
   }

   return $ans ;

}



=head2 SUB sql_table_spec(@tablelist)

   Title:       sql_table_spec()
   Function:    Return mysql DDL format table specs.
   Args:        $_ hashref pointing from tablename to specs
                $_->{i} = specs for ith table

=cut

sub sql_table_spec {

   my @tables = @_ ;
   my $all_tables = full_table_specs() ;

   my $ans ;

   foreach my $j ( 0 ..$#tables) {
      if (exists $all_tables->{$tables[$j]}) {
         my $cursql ;
         $cursql = '('."\n" ;

	 foreach my $k ( 0 .. $#{$all_tables->{$tables[$j]}->{field_name}}) {
	    $cursql .= "   ".$all_tables->{$tables[$j]}->{field_name}->[$k]." ".$all_tables->{$tables[$j]}->{field_spec}->[$k] ;
	    if ( $k < $#{$all_tables->{$tables[$j]}->{field_name}} ||
	         exists $all_tables->{$tables[$j]}->{prikey} ) {
	       $cursql .= ",\n" ;
	    } else {
	       $cursql .= "\n" ;
	    }
	 }

	 if (exists $all_tables->{$tables[$j]}->{prikey}) {
	    $cursql .= '   PRIMARY KEY '.$all_tables->{$tables[$j]}->{prikey}."\n" ; }

	 $cursql .= ")" ;
         $ans->{$tables[$j]} = $cursql ;
      }
   }

   return $ans ;

}


=head2 rawselect_tod()

   Title:       rawselect_tod()
   Function:    performs basic SELECT statement queries on a table on disk
   Args:        $_[0] = SQL SELECT-like statement
                $_[1] = fullfile

=cut

sub rawselect_tod {

   my $select_sql = shift ;
   my $fullfile ;
   if ($#_ >= 0) {
      $fullfile = shift ; }

   my $modtie_specs = modtie::set_modtie_specs() ;
   my $specs = $modtie_specs->{pibase_specs} ;

   my $bins = modtie::locate_binaries() ;

   my ($fields, $table) =
      ($select_sql =~ /(?:SELECT|select) (.+) (?:FROM|from) (\w+)$/) ;

   $fields =~ s/ //g ;
   my @fields = split(/\,/, $fields) ;

   my $tablespecs = table_spec($table) ;
   my $field_rev ;
   foreach my $j (0 .. $#{$tablespecs->[0]->{specs}->{field_name}}) {
      $field_rev->{$tablespecs->[0]->{specs}->{field_name}->[$j]} = $j ; }

   my @fields_id ;
   foreach my $j ( 0 .. $#fields) {
      if (!exists $field_rev->{$fields[$j]}) {
         print STDERR "ERROR rawselect_tod(): field $fields[$j] doesn't exist\n";
	 return "ERROR" ;
      }
      $fields_id[$j] = $field_rev->{$fields[$j]} ;
   }

   if (!defined $fullfile) {
      $fullfile = $specs->{tod_dir}.'/'.$table ;
      if (!-s $fullfile) {
         print STDERR "ERROR rawselect_tod(): $fullfile not found\n" ;
         return "ERROR" ;
      }
   }

   if ($fullfile =~ /gz$/) { $fullfile = "$bins->{zcat} $fullfile |" ; }
   my @results ;
   open(TABLEF, $fullfile) ;
   while (my $line = <TABLEF>) {
      chomp $line;
      my @curf = split(/\t/, $line) ;
      foreach my $j ( 0 .. $#fields_id) {
         push @{$results[$j]}, $curf[$fields_id[$j]] ; }
   }
   close(TABLEF) ;

   return @results ;
}


=head2 timestamp()

   Function:    Returns a timestamp
   Args:        none
   Return:      $_[0] = timestamp: <4-digit YEAR><2-digit MONTH><2-digit DAY>_
                  <2-digit HOUR><2-digit MINUTE>

=cut

sub timestamp {

   my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time);

   my $f_day = '0'x(2 - length($mday)).$mday;
   my $f_mon = '0'x(2 - length(($mon + 1))).($mon + 1);
   if ($year > 100) {
      $year = $year - 100 ; }
   my $f_year = '0'x(2-length($year)).$year;

   my $f_hour = '0'x(2 - length($hour)).$hour;
   my $f_min = '0'x(2 - length($min)).$min;


   my $time_stamp = $f_year.$f_mon.$f_day.'_'.$f_hour.$f_min ;
   return $time_stamp ;

}


=head2 safe_move()

   Function:    Safely move a file to a directory (using File::Copy::move),
                  retries 14 times, and prints an error if it didnt work
   Return:	nothing
   Args:        $_[0] = source filename
                $_[1] = target directory

=cut

sub safe_move {

   my $file = shift ;
   my $dir = shift ;
   my $tries = 15 ;

   my $res = 0 ;
   if (!-s $dir) { File::Path::mkpath($dir) ; }

   while (($tries > 0) && ($res == 0 )) {
      $res = File::Copy::move($file, $dir) ;
      $tries-- ;
   }

   if (-s $file) {
      print STDERR "ERROR: couldnt move $file to $dir\t$!\n";
   }

   return ;

}


=head2 safe_copy()

   Function:    Safely copy a file to a directory (using File::Copy::copy),
                  retries 14 times, and prints an error if it didnt work
   Return:	nothing
   Args:        $_[0] = source filename
                $_[1] = target directory

=cut

sub safe_copy {

   my $file = shift ;
   my $dest = shift ;
   my $tries = 15 ;

   if (!-s $file) {
      print STDERR "ERROR: couldnt find $file\n"; return}

   my $res = 0 ;
   while (($tries > 0) && ($res == 0 )) {
      $res = File::Copy::copy($file, $dest) ;
      $tries-- ;
   }

   if (!-s $dest) {
      print STDERR "ERROR: couldnt copy $file to $dest\t$!\n"; }

   return ;

}



=head2 replace_undefs(arrayref, undef_sub)

   Name:        replace_undefs() ;
   Function:    Takes an array reference and replaces (inplace) undefined
                  values to a specified substitution value.
   Args:        $_[0] = array references
                $_[1] = undefined substitution - if the cell contains an
                  undefined value, replace with this value
   Returns:	$_ = array reference

=cut

sub replace_undefs {

   my $array = shift ;
   my $undef_sub = shift;

   foreach my $j ( 0 .. $#{$array}) {
      if (! defined $array->[$j]) { $array->[$j] = $undef_sub; }}

   return $array ;

}


=head2 rawselect_metatod()

   Title:       rawselect_metatod()
   Function:    selects specified fields from a table-on-disk.
      Note: WHERE clause does not work, this command just recognizes
            the field names and returns the appropriate columns.
   Args:        $_[0] = filename
                $_[1] = SELECT sql command
   Returns:     array of query results

=cut


sub rawselect_metatod {

   my $filename = shift ;
   my $select_sql = shift ;

   my $bins = modtie::locate_binaries() ;
   my $modtie_specs = modtie::set_modtie_specs() ;
   my $specs = $modtie_specs->{pibase_specs} ;

   my $proto = $filename;
   $proto =~ s/.*\/// ;
   $proto =~ s/_[0-9]+.*// ;
   $proto =~ s/\.[0-9]+.*// ;
   $proto .= '_prototype' ;

   if ($proto eq 'secstrx_prototype') {
      $proto = 'bdp_secstrx_prototype' ; }

   if ($proto eq 'secstrx_basic_prototype') {
      $proto = 'bdp_secstrx_basic_prototype' ; }

   my $tablespecs = table_spec($proto) ;

   my $fields = $select_sql ;
   $fields =~ s/(FROM|from).*$// ;
   $fields =~ s/^(SELECT|select)// ;
   $fields =~ s/ //g ;
   my @fields = split(/\,/, $fields) ;

   my $field_rev ;
   foreach my $j (0 .. $#{$tablespecs->[0]->{specs}->{field_name}}) {
      $field_rev->{$tablespecs->[0]->{specs}->{field_name}->[$j]} = $j ; }

   my @fields_id ;
   foreach my $j ( 0 .. $#fields) {
      if (!exists $field_rev->{$fields[$j]}) {
         print STDERR "rawselect_metatod() error: field $fields[$j] doesn't exist\n" ;
	 return "ERROR" ;
      }
      $fields_id[$j] = $field_rev->{$fields[$j]} ;
   }

   my @results ;
   if (!-s $filename) {
      print STDERR "ERROR rawselect_metatod(): $filename not found\n" ;
      return "ERROR" ;
   }

   if ($filename =~ /gz$/) {
      $filename = "$bins->{zcat} $filename |" ; }

   open(TABLEF, $filename) ;
   while (my $line = <TABLEF>) {
      chomp $line;
      my @curf = split(/\t/, $line) ;
      foreach my $j ( 0 .. $#fields_id) {
         push @{$results[$j]}, $curf[$fields_id[$j]] ; }
   }
   close(TABLEF) ;

   return @results ;

}


=head2 mysql_hashload(dbh, query)

   Name:        mysql_hashload()
   Function:    Processes a 2 column query and returns a hash pointing from
                  column1 values to column2 values (use for 1:1 relationships)
   Args:        $_[0] = dbh - DBI database handle
                $_[1] = SQL SELECT command
   Returns:     $a - hashref
                $a->{col1} = col2

=cut

sub mysql_hashload {

   my ($dbh, $query) = @_ ;

   $query = qq{$query} ;
   my $sth = $dbh->prepare($query) ;
   $sth->execute() ;

   my $ans ;
   while ( my @currow = $sth->fetchrow() ) {
      $ans->{$currow[0]} = $currow[1] ; }

   return $ans ;

}


=head2 residue_int(resno)

   Title:       residue_int()
   Function:    Seperates the integer and insertion code of a residue number
      NOTE: assumes insertion code is always alphanumeric
   Args:        residue number (5 character - number and insertion code)
   Returns:     $_[0] integer portion of the residue number
                $_[1] insertion code of the residue number

=cut

sub residue_int {

   my $res_no = shift ;

   if (!defined $res_no) {
      die "res_no was undefined\n"; }

# Extract the integer and insertion code from the residue number.

   my ($res_no_int, $ins_code) = ($res_no =~ /(-?[0-9]+)(.?)/) ;

# If the insertion code is undefined, set to ''.

   if (!defined $ins_code) {
      $ins_code = ''; }

   return ($res_no_int, $ins_code) ;

}



1 ;
