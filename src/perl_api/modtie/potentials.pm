=head1 NAME

potentials.pm - MODTIE routines to build and use statistical potentials

=head1 DESCRIPTION

The potentials.pm perl library contains subroutines to calculate, benchmark,
and use statistical potentials to score putative protein interfaces.

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


package modtie::potentials ;
use strict;
use Exporter ;
our @ISA = qw/Exporter/ ;
our @EXPORT_OK = qw/maxval_arrref buildpotential_count buildpotential_postcalc/;

use modtie ;
use modtie::pibase ;

sub buildpotential_count {

   my $in = shift ;
   my $specs ;
   if (exists $in->{specs}) { $specs = $in->{specs} ; }
   else { $specs = modtie::set_modtie_specs() ; }

   my $aa3to1 = $specs->{gconst}->{aa3to1};

# get db specs, binary locations and preprocess subsets_residues_tables

   my $dbspecs = $specs->{pibase_specs} ;
   my $pibase = $dbspecs->{db} ;

   my $host = hostname() ;

# Load subsets_files locations into memory

   my $meta_tables;
   my $time = modtie::pibase::timestamp() ;
   print STDERR "reading in subsets_residues_tables: " ;
   {
      my @tod_res = modtie::pibase::rawselect_tod(
         "SELECT subset_id, file_path FROM subsets_files") ;
      foreach my $j ( 0 .. $#{$tod_res[0]} ) {
         $meta_tables->{subsets_files}->{$tod_res[0]->[$j]} =
            $tod_res[1]->[$j];
      }
   }
   print STDERR "X\n" ;

# Load subsets_residues_table into memory indexed with bdp_id
   print STDERR "reading in subsets_residues_tables: " ;
   {
      my @tod_res = modtie::pibase::rawselect_tod(
         "SELECT bdp_id, source_file FROM subsets_residues_tables") ;
      foreach my $j ( 0 .. $#{$tod_res[0]}) {
         $meta_tables->{subsets_residues}->{$tod_res[0]->[$j]} =
            $tod_res[1]->[$j];
      }
   }
   print STDERR "X\n" ;


# Load bdppaths into memory indexed with bdp_id

   my ($bdp2path) = modtie::pibase::todload_bdp_ids('bdp_id_2_path') ;

# Load the entire subsets table into memory and index with bdp_id

   print STDERR "reading in subsets:" ;
   my $allsubs ;
   my $allsubs_point ;
   my $sid2bdp ;
   {
      my @tod_res = modtie::pibase::rawselect_tod(
         "SELECT bdp_id, subset_id, description, subset_source_id, class ".
         "FROM subsets") ;

      foreach my $j ( 0 .. $#{$tod_res[0]}) {
         if (defined $tod_res[0]->[$j]) { # if bdp_id IS NOT NULL
            push @{$allsubs->{bdp_id}}, $tod_res[0]->[$j] ;
            push @{$allsubs->{subset_id}}, $tod_res[1]->[$j] ;
            push @{$allsubs->{description}}, $tod_res[2]->[$j] ;
            push @{$allsubs->{source_id}}, $tod_res[3]->[$j] ;
            push @{$allsubs->{class}}, $tod_res[4]->[$j] ;
         }
      }

      foreach my $j (0 .. $#{$allsubs->{bdp_id}}) {
         $sid2bdp->{$allsubs->{subset_id}->[$j]} = $allsubs->{bdp_id}->[$j] ;
         push @{$allsubs_point->{$allsubs->{bdp_id}->[$j]}}, $j ;
      }
   }
   print STDERR "X\n" ;

   my $tempdir = tempdir(CLEANUP=>1) ;
   chdir $tempdir ;

# Iterate through the bdp_ids.

   my $counts;
   my $types;
   $types->{sc}= ['mm', 'ms', 'ss', 'all'] ;
   $types->{rad}= [4,6,8] ;
   $types->{inter}= ['inter', 'intra'] ;

   foreach my $sc_type (@{$types->{sc}}) {
      foreach my $radius (@{$types->{rad}}) {
         foreach my $restype1 (keys %{$aa3to1}) {
            $counts->{$sc_type}->{$radius}->{n_P} = 0 ;
            foreach my $restype2 (keys %{$aa3to1}) {
               $counts->{$sc_type}->{$radius}->{maxafp_ij}->{$restype1}->{$restype2} = 0 ;
               $counts->{$sc_type}->{$radius}->{nij_P}->{$restype1}->{$restype2} = 0 ;
               $counts->{$sc_type}->{$radius}->{dnij_P__n_P}->{$restype1}->{$restype2} = 0 ;
            }
         }
      }
   }


   foreach my $entrynum (0 .. $#{$in->{entries}}) {
      print STDERR "now on $entrynum: " ;

      my $params ;
      my $estring = '' ;
      foreach my $key (sort keys %{$in->{entries}->[$entrynum]}) {
         $estring .= ' '.$in->{entries}->[$entrynum]->{$key} ;
         $params->{$key} = $in->{entries}->[$entrynum]->{$key} ; }
      print STDERR " $estring\n" ;

      $params->{counts} = $counts ;
      $params->{sid2bdp} = $sid2bdp;
      $params->{meta_tables} = $meta_tables ;
      $params->{allsubs_point} = $allsubs_point ;
      $params->{allsubs} = $allsubs ;
      $params->{R} = $types->{rad} ;
      $params->{R_max} = maxval_arrref($types->{rad}) ;
      $params->{sid2type} = "SCOP" ;

      _buildpotentialcount_countup_asshole($params) ;
   }

#for each type:
# {4 | 6 | 8 } {inter | intra | all} | {mm | ms | ss | all }

# global counts:
#       g_ij_num
#       g_ij_denom

# local counts:
#       dn_ij
#       n_p
#       n__p_i
#       n_j

   my @outvals = () ;
   foreach my $sc_type (@{$types->{sc}}) {
      foreach my $radius (@{$types->{rad}}) {
         @outvals = ($sc_type, $radius, 'n_P', $counts->{$sc_type}->{$radius}->{n_P}) ;
         print join("\t", @outvals)."\n" ;
         foreach my $restype1 (sort keys %{$aa3to1}) {
            foreach my $restype2 (sort keys %{$aa3to1}) {
               @outvals = ($sc_type, $radius,'maxafp_ij', $restype1, $restype2,
                  sprintf("%.3f",
                     $counts->{$sc_type}->{$radius}->{maxafp_ij}->{$restype1}->{$restype2})) ;
               print join("\t", @outvals)."\n" ;

               @outvals = ($sc_type, $radius,'nij_P', $restype1, $restype2,
                  $counts->{$sc_type}->{$radius}->{nij_P}->{$restype1}->{$restype2}) ;
               print join("\t", @outvals)."\n" ;

               @outvals = ($sc_type, $radius,'dnij_P__n_P', $restype1, $restype2,
                  sprintf("%.3f",
                     $counts->{$sc_type}->{$radius}->{dnij_P__n_P}->{$restype1}->{$restype2}));
               print join("\t", @outvals)."\n" ;
            }
         }
      }
   }
   chdir ;

}


sub _buildpotentialcount_countup_asshole {

   my $in = shift ; my $specs ;
   if (exists $in->{specs}) { $specs = $in->{specs} ; }
   else { $specs = modtie::set_modtie_specs() ; }

   my $aa3to1 = $specs->{gconst}->{aa3to1} ; 
   my $mcatoms = $specs->{gconst}->{mcatoms};
   my $numaaatoms = $specs->{gconst}->{numaaatoms};

   my $params = shift ;
   my $in = $params ;
   my $meta_tables = $in->{meta_tables} ;
   my $allsubs_point = $in->{allsubs_point} ;
   my $allsubs = $in->{allsubs} ;

# catpdb into kdcontacts

   my ($t_pdbfh, $t_pdbfn) = tempfile("temppdb.XXXXX", SUFFIX => ".pdb") ;
   close($t_pdbfh) ;


   my $sidtype ;
   my $strx_path;
   my $bdp_id ;
   if (exists $params->{sid_1}) {
      $sidtype = $params->{sid2type} ;
      my $tcom = "cat $meta_tables->{subsets_files}->{$params->{sid_1}} $meta_tables->{subsets_files}->{$params->{sid_2}} > $t_pdbfn" ;
      system($tcom) ;
      $strx_path = $t_pdbfn ;
      $bdp_id = $params->{sid2bdp}->{$params->{sid_1}} ;
   } else {
      unlink $t_pdbfn ;
      $sidtype = $params->{sid2type} ;
      $strx_path = $meta_tables->{subsets_files}->{$params->{sid}} ;
      $bdp_id = $params->{sid2bdp}->{$params->{sid}} ;
   }


   my $tables ;

# Get the file name that holds subsets_residues

   $tables->{subsets_residues}->{sourcefile} =
      $meta_tables->{subsets_residues}->{$bdp_id} ;

# mysql(SELECT subset_id, descsription, subset_source_id, class FROM subsets WHERE bdp_id = <bdp_id>)

   my ( $subset_id, $subset_description, $subset_source_id, $subset_class );
   {
      my @tind = @{$allsubs_point->{$bdp_id}} ;
      push @{$subset_id}, @{$allsubs->{subset_id}}[@tind] ;
      push @{$subset_description}, @{$allsubs->{description}}[@tind] ;
      push @{$subset_source_id}, @{$allsubs->{source_id}}[@tind] ;
      push @{$subset_class}, @{$allsubs->{class}}[@tind] ;
   }

# Iterate through subsets and creating hashes pointing from subset_id to subset_class and subset_id to subset_source_id

   my $subset_id_2_index ;
   my @subsets_source_id ;
   my @subsets_id ;
   my @subsets_class ;
   foreach my $j ( 0 .. $#{$subset_id}) {
      $subsets_id[$j] = $subset_id->[$j] ;
      $subsets_class[$j] = $subset_class->[$j] ;
      $subsets_source_id[$j] = $subset_source_id->[$j] ;
      $subset_id_2_index->{$subset_id->[$j]} = $j ;
   }

# Read in residue subset assignments from subsets_residues

#      print STDERR "   reading in subsets_residues" ;

   my ( $subres_chain_no, $subres_chain_id, $subres_resno,
           $subres_resno_int, $subres_subset_id ) =
      modtie::pibase::rawselect_metatod(
         $tables->{subsets_residues}->{sourcefile},
         "SELECT chain_no, chain_id, ".
         "resno, resno_int, subset_id FROM subsets_residues" ) ;

# If the residue subset assignments is empty, display to STDERR and abort this bdp_id.

   if ( $#{$subres_subset_id} < 0 ) {
      print STDERR "\nSKIPPED bdp $bdp_id: no residue subset assignments found\n";
      print STDERR "X\n" ;
      next;
   }

# Load subsets residue assignments into a hash pointing from residue to subsets

   my $subset_assign ;
   foreach my $j (0 .. $#{$subres_chain_no}) {

    if ($subres_subset_id->[$j] =~ /$sidtype/) {
      if (!defined $subres_chain_id->[$j]) {
         $subres_chain_id->[$j] = ' '; }

# BUG in tables holding subset residue number: blank space issues.

      $subres_resno->[$j] =~ s/ //g ;
#      print STDERR "sid assign: $subset_id_2_index->{$subres_subset_id->[$j]} for $subres_resno->[$j] on $subres_chain_id->[$j]\n" ;

      my $sig = $subres_resno->[$j]."\n".$subres_chain_id->[$j] ;

# assuming only one domain participation per residue (within a system - eg SCOP)
      $subset_assign->{$sig} = $subset_id_2_index->{$subres_subset_id->[$j]} ;
    }
   }


# Laod bdp_chains info into a reverse hash  pointing from chain_id to chain_no ;

# mysql(SELECT chain_id_1, resno_1, resno_1_int, resna_1, chain_id_2, resno_2, resno_2_int, resna_2, count(distance) from interatomic_contacts_<n> WHERE bdp_id = <bdp_id>) group by chain_id_1, resno_1, chain_id_2, resno_2

   my @cont_fields = ('resna1', 'resno1', 'inscode1', 'chain_id1', 'atomno1', 'atomna1', 'resna2', 'resno2', 'inscode2', 'chain_id2', 'atomno2', 'atomna2', 'dist' ) ;
   my $kdfield2no ;
   foreach my $j ( 0 .. $#cont_fields) { $kdfield2no->{$cont_fields[$j]} = $j;}

# Iterate through the contacts, build all res - res contact info vectors

   my $kdcont_out  = _buildpotentialcount_calc_res_pairs({
         radius => $params->{R_max},
         bdp_path => $strx_path
   }) ;

   if (exists $kdcont_out->{error_fl})  {
      print STDERR "ERROR (bdp_id $bdp_id): $kdcont_out->{error_fl}\n" ;
      next;
   }
   my $kdcontacts_out = $kdcont_out->{contacts_fh} ;

   my $respairs ;
   my $respair_atoms ;
   my $resnames;

   while (my $line = <$kdcontacts_out>) {
      chomp $line;
      if ($line =~ /^#/) {next;}

      my ($resna_1, $resno_1, $inscode1, $chain_id_1, undef, $atomna_1,
             $resna_2, $resno_2, $inscode2, $chain_id_2, undef, $atomna_2,
             $dist) = split(/\t/, $line) ;

      if ($resna_1 eq 'HSD') {$resna_1 = 'HIS'; }
      if ($resna_2 eq 'HSD') {$resna_2 = 'HIS'; }
      if (!exists $specs->{gconst}->{aa3to1}->{$resna_1}||
          !exists $specs->{gconst}->{aa3to1}->{$resna_2}) {next;}

      if (    (substr($atomna_1, 1, 1) eq 'H') ||
              (substr($atomna_2, 1, 1) eq 'H') ||
              (substr($atomna_1, 1, 1) eq 'Q') ||
              (substr($atomna_2, 1, 1) eq 'Q') ) {next;}

      if ((!defined $chain_id_1) || ($chain_id_1 eq '')) {$chain_id_1 = ' ';}
      if ((!defined $chain_id_2) || ($chain_id_2 eq '')) {$chain_id_2 = ' ';}

      $inscode1 =~ s/ //g ; $inscode2 =~ s/ //g ;
      $resno_1 .= $inscode1 ; $resno_2 .= $inscode2 ;
      my $sig1 = $resno_1."\n".$chain_id_1 ;
      my $sig2 = $resno_2."\n".$chain_id_2 ;

# take care of double counting and intra-residue contacts
      if ($sig1 le $sig2) {next;}

      my $ressig = $sig1."\n".$sig2 ;

      my $sc_type1 = 's'; my $sc_type2 = 's';
      if (exists $mcatoms->{$atomna_1}) { $sc_type1 = 'm' ; }
      if (exists $mcatoms->{$atomna_2}) { $sc_type2 = 'm' ; }
      my $sc_type = $sc_type1.$sc_type2 ;

      $resnames->{$sig1} = $resna_1 ; $resnames->{$sig2} = $resna_2 ;

      if (!exists $respair_atoms->{$ressig}->{$sc_type}->[0]->{$atomna_1}->{min_dist} ||
          ($dist<=$respair_atoms->{$ressig}->{$sc_type}->[0]->{$atomna_1}->{min_dist})) {
                  $respair_atoms->{$ressig}->{$sc_type}->[0]->{$atomna_1}->{min_dist} = $dist; }

      if (!exists $respair_atoms->{$ressig}->{$sc_type}->[1]->{$atomna_2}->{min_dist} ||
          ($dist<=$respair_atoms->{$ressig}->{$sc_type}->[1]->{$atomna_2}->{min_dist})) {
                  $respair_atoms->{$ressig}->{$sc_type}->[1]->{$atomna_2}->{min_dist} = $dist; }

      if (!exists $respair_atoms->{$ressig}->{all}->[0]->{$atomna_1}->{min_dist} ||
          ($dist<=$respair_atoms->{$ressig}->{all}->[0]->{$atomna_1}->{min_dist})) {
                  $respair_atoms->{$ressig}->{all}->[0]->{$atomna_1}->{min_dist} = $dist; }

      if (!exists $respair_atoms->{$ressig}->{all}->[1]->{$atomna_2}->{min_dist} ||
          ($dist<=$respair_atoms->{$ressig}->{all}->[1]->{$atomna_2}->{min_dist})) {
                  $respair_atoms->{$ressig}->{all}->[1]->{$atomna_2}->{min_dist} = $dist; }

      if ((!exists $respairs->{$ressig}->{all}->{min_dist}) ||
         ($dist <= $respairs->{$ressig}->{all}->{min_dist})) {
                   $respairs->{$ressig}->{all}->{min_dist} = $dist ; }

      if ((!exists $respairs->{$ressig}->{$sc_type}->{min_dist}) ||
         ($dist <= $respairs->{$ressig}->{$sc_type}->{min_dist})) {
                   $respairs->{$ressig}->{$sc_type}->{min_dist} = $dist ; }

#CHANGEHERE - track what atoms are making the contacts for each sc_type; so can count later

   }
   close $kdcontacts_out ;
#   print STDERR "check out $kdcontacts_out\n" ;

# first assign inter/intra flags, then assing 4 / 6 / 8 cuts

# Iterate through the residue pairs

# For each new interacting atom (residue? fpd030714_0531), iterate through BDP subsets, and find ALL correspding domains, then output for import into subsets_residue.


   my $counts ;

   if (exists $params->{sid_2}) {
# inter counts

      foreach my $respair (keys %{$respairs}) {
         foreach my $sc_type (keys %{$respairs->{$respair}}) {
            my ($resno_1, $chain_id_1, $resno_2, $chain_id_2) =
               split(/\n/, $respair) ;

            my $sig1 = $resno_1."\n".$chain_id_1 ;
            my $sig2 = $resno_2."\n".$chain_id_2 ;
            my $resna1 = $resnames->{$sig1} ;
            my $resna2 = $resnames->{$sig2} ;


# Go looking for subset membership of res 1 and 2.  

            if (!exists $subset_assign->{$sig1} ||
                !exists $subset_assign->{$sig2}) { next; }
            my $sid1 = $subset_assign->{$sig1} ;
            my $sid2 = $subset_assign->{$sig2} ;

            if ($sid1 eq $sid2) {next;}
#            print STDERR " counting up on a $sid1 -- $sid2: $resna1 $resna2 on a $sc_type\n" ;
            
            foreach my $radius (@{$params->{R}}) {
               if ($respairs->{$respair}->{$sc_type}->{min_dist} <= $radius) {

#COUNTCHANGE to fractional increment instead of ++ ;
# afp = atom fractional participation = (interacting_atoms->{sc_type} / total_atoms->{sc_type})
# fraction =  min( afp_res1 , afp_res2)

                  my $sc1 = 'all'; my $sc2 = 'all' ;
                  if ($sc_type ne 'all') {
                     $sc1 = substr($sc_type, 0, 1) ;
                     $sc2 = substr($sc_type, 1, 1) ;
                  }

                  my $t_intatoms1 = 0 ;
                  foreach my $atom (keys %{$respair_atoms->{$respair}->{$sc_type}->[0]}) {
                     if ($respair_atoms->{$respair}->{$sc_type}->[0]->{$atom}->{min_dist} <= $radius) {
                        $t_intatoms1++ ; } }

                  my $t_intatoms2 = 0 ;
                  foreach my $atom (keys %{$respair_atoms->{$respair}->{$sc_type}->[1]}) {
                     if ($respair_atoms->{$respair}->{$sc_type}->[1]->{$atom}->{min_dist} <= $radius) {
                        $t_intatoms2++ ; } }


                  my $afp1 = $t_intatoms1 / $numaaatoms->{$sc1}->{$resna1} ;
                  my $afp2 = $t_intatoms2 / $numaaatoms->{$sc2}->{$resna2} ;

                  my $minafp = 0;
                  if ($afp1 < $afp2) {$minafp = $afp1;} else {$minafp = $afp2;}
                  if ($minafp > 1) {
                     my $t = "WARNING: minafp=$minafp, $sc_type on $respair" ;
                     $t =~ s/\n/:/g ;
                     print STDERR $t."\n" ;
                  }

                  my $rev_fl ;
                  my $rsc_type = $sc_type ;
                  if ($sc_type eq 'sm') { $rev_fl = 1; $rsc_type = 'ms';}


#COUNTCHANGE: 1,2,3
                  if ($rev_fl) {
                     if (!exists $counts->{$rsc_type}->{$radius}->{maxafp_ij}->{$resna2}->{$resna1} ||
                         ($minafp > $counts->{$rsc_type}->{$radius}->{maxafp_ij}->{$resna2}->{$resna1})) {
                        $counts->{$rsc_type}->{$radius}->{maxafp_ij}->{$resna2}->{$resna1} = $minafp ; }

                     $counts->{$rsc_type}->{$radius}->{dnij_P}->{$resna2}->{$resna1} += $minafp ;
                  } else {

                     if (!exists $counts->{$rsc_type}->{$radius}->{maxafp_ij}->{$resna1}->{$resna2} ||
                         ($minafp > $counts->{$rsc_type}->{$radius}->{maxafp_ij}->{$resna1}->{$resna2})) {
                        $counts->{$rsc_type}->{$radius}->{maxafp_ij}->{$resna1}->{$resna2} = $minafp ; }

                     $counts->{$rsc_type}->{$radius}->{dnij_P}->{$resna1}->{$resna2} += $minafp ;

                     if (($resna1 ne $resna2) && ($rsc_type ne 'ms')) {

                        if (!exists $counts->{$rsc_type}->{$radius}->{maxafp_ij}->{$resna2}->{$resna1} ||
                            ($minafp > $counts->{$rsc_type}->{$radius}->{maxafp_ij}->{$resna2}->{$resna1})) {
                           $counts->{$rsc_type}->{$radius}->{maxafp_ij}->{$resna2}->{$resna1} = $minafp ; }

                        $counts->{$rsc_type}->{$radius}->{dnij_P}->{$resna2}->{$resna1} += $minafp; }
                  }

                  $counts->{$rsc_type}->{$radius}->{ni_P}->{$sid1}->{$resna1}->{$sig1}++ ;
                  $counts->{$rsc_type}->{$radius}->{ni_P}->{$sid2}->{$resna2}->{$sig2}++ ;

                  $counts->{$rsc_type}->{$radius}->{n_P}->{$sid1}->{$sig1}++ ;
                  $counts->{$rsc_type}->{$radius}->{n_P}->{$sid2}->{$sig2}++ ;
               }
            }
         }
      }

      my $truecounts = $params->{counts}; #hashref: maniuplates global counts
      foreach my $sc_type (keys %{$counts}) {
         foreach my $radius (keys %{$counts->{$sc_type}}) {

            my @sids = keys %{$counts->{$sc_type}->{$radius}->{n_P}} ;
            my $numall1= keys %{$counts->{$sc_type}->{$radius}->{n_P}->{$sids[0]}};
            my $numall2= keys %{$counts->{$sc_type}->{$radius}->{n_P}->{$sids[1]}};
            my $numall12 = $numall1 + $numall2 ;
            $truecounts->{$sc_type}->{$radius}->{n_P} += $numall1 + $numall2 ;

            foreach my $restype1
               (keys %{$counts->{$sc_type}->{$radius}->{dnij_P}}) {

               if (($sc_type =~ /^s/) && $restype1 eq 'GLY') {next;}

               foreach my $restype2 
                  (keys %{$counts->{$sc_type}->{$radius}->{dnij_P}->{$restype1}}) {

                  if (($sc_type =~ /s$/) && $restype2 eq 'GLY') {next;}

                  if ($counts->{$sc_type}->{$radius}->{maxafp_ij}->{$restype1}->{$restype2} >
                  $truecounts->{$sc_type}->{$radius}->{maxafp_ij}->{$restype1}->{$restype2}) {
                     $truecounts->{$sc_type}->{$radius}->{maxafp_ij}->{$restype1}->{$restype2} =
                        $counts->{$sc_type}->{$radius}->{maxafp_ij}->{$restype1}->{$restype2} ; }

                  $truecounts->{$sc_type}->{$radius}->{dnij_P__n_P}->{$restype1}->{$restype2} +=
                     ($counts->{$sc_type}->{$radius}->{dnij_P}->{$restype1}->{$restype2} * $numall12) ;
               }

            }


            foreach my $restype1 (keys %{$aa3to1}) {
               if (($sc_type =~ /^s/) && $restype1 eq 'GLY') {next;}

               my $t_num1_a = keys %{$counts->{$sc_type}->{$radius}->{ni_P}->{$sids[0]}->{$restype1}} ;
               my $t_num1_b = keys %{$counts->{$sc_type}->{$radius}->{ni_P}->{$sids[1]}->{$restype1}} ;

               foreach my $restype2 (keys %{$aa3to1}) {

                  if (($sc_type =~ /s$/) && $restype2 eq 'GLY') {next;}

                  my $t_num2_a = 0 ;
                  if (exists $counts->{$sc_type}->{$radius}->{ni_P}->{$sids[0]}->{$restype2}) {
                     $t_num2_a = keys %{$counts->{$sc_type}->{$radius}->{ni_P}->{$sids[0]}->{$restype2}} ; }

                  my $t_num2_b = 0 ;
                  if (exists $counts->{$sc_type}->{$radius}->{ni_P}->{$sids[1]}->{$restype2}) {
                     $t_num2_b = keys %{$counts->{$sc_type}->{$radius}->{ni_P}->{$sids[1]}->{$restype2}} ;}

                  $truecounts->{$sc_type}->{$radius}->{nij_P}->{$restype1}->{$restype2} +=
                     ($t_num1_a * $t_num2_b)  ;

# fix ms diagonal bug
                  if ((($restype1 eq $restype2) && $sc_type eq 'ms') ||
                      ($restype1 ne $restype2) ) {
                     $truecounts->{$sc_type}->{$radius}->{nij_P}->{$restype1}->{$restype2} +=
                        ($t_num1_b * $t_num2_a) ; }
               }
            }
         }
      }

      unlink $t_pdbfn ;

   } else {

# intra counts

      foreach my $respair (keys %{$respairs}) {
         foreach my $sc_type (keys %{$respairs->{$respair}}) {
            my ($resno_1, $chain_id_1, $resno_2, $chain_id_2) =
               split(/\n/, $respair) ;

            my $sig1 = $resno_1."\n".$chain_id_1 ;
            my $sig2 = $resno_2."\n".$chain_id_2 ;
            my $resna1 = $resnames->{$sig1} ;
            my $resna2 = $resnames->{$sig2} ;

            foreach my $radius (@{$params->{R}}) {
               if ($respairs->{$respair}->{$sc_type}->{min_dist} <= $radius) {

#COUNTCHANGE to fractional increment instead of ++ ;
# afp = atom fractional participation = (interacting_atoms->{sc_type} / total_atoms->{sc_type})
# fraction =  min( afp_res1 , afp_res2)

                  my $sc1 = 'all'; my $sc2 = 'all' ;
                  if ($sc_type ne 'all') {
                     $sc1 = substr($sc_type, 0, 1) ;
                     $sc2 = substr($sc_type, 1, 1) ;
                  }

                  my $t_intatoms1 = 0 ;
                  foreach my $atom (keys %{$respair_atoms->{$respair}->{$sc_type}->[0]}) {
                     if ($respair_atoms->{$respair}->{$sc_type}->[0]->{$atom}->{min_dist} <= $radius) {
                        $t_intatoms1++ ; } }

                  my $t_intatoms2 = 0 ;
                  foreach my $atom (keys %{$respair_atoms->{$respair}->{$sc_type}->[1]}) {
                     if ($respair_atoms->{$respair}->{$sc_type}->[1]->{$atom}->{min_dist} <= $radius) {
                        $t_intatoms2++ ; } }

                  my $afp1 = $t_intatoms1 / $numaaatoms->{$sc1}->{$resna1} ;
                  my $afp2 = $t_intatoms2 / $numaaatoms->{$sc2}->{$resna2} ;
                  my $minafp = 0;
                  if ($afp1 < $afp2) {$minafp = $afp1;} else {$minafp = $afp2;}


                  my $rev_fl ;
                  my $rsc_type = $sc_type ;
                  if ($sc_type eq 'sm') { $rev_fl = 1; $rsc_type = 'ms';}

#COUNTCHANGE: 1,2,3
                  if ($rev_fl) {

                     if (!exists $counts->{$rsc_type}->{$radius}->{maxafp_ij}->{$resna2}->{$resna1} ||
                         $minafp > $counts->{$rsc_type}->{$radius}->{maxafp_ij}->{$resna2}->{$resna1}) {
                        $counts->{$rsc_type}->{$radius}->{maxafp_ij}->{$resna2}->{$resna1} = $minafp ; }

                     $counts->{$rsc_type}->{$radius}->{dnij_P}->{$resna2}->{$resna1} += $minafp ;

                  } else {

                     if (!exists $counts->{$rsc_type}->{$radius}->{maxafp_ij}->{$resna1}->{$resna2} ||
                         $minafp > $counts->{$rsc_type}->{$radius}->{maxafp_ij}->{$resna1}->{$resna2}) {
                        $counts->{$rsc_type}->{$radius}->{maxafp_ij}->{$resna1}->{$resna2} = $minafp ; }

                     $counts->{$rsc_type}->{$radius}->{dnij_P}->{$resna1}->{$resna2} += $minafp ;

                     if (($resna1 ne $resna2) && $sc_type ne 'ms') {

                        if (!exists $counts->{$rsc_type}->{$radius}->{maxafp_ij}->{$resna2}->{$resna1} ||
                            $minafp > $counts->{$rsc_type}->{$radius}->{maxafp_ij}->{$resna2}->{$resna1}) {
                           $counts->{$rsc_type}->{$radius}->{maxafp_ij}->{$resna2}->{$resna1} = $minafp ; }

                        $counts->{$rsc_type}->{$radius}->{dnij_P}->{$resna2}->{$resna1} += $minafp; }

                  }


                  $counts->{$rsc_type}->{$radius}->{ni_P}->{$resna1}->{$sig1}++ ;
                  $counts->{$rsc_type}->{$radius}->{ni_P}->{$resna2}->{$sig2}++ ;

                  $counts->{$rsc_type}->{$radius}->{n_P}->{$sig1}++ ;
                  $counts->{$rsc_type}->{$radius}->{n_P}->{$sig2}++ ;

               }
            }
         }
      }

      my $truecounts = $params->{counts}; #hashref: maniuplates global counts
      foreach my $sc_type (keys %{$counts}) {
         foreach my $radius (keys %{$counts->{$sc_type}}) {

            my $numall = keys %{$counts->{$sc_type}->{$radius}->{n_P}} ;
            $truecounts->{$sc_type}->{$radius}->{n_P} += $numall ;

            foreach my $restype1 (keys %{$counts->{$sc_type}->{$radius}->{dnij_P}}) {
               if (($sc_type =~ /^s/) && $restype1 eq 'GLY') {next;}

               foreach my $restype2 (keys %{$counts->{$sc_type}->{$radius}->{dnij_P}->{$restype1}}) {

                  if (($sc_type =~ /s$/) && $restype2 eq 'GLY') {next;}

#COUNTCHANGE
                  if ($counts->{$sc_type}->{$radius}->{maxafp_ij}->{$restype1}->{$restype2} >
                  $truecounts->{$sc_type}->{$radius}->{maxafp_ij}->{$restype1}->{$restype2}) {
                     $truecounts->{$sc_type}->{$radius}->{maxafp_ij}->{$restype1}->{$restype2} =
                        $counts->{$sc_type}->{$radius}->{maxafp_ij}->{$restype1}->{$restype2} ; }

                  $truecounts->{$sc_type}->{$radius}->{dnij_P__n_P}->{$restype1}->{$restype2} +=
                    $counts->{$sc_type}->{$radius}->{dnij_P}->{$restype1}->{$restype2} * $numall ;
               }
            }


            foreach my $restype1 (keys %{$counts->{$sc_type}->{$radius}->{ni_P}}) {
               if (($sc_type =~ /^s/) && $restype1 eq 'GLY') {next;}
               my $t_num1 = keys %{$counts->{$sc_type}->{$radius}->{ni_P}->{$restype1}} ;
               $truecounts->{$sc_type}->{$radius}->{ni_P}->{$restype1} += $t_num1 ;

               foreach my $restype2 (keys %{$counts->{$sc_type}->{$radius}->{ni_P}}) {
                  if (($sc_type =~ /s$/) && $restype2 eq 'GLY') {next;}

                  if ($restype1 eq $restype2) {

#fix ms diagonal bug
                     if ($sc_type eq 'ms') {
                        $truecounts->{$sc_type}->{$radius}->{nij_P}->{$restype1}->{$restype2} +=
                           ($t_num1) * ($t_num1 - 1) ;
                     } else {
                        $truecounts->{$sc_type}->{$radius}->{nij_P}->{$restype1}->{$restype2} +=
                           ($t_num1) * ($t_num1 - 1) / 2 ;
                     }

                  } else {
                     my $t_num2 = keys %{$counts->{$sc_type}->{$radius}->{ni_P}->{$restype2}} ;
                     $truecounts->{$sc_type}->{$radius}->{nij_P}->{$restype1}->{$restype2} +=
                        $t_num1 * $t_num2 ;
                  }
               }
            }
         }
      }
   }
   
}


sub _buildpotentialcount_calc_res_pairs {
   my $params = shift ;

   my $kdcontacts_radius = $params->{radius} || "6.6" ;

   my $binaries = modtie::locate_binaries() ;
   my $kdcontacts_bin = $binaries->{'kdcontacts'}." $kdcontacts_radius" ;
   my $altloc_check = $binaries->{'altloc_check'} ;
   my $altloc_filter = $binaries->{'altloc_filter'} ;

   my $bdp_file_path = $params->{bdp_path} ;

   my $host = hostname() ;

# Check if the pdb file contains altloc identifiers. If so, first filter with altloc_filter.pl

   my $altloc_fl = `$altloc_check < $bdp_file_path` ; chomp $altloc_fl ;
   my $tcom ;
   if ($altloc_fl) {
      $tcom = "$altloc_filter $bdp_file_path" ;
   } else {
      $tcom = "cat $bdp_file_path" ;
   }

# system(cat <bdp_file_path> | kdcontacts_bin 2><kdcontacts_err> ><kdcontacts_out>)

   my $out_fh ;
   $tcom .= " | $kdcontacts_bin | " ;
   open ($out_fh, $tcom) ;

#   $tcom .= " | $kdcontacts_bin 2>$kdcontacts_err >$kdcontacts_out" ;
#   system($tcom) ;

# send back filehandle to kdcontacts output

# read in kdcontacts output into an array

   my $results = { contacts_fh => $out_fh };

   return $results ;
}


sub buildpotential_postcalc {

#original form: buildpotential.count.pl

   my $in = shift ;
   my $specs ;
   if (exists $in->{specs}) { $specs = $in->{specs} ; }
   else { $specs = modtie::set_modtie_specs() ; }

   my $aa3to1 = $specs->{gconst}->{aa3to1} ;

# get db specs, binary locations and preprocess subsets_residues_tables

   my $types;
   $types->{sc}= ['mm', 'ms', 'ss', 'all'] ;
   $types->{rad}= [4,6,8] ;
   $types->{inter}= ['inter', 'intra'] ;

   my $counts = $in->{counts};
# read in (and sum up) g_ij

   my ($g_ij, $G_ij, $w_ij) ;
   my $G_ij_denom ;
   foreach my $sc_type (@{$types->{sc}}) {
      foreach my $radius (@{$types->{rad}}) {
         $G_ij_denom->{$sc_type}->{$radius} = 0 ;
         foreach my $restype1 (sort keys %{$aa3to1}) {
            $counts->{$sc_type}->{$radius}->{n_P} = 0 ;
            foreach my $restype2 (sort keys %{$aa3to1}) {

               if ($counts->{$sc_type}->{$radius}->{dnij_P__n_P}->{$restype1}->{$restype2} == 1) {
                  $g_ij->{$sc_type}->{$radius}->{$restype1}->{$restype2} = 
 1 / $counts->{$sc_type}->{$radius}->{nij_P}->{$restype1}->{$restype2} ;
               } else {
                  $g_ij->{$sc_type}->{$radius}->{$restype1}->{$restype2} =
                     $counts->{$sc_type}->{$radius}->{dnij_P__n_P}->{$restype1}->{$restype2} / 
                     ($counts->{$sc_type}->{$radius}->{nij_P}->{$restype1}->{$restype2} *
                     $counts->{$sc_type}->{$radius}->{maxafp_ij}->{$restype1}->{$restype2}) ;
               } 

               $G_ij_denom->{$sc_type}->{$radius} +=
                  $g_ij->{$sc_type}->{$radius}->{$restype1}->{$restype2} ;
            }
         }
      }
   }

   my ($out_fh, $out_fn) ;

   foreach my $sc_type (@{$types->{sc}}) {
      foreach my $radius (@{$types->{rad}}) {
         my $tt = $sc_type.'.'.$radius ;
         if (exists $in->{outfn_prefix}) {
            $out_fn->{$tt} = $in->{outfn_prefix}.'.'.$sc_type.'.'.$radius.'.pot' ;
            open($out_fh->{$tt}, ">$out_fn->{$tt}") ;
         } else {
            my $outprefix = "modtiepotential";
            ($out_fh->{$tt}, $out_fn->{$tt}) =
               tempfile($outprefix.'.'.$sc_type.'.'.$radius.".XXXXX",
                     SUFFIX => ".modtie.pot") ;
         }
      }
   }

   foreach my $sc_type (@{$types->{sc}}) {
      foreach my $radius (@{$types->{rad}}) {
         foreach my $restype1 (sort keys %{$aa3to1}) {
            $counts->{$sc_type}->{$radius}->{n_P} = 0 ;
            foreach my $restype2 (sort keys %{$aa3to1}) {
               $G_ij->{$sc_type}->{$radius}->{$restype1}->{$restype2} =

                 (400 * $g_ij->{$sc_type}->{$radius}->{$restype1}->{$restype2})/
                  $G_ij_denom->{$sc_type}->{$radius} ;

               $w_ij->{$sc_type}->{$radius}->{$restype1}->{$restype2} =
                 -1 * POSIX::log ( $G_ij->{$sc_type}->{$radius}->{$restype1}->{$restype2} ) ;

               my @outvals = ('w_ij', $sc_type, $radius, $restype1, $restype2,
                  sprintf("%.3f",
                  $w_ij->{$sc_type}->{$radius}->{$restype1}->{$restype2})) ;
               print join("\t", @outvals)."\n" ;
               print {$out_fh->{$sc_type.'.'.$radius}} join("\t", @outvals)."\n" ;

               @outvals = ('maxafp_ij', $sc_type, $radius, $restype1, $restype2,
                  sprintf("%.3f",
                  $counts->{$sc_type}->{$radius}->{maxafp_ij}->{$restype1}->{$restype2})) ;
               print {$out_fh->{$sc_type.'.'.$radius}} join("\t", @outvals)."\n" ;
               print join("\t", @outvals)."\n" ;
            }
         }
      }
   }


#for each type:
# {4 | 6 | 8 } {inter | intra | all} | {mm | ms | ss | all }

# global counts:
#       g_ij_num
#       g_ij_denom

# local counts:
#       dn_ij
#       n_p
#       n__p_i
#       n_j

}


sub build_roc_jackknife {

   my $in = shift;
   my $results ;


   my $num_neg = 0 ;
   my $num_pos = 0 ;
   my @theseinds ;
   if (exists $in->{indices}) { @theseinds = @{$in->{indices}} ;}
   else { @theseinds = (0 .. $#{$in->{scores}->{p}}) ;}

   my (@allscores, @labels) ; 
   foreach my $ind ( @theseinds) {
      push @allscores, @{$in->{scores}->{p}->[$ind]} ;
      push @labels, (('P') x ($#{$in->{scores}->{p}->[$ind]} + 1)) ;
      $num_pos += $#{$in->{scores}->{p}->[$ind]} + 1 ;

      push @allscores, @{$in->{scores}->{n}->[$ind]} ;
      push @labels, (('N') x ($#{$in->{scores}->{n}->[$ind]} + 1)) ;
      $num_neg += $#{$in->{scores}->{n}->[$ind]} + 1 ;
   }

   my @sort_order = sort {$allscores[$a] <=> $allscores[$b]} (0 .. $#allscores);

   my $last_score ;
   my $fp = 0; my $tp = 0 ;
   my $last_fp = 0; my $last_tp = 0 ;
   my $auc = 0 ;

   my $cump =  0; my $cumn = 0 ;
   my $optimal_thresh ;

   foreach my $ind (@sort_order) {
      if (!defined $last_score ||
          $allscores[$ind] != $last_score) {

         if (!defined $last_score) {$last_score = 'undef';}

         my $pbz = 0 ;
         if ($cump + $cumn > 0) {
            $pbz = sprintf("%.3f",
                   ($cump / $num_pos) /
                   ( ($cump / $num_pos) + ($cumn / $num_neg))) ;
         }

         my $curfpr = $fp / $num_neg ;
         my $curtpr = $tp / $num_pos ;

         push @{$results->{rocpoints}}, {
            x => sprintf("%.3f", $curfpr),
            y => sprintf("%.3f", $curtpr),
            last_score => $last_score,
            pbz => $pbz
         } ;

         $auc += _build_roc_trap_area($fp, $last_fp, $tp, $last_tp) ;

         my $curdist =
            ($results->{rocpoints}->[$#{$results->{rocpoints}}]->{x}) ** 2 +
            (1 - $results->{rocpoints}->[$#{$results->{rocpoints}}]->{y}) ** 2 ;

         if (!exists $optimal_thresh->{dist} ||
             $curdist < $optimal_thresh->{dist} || 
             ($optimal_thresh->{dist} == $curdist && 
              $optimal_thresh->{tpr} < $curtpr )) {

            $optimal_thresh = {
               tpr => $curtpr,
               fpr => $curfpr,
               dist => $curdist,
               thresh => $last_score } ;
         }

         $last_score = $allscores[$ind] ;
         $last_fp = $fp ;
         $last_tp = $tp ;

      }

      if ($labels[$ind] eq 'P') {
         $cump++ ;
         $tp++;
      } else {
         $cumn++ ;
         $fp++;
      }
   }

# do variance on the calculated z-score - jacknife the negative set
# see how the ROC's vary?

   my $curfpr = $fp / $num_neg ;
   my $curtpr = $tp / $num_pos ;

      push @{$results->{rocpoints}}, {
         x => sprintf("%.3f", $curfpr),
         y => sprintf("%.3f", $curtpr),
         pbz => sprintf("%.3f", ($cump / $num_pos) /
                      ( ($cump / $num_pos) + ($cumn / $num_neg))),
      last_score => $last_score
   } ;
   my $curdist= ($results->{rocpoints}->[$#{$results->{rocpoints}}]->{x}) ** 2 +
            (1 - $results->{rocpoints}->[$#{$results->{rocpoints}}]->{y}) ** 2 ;

   if (!exists $optimal_thresh->{dist} ||
       $curdist < $optimal_thresh->{dist} || 
       ( $optimal_thresh->{dist} == $curdist &&
         $optimal_thresh->{tpr} < $curtpr )) {
       $optimal_thresh = {
               tpr => $curtpr,
               fpr => $curfpr,
               dist => $curdist,
               thresh => $last_score } ;
   }

   $auc += _build_roc_trap_area($fp, $last_fp, $tp, $last_tp) ;
   $auc /= ($num_pos * $num_neg) ;

   $results->{auc} = $auc ;
   $results->{optimal} = $optimal_thresh ;

   return $results ;

#The receiver operating characteristic (ROC) curve was analyzed, the area under the curve was calculated, and the best threshold was defined as the VAS value that minimized the distance to the ideal point of the ROC curve (i.e., sensitivity = specificity = 1). When no significant difference was observed between the distances of two thresholds, the one associated with the highest sensitivity was retained.
}


sub build_roc {

   my $in = shift;
   my $results ;

   my @theseinds = (0 .. $#{$in->{scores}->{p}}) ;

   my $num_neg = 0 ;
   my $num_pos = 0 ;
   my (@allscores, @labels) ; 
   foreach my $ind ( @theseinds) {
      push @allscores, @{$in->{scores}->{p}->[$ind]} ;
      push @labels, (('P') x ($#{$in->{scores}->{p}->[$ind]} + 1)) ;
      $num_pos += $#{$in->{scores}->{p}->[$ind]} + 1 ;

      push @allscores, @{$in->{scores}->{n}->[$ind]} ;
      push @labels, (('N') x ($#{$in->{scores}->{n}->[$ind]} + 1)) ;
      $num_neg += $#{$in->{scores}->{n}->[$ind]} + 1 ;
   }

   my @sort_order = sort {$allscores[$a] <=> $allscores[$b]} (0 .. $#allscores);

   my $last_score ;
   my $fp = 0; my $tp = 0 ;
   my $last_fp = 0; my $last_tp = 0 ;
   my $auc = 0 ;

   my $cump =  0; my $cumn = 0 ;
   my $optimal_thresh ;

   foreach my $ind (@sort_order) {
      if (!defined $last_score ||
          $allscores[$ind] != $last_score) {

         if (!defined $last_score) {$last_score = 'undef';}

         my $pbz = 0 ;
         if ($cump + $cumn > 0) {
            $pbz = sprintf("%.3f",
                   ($cump / $num_pos) /
                   ( ($cump / $num_pos) + ($cumn / $num_neg))) ;
         }

         my $curfpr = $fp / $num_neg ;
         my $curtpr = $tp / $num_pos ;

         push @{$results->{rocpoints}}, {
            x => sprintf("%.3f", $curfpr),
            y => sprintf("%.3f", $curtpr),
            last_score => $last_score,
            pbz => $pbz
         } ;

         $auc += _build_roc_trap_area($fp, $last_fp, $tp, $last_tp) ;

         my $curdist =
            ($results->{rocpoints}->[$#{$results->{rocpoints}}]->{x}) ** 2 +
            (1 - $results->{rocpoints}->[$#{$results->{rocpoints}}]->{y}) ** 2 ;

         if (!exists $optimal_thresh->{dist} ||
             $curdist < $optimal_thresh->{dist} || 
             ($optimal_thresh->{dist} == $curdist && 
              $optimal_thresh->{tpr} < $curtpr )) {

            $optimal_thresh = {
               tpr => $curtpr,
               fpr => $curfpr,
               dist => $curdist,
               thresh => $last_score } ;
         }

         $last_score = $allscores[$ind] ;
         $last_fp = $fp ;
         $last_tp = $tp ;

      }

      if ($labels[$ind] eq 'P') {
         $cump++ ;
         $tp++;
      } else {
         $cumn++ ;
         $fp++;
      }
   }

# do variance on the calculated z-score - jacknife the negative set
# see how the ROC's vary?

   my $curfpr = $fp / $num_neg ;
   my $curtpr = $tp / $num_pos ;

      push @{$results->{rocpoints}}, {
         x => sprintf("%.3f", $curfpr),
         y => sprintf("%.3f", $curtpr),
         pbz => sprintf("%.3f", ($cump / $num_pos) /
                      ( ($cump / $num_pos) + ($cumn / $num_neg))),
      last_score => $last_score
   } ;
   my $curdist= ($results->{rocpoints}->[$#{$results->{rocpoints}}]->{x}) ** 2 +
            (1 - $results->{rocpoints}->[$#{$results->{rocpoints}}]->{y}) ** 2 ;

   if (!exists $optimal_thresh->{dist} ||
       $curdist < $optimal_thresh->{dist} || 
       ( $optimal_thresh->{dist} == $curdist &&
         $optimal_thresh->{tpr} < $curtpr )) {
       $optimal_thresh = {
               tpr => $curtpr,
               fpr => $curfpr,
               dist => $curdist,
               thresh => $last_score } ;
   }

   $auc += _build_roc_trap_area($fp, $last_fp, $tp, $last_tp) ;
   $auc /= ($num_pos * $num_neg) ;

   $results->{auc} = $auc ;
   $results->{optimal} = $optimal_thresh ;

   return $results ;

#The receiver operating characteristic (ROC) curve was analyzed, the area under the curve was calculated, and the best threshold was defined as the value that minimized the distance to the ideal point of the ROC curve (i.e., sensitivity = specificity = 1). When no significant difference was observed between the distances of two thresholds, the one associated with the highest sensitivity was retained.
}


sub _build_roc_trap_area {

   my ($x1, $x2, $y1, $y2) = @_ ;

   my $area = abs($x1 - $x2) * ($y1 + $y2) / 2 ;
   return $area ;

}


sub maxval_arrref {
   my $in = shift;
   my $max = $in->[0] ;
   foreach my $j ( 1 .. $#{$in}) {
      if ($in->[$j] > $max) {$max = $in->[$j];} }
   return $max ;
}


1 ;
