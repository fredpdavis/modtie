=head1 NAME

modtie.pm - routines for structure-based protein interaction prediction

=head1 DESCRIPTION

The modtie.pm perl library contains subroutines to identify and score
candidate protein--protein interactions based on similarity to template
protein complexes stored in PIBASE.

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


package modtie ;
require Exporter ;

use strict ;
use warnings ;
use vars qw($VERSION @ISA @EXPORT @EXPORT_OK) ;
our $VERSION = "1.11" ;

use Carp qw/croak/ ;

use Sys::Hostname qw/hostname/ ;
use File::Copy qw/move/ ;
use File::Path qw/mkpath/ ;
use Cwd ;

use POSIX qw/ceil floor/ ;
use File::Temp qw/tempfile tempdir/;
use Sys::Hostname qw/hostname/ ;

use modtie::complexes ;
use modtie::SGE ;
use modtie::yeast ;
use modtie::potentials ;

use modtie::pibase qw/safe_move safe_copy sid_2_domdir replace_undefs subset_extract get_salign parse_modeller_ali todload_bdp_ids/;


@ISA = qw/Exporter/ ;
@EXPORT_OK = qw/$modtie_specs/;
push @EXPORT_OK, qw/getstdinput/ ;
push @EXPORT_OK, qw/databaseid_2_seqid/ ;

push @EXPORT_OK, qw/build_roc/ ;

push @EXPORT_OK, qw/_getstdinput_format_modbaseoutput format_modbaseoutput/ ;
push @EXPORT_OK, qw/_getstdinput_calc_p_b_prior _getstdinput_assess_yeast_results/ ;
push @EXPORT_OK, qw/_getstdinput_runmodtie_modbase runmodtie_modbase/;
push @EXPORT_OK, qw/_getstdinput_runmodtie_scorecomplex runmodtie_scorecomplex/;
push @EXPORT_OK, qw/runmodtie_scorecomplex_alascan/;
push @EXPORT_OK, qw/_getstdinput_runmodtie_targetstrxs_template runmodtie_targetstrxs_template/ ;
push @EXPORT_OK, qw/_getstdinput_model_2_domains model_2_domains/;
push @EXPORT_OK, qw/_getstdinput_seqid_2_domains seqid_2_domains/;
push @EXPORT_OK, qw/_getstdinput_model_2_pdbstart model_2_pdbstart/;
push @EXPORT_OK, qw/_getstdinput_seqid_2_domainarch seqid_2_domainarch/;
push @EXPORT_OK, qw/_getstdinput_buildpotential_count /;
push @EXPORT_OK, qw/_getstdinput_buildpotential_postcalc/;
push @EXPORT_OK, qw/_getstdinput_runmodtie_benchmark runmodtie_benchmark/;
push @EXPORT_OK, qw/_getstdinput_runmodtie_roc runmodtie_roc/;
push @EXPORT_OK, qw/scoring_main salign_targ_tmpl_domains cut_domains/ ;

push @EXPORT_OK, qw/generate_interface_list extract_required_pibase_datafiles/ ;

=head2 set_modtie_specs()

   Title:       set_modtie_specs()
   Function:    gives modtie specifications
   Args:        NOTHING
   Returns:     returns completed modtie specifications

=cut

sub set_modtie_specs {

   my $modtie_specs = {};

#----- LOCAL CONFIGURATION DETAILS ------------------------------------------
#
# MODTIE installation directory
   $modtie_specs->{root} =
      '/groups/eddy/home/davisf/work/modtie/modtie_work/trunk/' ;

# Directory to store retrieved and generated files
   $modtie_specs->{runroot} =
      '/groups/eddy/home/davisf/work/modtie/modtie_201003/test_runs' ;

# PDB mirror location: expects files in pdb_dir/substr(CODE,1,2)/pdbCODE.ent.gz
   $modtie_specs->{pdb_dir}=
      '/groups/eddy/home/davisf/work/databases/pdb/data/structures' ;

# MODELLER binary location
   $modtie_specs->{modeller_bin} =
      "/groups/eddy/home/davisf/bin/modeller9v4/bin/mod9v4" ;

# SGE compute cluster configuration
   $modtie_specs->{cluster}->{head_node}= 'login-eddy'; # machine to qsub from
   $modtie_specs->{cluster}->{cluster_mode} = 1 ;       # cluster_mode (1=on)
   $modtie_specs->{cluster}->{qstat_sleep} = 60 ;       # poll frequency (sec)
   $modtie_specs->{cluster}->{priority} = 0;            # SGE job priority
   $modtie_specs->{cluster}->{numjobs} = 100 ;           # max # SGE jobs 
   $modtie_specs->{cluster}->{nodespecs} = '';          # Node requirements
#      "#\$ -l cpu500=false,cpu600=false,cpu933=false,cpu1500=false
#\$ -l pansali=1G,alto1=1G,alto2=1G,alto3=1G,diva1=1G,diva3=1G,scratch=1G" ;
#
#
#----- LOCAL CONFIGURATION DETAILS ---DONE-----------------------------------


   $modtie_specs->{cluster}->{minlines_salign} = 500 ;  # minimum SALIGN job
   $modtie_specs->{cluster}->{minlines_domcut} = 500 ;  # minimum domain cut job


   $modtie_specs->{dataroot} = $modtie_specs->{root}.'/modtie_data/' ;
   $modtie_specs->{binary_root} = $modtie_specs->{root}.'/src/auxil/' ;
   $modtie_specs->{libpath} = $modtie_specs->{root}.'/src/perl_api' ;
   $modtie_specs->{local_modbase_ali_dir} = $modtie_specs->{runroot}.
                                            '/modbase_ali/' ;
   $modtie_specs->{local_modbase_models_dir} = $modtie_specs->{runroot}.
                                               '/modbase_models/' ;



# DATABASE SPECS --------------------------------------------------
# PIBASE mysql db specs (used in local access run mode)
   $modtie_specs->{wget_specs}->{connect_timeout} = 30 ; #sec to kill wget
   $modtie_specs->{wget_specs}->{waittime}= 0.25 ; #seconds to wait to next wget
   $modtie_specs->{wget_specs}->{retries}= 5 ; #try 5 times to call wget

   $modtie_specs->{pibase_specs} = {
      tod_dir => $modtie_specs->{dataroot}."/pibase_tod",
      metatod_dir => $modtie_specs->{dataroot}."/pibase_metatod",
      subsets_dir => $modtie_specs->{dataroot}."/pibase_data/subsets_files",
      db => 'pibase_dbname',
      host => 'pibase_host',
      user => 'pibase_user' ,
      pass => 'pibase_password',
      old_root => '/groups/eddy/home/davisf/work/pibase/pibase2010', } ;

# MODBASE mysql db specs (used in local access run mode)
   $modtie_specs->{modbase_specs} = {
      db => 'modbase_dbname',
      host => 'modbase_host',
      user => 'modbase_user' ,
      pass => 'modbase_password' } ;

# PIBASE details necessary for accessing tod and metatod files properly.
# Sets up basedirectory specified in the table files
#  - not reflecting actual physical location on user's disk
#    (substituted when needed)
   $modtie_specs->{pibase_specs}->{old_dataroot} =
      $modtie_specs->{pibase_specs}->{old_root}.'/data' ;
   $modtie_specs->{pibase_specs}->{old_subsets_dir} =
      $modtie_specs->{pibase_specs}->{old_dataroot}.'/subsets_files' ;
   $modtie_specs->{pibase_specs}->{old_metatod_dir} =
      $modtie_specs->{pibase_specs}->{old_dataroot}.'/metatod' ;
   $modtie_specs->{pibase_specs}->{old_tod_dir} =
      $modtie_specs->{pibase_specs}->{old_dataroot}.'/tod' ;


# DEFAULTS ---------------------------------------------------------
   $modtie_specs->{runmode}->{default} = 'alilist_fl' ;
   $modtie_specs->{scoplevel}->{default} = 'superfamily' ;

# Default MODBASE access mode
   $modtie_specs->{modbase_access}->{default} = 'remote' ;

# THRESHOLDS -------------------------------------------------------
# if an interaction fails both then not printed;
   $modtie_specs->{printthresh_zscore} = -1.0;
   $modtie_specs->{printthresh_score} = 0; 

   $modtie_specs->{salign_quality_score_thresh} = 60 ;
   $modtie_specs->{complex_aln_tmpl_contacts_thresh} = 0.5 ;
   $modtie_specs->{complex_z_thresh} = -2.0 ;

# length ratio of domains (shortlen / longlen <= targ_tmpl_relative_size)
   $modtie_specs->{targ_tmpl_relative_size} = 0.60 ;


# MODBASE URL to retrieve all models for a seq_id
   $modtie_specs->{modbase_pdb_url} =
      'http://modbase.compbio.ucsf.edu/modbase-cgi/retrieve/modbase/'.
      '?type=model&seqID=' ; 

# MODBASE URL to retrieve all alignments for a seq_id
   $modtie_specs->{modbase_ali_url} =
      'http://modbase.compbio.ucsf.edu/modbase-cgi/retrieve/modbase/'.
      '?type=alignment&seqID=' ;

   @{$modtie_specs->{file_format}->{model_list}} =
      qw/seq_id align_id model_id seq_ident_global evalue target_length target_beg modpipe_style run mfu pdb_code baseali_dir basemodels_dir/ ;


   $modtie_specs->{target_domains_dir} = $modtie_specs->{runroot}.
                                         "/target_domains" ;
   $modtie_specs->{salign_ali_dir} = $modtie_specs->{runroot}."/salign_aln" ;


# By default, don't merge complexes through bridging sequences:
   $modtie_specs->{covalent_merge_fl} = 0 ;

# Set statistical potential locations
   $modtie_specs->{pot_dir}= $modtie_specs->{dataroot}."/potentials/";
   $modtie_specs->{pots} = [] ;
   $modtie_specs->{potsspec2fn} = {} ;
# oldpot file naming my $tfn_base = "pot.$type.allbut100.050917.$sc.$r.pot";
#   foreach my $r (qw/4 6 8/)
#      foreach my $type (qw/intra inter/)
#         foreach my $sc (qw/mm ms ss all/)
   foreach my $r (qw/8/) {
      foreach my $type (qw/inter/) {
         foreach my $sc (qw/ss/) {
            my $tfn_base = "modtie_potential.$type.060119.$sc.$r.pot";
            my $tfn = $modtie_specs->{pot_dir}.'/'.$tfn_base ;
            push @{$modtie_specs->{pots}}, {
               r => $r,
               type => 'pi',
               details => $sc,
               fn => $tfn,
               fn_base => $tfn_base
            } ;
            $modtie_specs->{potsspec2fn}->{'pi'}->{$type}->{$sc}->{$r} = $tfn ;
         }
      }
   }


# INTERNALLY SET FILE LOCATIONS
   $modtie_specs->{scripts}= $modtie_specs->{root}."/src/scripts" ;
   $modtie_specs->{templates_fn}->{default} = $modtie_specs->{dataroot}.
      "/templates/template_interfaces_cont1000.pibase2010.txt" ;
   $modtie_specs->{interface_clusters_fn} = $modtie_specs->{dataroot}.
      "/templates/interface_clusters_all.pibase2010.txt" ;

# DATA FILES FOR YEAST ASSESSMENT
   $modtie_specs->{yeastgfp_fn} = $modtie_specs->{dataroot}.
                                  "/yeastgfp/allOrfData.txt.yeastgfp" ;
   $modtie_specs->{seqid_2_sgd_fn} = $modtie_specs->{dataroot}.
                                     "/sgd/seqid_2_sgd.txt" ;
   $modtie_specs->{bind_complex2subunits} = $modtie_specs->{dataroot}.
                                       "/bind/20050811.complex2subunits.txt" ;
   $modtie_specs->{bind_nrints} = $modtie_specs->{dataroot}.
                                  "/bind/20050811.nrints.txt" ;
   $modtie_specs->{bind_division} = $modtie_specs->{dataroot}.
                                    "/bind/division_rgid.txt" ;
   $modtie_specs->{sgd_protein_properties_fn} = $modtie_specs->{dataroot}.
                                                "/sgd/protein_properties.tab" ;
   $modtie_specs->{sgd_features_fn} = $modtie_specs->{dataroot}.
                                      "/sgd/SGD_features.tab" ;
   $modtie_specs->{sgd_go_fn} = $modtie_specs->{dataroot}.
                                "/sgd/gene_association.sgd" ;

   $modtie_specs->{mips_complexcat_fn} = $modtie_specs->{dataroot}.
                                         "/mips/complexcat_data_14112005" ;
   $modtie_specs->{mips_ppi_fn} = $modtie_specs->{dataroot}.
                                  "/mips/PPI_221205.tab" ;

   $modtie_specs->{cellzome_complexes_fn} = $modtie_specs->{dataroot}.
                                          "/cellzome/cellzome06_complexes.txt" ;
   $modtie_specs->{cellzome_binary_fn} = $modtie_specs->{dataroot}.
                                          "/cellzome/socio-affinities.dat" ;

   $modtie_specs->{gconst}->{aa3to1} = {
      'ALA' => 'A' ,
      'ARG' => 'R' ,
      'ASN' => 'N' ,
      'ASP' => 'D' ,
      'CYS' => 'C' ,
      'GLN' => 'Q' ,
      'GLU' => 'E' ,
      'GLY' => 'G' ,
      'HIS' => 'H' ,
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
      'VAL' => 'V'
   } ;

   $modtie_specs->{gconst}->{mcatoms}->{" C  "} = 1 ;
   $modtie_specs->{gconst}->{mcatoms}->{" CA "} = 1 ;
   $modtie_specs->{gconst}->{mcatoms}->{" N  "} = 1 ;
   $modtie_specs->{gconst}->{mcatoms}->{" O  "} = 1 ;
   $modtie_specs->{gconst}->{mcatoms}->{" OXT"} = 1 ;

   $modtie_specs->{gconst}->{numaaatoms} = {};
   $modtie_specs->{gconst}->{numaaatoms}->{s} = {
      'ALA' => 1,
      'CYS' => 2,
      'ASP' => 4,
      'GLU' => 5,
      'PHE' => 7,
      'GLY' => 0,
      'HIS' => 6,
      'ILE' => 4,
      'LYS' => 5,
      'LEU' => 4,
      'MET' => 4,
      'ASN' => 4,
      'PRO' => 5,
      'GLN' => 5,
      'ARG' => 7,
      'SER' => 2,
      'THR' => 3,
      'VAL' => 3,
      'TYR' => 8,
      'TRP' => 10
   } ;

   foreach my $tres (keys %{$modtie_specs->{gconst}->{numaaatoms}->{s}}) {
      $modtie_specs->{gconst}->{numaaatoms}->{m}->{$tres} = 4 ;
      $modtie_specs->{gconst}->{numaaatoms}->{all}->{$tres} =
         4 + $modtie_specs->{gconst}->{numaaatoms}->{s}->{$tres};
   }


   $modtie_specs->{runmode_reqs} = {
      model_2_domains => {
         model_id => 1,
         seq_id => 1,
         align_id => 1,
         run_id => 1,
         seq_ident_global => 1,
         mfu => 1,
      },

      seqid_2_domains => {
         model_id => 1,
         seq_ident_global => 1,
         mfu => 1,
         start_resno => 1,
         end_resno => 1,
         subset_id => 1,
         subset_class => 1,
         numres_target_seg => 1,
         numres_template_seg => 1,
         numres_aln_seg => 1,
         numres_match_seg => 1,
         numres_target_full => 1,
         numres_template_full => 1,
         numres_aln_full => 1,
         numres_match_full => 1,
         subset_size => 1
      },

      seqid_2_domainarch => {
         model_id => 1,
         start_resno => 1,
         end_resno => 1,
         tmpl_subset_id => 1,
         subset_class => 1,
         subset_length => 1,
         segment_no => 1,
         num_segments => 1,
         seq_length => 1,
         scop_quality => 1,
      },

      seq_seq_score => {
         seq_1 => 1,
         seq_2 => 1
      },

      strx_strx_score => {
         fn => 1,
         subset_def_1 => 1,
         subset_def_2 => 1,
      },

      domarch_2_candidates => {
      }
   } ;

   $modtie_specs->{binaries} = locate_binaries({specs => $modtie_specs}) ;

   return $modtie_specs ;
}


=head2 generate_interface_list()

   Title:       generate_interface_list
   Function:    Generates template interface list by qeurying PIBASE to 
   Input:       $_ = specs hashref
   Return:      specs - hashref

=cut

sub generate_interface_list {
   require DBI ;

   my $in = shift ; my $specs ;
   if (exists $in->{specs}) { $specs = $in->{specs} ; }
   else { $specs = set_modtie_specs() ; }

   my ($dbh,) = modtie::pibase::connect_pibase() ;

# Make list of all interface cluster assignments
   {
   my ($bdp_id, $sid1, $sid2, $fam_pair, $cluster_no, $member_no) =
      modtie::pibase::mysql_fetchcols($dbh, "SELECT bdp_id, subset_id_1, ".
      " subset_id_2, scopclass_pair, cluster_no, member_no ".
      " FROM scop_interface_clusters WHERE cluster_level = \"fam\"") ;
   open(OUTF, ">interface_clusters_all.$$.txt") ;
   foreach my $j ( 0 .. $#{$bdp_id}) {
      my ($class1, $class2) = split(/\_/, $fam_pair->[$j]) ;
      my @outvals = ($class1, $class2, $cluster_no->[$j],
                     $sid1->[$j], $sid2->[$j], $member_no->[$j]) ;
      print OUTF join("\t", @outvals)."\n" ;
   }
   close(OUTF) ;
   }

# Make list of template interfaces (minimum 1000 num_contats) 
   {
   my ($bdp_id, $sid1, $sid2) = modtie::pibase::mysql_fetchcols($dbh,
      "SELECT a.bdp_id, a.subset_id_1, a.subset_id_2 ".
      "FROM scop_interface_clusters as a, intersubset_contacts as b ".
      "WHERE a.cluster_level = \"fam\" AND a.bdp_id = b.bdp_id AND ".
      "a.subset_id_1 = b.subset_id_1 AND a.subset_id_2 = b.subset_id_2 AND ".
      " member_no = 1 AND num_contacts >= 1000") ;

   open(OUTF, ">template_interfaces_cont1000.$$.txt") ;
   foreach my $j ( 0 .. $#{$bdp_id}) {
      print OUTF join("\t", $bdp_id->[$j], $sid1->[$j], $sid2->[$j])."\n" ; }
   close(OUTF) ;
   }

}


=head2 complete_modtie_specs(specs)

   Title:       complete_modtie_specs
   Function:    Fills in blanks in specs with default values.
   Input:       $_ = specs hashref
   Return:      specs - hashref

=cut

sub complete_modtie_specs {

   my $specs = shift  ;

   my $clean_specs = set_modtie_specs() ;
   my @keys = keys %{$clean_specs} ;

   foreach my $key (@keys) {
      if (!exists $specs->{$key}) {
         $specs->{$key} = $clean_specs->{$key} ; } }

   return $specs ;
}


=head2 locate_binaries()

   Function:    Returns location of necessary program binaries
   Return:      $_->{program} = program location.
      perl, zcat, rigor, subset_extractor, altloc_check
   Args:        none

=cut

sub locate_binaries {

   use Sys::Hostname qw/hostname/ ;

   my $hostname = Sys::Hostname::hostname ;

   my $in = shift ;
   my $specs ;
   if (exists $in->{specs}) { $specs = $in->{specs} ; }
   else { $specs = set_modtie_specs() ; }

   my $bin_rootdir = $specs->{'binary_root'};

   my $binaries ;
   $binaries->{'modeller'} = $specs->{modeller_bin} ;

   $binaries->{'subset_extractor'} = "$bin_rootdir/".
                                     "subset_extractor/subset_extractor" ;
   $binaries->{'altloc_check'} = "$bin_rootdir/altloc_check/altloc_check";
   $binaries->{'altloc_filter'} = "$bin_rootdir/altloc_filter/altloc_filter.pl";
   $binaries->{'kdcontacts'} = "$bin_rootdir/kdcontacts/kdcontacts" ;

   $binaries->{'perl'} = "perl" ;
   $binaries->{'zcat'} = "zcat" ;

# Old binary lcocation logic for multiple architectures
   if (0) {
      my $mach = 'i386' ;
   
   # GET MACHINE TYPE from UNAME
      my $proc_type = `uname -p` ; chomp $proc_type;
      if (!defined $proc_type ) {$proc_type = 'o64' ;}
   
      if ( ($hostname eq 'alto') || ($hostname eq 'diva') ) {
         $binaries->{'perl'} = 'perl5.6.1' ;
         $binaries->{'zcat'} = 'gzcat' ;
         $mach = 'sun4u' ;
      }
   
      if ($proc_type eq 'x86_64'||
          $hostname =~ /^o64/ || $hostname =~ /marimba/ ||
          $hostname =~ /^opt/ || $hostname =~ /lyre/) {
         $mach = 'o64' ;
      }
   
      if ($hostname =~ /^intel/) { $mach = 'ia64' ; }
      $binaries->{'subset_extractor'} .= ".$mach" ;
      $binaries->{'altloc_check'} .= ".$mach" ;
      $binaries->{'kdcontacts'} .= ".$mach" ;
   }


# SET BINARY LOCATIONS DEPENDING ON MACHINE TYPE
   if (! -e $binaries->{'subset_extractor'}) {
      $binaries->{'subset_extractor'} = "ERROR" ; }

   if (! -e $binaries->{'altloc_check'}) {
      $binaries->{'altloc_check'} = "ERROR";}

   if (! -e $binaries->{'altloc_filter'}) {
      $binaries->{'altloc_filter'} = "ERROR" ;
   } else {
      $binaries->{'altloc_filter'} = $binaries->{perl}." ".
         $binaries->{'altloc_filter'} ;
   }

   if (! -e $binaries->{'kdcontacts'}) {
      $binaries->{'kdcontacts'} = "ERROR" ; }

   return $binaries ;

}


=head2 sid_modelid_2_alidir()

   Title:       sid_modelid_2_alidir()
   Function:    given template subset_id and target model_id, returns
                directory holding alignment

   Args:        $_->{sid} = subset_id
                $_->{model_id} = target model_id
   Returns:     returns directory holding template-target SALIGN alignment

=cut

sub sid_modelid_2_alidir {

   my $in = shift ; my $specs ;
   if (exists $in->{specs}) { $specs = $in->{specs} ; }
   else { $specs = set_modtie_specs() ; }

   my ($bdp) = ($in->{sid} =~ /BDP([0-9]+)/) ;
   my $siddirnum = POSIX::floor($bdp / 10) ;

   my $dir = $specs->{salign_ali_dir}."/$siddirnum/$bdp/".
             substr($in->{modelid},0,3).'/'.$in->{modelid} ;
   return $dir ;
}



=head2 modelid_2_domdir()

   Title:       modelid_2_domdir()
   Function:    given model_id, returns directory holding PDB file
   Args:        $_ = model id
   Returns:     returns directory holding model PDB file

=cut

sub modelid_2_domdir {

   my $in = shift ; my $specs ;
   if (exists $in->{specs}) { $specs = $in->{specs} ; }
   else { $specs = set_modtie_specs() ; }

   my $modelid = $in->{modelid};
   return $specs->{target_domains_dir}.'/'.substr($modelid,0,3).'/'.$modelid ;

}


=head2 runmodtie_benchmark()

   Title:       runmodtie_benchmark()
   Function:    Benchmark the statistical potentials using known complex strx

=cut

sub runmodtie_benchmark {

#in: set of domain pairs
# potential locations  - if not the default modtie_specs pots

   my $in = shift;
   # use the scorecomplex() routine
   # - make sure to mode to new minafp fractional contacts scheme

   my $specs ;
   if (exists $in->{specs}) { $specs = $in->{specs} ; }
   else { $specs = set_modtie_specs() ; }

   my ($fn, $fh) ;

   if (!exists $in->{scores_fn}) {
      ($fh->{scores}, $fn->{scores}) = tempfile("benchmarkmodtie.scores.XXXXX",
                                    SUFFIX => '.modtie.out') ;
      close($fh->{scores}) ;
   } else {
      $fn->{scores} = $in->{scores_fn} ;
   }


   my $pots ;
   if (!exists $in->{pots}) {
      $pots = $specs->{pots} ;
   } else {
      $pots = $in->{pots} ;
   }
   readin_potentials($pots) ;

   my ($t_rmax, @trs) ;
   {
         foreach my $pot (@{$pots}) { push @trs, $pot->{r} ; }
         ($t_rmax, undef) = sort {$b <=> $a} @trs ;
   }

   my $pb_pl = scoring_pibase_preload() ;

   if (-s $fn->{scores}) {
         die "ERROR: $fn->{scores} is not empty" ; }

   open ($fh->{scores}, ">$fn->{scores}") ;

   foreach my $complex (@{$in->{complexes}}) {
      my $bdp_id = $complex->{bdp_id};
      my $sids = $complex->{subset_ids} ;
      print STDERR "now on bdp $bdp_id, sids ".join(',', @{$sids})."\n" ;

#make a file with the domain cuts in it before sending to get_domain_contacts;

      my $subsres = get_subsres($pb_pl->{tn}->{subsets_res}->{$bdp_id}) ;
      my ($temppdb_fh, $temppdb_fn) = tempfile() ;
      close($temppdb_fh) ;

      foreach my $sid (@{$sids}) {
         my $orig_path = modtie::pibase::sid_2_domdir({
                           sid => $sid,
                           specs => $specs,
                           })."/$sid.pdb.gz" ;
         my $t_com = $specs->{binaries}->{zcat}." $orig_path >> $temppdb_fn" ;
         system($t_com) ;
      }

      my $contacts = get_domain_contacts_complex({
         R_max => $t_rmax,
         R_list => \@trs,
         pdb_fn => $temppdb_fn,
         subset_residues => $subsres,
      }) ;

      foreach my $pot (@{$pots}) {
         my $rawscore = score_potential_complex({
            label => $temppdb_fn,
            contacts => $contacts->{respairs},
            pot => $pot,
            resnames => $contacts->{resnames},
            dompair2respair => $contacts->{dompair2respair}
         }) ;

#         print STDERR "rawscore is ".$rawscore->{whole}."\n" ;

         my $ranscores = get_rand_score_complex({
            contacts => $contacts->{respairs},
            pot => $pot,
            resnames => $contacts->{resnames},
            dompair2respair => $contacts->{dompair2respair}
         }) ;

         my $zscore = {};
         $zscore->{whole}= calc_zscore($rawscore->{whole}, $ranscores->{whole});
         if (exists $zscore->{whole}->{error}) {
            my @outvals = ('COMPLEX', join(', ', @{$sids}),
                        $pot->{fn}, $pot->{type}, $pot->{details}, $pot->{r},
                        $rawscore->{whole}, 'ERROR', $zscore->{whole}->{error});

            print {$fh->{scores}} join("\t", @outvals)."\n" ;
            print STDERR "ERROR: ".$zscore->{whole}->{error}."\n";next;
         }

         my @outvals = ('COMPLEX', join(', ', @{$sids}),
                     $pot->{fn}, $pot->{type}, $pot->{details}, $pot->{r},
                     $rawscore->{whole}, $zscore->{whole}->{z_score},
                     $zscore->{whole}->{z_prime}, $zscore->{whole}->{z_2},
                     $zscore->{whole}->{avg}, $zscore->{whole}->{min},
                     $zscore->{whole}->{min_tn}, $zscore->{whole}->{stdev},
                     $zscore->{whole}->{z_min}, $zscore->{whole}->{z_min_tn},
                     $zscore->{whole}->{false_pos},
                     join(', ', @{$zscore->{whole}->{z_bg}})) ;

         print {$fh->{scores}} join("\t", @outvals)."\n" ;

         foreach my $dom12sig ( keys %{$rawscore->{interface}}) {
            $zscore->{interface}->{$dom12sig} =
                calc_zscore($rawscore->{interface}->{$dom12sig},
                            $ranscores->{interface}->{$dom12sig});

            my $outdom12 = $dom12sig; $outdom12 =~ s/\n/\t/g ;
            if (exists $zscore->{whole}->{error}) {
               my @outvals = ('INTERFACE', $outdom12,
                           $pot->{fn}, $pot->{type}, $pot->{details}, $pot->{r},
                           $rawscore->{interface}->{$dom12sig}, 'ERROR',
                           $zscore->{interface}->{$dom12sig}->{error}) ;

               print {$fh->{scores}} join("\t", @outvals)."\n" ;
               print STDERR "ERROR: ".
                  $zscore->{interface}->{$dom12sig}->{error}."\n";
               next;
            }

            my @outvals = ('INTERFACE', $outdom12,
                        $pot->{fn}, $pot->{type}, $pot->{details}, $pot->{r},
                        $rawscore->{interface}->{$dom12sig},
                        $zscore->{interface}->{$dom12sig}->{z_score},
                        $zscore->{interface}->{$dom12sig}->{z_prime},
                        $zscore->{interface}->{$dom12sig}->{z_2},
                        $zscore->{interface}->{$dom12sig}->{avg},
                        $zscore->{interface}->{$dom12sig}->{min},
                        $zscore->{interface}->{$dom12sig}->{min_tn},
                        $zscore->{interface}->{$dom12sig}->{stdev},
                        $zscore->{interface}->{$dom12sig}->{z_min},
                        $zscore->{interface}->{$dom12sig}->{z_min_tn},
                        $zscore->{interface}->{$dom12sig}->{false_pos},
                        join(', ', @{$zscore->{interface}->{$dom12sig}->{z_bg}})) ;

            print {$fh->{scores}} join("\t", @outvals)."\n" ;
         }
      }
   }
}


=head2 modelid_2_domdir()

   Title:       modelid_2_domdir()
   Function:    Given calculated z-scores for positive and negative examples
                 calculate jack-knifed ROC curve.
   Args:        $_ = model id
   Returns:     prints out two files:
                File 1: statistical potential ROC performance
                  1. potential type
                  2. AUC
                  3. optimal z-score threshold
                  4. true positive rate at optimal z-score threshold
                  5. false positive rate at optimal z-score threshold
                  6. File containing ROC points
                  7. jack-knifed AUC mean
                  8. jack-knifed AUC stdev
                  9. jack-knifed optimal_thresh mean
                  10. jack-knifed optimal_thresh stdev
                  11. jack-knifed optimal_tpr mean
                  12. jack-knifed optimal_tpr stdev
                  13. jack-knifed optimal_fpr mean
                  14. jack-knifed optimal_fpr stdev

                File 2: ROC points for each potential:
                  1. x
                  2. y
                  3. z-score

=cut

sub runmodtie_roc {

   my $in = shift ;
   my $potinfo = $in->{potinfo} ;
   my $z_scores = $in->{z_scores} ;
   my @headers = qw/pottype auc optimal_thresh optimal_tpr optimal_fpr roc_fn/ ;
   print '#'.join("\t", @headers)."\n" ;

   my $rocpoints ;
   foreach my $pottype (keys %{$z_scores}) {
      my ($roc_fh, $roc_fn) ;
      if (!exists $in->{rocfile_prefix}) {
         ($roc_fh, $roc_fn) = tempfile("modtie_benchmark.$pottype.XXXXX",
                                       SUFFIX => '.rocpoints') ;
      } else {
         $roc_fn = $in->{rocfile_prefix}.".".$pottype.".rocpoints" ;
         if (-s $roc_fn) {
            die "FATAL ERROR: $roc_fn already exists"; }
         open ($roc_fh,">$roc_fn") ;
      }
      my $roccurve = build_roc({
         scores => $z_scores->{$pottype},
      });

      my $jacks ;
      foreach my $j ( 1 .. 20) {
         my @tind = (0 .. $#{$z_scores->{$pottype}->{p}}) ;
         fy_shuffle(\@tind) ;
         $#tind =($#tind - 20) ;
         my $roc_jack = build_roc_jackknife({
            indices => \@tind,
            scores => $z_scores->{$pottype},
         }) ;
         push @{$jacks->{auc}}, $roc_jack->{auc} ;
         push @{$jacks->{optimal_thresh}}, $roc_jack->{optimal}->{thresh} ;
         push @{$jacks->{optimal_tpr}}, $roc_jack->{optimal}->{tpr} ;
         push @{$jacks->{optimal_fpr}}, $roc_jack->{optimal}->{fpr} ;
      }

      my @jack_stats ;
      foreach my $key ( qw/auc optimal_thresh optimal_tpr optimal_fpr/) {
         my $stats = array_stats($jacks->{$key}) ;
         push @jack_stats, ($stats->{mean}, $stats->{stdev}) ;
      }

      my @outvals = ($pottype,
                     sprintf("%.3f", $roccurve->{auc}),
                     $roccurve->{optimal}->{thresh},
                     sprintf("%.3f", $roccurve->{optimal}->{tpr}),
                     sprintf("%.3f", $roccurve->{optimal}->{fpr}),
                     $roc_fn, @jack_stats) ;

      print join("\t", @outvals)."\n" ;

      foreach my $j ( 0 .. $#{$roccurve->{rocpoints}}) {
         print {$roc_fh} $roccurve->{rocpoints}->[$j]->{x}.' '.
                         $roccurve->{rocpoints}->[$j]->{y}.' '.
                         $roccurve->{rocpoints}->[$j]->{last_score}."\n" ;
      }
   }

}


=head2 array_stats()

   Title:       array_stats()
   Function:    given array of numbers, returns mean,stdev

   Args:        $_->[] = number
   Returns:     $->{mean} = mean
                $->{stdev} = standard deviation

=cut

sub array_stats {

   my $in = shift ;
   my $num = $#{$in} + 1 ;

   my $total = 0 ;
   foreach my $j ( 0 .. $#{$in}) {
      $total += $in->[$j] ; }
   my $mean = $total / $num ;

   my $stdevsq = 0 ;
   foreach my $j ( 0 .. $#{$in}) {
      $stdevsq += ($in->[$j] - $mean) * ($in->[$j] - $mean) ; }
   $stdevsq = $stdevsq / $num;
   my $stdev = sqrt($stdevsq) ;

   return {
      mean => $mean,
      stdev => $stdev
   } ;

}


=head2 runmodtie_scorecomplex_alascan()

   Title:       runmodtie_scorecomplex_alascan()
   Function:    Given list of complexes, ALA scan using MODTIE potentials

   Args:        $_->{pots} = statistical potential(s) to use for scoringj
                $_->{complexes}->{pdb_filename}->{domain_id} = {
                  chain => CHAIN_IDENTIFIER
                  start => START_RESIDUE
                  end => END_RESIDUE
                }

   Returns:
   Output:

=cut

sub runmodtie_scorecomplex_alascan {

#in:
# 1. target complex pdb file
# 2-n. domain definitions (get_subsres() format: hash->{ressig} = subset_id
#    use the same input format as in the candidates file list...

#out:
# - z-score, raw score, and bg scores

   my $in = shift;
   my $specs ;
   if (exists $in->{specs}) { $specs = $in->{specs} ; }
   else { $specs = set_modtie_specs() ; }

   my $pots ;
   if (!exists $in->{pots}) {
      $pots = $specs->{pots} ;
   } else {
      $pots = $in->{pots} ;
   }
   readin_potentials($pots) ;

   my ($fh, $fn) ;
   if (!exists $in->{out_fn}->{score_fn} &&
       !exists $in->{out_scores_fn}) {
      ($fh->{scores}, $fn->{scores}) =
         tempfile("complexscores.XXXXX", SUFFIX => ".out.modtie");
   } else {
      if (exists $in->{out_fn}->{score_fn}) {
         open($fh->{scores}, ">".$in->{out_fn}->{score_fn}) ;
      } elsif (exists $in->{out_scores_fn}) {
         open($fh->{scores}, ">".$in->{out_scores_fn}) ;
      }
   }


   foreach my $complex_pdbfn ( keys %{$in->{complexes}} ) {

      my $subsres = {} ;
      foreach my $domain_id ( keys %{$in->{complexes}->{$complex_pdbfn}}) {
         my ($tempcut_fh, $tempcut_fn) = tempfile(SUFFIX => ".pdb") ;
         close($tempcut_fh) ;

# 0. foreach domain: create a pdb file

         my $errfl = modtie::pibase::subset_extract({
            in_fn => $complex_pdbfn,
            out_fn => $tempcut_fn,
            chain => $in->{complexes}->{$complex_pdbfn}->{$domain_id}->{chain},
            start => $in->{complexes}->{$complex_pdbfn}->{$domain_id}->{start},
            end => $in->{complexes}->{$complex_pdbfn}->{$domain_id}->{end},
         }) ;

# 1. foreach domain: read in the created pdb file and populate subsres

         my $reslist = get_pdbfile_residuelist({pdb_fn => $tempcut_fn}) ;
         unlink $tempcut_fn;
         foreach my $curres ( @{$reslist->{residue_list}}) {
            $subsres->{$curres} = $domain_id ; }
      }

# 2. foreach complex: calculate interface contacts

      my ($t_rmax, @trs) ;
      {
         foreach my $pot (@{$pots}) { push @trs, $pot->{r} ; }
         ($t_rmax, undef) = sort {$b <=> $a} @trs ;
      }

      my $contacts = get_domain_contacts_complex({
         R_max => $t_rmax,
         R_list => \@trs,
         subset_residues => $subsres,
         pdb_fn => $complex_pdbfn
      }) ;

# this scores all interfaces together:
# if you want can restrict to pairwise contacts - to test wether  or not all pairwise contacts in actual pdb complexes score higher than the threshold
# as such will randomize and score as if one big pairwise interface

# 3. foreach complex: score the interface contacts

      foreach my $pot (@{$pots}) {
         my $rawscore = score_potential_complex_alascan({
            label => $complex_pdbfn,
            contacts => $contacts->{respairs},
            pot => $pot,
            resnames => $contacts->{resnames},
            dompair2respair => $contacts->{dompair2respair}
         }) ;

         my $ranscores = get_rand_score_complex_alascan({
            contacts => $contacts->{respairs},
            pot => $pot,
            resnames => $contacts->{resnames},
            dompair2respair => $contacts->{dompair2respair}
         }) ;

         my $zscore = {};
         $zscore->{whole}= calc_zscore($rawscore->{whole}, $ranscores->{whole});
         if (exists $zscore->{whole}->{error}) {
            my @outvals = ('COMPLEX', $complex_pdbfn,
                        join(", ", keys %{$in->{complexes}->{$complex_pdbfn}}),
                        $pot->{fn}, $pot->{type}, $pot->{details}, $pot->{r},
                        $rawscore->{whole}, 'ERROR', $zscore->{whole}->{error});

            print {$fh->{scores}} join("\t", @outvals)."\n" ;
            print STDERR "ERROR: ".$zscore->{whole}->{error}."\n";next;
         }

         my @outvals = ('COMPLEX', $complex_pdbfn,
                     join(", ", keys %{$in->{complexes}->{$complex_pdbfn}}),
                     $pot->{fn}, $pot->{type}, $pot->{details}, $pot->{r},
                     $rawscore->{whole}, $zscore->{whole}->{z_score},
                     $zscore->{whole}->{z_prime}, $zscore->{whole}->{z_2},
                     $zscore->{whole}->{avg}, $zscore->{whole}->{min},
                     $zscore->{whole}->{min_tn}, $zscore->{whole}->{stdev},
                     $zscore->{whole}->{z_min}, $zscore->{whole}->{z_min_tn},
                     $zscore->{whole}->{false_pos}) ;

         if (exists $in->{benchmark_fl}) {
            push @outvals, join(', ', @{$zscore->{whole}->{z_bg}}); }

         print {$fh->{scores}} join("\t", @outvals)."\n" ;

         foreach my $dom12sig ( keys %{$rawscore->{interface}}) {
            $zscore->{interface}->{$dom12sig} =
                calc_zscore($rawscore->{interface}->{$dom12sig},
                            $ranscores->{interface}->{$dom12sig});

            my $outdom12 = $dom12sig; $outdom12 =~ s/\n/\t/g ;
            my ($dom1, $dom2) = split(/\t/, $outdom12) ;
            if (exists $zscore->{whole}->{error}) {
               my @outvals = ('INTERFACE', $complex_pdbfn, $outdom12,
                           $pot->{fn}, $pot->{type}, $pot->{details}, $pot->{r},
                           $rawscore->{interface}->{$dom12sig}, 'ERROR',
                           $zscore->{interface}->{$dom12sig}->{error}) ;

               print {$fh->{scores}} join("\t", @outvals)."\n" ;
               print STDERR "ERROR: ".
                  $zscore->{interface}->{$dom12sig}->{error}."\n";
               next;
            }
            my $numcontacts = keys %{$contacts->{respairs}->{$dom12sig}} ;

            my @outvals = ('INTERFACE', $complex_pdbfn,
                        $outdom12,
                        $pot->{fn}, $pot->{type}, $pot->{details}, $pot->{r},
                        $rawscore->{interface}->{$dom12sig},
                        $zscore->{interface}->{$dom12sig}->{z_score},
                        $numcontacts,
                        $zscore->{interface}->{$dom12sig}->{z_prime},
                        $zscore->{interface}->{$dom12sig}->{z_2},
                        $zscore->{interface}->{$dom12sig}->{avg},
                        $zscore->{interface}->{$dom12sig}->{min},
                        $zscore->{interface}->{$dom12sig}->{min_tn},
                        $zscore->{interface}->{$dom12sig}->{stdev},
                        $zscore->{interface}->{$dom12sig}->{z_min},
                        $zscore->{interface}->{$dom12sig}->{z_min_tn},
                        $zscore->{interface}->{$dom12sig}->{false_pos}) ;

            if (exists $in->{benchmark_fl}) {
               push @outvals,
                  join(', ', @{$zscore->{interface}->{$dom12sig}->{z_bg}}) ; }

            print {$fh->{scores}} join("\t", @outvals)."\n" ;

            foreach my $resno (keys %{$rawscore->{interface_perres}->{$dom12sig}}) {

               my $resno_n = $resno ; $resno_n =~ s/\t/\n/ ;
               my $curresna ;
               if (exists $contacts->{resnames}->{$dom1}->{$resno_n}) {
                  $curresna = $contacts->{resnames}->{$dom1}->{$resno_n} ;
               } else {
                  $curresna = $contacts->{resnames}->{$dom2}->{$resno_n} ;
               }

               my $newrawscore = $rawscore->{interface}->{$dom12sig} -
                  $rawscore->{interface_perres}->{$dom12sig}->{$resno} + 
                  $rawscore->{interface_perres_ala}->{$dom12sig}->{$resno} ;

               my @tresnapart = ();
               my @tresnopart = ();
               my @toldraw = ();
               my @tnewraw = ();
               my @tdiffraw = ();
               foreach my $res2 (keys %{$rawscore->{interface_perres_partner}->{$dom12sig}->{$resno}}) {
                  my $res2_n = $res2 ; $res2_n =~ s/\t/\n/ ;
                  if (exists $contacts->{resnames}->{$dom1}->{$res2_n}) {
                     push @tresnapart, $contacts->{resnames}->{$dom1}->{$res2_n} ;
                  } else {
                     push @tresnapart, $contacts->{resnames}->{$dom2}->{$res2_n} ;
                  }
                  my $t = $res2 ; $t =~ s/\t/:/g ; push @tresnopart, $t ;
                  push @toldraw, sprintf("%.3f", $rawscore->{interface_perres_partner}->{$dom12sig}->{$resno}->{$res2}) ;
                  push @tnewraw, sprintf("%.3f", $rawscore->{interface_perres_partner_ala}->{$dom12sig}->{$resno}->{$res2}) ;
                  push @tdiffraw, sprintf("%.3f", $rawscore->{interface_perres_partner_ala}->{$dom12sig}->{$resno}->{$res2} - $rawscore->{interface_perres_partner}->{$dom12sig}->{$resno}->{$res2}) ;
               }

               $newrawscore = sprintf("%.3f", $newrawscore) ;
               my $newzscore =
                  calc_zscore($newrawscore,
                              $ranscores->{interface}->{$dom12sig});
               my $tres = $resno ; $tres =~ s/\n/:/g ; $tres =~ s/\t/:/g ;
               my @outvals2 = ('INTERFACE_ALA', $complex_pdbfn,
                     $outdom12, $tres, $curresna,
                     $pot->{fn}, $pot->{type}, $pot->{details}, $pot->{r},
                     $newrawscore,
                     sprintf("%.3f", ($newrawscore - $rawscore->{interface}->{$dom12sig})),
                     $newzscore->{z_score},
                     sprintf("%.3f",($newzscore->{z_score} -
                        $zscore->{interface}->{$dom12sig}->{z_score})),
                     join(',', @tresnapart),
                     join(',', @tresnopart),
                     join(',', @toldraw),
                     join(',', @tnewraw),
                     join(',', @tdiffraw),
                        ) ;
               print {$fh->{scores}} join("\t", @outvals2)."\n" ;
            }

         }
      }
   }

}


=head2 runmodtie_scorecomplex()

   Title:       runmodtie_scorecomplex()
   Function:    Given list of complexes, score using MODTIE potentials

   Args:        $_->{pots} = statistical potential(s) to use for scoring
                $_->{complexes}->{pdb_filename}->{domain_id} = {
                  chain => CHAIN_IDENTIFIER
                  start => START_RESIDUE
                  end => END_RESIDUE
                }

=cut

sub runmodtie_scorecomplex {

#in:
# 1. target complex pdb file
# 2-n. domain definitions (get_subsres() format: hash->{ressig} = subset_id
#    use the same input format as in the candidates file list...

#out:
# - z-score, raw score, and bg scores

   my $in = shift;
   my $specs ;
   if (exists $in->{specs}) { $specs = $in->{specs} ; }
   else { $specs = set_modtie_specs() ; }

   my $pots ;
   if (!exists $in->{pots}) {
      $pots = $specs->{pots} ;
   } else {
      $pots = $in->{pots} ;
   }
   readin_potentials($pots) ;

   my ($fh, $fn) ;
   if (!exists $in->{out_fn}->{score_fn} &&
       !exists $in->{out_scores_fn}) {
      ($fh->{scores}, $fn->{scores}) =
         tempfile("complexscores.XXXXX", SUFFIX => ".out.modtie");
   } else {
      if (exists $in->{out_fn}->{score_fn}) {
         open($fh->{scores}, ">".$in->{out_fn}->{score_fn}) ;
      } elsif (exists $in->{out_scores_fn}) {
         open($fh->{scores}, ">".$in->{out_scores_fn}) ;
      }
   }


   my $results ;
   foreach my $complex_pdbfn ( keys %{$in->{complexes}} ) {
      my $subsres = {} ;
      foreach my $domain_id ( keys %{$in->{complexes}->{$complex_pdbfn}}) {
         my ($tempcut_fh, $tempcut_fn) = tempfile(SUFFIX => ".pdb") ;
         close($tempcut_fh) ;

# 0. foreach domain: create a pdb file

         my $errfl = modtie::pibase::subset_extract({
            in_fn => $complex_pdbfn,
            out_fn => $tempcut_fn,
            chain => $in->{complexes}->{$complex_pdbfn}->{$domain_id}->{chain},
            start => $in->{complexes}->{$complex_pdbfn}->{$domain_id}->{start},
            end => $in->{complexes}->{$complex_pdbfn}->{$domain_id}->{end},
         }) ;

# 1. foreach domain: read in the created pdb file and populate subsres

         my $reslist = get_pdbfile_residuelist({pdb_fn => $tempcut_fn}) ;
         unlink $tempcut_fn;
         foreach my $curres ( @{$reslist->{residue_list}}) {
            $subsres->{$curres} = $domain_id ; }
      }

# 2. foreach complex: calculate interface contacts

      my ($t_rmax, @trs) ;
      {
         foreach my $pot (@{$pots}) { push @trs, $pot->{r} ; }
         ($t_rmax, undef) = sort {$b <=> $a} @trs ;
      }

      my $contacts = get_domain_contacts_complex({
         R_max => $t_rmax,
         R_list => \@trs,
         subset_residues => $subsres,
         pdb_fn => $complex_pdbfn
      }) ;
      $results->{contacts} = $contacts ;

# this scores all interfaces together:
# if you want can restrict to pairwise contacts - to test wether  or not all pairwise contacts in actual pdb complexes score higher than the threshold
# as such will randomize and score as if one big pairwise interface

# 3. foreach complex: score the interface contacts

      my @headers_complex= qw/COMPLEX pdb_fn domains
                     pot_fn pot_type pot_details pot_dist_threshold
                     rawscore_whole zscore_whole
                     zprime_whole z_2_whole avg_whole min_whole
                     min_tn_whole stdev_whole z_min_whole z_min_tn_whole
                     false_pos_whole/ ;
      print {$fh->{scores}} '#'.join("\t", @headers_complex)."\n";

      my @headers_interface = qw/INTERFACE pdb_fn dom1 dom2
                        pot_fn pot_type pot_details pot_dist_threshold
                        rawscore zscore numcontacts
                        z_prime z_2 avg min min_tn stdev z_min z_min_tn 
                        false_pos/ ;
      print {$fh->{scores}} '#'.join("\t", @headers_interface)."\n" ;

      foreach my $pot (@{$pots}) {
         my $rawscore = score_potential_complex({
            label => $complex_pdbfn,
            contacts => $contacts->{respairs},
            pot => $pot,
            resnames => $contacts->{resnames},
            dompair2respair => $contacts->{dompair2respair}
         }) ;
         $results->{rawscore}->{$pot} = $rawscore ;

         my $ranscores = get_rand_score_complex({
            contacts => $contacts->{respairs},
            pot => $pot,
            resnames => $contacts->{resnames},
            dompair2respair => $contacts->{dompair2respair}
         }) ;

         my $zscore = {};
         $zscore->{whole}= calc_zscore($rawscore->{whole}, $ranscores->{whole}) ;
         if (exists $zscore->{whole}->{error}) {
            my @outvals = ('COMPLEX', $complex_pdbfn,
                        join(", ", keys %{$in->{complexes}->{$complex_pdbfn}}),
                        $pot->{fn}, $pot->{type}, $pot->{details}, $pot->{r},
                        $rawscore->{whole}, 'ERROR', $zscore->{whole}->{error});

            print {$fh->{scores}} join("\t", @outvals)."\n" ;
            print STDERR "ERROR: ".$zscore->{whole}->{error}."\n";next;
         }

         my @outvals = ('COMPLEX', $complex_pdbfn,
                     join(", ", keys %{$in->{complexes}->{$complex_pdbfn}}),
                     $pot->{fn}, $pot->{type}, $pot->{details}, $pot->{r},
                     $rawscore->{whole}, $zscore->{whole}->{z_score},
                     $zscore->{whole}->{z_prime}, $zscore->{whole}->{z_2},
                     $zscore->{whole}->{avg}, $zscore->{whole}->{min},
                     $zscore->{whole}->{min_tn}, $zscore->{whole}->{stdev},
                     $zscore->{whole}->{z_min}, $zscore->{whole}->{z_min_tn},
                     $zscore->{whole}->{false_pos}) ;

         if (exists $in->{benchmark_fl}) {
            push @outvals, join(', ', @{$zscore->{whole}->{z_bg}}); }

         print {$fh->{scores}} join("\t", @outvals)."\n" ;

         foreach my $dom12sig ( keys %{$rawscore->{interface}}) {
            $zscore->{interface}->{$dom12sig} =
                calc_zscore($rawscore->{interface}->{$dom12sig},
                            $ranscores->{interface}->{$dom12sig});

            my $outdom12 = $dom12sig; $outdom12 =~ s/\n/\t/g ;
            if (exists $zscore->{whole}->{error}) {
               my @outvals = ('INTERFACE', $complex_pdbfn, $outdom12,
                           $pot->{fn}, $pot->{type}, $pot->{details}, $pot->{r},
                           $rawscore->{interface}->{$dom12sig}, 'ERROR',
                           $zscore->{interface}->{$dom12sig}->{error}) ;

               print {$fh->{scores}} join("\t", @outvals)."\n" ;
               print STDERR "ERROR: ".
                  $zscore->{interface}->{$dom12sig}->{error}."\n";
               next;
            }
            my $numcontacts = keys %{$contacts->{respairs}->{$dom12sig}} ;

            my @outvals = ('INTERFACE', $complex_pdbfn,
                        $outdom12,
                        $pot->{fn}, $pot->{type}, $pot->{details}, $pot->{r},
                        $rawscore->{interface}->{$dom12sig},
                        $zscore->{interface}->{$dom12sig}->{z_score},
                        $numcontacts,
                        $zscore->{interface}->{$dom12sig}->{z_prime},
                        $zscore->{interface}->{$dom12sig}->{z_2},
                        $zscore->{interface}->{$dom12sig}->{avg},
                        $zscore->{interface}->{$dom12sig}->{min},
                        $zscore->{interface}->{$dom12sig}->{min_tn},
                        $zscore->{interface}->{$dom12sig}->{stdev},
                        $zscore->{interface}->{$dom12sig}->{z_min},
                        $zscore->{interface}->{$dom12sig}->{z_min_tn},
                        $zscore->{interface}->{$dom12sig}->{false_pos}) ;

            if (exists $in->{benchmark_fl}) {
               push @outvals,
                  join(', ', @{$zscore->{interface}->{$dom12sig}->{z_bg}}) ; }

            print {$fh->{scores}} join("\t", @outvals)."\n" ;

         }
      }
   }

   return $results ;

}


=head2 get_pdbfile_residuelist()

   Title:       get_pdbfile_residuelist()
   Function:    Parses PDB file for list of residues
   Args:        $_->{pdb_fn} = name of PDB file
                $_->{pdb_fh} = filehandle of opened PDB file(optional)
   Returns:     $->{residue_list} = [residue_number."\n".chain_id , ...]

=cut

sub get_pdbfile_residuelist {

   my $in = shift;
   my $fh ;
   if (exists $in->{pdb_fh}) { $fh = $in->{pdb_fh};}
   else { open($fh, $in->{pdb_fn}) ; }

   my $reslist ;
   while (my $line = <$fh>) {
      chomp $line;
      my $lastsig = '' ;
      if ($line =~ /^ATOM/) {
         my $curresno = substr($line, 22, 5) ; $curresno =~ s/ //g ;
         my $curchain = substr($line, 21, 1) ;

         my $ressig = $curresno."\n".$curchain ;
         if ($ressig ne $lastsig) { push @{$reslist}, $ressig ; }
      }
   }

   return {
      residue_list => $reslist
   } ;

}


=head2 runmodtie_target_template_alignment()

   Title:       runmodtie_target_template_alignment()
   Function:    Given target-template alignment, score a candidate complex
   Args:        $_->{pot} = description of potential files (optional)
                $_->{complexes}->{pdb_filename}->{domain_id} = {
                  chain => CHAIN_IDENTIFIER
                  start => START_RESIDUE
                  end => END_RESIDUE
                }
   Returns:     nothing
   Output:      candidate complex scores

=cut

sub runmodtie_target_template_alignment {

   my $in = shift;
   my $specs ;
   if (exists $in->{specs}) { $specs = $in->{specs} ; }
   else { $specs = set_modtie_specs() ; }

   my $pots ;
   if (!exists $in->{pots}) {
      $pots = $specs->{pots} ;
   } else {
      $pots = $in->{pots} ;
   }
   readin_potentials($pots) ;

   my ($fh, $fn) ;
   if (!exists $in->{out_fn}->{score_fn}) {
      ($fh->{scores}, $fn->{scores}) =
         tempfile("complexscores.XXXXX", SUFFIX => ".out.modtie");
   } else {
      open($fh->{scores}, ">".$in->{out_fn}->{score_fn}) ;
   }


# out:
# z-score, raw score, and bg scores

   my @cheaders = ('COMPLEX', 'complex_id', 'domains', 'potential_path', 'potential_type', 'potential_details', 'potential_cutoff', 'rawscore', 'z-score', 'z_prime', 'z_2', 'avg_raw', 'min_raw', 'min_tn', 'stdev', 'z_min', 'z_min_tn', 'false_pos') ;
   print {$fh->{scores}} '#'.join("\t", @cheaders)."\n" ;

   my @iheaders = ('INTERFACE', 'complex_id', 'domain1', 'domain2', 'potential_path', 'potential_type', 'potential_details', 'potential_cutoff', 'rawscore', 'z-score', 'z_prime', 'z_2', 'avg_raw', 'min_raw', 'min_tn', 'stdev', 'z_min', 'z_min_tn', 'false_pos') ;
   print {$fh->{scores}} '#'.join("\t", @iheaders)."\n" ;



   foreach my $complex_id ( keys %{$in->{complexes}}) {

      print STDERR "now on complex $complex_id\n";

# 0. populate template subsres files - cut up template pdb files and read in?

      my $template_pdbfn = $in->{complexes}->{$complex_id}->{template_pdbfn} ;
      my $tmpl_subsres ;
      my ($resmap, $resnames_t) ;

      foreach my $domain_id (keys %{$in->{complexes}->{$complex_id}->{domains}}) {
         my ($tmpl_cutfh, $tmpl_cutfn) = tempfile(SUFFIX => ".pdb") ;
         close($tmpl_cutfh) ;

         my $errfl = modtie::pibase::subset_extract({
            in_fn => $template_pdbfn,
            out_fn => $tmpl_cutfn,
            chain =>
              $in->{complexes}->{$complex_id}->{domdef}->{$domain_id}->{template}->{chain},
            start =>
              $in->{complexes}->{$complex_id}->{domdef}->{$domain_id}->{template}->{start},
            end =>
              $in->{complexes}->{$complex_id}->{domdef}->{$domain_id}->{template}->{end}});

         my $reslist = get_pdbfile_residuelist({pdb_fn => $tmpl_cutfn}) ;
         foreach my $curres ( @{$reslist->{residue_list}}) {
            $tmpl_subsres->{$curres} = $domain_id ; }

         my $targ_cutfn = $in->{complexes}->{$complex_id}->{domains}->{$domain_id}->{target}->{pdb_fn};
         my $targetwascut_fl = 0;
         if (exists
             $in->{complexes}->{$complex_id}->{domdef}->{$domain_id}->{target}){
            my $tempfh ;
            ($tempfh, $targ_cutfn) = tempfile(SUFFIX => ".pdb") ;
            my $errfl = modtie::pibase::subset_extract({
            in_fn =>
              $in->{complexes}->{$complex_id}->{domains}->{$domain_id}->{target}->{pdb_fn},
            out_fn => $targ_cutfn,
            chain =>
               $in->{complexes}->{$complex_id}->{domdef}->{$domain_id}->{target}->{chain},
            start =>
               $in->{complexes}->{$complex_id}->{domdef}->{$domain_id}->{target}->{start},
            end => 
               $in->{complexes}->{$complex_id}->{domdef}->{$domain_id}->{target}->{end},
            }) ;

            $targetwascut_fl = 1;
         }

# 1. populate target subsres files (read in target_pdb_1 .. target_pdb_n)
# why? - leave for now; but if necessary copy same as template subsres populator
# excpet for an additional check to see if a domain definition is necessary (exists domdef?) - otherwise just read int he whole target_pdbfn as specified in the input
#
# need cuts anyway - to do alignments

# 2. get target - template alignments
# keep in mind the residue numbering in the individual target domain files may be redundant - dont assume that target resno are able to uniquely identify target residues; fixed it so contacts are now organized by pairs of domain_ids; then reisude names are retrieved from domain specific residue listings

         my $modeller_bin = $specs->{binaries}->{modeller} ;

         my $out_salign= modtie::pibase::get_salign({
            modeller_bin => $modeller_bin,
            pdb_fn_1 => $tmpl_cutfn,
            pdb_fn_2 => $targ_cutfn,
         }) ;


         if (exists $out_salign->{error}) {
            print STDERR "ERROR (complex $complex_id, domain $domain_id): ".
                         "SALIGN error\n" ;
            unlink $out_salign->{ali_log} ;
            next;
         }

         my ($ali_fh, $ali_fn)  = tempfile("align.XXXXX", SUFFIX => ".salign") ;
         open (SALIGNLOG, $out_salign->{ali_log}) ;
         my $gotin = 0 ;
         while (my $logline = <SALIGNLOG>) {
            if ($logline =~ /^Current/) {$gotin++;}
            if ($gotin == 3) {print $ali_fh $logline ;}
         }
         close(SALIGNLOG) ;
         close($ali_fh) ;
         unlink $out_salign->{ali_log};  unlink $out_salign->{ali_top}; 

         my $t_resmap = salign_resequiv_parse({fn => $ali_fn}) ;
         $resmap->{$domain_id} = $t_resmap->{resequiv} ;
         $resnames_t->{$domain_id} = $t_resmap->{resnames}->[1] ;

         unlink $tmpl_cutfn;
         unlink $ali_fn ;
         if ($targetwascut_fl) { unlink $targ_cutfn ; }
      }

      # 3. calculate template interdomain contacts
      my ($t_rmax, @trs) ;
      {
         foreach my $pot (@{$pots}) { push @trs, $pot->{r} ; }
         ($t_rmax, undef) = sort {$b <=> $a} @trs ;
      }

      my $tmpl_contacts = get_domain_contacts_complex({
         R_max => $t_rmax,
         R_list => \@trs,
         subset_residues => $tmpl_subsres,
         pdb_fn => $template_pdbfn
      }) ;

      my $targ_contacts= contact_equiv_complex({
         contacts => $tmpl_contacts->{respairs},
         seq => $resnames_t,
         resmap => $resmap,
         tmpl_subsres => $tmpl_subsres,
         dompair2respair => $tmpl_contacts->{dompair2respair},
      }) ;

      # 4. binscore
      foreach my $pot (@{$pots}) {
         my $rawscore = score_potential_complex({
            contacts => $targ_contacts->{respairs},
            pot => $pot,
            resnames => $targ_contacts->{resnames},
            dompair2respair => $targ_contacts->{dompair2respair},
         }) ;

         my $ranscores = get_rand_score_complex({
            contacts => $targ_contacts->{respairs},
            pot => $pot,
            resnames => $targ_contacts->{resnames},
            dompair2respair => $targ_contacts->{dompair2respair},
         }) ;

         my $zscore ={};
         $zscore->{whole}= calc_zscore($rawscore->{whole}, $ranscores->{whole}) ;
         if (exists $zscore->{whole}->{error}) {
            my @outvals = ('COMPLEX', $complex_id,
               join(", ", keys %{$in->{complexes}->{$complex_id}->{domains}}),
               $pot->{fn}, $pot->{type}, $pot->{details}, $pot->{r},
               $rawscore->{whole}, 'ERROR', $zscore->{whole}->{error});

            print {$fh->{scores}} join("\t", @outvals)."\n" ;
            print STDERR "ERROR: ".$zscore->{whole}->{error}."\n";next;
         }

         my @outvals = ('COMPLEX', $complex_id,
               join(", ", keys %{$in->{complexes}->{$complex_id}->{domains}}),
               $pot->{fn}, $pot->{type}, $pot->{details}, $pot->{r},
               $rawscore->{whole}, $zscore->{whole}->{z_score},
               $zscore->{whole}->{z_prime}, $zscore->{whole}->{z_2},
               $zscore->{whole}->{avg}, $zscore->{whole}->{min},
               $zscore->{whole}->{min_tn}, $zscore->{whole}->{stdev},
               $zscore->{whole}->{z_min}, $zscore->{whole}->{z_min_tn},
               $zscore->{whole}->{false_pos}) ;
#              join(',', @{$ranscores->{whole}})) ;

         print {$fh->{scores}} join("\t", @outvals)."\n" ;

         foreach my $dom12sig ( keys %{$rawscore->{interface}}) {
            $zscore->{interface}->{$dom12sig} =
                calc_zscore($rawscore->{interface}->{$dom12sig},
                            $ranscores->{interface}->{$dom12sig});

            my $outdom12 = $dom12sig; $outdom12 =~ s/\n/\t/g ;
            if (exists $zscore->{whole}->{error}) {
               my @outvals = ('INTERFACE', $complex_id, $outdom12,
                           $pot->{fn}, $pot->{type}, $pot->{details}, $pot->{r},
                           $rawscore->{interface}->{$dom12sig}, 'ERROR',
                           $zscore->{interface}->{$dom12sig}->{error}) ;

               print {$fh->{scores}} join("\t", @outvals)."\n" ;
               print STDERR "ERROR: ".
                  $zscore->{interface}->{$dom12sig}->{error}."\n";
               next;
            }

            my @outvals = ('INTERFACE', $complex_id, $outdom12,
                        $pot->{fn}, $pot->{type}, $pot->{details}, $pot->{r},
                        $rawscore->{interface}->{$dom12sig},
                        $zscore->{interface}->{$dom12sig}->{z_score},
                        $zscore->{interface}->{$dom12sig}->{z_prime},
                        $zscore->{interface}->{$dom12sig}->{z_2},
                        $zscore->{interface}->{$dom12sig}->{avg},
                        $zscore->{interface}->{$dom12sig}->{min},
                        $zscore->{interface}->{$dom12sig}->{min_tn},
                        $zscore->{interface}->{$dom12sig}->{stdev},
                        $zscore->{interface}->{$dom12sig}->{z_min},
                        $zscore->{interface}->{$dom12sig}->{z_min_tn},
                        $zscore->{interface}->{$dom12sig}->{false_pos}) ;
            print {$fh->{scores}} join("\t", @outvals)."\n" ;

         }
      }
   }

}


=head2 runmodtie_targetstrxs_template()

   Title:       runmodtie_targetstrxs_temaplate()
   Function:    Given target and template structures, align and score
                candidate complex
   Args:        $_->{pot} = description of potential files (optional)
                $_->{complexes}->{pdb_filename}->{domain_id} = {
                  chain => CHAIN_IDENTIFIER
                  start => START_RESIDUE
                  end => END_RESIDUE
                }
   Returns:     nothing
   Output:      candidate complex scores

=cut

sub runmodtie_targetstrxs_template {

   my $in = shift;
   my $specs ;
   if (exists $in->{specs}) { $specs = $in->{specs} ; }
   else { $specs = set_modtie_specs() ; }

   my $pots ;
   if (!exists $in->{pots}) {
      $pots = $specs->{pots} ;
   } else {
      $pots = $in->{pots} ;
   }
   readin_potentials($pots) ;

   my ($fh, $fn) ;
   if (!exists $in->{out_fn}->{score_fn}) {
j     ($fh->{scores}, $fn->{scores}) =
         tempfile("complexscores.XXXXX", SUFFIX => ".out.modtie");
   } else {
      open($fh->{scores}, ">".$in->{out_fn}->{score_fn}) ;
   }


# out:
# z-score, raw score, and bg scores

   my @cheaders = ('COMPLEX', 'complex_id', 'domains', 'potential_path', 'potential_type', 'potential_details', 'potential_cutoff', 'rawscore', 'z-score', 'z_prime', 'z_2', 'avg_raw', 'min_raw', 'min_tn', 'stdev', 'z_min', 'z_min_tn', 'false_pos') ;
   print {$fh->{scores}} '#'.join("\t", @cheaders)."\n" ;

   my @iheaders = ('INTERFACE', 'complex_id', 'domain1', 'domain2', 'potential_path', 'potential_type', 'potential_details', 'potential_cutoff', 'rawscore', 'z-score', 'z_prime', 'z_2', 'avg_raw', 'min_raw', 'min_tn', 'stdev', 'z_min', 'z_min_tn', 'false_pos') ;
   print {$fh->{scores}} '#'.join("\t", @iheaders)."\n" ;



   foreach my $complex_id ( keys %{$in->{complexes}}) {

      print STDERR "now on complex $complex_id\n";

# 0. populate template subsres files - cut up template pdb files and read in?

      my $template_pdbfn = $in->{complexes}->{$complex_id}->{template_pdbfn} ;
      my $tmpl_subsres ;
      my ($resmap, $resnames_t) ;

      foreach my $domain_id (keys %{$in->{complexes}->{$complex_id}->{domains}}) {
         my ($tmpl_cutfh, $tmpl_cutfn) = tempfile(SUFFIX => ".pdb") ;
         close($tmpl_cutfh) ;

         my $errfl = modtie::pibase::subset_extract({
            in_fn => $template_pdbfn,
            out_fn => $tmpl_cutfn,
            chain =>
              $in->{complexes}->{$complex_id}->{domdef}->{$domain_id}->{template}->{chain},
            start =>
              $in->{complexes}->{$complex_id}->{domdef}->{$domain_id}->{template}->{start},
            end =>
              $in->{complexes}->{$complex_id}->{domdef}->{$domain_id}->{template}->{end},
         }) ;

         my $reslist = get_pdbfile_residuelist({pdb_fn => $tmpl_cutfn}) ;
         foreach my $curres ( @{$reslist->{residue_list}}) {
            $tmpl_subsres->{$curres} = $domain_id ; }

         my $targ_cutfn = $in->{complexes}->{$complex_id}->{domains}->{$domain_id}->{target}->{pdb_fn};
         my $targetwascut_fl = 0;
         if (exists $in->{complexes}->{$complex_id}->{domdef}->{$domain_id}->{target}) {
            my $tempfh ;
            ($tempfh, $targ_cutfn) = tempfile(SUFFIX => ".pdb") ;
            my $errfl = modtie::pibase::subset_extract({
   in_fn => $in->{complexes}->{$complex_id}->{domains}->{$domain_id}->{target}->{pdb_fn},
   out_fn => $targ_cutfn,
   chain => $in->{complexes}->{$complex_id}->{domdef}->{$domain_id}->{target}->{chain},
   start => $in->{complexes}->{$complex_id}->{domdef}->{$domain_id}->{target}->{start},
   end => $in->{complexes}->{$complex_id}->{domdef}->{$domain_id}->{target}->{end}
            }) ;

            $targetwascut_fl = 1;
         }

# 1. populate target subsres files (read in target_pdb_1 .. target_pdb_n)
# why? - leave for now; but if necessary copy same as template subsres populator
# excpet for an additional check to see if a domain definition is necessary (exists domdef?) - otherwise just read int he whole target_pdbfn as specified in the input
#
# need cuts anyway - to do alignments

# 2. get target - template alignments
# keep in mind the residue numbering in the individual target domain files may be redundant - dont assume that target resno are able to uniquely identify target residues; fixed it so contacts are now organized by pairs of domain_ids; then reisude names are retrieved from domain specific residue listings

         my $modeller_bin = $specs->{binaries}->{modeller} ;

         my $out_salign= modtie::pibase::get_salign({
            modeller_bin => $modeller_bin,
            pdb_fn_1 => $tmpl_cutfn,
            pdb_fn_2 => $targ_cutfn,
         }) ;


         if (exists $out_salign->{error}) {
            print STDERR "ERROR (complex $complex_id, domain $domain_id): SALIGN error\n" ;
            unlink $out_salign->{ali_log} ;
            next;
         }

         my ($ali_fh, $ali_fn)  = tempfile("align.XXXXX", SUFFIX => ".salign") ;
         open (SALIGNLOG, $out_salign->{ali_log}) ;
         my $gotin = 0 ;
         while (my $logline = <SALIGNLOG>) {
            if ($logline =~ /^Current/) {$gotin++;}
            if ($gotin == 3) {print $ali_fh $logline ;}
         }
         close(SALIGNLOG) ;
         close($ali_fh) ;
         unlink $out_salign->{ali_log};  unlink $out_salign->{ali_top}; 

         my $t_resmap = salign_resequiv_parse({fn => $ali_fn}) ;
         $resmap->{$domain_id} = $t_resmap->{resequiv} ;
         $resnames_t->{$domain_id} = $t_resmap->{resnames}->[1] ;

         unlink $tmpl_cutfn;
         unlink $ali_fn ;
         if ($targetwascut_fl) { unlink $targ_cutfn ; }
      }

# 3. calculate template interdomain contacts
      my ($t_rmax, @trs) ;
      {
         foreach my $pot (@{$pots}) { push @trs, $pot->{r} ; }
         ($t_rmax, undef) = sort {$b <=> $a} @trs ;
      }

      my $tmpl_contacts = get_domain_contacts_complex({
         R_max => $t_rmax,
         R_list => \@trs,
         subset_residues => $tmpl_subsres,
         pdb_fn => $template_pdbfn
      }) ;

      my $targ_contacts= contact_equiv_complex({
         contacts => $tmpl_contacts->{respairs},
         seq => $resnames_t,
         resmap => $resmap,
         tmpl_subsres => $tmpl_subsres,
         dompair2respair => $tmpl_contacts->{dompair2respair},
      }) ;

# 4. binscore
      foreach my $pot (@{$pots}) {
         my $rawscore = score_potential_complex({
            contacts => $targ_contacts->{respairs},
            pot => $pot,
            resnames => $targ_contacts->{resnames},
            dompair2respair => $targ_contacts->{dompair2respair},
         }) ;

         my $ranscores = get_rand_score_complex({
            contacts => $targ_contacts->{respairs},
            pot => $pot,
            resnames => $targ_contacts->{resnames},
            dompair2respair => $targ_contacts->{dompair2respair},
         }) ;

         my $zscore ={};
         $zscore->{whole}= calc_zscore($rawscore->{whole}, $ranscores->{whole});
         if (exists $zscore->{whole}->{error}) {
            my @outvals = ('COMPLEX', $complex_id,
               join(", ", keys %{$in->{complexes}->{$complex_id}->{domains}}),
               $pot->{fn}, $pot->{type}, $pot->{details}, $pot->{r},
               $rawscore->{whole}, 'ERROR', $zscore->{whole}->{error});

            print {$fh->{scores}} join("\t", @outvals)."\n" ;
            print STDERR "ERROR: ".$zscore->{whole}->{error}."\n";next;
         }

         my @outvals = ('COMPLEX', $complex_id,
                join(", ", keys %{$in->{complexes}->{$complex_id}->{domains}}),
                $pot->{fn}, $pot->{type}, $pot->{details}, $pot->{r},
                $rawscore->{whole}, $zscore->{whole}->{z_score},
                $zscore->{whole}->{z_prime}, $zscore->{whole}->{z_2},
                $zscore->{whole}->{avg}, $zscore->{whole}->{min},
                $zscore->{whole}->{min_tn}, $zscore->{whole}->{stdev},
                $zscore->{whole}->{z_min}, $zscore->{whole}->{z_min_tn},
                $zscore->{whole}->{false_pos}) ;
#                     join(',', @{$ranscores->{whole}})) ;

         print {$fh->{scores}} join("\t", @outvals)."\n" ;

         foreach my $dom12sig ( keys %{$rawscore->{interface}}) {
            $zscore->{interface}->{$dom12sig} =
                calc_zscore($rawscore->{interface}->{$dom12sig},
                            $ranscores->{interface}->{$dom12sig});

            my $outdom12 = $dom12sig; $outdom12 =~ s/\n/\t/g ;
            if (exists $zscore->{whole}->{error}) {
               my @outvals = ('INTERFACE', $complex_id, $outdom12,
                           $pot->{fn}, $pot->{type}, $pot->{details}, $pot->{r},
                           $rawscore->{interface}->{$dom12sig}, 'ERROR',
                           $zscore->{interface}->{$dom12sig}->{error}) ;

               print {$fh->{scores}} join("\t", @outvals)."\n" ;
               print STDERR "ERROR: ".
                  $zscore->{interface}->{$dom12sig}->{error}."\n";
               next;
            }

            my @outvals = ('INTERFACE', $complex_id, $outdom12,
                        $pot->{fn}, $pot->{type}, $pot->{details}, $pot->{r},
                        $rawscore->{interface}->{$dom12sig},
                        $zscore->{interface}->{$dom12sig}->{z_score},
                        $zscore->{interface}->{$dom12sig}->{z_prime},
                        $zscore->{interface}->{$dom12sig}->{z_2},
                        $zscore->{interface}->{$dom12sig}->{avg},
                        $zscore->{interface}->{$dom12sig}->{min},
                        $zscore->{interface}->{$dom12sig}->{min_tn},
                        $zscore->{interface}->{$dom12sig}->{stdev},
                        $zscore->{interface}->{$dom12sig}->{z_min},
                        $zscore->{interface}->{$dom12sig}->{z_min_tn},
                        $zscore->{interface}->{$dom12sig}->{false_pos}) ;
            print {$fh->{scores}} join("\t", @outvals)."\n" ;

         }
      }

   }

}


=head2 runmodtie_modbase()

   Title:       runmodtie_modbase()
   Function:    Find and score candidate interactions within 1 set or between 2
                sets of proteins, specified as a list of seq_id or MODPIPE run
                identifiers

   Args:        $_->{run} = MODPIPE run (intra run)
                $_->{seqid_set} = file containing list of seq id (intra run)

                $_->{run1} = MODPIPE run identifier (inter run)
                $_->{seqid_set} = file containing list of seq id (set 1)
                $_->{run2} = MODPIPE run identifier (inter run)
                $_->{seqid_set} = file containing list of seq id (set 2)
                $_->{}...

   Tables:      modbase.newmodels:      run,seq_id, align_id, model_id,
                                        maxseq_ident, evalue, target_length,
                                        target_beg, run, model_mfu, template_id

                modbase.modpipe_runs:   align_path, model_path, modpipe_style

                modbase.templates:      pdb_code, template_id, align_id

                modbase.aasequences:    seq_id, sequence

=cut

sub runmodtie_modbase {

   my $bins = locate_binaries();

   my $in = shift;
   my $specs ;
   if (exists $in->{specs}) { $specs = $in->{specs} ; }
   else { $specs = set_modtie_specs() ; }

   my ($temp_fn, $temp_fh, $tcom) ;

   my @optionkeys = keys %{$in} ;
   print STDERR "Started runmodtie_modbase() with options: ".
      join(", ", @optionkeys)."\n";

# Make sure run or seqid_set specified for either intra- or inter-set run
   if (exists $in->{run} || exists $in->{seqid_set}) {
      $in->{xset_fl} = 0 ;
   } elsif ((exists $in->{seqid_set1} || exists $in->{run1}) &&
            (exists $in->{seqid_set2} || exists $in->{run2})) {
      $in->{xset_fl} = 1 ;
   } else {
      croak "ERROR: runtype not recognized: either specify run or seqid_set to predict complexes within a set of proteins, or run1/seqid_set1 and run2/seqid_set2 to predict interactions between two sets of proteins" ;
   }

   if (!exists $in->{local_modbase_models_dir}) {
      $in->{local_modbase_models_dir} = $specs->{local_modbase_models_dir} ; }

   if (!exists $in->{local_modbase_ali_dir}) {
      $in->{local_modbase_ali_dir} = $specs->{local_modbase_ali_dir} ; }

# If not specified, set scoplevel to default 
   if (!exists $in->{scoplevel}) {
      $in->{scoplevel} = $specs->{scoplevel}->{default} ; }

# If not specified, set MODBASE access to default remote
   if (!exists $in->{modbase_access}) {
      $in->{modbase_access} = $specs->{modbase_access}->{default} ; }

# Seq-id must be specified if modbase access is remote
   if ($in->{modbase_access} eq 'remote' &&
       !exists $in->{seqid_set} && !exists $in->{seqid_set1}) {
      croak "ERROR: if modbase_access is remote, must provide seq_id input" ;
   }

   my $cookie_fn = '' ;
   if (exists $in->{modweb_cookie}) { $cookie_fn = $in->{modweb_cookie}; }

# Initialize sequence set(s)
   my @tin_seqset ;
   if (exists $in->{seqid_set})  { push @tin_seqset, [0, $in->{seqid_set} ] ; }
   if (exists $in->{seqid_set1}) { push @tin_seqset, [1, $in->{seqid_set1}] ; }
   if (exists $in->{seqid_set2}) { push @tin_seqset, [2, $in->{seqid_set2}] ; }

# Initialize MODBASE run identifiers
   my @in_runs;
   if (exists $in->{run})       { push @in_runs, [0, $in->{run} ] ;}
   if (exists $in->{run1})      { push @in_runs, [1, $in->{run1}] ;}
   if (exists $in->{run2})      { push @in_runs, [2, $in->{run2}] ;}

# If not specified, set to default interface cluster file
   if (!exists $in->{templates_fn}) {
      $in->{templates_fn} = $specs->{templates_fn}->{default} ; }

# If not specified, set to default cluster mode
   if (!exists $in->{cluster_fl}) {
      $in->{cluster_fl} = $specs->{cluster}->{cluster_mode} ; }

# If cluster flag set, set specs
   if ($in->{cluster_fl}) {
      foreach my $t_key (keys %{$specs->{cluster}}) {
         if ($t_key eq 'cluster_mode') {next;}
         if (!exists $in->{cluster}->{$t_key}) {
            $in->{cluster}->{$t_key} = $specs->{cluster}->{$t_key}; }
      }
   }

# Initialize model entries, sequence set information
   my ($models, $model_entries, $seqset) ;

# Step 1. Input Prep: (i) get necessary model PDB/ALI files,
#                    (ii) get model details
#                   (iii) set seqid_setinfo (esp for xset_fl runs)

   if (!exists $in->{model_list} || !exists $in->{seqid_setinfo}) {

# Step 1a. Input Prep with local DBI access to MODBASE
    if ($in->{modbase_access} eq 'local') { #Setup model info/datafiles via DBI
      my $dbh ; $dbh->{mb} = connect_modbase() ;

# If sequence set specified, retrieve MODBASE run identifiers for each sequence
      my ($seqid2run, $runs) ;
      if ($#tin_seqset >= 0) {
         foreach my $t_seqset (@tin_seqset) {
            my ($set_no, $set_fn) = @$t_seqset ;
            $seqset->{$set_no} = readin_1col_arr({fn =>$set_fn}) ;
            if (exists $in->{model_details}) {next;}

            foreach my $tseqid ( @{$seqset->{$set_no}}) {
               my ($truns) = modtie::pibase::mysql_fetchcols($dbh->{mb},
                  "SELECT run FROM newmodels where seq_id = \"$tseqid\"");
               foreach my $trun (@{$truns}) {
                  $seqid2run->{$tseqid}->{$trun}++ ;
                  $runs->{$trun}->{$tseqid}++ ;
               }
            }
         }
      }

# If MODBASE run identifier specified, retrieve all sequence ids in the run
      if ($#in_runs >= 0) {
         foreach my $t_run (@in_runs) {
            my ($set_no,$set_run) = @$t_run ;
            my ($tseqids) = modtie::pibase::mysql_fetchcols($dbh->{mb},
              "SELECT distinct seq_id FROM newmodels where run = \"$set_run\"");

            foreach my $tseqid (@{$tseqids}) {
               push @{$seqset->{$set_no}}, $tseqid ;
               if (exists $in->{model_details}) {next;}
               $runs->{$set_run}->{$tseqid}++ ;
               $seqid2run->{$tseqid}->{$set_run}++ ;
            }
         }
      }

# Set up $seqset hash to point from seqid to set membership
      if ($in->{xset_fl} == 1 ) {
         foreach my $set_no (1 ..2) {
            foreach my $tseqid (@{$seqset->{$set_no}}) { # set to 1,2,or 3
               $seqset->{union}->{$tseqid} += $set_no ; } }
      } else {
         foreach my $tseqid (@{$seqset->{0}}) {
            $seqset->{union}->{$tseqid} = 0; }
      }


#Get MODPIPE run alignment and model paths
      my $rundata_copy = {} ; #copies over data from each run whose data
                              # is not in the local_modbase_(ali|models)_dir
      my $run2align_path ;
      my $run2model_path ;
      my $run2modpipe_style;
      foreach my $trun (keys %{$runs}) {
         print STDERR "STATUS: ".
                      "getting modpipe run info FOR $trun (".__LINE__."\n";

         if (exists $specs->{ali_dir}->{$trun} &&
             exists $specs->{model_dir}->{$trun}) {
            $run2align_path->{$trun} = $specs->{ali_dir}->{$trun};
            $run2model_path->{$trun} = $specs->{model_dir}->{$trun};
         } else {
            my @t = modtie::pibase::mysql_fetchcols($dbh->{mb},
              "SELECT align_path, model_path, modpipe_style FROM modpipe_runs ".
              "WHERE run=\"$trun\"");
            if ($#{$t[0]} >= 0) {
               $run2align_path->{$trun} = $t[0]->[0] ;
               $run2model_path->{$trun} = $t[1]->[0] ;
               $run2modpipe_style->{$trun} = $t[2]->[0] ;
            }
         }

         if (!defined $run2align_path->{$trun} ||
             !defined $run2model_path->{$trun}) {
           print STDERR "ERROR: align_path or model_path undefined for $trun\n";
            next;
         }

         if ($run2align_path->{$trun} ne $in->{local_modbase_ali_dir}) {
            $rundata_copy->{$trun}->{ali}++ ; }

         if ($run2align_path->{$trun} ne $in->{local_modbase_models_dir}) {
            $rundata_copy->{$trun}->{mod}++ ; }
      }

# Load model details from MODBASE tables: newmodels, templates, aasequences
      if (!exists $in->{binscores}) {
         print STDERR "STATUS: loading model info from modbase (".__LINE__."\n";

# If MODPIPE run name specified, load by run
         if ($#in_runs >= 0) {
            foreach my $t_runinfo (@in_runs) {
               my ($t_set_no, $trun) = @$t_runinfo ;
               print STDERR "STATUS: loading model info from modbase ".
                            "for run $trun (".__LINE__."\n";
               my $tselect = "SELECT a.seq_id, a.align_id, a.model_id, a.maxseq_ident, a.evalue, a.target_length, a.target_beg, a.run, a.model_mfu, b.pdb_code FROM newmodels as a, templates as b WHERE a.run = \"$trun\" and a.template_id = b.template_id and a.align_id = b.align_id ORDER by pdb_code" ;
               my $tmodel_entries ;
               ( $tmodel_entries->{seq_id},
                 $tmodel_entries->{align_id},
                 $tmodel_entries->{model_id},
                 $tmodel_entries->{seq_ident_global},
                 $tmodel_entries->{evalue},
                 $tmodel_entries->{target_length},
                 $tmodel_entries->{target_beg},
                 $tmodel_entries->{run},
                 $tmodel_entries->{mfu},
                 $tmodel_entries->{pdb_code} ) = 
               modtie::pibase::mysql_fetchcols($dbh->{mb},$tselect) ;

               if ($#{$tmodel_entries->{seq_id}} < 0) {
                  print STDERR "WARNING: no models found for run $trun\n" ;
                  next;
               }

               foreach my $tfield ( keys %{$tmodel_entries}) {
                  push @{$model_entries->{$tfield}},
                       @{$tmodel_entries->{$tfield}}; }
            }

         } else {

# If sequence sets are specified, query by seq_id
            foreach my $tseqid (keys %{$seqset->{union}}) {

               my $tselect = "SELECT a.seq_id, a.align_id, a.model_id, a.maxseq_ident, a.evalue, a.target_length, a.target_beg, a.run, a.model_mfu, b.pdb_code FROM newmodels as a, templates as b WHERE a.seq_id = \"$tseqid\" and a.template_id = b.template_id and a.align_id = b.align_id ORDER by pdb_code" ;

               my $tmodel_entries ;
               ( $tmodel_entries->{seq_id},
                 $tmodel_entries->{align_id},
                 $tmodel_entries->{model_id},
                 $tmodel_entries->{seq_ident_global},
                 $tmodel_entries->{evalue},
                 $tmodel_entries->{target_length},
                 $tmodel_entries->{target_beg},
                 $tmodel_entries->{run},
                 $tmodel_entries->{mfu},
                 $tmodel_entries->{pdb_code} ) = 
               modtie::pibase::mysql_fetchcols($dbh->{mb},$tselect) ;


               if ($#{$tmodel_entries->{seq_id}} < 0) {
                  print STDERR "WARNING: no model found for ".
                               "$tseqid ($seqset->{union}->{$tseqid})\n" ;
                  next;
               }
               foreach my $tfield ( keys %{$tmodel_entries}) {
                  push @{$model_entries->{$tfield}},
                       @{$tmodel_entries->{$tfield}};}
            }
         }
      }

# supplement model_entries with modpipe_style entries
      foreach my $j ( 0 .. $#{$model_entries->{run}}) {
         my $trun = $model_entries->{run}->[$j] ;
         if ($trun eq 'new') {push @{$model_entries->{modpipe_style}}, 'new';}
                        else {push @{$model_entries->{modpipe_style}}, 'old';}
      }

# Copy over files to cluster-accessible directory if necessary.
      print STDERR "STATUS: copying over to cluster run directory (".
                   __LINE__."\n";
      my $numrun_copy = keys %{$rundata_copy} ;
      if ($numrun_copy > 0 && !exists $in->{copyover_done}) {
         my $copythreshold = 200 ;
         my @copythese ; my @copydest ; my @copydest_dir ;

         foreach my $j (0 .. $#{$model_entries->{seq_id}}) {
            my $tseq_id = $model_entries->{seq_id}->[$j] ;
            my $tmodel_id = $model_entries->{model_id}->[$j] ;
            my $talign_id = $model_entries->{align_id}->[$j] ;
            my $trun = $model_entries->{run}->[$j] ;

            if (!exists $rundata_copy->{$trun}) {next;}

            if (exists $rundata_copy->{$trun}->{ali}) {
               my $src_alifn ;
               if ($run2modpipe_style->{$trun} eq 'new') {
                  $src_alifn = $run2align_path->{$trun}."/".
                     substr($tseq_id,0,3)."/$tseq_id/alignments".
                     "/$talign_id.ali"; ;
               } else {
                  $src_alifn = $run2align_path->{$trun}."/".
                     substr($talign_id,0,2)."/$talign_id/$talign_id.ali";
               }

               my $new_alidir = $in->{local_modbase_ali_dir}."/".
                  substr($tseq_id, 0, 3)."/$tseq_id/alignments" ;
               my $new_alifn = $new_alidir."/$talign_id.ali" ;

               push @copythese, $src_alifn ;
               push @copydest, $new_alifn ;
               push @copydest_dir, $new_alidir ;
            }

            if (exists $rundata_copy->{$trun}->{mod}) {
               my $src_modfn ;
               if ($run2modpipe_style->{$trun} eq 'new') {
                  $src_modfn = $run2model_path->{$trun}."/".
                               substr($tseq_id,0,3)."/$tseq_id/models".
                               "/$tmodel_id.pdb.gz" ;
               } else {
                  $src_modfn = $run2model_path->{$trun}."/".
                               substr($tmodel_id, 0, 2).'/'.$tmodel_id.'/'.
                               $tmodel_id.'.pdb.gz' ;
               }
               my $new_moddir = $in->{local_modbase_models_dir}."/".
                  substr($tseq_id, 0, 3)."/$tseq_id/models" ;
               my $new_modfn = $new_moddir."/$tmodel_id.pdb.gz" ;

               push @copythese, $src_modfn ;
               push @copydest, $new_modfn ;
               push @copydest_dir, $new_moddir ;
            }

            if ($#copythese >= $copythreshold) {
               foreach my $j ( 0 .. $#copythese) {
                  print STDERR "PROGRESS: copying to $copydest[$j]\n" ;
                  if (! -d $copydest_dir[$j]) { mkpath($copydest_dir[$j]) ; }
                  if (-s $copydest[$j]) {next;}
                  modtie::pibase::safe_copy($copythese[$j], $copydest[$j]) ;}
               @copythese = () ;
               @copydest = () ;
               @copydest_dir = () ;
            }
         }

         foreach my $j ( 0 .. $#copythese) {
            print STDERR "PROGRESS: copying to $copydest[$j]\n" ;
               if (! -d $copydest_dir[$j]) { mkpath($copydest_dir[$j]) ; }
               if (-s $copydest[$j]) {next;}
               modtie::pibase::safe_copy($copythese[$j], $copydest[$j]) ;}
         @copythese = () ;
         @copydest = () ;
         @copydest_dir = () ;
      }

# Step 1b. Input Prep with remote web access to MODBASE
    } else { #Setup model info/datafiles via MODBASE web access

# Load sequence_id sets
      foreach my $t_seqset (@tin_seqset) {
         my ($set_no, $set_fn) = @$t_seqset ;
         $seqset->{$set_no} = readin_1col_arr({fn =>$set_fn}) ; }

# Set up $seqset hash to point from seqid to set membership
      if ($in->{xset_fl} == 1 ) {
         foreach my $set_no (1 ..2) {
            foreach my $tseqid (@{$seqset->{$set_no}}) { # set to 1,2,or 3
               $seqset->{union}->{$tseqid} += $set_no ; } }
      } else {
         foreach my $tseqid (@{$seqset->{0}}) {
            $seqset->{union}->{$tseqid} = 0; }
      }


# Retrieve necessary files from MODBASE and parse model information
      my $run_specs = {} ;
      if (exists $in->{run})    { $run_specs->{run} = $in->{run}; } 
      if (exists $in->{run1})   { $run_specs->{run1} = $in->{run1}; } 
      if (exists $in->{run2})   { $run_specs->{run2} = $in->{run2}; } 

      $model_entries = webget_modbase_files({ seq_id => $seqset->{union},
                                              specs  => $specs,
                                              cookie => $cookie_fn,
                                              run_specs => $run_specs}) ;
      my $nummodels_read = $#{$model_entries->{model_id}} ;
      print STDERR "READ IN $nummodels_read models\n" ;
    }

# Print seqid set information to disk
    ($temp_fh->{seqid_setinfo}, $temp_fn->{seqid_setinfo}) =
         tempfile("seqid_setinfo_XXXXX", SUFFIX => '.modtie') ;
    write_seqid_setinfo({
      seqset => $seqset,
      out_fh => $temp_fh->{seqid_setinfo}
    }) ;
    close $temp_fh->{seqid_setinfo} ;

# Print model list to disk
    ($temp_fh->{model_list}, $temp_fn->{model_list}) =
      tempfile("model_list_XXXXX", SUFFIX => '.modtie') ;
    write_model_list({
      out_fh => $temp_fh->{model_list},
      model_entries => $model_entries,
      ali_dir => $in->{local_modbase_ali_dir},
      models_dir => $in->{local_modbase_models_dir},
    }) ;
    close($temp_fh->{model_list}) ;

   } else {
      $models = read_model_list({fn => $in->{model_list}}) ;
      $model_entries = $models->{model_entries} ;
      $seqset = read_seqid_setinfo({fn => $in->{seqid_setinfo}}) ;

      $temp_fn->{model_list} = $in->{model_list} ;
      $temp_fn->{seqid_setinfo} = $in->{seqid_setinfo} ;
   }

   if (exists $in->{prepare_models_only}) {
      print STDERR "Models have been prepared\n"; exit(1); }

# Steps 2-4. Domain assignments

# Step 2. assign model domains (if assignment file not specified)
   if (!exists $in->{model_domains}) {
      ($temp_fh->{model_domains_out}, $temp_fn->{model_domains_out}) =
         tempfile("model_domains_out_XXXXX", SUFFIX => '.modtie') ;
      ($temp_fh->{model_domains_err}, $temp_fn->{model_domains_err}) =
         tempfile("model_domains_err_XXXXX", SUFFIX => '.modtie') ;
      close($temp_fh->{model_domains_out}) ;
      close($temp_fh->{model_domains_err}) ;

      print STDERR "Running assign_model_domains\n" ;
      if ($in->{cluster_fl} == 0) {
         $tcom = $specs->{scripts}."/assign_model_domains.mti.pl < $temp_fn->{model_list} 2> $temp_fn->{model_domains_err} > $temp_fn->{model_domains_out}" ;
         print STDERR "local run: $tcom\n" ;
         system($tcom) ;
         if (-z $temp_fn->{model_domains_err}) {
            unlink $temp_fn->{model_domains_err};}

      } else {

         my $split_dir = tempdir("splits_assign_model_domains.XXXXX") ;
         my $splits = modtie::SGE::_clust_split_ins({
            fn => $temp_fn->{model_list},
            dir => $split_dir,
            numjobs => $in->{cluster}->{numjobs}
         });

         my ($sgescript_fh, $sgescript_fn) =
            tempfile("mt.assign_model_domains.XXXXX", SUFFIX => ".SGE.sh") ;
         my $sge_outdir = tempdir("SGEOUT.assign_model_domains.XXXXX") ;

         print {$sgescript_fh} "#!/bin/csh
#\$ -S /bin/csh
#\$ -cwd
#\$ -o $sge_outdir
#\$ -e $sge_outdir
#\$ -r y
$in->{cluster}->{nodespecs}
#\$ -p $in->{cluster}->{priority}
#\$ -t 1-$splits->{numjobs}

set tasks1=( $splits->{tasklist} )

set input1=\$tasks1[\$SGE_TASK_ID\]

set curdir=`pwd`
set curhost=`hostname`
set curtime=`date`
set scratchdir=/scratch/fred/\$input1.\$\$\

rm -rf \$scratchdir
mkdir -p \$scratchdir

cp $specs->{scripts}/assign_model_domains.mti.pl \$scratchdir
cp $split_dir/\$input1 \$scratchdir

cd \$scratchdir


echo \"#sgejob run started on \$curhost at \$curtime\"
perl assign_model_domains.mti.pl < \$input1

set curtime=`date`
echo \"#sgejob run finished on \$curhost at \$curtime\"

rm -f \$scratchdir/\$input1 \$scratchdir/assign_model_domains.mti.pl
cd \$curdir
rmdir \$scratchdir\n" ;
         close($sgescript_fh) ;

         my $qsub_job_id = modtie::SGE::_clust_qsub({
            host => $in->{cluster}->{head_node},
            sgescript_fn => $sgescript_fn,
         }) ;

         while (1) {
            sleep $specs->{cluster}->{qstat_sleep} ;
            my $job_status= modtie::SGE::_clust_qstat({
               host => $in->{cluster}->{head_node},
               job_id => $qsub_job_id
            });
            if ($job_status) {last;}
         }

         modtie::SGE::_clust_merge_outs({
            script_fn => $sgescript_fn,
            out_fn => $temp_fn->{model_domains_out},
            err_fn => $temp_fn->{model_domains_err},
            job_id => $qsub_job_id,
            outdir => $sge_outdir,
            numjobs => $splits->{numjobs}
         }) ;
      }
   } else {
      $temp_fn->{model_domains_out} = $in->{model_domains} ;
   }

# Step 3. assign sequence domains (if assignment file not specified)
   if (!exists $in->{seqid_domains}) {
      ($temp_fh->{seqid_domains_out}, $temp_fn->{seqid_domains_out}) =
         tempfile("seqid_domains_out_XXXXX", SUFFIX => '.modtie') ;
      ($temp_fh->{seqid_domains_err}, $temp_fn->{seqid_domains_err}) =
         tempfile("seqid_domains_err_XXXXX", SUFFIX => '.modtie') ;
      close($temp_fh->{seqid_domains_err}) ;
      close($temp_fh->{seqid_domains_out}) ;
      print STDERR "set seqid_domains_out: $temp_fn->{seqid_domains_out} ".
                   "(err $temp_fn->{seqid_domains_err})\n" ;

      $tcom = $specs->{scripts}."/assign_seqid_domains.mti.pl < $temp_fn->{model_domains_out} 2>$temp_fn->{seqid_domains_err} >$temp_fn->{seqid_domains_out}";
      print STDERR "Running assign_seqid_domains\n" ;
      system($tcom) ;
      if (-z $temp_fn->{seqid_domains_err}) {
         unlink $temp_fn->{seqid_domains_err};}
   } else {
      $temp_fn->{seqid_domains_out} = $in->{seqid_domains} ;
   }

# Step 3b. grep first PDB ATOM record residue number for each model
# added (100925_2013) to adjust for new residue numbering in ModPipe > svn-r533.
   if (!exists $in->{model_pdbstart}) {
      ($temp_fh->{model_pdbstart_out}, $temp_fn->{model_pdbstart_out}) =
         tempfile("model_pdbstart_out_XXXXX", SUFFIX => '.modtie') ;
      ($temp_fh->{model_pdbstart_err}, $temp_fn->{model_pdbstart_err}) =
         tempfile("model_pdbstart_err_XXXXX", SUFFIX => '.modtie') ;
      close($temp_fh->{model_pdbstart_err}) ;
      close($temp_fh->{model_pdbstart_out}) ;

      $tcom = $specs->{scripts}."/calc_model_pdbstart.mti.pl < $temp_fn->{seqid_domains_out} 2>$temp_fn->{model_pdbstart_err} >$temp_fn->{model_pdbstart_out}";
      print STDERR "Running calc_model_pdbstart\n" ;
      system($tcom) ;
      if (-z $temp_fn->{model_pdbstart_err}) {
         unlink $temp_fn->{model_pdbstart_err};}
   } else {
      $temp_fn->{model_pdbstart_out} = $in->{model_pdbstart} ;
   }


# Step 4. build sequence domain architectures (if file not specified)
   if (!exists $in->{seqid_domainarch}) {
      ($temp_fh->{seqid_domainarch_out}, $temp_fn->{seqid_domainarch_out}) =
         tempfile("seqid_domainarch_out_XXXXX", SUFFIX => '.modtie') ;
      ($temp_fh->{seqid_domainarch_err}, $temp_fn->{seqid_domainarch_err}) =
         tempfile("seqid_domainarch_err_XXXXX", SUFFIX => '.modtie') ;
      close($temp_fh->{seqid_domainarch_err}) ;
      close($temp_fh->{seqid_domainarch_out}) ;

      $tcom = $specs->{scripts}."/assign_seqid_domainarch.mti.pl < $temp_fn->{seqid_domains_out} 2>$temp_fn->{seqid_domainarch_err} >$temp_fn->{seqid_domainarch_out}";
      print STDERR "Running assign_seqid_domarch\n" ;
      system($tcom) ;
      if (-z $temp_fn->{seqid_domainarch_err}) {
         unlink $temp_fn->{seqid_domainarch_err};}
   } else {
      $temp_fn->{seqid_domainarch_out} = $in->{seqid_domainarch} ;
   }

#100326_1127 THIS STEP IS UNNECESSARY - don't really need aa sequence info
#  changes: 1. edit out aaseq expectations through scoring_main(); 
#  changes: 2. just read in the domainarch output from last step.
# 1d. list seq_id, aaseq, domarch (if file not specified)
#   if (!exists $in->{seqid_aaseq_arch}) {
#      ($temp_fh->{seqid_aaseq_arch_out}, $temp_fn->{seqid_aaseq_arch_out}) =
#         tempfile("seqid_aaseq_arch_out_XXXXX", SUFFIX => '.modtie') ;
#      open(DOMARCH, $temp_fn->{seqid_domainarch_out}) ;
#      while (my $line = <DOMARCH>) {
#         chomp $line;
#         my ($t_seqid, $t_length, $t_cover_len, $t_num_domains, $t_domarch_ss,
#             $t_domarch_s, $t_domarch) = split(/\t/, $line) ;
#         print {$temp_fh->{seqid_aaseq_arch_out}}
#            join("\t", ($t_seqid, $aaseq->{$t_seqid}, $t_domarch))."\n";
#      }
#      close(DOMARCH) ;
#      close($temp_fh->{seqid_aaseq_arch_out});
#   } else {
#      $temp_fn->{seqid_aaseq_arch_out} = $in->{seqid_aaseq_arch} ;
#   }


# Step 5. determine candidate domain-domain interactions, if alilist or
#            ali_done not specified
#
   if (!exists $in->{alilist} && !exists $in->{ali_done}) {
      print STDERR "Running candidate generator (ali, cut lister)\n" ;
      ($temp_fh->{alilist_out}, $temp_fn->{alilist_out}) =
         tempfile("alilist_out_XXXXX", SUFFIX => '.modtie') ;
      ($temp_fh->{cutlist_out}, $temp_fn->{cutlist_out}) =
         tempfile("cutlist_out_XXXXX", SUFFIX => '.modtie') ;
      close($temp_fh->{alilist_out}) ;
      close($temp_fh->{cutlist_out}) ;

      scoring_main({
         xset_fl => $in->{xset_fl},
         scoplevel => $in->{scoplevel},
         alilist_fl => 1,
         seqset => $seqset,
         in_fn => {
            partial_alipaths_err_fn => $in->{partial_alipaths_err},
            seqid_setinfo_fn => $temp_fn->{seqid_setinfo},
            seqid_domainarch_fn => $temp_fn->{seqid_domainarch_out},
            seqdomains_fn => $temp_fn->{seqid_domains_out},
            model_pdbstart_fn => $temp_fn->{model_pdbstart_out},
            templates_fn => $in->{templates_fn}
         },
         out_fn => {
            alilist_fn => $temp_fn->{alilist_out},
            cutlist_fn => $temp_fn->{cutlist_out},
         }
      }) ;

   } else {
      if (exists $in->{alilist}) {
         $temp_fn->{alilist_out} = $in->{alilist} ;}
   }


# Step 6. Extract domain PDB files for each putatively interacting target domain
   if (exists $in->{cutlist}) {
      $temp_fn->{cutlist_out} = $in->{cutlist} ; }

   if (!exists $in->{cuts_done} && -s $temp_fn->{cutlist_out}) {
      print STDERR "Cutting up target domains\n" ;
      ($temp_fh->{cutpaths_out}, $temp_fn->{cutpaths_out}) =
         tempfile("cutpaths_out_XXXXX", SUFFIX => '.modtie') ;
      ($temp_fh->{cutpaths_err}, $temp_fn->{cutpaths_err}) =
         tempfile("cutpaths_err_XXXXX", SUFFIX => '.modtie') ;

      if ($in->{cluster_fl} == 0) {

         $tcom = $specs->{scripts}.
            "/cut_target_domains.mti.pl < $temp_fn->{cutlist_out} ".
            "2>$temp_fn->{cutpaths_err} >$temp_fn->{cutpaths_out}";
         print STDERR "RUN: $tcom\n" ;
         system($tcom) ;
         if (-z $temp_fn->{cutpaths_err}) {
            unlink $temp_fn->{cutpaths_err};}

      } else {
         my $split_dir = tempdir("splits_cut_target_domains.XXXXX") ;
         my $splits = modtie::SGE::_clust_split_ins({
            fn => $temp_fn->{cutlist_out},
            dir => $split_dir,
            minlines => $specs->{cluster}->{minlines_domcut},
            numjobs => $in->{cluster}->{numjobs}
         });

         my ($sgescript_fh, $sgescript_fn) =
            tempfile("mt.cut_target_domains.XXXXX", SUFFIX => ".SGE.sh") ;
         my $sge_outdir = tempdir("SGEOUT.cut_target_domains.XXXXX") ;

         print {$sgescript_fh} "#!/bin/csh
#\$ -S /bin/csh
#\$ -cwd
#\$ -o $sge_outdir
#\$ -e $sge_outdir
#\$ -r y
$in->{cluster}->{nodespecs}
#\$ -p $in->{cluster}->{priority}
#\$ -t 1-$splits->{numjobs}

set tasks1=( $splits->{tasklist} )

set input1=\$tasks1[\$SGE_TASK_ID\]

set curdir=`pwd`
set curhost=`hostname`
set curtime=`date`
set scratchdir=/scratch/fred/\$input1.\$\$\

rm -rf \$scratchdir
mkdir -p \$scratchdir

cp $specs->{scripts}/cut_target_domains.mti.pl \$scratchdir
cp $split_dir/\$input1 \$scratchdir

cd \$scratchdir


echo \"#sgejob run started on \$curhost at \$curtime\"
perl cut_target_domains.mti.pl < \$input1

set curtime=`date`
echo \"#sgejob run finished on \$curhost at \$curtime\"

rm -f \$scratchdir/\$input1 \$scratchdir/cut_target_domains.mti.pl
cd \$curdir
rmdir \$scratchdir\n" ;
         close($sgescript_fh) ;

         my $qsub_job_id = modtie::SGE::_clust_qsub({
            host => $in->{cluster}->{head_node},
            sgescript_fn => $sgescript_fn,
         }) ;

         while (1) {
            sleep $specs->{cluster}->{qstat_sleep} ;
            my $job_status= modtie::SGE::_clust_qstat({
               host => $in->{cluster}->{head_node},
               job_id => $qsub_job_id
            });
            if ($job_status) {last;}
         }

         modtie::SGE::_clust_merge_outs({
            script_fn => $sgescript_fn,
            out_fn => $temp_fn->{cutpaths_out},
            err_fn => $temp_fn->{cutpaths_err},
            job_id => $qsub_job_id,
            outdir => $sge_outdir,
            numjobs => $splits->{numjobs}
         }) ;
      }
   }


# Step 7. align candidate target - template domains: SALIGN
   if (!exists $in->{ali_done}) {
      print STDERR "Running alignments\n" ;
      ($temp_fh->{alipaths_out}, $temp_fn->{alipaths_out}) =
         tempfile("alipaths_out_XXXXX", SUFFIX => '.modtie') ;
      ($temp_fh->{alipaths_err}, $temp_fn->{alipaths_err}) =
         tempfile("alipaths_err_XXXXX", SUFFIX => '.modtie') ;

      if ($in->{cluster_fl} == 0) {
         $tcom = $specs->{scripts}."/salign_tmpl_targ_domains.mti.pl < $temp_fn->{alilist_out} 2>$temp_fn->{alipaths_err} >$temp_fn->{alipaths_out}";
         print STDERR "RUN: $tcom\n" ;
         system($tcom) ;
         if (-z $temp_fn->{alipaths_err}) {
            unlink $temp_fn->{alipaths_err};}

      } else {

         my $split_dir = tempdir("splits_salign_tmpl_targdomains.XXXXX") ;
         my $splits = modtie::SGE::_clust_split_ins({
            fn => $temp_fn->{alilist_out},
            dir => $split_dir,
            minlines => $specs->{cluster}->{minlines_salign},
            numjobs => $in->{cluster}->{numjobs}
         });

         my ($sgescript_fh, $sgescript_fn) = tempfile(
            "mt.salign_tmpl_targ_domains.XXXXX", SUFFIX => ".SGE.sh") ;
         my $sge_outdir = tempdir("SGEOUT.salign_tmpl_targ_domains.XXXXX") ;

         print {$sgescript_fh} "#!/bin/csh
#\$ -S /bin/csh
#\$ -cwd
#\$ -o $sge_outdir
#\$ -e $sge_outdir
#\$ -r y
$in->{cluster}->{nodespecs}
#\$ -p $in->{cluster}->{priority}
#\$ -t 1-$splits->{numjobs}

set tasks1=( $splits->{tasklist} )

set input1=\$tasks1[\$SGE_TASK_ID\]

set curdir=`pwd`
set curhost=`hostname`
set curtime=`date`
set scratchdir=/scratch/fred/\$input1.\$\$\

rm -rf \$scratchdir
mkdir -p \$scratchdir

cp $specs->{scripts}/salign_tmpl_targ_domains.mti.pl \$scratchdir
cp $split_dir/\$input1 \$scratchdir

cd \$scratchdir


echo \"#sgejob run started on \$curhost at \$curtime\"
perl salign_tmpl_targ_domains.mti.pl < \$input1

set curtime=`date`
echo \"#sgejob run finished on \$curhost at \$curtime\"

rm -f \$scratchdir/\$input1 \$scratchdir/salign_tmpl_targ_domains.mti.pl
cd \$curdir
rmdir \$scratchdir\n" ;
         close($sgescript_fh) ;

         my $qsub_job_id = modtie::SGE::_clust_qsub({
            host => $in->{cluster}->{head_node},
            sgescript_fn => $sgescript_fn,
         }) ;

         while (1) {
            sleep $specs->{cluster}->{qstat_sleep} ;
            my $job_status= modtie::SGE::_clust_qstat({
               host => $in->{cluster}->{head_node},
               job_id => $qsub_job_id
            });
            if ($job_status) {last;}
         }

         modtie::SGE::_clust_merge_outs({
            script_fn => $sgescript_fn,
            out_fn => $temp_fn->{alipaths_out},
            err_fn => $temp_fn->{alipaths_err},
            job_id => $qsub_job_id,
            outdir => $sge_outdir,
            numjobs => $splits->{numjobs}
         }) ;
      }

      if (exists $in->{partial_alipaths_err}) {
         if ($in->{partial_alipaths_err} =~ /\.gz$/) {
            my $t_com = $bins->{'zcat'}.
               " $in->{partial_alipaths_err} >> $temp_fn->{alipaths_err}" ;
            system($t_com) ;
         } else {
            system("cat $in->{partial_alipaths_err} >> ".
                   "$temp_fn->{alipaths_err}") ;
         }
      }
   } else {
      $temp_fn->{alipaths_err} = $in->{alipaths_err} ;
   }

# Step 8. Regenerate candidate list taking into account target/template
#         alignments that failed or scored poorly and
#
   if (!exists $in->{postali_candilist}) {

      if (!exists $temp_fn->{alipaths_err} ||
          !-s $temp_fn->{alipaths_err}) {
         croak "ERROR: postali_candilist requires alipaths_err file\n" ;
      }

      print STDERR "Re-running candidate generator\n" ;
      ($temp_fh->{postali_candilist_out}, $temp_fn->{postali_candilist_out}) =
         tempfile("postali_candilist_out_XXXXX", SUFFIX => '.modtie') ;
      close($temp_fh->{postali_candilist_out}) ;
      scoring_main({
         xset_fl => $in->{xset_fl},
         scoplevel => $in->{scoplevel},
         postali_candilist_fl => 1,
         seqset => $seqset,
         in_fn => {
            alipaths_err_fn => $temp_fn->{alipaths_err},
            seqid_setinfo_fn => $temp_fn->{seqid_setinfo},
            seqid_domainarch_fn => $temp_fn->{seqid_domainarch_out},
            seqdomains_fn => $temp_fn->{seqid_domains_out},
            model_pdbstart_fn => $temp_fn->{model_pdbstart_out},
            templates_fn => $in->{templates_fn}
         },
         out_fn => {
            postali_candilist_fn => $temp_fn->{postali_candilist_out},
         }
      }) ;

   } else {
      $temp_fn->{postali_candilist_out} = $in->{postali_candilist} ;
   }

# Step 9. Score candidate binary interactions
#
   if (!exists $in->{binscores}) {
      print STDERR "Assessing binary interactions\n" ;
      ($temp_fh->{binscores_out}, $temp_fn->{binscores_out}) =
         tempfile("binscores_out_XXXXX", SUFFIX => '.modtie') ;
      ($temp_fh->{binscores_err}, $temp_fn->{binscores_err}) =
         tempfile("binscores_err_XXXXX", SUFFIX => '.modtie') ;

      if ($in->{cluster_fl} == 0) {

         scoring_main({
            scoplevel => $in->{scoplevel},
            assess_fl => 1,
            in_fn => {
               postali_candilist_fn => $temp_fn->{postali_candilist_out},
               seqid_domainarch_fn => $temp_fn->{seqid_domainarch_out},
               model_pdbstart_fn => $temp_fn->{model_pdbstart_out},
               seqdomains_fn => $temp_fn->{seqid_domains_out}
            },
            out_fn => {
               interactions_fn => $temp_fn->{binscores_out},
            }
         }) ;

      } else {
         my $split_dir = tempdir("splits_binscores.XXXXX") ;
         my $splits = modtie::SGE::_clust_split_ins({
            fn => $temp_fn->{postali_candilist_out},
            dir => $split_dir,
            numjobs => $in->{cluster}->{numjobs}
         });

         my ($perlscript_fh, $perlscript_fn) =
            tempfile("mt.binscores.XXXXX", SUFFIX => ".mti.pl") ;
         my ($sgescript_fh, $sgescript_fn) =
            tempfile("mt.binscores.XXXXX", SUFFIX => ".SGE.sh") ;
         my $sge_outdir = tempdir("SGEOUT.binscores.XXXXX") ;

         print {$perlscript_fh} '#!/usr/local/bin/perl'."
use strict;
use warnings;
use lib \'".$specs->{libpath}."\';
use modtie qw/scoring_main/ ;

main() ;

sub main {

         scoring_main({
            scoplevel => \"".$in->{scoplevel}."\",
            assess_fl => 1,
            in_fn => {
               postali_candilist_fn => \$ARGV[0],
               seqid_domainarch_fn => \"".$temp_fn->{seqid_domainarch_out}."\",
               model_pdbstart_fn => \"".$temp_fn->{model_pdbstart_out}."\",
               seqdomains_fn => \"".$temp_fn->{seqid_domains_out}."\"
            },
         }) ;

}\n" ;

         print {$sgescript_fh} "#!/bin/csh
#\$ -S /bin/csh
#\$ -cwd
#\$ -o $sge_outdir
#\$ -e $sge_outdir
#\$ -r y
$in->{cluster}->{nodespecs}
#\$ -p $in->{cluster}->{priority}
#\$ -t 1-$splits->{numjobs}

set tasks1=( $splits->{tasklist} )

set input1=\$tasks1[\$SGE_TASK_ID\]

set curdir=`pwd`
set curhost=`hostname`
set curtime=`date`
set scratchdir=/scratch/fred/\$input1.\$\$\

rm -rf \$scratchdir
mkdir -p \$scratchdir

cp $perlscript_fn \$scratchdir
cp $split_dir/\$input1 \$scratchdir
cp $temp_fn->{seqid_domainarch_out} \$scratchdir
cp $temp_fn->{seqid_domains_out} \$scratchdir
cp $temp_fn->{model_pdbstart_out} \$scratchdir

cd \$scratchdir

echo \"#sgejob run started on \$curhost at \$curtime\"
perl $perlscript_fn \$input1

set curtime=`date`
echo \"#sgejob run finished on \$curhost at \$curtime\"

rm -f \$scratchdir/\$input1 \$scratchdir/$perlscript_fn \$scratchdir/$temp_fn->{seqid_domainarch_out} \$scratchdir/$temp_fn->{seqid_domains_out} \$scratchdir/$temp_fn->{model_pdbstart_out}
cd \$curdir
rmdir \$scratchdir\n" ;

         close($perlscript_fh) ;
         close($sgescript_fh) ;

         my $qsub_job_id = modtie::SGE::_clust_qsub({
            host => $in->{cluster}->{head_node},
            sgescript_fn => $sgescript_fn,
         }) ;

         while (1) {
            sleep $specs->{cluster}->{qstat_sleep} ;
            my $job_status= modtie::SGE::_clust_qstat({
               host => $in->{cluster}->{head_node},
               job_id => $qsub_job_id
            });
            if ($job_status) {last;}
         }

         my $binscore_headers = scoring_main({
            get_assess_headers => 1 });
         modtie::SGE::_clust_merge_outs({
            script_fn => $sgescript_fn,
            headers => $binscore_headers,
            out_fn => $temp_fn->{binscores_out},
            err_fn => $temp_fn->{binscores_err},
            job_id => $qsub_job_id,
            outdir => $sge_outdir,
            numjobs => $splits->{numjobs}
         }) ;
      }
      $temp_fn->{interactions_out} = $temp_fn->{binscores_out} ;
   } else {
      $temp_fn->{interactions_out} = $in->{binscores} ;
   }

# Step 10. Generate candidate list of higher-order complexes
#
   if (!exists $in->{complexes_done}) {
      print STDERR "Building higher-order complexes\n" ;
      modtie::complexes::list_complexes({
         in_fn => {
            interactions_fn => $temp_fn->{interactions_out},
            seqid_domainarch_fn => $temp_fn->{seqid_domainarch_out},
         },
      }) ;
   }
}


=head2 webget_modbase_files()

   Title:       webget_modbase_files()
   Function:    Retrive/setup necessary files from MODBASE to process a
                 given set of sequences.
   Args:        $_->{seqid_set} = hash list of seqids
   Returns:     $_->{error_fl} = if an error is encountered

=cut

sub webget_modbase_files {

   my $in = shift ;

   my $specs ;
   if (exists $in->{specs}) { $specs = $in->{specs} ; }
   else { $specs = set_modtie_specs() ; }

# make a temporary directory to play in
   my $curdir = getcwd;
   my ($tempdir) = File::Temp::tempdir(CLEANUP => 1) ;
   my $cookie_com = '';
   if (exists $in->{cookie_fn} && $in->{cookie_fn} ne '') {
      if (!-s $in->{cookie_fn}) {
         print STDERR "WARNING: cookie file not found: ".$in->{cookie_fn}."\n";}
      else {
         my $basename = basename($in->{cookie_fn}) ;
         my $new_loc = $tempdir.'/'.$basename ;
         modtie::pibase::safe_copy($in->{cookie_fn}, $new_loc) ;
         $cookie_com = " --load-cookies $basename ";
      }
   }
   chdir $tempdir ;

   my $model_info = {};
   foreach my $seq_id (sort keys %{$in->{seq_id}}) {

# If a run name was specified, use it to restrict MODBASE results
      my $dataset_spec = '';
      my $cur_setnum = $in->{seq_id}->{$seq_id};
      #0=intra; 1=inter-set1; 2=inter-set2; 3=inter-sets1+2
      if ($cur_setnum == 0 && exists $in->{run_specs}->{run}) {
         $dataset_spec = '&dataset='.$in->{run_specs}->{run};
      } elsif ($cur_setnum > 0 && $cur_setnum < 3 &&
               exists $in->{run_specs}->{"run".$cur_setnum}) {
         $dataset_spec = '&dataset='.$in->{run_specs}->{"run".$cur_setnum};
      }

# Retrieve MODBASE alignment XML file
      sleep($specs->{wget_specs}->{waittime}) ;
      my $modbase_ali_url = $specs->{modbase_ali_url}.$seq_id.$dataset_spec ;
      my $ali_out_fn = "$seq_id.modbase_ali" ;
      {
         my $wget_done = 0 ;
         my $num_wget_tries = $specs->{wget_specs}->{retries} ;
         my $wget_com = "wget --connect-timeout=".
            $specs->{wget_specs}->{connect_timeout}." ".
            "\'$modbase_ali_url\' ".$cookie_com."-O $ali_out_fn 2>/dev/null";
         print STDERR $wget_com."\n";
         if (!$wget_done && $num_wget_tries > 0) {
            system($wget_com) ;
            $num_wget_tries-- ;
            if (-s $ali_out_fn) {$wget_done = 1;}
         }
      }

# Make sure file completely downloaded
      if (!-s $ali_out_fn) {
         print STDERR "ERROR: did not retrieve alignments for seq $seq_id\n";
         next; }

      system("tail -1 $ali_out_fn > /tmp/lastline.$$.out") ;
      my $lastline = `cat /tmp/lastline.$$.out`; chomp $lastline ;
      system("rm /tmp/lastline.$$.out") ;
      if ($lastline !~ /<\/files/) {
         print STDERR "ERROR: incomplete download of alignments ".
                      "for seq $seq_id\n" ;
         next;
      }

# Parse ALI XML file and extract to individual .ali files
      my $ali_info = parse_modbase_ali_xml({
         seq_id => $seq_id,
         fn => $ali_out_fn,
         dest_dir => $specs->{local_modbase_ali_dir}
      }) ;
      unlink $ali_out_fn ;

# Retrieve MODBASE model PDB XML file
      sleep($specs->{wget_specs}->{waittime}) ;
      my $modbase_pdb_url = $specs->{modbase_pdb_url}.$seq_id.$dataset_spec ;
      my $models_out_fn = "$seq_id.modbase_pdb" ;
      {
         my $wget_done = 0 ;
         my $num_wget_tries = $specs->{wget_specs}->{retries} ;
         my $wget_com = "wget --connect-timeout=".
            $specs->{wget_specs}->{connect_timeout}." ".
            "\'$modbase_pdb_url\' ".$cookie_com."-O $models_out_fn 2>/dev/null" ;
         print STDERR $wget_com."\n";
         if (!$wget_done && $num_wget_tries > 0) {
            system($wget_com) ;
            $num_wget_tries-- ;
            if (-s $models_out_fn) {$wget_done = 1;}
         }
      }

# Make sure file completely downloaded
      if (!-s $models_out_fn) {
         print STDERR "ERROR: did not retrieve models for seq $seq_id\n";
         next; }

      system("tail -1 $models_out_fn > /tmp/lastline.$$.out") ;
      $lastline = `cat /tmp/lastline.$$.out`; chomp $lastline ;
      system("rm /tmp/lastline.$$.out") ;
      if ($lastline !~ /<\/files/) {
         print STDERR "ERROR: incomplete download of models for seq $seq_id\n";
         next;
      }

# Parse PDB XML file and extract to individual .ali files
      my $pdb_info = parse_modbase_pdb_xml({
         seq_id => $seq_id,
         fn => $models_out_fn,
         dest_dir => $specs->{local_modbase_models_dir}
      }) ;
      unlink $models_out_fn ;


# Map from alignment ID to PDB info entries
      my $aln2pdb_info_ind = {};
      foreach my $j ( 0.. $#{$pdb_info}) {
         $aln2pdb_info_ind->{$pdb_info->[$j]->{"MODPIPE ALIGNMENT ID"}} = $j; }

# Load information into the model_entries format required by internal routines
      foreach my $j ( 0 .. $#{$ali_info}) {
         my $pdbinfo_ind = $aln2pdb_info_ind->{$ali_info->[$j]->{align_id}} ;

         push @{$model_info->{model_id}}, $ali_info->[$j]->{model_id} ;
         push @{$model_info->{seq_id}}, $ali_info->[$j]->{seq_id} ;
         push @{$model_info->{align_id}}, $ali_info->[$j]->{align_id} ;
         push @{$model_info->{seq_ident_global}},
            $ali_info->[$j]->{sequence_identity} ;
         push @{$model_info->{evalue}}, $ali_info->[$j]->{evalue} ;
         push @{$model_info->{target_length}},
            $pdb_info->[$pdbinfo_ind]->{"TARGET LENGTH"} ;
         push @{$model_info->{target_beg}}, $ali_info->[$j]->{target_begin} ;
         if ($ali_info->[$j]->{software} eq 'ModPipe1.0') {
            push @{$model_info->{modpipe_style}}, 'old' ;
         } else {
            push @{$model_info->{modpipe_style}}, 'new' ;
         }
         push @{$model_info->{run}}, $ali_info->[$j]->{dataset} ;
         push @{$model_info->{pdb_code}}, $ali_info->[$j]->{template_pdb_code};
#         push @{$model_info->{zdope}},
#            $pdb_info->[$pdbinfo_ind]->{zDOPE};
#         push @{$model_info->{modelscore}},
#            $pdb_info->[$pdbinfo_ind]->{"MODEL SCORE"};
#         push @{$model_info->{mpqs}},
#            $pdb_info->[$pdbinfo_ind]->{"ModPipe Quality Score"};


# MFU criteria from modbase webpage:
# http://modbase.compbio.ucsf.edu/modbase/modbase_help.html#filter
         my $model_score = $pdb_info->[$pdbinfo_ind]->{"MODEL SCORE"}  ;
         my $evalue = $ali_info->[$j]->{evalue} ;

         my $cur_mfu = '' ;
         if ($evalue <= 0.0001) {$cur_mfu .= 'F';}
         if ($model_score >= 0.7) {$cur_mfu .= 'M';}
         if ($cur_mfu eq '') {$cur_mfu = ' ';}
         push @{$model_info->{mfu}}, $cur_mfu ;
      }
   }

   chdir $curdir ;
   unlink $tempdir ;

   my $nummodels_read = $#{$model_info->{model_id}} ;

   return $model_info ;

}


=head2 read_seqset_xset()

   Title:       read_seqset_xset()
   Function:    Builds a hash pointing from seqid to protein set
   Args:        $_->{seqid_set1} = name of file with list of seqid in set 1
                $_->{seqid_set2} = name of file with list of seqid in set 2
   Returns:     $_->{union}->{SEQID} = 1 (set 1)
                                       2 (set 2)
                                       3 (set 1 and 2)

=cut

sub read_seqset_xset {

   my $in = shift;

   my $seqset;
   $seqset->{1} = readin_1col_arr({fn => $in->{seqid_set1}}) ;
   $seqset->{2} = readin_1col_arr({fn => $in->{seqid_set2}}) ;

   foreach my $setno(1 ..2) {
      foreach my $tseqid (@{$seqset->{$setno}}) {
         $seqset->{union}->{$tseqid} += $setno ; #1,2,or 3
      }
   }

   return $seqset ;

}


=head2 getstdinput()

   Title:       getstdinput()
   Function:    Retrieves input from ARGV, STDIN, or files depending on run mode
   Args:        $_->{mode} = MODTIE run mode
   Returns:     $_ = INPUT information to pass to actual run routine

=cut

sub getstdinput {

   my $in = shift ;
   my $input ;
   if ($in->{mode} eq 'model_2_domains') {
      $input = _getstdinput_model_2_domains({ARGV => $in->{ARGV}}) ;
   } elsif ($in->{mode} eq 'format_modbaseoutput') {
      $input = _getstdinput_format_modbaseoutput({ARGV => $in->{ARGV}}) ;
   } elsif ($in->{mode} eq 'calc_p_b_prior') {
      $input = _getstdinput_calc_p_b_prior({ARGV => $in->{ARGV}}) ;
   } elsif ($in->{mode} eq 'assess_yeast_results') {
      $input = _getstdinput_assess_yeast_results({ARGV => $in->{ARGV}}) ;
   } elsif ($in->{mode} eq 'runmodtie_roc') {
      $input = _getstdinput_runmodtie_roc({ARGV => $in->{ARGV}}) ;
   } elsif ($in->{mode} eq 'runmodtie_benchmark') {
      $input = _getstdinput_runmodtie_benchmark({ARGV => $in->{ARGV}}) ;
   } elsif ($in->{mode} eq 'buildpotential_count') {
      $input = _getstdinput_buildpotential_count({ARGV => $in->{ARGV}}) ;
   } elsif ($in->{mode} eq 'buildpotential_postcalc') {
      $input = _getstdinput_buildpotential_postcalc({ARGV => $in->{ARGV}}) ;
   } elsif ($in->{mode} eq 'seqid_2_domains') {
      $input = _getstdinput_seqid_2_domains({ARGV => $in->{ARGV}}) ;
   } elsif ($in->{mode} eq 'runmodtie_modbase') {
      $input = _getstdinput_runmodtie_modbase({ARGV => $in->{ARGV}}) ;
   } elsif ($in->{mode} eq 'runmodtie_scorecomplex') {
      $input = _getstdinput_runmodtie_scorecomplex({ARGV => $in->{ARGV}}) ;
   } elsif ($in->{mode} eq 'runmodtie_targetstrxs_template') {
      $input = _getstdinput_runmodtie_targetstrxs_template({ARGV => $in->{ARGV}}) ;
   } elsif ($in->{mode} eq 'seqid_2_domainarch') {
      $input = _getstdinput_seqid_2_domainarch({ARGV => $in->{ARGV}}) ;
   } elsif ($in->{mode} eq 'model_2_pdbstart') {
      $input = _getstdinput_model_2_pdbstart({ARGV => $in->{ARGV}}) ;
   }

   return $input ;

}


sub _getstdinput_format_modbaseoutput {

   my $in = shift ;
   my $options ;
   my $fields = {
      binary_fn => 0,
      complexes_fn => 0,
      out_fn => 0,
      project => 0,
   } ;

   my $j = 0 ;
   while ($j <= $#{$in->{ARGV}}) {
      my $type = $in->{ARGV}->[$j] ; $type =~ s/^-// ;
      my $val = $in->{ARGV}->[($j + 1)] ;
      if (!exists $fields->{$type}) {
         die "ERROR: option $type not recognized\n" ; }
      print STDERR "set $type: $val\n" ;
      $options->{$type} = $val ;
      $j += 2 ;
   }

   return $options ;

}


sub _getstdinput_assess_yeast_results {

   my $in = shift ;
   my $options ;
   my $fields = {
      binary_fn => 0,
      complexes_fn => 0,
      binary_fl => 0,
      complex_fl => 0,
      seqid_aaseq_arch_fn => 0,
   } ;

   my $j = 0 ;
   while ($j <= $#{$in->{ARGV}}) {
      my $type = $in->{ARGV}->[$j] ; $type =~ s/^-// ;
      my $val = $in->{ARGV}->[($j + 1)] ;
      if (!exists $fields->{$type}) {
         die "ERROR: option $type not recognized\n" ; }
      print STDERR "set $type: $val\n" ;
      $options->{$type} = $val ;
      $j += 2 ;
   }

   return $options ;

}


sub _getstdinput_calc_p_b_prior {

   my $in = shift ;
   my $options ;
   my $fields = {
      seqid_aaseq_arch => 0,
   } ;

   my $j = 0 ;
   while ($j <= $#{$in->{ARGV}}) {
      my $type = $in->{ARGV}->[$j] ; $type =~ s/^-// ;
      my $val = $in->{ARGV}->[($j + 1)] ;
      if (!exists $fields->{$type}) {
         die "ERROR: option $type not recognized\n" ; }
      print STDERR "set $type: $val\n" ;
      $options->{$type} = $val ;
      $j += 2 ;
   }

   return $options ;

}


sub _getstdinput_buildpotential_postcalc {

   my $in = shift;
   my $options;

   if (exists $in->{ARGV} && exists $in->{ARGV}->[0]) {
      $options->{outfn_prefix} = shift @{$in->{ARGV}} ; }

   my $specs ;
   if (exists $in->{specs}) { $specs = $in->{specs} ; }
   else { $specs = set_modtie_specs() ; }

   my $types;
   $types->{sc}= ['mm', 'ms', 'ss', 'all'] ;
   $types->{rad}= [4,6,8] ;
   $types->{inter}= ['inter', 'intra'] ;

   my $counts ;
   foreach my $sc_type (@{$types->{sc}}) {
      foreach my $radius (@{$types->{rad}}) {
         foreach my $restype1 (keys %{$specs->{gconst}->{aa3to1}}) {
            $counts->{$sc_type}->{$radius}->{n_P} = 0 ;
            foreach my $restype2 (keys %{$specs->{gconst}->{aa3to1}}) {

               $counts->{$sc_type}->{$radius}->{maxafp_ij}->{$restype1}->{$restype2} = 0 ;
               $counts->{$sc_type}->{$radius}->{nij_P}->{$restype1}->{$restype2} = 1 ;
               $counts->{$sc_type}->{$radius}->{dnij_P__n_P}->{$restype1}->{$restype2} = 1 ;
            }
         }
      }
   }

   while (my $line = <STDIN>) {
      if ($line =~ /^\#/ || $line =~ /^$/) {next;}
      chomp $line;

      my ($sc_type, $radius, $type, @vals)= split("\t", $line) ;

      if ($type eq 'n_P') {
         $counts->{$sc_type}->{$radius}->{$type} += $vals[0] ;
      } elsif ($type eq 'maxafp_ij') {
         my ($res1, $res2, $val) = @vals ;
         if ($val > $counts->{$sc_type}->{$radius}->{maxafp_ij}->{$res1}->{$res2}) { $counts->{$sc_type}->{$radius}->{maxafp_ij}->{$res1}->{$res2} = $val;}
      } else {
         my ($res1, $res2, $val) = @vals ;
         $counts->{$sc_type}->{$radius}->{$type}->{$res1}->{$res2} += $val ;
      }
   }

   $options->{counts} = $counts;
   return $options ;
}


sub _getstdinput_buildpotential_count {

   my $in = shift;
   my $options;

   my $entries;
   while (my $line = <STDIN>) {
      chomp $line;
      my @vals = split("\t", $line) ;
      if ($#vals == 0 ) {
         push @{$entries}, {sid => $vals[0]} ;
      } elsif ($#vals == 1) {
         push @{$entries}, {sid_1 => $vals[0], sid_2 => $vals[1]} ;
      }
   }

   $options->{entries} = $entries ;
   return $options ;
}


sub _getstdinput_runmodtie_targetstrxs_template {

   my $in = shift;
   my $options;

# format: jobid (target|template) pdb_fn domain_id domdef if necessary
# if domdef not provided (only optional for target: assumes entire pdb_fn

   my $complexes;
   while (my $line = <STDIN>) {
    if ($line !~ /^$/ && $line !~ /^\#/) {


      chomp $line;
      my @t = split(/\t/, $line) ;
      my $complex_id = shift @t ;
      my $type = shift @t ;
      my $domid = shift @t ;
      my $pdb_fn = shift @t ;

      $complexes->{$complex_id}->{domains}->{$domid}->{$type}->{pdb_fn} = $pdb_fn ;
      if ($type eq 'template') {
         $complexes->{$complex_id}->{template_pdbfn} = $pdb_fn ; }

      if ($#t > 0) {
         my $j = 0 ;
         while ($j <= $#t) {
            my $cur_start = $t[$j] ;
            my $cur_end = $t[($j + 1)] ;
            my $cur_chain = $t[($j + 2)] ;
            if (!defined $cur_start) {$cur_start = '';}
            if (!defined $cur_end) {$cur_end= '';}
            if (!defined $cur_chain) {$cur_chain= ' ';}
            push @{$complexes->{$complex_id}->{domdef}->{$domid}->{$type}->{start}}, $cur_start ;
            push @{$complexes->{$complex_id}->{domdef}->{$domid}->{$type}->{end}}, $cur_end ;
            push @{$complexes->{$complex_id}->{domdef}->{$domid}->{$type}->{chain}}, $cur_chain ;
            $j+=3 ;
         }
      }
    }
   }

   $options->{complexes} = $complexes ;
   return $options ;

}


sub _getstdinput_runmodtie_scorecomplex {

   my $in = shift;
   my $options;

   my $fields = {
     pot => 0,
     out_scores_fn => 0,
   } ;

   my $j = 0 ;
   while ($j <= $#{$in->{ARGV}}) {
      my $type = $in->{ARGV}->[$j] ; $type =~ s/^-// ;
      my $val = $in->{ARGV}->[($j + 1)] ;
      if (!exists $fields->{$type}) {
         die "ERROR: option $type not recognized\n" ; }
      print STDERR "set $type: $val\n" ;
      $options->{$type} = $val ;
      $j += 2 ;
   }

   my $complexes ;
   while (my $line = <STDIN>) {
      chomp $line;
      my @t  = split(/\t/, $line) ;
      my $curfn = shift @t ;
      my $j = 0 ;
      while ($j <= $#t) {
         my $cur_domid = $t[$j] ;
         my $cur_start = $t[($j + 1)] ;
         my $cur_end = $t[($j + 2)] ;
         my $cur_chain = $t[($j + 3)] ;
         push @{$complexes->{$curfn}->{$cur_domid}->{start}}, $cur_start ;
         push @{$complexes->{$curfn}->{$cur_domid}->{end}}, $cur_end ;
         push @{$complexes->{$curfn}->{$cur_domid}->{chain}}, $cur_chain ;
         $j+=4 ;
      }
   }

   $options->{complexes} = $complexes ;
   return $options ;
}


sub _getstdinput_runmodtie_roc {

   my $in = shift ;
   my $options ;
   my $fields = {
      rocfile_prefix => 0,
   } ;

   my $j = 0 ;
   while ($j <= $#{$in->{ARGV}}) {
      my $type = $in->{ARGV}->[$j] ; $type =~ s/^-// ;
      my $val = $in->{ARGV}->[($j + 1)] ;
      if (!exists $fields->{$type}) {
         die "ERROR: option $type not recognized\n" ; }
      print STDERR "set $type: $val\n" ;
      $options->{$type} = $val ;
      $j += 2 ;
   }

   my $z_scores ;
   my $potinfo ;
   while (my $line = <STDIN>) {
      if ($line =~ /^INTERFACE/) {
         my @t = split(/\t/, $line) ;
         my $inttype = 'inter';
         if ($t[3] =~ /\.intra\./) { $inttype = 'intra' ; }
         my $pottype = join('.', ($t[4], $inttype, @t[5,6]));

         $potinfo->{$pottype} = $t[3] ;

         my $t_zsc = $t[8] ;
         my $t_zsc_bg = $t[18] ;

# changed the structure of the zscores so a positive and
# its associated negatives can be removed together
         push @{$z_scores->{$pottype}->{p}},  [$t_zsc] ;
         push @{$z_scores->{$pottype}->{n}}, [split(', ', $t_zsc_bg)] ;
      }
   }

   $options->{potinfo} = $potinfo ;
   $options->{z_scores} = $z_scores;

   return $options ;

}


sub _getstdinput_runmodtie_benchmark {

   my $in = shift;
   my $options;

   my $fields = {
      scores_fn => 0,
   } ;

   my $j = 0 ;
   while ($j <= $#{$in->{ARGV}}) {
      my $type = $in->{ARGV}->[$j] ; $type =~ s/^-// ;
      my $val = $in->{ARGV}->[($j + 1)] ;
      if (!exists $fields->{$type}) {
         die "ERROR: option $type not recognized\n" ; }
      print STDERR "set $type: $val\n" ;
      $options->{$type} = $val ;
      $j += 2 ;
   }


   my $complexes ;
   while (my $line = <STDIN>) {
      chomp $line;
      my ($bdp_id, @sids) = split(/\t/, $line) ;
      push @{$complexes}, {
         bdp_id => $bdp_id,
         subset_ids => \@sids
      } ;
   }
   $options->{complexes} = $complexes ;

   return $options;

}


sub _getstdinput_runmodtie_modbase {

   my $in = shift;

   my $fields = {
      cluster_fl => 0,
      scoplevel => 0,
      candilist_fl => 0,
      local_modbase_models_dir => 0,
      local_modbase_ali_dir => 0,
      prepare_models_only => 0,
      model_pdbstart => 0,
      postali_candilist => 0,
      modweb_cookie => 0,
      model_list => 0,
      seqid_setinfo => 0,
      model_domains => 0,
      seqid_domains => 0,
      seqid_domainarch => 0,
      partial_alipaths_err => 0,
      alipaths_err => 0,
      copyover_done => 0,
      cuts_done => 0,
      ali_done => 0,
      bin_done => 0,
      alilist => 0,
      cutlist => 0,
      binscores => 0,
      templates_fn => 0,
      complexes_done => 0,
      seqid_set => 0,
      seqid_set1 => 0,
      seqid_set2 => 0,
      run => 0,
      run1 => 0,
      run2 => 0,
   } ;

   my $options;
   my $j = 0 ;
   while ($j <= $#{$in->{ARGV}}) {
      my $type = $in->{ARGV}->[$j] ; $type =~ s/^-// ;
      my $val = $in->{ARGV}->[($j + 1)] ;
      if (!exists $fields->{$type}) {
         die "ERROR: option $type not recognized\n" ; }
      print STDERR "set $type: $val\n" ;
      $options->{$type} = $val ;
      $j += 2 ;
   }
   print STDERR "\n";

   return $options ;
}


=head2 readin_1col_arr()

   Title:       readin_1col_arr()
   Function:    Reads in a single column file into an array
   Args:        $_->{fn} = filename
   Returns:     $_->[] = array of values in file

=cut

sub readin_1col_arr {

   my $in = shift;
   my $fn = $in->{fn} ;

   open (FH, $in->{fn}) ;
   my $results; 
   while (my $line = <FH>) {
      chomp $line;
      push @{$results}, $line ;
   }

   return $results;

}


sub _getstdinput_model_2_pdbstart {

   my $in = shift;
   my $models ;

# Format specified in assign_seqid_domains()
   while (my $line = <STDIN>) {
      chomp $line;
      if ($line =~ /^#/) {next;}
      my ($seq_id, $start_resno, $end_resno, $model_id, $target_beg, $subset_length, $seq_length, $tmpl_subset_id, $tmpl_class, $tmpl_subset_size, $subset_class, $subset_seq_ident, $mfu, $evalue, $scop_quality, $segment_no, $num_segments, $modpipe_style, $baseali_dir, $basemodels_dir) = split(/\t/, $line) ;

      $models->{$seq_id}->{$model_id} = {
         modpipe_style => $modpipe_style,
         basemodels_dir => $basemodels_dir
      };
   }

   return {models => $models, display => 1, datareturn => 0};
}


sub _getstdinput_seqid_2_domainarch {

   my $in = shift;
   my $doms ;

   my $j = 0 ;
# Format specified in assign_seqid_domains()
   while (my $line = <STDIN>) {
      chomp $line;
      if ($line =~ /^#/) {next;}
      my ($seq_id, $start_resno, $end_resno, $model_id, $target_beg, $subset_length, $seq_length, $tmpl_subset_id, $tmpl_class, $tmpl_subset_size, $subset_class, $subset_seq_ident, $mfu, $evalue, $scop_quality, $segment_no, $num_segments, $modpipe_style, $baseali_dir, $basemodels_dir) = split(/\t/, $line) ;

      push @{$doms->{seq_id}}, $seq_id ;
      push @{$doms->{model_id}}, $model_id ;
      push @{$doms->{target_beg}}, $target_beg ;
      push @{$doms->{start_resno}}, $start_resno;
      push @{$doms->{end_resno}}, $end_resno;
      push @{$doms->{susbset_id}}, $tmpl_subset_id;
      push @{$doms->{subset_class}}, $subset_class;
      push @{$doms->{subset_length}}, $subset_length;
      push @{$doms->{segment_no}}, $segment_no;
      push @{$doms->{num_segments}}, $num_segments;
      push @{$doms->{seq_length}}, $seq_length;
      push @{$doms->{scop_quality}}, $scop_quality;

      push @{$doms->{modpipe_style}}, $modpipe_style;
      push @{$doms->{baseali_dir}}, $baseali_dir;
      push @{$doms->{basemodels_dir}}, $basemodels_dir;

      push @{$doms->{seqid2doms}->{$seq_id}}, $j ;
      $j++ ;
   }

   return {doms => $doms, display => 1, datareturn => 0};
}


sub _getstdinput_seqid_2_domains {

   my $in = shift;
   my $doms ;

   my $readfh ;
   if (exists $in->{fn}) {
      open($readfh, $in->{fn}) ;
   } else {
      open($readfh, '-') ;
   }

# Format specified in assign_model_domains()
   my $j = 0 ;
   while (my $line = <$readfh>) {
      chomp $line;
      if ($line =~ /^#/) {next;}
      my ($seq_id, $target_length, $target_beg, $model_id, $modpipe_style, $baseali_dir, $basemodels_dir, $seq_ident_global, $evalue, $mfu, $seg_start, $seg_stop, $subset_id, $class, $numres_target_seg, $numres_tmpl_seg, $numres_aln_seg, $numres_match_seg, $numres_target_full, $numres_tmpl_full, $numres_aln_full, $numres_match_full, $subset_size) = split(/\t/, $line) ;

      push @{$doms->{seq_id}}, $seq_id ;
      push @{$doms->{target_length}}, $target_length;
      push @{$doms->{target_beg}}, $target_beg;
      push @{$doms->{model_id}}, $model_id ;

      push @{$doms->{modpipe_style}}, $modpipe_style;
      push @{$doms->{baseali_dir}}, $baseali_dir;
      push @{$doms->{basemodels_dir}}, $basemodels_dir;

      push @{$doms->{seq_ident_global}}, $seq_ident_global;
      push @{$doms->{evalue}}, $evalue;
      push @{$doms->{mfu}}, $mfu;

      push @{$doms->{seg_start}}, $seg_start ;
      push @{$doms->{seg_end}}, $seg_stop ;
      push @{$doms->{sid}}, $subset_id;
      push @{$doms->{class}}, $class;

      push @{$doms->{numres_targ_seg}}, $numres_target_seg ;
      push @{$doms->{numres_tmpl_seg}}, $numres_tmpl_seg ;
      push @{$doms->{numres_aln_seg}}, $numres_aln_seg ;
      push @{$doms->{numres_match_seg}}, $numres_match_seg ;

      push @{$doms->{numres_targ_full}}, $numres_target_full ;
      push @{$doms->{numres_tmpl_full}}, $numres_tmpl_full ;
      push @{$doms->{numres_aln_full}}, $numres_aln_full ;
      push @{$doms->{numres_match_full}}, $numres_match_full ;

      push @{$doms->{subset_size}}, $subset_size;

      push @{$doms->{seqid2doms}->{$seq_id}}, $j ;
      $j++ ;
   }
   close($readfh) ;

   return {doms => $doms, display => 1, datareturn => 0};
}



sub _getstdinput_model_2_domains {

   my $j = 0 ;
   my $models = read_model_list({ fn => '-' }); #readin from STDIN

   return {models => $models, display => 1, datareturn => 0};
}


=head2 model_2_domains()

   Title:       model_2_domains()
   Function:    Assign SCOP domains to models based on alignment
   Args:        $_->{display} = option to display domains to STDOUT
                $_->{datareturn} = option to send back domain assignments as
                                   variable

                $_->{models}->{align_path}->{RUN} = alignment directory for RUN
                $_->{models}->{model_path}->{RUN} = model PDB directory for RUN
                $_->{models}->{modpipe_style}->{RUN} = modpipe style (old/new)

   Returns:     $_->{model_id} = {      #if $in->{datareturn} specified
                  run
                  seq_ident_global
                  evalue
                  target_length
                  target_beg
                  mfu
                  seg_start
                  seg_stop
                  subset_id
                  class
                  seg_targ
                  seg_tmpl
                  seg_numres_aln
                  seg_numres_match
                  targ_subset_id
                  tmpl_subset_id
                  numres_aln
                  numres_match
                  subset_size
                }

   Displays:    model assignment information as tab-delimited values
                 of hash above.

=cut

sub model_2_domains {

   my $in = shift;
   my $models = $in->{models} ;
   my $model_entries = $models->{model_entries} ;
   my $pdb2models = $models->{pdb2models} ;

   my ($display, $datareturn) ;
   if (exists $in->{display}) { $display = $in->{display}; }
   if (exists $in->{datareturn}) { $datareturn = $in->{datareturn}; }

   my $pb_pl = assign_model_domains_preload() ;

   my $model2dom ;
# wrap up models in batches of same pdb_id then send to assign_model_domains
   foreach my $pdb_id (keys %{$pdb2models}) {
      my $model_arg ;

      map { push @{$model_arg->{$_}},
            @{$model_entries->{$_}}[@{$pdb2models->{$pdb_id}}] }
            (keys %{$model_entries}) ;

      assign_model_domains({
         results => $model2dom,
         models => $model_arg,
         display => $display,
         datareturn => $datareturn,
         pb_preload => $pb_pl
      }) ;
   }

   return $model2dom ;

}


=head2 seqid_2_domains()

   Title:       seqid_2_domains()
   Function:    Assign SCOP domains to target sequences
   Args:        $_->{display} = option to display domains to STDOUT
                $_->{datareturn} = option to send back domain assignments as
                                   variable

                $_->{doms}->{align_path}->{RUN} = alignment directory for RUN
                $_->{doms}->{model_path}->{RUN} = model directory for RUN
                $_->{doms}->{modpipe_style}->{RUN} = modpipe style for RUN

   Returns:     $_->{seq_id} = {      #if $in->{datareturn} specified
                  start
                  end
                  model_id
                  run
                  targetbeg
                  numres_tmpl
                  cur_length
                  subset_id
                  tmpl_class
                  subset_size
                  class
                  subset_seqident
                  mfu
                  evalue
                  scop_quality
                  frag_no
                  num_frags
                }

   Displays:    domain assignment information as tab-delimited values
                 of hash above.

=cut

sub seqid_2_domains {
   my $in = shift;

   my $seq2dom ;

   my ($display, $datareturn) ;
   if (exists $in->{display}) { $display = $in->{display}; }
   if (exists $in->{datareturn}) { $datareturn = $in->{datareturn}; }

   assign_seqid_domains({
      results => $seq2dom,
      doms => $in->{doms},
      display => $display,
      datareturn => $datareturn,
   }) ;

   return $seq2dom ;
}


=head2 seqid_2_domainarch()

   Title:       seqid_2_domainarch()
   Function:    Calculate SCOP domain architecture for each target sequence
   Args:        $_->{display} = option to display domains to STDOUT
                $_->{datareturn} = option to send back domain assignments as
                                   variable

                $_->{doms}->{align_path}->{RUN} = alignment directory for RUN
                $_->{doms}->{model_path}->{RUN} = model directory for RUN
                $_->{doms}->{modpipe_style}->{RUN} = modpipe style for RUN

   Returns:     $_->{seq_id} = {      #if $in->{datareturn} specified
                  seq_length
                  cover_length
                  num_domains
                  domarch_ss
                  domarch_s = 
                  domarch_f=[domain_no:class:fragno:fraglen:domlen:scop_quality]
                }

   Displays:    domain architecture information as tab-delimited values
                 of hash above.

=cut

sub seqid_2_domainarch {
   my $in = shift;

   my $seq2domarch ;

   my ($display, $datareturn) ;
   if (exists $in->{display}) { $display = $in->{display}; }
   if (exists $in->{datareturn}) { $datareturn = $in->{datareturn}; }

   assign_seqid_domainarch({
      results => $seq2domarch,
      doms => $in->{doms},
      display => $display,
      datareturn => $datareturn,
   }) ;

   return $seq2domarch ;
}



=head2 model_2_pdbstart()

   Title:       model_2_pdbstart()
   Function:    Get residue number of first PDB ATOM record for each model
   Args:        $_->{display} = option to display domains to STDOUT
                $_->{datareturn} = option to send back domain assignments as
                                   variable

                $_->{doms}->{align_path}->{RUN} = alignment directory for RUN
                $_->{doms}->{model_path}->{RUN} = model directory for RUN
                $_->{doms}->{modpipe_style}->{RUN} = modpipe style for RUN

   Returns:     $_->{model_id} = {      #if $in->{datareturn} specified
                  pdbstart_resno
                }

   Displays:    domain architecture information as tab-delimited values
                 of hash above.

=cut

sub model_2_pdbstart {
   my $in = shift;

   my $model2pdbstart = {};

   my ($display, $datareturn) ;
   if (exists $in->{display}) { $display = $in->{display}; }
   if (exists $in->{datareturn}) { $datareturn = $in->{datareturn}; }

   calc_model_pdbstart({
      results => $model2pdbstart,
      models => $in->{models},
      display => $display,
      datareturn => $datareturn,
   }) ;

   return $model2pdbstart;
}

=head2 assign_model_domains_preload()

   Title:       assign_model_domains_preload()
   Function:    Preload PIBASE data necessary for model domain assignment
   Args:        NONE

   Tables:      pibase.bdp_files (tod)
                pibase.bdp_residues_tables (tod)
                pibase.subsets_residues (tod)
                pibase.subsets (tod)

   Returns:     $->{subset_2_class}->{subset_id} = subsets_class
                $->{pdb_2_subsets}->{PDB_ID} = [subset_id_1, _2, ...]
                $->{pdbid_2_subres}->{PDB_ID} = file containing subsets residues
                $->{pdbid_2_residues}->{PDB_ID} = file containing all residues
                $->{raw_pdb}->{BDP_ID} = 1 if PDB, 0 otherwise (eg PQS)
                $->{pdb_2_bdp_id}->{PDB_ID} = BDP_ID

=cut

sub assign_model_domains_preload {

   my $specs = set_modtie_specs() ;

   my $raw_pdb ;
   my $pdb_2_bdp_id ;
   {
      my ($t_bdp_id, $t_pdb_id, $t_raw_pdb) =
         modtie::pibase::rawselect_tod(
            "SELECT bdp_id, pdb_id, raw_pdb FROM bdp_files");
      foreach my $j ( 0 .. $#{$t_bdp_id}) {
         if ($t_raw_pdb->[$j] == 1) {
            $pdb_2_bdp_id->{$t_pdb_id->[$j]} = $t_bdp_id->[$j] ;
            $raw_pdb->{$t_bdp_id->[$j]} = $t_pdb_id->[$j] ; } }
   }

# filepath to residue listing tables
   my $pdbid_2_residues ;
   {
      my ($t_bdp_id, $t_sourcefile) =
         modtie::pibase::rawselect_tod("SELECT bdp_id, source_file ".
                               "FROM bdp_residues_tables") ;
      foreach my $j ( 0 .. $#{$t_bdp_id}) {
         if (exists $raw_pdb->{$t_bdp_id->[$j]}) {
            $pdbid_2_residues->{$raw_pdb->{$t_bdp_id->[$j]}} =
               $t_sourcefile->[$j] ; } }
      replace_basedir_inplace({
         hash => $pdbid_2_residues,
         old_dir => $specs->{pibase_specs}->{old_metatod_dir}."/bdp_residues",
         new_dir => $specs->{pibase_specs}->{metatod_dir}."/bdp_residues"
      }) ;
   }

# subsets residues tables: residue subset assignments
   my $pdbid_2_subres ;
   {
      my ($t_bdp_id, $t_sourcefile) =
         modtie::pibase::rawselect_tod("SELECT bdp_id, source_file ".
                               "FROM subsets_residues_tables") ;
      foreach my $j ( 0 .. $#{$t_bdp_id}) {
         if (exists $raw_pdb->{$t_bdp_id->[$j]}) {
            $pdbid_2_subres->{$raw_pdb->{$t_bdp_id->[$j]}} =
               $t_sourcefile->[$j] ; } }
      replace_basedir_inplace({
         hash => $pdbid_2_subres,
         old_dir => $specs->{pibase_specs}->{old_metatod_dir}.
            "/subsets_residues",
         new_dir => $specs->{pibase_specs}->{metatod_dir}."/subsets_residues"
      }) ;
   }

# subset classes
   my $pdb_2_subsets;
   my $subsets_2_class ;
   {
      my ($t_bdp_id, $t_subset_id, $t_class) =
         modtie::pibase::rawselect_tod(
            "SELECT bdp_id, subset_id, class FROM subsets") ;
      foreach my $j ( 0 .. $#{$t_bdp_id}) {
         if (exists $raw_pdb->{$t_bdp_id->[$j]}) {
            $subsets_2_class->{$t_subset_id->[$j]} = $t_class->[$j] ;
            push @{$pdb_2_subsets->{$raw_pdb->{$t_bdp_id->[$j]}}},
               $t_subset_id->[$j]; } }
   }

   my $pb = {
      subsets_2_class => $subsets_2_class,
      pdb_2_subsets => $pdb_2_subsets,
      pdbid_2_subres => $pdbid_2_subres,
      pdbid_2_residues => $pdbid_2_residues,
      raw_pdb => $raw_pdb,
      pdb_2_bdp_id => $pdb_2_bdp_id
   } ;

   return $pb ;
}



=head2 assign_model_domains()

   Title:       assign_model_domains()
   Function:    Calculate domain assignments for each model by parsing
                 MODBASE alignment
   Args:        $_->{display} = option to display domains to STDOUT
                $_->{datareturn} = option to send back domain assignments as
                                   variable
                $_->{models}->{model_id} = [] list of model identifiers
                $_->{pb_preload} = preload of PIBASE
                $_->{results} = reference to send back results

   Tables:      pibase.bdp_residues (tod)
                pibase.bdp_residues (tod)
                pibase.subsets_residues (tod)

   Returns:     (populates $_->{results}
                ->{MODEL_ID} = {
                  run
                  seq_ident_global
                  evalue
                  target_length
                  target_beg
                  mfu
                  seg_start
                  seg_stop
                  subset_id
                  class
                  seg_targ
                  seg_tmpl
                  seg_numres_aln
                  seg_numres_match
                  targ_subset_id
                  tmpl_subset_id
                  numres_aln
                  numres_match
                  subset_size
                }

=cut

sub assign_model_domains {

   my $subinfo ;
   ($subinfo->{package},
    $subinfo->{filename},
    $subinfo->{line},
    $subinfo->{name},
    $subinfo->{has_args},
    $subinfo->{iwantarray})= caller(0);
   my $in = shift ;
   my $specs ;
   if (exists $in->{specs}) { $specs = $in->{specs} ; }
   else { $specs = set_modtie_specs() ; }

   my ($display, $datareturn) ;
   if (exists $in->{display}) { $display = $in->{display}; }
   if (exists $in->{datareturn}) { $datareturn = $in->{datareturn}; }

   my $models = $in->{models};
   my $pb_pl = $in->{pb_preload} ;
   my $results = $in->{results} ;

   my $subsets_2_class = $pb_pl->{subsets_2_class} ;
   my $pdbid_2_subres = $pb_pl->{pdbid_2_subres} ;
   my $pdbid_2_residues = $pb_pl->{pdbid_2_residues} ;
#100323_1057   my $pdb_2_subsets = $pb_pl->{pdb_2_subsets} ;
#100323_1057   my $raw_pdb = $pb_pl->{raw_pdb} ;
#100323_1057   my $pdb_2_bdp_id = $pb_pl->{pdb_2_bdp_id} ;

# given a set of ModBase models (pre-ordered by template_pdb to cut down on pdb read in) assigns domain boundaries

   my ($chain_id, $chain_no, $resno, $resno_ser) ;
   my $resnochain_2_resnoser;

   my $last_pdb= 'XXXX' ;
   my ($resno_2_ser,  $resnoser_2_subset, $subset_size) ;
   my ($cr_serial_2_resno, $cr_resno_2_serial, $cr_cid_2_no, $cr_cno_2_id) ;
   my ($raw_serial_2_resno, $raw_resno_2_serial) ;
   my ($weird_res_raw, $weird_res_ser) ;
   my ($weirdo_res, $weirdo_res_ser) ;

   foreach my $j (0 .. $#{$models->{model_id}}) {
      my $seq_id = $models->{seq_id}->[$j] ;
      my $align_id = $models->{align_id}->[$j] ;
      my $model_id = $models->{model_id}->[$j] ;
      my $seq_ident_global = $models->{seq_ident_global}->[$j] ;
      my $evalue = $models->{evalue}->[$j] ;
      my $target_length = $models->{target_length}->[$j] ;
      my $target_beg = $models->{target_beg}->[$j] ;
      my $modpipe_style = $models->{modpipe_style}->[$j] ;
      my $mfu = $models->{mfu}->[$j] ;
      my $pdb_code = $models->{pdb_code}->[$j] ;
      my $baseali_dir = $models->{baseali_dir}->[$j] ;
      my $basemodels_dir = $models->{basemodels_dir}->[$j] ;

      my $ali_fn = $baseali_dir."/".substr($seq_id,0,3).
                   "/$seq_id/alignments/$align_id.ali";

# Parse meta information from alignment
      my $ali_info = get_ali_meta($ali_fn) ;
      if (exists $ali_info->{error}) {
         print STDERR  "ERROR ($seq_id,$align_id,$model_id): $ali_info->{error} ($subinfo->{name} ".__LINE__.")\n" ;
         next ;
      }

      my $tmpl_pdb = $ali_info->{template}->{pdb_id} ;
      my $tmpl_start = $ali_info->{template}->{start_resno} ;
      my $tmpl_end = $ali_info->{template}->{end_resno} ;
      my $tmpl_chain = $ali_info->{template}->{start_chain} ;

#If the residue listing or residue subset assignment for the template pdb does not exist, goto next model

      if (!exists $pdbid_2_residues->{$tmpl_pdb}) {
         print STDERR "ERROR ($seq_id,$align_id,$model_id): no residues ".
            "table for pdb $tmpl_pdb ($subinfo->{name} ".__LINE__.")\n" ;
         next ;
      }

      if (!exists $pdbid_2_subres->{$tmpl_pdb}){
         print STDERR "ERROR ($seq_id,$align_id,$model_id): no subsets ".
         "residues table for pdb $tmpl_pdb ($subinfo->{name} ".__LINE__.")\n" ;
         next ;
      }


# when parsing MODELLER alignments, act as though the MSE residue does not
# exist, i.e. replace it with a '-', so as though the target has no alignment
# otherwise throws serial numbering off
#
# so do we need cr_ series as though  MSE was actually parsed?
#  NO, throws bdp_residues serial assignments off


# If the current template pdb is not the same as the previous models, load template info

      if ($tmpl_pdb ne $last_pdb) {

# reinitialize template residue info
         $resno_2_ser = {} ;
         $resnoser_2_subset = {} ;
         $subset_size = {} ;
         $cr_serial_2_resno = {} ;
         $cr_resno_2_serial = {} ;
         $cr_cid_2_no = {} ;
         $cr_cno_2_id = {} ;

# Read in residue, serial number mapping from bdp_residues files
         {
            my ($resno, $resno_ser, $chain_id, $chain_no) =
               modtie::pibase::rawselect_metatod($pdbid_2_residues->{$tmpl_pdb},
                  "SELECT resno, resno_serial, chain_id, chain_no ".
                  "FROM $pdbid_2_residues->{$tmpl_pdb}") ;
            foreach my $j (0 .. $#{$resno}) {
               my $sig = $resno->[$j]."\n".$chain_id->[$j] ;
               my $sig_ser = $resno_ser->[$j]."\n".$chain_no->[$j] ;
               $cr_serial_2_resno->{$sig_ser} = $sig ;
               $cr_resno_2_serial->{$sig} = $sig_ser ;
               $cr_cid_2_no->{$chain_id->[$j]} = $chain_no->[$j] ;
               $cr_cno_2_id->{$chain_no->[$j]} = $chain_id->[$j] ;

               $resno_2_ser->{$resno->[$j]."\n".$chain_id->[$j]} =
                  $resno_ser->[$j] ;
            }
         }

#Read in residue subset_id mapping from subsets_residues files
         {
            my ($resno, $resno_ser, $chain_id, $chain_no, $subset_id) =
               modtie::pibase::rawselect_metatod($pdbid_2_subres->{$tmpl_pdb},
                  "SELECT resno, resno_serial, chain_id, chain_no, subset_id ".
                  "FROM $pdbid_2_subres->{$tmpl_pdb}") ;
            foreach my $j ( 0 .. $#{$resno}) {
               if ($subset_id->[$j] =~ /SCOP/) {
                  my $ressig = $resno_ser->[$j]."\n".$chain_id->[$j] ;
                  $resnoser_2_subset->{$ressig} = $subset_id->[$j] ;
                  $subset_size->{$subset_id->[$j]}++ ;
               }
            }
         }

#Determine MSE positions in the template pdb (get_weirdo_resser())
         $raw_serial_2_resno = {};
         $raw_resno_2_serial = {};
# in theory required for chain no shifts when a chain has only one residue which turns uot to be a MSE HETATM - doubtful - fuckit for now
         {
            my $result_weirdo = modtie::pibase::get_weirdo_resser({
               pdb_fn => $specs->{pdb_dir}.'/'.substr($tmpl_pdb,1,2).
                         '/pdb'.$tmpl_pdb.'.ent.gz'}) ;

            foreach my $j (0 .. $#{$result_weirdo->{resno}}) {
               my $sig = $result_weirdo->{resno}->[$j]."\n".
                         $result_weirdo->{chain_id}->[$j] ;
               my $sig_ser = $result_weirdo->{resno_serial}->[$j]."\n".
                             $result_weirdo->{chain_no}->[$j] ;

               $raw_serial_2_resno->{$sig_ser} = $sig ;
               $raw_resno_2_serial->{$sig} = $sig_ser ;
            }
            $weird_res_raw = $result_weirdo->{weird_res_raw} ;
            $weird_res_ser = $result_weirdo->{weird_res_ser} ;
         }
      }

# If there is only one chain in the template pdb, reset the template chain extracted from the alignment to the actual chain in the pdb file. Fixes cases were RCSB has changed the chain id from ' ' to, in most cases, 'A', after the ModPipe run.

      my $singlechain = 0;
      my $numchains = keys %{$cr_cid_2_no} ;
      if ($numchains  == 1) {
         $singlechain = 1 ; }

      if ((!defined $cr_cid_2_no->{$tmpl_chain}) && $singlechain) {
         ($tmpl_chain) =  (keys %{$cr_cid_2_no}); }

# Determine template target residue equivalence: calcula_resequiv()
      my $modpipe_newstyle_orderswitch = 0 ;
      if ($modpipe_style eq 'new') {
         $modpipe_newstyle_orderswitch = 1 ; }

      my $resmap = modeldomains_calcula_resequiv({
         ali_fn => $ali_fn,
         tmpl_chain => $tmpl_chain,
         start => $ali_info->{target}->{start_resno},
         end => $ali_info->{target}->{end_resno},
         ser2resno => $cr_serial_2_resno,
         resno2ser => $cr_resno_2_serial,
         cid2no => $cr_cid_2_no,
         cno2id => $cr_cno_2_id,
         raw_resno2ser => $raw_resno_2_serial,
         raw_resser2no => $raw_serial_2_resno,
         weird_res_raw => $weird_res_raw,
         weird_res_ser => $weird_res_ser,
         modpipe_newstyle_orderswitch => $modpipe_newstyle_orderswitch
      }) ;

#If calcula_resequiv() returns an error, print it to STDOUT and goto next model

      if (exists $resmap->{error}) {
         print STDERR "ERROR ($seq_id,$align_id,$model_id): $resmap->{error} ".
            "($subinfo->{name} ".__LINE__.")\n";
         next ;
      }

      my $assigns ;
      my @assigns ;

# Reset template start and end residue to that returned by calcula_resequiv() - incase of a start or end MSE residue. If they are undefined, say so to STDOUT, then abort to next model.

      $tmpl_start = $resmap->{tmpl_start} ;
      $tmpl_end = $resmap->{tmpl_end} ;

      my $tmpl_startser = $resmap->{tmpl_start_ser} ;
      my $tmpl_endser = $resmap->{tmpl_end_ser} ;

      if (!defined $tmpl_startser ) {
         print STDERR "ERROR ($seq_id,$align_id,$model_id): template start ($tmpl_start : $tmpl_chain) not identified in residues table for pdb $tmpl_pdb ($pdbid_2_residues->{$tmpl_pdb}) ($subinfo->{name} ".__LINE__.")\n";
         next ;
      }

      if (!defined $tmpl_endser ) {
         print STDERR "ERROR ($seq_id,$align_id,$model_id): template end ($tmpl_end : $tmpl_chain) not identified in residues table for pdb $tmpl_pdb ($pdbid_2_residues->{$tmpl_pdb}) ($subinfo->{name} ".__LINE__.")\n" ;
         next ;
      }

      my @match_segs ;
      my $last_sid = '';
      my $last_target_resno = '';

# Iterate over the template residues from the alignment start to end residue

      foreach my $j ( $tmpl_startser .. $tmpl_endser ) {

# Determine aligned target residue. If undefined (ie gap), set curcover to 0 and reset to the last aligned target residue

         my $cur_target_resno = $resmap->{serialch}->[0]->{$j."\n".$tmpl_chain};
         my $cur_resmatch=$resmap->{serialch_match}->[0]->{$j."\n".$tmpl_chain};

         my $curcover = 1;

#YO: dont know which residue listing instance to attach residue names to, but need to determine subset_id-local sequence identity.

         if (!defined $cur_target_resno) {
            $curcover = 0 ;
            $cur_target_resno = $last_target_resno ; }

         my $sersig = $j."\n".$tmpl_chain ;
         my $cur_sid ;

#Determine domain assignment of current template residue

         if (exists $resnoser_2_subset->{$sersig}) {

            $cur_sid = $resnoser_2_subset->{$sersig} ;

#If not already on it, push the subset_id onto the list of assignments

            if (!exists $assigns->{$cur_sid}) {
               $assigns->{$cur_sid}->{numres_match} = 0 ;
               $assigns->{$cur_sid}->{numres_aln} = 0 ;
               $assigns->{$cur_sid}->{targ} = 0 ;
               $assigns->{$cur_sid}->{tmpl} = 0 ;
               push @assigns, $cur_sid ;
            }

#Increment the number of aligned template residues belonging to the current subset

            $assigns->{$cur_sid}->{tmpl}++ ;

#If curcover is 1, increment the number of aligned target residues.

            if ($curcover) {
               $assigns->{$cur_sid}->{targ}++ ;
               $assigns->{$cur_sid}->{numres_aln}++ ;
               if ($cur_resmatch) {
                  $assigns->{$cur_sid}->{numres_match}++ ; }
            }
            if (!exists $assigns->{$cur_sid}->{targ}) {
               $assigns->{$cur_sid}->{targ} = 0 ; }
         } else {
            $cur_sid = '' ;
         }

# If the current subset_id is different than the last one: set the stop residue on the previous matching segment, and set the start residue on a new one.

         if ($last_sid ne $cur_sid) {
            if ($last_sid ne '') {
               $match_segs[$#match_segs]->{stop} =$last_target_resno ; }
            if ($cur_sid ne '') {
               push @match_segs, {subset_id => $cur_sid,
                  start => $cur_target_resno} ;
               $match_segs[$#match_segs]->{numres_aln} = 0 ;
               $match_segs[$#match_segs]->{numres_match} = 0 ;
               $match_segs[$#match_segs]->{targ} = 0 ;
               $match_segs[$#match_segs]->{tmpl} = 0 ;
            }
         }

         if ($cur_sid ne '') {
            $match_segs[$#match_segs]->{tmpl}++ ;
            if ($curcover) {
               $match_segs[$#match_segs]->{numres_aln}++ ;
               $match_segs[$#match_segs]->{targ}++ ;
               if ($cur_resmatch) {
                  $match_segs[$#match_segs]->{numres_match}++ ; }
            }
         }

# Set the last residue and chain to the current one.

         if (defined $cur_target_resno) {
            $last_target_resno = $cur_target_resno ; }
         $last_sid = $cur_sid ;

      }

# Set the stop residue on the last matching segment.

      if ($#match_segs >= 0 ) {
       if (!exists $match_segs[$#match_segs]->{stop}) {
          $match_segs[$#match_segs]->{stop} = $last_target_resno; } }

# If no subsets were matched, print so to STDOUT.

      if (($#assigns < 0) && $display) {
         my @outvals = ($seq_id, $model_id,
                        $ali_info->{target}->{start_resno},
                        $ali_info->{target}->{end_resno}, 'NONE') ;
         print STDERR "NOTE: ".join("\t", @outvals)."\n" ;
      }

# Iterate over matching segments and display to STDOUT.

      foreach my $k (0 .. $#match_segs) {
         my $subset_id = $match_segs[$k]->{subset_id} ;

         if ($display) {
            my @outvals = (
               $seq_id,
               $target_length,
               $target_beg,
               $model_id,
               $modpipe_style,
               $baseali_dir,
               $basemodels_dir,
               $seq_ident_global,
               $evalue,
               $mfu,
               $match_segs[$k]->{start},
               $match_segs[$k]->{stop},
               $subset_id,
               $subsets_2_class->{$subset_id},
               $match_segs[$k]->{targ},
               $match_segs[$k]->{tmpl},
               $match_segs[$k]->{numres_aln},
               $match_segs[$k]->{numres_match},
               $assigns->{$subset_id}->{targ},
               $assigns->{$subset_id}->{tmpl},
               $assigns->{$subset_id}->{numres_aln},
               $assigns->{$subset_id}->{numres_match},
               $subset_size->{$subset_id},
            ) ;
            print join("\t", @outvals)."\n" ;
         }

         if ($datareturn) {
            push @{$results->{$model_id}}, {
               modpipe_style=> $modpipe_style,
               baseali_dir => $baseali_dir,
               basemodels_dir => $basemodels_dir,
               seq_ident_global => $seq_ident_global,
               evalue => $evalue,
               target_length => $target_length,
               target_beg => $target_beg,
               mfu => $mfu,
               seg_start => $match_segs[$k]->{start},
               seg_stop => $match_segs[$k]->{stop},
               subset_id => $subset_id,
               class => $subsets_2_class->{$subset_id},
               seg_targ => $match_segs[$k]->{targ},
               seg_tmpl => $match_segs[$k]->{tmpl},
               seg_numres_aln => $match_segs[$k]->{numres_aln},
               seg_numres_match => $match_segs[$k]->{numres_match},
               targ_subset_id => $assigns->{$subset_id}->{targ},
               tmpl_subset_id => $assigns->{$subset_id}->{tmpl},
               numres_aln => $assigns->{$subset_id}->{numres_aln},
               numres_match => $assigns->{$subset_id}->{numres_match},
               subset_size => $subset_size->{$subset_id},
            } ;
         }
      }

      $last_pdb = $ali_info->{template}->{pdb_id}  ;
   }

}



=head2 get_ali_meta()

   Title:       get_ali_meta()
   Function:    Extract meta information from MODBASE target-template alignment
   Args:        $_ = alignment filename
   Returns:     ->{target} = {
                  start_resno
                  end_resno
                }

                ->{template} = {
                  pdb_id
                  start_resno
                  start_chain
                  end_resno
                  end_chain
                }

=cut

sub get_ali_meta {

   my $ali_fn = shift ;

   my $target ;
   my $template ;

   open (ALI_FN, $ali_fn) || return {error => "couldn't open $ali_fn"} ;
   while (my $line = <ALI_FN>) {
      chomp $line ;
      if ($line =~ /^structure/) {
         my @t = split(/:/, $line) ;
         $template->{pdb_id} = $t[1] ;
         $t[2] =~ s/ //g ;
         $template->{start_resno} = $t[2] ;
         $t[3] =~ s/ //g ; if ($t[3] eq '') {$t[3] = ' '; }
         $template->{start_chain} = $t[3] ;
         $t[4] =~ s/ //g ;
         $template->{end_resno} = $t[4] ;
         $t[5] =~ s/ //g ; if ($t[5] eq '') {$t[5] = ' '; }
         $template->{end_chain} = $t[5] ;

      } elsif ($line =~ /^sequence/) {
         my @t = split(/:/, $line) ;
         $t[2] =~ s/ //g ; $t[4] =~ s/ //g ;

         if ($t[2] eq '') { #modpipe_style=old
            my $string = $t[6] ;
            ($target->{start_resno}, $target->{end_resno}) =
               ($string =~ /.* ([0-9]+)-([0-9]+)/) ;

         } else { #modpipe_style=new

            $target->{start_resno} = $t[2] ;
            $target->{end_resno} = $t[4] ;
         }
      }

   }
   close(ALI_FN) ;

   my $results ;
   $results->{target} = $target ;
   $results->{template} = $template ;

   return $results;

}


=head2 modeldomains_calcula_resequiv()

   Title:       modeldomains_calcula_resequiv()
   Function:    combine serial residue number/positioning from
                 get_resequiv_serial with residue_info()

   Args:        $_->{ser2resno}
                $_->{resno2ser}
                $_->{cid2no}
                $_->{cno2id}
                $_->{raw_resno2ser}
                $_->{raw_resser2no}
                $_->{weird_res_ser}
                $_->{weird_res_raw}
                $_->{modpipe_newstyle_orderswitch}
                $_->{ali_fn}
                $_->{tmpl_chain}
                $_->{start}
                $_->{end}

   Returns:     ->{real}->[0]->{resno sequence 0} = aligned resno sequence 1
                ->{real}->[1]->{resno sequence 1} = aligned resno sequence 0

                ->{serial}->[0|1]->{serial resno 0|1} = aligned resno seq 1|0
                ->{serialch}->[0|1]->{serial resno 0|1."\n".chain_id} =
                  aligned resno seq 1|0

                ->{serialch_match}->[0|1]->{serial resno 0|1."\n".chain_id} =
                  1 if identical; 0 if not

                ->{tmpl_end_ser}   =    template end serial residue number
                ->{tmpl_start_ser} =    template start serial residue number
                ->{tmpl_end} =          template end residue number
                ->{tmpl_start}          template start residue number

=cut

sub modeldomains_calcula_resequiv {
   my $params = shift ;

   my $serial_2_resno = $params->{ser2resno} ;
   my $resno_2_serial = $params->{resno2ser};
   my $chainid_2_no = $params->{cid2no} ;
   my $chainno_2_id = $params->{cno2id} ;

   my $raw_resno_2_serial = $params->{raw_resno2ser}; #MSE (MODELLER) numbering
   my $raw_serial_2_resno = $params->{raw_resser2no}; #MSE (MODELLER) numbering
   my $weird_res_ser = $params->{weird_res_ser} ; #serial residue number
   my $weird_res_raw = $params->{weird_res_raw} ; #absolute residue number

   my $modpipe_newstyle_orderswitch = 0 ;
   if (exists $params->{modpipe_newstyle_orderswitch}) {
      $modpipe_newstyle_orderswitch = $params->{modpipe_newstyle_orderswitch} ;}
   my $re_alipars = modtie::pibase::parse_modeller_ali({
      ali_fn => $params->{ali_fn},
      modpipe_newstyle_orderswitch => $modpipe_newstyle_orderswitch
   }) ;

   $re_alipars->{chain_start}->[0] = $params->{tmpl_chain} ;
   $re_alipars->{chain_end}->[0] = $params->{tmpl_chain} ;

# ModBase PDB changes: catch those PDB files where the file was changed from
#   chain ids of blank to chain ids of 'A'
#

# target start/end
   $re_alipars->{resno_start}->[1] = $params->{start} ;
   $re_alipars->{resno_end}->[1] = $params->{end} ;

# fixed this by using get_weridres(): what happens when the first template residue in the alignment is a weirdo residue ; e.g. alignment b07c08caf0008296b6e3d1d052185cfa tempalte end (110 : A)

   my $t_startres = $re_alipars->{resno_start}->[0]."\n".
                    $re_alipars->{chain_start}->[0];

   my $t_endres = $re_alipars->{resno_end}->[0]."\n".
                  $re_alipars->{chain_end}->[0];

#   my $t ="initbounds: $t_startres -- $t_endres" ; $t =~ s/\n/:/g;
#   print STDERR $t."\n" ;

   if (exists $weird_res_raw->{$t_startres}) {
      my $t_ser = $raw_resno_2_serial->{$t_startres} ;
      my ($t_res, $t_ch) = split(/\n/, $t_ser) ;
      my $start_inc = 1 ;
      my $stillweird = 1;
      while ($stillweird && ($start_inc <= 5)) {
         if (!exists $weird_res_ser->{($start_inc + $t_res)."\n".$t_ch}) {
            $stillweird = 0 ;
            $t_startres=$raw_serial_2_resno->{($start_inc + $t_res)."\n".$t_ch};
            ($re_alipars->{resno_start}->[0], $re_alipars->{chain_start}->[0]) =
               split(/\n/, $t_startres) ;
         } else {
            $start_inc++ ;
         }
      }

      if ($stillweird) {
         return {error => "couldnt find a regular template start residue"} ;
      }
   }

   my $t_startserch = $chainid_2_no->{$re_alipars->{chain_start}->[0]} ;
   my $t_startserres = $resno_2_serial->{$re_alipars->{resno_start}->[0]."\n".$re_alipars->{chain_start}->[0]} ;

#   my $t ="   new start: $t_startres ($t_startserres)" ; $t =~ s/\n/:/g;
#   print STDERR $t."\n" ;

   if ((!defined $t_startserres) ||
       (!defined $t_startserch)) {
      return {error => "template start residue ($re_alipars->{resno_start}->[0]:$re_alipars->{chain_start}->[0]) not recognized"} ;
   }

   $t_startserres =~ s/\n.*// ;


   if (exists $weird_res_raw->{$t_endres}) {
      my $t_ser = $raw_resno_2_serial->{$t_endres} ;
      my ($t_res, $t_ch) = split(/\n/, $t_ser) ;
      my $end_dec = 1 ;
      my $stillweird = 1;
      while ($stillweird && ($end_dec <= 5)) {
         if (!exists $weird_res_ser->{($t_res - $end_dec)."\n".$t_ch}) {
            $stillweird = 0 ;
            $t_endres = $raw_serial_2_resno->{($t_res - $end_dec)."\n".$t_ch} ;
            ($re_alipars->{resno_end}->[0], $re_alipars->{chain_end}->[0]) =
               split(/\n/, $t_endres) ;
         } else {
            $end_dec++ ;
         }
      }
      if ($stillweird) {
         return {error => "couldnt find a regular template end residue"} ;
      }
   }

   my $t_endserch = $chainid_2_no->{$re_alipars->{chain_end}->[0]} ;
   my $t_endserres = $resno_2_serial->{$re_alipars->{resno_end}->[0]."\n".$re_alipars->{chain_end}->[0]} ;

#   my $t ="   new end: $t_endres ($t_endserres)" ; $t =~ s/\n/:/g;
#   print STDERR $t."\n" ;


   if ((!defined $t_endserres) ||
       (!defined $t_endserch)) {
      return {error => "template end residue ($re_alipars->{resno_end}->[0]:$re_alipars->{chain_end}->[0]) not recognized"} ;
   }

   $t_endserres =~ s/\n.*// ;

# adjust all serial numbers to the presence of weirdores

   my $weirdo_dec = 0 ;
   foreach my $weirdo (keys %{$weird_res_ser}) {
      my ($rn, $cn) = split(/\n/, $weirdo) ;
      if (($cn == $t_startserch) && ($rn <= $t_startserres)) {
         $weirdo_dec++ ;
#         print STDERR " weirdres $rn < $t_startserres: $weirdo_dec\n" ;
      }
   }

   my $weirdores ;
   foreach my $weirdo (keys %{$weird_res_ser}) {
      my ($rn, $cn) = split(/\n/, $weirdo) ;
      if (($cn == $t_startserch) && ($rn >= $t_startserres)) {
#         print STDERR "decrementing $rn by $weirdo_dec\n" ;
         $weirdores->{($rn - $weirdo_dec)."\n".$cn}++ ;
      }
   }

# iterate over alignment positions, if serial residue corresponds to MSE, replace with gap (i.e. remove equivalency), and decrement all further equivalences

   my $decrement = 0 ;

   my $done  = 1 ;
   my $maxlen = $re_alipars->{maxlength} - 1 ;
   my $pos = 0 ;
   while ($pos <= $maxlen) {
#      print STDERR "pos $pos: " ;
      if ( (exists $re_alipars->{alipos_2_serial}->[0]->[$pos]) &&
           (exists $re_alipars->{alipos_2_serial}->[1]->[$pos]) ) {

         my $cur_serialresno = $re_alipars->{alipos_2_serial}->[0]->[$pos] +
                               $t_startserres - 1 ;

         my $cur_chainno = $re_alipars->{alipos_2_chainno}->[0]->[$pos] +
                           $t_startserch - 1 ;

#	 if ($pos == 0) {
#	    print STDERR "\n   start tmpl res_ser: $re_alipars->{alipos_2_serial}->[0]->[$pos] + $t_startserres - 1\n" ;
#	    print STDERR "   going to check $cur_serialresno:$cur_chainno " ;
#	    print STDERR "    wr: ".$weirdores->{($cur_serialresno."\n".$cur_chainno)}."\n" ;
#	 }
#
#	 print STDERR "   got in tmpl $re_alipars->{alipos_2_serial}->[0]->[$pos] ($cur_serialresno) -- targ $re_alipars->{alipos_2_serial}->[1]->[$pos]\t" ;

         if ( exists $weirdores->{$cur_serialresno."\n".$cur_chainno} &&
              $weirdores->{$cur_serialresno."\n".$cur_chainno} > 0 ) {
#	    print STDERR "   weirdores def:\t" ;
            $re_alipars->{alipos_2_serial}->[0]->[$pos] = undef ;
            $weirdores->{$cur_serialresno."\n".$cur_chainno} = 0 ;
            $decrement++ ;
         } else {
#	    print STDERR "   dec by $decrement\t" ;
            $re_alipars->{alipos_2_serial}->[0]->[$pos] -= $decrement ;
         }

#	 print STDERR "\t(post-tmpl $re_alipars->{alipos_2_serial}->[0]->[$pos])" ;
      }
#      print STDERR "\n" ;
      $pos++ ;
   }


   my $resequiv ;
   my $resequiv_ser ;
   my $resequiv_serch ;
   my $resmatch_serch ;

# iterate over alignment positions

#   print STDERR "pos\ttmpl\ttarg\n" ;
   my $numgaps = 0 ;
   foreach my $pos ( 0 .. ($re_alipars->{maxlength} - 1)) {
#      print STDERR "$pos " ;
#      if ((!exists $re_alipars->{alipos_2_serial}->[0]->[$pos]) ||
#          (!defined $re_alipars->{alipos_2_serial}->[0]->[$pos])) {
#         print STDERR "   tmpl serial not exists/def\t" ; }
#      if (!exists $re_alipars->{alipos_2_serial}->[1]->[$pos]) {
#         print STDERR "   targ serial not exists" ; }

      if ( (exists $re_alipars->{alipos_2_serial}->[0]->[$pos]) &&
           (defined $re_alipars->{alipos_2_serial}->[0]->[$pos]) &&
           (exists $re_alipars->{alipos_2_serial}->[1]->[$pos]) ) {


         my $real_chainid ;
         my $real_resno ;
         my $ser_resno ;

         my $cur_serialresno = $re_alipars->{alipos_2_serial}->[0]->[$pos] +
                               $t_startserres - 1;
         my $cur_chainno = $re_alipars->{alipos_2_chainno}->[0]->[$pos] +
                           $t_startserch - 1 ;

         my $identical_match = 0 ;
         if ($re_alipars->{alipos_2_resna}->[0]->[$pos] eq
             $re_alipars->{alipos_2_resna}->[1]->[$pos] ) {
            $identical_match = 1 ; }

         $ser_resno->[0] = $cur_serialresno."\n".$cur_chainno ;

         $real_chainid->[0] = $chainno_2_id->{$cur_chainno} ;
         $real_resno->[0] =
            $serial_2_resno->{$cur_serialresno."\n".$cur_chainno} ;

#         print STDERR "pos is $pos, numgaps is $numgaps, resno start is ".
#                      $re_alipars->{resno_start}->[1]."\n" ;
         $real_resno->[1] = $re_alipars->{resno_start}->[1] + $pos - $numgaps;

         if (!defined $real_resno->[0]) {
            return {error => "real resno 0 is undefined: ".
                              "$cur_serialresno:$cur_chainno"} ; }
         if (!defined $real_resno->[1]) {
            return {error => "real resno 1 is undefined"} ; }

         $resequiv->[0]->{$real_resno->[0]} = $real_resno->[1] ;
         $resequiv->[1]->{$real_resno->[1]} = $real_resno->[0] ;


# 0 = tmpl strx, 1 = target seq

         $resequiv_ser->[0]->{$ser_resno->[0]} = $real_resno->[1] ;
         $resequiv_serch->[0]->{$cur_serialresno."\n".$real_chainid->[0]} =
            $real_resno->[1] ;
         $resmatch_serch->[0]->{$cur_serialresno."\n".$real_chainid->[0]} =
            $identical_match ;

         $resequiv_ser->[1]->{$real_resno->[1]} = $ser_resno->[0] ;

#         my $t = "$cur_serialresno ($real_resno->[0])\t$real_resno->[1]\n" ;
#	 $t =~ s/\n/:/g ; print STDERR $t ;
      } elsif (!exists $re_alipars->{alipos_2_serial}->[1]->[$pos]) {
         $numgaps++ ;
      }
#      print STDERR "\n" ;
   }

   return {
      real =>$resequiv,
      serial => $resequiv_ser,
      serialch => $resequiv_serch,
      serialch_match => $resmatch_serch,
      tmpl_end_ser => $t_endserres,
      tmpl_start_ser => $t_startserres,
      tmpl_end => $re_alipars->{resno_end}->[0],
      tmpl_start => $re_alipars->{resno_start}->[0]
   } ;

}


=head2 assign_seqid_domains()

   Title:       assign_seqid_domains()
   Function:    Combine model assigments to assign domains to sequence id
   Args:        $_->{display} = option to display domains to STDOUT
                $_->{datareturn} = option to send back domain assignments as
                                   variable
                $_->{results} = reference to send back results
                $_->{doms}->{align_path}->{RUN} = alignment directory for RUN
                $_->{doms}->{model_path}->{RUN} = model directory for RUN
                $_->{doms}->{modpipe_style}->{RUN} = modpipe style for RUN

   Returns:     populates $->{results} if in->{datareturn}:
                ->{SEQID} = [
                  start,
                  end,
                  model_id,
                  modpipe_style,
                  targetbeg,
                  numres_tmpl,
                  cur_length,
                  subset_id,
                  tmpl_class,
                  subset_size,
                  class,
                  subset_seqident,
                  mfu,
                  evalue,
                  scop_quality,
                  frag_no,
                  num_frags,
                ]

                if in->{display}: prints above fields tab-delimited to STDOUT

=cut

sub assign_seqid_domains {

   my $in = shift ;
   my $results = $in->{results} ;

   my ($display, $datareturn) ;
   if (exists $in->{display}) { $display = $in->{display}; }
   if (exists $in->{datareturn}) { $datareturn = $in->{datareturn}; }

   my $thresh ;
   $thresh->{len_cover} = 0.7;
   $thresh->{max_overlap} = 0.3;

   my $modelid_2_targetbeg ;

   foreach my $cur_seqid (sort keys %{$in->{doms}->{seqid2doms}}) {
      my @rows = @{$in->{doms}->{seqid2doms}->{$cur_seqid}} ;

      my ($c_modid, $c_seqident_g, $c_mfu, $c_start, $c_end, $c_sid, $c_class,
          $c_numres_targ_seg, $c_numres_tmpl_seg, $c_numres_aln_seg,
          $c_numres_match_seg,
          $c_numres_targ_ful, $c_numres_tmpl_ful, $c_numres_aln_ful,
          $c_numres_match_ful,
          $c_subset_size, $c_length, $c_targetbeg, $c_eval,
          $c_modpipe_style, $c_baseali_dir, $c_basemodels_dir,
          ) = (
            [@{$in->{doms}->{model_id}}[@rows]],
            [@{$in->{doms}->{seq_ident_global}}[@rows]],
            [@{$in->{doms}->{mfu}}[@rows]],
            [@{$in->{doms}->{seg_start}}[@rows]],
            [@{$in->{doms}->{seg_end}}[@rows]],
            [@{$in->{doms}->{sid}}[@rows]],
            [@{$in->{doms}->{class}}[@rows]],
            [@{$in->{doms}->{numres_targ_seg}}[@rows]],
            [@{$in->{doms}->{numres_tmpl_seg}}[@rows]],
            [@{$in->{doms}->{numres_aln_seg}}[@rows]],
            [@{$in->{doms}->{numres_match_seg}}[@rows]],
            [@{$in->{doms}->{numres_targ_full}}[@rows]],
            [@{$in->{doms}->{numres_tmpl_full}}[@rows]],
            [@{$in->{doms}->{numres_aln_full}}[@rows]],
            [@{$in->{doms}->{numres_match_full}}[@rows]],
            [@{$in->{doms}->{subset_size}}[@rows]],
            [@{$in->{doms}->{target_length}}[@rows]],
            [@{$in->{doms}->{target_beg}}[@rows]],
            [@{$in->{doms}->{evalue}}[@rows]],
            [@{$in->{doms}->{modpipe_style}}[@rows]],
            [@{$in->{doms}->{baseali_dir}}[@rows]],
            [@{$in->{doms}->{basemodels_dir}}[@rows]],
        ) ;

      my $cur_length = $c_length->[0] ;

      my $subsets ;
      foreach my $k ( 0 .. $#{$c_modid}) {
       $modelid_2_targetbeg->{$c_modid->[$k]} = $c_targetbeg->[$k] ;
       if ($c_numres_aln_ful->[$k] > 0 ) {
         my $sig = $c_modid->[$k].$c_sid->[$k] ;

         my $cur_seqident = $c_numres_match_ful->[$k] / $c_numres_aln_ful->[$k];
         if ($cur_seqident == 0 ) {
            next;}

         my $cur_class ;
         my $scop_quality ;
         if ($cur_seqident >= 0.3 && $c_mfu->[$k] =~ /M/) {
            $cur_class = $c_class->[$k];
            $scop_quality = 'fa' ;
         } elsif (($cur_seqident < 0.3 && $c_mfu->[$k] =~ /M/) ||
                  ($cur_seqident >= 0.3 && $c_mfu->[$k] =~ /F/) ) {
            ($cur_class) = ($c_class->[$k] =~ /^(.*\.[0-9]+\.[0-9]+)\.[0-9]+$/);
            $scop_quality = 'sf' ;
         } else {
            ($cur_class) = ($c_class->[$k] =~ /^(.*\.[0-9]+)\.[0-9]+\.[0-9]+$/);
            $scop_quality = 'fo' ;
         }

         if (!exists $subsets->{$sig}->{fullsub}) {
            my ($cur_eval) = $c_eval->[$k] ;

            $subsets->{$sig}->{fullsub} = {
               tmpl_class => $c_class->[$k],
               scop_quality  => $scop_quality,
               class => $cur_class,
               mfu => $c_mfu->[$k],
               model_id => $c_modid->[$k],
               evalue => $cur_eval,
               subset_id => $c_sid->[$k],

               modpipe_style => $c_modpipe_style->[$k],
               baseali_dir => $c_baseali_dir->[$k],
               basemodels_dir  => $c_basemodels_dir->[$k],

               numres_targ => $c_numres_targ_ful->[$k],
               numres_tmpl => $c_numres_tmpl_ful->[$k],
               numres_aln => $c_numres_aln_ful->[$k],
               numres_match => $c_numres_match_ful->[$k],
               subset_seq_ident =>
                  ($c_numres_match_ful->[$k] / $c_numres_aln_ful->[$k]),
               global_seq_ident => $c_seqident_g->[$k],
               subset_size => $c_subset_size->[$k],
               num_frags => 0,
            } ;
         }

         $subsets->{$sig}->{fullsub}->{num_frags}++ ;

         push @{$subsets->{$sig}->{segment}}, {
            start => $c_start->[$k],
            end => $c_end->[$k],
            numres_targ => $c_numres_targ_seg->[$k],
            numres_tmpl => $c_numres_tmpl_seg->[$k],
            numres_aln => $c_numres_aln_seg->[$k],
            numres_match => $c_numres_match_seg->[$k],
            frag_no => $subsets->{$sig}->{fullsub}->{num_frags},
         } ;
       }
      }

      my @deletethese ;
      foreach my $dom ( keys %{$subsets}) {
         if ($subsets->{$dom}->{fullsub}->{numres_tmpl} <
             $thresh->{len_cover} * $subsets->{$dom}->{fullsub}->{subset_size}){
            push @deletethese, $dom ; } }

      foreach my $k ( 0 .. $#deletethese) {
         delete $subsets->{$deletethese[$k]} ; }

      my $fullseq ; #array of residues in the full target sequence
      foreach my $sig ( keys %{$subsets}) {
         foreach my $seg ( @{$subsets->{$sig}->{segment}}) {
            foreach my $l ( $seg->{start} .. $seg->{end} ) {
               $fullseq->[$l]->{$sig}++ ; } } }

      my $overlap_single ;
      my $overlaps ;
      foreach my $k (1 .. $#{$fullseq}) {
         my @sid = (keys %{$fullseq->[$k]}) ;
         if ($#sid > 0 ) {
            foreach my $l (0 .. ($#sid - 1)) {
               foreach my $m (($l+1) .. $#sid ) {
                  $overlaps->{$sid[$l]."\n".$sid[$m]}++ ;
                  $overlap_single->{$sid[$l]}++ ;
                  $overlap_single->{$sid[$m]}++ ;
               }
            }
         }
      }

      my ($bs_overlaps) ;
      my ($verdict, $verdict_loss, $wins, $loss) ;
      foreach my $sidpair (keys %{$overlaps}) {
         my ($sid1, $sid2) = split(/\n/, $sidpair) ;
         my $rev_sidpair = $sid2."\n".$sid1 ;
         my $length_1 = $subsets->{$sid1}->{fullsub}->{numres_targ} ;
         my $length_2 = $subsets->{$sid2}->{fullsub}->{numres_targ} ;
         my $seqident_1= $subsets->{$sid1}->{fullsub}->{subset_seq_ident} ;
         my $seqident_2= $subsets->{$sid2}->{fullsub}->{subset_seq_ident} ;

         my $olen = $overlaps->{$sidpair} ;
         if (($olen >= $length_1 * $thresh->{max_overlap}) ||
             ($olen >= $length_2 * $thresh->{max_overlap})) {
            my $score ;
            $score->{seq} = log($seqident_1) - log($seqident_2) ;
            $score->{len} = ($length_1 - $length_2) / $cur_length ;

            my ($c_wins, $c_loss) ;
            if ($score->{seq} > 0 )     { $c_wins->{seq} = $sid1;
                                          $c_loss->{seq} = $sid2; }
            elsif ($score->{seq} < 0)   { $c_wins->{seq} = $sid2;
                                          $c_loss->{seq} = $sid1; }

            if ($score->{len} > 0 )     { $c_wins->{len} = $sid1;
                                          $c_loss->{len} = $sid2; }
            elsif ($score->{len} < 0)   { $c_wins->{len} = $sid2;
                                          $c_loss->{len} = $sid1; }

            if ((!exists $c_wins->{seq}) && (!exists $c_wins->{len} )) {
               $verdict->{$sidpair} = $sid1 ;
               $verdict->{$rev_sidpair} = $sid1 ;
               $verdict_loss->{$sidpair} = $sid2 ;
               $verdict_loss->{$rev_sidpair} = $sid2 ;
               $wins->{$sid1}++ ; $loss->{$sid1}++ ;
            } elsif (exists $c_wins->{seq} && !exists $c_wins->{len}) {
               $verdict->{$sidpair} = $c_wins->{seq} ;
               $verdict->{$rev_sidpair} = $c_wins->{seq} ;
               $verdict_loss->{$sidpair} = $c_loss->{seq} ;
               $verdict_loss->{$rev_sidpair} = $c_loss->{seq} ;
               $wins->{$c_wins->{seq}}++ ;
               $loss->{$c_loss->{seq}}++ ;
            } elsif (exists $c_wins->{len} && !exists $c_wins->{seq}) {
               $verdict->{$sidpair} = $c_wins->{len} ;
               $verdict->{$rev_sidpair} = $c_wins->{len} ;
               $verdict_loss->{$sidpair} = $c_loss->{len} ;
               $verdict_loss->{$rev_sidpair} = $c_loss->{len} ;
               $wins->{$c_wins->{len}}++ ;
               $loss->{$c_loss->{len}}++ ;
            } else {
               if (abs($score->{seq}) >= 0.5)  {
                  $verdict->{$sidpair} = $c_wins->{seq} ;
                  $verdict->{$rev_sidpair} = $c_wins->{seq} ;
                  $verdict_loss->{$sidpair} = $c_loss->{seq} ;
                  $verdict_loss->{$rev_sidpair} = $c_loss->{seq} ;
                  $wins->{$c_wins->{seq}}++ ;
                  $loss->{$c_loss->{seq}}++ ;
               } else {
                  if (abs($score->{len}) >= 0.1) {
                     $verdict->{$sidpair} = $c_wins->{len} ;
                     $verdict->{$rev_sidpair} = $c_wins->{len} ;
                     $verdict_loss->{$sidpair} = $c_loss->{len} ;
                     $verdict_loss->{$rev_sidpair} = $c_loss->{len} ;
                     $wins->{$c_wins->{len}}++ ;
                     $loss->{$c_loss->{len}}++ ;
                  } else {
                     $verdict->{$sidpair} = $c_wins->{seq} ;
                     $verdict->{$rev_sidpair} = $c_wins->{seq} ;
                     $verdict_loss->{$sidpair} = $c_loss->{seq} ;
                     $verdict_loss->{$rev_sidpair} = $c_loss->{seq} ;
                     $wins->{$c_wins->{seq}}++ ;
                     $loss->{$c_loss->{seq}}++ ;
                  }
               }
            }
         } else {
            push @{$bs_overlaps}, $sidpair ;
         }
      }

      foreach my $sidpair (@{$bs_overlaps}) {
         delete $overlaps->{$sidpair} ; }

      my $final ;
      foreach my $dom ( keys %{$subsets}) {
         if (!exists $loss->{$dom} ||
             (exists $wins->{$dom} && $wins->{$dom} > $loss->{$dom})) {
            $final->{$dom}++ ; } }

      my @finals = keys %{$final} ;
      my $final_deletes ;
      foreach my $k ( 0 .. ($#finals - 1)) {
         foreach my $l ( ($k + 1) .. $#finals) {
            if ((exists $overlaps->{$finals[$k]."\n".$finals[$l]} ||
                 exists $overlaps->{$finals[$l]."\n".$finals[$k]}) &&
                 !(exists $final_deletes->{$finals[$k]} ||
                   exists $final_deletes->{$finals[$l]})) {
              $final_deletes->{$verdict_loss->{$finals[$k]."\n".$finals[$l]}}++;
            } } }

      foreach my $dom (@finals) {
         if (!exists $final_deletes->{$dom}) {
            my $nowon = $subsets->{$dom}->{fullsub} ;
            foreach my $k (0 .. $#{$subsets->{$dom}->{segment}}) {
               my $seg = $subsets->{$dom}->{segment}->[$k] ;
               if ($display) {
                  my @outvals = (
                     $cur_seqid,
                     $seg->{start},
                     $seg->{end},
                     $nowon->{model_id},
                     $modelid_2_targetbeg->{$nowon->{model_id}},
                     $nowon->{numres_tmpl},
                     $cur_length,
                     $nowon->{subset_id},
                     $nowon->{tmpl_class},
                     $nowon->{subset_size},
                     $nowon->{class},
                     sprintf("%.2f", $nowon->{subset_seq_ident}),
                     $nowon->{mfu},
                     $nowon->{evalue},
                     $nowon->{scop_quality},
                     $seg->{frag_no},
                     $nowon->{num_frags},
                     $nowon->{modpipe_style},
                     $nowon->{baseali_dir},
                     $nowon->{basemodels_dir},
                  ) ;
                  print join("\t", @outvals)."\n" ;
               }

               if ($datareturn) {
                  push @{$results->{$cur_seqid}}, {
                     start => $seg->{start},
                     end => $seg->{end},
                     model_id => $nowon->{model_id},
                     targetbeg => $modelid_2_targetbeg->{$nowon->{model_id}},
                     numres_tmpl => $nowon->{numres_tmpl},
                     cur_length => $cur_length,
                     subset_id => $nowon->{subset_id},
                     tmpl_class => $nowon->{tmpl_class},
                     subset_size => $nowon->{subset_size},
                     class => $nowon->{class},
                     subset_seqident =>
                        sprintf("%.2f", $nowon->{subset_seq_ident}),
                     mfu => $nowon->{mfu},
                     evalue => $nowon->{evalue},
                     scop_quality => $nowon->{scop_quality},
                     frag_no => $seg->{frag_no},
                     num_frags => $nowon->{num_frags},
                     modpipe_style => $nowon->{modpipe_style},
                     baseali_dir => $nowon->{baseali_dir},
                     basemodels_dir => $nowon->{basemodels_dir},
                  } ;
               }
            }
         }
      }
   }
}


=head2 assign_seqid_domainarch()

   Title:       assign_seqid_domainarch()
   Function:    Calculate a SCOP domain architecture for each target sequence
   Args:        $_->{display} = option to display domains to STDOUT
                $_->{datareturn} = option to send back domain assignments as
                                   variable
                $_->{results} = reference to send back results
                $_->{doms}->{align_path}->{RUN} = alignment directory for RUN
                $_->{doms}->{model_path}->{RUN} = model directory for RUN
                $_->{doms}->{modpipe_style}->{RUN} = modpipe style for RUN

   Returns:     $_->{seq_id} = {      #if $in->{datareturn} specified
                  seq_length
                  cover_length
                  num_domains
                  domarch_ss
                  domarch_s = 
                  domarch_f=[domain_no:class:fragno:fraglen:domlen:scop_quality]
                }

                if in->{display}: prints above fields tab-delimited to STDOUT

=cut

sub assign_seqid_domainarch {

   my $in = shift;
   my $results = $in->{results} ;

   my ($display, $datareturn) ;
   if (exists $in->{display}) { $display = $in->{display}; }
   if (exists $in->{datareturn}) { $datareturn = $in->{datareturn}; }

# stoped using subseT_length fromseq_domains because that is only the number of target residues that were covered by the alignment - excludes gapped residues; now recalculates the subset_length by adding up fragments

   foreach my $cur_seqid (sort keys %{$in->{doms}->{seqid2doms}}) {
      my @rows = @{$in->{doms}->{seqid2doms}->{$cur_seqid}} ;
      my @order_resno = sort { $in->{doms}->{start_resno}->[$a] <=> $in->{doms}->{start_resno}->[$b] } @rows ;

      my ($model_id,
          $start_resno, $end_resno,
          $subset_id, $subset_class, $subset_length,
          $segment_no, $num_segments, $seq_length, $scop_quality) = (
            [@{$in->{doms}->{model_id}}[@order_resno]],
            [@{$in->{doms}->{start_resno}}[@order_resno]],
            [@{$in->{doms}->{end_resno}}[@order_resno]],
            [@{$in->{doms}->{susbset_id}}[@order_resno]],
            [@{$in->{doms}->{subset_class}}[@order_resno]],
            [@{$in->{doms}->{subset_length}}[@order_resno]],
            [@{$in->{doms}->{segment_no}}[@order_resno]],
            [@{$in->{doms}->{num_segments}}[@order_resno]],
            [@{$in->{doms}->{seq_length}}[@order_resno]],
            [@{$in->{doms}->{scop_quality}}[@order_resno]],
        ) ;
# must be ordered by start_resno, ala:
#         "SELECT model_id, start_resno, end_resno, tmpl_subset_id, subset_class, subset_length, segment_no, num_segments, seq_length, scop_quality FROM seq_domains where seq_id = \"$cur_seq_id\" order by start_resno") ;

      my $doms ;
      my @cells ;
      my $cur_domno = 0;
      foreach my $k ( 0 .. $#{$model_id}) {
         my $sid = $model_id->[$k]."\n".$subset_id->[$k] ;
         if (!exists $doms->{$sid}) {
            $cur_domno++ ;
            my $t_scopqual ;
            if ((!defined $scop_quality->[$k]) ||
                ($scop_quality->[$k] eq ' ') ||
                ($scop_quality->[$k] eq '') ) {
               $t_scopqual = 'u' ;
            } else {
               $t_scopqual = $scop_quality->[$k] ;
            }

#	       len  => $subset_length->[$k],
            $doms->{$sid} = {
               class => $subset_class->[$k],
               num_seg => $num_segments->[$k],
               scop_quality => $t_scopqual,
               len  => 0,
               domno => $cur_domno,
            } ;
         }
         $doms->{$sid}->{len} += ($end_resno->[$k] - $start_resno->[$k] + 1);

         $doms->{$sid}->{frag}->{$segment_no->[$k]} = {
            start => $start_resno->[$k],
            end => $end_resno->[$k],
            len => ($end_resno->[$k] - $start_resno->[$k] + 1 )
         } ;

         push @cells, {
            subset_id => $sid,
            seg_no => $segment_no->[$k],
            start => $start_resno->[$k],
            end => $end_resno->[$k]
         } ;
      }
      my $num_domains = $cur_domno ;

      my $nowon = 0;
      my $domarch ;
      my $unannot = 0 ;
      foreach my $k ( 0 .. $#cells ) {
         my $sid = $cells[$k]->{subset_id} ;
         my $fragno = $cells[$k]->{seg_no} ;
         if ($cells[$k]->{start} > ($nowon + 1)) {
            my $blank = $cells[$k]->{start} - $nowon - 1 ;
            $unannot += $blank ;
            $domarch->{s} .= "<$blank>" ;
            $domarch->{f} .= "<$blank>" ;
         } elsif ($cells[$k]->{start} < ($nowon + 1)) {
            my $overlap = $nowon + 1 - $cells[$k]->{start} ;
            $domarch->{s} .= "($overlap)" ;
            $domarch->{f} .= "($overlap)" ;
         }

         my $sep_s = '['; my $sep_e = ']' ;
         if ($cells[$k]->{seg_no} > 1) { $sep_s = '{'; $sep_e = '}' ; }

         $domarch->{ss} .= $sep_s.$doms->{$sid}->{class}.$sep_e ;
         $domarch->{s} .= $sep_s.$doms->{$sid}->{class}.$sep_e ;
         $domarch->{f} .= "[".join(':', (
            $doms->{$sid}->{domno},
            $doms->{$sid}->{class},
            $fragno,
            $doms->{$sid}->{frag}->{$fragno}->{len},
            $doms->{$sid}->{len},
            $doms->{$sid}->{scop_quality}, 
         ))."]" ;

         $nowon = $cells[$k]->{end} ;
      }

      if ($seq_length->[0] > $nowon) {
         my $blank = $seq_length->[0] - $nowon ;
         $unannot += $blank ;
         $domarch->{s} .= "<$blank>" ;
         $domarch->{f} .= "<$blank>" ;
      }

      my $cover_len = $seq_length->[0] - $unannot ;

      if ($display) {
         my @outvals = (
            $cur_seqid,
            $seq_length->[0],
            $cover_len,
            $num_domains,
            $domarch->{ss},
            $domarch->{s},
            $domarch->{f}
         ) ;
         print join("\t",@outvals)."\n" ;
      }

      if ($datareturn) {
         $results->{$cur_seqid} = {
            seq_length => $seq_length->[0],
            cover_length => $cover_len,
            num_domains => $num_domains,
            domarch_ss => $domarch->{ss},
            domarch_s => $domarch->{s},
            domarch_f => $domarch->{f}
         } ;
      }
   }

   return $results ;
}


=head2 read_seqid_domainarch()

   Title:       read_seqid_domainarch()
   Function:    Read in domain architecture assignments
   Args:        $_->{fn} = output filename
   Returns:     $->{seq_id} = DOMARCH

=cut

sub read_seqid_domainarch {

   my $in = shift ;
   my $fn = $in->{fn} ;
   if (!-s $fn) {
      die "ERROR: couldn't find seqid_domainarch file $fn\n"; }

   open(INF, $fn) ;
   my $domarch ;
   while (my $line = <INF>) {
      chomp $line;
      my ($t_seqid, $t_length, $t_cover_len, $t_num_domains, $t_domarch_ss,
          $t_domarch_s, $t_domarch) = split(/\t/, $line) ;
      $domarch->{$t_seqid} = $t_domarch ;
   }

   return $domarch ;
}



=head2 calc_model_pdbstart()

   Title:       calc_model_pdbstart()
   Function:    Greps PDB file for residue number of first ^ATOM record
                 for all models specifid in a seqid_domains file
   Args:        $_->{display} = option to display domains to STDOUT
                $_->{datareturn} = option to send back domain assignments as
                                   variable
                $_->{results} = reference to send back results

   Returns:     $_->{model_id} = {      #if $in->{datareturn} specified
                  pdbstart_resno
                }

                if in->{display}: prints above fields tab-delimited to STDOUT

=cut

sub calc_model_pdbstart {

#************************************************************************
#              WARNING - WARNING - WARNING - WARNING                    *
#                                                                       *
#                NO INSCODES!!! RESIDUE NUMBERS ONLY!!!                 *
#                    ONLY USE FOR MODPIPE MODELS                        *
#                                                                       *
#************************************************************************

   my $in = shift;
   my $results = $in->{results} ;
   my $bins = locate_binaries();

   my ($display, $datareturn) ;
   if (exists $in->{display}) { $display = $in->{display}; }
   if (exists $in->{datareturn}) { $datareturn = $in->{datareturn}; }

   foreach my $seq_id (sort keys %{$in->{models}}) {
    foreach my $model_id (sort keys %{$in->{models}->{$seq_id}}) {

      my $model_fn =  $in->{models}->{$seq_id}->{$model_id}->{basemodels_dir}.'/'.
                           substr($seq_id,0,3)."/$seq_id/models/$model_id.pdb.gz" ;

      my $pdbstart_resno = '' ;
      open(PDBF, $bins->{'zcat'}." $model_fn |") ;
      while (my $pdbline = <PDBF>) {
         if ($pdbline !~ /^ATOM/) {next;}
         $pdbstart_resno = substr($pdbline, 22, 4) ;
         last;
      }
      close(PDBF) ;

      if ($display) {
         my @outvals = (
            $model_id,
            $pdbstart_resno,
         ) ;
         print join("\t",@outvals)."\n" ;
      }

      if ($datareturn) {
         $results->{$model_id} = {
            pdbstart_resno => $pdbstart_resno,
         } ;
      }
    }
   }

   return $results ;

}


=head2 connect_modbase()

   Title:       connect_modbase()
   Function:    Connects to mysql database for MODBASE
   Args:        none
   Returns:     DBI database handle

=cut

sub connect_modbase {

   my $in = shift ; my $specs ;
   if (exists $in->{specs}) { $specs = $in->{specs} ; }
   else { $specs = set_modtie_specs() ; }

   my ($dbh,) = modtie::pibase::connect_pibase({
                  host => $specs->{modbase_specs}->{host},
                  db => $specs->{modbase_specs}->{db},
                  user => $specs->{modbase_specs}->{user},
                  pass => $specs->{modbase_specs}->{pass} }) ;

   return $dbh ;
}



# Scoring/Candidate listing module

#BUG HERE scoring: 1. now need to calculate template interface contacts on the fly, as the interatomic_contacts_tables are no longer stored in pibase -> NOT TRUE, intersubset_cotnacts are still gzipped on disk

sub readin_candilist {

   my $in = shift ;
   my $fn = $in->{fn} ;
   if (!-s $fn) {
      die "FATAL ERROR: candidate binary interactions list $fn not found\n"; }

   my $results ;
   open(FH, $fn) ;
   while (my $line = <FH>) {
      chomp $line;
      my ($sid1, $seq1, $model_id_1, $resrange1,
          $sid2, $seq2, $model_id_2, $resrange2) = split(/\t/, $line) ;
      my $candint= $seq1."\t".$model_id_1."\t".$resrange1."\t".$seq2."\t".$model_id_2."\t".$resrange2;
      my $tmpl = $sid1."\t".$sid2 ;
      $results->{$tmpl}->{candidate_interaction}->{$candint}++ ;
      $results->{$tmpl}->{domain_tmpltarg}->{$sid1}->{$seq1."\t".$model_id_1."\t".$resrange1}++ ;
      $results->{$tmpl}->{domain_tmpltarg}->{$sid2}->{$seq2."\t".$model_id_2."\t".$resrange2}++ ;
   }
   return $results ;

}


sub readin_tmpllist {

   my $in = shift ;
   my $fn = $in->{fn} ;

   if (! -s $fn) {
      die "FATAL ERROR: template interface file $fn not found\n" ; }

   my $results ;
   my $seq ;

   open(FH, $fn) ;
   while (my $line = <FH>) {
      chomp $line;
      my ($bdp_id, $sid1, $sid2) = split(/\t/, $line) ;
      push @{$results->{bdp_id}},  $bdp_id;
      push @{$results->{sid1}},  $sid1 ;
      push @{$results->{sid2}},  $sid2 ;
   }

   return $results;
}


sub readin_badalis_alipaths_err {

   my $in = shift ;
   my $fn = $in->{fn} ;
   if (! -s $fn) {
      die "FATAL ERROR: badali paths file $fn not found\n" ; }

   my $com ;
   my $bins = locate_binaries();
   if ($fn =~ /gz$/) { $com = $bins->{'zcat'}." $fn" ;
   } else { $com = "cat $fn" ; }

   open(FH, "$com |") ;

   my $badalis ;
   while (my $line = <FH>) {
      chomp $line;
      if ($line !~ /^ERROR/ && $line !~ /WARNING/) {next;}
      $line =~ s/^ERROR //g ; $line =~ s/^WARNING //g ;
      $line =~ s/\:.*$//g ;
      $line =~ s/\(// ; $line =~ s/\)// ;
      $line =~ s/,//g ;
      my ($sid, $model_id, $res_range) = split(' ', $line) ;
      $badalis->{$sid."\t".$model_id."\t".$res_range}++ ;
   }

   return $badalis ;
}


=head2 scoring_pibase_preload()

   Title:       scoring_pibase_preload()
   Function:    Preload PIBASE data necessary for scoring routine
   Args:        NONE
   Returns:     $->{meta_tables}->{TABLE_TYPE}->{BDP_ID} = FILE_PATH
                $->{tn}->{subsets_res} = $->{meta_tables}->{subsets_residues}
                $->{sid2scop}->{SUBSET_ID} = CLASS
                $->{sid2bdp}->{SUBSET_ID} = BDP_ID
                $->{sid2length}->{SUBSET_ID} = sequence_length
                $->{sid2aaseq}->{SUBSET_ID} = aa_seq
                $->{all_subsdef}->{SUBSET_ID}->{ch} = chain
                $->{all_subsdef}->{SUBSET_ID}->{s} = start residue
                $->{all_subsdef}->{SUBSET_ID}->{e} = end residue

=cut

sub scoring_pibase_preload {

   my $tn ;
   my $meta_tables;

   my $specs = set_modtie_specs() ;

   print STDERR "PIBASE data load\n" ;
   print STDERR "\tsubsets_residues_tables: " ;
   {
      my @tod_res = modtie::pibase::rawselect_tod(
         "SELECT bdp_id, source_file FROM subsets_residues_tables") ;
      foreach my $j ( 0 .. $#{$tod_res[0]}) {
         $meta_tables->{subsets_residues}->{$tod_res[0]->[$j]} =
            $tod_res[1]->[$j];
      }
      replace_basedir_inplace({
         hash => $meta_tables->{subsets_residues},
         old_dir => $specs->{pibase_specs}->{old_metatod_dir}.
            "/subsets_residues",
         new_dir => $specs->{pibase_specs}->{metatod_dir}."/subsets_residues"
      }) ;
   }
   $tn->{subsets_res} = $meta_tables->{subsets_residues};
   print STDERR "X\n" ;


   print STDERR "\tsubsets_files: " ;
   {
      my @tod_res = modtie::pibase::rawselect_tod(
         "SELECT subset_id, file_path FROM subsets_files") ;
      foreach my $j ( 0 .. $#{$tod_res[0]}) {
         $meta_tables->{subsets_files}->{$tod_res[0]->[$j]} =
            $tod_res[1]->[$j];
      }
      replace_basedir_inplace({
         hash => $meta_tables->{subsets_files},
         old_dir => $specs->{pibase_specs}->{old_subsets_dir},
         new_dir => $specs->{pibase_specs}->{subsets_dir}
      }) ;
   }
   print STDERR "X\n" ;


   my ($sid2scop, $sid2bdp) ;
   print STDERR "\tsubsets: " ;
   {
      my @tod_res = modtie::pibase::rawselect_tod(
         "SELECT bdp_id, subset_id, class FROM subsets") ;
      foreach my $j ( 0 .. $#{$tod_res[0]}) {
         if ($tod_res[0]->[$j] ne "") {
            $sid2bdp->{$tod_res[1]->[$j]} = $tod_res[0]->[$j] ; }
#fpd100921_1437  - Used to check for subset_source_id=1, but new
# pibase2010 is different; doesn't really matter
         $sid2scop->{$tod_res[1]->[$j]} = $tod_res[2]->[$j];
      }
   }
   print STDERR "X\n" ;

   my $sid2length ;
   print STDERR "\tsubsets_sequence: " ;
   {
      my @tod_res = modtie::pibase::rawselect_tod(
         "SELECT subset_id, num_res, sequence FROM subsets_sequence") ;
      foreach my $j ( 0 .. $#{$tod_res[0]}) {
         $sid2length->{$tod_res[0]->[$j]} = $tod_res[1]->[$j];
      }
   }
   print STDERR "X\n" ;


   my $pb = {
      meta_tables => $meta_tables,
      tn => $tn,
      sid2scop => $sid2scop,
      sid2bdp => $sid2bdp,
      sid2length => $sid2length,
   } ;

   return $pb ;
}


sub setoutfiles_scoring_main {
   my $in = shift ;
   my $fh = $in->{fh} ;
   my $fn = $in->{fn} ;
   my $params = $in->{params} ;

   if ($params->{alilist_fl}) {
      if (!exists $params->{out_fn}->{alilist_fn}) {
         ($fh->{ali},$fn->{ali})=tempfile("alignments.XXXXX", SUFFIX=>".out");
         print STDERR "\talignments: $fn->{ali}\n" ;
      } else {
         $fn->{ali} = $params->{out_fn}->{alilist_fn} ;
         open($fh->{ali}, ">$fn->{ali}") ;
      }

      if (!exists $params->{out_fn}->{cutlist_fn}) {
         ($fh->{cutlist},$fn->{cutlist})=tempfile("cutlist.XXXXX", SUFFIX=>".out");
         print STDERR "\tcutlist: $fn->{cutlist}\n" ;
      } else {
         $fn->{cutlist} = $params->{out_fn}->{cutlist_fn} ;
         open($fh->{cutlist}, ">$fn->{cutlist}") ;
      }
   }

   if ($params->{postali_candilist_fl}) {
      if (!exists $params->{out_fn}->{postali_candilist_fn}) {
         ($fh->{postali_cand},$fn->{postali_cand}) =
            tempfile("postali_candidates.XXXXX", SUFFIX=>".out");
         print STDERR "\tpostali candidates: $fn->{postali_cand}\n" ;
      } else {
         $fn->{postali_cand} = $params->{out_fn}->{postali_candilist_fn} ;
         open($fh->{postali_cand}, ">$fn->{postali_cand}") ;
      }
   }

#($fh->{binscores},$fn->{binscores})=tempfile("binscores.XXXXX", SUFFIX=>".out")
   if ($params->{assess_fl} || $params->{benchmark_fl}) {
      if (!exists $params->{out_fn}->{interactions_fn}) {
         $params->{out_fn}->{interactions_fn} = '-'; }

      $fn->{binscores} = $params->{out_fn}->{interactions_fn} ;
      open($fh->{binscores}, ">$fn->{binscores}") ;
   }

}


=head2 scoring_main()

   Title:       scoring_main()
   Function:    Main scoring routine
   Args:        $_->{scoplevel}
                $_->{alilist_fl}           #run mode: list alignments
                $_->{postali_candilist_fl} #run mode: make candidate list
                $_->{assess_fl}            #run mode: assess candidates

                $_->{in_fn}->{seqid_domainarch_fn}
                $_->{in_fn}->{seqdomains_fn}

   Returns:     NOTHING
   Displays:    to fh->{binscores_fn} (as set by setoutfiles_scoring_main()):
                Lots of stuff including sequence identifiers, model identifiers,
                  domain identifiers, raw score, z-score, alignment info,
                  sequence identities, etc...

=cut

sub scoring_main {

   my $params = shift ;
   my $specs ;
   if (exists $params->{specs}) { $specs = $params->{specs} ; }
   else { $specs = set_modtie_specs() ; }

   my $headers;
   @{$headers->{assess}} = qw/INTERFACE subset_id_1 subset_id_2 seq1 model_id_1 t_def1 class_1 seq2 model_id_2 t_def2 class_2 tmpl_contacts aln_contacts pot_file pot_source pot_type pot_dist raw z z_p1 z_p2 avg min min_tn stdev z_min z_min_tn false_pos rmsd_dom1 equivpos_dom1 numres_tmpl_dom1 numres_targ_dom1 numres_ident_dom1 numres_tmpl_bs_dom1 numres_targ_bs_dom1 numres_ident_bs_dom1 rmsd_dom2 equivpos_dom2 numres_tmpl_dom2 numres_targ_dom2 numres_ident_dom2 numres_tmpl_bs_dom2 numres_targ_bs_dom2 numres_ident_bs_dom2/ ;

   if (exists $params->{get_assess_headers}) {
      return $headers->{assess} ; }

# 0. figure out runmode
   if (!exists $params->{scoplevel}) {
      $params->{scoplevel} = $specs->{scoplevel}->{default} ; }

   if (!exists $params->{alilist_fl} &&
       !exists $params->{postali_candilist_fl} &&
       !exists $params->{assess_fl}) {
      $params->{$specs->{runmode}->{default}}++ ;
   }

   if (!exists $params->{in_fn}->{seqid_domainarch_fn}) {
      croak "ERROR: sequence id / arch not specified\n" ; }

   if (!exists $params->{in_fn}->{model_pdbstart_fn}) {
      croak "ERROR: model_pdbstart file not specified\n" ; }

#1a. read in target protein sequences and domain architectures
   my $seqdomfraginfo = readin_seqdomains_start_modelid({
      fn => $params->{in_fn}->{seqdomains_fn}}) ;
   my $modelid_2_pdbstart = readin_model_pdbstart({
      fn => $params->{in_fn}->{model_pdbstart_fn}}) ;

   my $seqfrag2model = $seqdomfraginfo->{seqfrag2model} ;
   my $seqid_2_run = $seqdomfraginfo->{seqid_2_run} ;
   my $modelid_2_target_beg = $seqdomfraginfo->{modelid_2_target_beg} ;
   my $basemodels_dir = $seqdomfraginfo->{basemodels_dir} ;
   my $baseali_dir = $seqdomfraginfo->{baseali_dir} ;
   my $modpipe_style = $seqdomfraginfo->{modpipe_style} ;

   $params->{seqinfo}->{arch} = read_seqid_domainarch({
      fn => $params->{in_fn}->{seqid_domainarch_fn}});
   map {$params->{seqs}->{$_}++;} (keys %{$params->{seqinfo}->{arch}});

   my $fh = {}; my $fn = {};
   setoutfiles_scoring_main({params => $params, fh => $fh, fn => $fn}) ;
   my $pb_pl = scoring_pibase_preload() ;

   if ($params->{alilist_fl} || $params->{postali_candilist_fl}) {
      my $partial_badalis = {} ;
      if (defined $params->{in_fn}->{partial_alipaths_err_fn} &&
          -s $params->{in_fn}->{partial_alipaths_err_fn}) {
         $partial_badalis = readin_badalis_alipaths_err(
            {fn => $params->{in_fn}->{partial_alipaths_err_fn}}) ;
      }

      $params->{seqset} = read_seqid_setinfo({
         fn => $params->{in_fn}->{seqid_setinfo_fn}}) ;

      if (!exists $params->{in_fn}->{templates_fn}) {
         $params->{in_fn}->{templates_fn} = $specs->{templates_fn}->{default} ;}

      my $pbtmpl = readin_tmpllist({fn => $params->{in_fn}->{templates_fn}});

      my $mbinfo = _parse_domarchs($params) ;
      my ( $scop2seq, $scop2seqset, $seq2dom, $seq2domfrag, $scop2num ) =
         ($mbinfo->{scop2seq},
          $mbinfo->{scop2seqset},
          $mbinfo->{seq2dom},
          $mbinfo->{seq2domfrag},
          $mbinfo->{scop2num}) ;


      if (!exists $params->{seqset}) {
         @{$params->{seqset}->{1}} = keys %{$params->{seqs}} ;
         foreach my $tseq (@{$params->{seqset}->{1}}) {
            $params->{seqset}->{union}->{$tseq}++ ;
         }
      }

      my $tmplcounts ;
      foreach my $tmplno ( 0 .. $#{$pbtmpl->{sid1}} ) {
         my ($bdp_id, $sid) ;
         ($bdp_id, $sid->[0], $sid->[1]) =
            ($pbtmpl->{bdp_id}->[$tmplno],
             $pbtmpl->{sid1}->[$tmplno],
             $pbtmpl->{sid2}->[$tmplno]) ;

         if (!exists $pb_pl->{tn}->{subsets_res}->{$bdp_id}) {
            print STDERR "ERROR: couldnt find subsets_res for $bdp_id\n" ;
            next;
         }


         my $scop_fa ;
         $scop_fa->[0] = $pb_pl->{sid2scop}->{$sid->[0]} ;
         $scop_fa->[1] = $pb_pl->{sid2scop}->{$sid->[1]} ;
#         print STDERR "now on interface template $tmplno ($sid->[0] $scop_fa->[0] -- $sid->[1] $scop_fa->[1]) ".__LINE__."\n";

         my $dom_types ;

         if ($params->{scoplevel} eq 'superfamily') {
            ($dom_types->[0]) = ($scop_fa->[0] =~ /(.+)\.[0-9]+$/) ;
            ($dom_types->[1]) = ($scop_fa->[1] =~ /(.+)\.[0-9]+$/) ;
         } elsif ($params->{scoplevel} eq 'family') {
            $dom_types->[0] = $scop_fa->[0] ;
            $dom_types->[1] = $scop_fa->[1] ;
         }

         if (!exists $scop2seq->{$dom_types->[0]} || 
             !exists $scop2seq->{$dom_types->[1]}) { next;}

         my $crosstype ;
         if ($params->{xset_fl} == 1) {
            if (exists $scop2seqset->{1}->{$dom_types->[0]} &&
                exists $scop2seqset->{2}->{$dom_types->[1]}) {
                  $crosstype->{12}++;}

            if (exists $scop2seqset->{2}->{$dom_types->[0]} &&
                exists $scop2seqset->{1}->{$dom_types->[1]}) {
                  $crosstype->{21}++;}

            if (!exists $crosstype->{12} && !exists $crosstype->{21}) {
               next; }
         } elsif ($params->{xset_fl} == 0) {
            $crosstype->{1}++ ;
         }

         my @t = sort ($dom_types->[0], $dom_types->[1]) ;
         $tmplcounts->{inttype}->{$t[0]."\n".$t[1]}->{$sid->[0]."\n".$sid->[1]}++ ;
         $tmplcounts->{inttype2tmpltype}->{$t[0]."\n".$t[1]}->{$sid->[0]."\n".$sid->[1]}->{$sid->[0]} = $dom_types->[0] ;
         $tmplcounts->{inttype2tmpltype}->{$t[0]."\n".$t[1]}->{$sid->[0]."\n".$sid->[1]}->{$sid->[1]} = $dom_types->[1] ;
         $tmplcounts->{inttmpl}++ ;

         $tmplcounts->{x}->{$dom_types->[0]}->{$sid->[0]}++ ;
         $tmplcounts->{x}->{$dom_types->[1]}->{$sid->[1]}++ ;
      }
      my $numtypes = keys %{$tmplcounts->{inttype}} ;
      print STDERR "NOTE: $tmplcounts->{inttmpl} interface templates ".
         "($numtypes types) could be used (no alignment considerations yet)\n";


#YO 061206_1134: NEED PROPER XSET CHECKING
      if ($params->{alilist_fl}) {
         print STDERR "NOW: making list of alignments necessary: " ;
         my $numali = 0 ;
         foreach my $domtype (keys %{$tmplcounts->{x}}) {
            my $numx = keys %{$tmplcounts->{x}->{$domtype}} ;
            if (!exists $scop2seq->{$domtype}) {next;}
            my $numt = $scop2num->{$domtype} ;
            $numali +=  $numx * $numt ;
# note number alignments needed is going to be less than specified in numali
# becuase of the size mismatches; numali is an upper bound.

            my $cuthash ;
            foreach my $x ( keys %{$tmplcounts->{x}->{$domtype}}) {
               foreach my $seq ( keys %{$scop2seq->{$domtype}}) {
                  foreach my $dom (@{$scop2seq->{$domtype}->{$seq}}) {

                     my $relsize = $pb_pl->{sid2length}->{$x} /
                                    $seq2dom->{$seq}->[$dom]->{domsize} ;
                     if ($relsize > 1) {$relsize = 1 / $relsize;}
                     if ($relsize < $specs->{targ_tmpl_relative_size}) { next; }

                     my @defs ;
                     my @cutdefs ;
                     foreach my $j ( 1 .. $#{$seq2domfrag->{$seq}->[$dom]}) {
                        push @defs, $seq2domfrag->{$seq}->[$dom]->[$j]->{start}.
                                 '-'.$seq2domfrag->{$seq}->[$dom]->[$j]->{end} ;
                        push @cutdefs,
                           ('', $seq2domfrag->{$seq}->[$dom]->[$j]->{start},
                                  $seq2domfrag->{$seq}->[$dom]->[$j]->{end}) ;
                     }
                     my $def = join('_', @defs) ;


                     my $model_id =
         $seqfrag2model->{$seq}->{$seq2domfrag->{$seq}->[$dom]->[1]->{start}};

                     if (exists $partial_badalis->{$x."\t".$model_id."\t".$def})
                        { next; }

                     my $t_alioutdir = sid_modelid_2_alidir({
                                          specs => $specs,
                                          sid => $x,
                                          modelid => $model_id});

                     my $t_alipath =
                        "$t_alioutdir/align.$x.$model_id.$def.salign.gz";
                     if (-s $t_alipath) {next;}

                     print {$fh->{ali}} $x."\t".$model_id."\t".$def."\n" ;

                     if (!exists $cuthash->{$seq.$def}) {
                        my $t_fn ;

                        my $t_cutdir = modelid_2_domdir({
                                          modelid => $model_id,
                                          specs => $specs }) ;
                        my $t_cutpdb = $t_cutdir.'/'.$model_id.'.'.$def.'.pdb' ;
                        if (-s $t_cutpdb) {
                           $cuthash->{$seq.$def}++ ;
                           next;
                        }

                        $t_fn = $basemodels_dir->{$model_id}.'/'.
                           substr($seq,0,3)."/$seq/models/$model_id.pdb.gz" ;

# fpd100925_2104 - changing to adjust to new MODPIPE residue numbering
# old logic: if modpipe_style='new': provide target_beg number; eg ModPipe2.0
#                         otherwise: PDB number = real numbering; eg ModPipe1.0
#
# new logic: if $modelid_2_pdbstart != $modelid_2_target_beg: provide number
#                                                  otherwise: PDB number = real
#  *modpipe_style was a proxy for this behavior, but no longer holds
#   eg new modpipe runs, eg 2.0.2 are 'new' style, but numbering is old

                        if ($modelid_2_pdbstart->{$model_id} !=
                            $modelid_2_target_beg->{$model_id}) {
                           print {$fh->{cutlist}} join("\t",
                              $t_fn,
                              'modelnew'.$modelid_2_target_beg->{$model_id},
                              $model_id,
                              @cutdefs)."\n";
                        } else {
                           print {$fh->{cutlist}} join("\t",
                              $t_fn,
                              'modelold',
                              $model_id,
                              @cutdefs)."\n";
                        }

                        $cuthash->{$seq.$def}++ ;
                     }
                  }
               }
            }
         }
         print STDERR "$numali alignments needed\n" ;
      }


      if ($params->{postali_candilist_fl}) {
         my $badalis ;
         if (exists $params->{in_fn}->{alipaths_err_fn}) {
            $badalis = readin_badalis_alipaths_err({
               fn => $params->{in_fn}->{alipaths_err_fn}}) ; }

         my $numcand = 0 ;
         foreach my $inttype (keys %{$tmplcounts->{inttype}}) {
            my ($d1, $d2) = split(/\n/, $inttype) ;

            my $candset ;
            my $candset2len ;
            my @domtypes = ($d1);
            if ($d1 ne $d2) {push @domtypes, $d2;}

            foreach my $domtype (@domtypes) {
               foreach my $seq ( keys %{$scop2seq->{$domtype}}) {
                  foreach my $dom (@{$scop2seq->{$domtype}->{$seq}}) {
                     my @defs ;
                     foreach my $j ( 1 .. $#{$seq2domfrag->{$seq}->[$dom]}) {
                        push @defs, $seq2domfrag->{$seq}->[$dom]->[$j]->{start}.
                                    '-'.$seq2domfrag->{$seq}->[$dom]->[$j]->{end} ;
                     }
                     my $def = join('_', @defs) ;
                     my $model_id =
          $seqfrag2model->{$seq}->{$seq2domfrag->{$seq}->[$dom]->[1]->{start}};

                     $candset->{$domtype}->{$seq."\t".$model_id."\t".$def}++ ;
                     $candset2len->{$seq."\t".$model_id."\t".$def} =
                        $seq2dom->{$seq}->[$dom]->{domsize};
                  }
               }
            }

            my $seen ;
            foreach my $cur_tmpl (sort keys
               %{$tmplcounts->{inttype}->{$inttype}}) {

               my ($x1, $x2) = split("\n", $cur_tmpl) ;

               my $sids = [$x1, $x2] ;
               my $domtypes = [
         $tmplcounts->{inttype2tmpltype}->{$inttype}->{$cur_tmpl}->{$x1},
         $tmplcounts->{inttype2tmpltype}->{$inttype}->{$cur_tmpl}->{$x2} ] ;

               print STDERR "STATUS: template $x1 -- $x2 ".
                  "($domtypes->[0] -- $domtypes->[1]), " ;

               my $sizecheck_cands ;
               foreach my $side (0 .. 1) {
                  foreach my $cand (sort keys %{$candset->{$domtypes->[$side]}}) {
                     my ($seq, $modelid, $def) = split(/\t/, $cand) ;

                     my $relsize = $pb_pl->{sid2length}->{$sids->[$side]} /
                        $candset2len->{$cand} ;

                     if ($relsize > 1) {$relsize = 1 / $relsize;}
                     if ($relsize < $specs->{targ_tmpl_relative_size}) { next; }

                     if (exists $badalis->{$sids->[$side].
                                 "\t".$modelid."\t".$def}) {
                        next; }

                     push @{$sizecheck_cands->[$side]}, $cand ;
                  }
               }

               if ($#{$sizecheck_cands->[0]} < 0 ||
                   $#{$sizecheck_cands->[1]} < 0) {
                  print STDERR " candidates: none\n" ;
                  next;
               }

               my $numside1 = $#{$sizecheck_cands->[0]} + 1;
               my $numside2 = $#{$sizecheck_cands->[1]} + 1;
               print STDERR " candidates: $numside1 vs $numside2\n" ;

               foreach my $cand1 (@{$sizecheck_cands->[0]}) {

                  my ($seq1, $cur_modelid_1, $cur_def_1) = split(/\t/, $cand1) ;

                  my $set1 ;
                  if ($params->{xset_fl} == 1 && defined $scop2seqset) {
                     $set1 = $params->{seqset}->{union}->{$seq1} ; }

                  foreach my $cand2 (@{$sizecheck_cands->[1]}) {

                     my ($seq2, $cur_modelid_2, $cur_def_2) =
                        split(/\t/, $cand2) ;

                     my $set2 ;
                     if ($params->{xset_fl} == 1 && defined $scop2seqset) {
                        $set2 = $params->{seqset}->{union}->{$seq2} ;
                        if ($set1 != 3 && $set1 == $set2) {next;}
                     }

                     print {$fh->{postali_cand}} join("\t", ($x1, $cand1,
                           $x2, $cand2))."\n" ;
                     $numcand++ ;
                  }
               }
            }
         }
         print STDERR "$numcand binary interaction candidates (target/tmpl interaction entries) \n" ;
      }

      if ($params->{alilist_fl}) {close $fh->{ali} ;}
      if ($params->{postali_candilist_fl}) {close $fh->{postali_cand} ;}
   }

   if ($params->{assess_fl}) {

#100326_1131      my $aaseq = $params->{seqinfo}->{aaseq} ;

      my ($pots, $const_rmax, @const_trs) ;
      if (!exists $params->{pots}) {
         $pots = $specs->{pots} ;
      } else {
         $pots = $params->{pots} ;
      }
      readin_potentials($pots) ;
      foreach my $pot (@{$pots}) { push @const_trs, $pot->{r} ; }
      ($const_rmax, undef) = sort {$b <=> $a} @const_trs ;


      my $candidates = readin_candilist(
         {fn => $params->{in_fn}->{postali_candilist_fn}}) ;

      foreach my $tmpl ( sort keys %{$candidates}) {
         my $sid ;
         ($sid->[0], $sid->[1]) = split(/\t/, $tmpl) ;
         my $scop_fa;
         my $bdp_id = $pb_pl->{sid2bdp}->{$sid->[0]} ;
         print STDERR "now on template $sid->[0] -- $sid->[1] bdp $bdp_id\n" ;

         my ($t_pdbfh, $t_pdbfn) = tempfile("temppdb.XXXXX", SUFFIX => ".pdb",
                                             CLEANUP => 0); close($t_pdbfh) ;

         foreach my $sid_pdbfn (
          $pb_pl->{meta_tables}->{subsets_files}->{$sid->[0]},
          $pb_pl->{meta_tables}->{subsets_files}->{$sid->[1]}) {
            my $catcom = 'cat' ;
            if ($sid_pdbfn =~ /\.gz$/) {$catcom = $specs->{binaries}->{zcat};}
            my $tcom = "$catcom $sid_pdbfn >> $t_pdbfn" ;
            system($tcom) ;
         }

         my $strx_path = $t_pdbfn ;

         my $tmpl_subsres = get_subsres($pb_pl->{tn}->{subsets_res}->{$bdp_id});

         my $tmpl_con = get_domain_contacts_complex({
            subset_residues => $tmpl_subsres,
            pdb_fn => $t_pdbfn,
            R_max => $const_rmax,
            R_list => \@const_trs
         }) ;
         unlink $t_pdbfn ;

         my $aln_cache ;
         my $targetres_cache ;
         my $def_cache ;
         foreach my $side (0 .. 1) {
            $scop_fa->[$side] = $pb_pl->{sid2scop}->{$sid->[$side]} ;

            foreach my $target (keys
                  %{$candidates->{$tmpl}->{domain_tmpltarg}->{$sid->[$side]}}) {
               my ($seq, $model_id, $def) = split(/\t/, $target) ;
               my $model_id_def = $model_id."\t".$def ;

               my $target_beg = $modelid_2_target_beg->{$model_id} ;

               my $cur_model_path = $basemodels_dir->{$model_id} ;
               if (!defined $cur_model_path) {
                  if (exists $specs->{local_modbase_models_dir}) {
                     $cur_model_path = $specs->{local_modbase_models_dir} ;
                  } else {
                     print STDERR "WARNING: model_path not found for ".
                           "model $model_id, seq $seq; skipping this model\n";
                     next;
                  }
               }

               {
                  my $ali_fn = sid_modelid_2_alidir({
                     specs => $specs,
                     modelid => $model_id,
                     sid => $sid->[$side]
                  })."/align.".$sid->[$side].".$model_id.$def.salign.gz";

                  if (!-s $ali_fn) {
                     print STDERR "ERROR: $ali_fn not found ".
                           "(seq $seq model $model_id range $def ) ".
                           "aligned to $sid->[$side])\n" ;
                     next;
                  }
#100329_2341 brought back in def_cache
                  $def_cache->{$seq}->{$model_id_def} = $def ;

# fpd100925_2111 - changing to adjust to new MODPIPE residue numbering
# old logic: if modpipe_style='new': provide target_beg number; eg ModPipe2.0
#                         otherwise: PDB number = real numbering; eg ModPipe1.0
#
# new logic: if $modelid_2_pdbstart != $modelid_2_target_beg: provide number
#                                                  otherwise: PDB number = real
#  *modpipe_style was a proxy for this behavior, but no longer holds
#   eg new modpipe runs, eg 2.0.2 are 'new' style, but numbering is old

                  my $t_resmap ;
                  if ($modelid_2_pdbstart->{$model_id} !=
                      $modelid_2_target_beg->{$model_id}) {
                     $t_resmap = salign_resequiv_parse({
                        target_beg => $target_beg,
                        fn => $ali_fn
                     }) ;
                  } else {
                     $t_resmap = salign_resequiv_parse({ fn => $ali_fn }) ;
                  }

                  if (!defined $t_resmap->{alistats}->{quality_score}){
                     print STDERR "ERROR: no quality score computed, ".
                        "bad alignment ".
                        "seq $seq model $model_id range $def ".
                        "aligned to $sid->[$side])\n" ;
                     next;
                  }


                  if ((exists $t_resmap->{alistats}->{quality_score}) &&
                      ($t_resmap->{alistats}->{quality_score} < 
                          $specs->{salign_quality_score_thresh})) {
                        print STDERR "WARNING: skipping low quality alignment ".
                           " (quality_score = ".
                           $t_resmap->{alistats}->{quality_score}."):".
                           "seq $seq model $model_id range $def ".
                           "aligned to $sid->[$side])\n" ;
                        next;
                  }

# also parse in the RMSD and sequence identity values here?
                  if (!defined $t_resmap->{resequiv}) {next;}
                  $aln_cache->{$side}->{resmap}->{$seq}->{$model_id_def} =
                     $t_resmap->{resequiv} ;
                  $aln_cache->{$side}->{resnames}->{$seq}->{$model_id_def} =
                     $t_resmap->{resnames}->[1] ;
                  $aln_cache->{$side}->{alistats}->{$seq}->{$model_id_def} =
                     $t_resmap->{alistats} ;
               }
            }
         }

         my $binary_counter = 1;
         foreach my $candint (keys
               %{$candidates->{$tmpl}->{candidate_interaction}}) {
            print STDERR "now on candidate $binary_counter " ;
            my ($seq1, $model_id_1, $resrange1, $seq2, $model_id_2, $resrange2)=
               split(/\t/, $candint) ;
            my $dom1 = $model_id_1."\t".$resrange1 ;
            my $dom2 = $model_id_2."\t".$resrange2 ;

#100326_1131            my $targetres_1 = $targetres_cache->{$seq1} ;
#100326_1131            my $targetres_2 = $targetres_cache->{$seq2} ;

            my ($resmap, $resnames) ;
            $resmap->{$sid->[0]} = $aln_cache->{0}->{resmap}->{$seq1}->{$dom1} ;
            $resnames->{$sid->[0]} =
               $aln_cache->{0}->{resnames}->{$seq1}->{$dom1};

            $resmap->{$sid->[1]}= $aln_cache->{1}->{resmap}->{$seq2}->{$dom2} ;
            $resnames->{$sid->[1]} =
               $aln_cache->{1}->{resnames}->{$seq2}->{$dom2};

            my $model_con = contact_equiv_complex({
               contacts => $tmpl_con->{respairs},
               tmpl_resnames => $tmpl_con->{resnames},
               seq => $resnames,
               resmap => $resmap,
               tmpl_subsres => $tmpl_subsres,
               dompair2respair => $tmpl_con->{dompair2respair},
               labelnames => {
                  targ1 => $model_id_1,
                  tmpl1 => $sid->[0],
                  targ2 => $model_id_2,
                  tmpl2 => $sid->[1],
               }
            });

            my $num_tmplcon =
                     keys %{$tmpl_con->{respairs}->{$sid->[0]."\n".$sid->[1]}} ;
            my $num_alncon =
                     keys %{$model_con->{respairs}->{$sid->[0]."\n".$sid->[1]}};

            if ($num_alncon == 0) {
               print STDERR " ERROR: no aligned contacts in candidate ".
                  "interaction ($sid->[0] -- $sid->[1]  aligned to ".
                  "$model_id_1 -- $model_id_2 )\n";
               next;
            }

            if (($num_alncon / $num_tmplcon) <
                 $specs->{complex_aln_tmpl_contacts_thresh}) {
               print STDERR " WARNING: very low number of aligned contacts ".
                  "in candidate interaction ($sid->[0] -- $sid->[1]  ".
                  "aligned to $model_id_1 -- $model_id_2 ), skipping\n";
               next;
            }

            foreach my $pot ( @{$pots}) {
               my $rawscore = score_potential_complex({
                  contacts => $model_con->{respairs},
                  pot => $pot,
                  resnames => $model_con->{resnames},
                  dompair2respair => $model_con->{dompair2respair}
               }) ;

               my $ranscores = get_rand_score_complex({
                  contacts => $model_con->{respairs},
                  pot => $pot,
                  resnames => $model_con->{resnames},
                  dompair2respair => $model_con->{dompair2respair}
               }) ;

               my $zscore = {};
# FOR NOW COMMENTED OUT BECAUSE REDUNDANT - but code is correct
#                     $zscore->{whole}= calc_zscore($rawscore->{whole}, $ranscores->{whole}) ;
#                     if (exists $zscore->{whole}->{error}) {
#                        my @outvals = ('COMPLEX', $sid->[0], $sid->[1],
#    $seq1, $model_id_1, $def_cache->{$seq1}->{$dom1}, $scop_fa->[0],
#    $seq2, $model_id_2, $def_cache->{$seq2}->{$dom2}, $scop_fa->[1],
#    $pot->{fn}, $pot->{type}, $pot->{details}, $pot->{r},
#    $rawscore->{whole}, 'ERROR', $zscore->{whole}->{error});
#
#                        print {$fh->{binscores}} join("\t", @outvals)."\n" ;
#                        print STDERR "ERROR ($sid->[0] -- $sid->[1]): ".
#                           $zscore->{whole}->{error}."\n";
#                        next;
#                     }
#
#                     my @outvals = ('COMPLEX',
#    $seq1, $model_id_1, $def_cache->{$seq1}->{$dom1}, $sid->[0], $scop_fa->[0],
#    $seq2, $model_id_2, $def_cache->{$seq2}->{$dom2}, $sid->[1], $scop_fa->[1],
#    $num_tmplcon, $num_alncon,
#    $pot->{fn_base}, $pot->{type}, $pot->{details}, $pot->{r},
#    $rawscore->{whole}, $zscore->{whole}->{z_score},
#    $zscore->{whole}->{z_prime}, $zscore->{whole}->{z_2},
#    $zscore->{whole}->{avg}, $zscore->{whole}->{min},
#    $zscore->{whole}->{min_tn}, $zscore->{whole}->{stdev},
#    $zscore->{whole}->{z_min}, $zscore->{whole}->{z_min_tn},
#    $zscore->{whole}->{false_pos}) ;
#
#                     print {$fh->{binscores}} join("\t", @outvals)."\n" ;
#

# in general display more alignment information - eg num aligned res, num ident res

               my $dom12sig = $sid->[0]."\n".$sid->[1] ;
               $zscore->{interface}->{$dom12sig} =
                      calc_zscore($rawscore->{interface}->{$dom12sig},
                                  $ranscores->{interface}->{$dom12sig});

               my $outdom12 = $dom12sig; $outdom12 =~ s/\n/\t/g ;
               if (exists $zscore->{whole}->{error}) {

                  my @outvals = ('INTERFACE', $outdom12,
      $seq1, $model_id_1, $def_cache->{$seq1}->{$dom1}, $scop_fa->[0],
      $seq2, $model_id_2, $def_cache->{$seq2}->{$dom2}, $scop_fa->[1],
      $num_tmplcon, $num_alncon,
      $pot->{fn}, $pot->{type}, $pot->{details}, $pot->{r},
      $rawscore->{interface}->{$dom12sig}, 'ERROR',
      $zscore->{interface}->{$dom12sig}->{error}) ;

                  print {$fh->{binscores}} join("\t", @outvals)."\n" ;
                  print STDERR " ERROR: ".
                  $zscore->{interface}->{$dom12sig}->{error}."\n";
                     next;

               }


               if (($rawscore->{interface}->{$dom12sig} >
                    $specs->{printthresh_score}) &&
                   ($zscore->{interface}->{$dom12sig}->{z_score} > 
                    $specs->{printthresh_zscore})) {next;}

               my @outvals = ('INTERFACE',
                   $outdom12,
                   $seq1,
                   $model_id_1,
                   $def_cache->{$seq1}->{$dom1},
                   $scop_fa->[0],
                   $seq2,
                   $model_id_2,
                   $def_cache->{$seq2}->{$dom2},
                   $scop_fa->[1],
                   $num_tmplcon,
                   $num_alncon,
                   $pot->{fn_base},
                   $pot->{type},
                   $pot->{details},
                   $pot->{r},
                   $rawscore->{interface}->{$dom12sig},
                   $zscore->{interface}->{$dom12sig}->{z_score},
                   $zscore->{interface}->{$dom12sig}->{z_prime},
                   $zscore->{interface}->{$dom12sig}->{z_2},
                   $zscore->{interface}->{$dom12sig}->{avg},
                   $zscore->{interface}->{$dom12sig}->{min},
                   $zscore->{interface}->{$dom12sig}->{min_tn},
                   $zscore->{interface}->{$dom12sig}->{stdev},
                   $zscore->{interface}->{$dom12sig}->{z_min},
                   $zscore->{interface}->{$dom12sig}->{z_min_tn},
                   $zscore->{interface}->{$dom12sig}->{false_pos},
                   $aln_cache->{0}->{alistats}->{$seq1}->{$dom1}->{rmsd},
                   $aln_cache->{0}->{alistats}->{$seq1}->{$dom1}->{numequivpos},
                   $aln_cache->{0}->{alistats}->{$seq1}->{$dom1}->{numres1},
                   $aln_cache->{0}->{alistats}->{$seq1}->{$dom1}->{numres2},
                   $aln_cache->{0}->{alistats}->{$seq1}->{$dom1}->{numident},
                   $model_con->{bs_seqstats}->{tmpl_numres_bs1},
                   $model_con->{bs_seqstats}->{targ_numres_bs1},
                   $model_con->{bs_seqstats}->{numident_bs1},
                   $aln_cache->{1}->{alistats}->{$seq2}->{$dom2}->{rmsd},
                   $aln_cache->{1}->{alistats}->{$seq2}->{$dom2}->{numequivpos},
                   $aln_cache->{1}->{alistats}->{$seq2}->{$dom2}->{numres1},
                   $aln_cache->{1}->{alistats}->{$seq2}->{$dom2}->{numres2},
                   $aln_cache->{1}->{alistats}->{$seq2}->{$dom2}->{numident},
                   $model_con->{bs_seqstats}->{tmpl_numres_bs2},
                   $model_con->{bs_seqstats}->{targ_numres_bs2},
                   $model_con->{bs_seqstats}->{numident_bs2},
               ) ;

               modtie::pibase::replace_undefs(\@outvals, '') ;

               print {$fh->{binscores}} join("\t", @outvals)."\n" ;
#               print STDERR  " passed " ;
            }
            print STDERR "\n" ;

            $binary_counter++ ;
         }
      }
      close $fh->{binscores} ;
   }

}


=head2 seq_2_arr()

   Title:       seq_2_arr()
   Function:    Gien an aa sequence, returns array of three letter symbols
   Args:        $_ = aa sequence
   Returns:     $_->[] - aa 3-letter abbreviation

=cut

sub seq_2_arr {

   my $seq = shift ;
   my @aa = split('', $seq) ;
   my $one2three = {
      'G' => 'GLY',
      'A' => 'ALA',
      'V' => 'VAL',
      'L' => 'LEU',
      'I' => 'ILE',
      'P' => 'PRO',
      'F' => 'PHE',
      'Y' => 'TYR',
      'W' => 'TRP',
      'C' => 'CYS',
      'M' => 'MET',
      'S' => 'SER',
      'T' => 'THR',
      'K' => 'LYS',
      'R' => 'ARG',
      'H' => 'HIS',
      'D' => 'ASP',
      'E' => 'GLU',
      'N' => 'ASN',
      'Q' => 'GLN'
   } ;

   my $aa3 ;
   foreach my $j ( 0 .. $#aa) { $aa3->[$j] = $one2three->{$aa[$j]} ; }

   return $aa3 ;
}


sub readin_seqdomains_start_modelid {

   my $in = shift ;

   print STDERR "reading in modbase domain definitions: " ;
   my $seqfrag2model ;
   my $modelid_2_target_beg ;
   my ($basemodels_dir, $baseali_dir, $modpipe_style) ;
   open(SEQDOMAINS, $in->{fn}) ;

# See file format in the write routine: assign_seqid_domains()
   while (my $line = <SEQDOMAINS>) {
      chomp $line;
      if ($line =~ /^#/) {next;}
      my @t = split(/\t/, $line) ;
      $seqfrag2model->{$t[0]}->{$t[1]} = $t[3] ;

      $modpipe_style->{$t[3]} = $t[17] ;
      $baseali_dir->{$t[3]} = $t[18] ;
      $basemodels_dir->{$t[3]} = $t[19] ;

      $modelid_2_target_beg->{$t[3]} = $t[4] ;
   }
   close(SEQDOMAINS) ;
   print STDERR "X\n" ;

   return {
      seqfrag2model => $seqfrag2model,
      basemodels_dir => $basemodels_dir,
      baseali_dir => $baseali_dir,
      modpipe_style => $modpipe_style,
      modelid_2_target_beg => $modelid_2_target_beg
   } ;

}


sub readin_model_pdbstart {

   my $in = shift ;

   print STDERR "reading in model PDB start resno: " ;
   my $modelid_2_pdbstart;
   open(MODELPDBSTARTS, $in->{fn}) ;

# See file format in the write routine: assign_seqid_domains()
   while (my $line = <MODELPDBSTARTS>) {
      chomp $line;
      if ($line =~ /^#/) {next;}
      my @t = split(/\t/, $line) ;
      $modelid_2_pdbstart->{$t[0]} = $t[1] ;
   }
   close(MODELPDBSTARTS) ;
   print STDERR "X\n" ;

   return $modelid_2_pdbstart ;

}


sub _parse_domarchs {

   my $params = shift ;

   my $scop2seq ;
   my $mb2dom ;
   my $mbdomfrag ;
   my $type2num ;

   foreach my $seq_id ( keys %{$params->{seqinfo}->{arch}} ) {
      my $cur_arch = $params->{seqinfo}->{arch}->{$seq_id} ;
      my @segments = ($cur_arch =~ /([\<\(\[][^\(\)\[\]\<\>]+[\>\)\]])/g) ;
      my $curres = 1 ;
      foreach my $k ( 0 .. $#segments) {
         if ($segments[$k] =~ /^\(/) {
            my ($num) = ($segments[$k] =~ /([0-9]+)/) ;
            $curres -= $num ;

         } elsif ($segments[$k] =~ /^\</) {

            my ($num) = ($segments[$k] =~ /([0-9]+)/) ;
            $curres += $num ;

         } else {

            $segments[$k] =~ s/^\[// ;
            $segments[$k] =~ s/\]$// ;

            my @t = split(':', $segments[$k]) ;
            my ($domno,$scopcl,$fragno,$segsize,$domsize,$conf) = @t ;

            my ($scopfa, $scopsf) ;
            my @pieces = split(/\./, $scopcl) ;
            if ($#pieces == 3) { $scopfa = $scopcl ; }
            if ($#pieces >= 2) { $scopsf = $pieces[0].'.'.$pieces[1].
                                           '.'.$pieces[2] ; }

            if ($fragno == 1) {
               if (defined $scopfa) {
                  $type2num->{$scopfa}++ ;
                  push @{$scop2seq->{$scopfa}->{$seq_id}}, $domno;
               }

               if (defined $scopsf) {
                  $type2num->{$scopsf}++ ;
                  push @{$scop2seq->{$scopsf}->{$seq_id}}, $domno;
               }

               $mb2dom->{$seq_id}->[$domno] = {
                  scop => $scopcl,
                  numfrag => $fragno,
                  segsize => $segsize,
                  domsize => $domsize
               } ;
            }

            $mbdomfrag->{$seq_id}->[$domno]->[$fragno] = {
               start => $curres,
               end => ($curres + $segsize - 1)
            } ;
            $curres += $segsize ;
         }
      }

      foreach my $curdom ( 1 .. $#{$mbdomfrag->{$seq_id}}) {
         my @cdefs;
         foreach my $curfrag ( 1 .. $#{$mbdomfrag->{$seq_id}->[$curdom]}) {
            push @cdefs, $mbdomfrag->{$seq_id}->[$curdom]->[$curfrag]->{start}.
               '-'.$mbdomfrag->{$seq_id}->[$curdom]->[$curfrag]->{end}
         }
         $mb2dom->{$seq_id}->[$curdom]->{def} = join('_', @cdefs) ;
      }
   }

   my $results = { scop2seq => $scop2seq,
                   seq2dom => $mb2dom,
                   seq2domfrag => $mbdomfrag,
                   scop2num => $type2num } ;

   if (defined $params->{seqset}) {
      my $scop2seqset ;
      foreach my $type (keys %{$scop2seq}) {
         foreach my $seq_id (keys %{$scop2seq->{$type}}) {
            my $setno = $params->{seqset}->{union}->{$seq_id} ;
            if ($setno == 1) {
               push @{$scop2seqset->{1}->{$type}->{$seq_id}},
                  @{$scop2seq->{$type}->{$seq_id}} ;
            } elsif ($setno == 2) {
               push @{$scop2seqset->{2}->{$type}->{$seq_id}},
                  @{$scop2seq->{$type}->{$seq_id}} ;
            } elsif ($setno == 3) {
               push @{$scop2seqset->{1}->{$type}->{$seq_id}},
                  @{$scop2seq->{$type}->{$seq_id}} ;
               push @{$scop2seqset->{2}->{$type}->{$seq_id}},
                  @{$scop2seq->{$type}->{$seq_id}} ;
            }
         }
      }
      $results->{scop2seqset} = $scop2seqset ;
   }

   return $results ;
}


=head2 fy_shuffle()

   Title:       fy_shuffle()
   Function:    Fisher-yates shuffle an array (in-place)
   Args:        $_ = reference to array
   Returns:     NOTHING - shuffles array in place

=cut

sub fy_shuffle {
   my $deck = shift;  # $deck is a reference to an array
   my $i = @$deck;
   while ($i--) {
      my $j = int rand ($i+1);
      @$deck[$i,$j] = @$deck[$j,$i];
   }
}


=head2 calc_zscore()

   Title:       calc_zscore()
   Function:    Fisher-yates shuffle an array (in-place)
   Args:        $_->[0] = target score (scalar)
                $_->[bg] = list of background scores (arrayref)
   Returns:     $->{z_score}
                $->{z_min_tn}
                $->{z_prime}
                $->{z_2}

=cut

sub calc_zscore {

   my $score = shift ;
   my $bg = shift ;

   my $total ;
   my $min ;
   my $mintrueneg ;

   my $results ;

   my $false_pos = 0 ;
   foreach my $j ( 0 .. $#{$bg}) {
      $total += $bg->[$j] ;
      if ($bg->[$j] <= $score) {
         $false_pos++ ; }

      if ( ($bg->[$j] > $score) &&
           ((!exists $results->{min_tn}) || ($bg->[$j] <= $results->{min_tn}))) {
         $results->{min_tn}= $bg->[$j] ; }

      if ((!exists $results->{min}) || ($bg->[$j] < $results->{min})) {
         $results->{min} = $bg->[$j] ; }
   }

   $results->{false_pos} = $false_pos / ($#{$bg} + 1 ) ;
   $results->{avg} = $total / ($#{$bg} + 1) ;

   $results->{stdev}= 0 ;
   foreach my $j ( 0 .. $#{$bg}) {
      $results->{stdev} += (($bg->[$j] - $results->{avg}) *
                           ($bg->[$j] - $results->{avg})) ; }

   $results->{stdev} = $results->{stdev} / ($#{$bg} + 1);
   if ($results->{stdev} == 0 ) { $results->{error} = "stdev of 0" ; return $results; }
   $results->{stdev} = sqrt($results->{stdev});

   $results->{z_bg} = [ ] ;
   foreach my $j ( 0 .. $#{$bg}) {
      push @{$results->{z_bg}}, sprintf("%.3f", (($bg->[$j] - $results->{avg}) / $results->{stdev})) ; }

   $results->{z_score} = ($score - $results->{avg}) / $results->{stdev} ;
   if (exists $results->{min_tn}) {
      $results->{z_min_tn} = ($results->{min_tn} - $results->{avg}) / $results->{stdev} ;
      $results->{z_prime} = $results->{z_score} - $results->{z_min_tn} ;
   } else {
      $results->{z_min_tn} = 'undef' ;
      $results->{z_prime} = 'undef' ;
   }

   $results->{z_min} = ($results->{min} - $results->{avg}) / $results->{stdev} ;
   $results->{z_2} = $results->{z_score} - $results->{z_min} ;

   foreach my $k (keys %{$results}) {
      if ($k ne 'z_bg' && $results->{$k} ne 'undef') {
         $results->{$k} = sprintf("%.3f", $results->{$k});}}

   return $results ;

}


=head2 contact_equiv_complex()

   Title:       contact_equiv_complex()
   Function:    Determine target contacts implied by an alignment to template
                 complex
   Args:        $_->{contacts} = template contacts
                $_->{seq} = target residues
                $_->{resmap} = target residues
                $_->{tmpl_resnames}
                $_->{tmpl_subsers}
                $_->{dompair2respair}
                $_->{labelnames}

   Returns:     $->{bs_seqstats} = {
                  tmpl_numres_bs1
                  tmpl_numres_bs2
                  tmpl_numres_bs12
                  targ_numres_bs1
                  targ_numres_bs2
                  targ_numres_bs12
                  numident_numres_bs1
                  numident_numres_bs2
                  numident_numres_bs12
                }
                $->{respairs}->{dom1."\t".dom2}->{res1sig."\n".res2sig}->
                  {all|SS|MM|MS}->{min_afp|min_dist}
                  min_afp: minimum fraction of atoms participating in contact
                  min_dist: minimum distance between atoms in residue pair

                $->{resnames}->{DOMAIN_ID}->{TARGET_RESSIG} = RESIDUE_NAME
                $->{dompair2respair}->{dom1."\t".dom2}->{res1sig."\n".res2sig}

=cut

sub contact_equiv_complex {

   my $in = shift ;

   my $contacts = $in->{contacts};
   my $targetres = $in->{seq} ;
   my $resmap = $in->{resmap} ;
   my $tmpl_resnames = $in->{tmpl_resnames} ;
   my $tmpl_subsres= $in->{tmpl_subsres} ;
   my $dompair2respair = $in->{dompair2respair};
   my $labelnames = $in->{labelnames} ;

   my $newcon ;
   my $newnames ;
   my $newdp2rp;

   my ($tmpl_bs, $targ_bs) ;
   $tmpl_bs = {} ;
   $tmpl_bs->{numres1} = {} ;
   $tmpl_bs->{numres2} = {} ;
   $targ_bs->{numres1} = {} ;
   $targ_bs->{numres2} = {} ;
   $targ_bs->{numident1} = {} ;
   $targ_bs->{numident2} = {} ;

   $targ_bs = {} ;

   foreach my $dom12sig ( keys %{$contacts}) {
    foreach my $ressig ( keys %{$contacts->{$dom12sig}}) {
      my ($domain_id1, $domain_id2) = split(/\n/, $dom12sig) ;
      if ($dompair2respair->{$dom12sig}->{$ressig} == 2) {
         my $t = $domain_id2; $domain_id2 = $domain_id1; $domain_id1 = $t ; }

      my ($res1, $ch1, $res2, $ch2) = split(/\n/, $ressig) ;
#      print STDERR "looking to align to template $res1 ($ch1) -- $res2 ($ch2)\n" ;

      my $sig1 = $res1."\n".$ch1 ;
      my $sig2 = $res2."\n".$ch2 ;

# go from template residue -> domain -> target residue equivalency (Resmap)

#      print STDERR "looking for equivalent to template $res1 $ch1 -- $res2 $ch2\n";
#      print STDERR "resmap1: ".$resmap->{$domain_id1}->[1]->{$sig1}."\n" ;
#      print STDERR "resmap2: ".$resmap->{$domain_id2}->[1]->{$sig2}."\n" ;

#      my $domain_id1 = $tmpl_subsres->{$sig1} ;
#      my $domain_id2 = $tmpl_subsres->{$sig2} ;

      $tmpl_bs->{numres1}->{$sig1}++ ;
      $tmpl_bs->{numres2}->{$sig2}++ ;

      if (defined $resmap->{$domain_id1}->[1]->{$sig1} &&
          defined $resmap->{$domain_id2}->[1]->{$sig2}) {

         my $trg_ressig_1 = $resmap->{$domain_id1}->[1]->{$sig1} ;
         my $trg_ressig_2 = $resmap->{$domain_id2}->[1]->{$sig2} ;

         my $trg_resna_1 = $targetres->{$domain_id1}->{$trg_ressig_1} ;
         if (!defined $trg_resna_1) {
            my $tt = "WARNING: contact_equiv_complex(): could not find resname for target (".$labelnames->{targ1}.") residue ".$trg_ressig_1." aligned to template (".$labelnames->{tmpl1}." residue $sig1/\n" ;
            $tt =~ s/\n/:/g ;
            print STDERR $tt."\n" ;
            next;
         }

         my $trg_resna_2 = $targetres->{$domain_id2}->{$trg_ressig_2} ;
         if (!defined $trg_resna_2) {
            my $tt = "WARNING: contact_equiv_complex(): could not find resname for target (".$labelnames->{targ2}.") residue ".$trg_ressig_2." aligned to template (".$labelnames->{tmpl2}." residue $sig2/\n" ;
            $tt =~ s/\n/:/g ;
            print STDERR $tt."\n" ;
            next;
         }

         $targ_bs->{numres1}->{$trg_ressig_1}++ ;
         $targ_bs->{numres2}->{$trg_ressig_2}++ ;

         if (exists $tmpl_resnames->{$domain_id1}->{$sig1} &&
             $tmpl_resnames->{$domain_id1}->{$sig1} eq
             $trg_resna_1) { $targ_bs->{numident1}->{$trg_ressig_1}++;}

         if (exists $tmpl_resnames->{$domain_id2}->{$sig2} &&
             $tmpl_resnames->{$domain_id2}->{$sig2} eq
             $trg_resna_2) { $targ_bs->{numident2}->{$trg_ressig_2}++;}

         $newnames->{$domain_id1}->{$trg_ressig_1} = $trg_resna_1 ;
         $newnames->{$domain_id2}->{$trg_ressig_2} = $trg_resna_2 ;

         my $trg_pair = join("\n", ($trg_ressig_1, $trg_ressig_2)) ;
         $newdp2rp->{$dom12sig}->{$trg_pair} =
            $dompair2respair->{$dom12sig}->{$ressig} ;

         my @t_dists ;
         foreach my $dets (keys %{$contacts->{$dom12sig}->{$ressig}}) {

            if ($dets =~ /^S/ && $trg_resna_1 eq 'GLY') {next;}
            if ($dets =~ /S$/ && $trg_resna_2 eq 'GLY') {next;}

            $newcon->{$dom12sig}->{$trg_pair}->{$dets}->{min_afp} =
               $contacts->{$dom12sig}->{$ressig}->{$dets}->{min_afp} ;

            $newcon->{$dom12sig}->{$trg_pair}->{$dets}->{min_dist} =
               $contacts->{$dom12sig}->{$ressig}->{$dets}->{min_dist} ;

            push @t_dists,
               $newcon->{$dom12sig}->{$trg_pair}->{$dets}->{min_dist} ;
         }
         my $min_all = min(\@t_dists) ;
         $newcon->{$dom12sig}->{$trg_pair}->{'all'}->{'min_dist'} = $min_all ;
      }
    }
   }

   my $bs_seqstats ;
   $bs_seqstats->{tmpl_numres_bs1} = keys %{$tmpl_bs->{numres1}} ;
   $bs_seqstats->{tmpl_numres_bs2} = keys %{$tmpl_bs->{numres2}} ;
   $bs_seqstats->{tmpl_numres_bs12} = $bs_seqstats->{tmpl_numres_bs1} +
                                      $bs_seqstats->{tmpl_numres_bs2} ;

   $bs_seqstats->{targ_numres_bs1} = keys %{$targ_bs->{numres1}} ;
   $bs_seqstats->{targ_numres_bs2} = keys %{$targ_bs->{numres2}} ;
   $bs_seqstats->{targ_numres_bs12} = $bs_seqstats->{targ_numres_bs1} +
                                      $bs_seqstats->{targ_numres_bs2} ;

   $bs_seqstats->{numident_bs1} = keys %{$targ_bs->{numident1}} ; 
   $bs_seqstats->{numident_bs2} = keys %{$targ_bs->{numident2}} ; 
   $bs_seqstats->{numident_bs12} = $bs_seqstats->{numident_bs1} +
                                   $bs_seqstats->{numident_bs2} ;

#   foreach my $key (keys %{$bs_seqstats} ) {
#      print STDERR "   contact_equiv $key: $bs_seqstats->{$key}\n" ; }

   my $results ;
   $results = {
      bs_seqstats => $bs_seqstats,
      respairs => $newcon,
      resnames => $newnames,
      dompair2respair => $newdp2rp
   } ;

   return $results ;

# need to change resnames format to allow for duplicate residue numbering

}


=head2 get_subsres()

   Title:       get_subsres()
   Function:    Read in subsets_residues file to get residue2subset_id map
   Args:        $_ = PIBASE subsets_residues table-on-disk filename
   Returns:     $_->{resno."\n".chain_id} = subset_id

=cut

sub get_subsres {

   my $fn = shift ;

   my $assign ;
   ($assign->{subset_id}, $assign->{chain_id}, $assign->{resno_serial},
    $assign->{resno}) = modtie::pibase::rawselect_metatod($fn,
      "SELECT subset_id, chain_id, resno_serial, resno ".
      "FROM $fn") ;

   my $res2sub ;
   foreach my $j ( 0 .. $#{$assign->{subset_id}}) {
      if ($assign->{subset_id}->[$j] =~ /SCOP/) {
         my $sig = $assign->{resno}->[$j]."\n".$assign->{chain_id}->[$j] ;
         $res2sub->{$sig} = $assign->{subset_id}->[$j] ;
      }
   }

   return $res2sub ;

}


=head2 get_domain_contacts_complex()

   Title:       get_domain_contacts_complex()
   Function:    Compute inter-domain residue contacts given a PDB file and defs
   Args:        $_->{subsets_residues}->{RESNO."\n".CHAIN_ID} = SUBSET_ID
                $_->{pdb_fn} = PDB filename
                $_->{R_max} = maximum distance
                $_->{R_list} = list of max dists at which to compute contacts
   Returns:     $_->{respairs}->{sid12sig}->{res12sig}
                     ->{all}->{min_afp|min_dist}
                     ->{MM|MS|SS}->{min_afp}->{DIST_THRESH} = value
                $->{resnames}->{SUBSET_ID}->{RESNO."\n".CHAIN_ID} = RESIDUE NAME
                $->{dompair2respair}->{sid12sig}->{res12sig}
                  = 1 if sid1 .lt. sid2
                  = 2 if sid1 .ge. sid2

=cut

sub get_domain_contacts_complex {
# out of pibase framework - take pdbfn / boundaires / subsres from input
#  if not exists fuckoff
#
# differences from get_domain_Contactxs() :
# 1. takes in full complex structure versus individual domain pdb files

# have one full respair list and one that assigns each respair to a pair of domains

   my $in = shift ;
   my $subsres = $in->{subset_residues} ;
   my $strx_path = $in->{pdb_fn} ;

   my $specs ;
   if (exists $in->{specs}) { $specs = $in->{specs} ;
   } else { $specs = set_modtie_specs() ; }

   my $mcatoms = $specs->{gconst}->{mcatoms};

# number of side chain atoms
   my $numaaatoms = $specs->{gconst}->{numaaatoms};

   my $kdcont_out  = calc_res_pairs({
         radius => $in->{R_max},
         bdp_path => $strx_path
   }) ;

   if (exists $kdcont_out->{error_fl})  {
      print STDERR "ERROR ($strx_path): $kdcont_out->{error_fl}\n" ;
      next;
   }
   my $kdcontacts_out = $kdcont_out->{contacts_fh} ;

   my $respairs ;
   my $respair_atoms ;
   my $resnames;
   my $sidpair2respair;
   while (my $line = <$kdcontacts_out>) {
      chomp $line;
      if ($line =~ /^#/) {next;}

      my ($resna_1, $resno_1, $inscode1, $chain_id_1, undef, $atomna_1,
             $resna_2, $resno_2, $inscode2, $chain_id_2, undef, $atomna_2,
             $dist) = split(/\t/, $line) ;

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

      if (!exists $subsres->{$sig1} ||
          !exists $subsres->{$sig2} ||
          ($subsres->{$sig1} eq $subsres->{$sig2})) { next; }


      my $sid1 = $subsres->{$sig1} ;
      my $sid2 = $subsres->{$sig2} ;

#      print STDERR " $sid1 -- $sid2: $resno_1 $chain_id_1 -- $resno_2 $chain_id_2\n" ;

      my $ressig = $sig1."\n".$sig2 ;

      my $sc_type1 = 's'; my $sc_type2 = 's';
      if (exists $mcatoms->{$atomna_1}) { $sc_type1 = 'm' ; }
      if (exists $mcatoms->{$atomna_2}) { $sc_type2 = 'm' ; }

      my $sc_type = $sc_type1.$sc_type2 ;

      if ($resna_1 eq 'HSD') {$resna_1 = 'HIS'; }
      if ($resna_2 eq 'HSD') {$resna_2 = 'HIS'; }

      $resnames->{$sid1}->{$sig1} = $resna_1 ;
      $resnames->{$sid2}->{$sig2} = $resna_2 ;
      my $sid12sig ;
      if ($sid1 lt $sid2) {
         $sid12sig = $sid1."\n".$sid2 ;
         $sidpair2respair->{$sid12sig}->{$ressig} = 1 ;
      } else {
         $sid12sig = $sid2."\n".$sid1 ;
         $sidpair2respair->{$sid12sig}->{$ressig} = 2 ;
      }

      if (!exists $respair_atoms->{$sid12sig}->{$ressig}->{$sc_type}->[0]->{$atomna_1}->{min_dist} ||
         ($dist <= $respair_atoms->{$sid12sig}->{$ressig}->{$sc_type}->[0]->{$atomna_1}->{min_dist})) {
         $respair_atoms->{$sid12sig}->{$ressig}->{$sc_type}->[0]->{$atomna_1}->{min_dist} = $dist; }

      if (!exists $respair_atoms->{$sid12sig}->{$ressig}->{$sc_type}->[1]->{$atomna_2}->{min_dist} ||
         ($dist <= $respair_atoms->{$sid12sig}->{$ressig}->{$sc_type}->[1]->{$atomna_2}->{min_dist})) {
         $respair_atoms->{$sid12sig}->{$ressig}->{$sc_type}->[1]->{$atomna_2}->{min_dist} = $dist; }


      if (!exists $respair_atoms->{$sid12sig}->{$ressig}->{all}->[0]->{$atomna_1}->{min_dist} ||
         ($dist <= $respair_atoms->{$sid12sig}->{$ressig}->{all}->[0]->{$atomna_1}->{min_dist})) {
         $respair_atoms->{$sid12sig}->{$ressig}->{all}->[0]->{$atomna_1}->{min_dist} = $dist; }

      if (!exists $respair_atoms->{$sid12sig}->{$ressig}->{all}->[1]->{$atomna_2}->{min_dist} ||
          ($dist <= $respair_atoms->{$sid12sig}->{$ressig}->{all}->[1]->{$atomna_2}->{min_dist})) {
         $respair_atoms->{$sid12sig}->{$ressig}->{all}->[1]->{$atomna_2}->{min_dist} = $dist; }


      if ((!exists $respairs->{$sid12sig}->{$ressig}->{all}->{min_dist}) ||
         ($dist <= $respairs->{$sid12sig}->{$ressig}->{all}->{min_dist})) {
                   $respairs->{$sid12sig}->{$ressig}->{all}->{min_dist} = $dist ; }

      if ((!exists $respairs->{$sid12sig}->{$ressig}->{$sc_type}->{min_dist}) ||
         ($dist <= $respairs->{$sid12sig}->{$ressig}->{$sc_type}->{min_dist})) {
                     $respairs->{$sid12sig}->{$ressig}->{$sc_type}->{min_dist} = $dist ; }
   }
   close $kdcontacts_out ;

   my $f_respairs ;
   my $f_resnames;
   my $f_sidpair2respair;
   my $f_sidpair2respair_afp;
   foreach my $sid12sig (keys %{$sidpair2respair}) {
      my @t_doms = split(/\n/, $sid12sig) ;
      foreach my $ressig (keys %{$sidpair2respair->{$sid12sig}}) {

         my ($tr1, $tc1, $tr2, $tc2) = split(/\n/, $ressig) ;
         my $curres = [$tr1."\n".$tc1, $tr2."\n".$tc2] ;


         my $domsigorder = $sidpair2respair->{$sid12sig}->{$ressig} ;
         my $curdoms ;
         if ($domsigorder == 2) {
            $curdoms->[0] = $t_doms[1] ;
            $curdoms->[1] = $t_doms[0] ;
         } else {
            $curdoms->[0] = $t_doms[0] ;
            $curdoms->[1] = $t_doms[1] ;
         }


         foreach my $sc_type (keys %{$respairs->{$sid12sig}->{$ressig}}) {

            my $sc1 = 'all'; my $sc2 = 'all';
            if ($sc_type ne 'all') {
               $sc1 = substr($sc_type, 0, 1) ;
               $sc2 = substr($sc_type, 1, 1) ;
            }
            my @scs = ($sc1, $sc2) ;

            my $t_resna1 = $resnames->{$curdoms->[0]}->{$curres->[0]} ;
            my $t_resna2 = $resnames->{$curdoms->[1]}->{$curres->[1]} ;

            foreach my $cur_r (@{$in->{R_list}})  {
               my $t_intatoms1 = 0 ;
               foreach my $atom (keys %{$respair_atoms->{$sid12sig}->{$ressig}->{$sc_type}->[0]}) {
                  if ($respair_atoms->{$sid12sig}->{$ressig}->{$sc_type}->[0]->{$atom}->{min_dist} <= $cur_r) {
                        $t_intatoms1++ ; } }
               my $t_afp1 = $t_intatoms1 / $numaaatoms->{$sc1}->{$t_resna1} ;

               my $t_intatoms2 = 0 ;
               foreach my $atom (keys %{$respair_atoms->{$sid12sig}->{$ressig}->{$sc_type}->[1]}) {
                  if ($respair_atoms->{$sid12sig}->{$ressig}->{$sc_type}->[1]->{$atom}->{min_dist} <= $cur_r) {
                        $t_intatoms2++ ; } }
               my $t_afp2 = $t_intatoms2 / $numaaatoms->{$sc2}->{$t_resna2} ;

               my $minafp = $t_afp2; if ($t_afp1 < $t_afp2) {$minafp=$t_afp1;}

               $respairs->{$sid12sig}->{$ressig}->{$sc_type}->{min_afp}->{$cur_r} = 
                $minafp ;
            }
         }
      }
   }


   my $results;
   $results = {
      respairs => $respairs,
      resnames => $resnames,
      dompair2respair => $sidpair2respair
   } ;

   return $results;
}


=head2 SUB calc_res_pairs()

   Title:       calc_res_pairs()
   Function:    Given a PDB file, calculate all residue pairs
   Args:        $_->{radius} - maximum distance for kdcontacts run
                $_->{bdp_path} - location of PDB file
   Returns:     $->{contacts_fh} = filehandle to kdcontacts output

=cut

sub calc_res_pairs {

   my $params = shift ;

   my $kdcontacts_radius = $params->{radius} || "6.6" ;

   my $binaries = locate_binaries() ;
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

# send back filehandle to kdcontacts output
# read in kdcontacts output into an array

   my $results = { contacts_fh => $out_fh };

   return $results ;
}



=head2 SUB readin_potentials()

   Title:       readin_potentials()
   Function:    Read in MODTIE or SMY potential
   Args:        $_->{pots} = [{
                  type => 'smy' or 'pi'
                  fn => potential file path
                }]

   Returns:     populates an additional hash entry above:
                $_->{pots}->[i]->{val}->{potential_entry_type}->{aa1.aa2} = val

=cut

sub readin_potentials {

   my $pots = shift ;

   foreach my $pot ( @{$pots}) {

      if ($pot->{type} eq 'smy') {
         _readpot_smy($pot) ;
      } elsif ($pot->{type} =~ /^pi/) {
         _readpot_pi($pot) ;
      } else {
         $pot->{error} = 'unrecognized type' ;
      }
   }

   return $pots ;

}


=head2 get_rand_score_complex_alascan()

   Title:       get_rand_score_complex_alascan()
   Function:    randomize residues and recompute raw MODTIE potential scores
                to calculate background for ALASCAN scoring

   Args:        $_->{contacts}
                $_->{pot}
                $_->{resnames}
                $_->{dompair2respair}

   Returns:     $_->{whole} = [] whole complex raw scores
                $_->{interface}->{SID12SIG} = [] raw score for each interface

=cut

sub get_rand_score_complex_alascan {
   my $in = shift ;

   my $contacts = $in->{contacts} ;
   my $pot = $in->{pot} ;
   my $resna = $in->{resnames} ;
   my $dompair2respair= $in->{dompair2respair} ;

   my $ranscores ;
   foreach my $j ( 0 .. 999)  {
      my $ranresna ;
      foreach my $domain ( keys %{$resna}) {
         my @ranresna = values(%{$resna->{$domain}}) ;
         my @ressig = keys(%{$resna->{$domain}});
         fy_shuffle(\@ranresna) ;
         @{$ranresna->{$domain}}{@ressig} = @ranresna ;
      }

      my $ranscore = score_potential_complex_alascan({
         contacts => $contacts,
         pot => $pot,
         resnames => $ranresna,
         dompair2respair => $dompair2respair
      }) ;
      push @{$ranscores->{whole}}, $ranscore->{whole} ;

      foreach my $int (keys %{$ranscore->{interface}}) {
         push @{$ranscores->{interface}->{$int}},
            $ranscore->{interface}->{$int}; }
   }

   return $ranscores ;
}


=head2 get_rand_score_complex()

   Title:       get_rand_score_complex()
   Function:    randomize residues and recompute raw MODTIE potential scores
                for z-scores scoring

   Args:        $_->{contacts}
                $_->{pot}
                $_->{resnames}
                $_->{dompair2respair}

   Returns:     $_->{whole} = [] whole complex raw scores
                $_->{interface}->{SID12SIG} = [] raw score for each interface

=cut

sub get_rand_score_complex {
   my $in = shift ;

   my $contacts = $in->{contacts} ;
   my $pot = $in->{pot} ;
   my $resna = $in->{resnames} ;
   my $dompair2respair= $in->{dompair2respair} ;

   my $ranscores ;
   foreach my $j ( 0 .. 999)  {
      my $ranresna ;
      foreach my $domain ( keys %{$resna}) {
         my @ranresna = values(%{$resna->{$domain}}) ;
         my @ressig = keys(%{$resna->{$domain}});
         fy_shuffle(\@ranresna) ;
         @{$ranresna->{$domain}}{@ressig} = @ranresna ;
      }

      my $ranscore = score_potential_complex({
         contacts => $contacts,
         pot => $pot,
         resnames => $ranresna,
         dompair2respair => $dompair2respair
      }) ;
      push @{$ranscores->{whole}}, $ranscore->{whole} ;

      foreach my $int (keys %{$ranscore->{interface}}) {
         push @{$ranscores->{interface}->{$int}},
            $ranscore->{interface}->{$int}; }
   }

   return $ranscores ;
}


=head2 score_potential_complex_alascan()

   Title:       score_potential_complex_alascan()
   Function:    ALASCAN scores a complex given pairwise residue contacts and
                 a potential. WARNING NOT BENCHMARKED

   Args:        $_->{contacts}
                $_->{pot}
                $_->{resnames}
                $_->{dompair2respair}

   Returns:     $_->{whole} = raw score for whole complex
                $_->{interface}->{SID12SIG} = raw score for each interface

=cut

sub score_potential_complex_alascan {

   my $in = shift ;
   my $contacts = $in->{contacts} ;
   my $dompair2respair = $in->{dompair2respair} ;
   my $resna = $in->{resnames} ;
   my $pot = $in->{pot} ;

   my $score = {} ;

   if ($pot->{type} eq 'pi') {

      my $dets = $pot->{details} ;
      my $r_thresh = $pot->{r} ;

      $score->{whole} = 0 ;

      foreach my $dom12sig ( keys %{$contacts}) {
         $score->{interface}->{$dom12sig} = 0 ;
         my ($dom1, $dom2) = split(/\n/, $dom12sig) ;

         foreach my $ressig ( keys %{$contacts->{$dom12sig}}) {
            my ($dom1, $dom2) = split(/\n/, $dom12sig) ;
            if ($dompair2respair->{$dom12sig}->{$ressig} == 2) {
               my $t = $dom2; $dom2 = $dom1; $dom1 = $t ; }

            my $revdets = 'ohshit' ;
            if ($dets eq 'ms') { $revdets = 'sm'; }
            if ((exists $contacts->{$dom12sig}->{$ressig}->{$dets} &&
                 $contacts->{$dom12sig}->{$ressig}->{$dets}->{min_dist} <=
                  $r_thresh) || 
                (exists $contacts->{$dom12sig}->{$ressig}->{$revdets} &&
                $contacts->{$dom12sig}->{$ressig}->{$revdets}->{min_dist} <=
                  $r_thresh)) {
# HAVE TO COUNT MS and SM
               my ($res1, $ch1, $res2, $ch2 ) = split(/\n/, $ressig) ;

               my $resna1 = $resna->{$dom1}->{$res1."\n".$ch1} ;
               my $resna2 = $resna->{$dom2}->{$res2."\n".$ch2} ;
#               print STDERR " the residue names are: $resna1 -- $resna2\n" ;
               if (!defined $resna1) {
                  print STDERR " *** undefined resna1 for ".
                     "$res1 on chain $ch1 (domain $dom1)/\n" ;}
               if (!defined $resna2) {
                  print STDERR " *** undefined resna2 for ".
                     "$res2 on chain $ch2 (domain $dom2)\n" ;}

               if ($resna1 eq 'HSD') { $resna1 = 'HIS';}
               if ($resna2 eq 'HSD') { $resna2 = 'HIS';}

               if ($dets =~ /^S/ && $resna1 eq 'GLY') {next;}
               if ($dets =~ /S$/ && $resna2 eq 'GLY') {next;}

               if (exists $contacts->{$dom12sig}->{$ressig}->{$dets} &&
                   exists $pot->{val}->{w_ij}->{$resna1.$resna2}) {

                  my $thisscore =$contacts->{$dom12sig}->{$ressig}->{$dets}->{min_afp}->{$r_thresh} * $pot->{val}->{w_ij}->{$resna1.$resna2} ; 
                  my $thisscore_ala1 = $contacts->{$dom12sig}->{$ressig}->{$dets}->{min_afp}->{$r_thresh} * $pot->{val}->{w_ij}->{"ALA".$resna2} ; 
                  my $thisscore_ala2 = $contacts->{$dom12sig}->{$ressig}->{$dets}->{min_afp}->{$r_thresh} * $pot->{val}->{w_ij}->{$resna1."ALA"} ; 

                  $score->{interface}->{$dom12sig} += $thisscore;
                  $score->{interface_perres_partner}->{$dom12sig}->{$res1."\t".$ch1}->{$res2."\t".$ch2} += $thisscore ;
                  $score->{interface_perres_partner}->{$dom12sig}->{$res2."\t".$ch2}->{$res1."\t".$ch1} += $thisscore ;
                  $score->{interface_perres_partner_ala}->{$dom12sig}->{$res1."\t".$ch1}->{$res2."\t".$ch2} += $thisscore_ala1;
                  $score->{interface_perres_partner_ala}->{$dom12sig}->{$res2."\t".$ch2}->{$res1."\t".$ch1} += $thisscore_ala2 ;

                  $score->{interface_perres}->{$dom12sig}->{$res1."\t".$ch1} += $thisscore ;
                  $score->{interface_perres}->{$dom12sig}->{$res2."\t".$ch2} += $thisscore ;
                  $score->{interface_perres_ala}->{$dom12sig}->{$res1."\t".$ch1} += $thisscore_ala1 ;
                  $score->{interface_perres_ala}->{$dom12sig}->{$res2."\t".$ch2} += $thisscore_ala2 ;
               }

               if (exists $contacts->{$dom12sig}->{$ressig}->{$revdets} &&
                   exists $pot->{val}->{w_ij}->{$resna2.$resna1}) {

                   my $thisscore =
      $contacts->{$dom12sig}->{$ressig}->{$revdets}->{min_afp}->{$r_thresh} *
      $pot->{val}->{w_ij}->{$resna2.$resna1} ; 

                   my $thisscore_ala1 =
      $contacts->{$dom12sig}->{$ressig}->{$revdets}->{min_afp}->{$r_thresh} *
      $pot->{val}->{w_ij}->{$resna2."ALA"} ; 

                   my $thisscore_ala2 =
      $contacts->{$dom12sig}->{$ressig}->{$revdets}->{min_afp}->{$r_thresh} *
      $pot->{val}->{w_ij}->{"ALA".$resna1} ; 

                  $score->{interface}->{$dom12sig} += $thisscore ;

                  $score->{interface_perres_partner}->{$dom12sig}->{$res1."\t".$ch1}->{$res2."\t".$ch2} += $thisscore ;
                  $score->{interface_perres_partner}->{$dom12sig}->{$res2."\t".$ch2}->{$res1."\t".$ch1} += $thisscore ;
                  $score->{interface_perres_partner_ala}->{$dom12sig}->{$res1."\t".$ch1}->{$res2."\t".$ch2} += $thisscore_ala1;
                  $score->{interface_perres_partner_ala}->{$dom12sig}->{$res2."\t".$ch2}->{$res1."\t".$ch1} += $thisscore_ala2 ;

                  $score->{interface_perres}->{$dom12sig}->{$res1."\t".$ch1} +=
                     $thisscore ;
                  $score->{interface_perres}->{$dom12sig}->{$res2."\t".$ch2} +=
                     $thisscore ;
               $score->{interface_perres_ala}->{$dom12sig}->{$res1."\t".$ch1} +=
                     $thisscore_ala1 ;
               $score->{interface_perres_ala}->{$dom12sig}->{$res2."\t".$ch2} +=
                     $thisscore_ala2 ;
               }
            }
         }

         $score->{whole} += $score->{interface}->{$dom12sig} ;

         $score->{interface}->{$dom12sig} =
            sprintf("%.3f", $score->{interface}->{$dom12sig}) ;
      }

      $score->{whole} = sprintf("%.3f", $score->{whole}) ;

   } else {
      print STDERR "ERROR: potential type $pot->{type} not recognized\n" ;
   }

   return $score ;
}


=head2 score_potential_complex()

   Title:       score_potential_complex()
   Function:    Scores a complex using pairwise residue contacts and
                 a potential

   Args:        $_->{contacts}
                $_->{pot}
                $_->{resnames}
                $_->{dompair2respair}

   Returns:     $_->{whole} = raw score for whole complex
                $_->{interface}->{SID12SIG} = raw score for each interface

=cut

sub score_potential_complex {

   my $in = shift ;
   my $contacts = $in->{contacts} ;
   my $dompair2respair = $in->{dompair2respair} ;
   my $resna = $in->{resnames} ;
   my $pot = $in->{pot} ;

   my $score = {} ;

   if ($pot->{type} eq 'pi') {

      my $dets = $pot->{details} ;
      my $r_thresh = $pot->{r} ;

      $score->{whole} = 0 ;

      foreach my $dom12sig ( keys %{$contacts}) {
         $score->{interface}->{$dom12sig} = 0 ;
         $score->{interface_percontact}->{$dom12sig} = {} ;
         my ($dom1, $dom2) = split(/\n/, $dom12sig) ;

         foreach my $ressig ( keys %{$contacts->{$dom12sig}}) {
            my ($dom1, $dom2) = split(/\n/, $dom12sig) ;
            if ($dompair2respair->{$dom12sig}->{$ressig} == 2) {
               my $t = $dom2; $dom2 = $dom1; $dom1 = $t ; }

            my $revdets = 'ohshit' ;
            if ($dets eq 'ms') { $revdets = 'sm'; }
            if ((exists $contacts->{$dom12sig}->{$ressig}->{$dets} &&
                 $contacts->{$dom12sig}->{$ressig}->{$dets}->{min_dist} <=
                  $r_thresh) || 
                (exists $contacts->{$dom12sig}->{$ressig}->{$revdets} &&
                $contacts->{$dom12sig}->{$ressig}->{$revdets}->{min_dist} <=
                  $r_thresh)) {
# HAVE TO COUNT MS and SM
               my ($res1, $ch1, $res2, $ch2 ) = split(/\n/, $ressig) ;

               my $resna1 = $resna->{$dom1}->{$res1."\n".$ch1} ;
               my $resna2 = $resna->{$dom2}->{$res2."\n".$ch2} ;
#               print STDERR " the residue names are: $resna1 -- $resna2\n" ;
               if (!defined $resna1) {
                  print STDERR " *** undefined resna1 for $res1 on chain $ch1 (domain $dom1)/\n" ;}
               if (!defined $resna2) {
                  print STDERR " *** undefined resna2 for $res2 on chain $ch2 (domain $dom2)\n" ;}

#               print STDERR "res1 = $res1 $ch1 on dom $dom1 ($resna1) \n" ;
#               print STDERR "res2 = $res2 $ch2 on dom $dom2 ($resna2) \n" ;

               if ($resna1 eq 'HSD') { $resna1 = 'HIS';}
               if ($resna2 eq 'HSD') { $resna2 = 'HIS';}

               if ($dets =~ /^S/ && $resna1 eq 'GLY') {next;}
               if ($dets =~ /S$/ && $resna2 eq 'GLY') {next;}

               if (exists $contacts->{$dom12sig}->{$ressig}->{$dets} &&
                   exists $pot->{val}->{w_ij}->{$resna1.$resna2}) {

                   $score->{interface_percontact}->{$dom12sig}->{$ressig."\n".$resna1."\n".$resna2} += 
                     $contacts->{$dom12sig}->{$ressig}->{$dets}->{min_afp}->{$r_thresh} * 
                     $pot->{val}->{w_ij}->{$resna1.$resna2} ; 

                   $score->{interface}->{$dom12sig} +=
                     $contacts->{$dom12sig}->{$ressig}->{$dets}->{min_afp}->{$r_thresh} * 
                     $pot->{val}->{w_ij}->{$resna1.$resna2} ; 
               }

               if (exists $contacts->{$dom12sig}->{$ressig}->{$revdets} &&
                   exists $pot->{val}->{w_ij}->{$resna2.$resna1}) {

                   $score->{interface_percontact}->{$dom12sig}->{$ressig."\n".$resna1."\n".$resna2} += 
                     $contacts->{$dom12sig}->{$ressig}->{$revdets}->{min_afp}->{$r_thresh} * 
                     $pot->{val}->{w_ij}->{$resna2.$resna1} ; 

                   $score->{interface}->{$dom12sig} +=
                     $contacts->{$dom12sig}->{$ressig}->{$revdets}->{min_afp}->{$r_thresh} * 
                     $pot->{val}->{w_ij}->{$resna2.$resna1} ; 
               }
            }
         }

         $score->{whole} += $score->{interface}->{$dom12sig} ;

         $score->{interface}->{$dom12sig} =
            sprintf("%.3f", $score->{interface}->{$dom12sig}) ;
      }

      $score->{whole} = sprintf("%.3f", $score->{whole}) ;

   } else {
      print STDERR "ERROR: potential type $pot->{type} not recognized\n" ;
   }

   return $score ;
}


sub _readpot_smy {

   my $pot = shift;
   $pot->{params}->{dmin} = 0.75 ;
   $pot->{params}->{dbin} = 0.5 ;
   $pot->{params}->{dmax} = 8.0 ;

   open(POTF, "statpotaa.smy") ;
   my $j = 0;
   my $skip = 0 ;
   my ($res1, $atm1, $res2, $atm2) ;
   my @restype = qw/X ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE/ ;
   push @restype,qw/LEU LYS MET PHE PRO SER THR TRP TYR VAL/ ;
   my @atmtype = ('X', ' N  ', ' CA ', ' C  ', ' O  ') ;
   while (my $line = <POTF>) {
      chomp $line;
      $line =~ s/^\s+// ;
      if (($j % 2) == 0) {
         ($res1, $atm1, $res2, $atm2) = split(' ', $line) ;
         if ($atm1 > 4 || $atm2 > 4) {
            $skip = 1;
         } elsif ( exists $pot->{val}->{$restype[$res1]."\n".$atmtype[$atm1]}->{$restype[$res2]."\n".$atmtype[$atm2]}) {
            $skip = 1 ;
         } else {
            $skip = 0;
         }
      } elsif (!$skip ) {
            my @t = split(' ', $line) ;
            $pot->{val}->{$restype[$res1]."\n".$atmtype[$atm1]}->{$restype[$res2]."\n".$atmtype[$atm2]} = \@t ;
            $pot->{val}->{$restype[$res2]."\n".$atmtype[$atm2]}->{$restype[$res1]."\n".$atmtype[$atm1]} = \@t ;
      }
      $j++ ;
   }
}


# adjusted to new fractional contribution
sub _readpot_pi {
   my $pot = shift;

   my $matrix_fn = $pot->{fn} ;
   if (! -e $matrix_fn) {
      die "can't open MATRIX $matrix_fn\n" ; }

   open(MAT, $matrix_fn) ;
   while (my $line = <MAT>) {
    if ($line !~ /^#/ && $line !~ /^$/) {
      chomp $line;
      my ($type,$sctype, $radius, $aa1, $aa2, $val)= split(' ', $line) ;
      $pot->{val}->{$type}->{$aa1.$aa2} = $val ;
    }
   }
   close(MAT) ;
}


=head2 salign_resequiv_parse()

   Title:       salign_resequiv_parse()
   Function:    Parse SALIGN output for target-template residue equivalences

   Args:        $_->{fn} = alignment file
   Returns:     $_->{resequiv}->[0|1]->{resno."\n".chain} = resno."\n".chain
                  points from sequence[0|1] residue to aligned residue

                $_->{resnames}->[0|1]->{resno."\n".chain} = residue name

                $_->{alistats}->{
                  rmsd
                  numequivpos
                  drms
                  numres1
                  numres2
                  numident
                  percseqid
                  quality_score
                }

=cut

sub salign_resequiv_parse {

   my $rescode = {
      'A' => 'ALA',
      'R' => 'ARG',
      'N' => 'ASN',
      'D' => 'ASP',
      'C' => 'CYS',
      'Q' => 'GLN',
      'E' => 'GLU',
      'G' => 'GLY',
      'H' => 'HIS',
      'I' => 'ILE',
      'L' => 'LEU',
      'K' => 'LYS',
      'M' => 'MET',
      'F' => 'PHE',
      'P' => 'PRO',
      'S' => 'SER',
      'T' => 'THR',
      'W' => 'TRP',
      'Y' => 'TYR',
      'V' => 'VAL'
   } ;


   my $in = shift ;
   my $fn = $in->{fn} ;

   my $com ;
   my $bins = locate_binaries();
   if ($fn =~ /gz$/) { $com = $bins->{'zcat'}." $fn" ;
   } else { $com = "cat $fn" ; }

   open(FILE, "$com |") ;

   my $imin = 0 ;
   my $resequiv ;
   my $resnames ;
   while (my $line = <FILE>) {
      chomp $line;

      if ($imin) {
         if ($line =~ /--------/) {$imin = 0 ; last;}
         my $linelength = length($line) ;
         if ($linelength < 41) {next;}

         my ($resna1) = substr($line, 39, 1) ;
         my ($resno1) = substr($line, 40, 4) ; $resno1 =~ s/ //g ;

         my $chain1 = ' ';
         if ($linelength >= 46) {
            ($chain1) = substr($line, 45, 1) ; }

#MODELLER output now seems to remove trailing blanks in log files.
# make sure line is long enough before try to slurp into a variable.

         my ($resna2, $resno2, $chain2) ;
         if ( length($line) <52 ||
              (substr($line, 47, 1) eq '-' )) {
            $resna2 = '-' ;
         } else {
            ($resna2) = substr($line, 47, 1) ;
            ($resno2) = substr($line, 48, 4) ; $resno2 =~ s/ //g ;

            if (exists $in->{target_beg}) {
               $resno2 = $resno2 + $in->{target_beg} - 1 ; }

            $chain2 = ''; #if strx has no chain id, salign output doesnt have that column
# what is the null convention here HERENOW - blank, undef, or ''?
            if (length($line) >= 54) {
               ($chain2) = substr($line, 53, 1) ; }
         }

         if ($resna2 ne '-' && $resna2 ne ' ') {
            my $full_resna2 =  $rescode->{$resna2} ;
            $resnames->[1]->{$resno2."\n".$chain2} = $full_resna2 ;
         }

         if ($resna1 ne '-' && $resna1 ne ' ') {
            my $full_resna1 =  $rescode->{$resna1} ;
            $resnames->[0]->{$resno1."\n".$chain1} = $full_resna1 ;
         }

         if (($resna1 ne '-') && ($resna2 ne '-') &&
             ($resna1 ne ' ') && ($resna2 ne ' ')) {
            $resequiv->[1]->{$resno1."\n".$chain1} = $resno2."\n".$chain2 ;
            $resequiv->[0]->{$resno2."\n".$chain2} = $resno1."\n".$chain1 ;
         }
      }

      if ($line =~ /------/) {$imin = 1; }
   }

   my ($rmsd, $numequivpos, $drms, $numres1,
       $numres2, $numident, $percseqid, $quality_score);
   while (my $line = <FILE>) {
      chomp $line;
      if ($line =~ /QUALITY_SCORE \(percentage\)/) {
         $line =~ s/.*\:\s+// ;
         $quality_score = $line ;
      } elsif ($line =~ /Upper = RMS, Lower = numb equiv positions/) {
         $line = <FILE> ; #blank line
         $line = <FILE> ; #header
         $line = <FILE> ;
         (undef, undef, $rmsd) = split(' ', $line) ;
         $line = <FILE> ;
         (undef, $numequivpos, undef) = split(' ', $line) ;
      } elsif ($line =~ /Upper = Distance RMS, Lower = numb equiv distances/){
         $line = <FILE> ; #blank line
         $line = <FILE> ; #header
         $line = <FILE> ;
         (undef, undef, $drms) = split(' ', $line) ;
      } elsif ($line =~ /Diag=numb res, Upper=numb equiv res, Lower = % seq/){
         $line = <FILE> ;#blank
         $line = <FILE> ;#blank
         $line = <FILE> ;#header
         $line = <FILE> ;#
         (undef, $numres1, $numident) = split(' ', $line) ;
         $line = <FILE> ;
         (undef, $percseqid, $numres2) = split(' ', $line) ;
         last;
      }
   }
   close(FILE) ;

   my $alistats = {
      rmsd => $rmsd,
      numequivpos => $numequivpos,
      drms => $drms,
      numres1 => $numres1,
      numres2 => $numres2,
      numident => $numident,
      percseqid => $percseqid,
      quality_score => $quality_score
   };

   return {
      resequiv => $resequiv,
      resnames => $resnames,
      alistats => $alistats,
   } ;
}


=head2 min()

   Title:       min()
   Function:    Finds minimum value in an array
   Args:        $_ = arrayref
   Returns:     $_ = lowest value in the array

=cut

sub min {
   my $arr = shift ;
   my $min = $arr->[0] ;
   foreach my $j ( 1 .. $#{$arr}) {if ($arr->[$j] < $min) {$min = $arr->[$j];}}
   return $min ;
}


=head2 salign_targ_tmpl_domains()

   Title:       salign_targ_tmpl_domains()
   Function:    dispatch template/model domain salign jobs on the cluster 
   Args:        Nothing
   Input:       tab-delimited STDIN:
                1. subset_id
                2. model_id
                3. domain definition

   Returns:     Nothing
   Output:      Calculated alignments are placed in alignment directory
                  as determined by sid_modelid_2_alidir()

=cut

sub salign_targ_tmpl_domains {

   my $in = shift ; my $specs ;
   if (exists $in->{specs}) { $specs = $in->{specs} ; }
   else { $specs = set_modtie_specs() ; }

   my $modeller_bin = $specs->{binaries}->{modeller} ;

   my $movethreshold = 200 ;
   my $tempdir = tempdir(CLEANUP => 1) ;
   chdir $tempdir ;
   my @movethese ;
   my @movedest ;


   while (my $line = <STDIN>) {
      chomp $line;
      my ($sid, $model, $def) = split(/\t/, $line) ;
      print STDERR "now on $sid - $model ($def)\n" ;

      my ($out_fh, $out_fn) ;
      $out_fn = "align.$sid.$model.$def.salign";

      my $xpath = modtie::pibase::sid_2_domdir({
                     sid => $sid,
                     specs => $specs
                  })."/$sid.pdb.gz";
      my $mpath = modelid_2_domdir({
                     modelid => $model,
                     specs => $specs,
                  })."/$model.$def.pdb.gz" ;
      my $outdir = sid_modelid_2_alidir({
                     sid => $sid,
                     modelid =>$model,
                     specs => $specs});
      my $finalpath = $outdir.'/'.$out_fn.'.gz' ;

      if (-s $finalpath && -s $finalpath >= 110) {
            print STDERR "done already\n"; next; }

      my ($fn, $fh) ;

      my $mbase =$mpath;
      $mbase =~ s/.*\///g ;
      modtie::pibase::safe_copy($mpath, "./$mbase") ;

      my $xbase =$xpath;
      $xbase =~ s/.*\///g ;
      modtie::pibase::safe_copy($xpath, "./$xbase") ;

      my $out_salign= modtie::pibase::get_salign({
         modeller_bin => $modeller_bin,
         pdb_fn_1 => $xbase,
         pdb_fn_2 => $mbase,
      }) ;


      if (exists $out_salign->{error}) {
         print STDERR "ERROR ($sid, $model, $def): SALIGN error\n" ;
         if (exists $out_salign->{ali_log} && -e $out_salign->{ali_log}) {
            unlink $out_salign->{ali_log} ; }
         next;
      }

      open($out_fh, ">$out_fn") ;

      open (SALIGNLOG, $out_salign->{ali_log}) ;
      my $gotin = 0 ;
      my $quality_score ;
      while (my $logline = <SALIGNLOG>) {
         if ($logline =~ /^Current/) {$gotin++;}
         if ($gotin == 3) {print $out_fh $logline ;}
         if ($logline =~ /QUALITY_SCORE \(percentage\)/) {
            $logline =~ s/.*\:\s+// ;
            $quality_score = $logline ;
         }
      }
      close(SALIGNLOG) ;
      close($out_fh) ;
      unlink $out_salign->{ali_log};  unlink $out_salign->{ali_top}; 

      if (!defined $quality_score ||
          $quality_score < $specs->{salign_quality_score_thresh}) {

         if (!defined $quality_score) {$quality_score = 'UNK';}
         print STDERR "WARNING ($sid, $model, $def): ".
                      "poor SALIGN quality $quality_score\n" ;
          unlink $out_fn ;
      } else {
         system("gzip $out_fn") ; $out_fn .= '.gz' ;
         push @movethese, $out_fn ;
         push @movedest, $outdir ;
         my @out = ($sid, $model, $def, "$outdir/$out_fn", $quality_score) ;
         print join("\t", @out)."\n" ;
      }


      if ($#movethese >= $movethreshold) {
         foreach my $j ( 0 .. $#movethese) {
            if (! -d $movedest[$j]) {
               mkpath($movedest[$j]) ; }
            modtie::pibase::safe_move($movethese[$j], $movedest[$j]) ;}
         @movethese = () ;
         @movedest = () ;
      }

      if ($xbase ne $xpath) {unlink $xbase;}
      if ($mbase ne $mpath) {unlink $mbase;}
   }

   foreach my $j ( 0 .. $#movethese) {
      if (! -d $movedest[$j]) {
         mkpath($movedest[$j]) ; }
      modtie::pibase::safe_move($movethese[$j], $movedest[$j]) ;}

}



=head2 cut_domains()

   Title:       cut_domains()
   Function:    Given (subset_id, model) domain pairs, will prepare input for
                  cutdom.calc.pl to extract the domains from source PDB file
   Args:        Nothing
   Input:       tab-delimited STDIN:
                1. BDP file path
                2. type (eg, modelnew)
                3. domain id
                4-n. domain definitions (tab-delimited chain,start,end for each
                                         segment)

   Returns:     Nothing
   Output:      Calculated alignments are placed in alignment directory
                  as determined by modelid_2_domdir()

=cut

sub cut_domains {

   my $binaries = locate_binaries() ;
   my $specs = modtie::set_modtie_specs();

#MAKESURE that Taret domains are deposited in right dir structure - a la modbase first 2 character prefix / modelid / model.pdb

   my $subsets ;
   my $movethreshold = 500 ;
   my $origdir = getcwd() ;
   my $tempdir = tempdir(CLEANUP => 1) ;
   chdir $tempdir ;
   my @movethese ;
   my @movedest ;

   my $lastsig ;
   while (my $line = <STDIN>) {
      chomp $line ;
      my ($file_path, $type, $id, @defs) = split(/\t/, $line) ;
      my $cur_fn = $file_path;

      my ($chain, $start, $end) ;
      my $j = 0 ;
      my @defstring;
      my (@chain, @start, @end) ;
      my $actualstart ; #deals with modpipe 2.0 1-n model residue numbering
      if ($type =~ /^modelnew/) { ($actualstart) = ($type =~ /([0-9]+)/) ; }

      while ($j <= $#defs) {
         if ($type =~ /^modelnew/) {
            push @{$chain}, $defs[$j] ; $j++ ;

            my $t_range = $defs[$j] ;
            push @{$start}, $defs[$j] - $actualstart + 1; $j++ ;

            $t_range .= "-$defs[$j]" ;
            push @{$end}, $defs[$j] -$actualstart + 1; $j++ ;

            push @defstring, $t_range ;
         } else {
            push @{$chain}, $defs[$j] ; $j++ ;

            my $t_range = $defs[$j] ;
            push @{$start}, $defs[$j]; $j++ ;

            $t_range .= "-$defs[$j]" ;
            push @{$end}, $defs[$j] ; $j++ ;

            push @defstring, $t_range ;
         }
      }
      my $defstring = ' ' ;
      if ($type =~ /^model/) { $defstring = join('_', @defstring) ;}

      my $tfinal_dir = modelid_2_domdir({ modelid => $id, specs => $specs}) ;
      my $tfinal_name = $id.'.'.$defstring.'.pdb' ;
      if (-s $tfinal_dir.'/'.$tfinal_name) {
         print STDERR "note: $tfinal_name already cut\n" ;
         next;
      }

      my ($compress_fl, $unc_fh, $unc_fn) ;
      if ($file_path =~ /\.gz$/) {
         $compress_fl = 1;
         ($unc_fh, $unc_fn) = tempfile() ; close($unc_fh) ;
         system("$binaries->{zcat} $cur_fn > $unc_fn") ;
         $cur_fn = $unc_fn ;
      }


      my $altloc_fl = `$binaries->{altloc_check} < $cur_fn` ;
      chomp $altloc_fl ;

      if ($altloc_fl) {
         my ($t_fh, $t_fn) =
            tempfile("tempsource.XXXXXX", SUFFIX=>"pdb") ; close($t_fh) ;

         print STDERR "note: $file_path has an altloc set, now filtering: ." ;
         system("$binaries->{altloc_filter} $cur_fn > $t_fn") ;
         print STDERR "\bx\n" ;
         if ($compress_fl) {unlink $unc_fn;}
         $cur_fn = $t_fn ;
      }

      my $fileid ;
      if ($type eq 'subset') {
         $fileid = $id ;
      } else {
         $fileid = $id.'.'.$defstring ;
      }

      my ($cutpdb_fn, $cutpdb_fh) ;
      if ($type =~ /^model/) {
         $cutpdb_fn = "$fileid.pdb.gz" ;
      } else {
         my ($cutpdb_fh, $cutpdb_fn) =
# added gz here - modtie::pibase::subset_extract will automatically gzip it
         tempfile( "$fileid.XXXXXX", SUFFIX => ".pdb.gz" ) ;
         close($cutpdb_fh) ;
      }

      my $extract_errors = modtie::pibase::subset_extract({
            in_fn => $cur_fn,
            out_fn => $cutpdb_fn,
            chain => $chain,
            start => $start,
            end => $end}) ;


      if ($#{$extract_errors} >= 0 ) {
         foreach my $j ( 0 .. $#{$extract_errors}) {
            print STDERR "ERROR: $file_path, extract $id: ".
            "modtie::pibase::subset_extract(): $extract_errors->[$j]\n" ;
         }
      } elsif (-z $cutpdb_fn) {
         print STDERR "ERROR: $file_path, extract $id: ".
                      "empty subset extract pdb file\n" ;
         unlink $cutpdb_fn ;
      } else {
         my $deposit_dir = modelid_2_domdir({modelid => $id, specs => $specs});
         push @movethese, $cutpdb_fn ;
         push @movedest, $deposit_dir ;
         print "$file_path\t$id\t$defstring\t$deposit_dir/$cutpdb_fn\n" ;
      }

      if ( ($altloc_fl) && ($file_path ne $cur_fn) ) {
         unlink $cur_fn;
      }

      if ($#movethese == $movethreshold) {
         foreach my $j ( 0 .. $#movethese) {
            if (! -d $movedest[$j]) {
               File::Path::mkpath($movedest[$j]) ; }
            modtie::pibase::safe_move($movethese[$j], $movedest[$j]) ;}
         @movethese = () ;
         @movedest = () ;
      }

   }

   foreach my $j ( 0 .. $#movethese) {
      if (! -d $movedest[$j]) {
         File::Path::mkpath($movedest[$j]) ; }
      modtie::pibase::safe_move($movethese[$j], $movedest[$j]) ;}

   chdir $origdir ;

}


=head2 format_modbaseoutput()

   Title:       format_modbaseoutput()
   Function:    Convert MODTIE output to format for MODBASE table import
   Args:        $_->{project}
                $_->{binary_fn}
                $_->{complexes_fn}

   Returns:     Nothing
   Output:      Tab-delimited format:
                1. complex number
                2. subunit number
                3. seq_id
                4. model_id
                5. resrange
                6. bdp_id
                7. SCOP sid
                8. score_max
                9. total number of subunits
               10. project name
               11. comment
=cut

sub format_modbaseoutput {

   my $in = shift ;

   if (!exists $in->{project} ||
       !exists $in->{binary_fn} ||
       !exists $in->{complexes_fn}) {
      die "format_modbaseoutput() requires: ".
          "project, binary_fn, complexes_fn\n" ;
   }

   if (!exists $in->{out_fn}) { $in->{out_fn} = '-' ; }
   my $count = 1 ;
   open(OUTF, ">$in->{out_fn}") ;


   open(BINF, $in->{binary_fn}) ;
   while (my $line = <BINF>) {
      if ($line !~ /^INTERFACE/) {next;}
      chomp $line;

      my (undef, $tmpl_sid1, $tmpl_sid2,
          $seq1, $model_id_1,  $resrange1, $scopfa1,
          $seq2, $model_id_2,  $resrange2, $scopfa2,
          $num_tmplcon, $num_alncon,
          $pot_fn_base, $pot_type, $pot_details, $pot_r,
          $rawscore, $zscore, $zprime, $z_2, $avgscore, $minscore, $minscore_tn,
          $stdevscore, $zmin, $zmin_tn, $falsepos,
          $rmsd1, $equivpos1, $numres1a, $numres1b, $numident1,
          $bs_numres1_tmpl, $bs_numres1_targ, $bs_numident1,
          $rmsd2, $equivpos2, $numres2a, $numres2b, $numident2,
          $bs_numres2_tmpl, $bs_numres2_targ, $bs_numident2) =
            split(/\t/, $line) ;

      my $seqsig = join("\n", sort ($seq1, $seq2)) ;

      my ($bdp_id) = ($tmpl_sid1 =~ /BDP([0-9]+)\-/);

      my $scopsid1 = $tmpl_sid1 ; $scopsid1 =~ s/BDP[0-9]+\-.*SCOP\.// ;
      my $scopsid2 = $tmpl_sid2 ; $scopsid2 =~ s/BDP[0-9]+\-.*SCOP\.// ;

      my @outvals1 = ($count, 1, $seq1, $model_id_1, $resrange1, $bdp_id,
            $scopsid1, $zscore,2, $in->{project});
      my @outvals2 = ( $count, 2, $seq2, $model_id_2, $resrange2, $bdp_id,
            $scopsid2, $zscore,2, $in->{project}) ;

      push @outvals1, join(' ', sort keys %{$in->{comments}->{binary}->{$seqsig}}) ;
      push @outvals2, join(' ', sort keys %{$in->{comments}->{binary}->{$seqsig}}) ;
      print OUTF join("\t", @outvals1)."\n" ;
      print OUTF join("\t", @outvals2)."\n" ;
      $count++ ;

   }
   close(BINF) ;

   my $lastcid ='';
   my $outlines = [] ;
   open(COMPLEXF, $in->{complexes_fn}) ;
   while (my $line = <COMPLEXF>) {
      if ($line =~ /^\#/) {next;}
      chomp $line;

      my ($cid, $subno, $seq_id, $model_id, $resrange, $bdp_id,
          $tmpl_sid, $score_avg, $score_max) = split(/\t/, $line) ;
      my $scopsid = $tmpl_sid ; $scopsid =~ s/BDP[0-9]+\-.*SCOP\.// ;

      if ($cid ne $lastcid && $lastcid ne '') {
         my $numsub = $#{$outlines} + 1 ;

         my $curcomment = join(' ',
               sort keys %{$in->{comments}->{complexes}->{$lastcid}}) ;

         foreach my $j ( 0 .. $#{$outlines}) {
               my @outvals ;
               push @outvals, @{$outlines->[$j]} ;
               push @outvals, ($numsub, $in->{project}) ;
               push @outvals, $curcomment ;
               print OUTF join("\t", @outvals)."\n" ;
         }
         $count++ ;
         $outlines = [] ;
      }

      push @{$outlines}, [ $count, $subno, $seq_id, $model_id,
                           $resrange, $bdp_id, $scopsid, $score_max ] ;
      $lastcid = $cid ;
   }
   close(COMPLEXF) ;

   {
      my $numsub = $#{$outlines} + 1 ;
      my $curcomment = join(' ',
         sort keys %{$in->{comments}->{complexes}->{$lastcid}}) ;

      foreach my $j ( 0 .. $#{$outlines}) {
         my @outvals ;
         push @outvals, @{$outlines->[$j]} ;
         push @outvals, ($numsub, $in->{project}) ;
         push @outvals, $curcomment ;
         print OUTF join("\t", @outvals)."\n" ;
      }
   }

   close(OUTF) ;

}


=head2 databaseid_2_seqid()

   Title:       databaseid_2_seqid()
   Function:    Converts database_id to MODBASE seq_id sequence identifiers via
                 MODBASE.synonyms DBI access
   Args:        Nothing
   Returns:     Nothing

   Input:       database_id

   Output:      Tab-delimited:
                1. database_id
                2. number of MODBASE hits
                3. MODBASE seq_id

=cut

sub databaseid_2_seqid {

   my $dbh ;
   ($dbh->{mb}) = modtie::connect_modbase() ;
   while (my $line = <STDIN>) {
      chomp $line;
      my ($seqid) = modtie::pibase::mysql_fetchcols($dbh->{mb},
         "SELECT distinct seq_id FROM synonyms where database_id = \"$line\"") ;
      if ($#{$seqid} >= 0) {
         my $numhits = $#{$seqid} + 1;
         foreach my $tseqid (@{$seqid}) {
            my @outvals = ($line, $numhits, $tseqid) ;
            print join("\t", @outvals)."\n" ;
         }
      } else {
         my @outvals = ($line, 0 ) ;
         print join("\t", @outvals)."\n" ;
      }
   }
   $dbh->{mb}->disconnect ;

}


=head2 parse_modbase_ali_xml()

   Title:       parse_modbase_ali_xml()
   Function:    Parse MODBASE alignment XML file for meta-info and actual aln
   Args:        $_->{fn} = alignment XML file name
   Returns:     $_->[i]->{field} = val; meta info about ith alignment

   Output:      MODBASE alignment files to :
                dest_dir/substr(seq_id,0,3)/seq_id/alignments/align_id.ali

=cut

sub parse_modbase_ali_xml {

   my $in = shift ;
   my $info ;
   my $in_content = 0 ;
   my $dest_dir = $in->{dest_dir} ;

   open(ALIF, $in->{fn}) ;
   my $aln; my $out_alifh ;
   while (my $line = <ALIF>) {
      chomp $line;
      if ($line =~ /\</) {                              # Non-alignment line
         if ($line =~ /^\<\?xml/ || $line =~ /\<files\>/
             || $line =~ /\<\/files>/) {
            next;
         } elsif ($line =~ /\<content/) {
            $in_content = 1;
         } elsif ($line =~ /\<alignmentfile/) {
            push @{$aln}, {} ;
            # ALSO OPEN FILEHANDLE TO ALI OUTPUT
         } elsif ($line =~ /\<\/alignmentfile/) {       # Store aln on disk
            # CLOSE CUR ALI FILEHANDLE
            close($out_alifh) ;
         } elsif ($line =~ /\<\/content/) {
            $in_content = 0;
         } else {                                       # Parse out meta fields
            my ($key, $val, $key2) =
               ($line =~ /\<([^\>]+)\>([^\>]*)\<\/([^\>]+)\>/) ;
            if ($key ne $key2) {
               print STDERR "ERROR: ALI XML parsing error $key != $key2\n".
                            "line: $line\n"; }
            if (defined $val) { #eg ignoring: <database_id></database_id>
               $val =~ s/^\s+// ; $val =~ s/\s+$// ;
               $aln->[$#{$aln}]->{$key} = $val ; }

            if ($key eq 'align_id') {                   # Open ALI filehandle
               my $out_dir = $in->{dest_dir}.'/'.
                  substr($in->{seq_id},0,3).'/'.$in->{seq_id}.'/alignments' ;
               if (!-d $out_dir) {mkpath($out_dir);}
               my $out_alifn = $out_dir.'/'.$val.'.ali' ;
               open($out_alifh, ">".$out_alifn) ;
               $aln->[$#{$aln}]->{ali_fn} = $out_alifn ;
            }
         }
      } elsif ($in_content) {                   # Store actual alignment lines
         print {$out_alifh} $line."\n" ;
      }
   }
   close(ALIF) ;

   return $aln ;
}



=head2 parse_modbase_pdb_xml()

   Title:       parse_modbase_pdb_xml()
   Function:    Parse MODBASE PDB XML file for meta-info and model PDB files
   Args:        $_->{fn} = PDB XML file name
   Returns:     $_->[i]->{meta_info}->{field} = val; info about ith model
   Output:      MODBASE model PDB files to :
                dest_dir/substr(seq_id,0,3)/seq_id/models/model_id.pdb.gz

=cut

sub parse_modbase_pdb_xml {

   my $in = shift ;
   my $info ;
   my $in_content = 0 ;
   my $dest_dir = $in->{dest_dir} ;

   open(PDBF, $in->{fn}) ;
   my $models; my $out_pdbfh ;
   while (my $line = <PDBF>) {
      chomp $line;
      if ($line =~ /\</) {             # Non-alignment line
         if ($line =~ /^\<\?xml/ || $line =~ /\<files\>/
             || $line =~ /\<\/files>/) {
            next;
         } elsif ($line =~ /\<content/) {
            $in_content = 1;
         } elsif ($line =~ /\<pdbfile/) {
            push @{$models}, {} ;
         } elsif ($line =~ /\<\/pdbfile/) {
            close($out_pdbfh) ;
            system("gzip -f ".$models->[$#{$models}]->{pdb_fn}) ;
         } elsif ($line =~ /\<\/content/) {
            $in_content = 0;
         } else {                      # Parse out meta fields
            my ($key, $val, $key2) =
               ($line =~ /\<([^\>]+)\>([^\>]+)\<\/([^\>]+)\>/) ;
            $val =~ s/^\s+// ; $val =~ s/\s+$// ;
            $models->[$#{$models}]->{$key} = $val ;
            if ($key eq 'model_id') {
               my $out_dir = $in->{dest_dir}.'/'.
                  substr($in->{seq_id},0,3).'/'.$in->{seq_id}.'/models' ;
               if (!-d $out_dir) {mkpath($out_dir);}
               my $out_pdbfn = $out_dir.'/'.$val.'.pdb' ;
               open($out_pdbfh, ">".$out_pdbfn) ;
               $models->[$#{$models}]->{pdb_fn} = $out_pdbfn ;
            }
         }
      } elsif ($in_content) {                           # Print PDB lines
         print {$out_pdbfh} $line."\n";
         if ($line =~ /^REMARK 220 .*\:/) { # store REMARK lines as meta info
            $line =~ s/^REMARK 220 // ;
            my ($key, $val) = ($line =~/(.+)\:(.+)/)  ;
            $val =~ s/^\s+// ; $val =~ s/\s+$// ;
            $models->[$#{$models}]->{$key} = $val ;
         }
      }
   }
   close(PDBF) ;

   return $models ;
}


=head2 write_model_list()

   Title:       write_model_list()
   Function:    Writes out model details for actual MODTIE run
   Args:        $_->{out_fn} = output filename
                $_->{model_entries}->{field} = [] ; array of values for
                                                    each model
   Returns:     Nothing
   Output:      to $_->{out_fn}: tab_delimited fields specified in:
                $specs->{file_format}->{model_list} ;

=cut

sub write_model_list {

   my $in = shift ;
   my $model_entries = $in->{model_entries} ;
   my $out_fh ;
   if (exists $in->{out_fh})  { $out_fh = $in->{out_fh} ; }
   else                       { open($out_fh, ">".$in->{out_fn}) ; }

   my $specs ;
   if (exists $in->{specs}) { $specs = $in->{specs} ; }
   else { $specs = set_modtie_specs() ; }

   my @out_fields = @{$specs->{file_format}->{model_list}} ;

   #pull off last 2: baseali_dir, basemodel_dir
   pop @out_fields; pop @out_fields ;

   foreach my $j ( 0 .. $#{$model_entries->{model_id}}) {
      my @outvals ;
      map {push @outvals, $model_entries->{$_}->[$j]} @out_fields ;
      push @outvals, $in->{ali_dir}; push @outvals, $in->{models_dir};
      print {$out_fh} join("\t", @outvals)."\n" ;
   }
   if (!exists $in->{out_fh}) {close($out_fh) ; }

   return ;
}


=head2 read_model_list()

   Title:       read_model_list()
   Function:    Reads model details for actual MODTIE run
   Args:        $_->{fn} = model list filename
   Returns:     $_->{field} = [] ; array of field value for each model
   Output:      Nothing

=cut

sub read_model_list {

   my $in = shift ;
   my $fn = $in->{fn} ;

   my $specs ;
   if (exists $in->{specs}) { $specs = $in->{specs} ; }
   else { $specs = set_modtie_specs() ; }

   my $model_entries = {} ; my $pdb2models = {} ;
   open(LISTF, $fn) ;
   while (my $line = <LISTF>) {
      chomp $line;
      if ($line =~ /^#/) {next;}
      my @t = split(/\t/, $line) ;

      map {push @{$model_entries->{$specs->{file_format}->{model_list}->[$_]}},
           $t[$_]} (0 .. $#{$specs->{file_format}->{model_list}});

      push @{$pdb2models->{$model_entries->{pdb_code}->[$#{$model_entries->{pdb_code}}]}}, $#{$model_entries->{pdb_code}} ;
   }
   close(LISTF) ;

   return {
      model_entries => $model_entries,
      pdb2models => $pdb2models
   };
}


=head2 write_seqid_setinfo()

   Title:       write_seqid_setinfo()
   Function:    Writes out sequence set membership for MODTIE run
   Args:        $_->{out_fn} = output filename
                $_->{seqset} = $seqset data structure:
                            ->{union}->{seq_id} = sum of set labels (1,2,or 3)
                            ->{1|2}->{seq_id} = hash list of seq_id in each set

   Returns:     Nothing
   Output:      to $_->{out_fn}: tab_delimited
                1. seq_id
                2. sum of set membership labels

=cut

sub write_seqid_setinfo {

   my $in = shift ;
   my $seqset = $in->{seqset} ;

   my $out_fh ;
   if (exists $in->{out_fh})  { $out_fh = $in->{out_fh} ; }
   else                       { open($out_fh, ">".$in->{out_fn}) ; }

   foreach my $seq_id ( sort keys %{$seqset->{union}}) {
      my @outvals = ($seq_id, $seqset->{union}->{$seq_id}) ;
      print {$out_fh} join("\t", @outvals)."\n" ;
   }
   if (!exists $in->{out_fh}) { close($out_fh) ; }

   return ;

}


=head2 read_seqid_setinfo()

   Title:       read_seqid_setinfo()
   Function:    Reads sequence set membership for MODTIE run
   Args:        $_->{fn} = filename

   Returns:     $_->{seqset} = $seqset data structure:
                            ->{union}->{seq_id} = sum of set labels (1,2,or 3)
                            ->{1|2}->{seq_id} = hash list of seq_id in each set

=cut

sub read_seqid_setinfo {

   my $in = shift ;
   my $fn = $in->{fn} ;

   my $seqset = {} ;
   if (! -s $fn) {
      die "FATAL ERROR: sequence file $fn not found\n" ; }

   open(SETF, $fn) ;
   while (my $line = <SETF>) {
      chomp $line;
      if ($line =~ /^#/) {next;}
      my ($seq_id, $setsum) = split(/\t/, $line) ;
      $seqset->{union}->{$seq_id} = $setsum ;
      if ($setsum < 3) {
         $seqset->{$setsum}->{$seq_id}++ ;
      } else {
         $seqset->{1}->{$seq_id}++ ;
         $seqset->{2}->{$seq_id}++ ;
      }
   }
   close(SETF) ;

   return $seqset;

}


=head2 extract_required_pibase_datafiles()

   Title:       extract_required_pibase_datafiles
   Function:    Gets necessary data files from a full PIBASE installation
                 to package for the MODTIE distribution.
   Input:       $_ = specs hashref
   Return:      specs - hashref

=cut

sub extract_required_pibase_datafiles {

   require DBI ;

   my $in = shift ;
   my $specs ;
   if (exists $in->{specs}) { $specs = $in->{specs} ; }
   else { $specs = set_modtie_specs() ; }

# Setup destination directories; create if necessary
# subsets_residues and bdp_residues metatod tars are just
#    manually copied from the PIBASE directories
#      subsets_residues => "subsets_residues",
#      bdp_residues => "bdp_residues",
   my $goal_dir = {
      subsets_files => "subsets_files",
      tod => "tod",
   } ;
   foreach my $type (keys %{$goal_dir}) {
      if (!-d $goal_dir) {mkpath($goal_dir->{$type});} }

   my $required_tod = {
      bdp_files => 1,
      bdp_residues_tables => 1,
      intersubset_contacts => 1,
      subsets => 1,
      subsets_files => 1,
      subsets_residues_tables => 1,
      subsets_sequence => 1
   } ;

# Set template interface list
   my $templates_fn ;
   if (!exists $in->{templates_fn}) {
      $in->{templates_fn} = $specs->{templates_fn}->{default}; }

# Read in the template interface list; get list of bdp_id and subset_id
   my $pbtmpl = readin_tmpllist({fn => $in->{templates_fn}});
   my $bdp_list; my $sid_list;
   foreach my $j ( 0 .. $#{$pbtmpl->{sid1}}) {
      $bdp_list->{$pbtmpl->{bdp_id}->[$j]}++ ;
      $sid_list->{$pbtmpl->{sid1}->[$j]}++ ;
      $sid_list->{$pbtmpl->{sid2}->[$j]}++ ;
   }

# connect to PIBASE and get all file locations
   my $dbh ;
   ($dbh->{pb}) = modtie::pibase::connect_pibase() ;
   my $subsets_residues_tables = modtie::pibase::mysql_hashload($dbh->{pb},
      "SELECT bdp_id, source_file FROM subsets_residues_tables ") ;

   my $bdp_residues_tables = modtie::pibase::mysql_hashload($dbh->{pb},
      "SELECT bdp_id, source_file FROM bdp_residues_tables ") ;

   my $subsets_files = modtie::pibase::mysql_hashload($dbh->{pb},
      "SELECT subset_id, file_path FROM subsets_files ".
      'WHERE subset_id LIKE "%SCOP%"') ;

# Copy necessary files into a local directory------------
# 1. bdp_residues
   if (0) {
   print STDERR "Getting bdp_residues tod\n";
   my $copy_these; my $copy_dest ;
   foreach my $bdp (sort {$a <=> $b} keys %{$bdp_list}) {
      my $orig_fn = $bdp_residues_tables->{$bdp} ;
      my $t_1 = $specs->{pibase_specs}->{metatod_dir}."/bdp_residues" ;
      my $t_2 = $goal_dir->{bdp_residues} ;

      my $new_fn = $orig_fn ;
      $new_fn =~ s/$t_1/$t_2/ ;
      if (-s $new_fn) {next;}

      my $dirname = File::Basename::dirname($new_fn) ;
      if (!-d $dirname) {mkpath($dirname)} ;
      push @{$copy_these}, $orig_fn ;
      push @{$copy_dest}, $new_fn ;
   }

   print STDERR "* copying: " ;
   foreach my $j ( 0 .. $#{$copy_these}) {
      print "Copying ".$copy_these->[$j]." to ".$copy_dest->[$j]."\n" ;
      modtie::pibase::safe_copy($copy_these->[$j], $copy_dest->[$j]) ; }
   print STDERR "X\n" ;
   }

# 2. subsets_residues
# Dump these into a local directory; while respecting dir structure.
   if (0) {
   print STDERR "Getting subsets_residues tod\n";
   my $copy_these; my $copy_dest ;
   foreach my $bdp (sort {$a <=> $b} keys %{$bdp_list}) {
      my $orig_fn = $subsets_residues_tables->{$bdp} ;
      my $t_1 = $specs->{pibase_specs}->{metatod_dir}."/subsets_residues" ;
      my $t_2 = $goal_dir->{subsets_residues} ;

      my $new_fn = $orig_fn ;
      $new_fn =~ s/$t_1/$t_2/ ;
      if (-s $new_fn) {next;}

      my $dirname = File::Basename::dirname($new_fn) ;
      if (!-d $dirname) {mkpath($dirname)} ;
      push @{$copy_these}, $orig_fn ;
      push @{$copy_dest}, $new_fn ;
   }

   print STDERR "* copying: " ;
   foreach my $j ( 0 .. $#{$copy_these}) {
      print "Copying ".$copy_these->[$j]." to ".$copy_dest->[$j]."\n" ;
      modtie::pibase::safe_copy($copy_these->[$j], $copy_dest->[$j]) ; }
   print STDERR "X\n" ;
   }


# 3. subsets_files
# Dump these into a local directory; while respecting dir structure.
   {
   print STDERR "Getting subsets_files tod\n";
   my $copy_these; my $copy_dest ;
   foreach my $sid (sort keys %{$sid_list}) {
      my $orig_fn = $subsets_files->{$sid} ;
      my $t_1 = $specs->{pibase_specs}->{old_subsets_dir} ;
      my $t_2 = $goal_dir->{subsets_files} ;

      my $new_fn = $orig_fn ;
      $new_fn =~ s/$t_1/$t_2/ ;
      if (-s $new_fn) {next;}

      my $dirname = File::Basename::dirname($new_fn) ;
      if (!-d $dirname) {mkpath($dirname)} ;
      push @{$copy_these}, $orig_fn ;
      push @{$copy_dest}, $new_fn ;
   }

   print STDERR "* copying: " ;
   foreach my $j ( 0 .. $#{$copy_these}) {
      print "Copying ".$copy_these->[$j]." to ".$copy_dest->[$j]."\n" ;
      modtie::pibase::safe_copy($copy_these->[$j], $copy_dest->[$j]) ; }
   print STDERR "X\n" ;
   }

# Setup TOD copies
   {
   print STDERR "Getting PIBASE tod\n";
   my $copy_these; my $copy_dest ;
   foreach my $table  (sort keys %{$required_tod}) {
      my $orig_fn = $specs->{pibase_specs}->{old_tod_dir}."/".$table ;
      my $t_1 = $specs->{pibase_specs}->{old_tod_dir} ;
      my $t_2 = $goal_dir->{tod} ;

      my $new_fn = $orig_fn ;
      $new_fn =~ s/$t_1/$t_2/ ;
      if (-s $new_fn) {next;}

      my $dirname = File::Basename::dirname($new_fn) ;
      if (!-d $dirname) {mkpath($dirname)} ;
      push @{$copy_these}, $orig_fn ;
      push @{$copy_dest}, $new_fn ;
   }

   print STDERR "* copying: " ;
   foreach my $j ( 0 .. $#{$copy_these}) {
      print "Copying ".$copy_these->[$j]." to ".$copy_dest->[$j]."\n" ;
      modtie::pibase::safe_copy($copy_these->[$j], $copy_dest->[$j]) ; }
   print STDERR "X\n" ;
   }

}


=head2 replace_basedir_inplace()

   Title:       replace_basedir_inplace()
   Function:    Changes file base directory in a key->filepath hash
   Input:       $_->{hash}->{key} = filepath
                $_->{old_dir} = old directory name
                $_->{new_dir} = new directory name
   Return:      specs - hashref

=cut

sub replace_basedir_inplace {

   my $in = shift ;
   my $hash = $in->{hash} ;
   my $old_dir = $in->{old_dir} ;
   my $new_dir = $in->{new_dir} ;

   foreach my $key (keys %{$hash}) {
      my $new_fn = $hash->{$key};
      $new_fn =~ s/$old_dir/$new_dir/ ;
      $hash->{$key} = $new_fn ;
   }
}

1 ;
