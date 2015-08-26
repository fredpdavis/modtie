=head1 NAME

complexes.pm - routines for structure-based prediction of protein complexes

=head1 DESCRIPTION

The complexes.pm perl library contains subroutines to predict higher
order protein complexes based on template complexes and structure-based
predictions of binary interactions.

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


package modtie::complexes ;
use strict;
use warnings;
use Exporter ;
our @ISA = qw/Exporter/ ;
our @EXPORT_OK = qw// ;

use modtie ;
use modtie::pibase ;
use modtie::potentials qw/maxval_arrref/ ;
use File::Temp qw/tempfile tempdir/;
use Digest::MD5 qw/md5/;
use Data::Dumper qw/Dumper/ ;


=head2 list_complexes()

   Title:       list_complexes()
   Function:    Determine possible higher order complexes given list of
                  binary interactions.

   Args:        $_->{in_fn}->{interactions_fn} - binary interactions file
   Returns:     Nothing
   Output:      

=cut

sub list_complexes {

   my $params = shift;

   my $specs ;
   if (exists $params->{specs}) { $specs = $params->{specs} ; }
                           else { $specs = modtie::set_modtie_specs() ; }

# Set file locations and thresholds if specified
   if (!exists $params->{in_fn}->{interactions_fn}) {
      print STDERR "ERROR: binary interaction file not specified\n" ;
      return ;}

   if (! -s $params->{in_fn}->{interactions_fn}) {
      print STDERR "ERROR: binary interaction file not found\n" ;
      return ;}

   if (!exists $params->{in_fn}->{seqid_domainarch_fn}) {
      print STDERR "ERROR: seqid_domainarch_fn not specified\n" ;
      return ;}

   if (!exists $params->{z_thresh}) {
      $params->{z_thresh} = $specs->{complex_z_thresh} ; }

   if (!exists $params->{complex_aln_tmpl_contacts_thresh}) {
      $params->{complex_aln_tmpl_contacts_thresh} =
         $specs->{complex_aln_tmpl_contacts_thresh} ;}

#   print "#INFO: z_score threshold: $params->{z_thresh}\n" ;

# Preload PIBASE Data
   print STDERR "Loading PIBASE data: " ;
   my $pb_pl = list_complexes_preload() ;
   print STDERR "X\n" ;

   my $sid2class = $pb_pl->{sid2class} ;
   my $sid2bdp = $pb_pl->{sid2bdp} ;
   my $sid12chains = $pb_pl->{sid12chains} ;
   my $bdp2sid12 = $pb_pl->{bdp2sid12} ;

# Load PIBASE interfaace cluster assignments
   my $iclusts = load_interface_clusters({
      specs => $specs }) ;

# Load sequence domain architecture strings
   my $domdata ; $domdata->{seqinfo} = {} ;
   $domdata->{seqinfo}->{arch} = modtie::read_seqid_domainarch({
      fn => $params->{in_fn}->{seqid_domainarch_fn}});

   my $mbinfo      = modtie::_parse_domarchs($domdata) ;
   my $scop2seq    = $mbinfo->{scop2seq} ;   #->{$scop_fa}->{$seq_id}++
                                             #->{$scop_sf}->{$seq_id}++
   my $seq2dom     = $mbinfo->{seq2dom} ;    #->{$seq_id}->[i]->{def} =
                                             # ith domain definition in seq_id
   my $seq2domfrag = $mbinfo->{seq2domfrag}; #->{$seq_id}->[i]->[j]->{start|end}
                                             # jth fragment of ith seq_id domain
   my $scop2num    = $mbinfo->{scop2num} ;   #->{$scop_fa} = #seq domains in fam

   my $headers = modtie::scoring_main({get_assess_headers => 1}) ;
   my $f; foreach my $j (0 .. $#{$headers}) { $f->{$headers->[$j]} = $j; }

# Load binary interactions
# Setup data structure to hold information on complexes:
   #->{summary}->{cid} = {} , ->{details}->{subunits}->{..}
   #->{summary} used to be $complexcompo
   #->{counter} used to be $complexcounter ;
   my $complexes = { counter => 0,
                     pb_pl => $pb_pl };
   my ($binarypred, $z_scores) ;

   print STDERR "Loading binary interactions: " ;
   open (BINTERACTIONS, $params->{in_fn}->{interactions_fn}) ;
   while (my $line = <BINTERACTIONS>) {
      if ($line !~ /^INTERFACE/ || $line =~ /ERROR\tstdev of 0/) {next;}

      chomp $line;
      my @fields = split(/\t/, $line) ;

# if the interaction doesnt pass z_score or fraction of template contacts
# threshold, skip it
      if (!defined $fields[$f->{z}] || $fields[$f->{z}] eq '' ||
          $fields[$f->{z}] > $params->{z_thresh} || 
          ($fields[$f->{aln_contacts}] / $fields[$f->{tmpl_contacts}]) <
            $params->{complex_aln_tmpl_contacts_thresh}) {
         next; }

      my $sid1 = $fields[$f->{subset_id_1}] ;
      my $sid2 = $fields[$f->{subset_id_2}] ;
      my $sid12= $sid1."\n".$sid2 ;
      my $bdp = $sid2bdp->{$sid1} ;


      my $sid12sig = $sid12 ;
      if ($sid2 lt $sid1) {$sid12sig = $sid2."\n".$sid1;}

# If the template interface is intra-chain, only keep the target interaction
#  if it involves different parts of the same sequence.
      if ($sid12chains->{$sid12sig} eq 'same' &&
          !($fields[$f->{seq1}] eq $fields[$f->{seq2}] &&
            $fields[$f->{t_def1}] ne $fields[$f->{t_def2}])) {next;}

      my $hit1 = $fields[$f->{seq1}]."^".$fields[$f->{model_id_1}].
                 "^".$fields[$f->{t_def1}] ;
      my $hit2 = $fields[$f->{seq2}]."^".$fields[$f->{model_id_2}].
                 "^".$fields[$f->{t_def2}] ;
      my $hit12= $hit1."\n".$hit2 ;

# tally up the individual domain (sid1, sid2) as well as interaction (sid12) that the entry describes
      $binarypred->{$bdp}->{sid}->{$sid1}->{$hit1}++ ;
      $binarypred->{$bdp}->{sid}->{$sid2}->{$hit2}++ ;
      $binarypred->{$bdp}->{sid12}->{$sid12}->{$hit12}++ ;

# iterate over all interfaces in the same cluster as the current template
# and populate those hits as well
# - those interfaces receive same z-score as the actual scored interface
      my $curiclus = $iclusts->{sidpair2clus}->{$sid1."\n".$sid2} ;
      foreach my $newsid12 (keys %{$iclusts->{fampairclus}->{$curiclus}}) {
         my ($newsid1, $newsid2) = split(/\n/, $newsid12) ;
         my $inherit_bdp = $sid2bdp->{$newsid1} ;
         $binarypred->{$inherit_bdp}->{sid}->{$newsid1}->{$hit1}++ ;
         $binarypred->{$inherit_bdp}->{sid}->{$newsid2}->{$hit2}++ ;
         $binarypred->{$inherit_bdp}->{sid12}->{$newsid12}->{$hit12}++ ;
         $z_scores->{$newsid12}->{$hit12} = $fields[$f->{z}];
      }
      $z_scores->{$sid12}->{$hit12} = $fields[$f->{z}] ;
      $complexes->{counter}++;

# Setup 'complex' structure for this interaction
# Going through the riggamarol to save on memory
# - sort template domains based on alphabetical order of hits
      my @sids = ($sid1, $sid2) ;
      my @hits = ($hit1, $hit2) ;

# - set topology string
      my $curtopostring_hits = '0-1' ;
      my $curtopostring_sids = '0-1' ;

# ORIGINAL DATA STRUCTURE:
# ->{compo} = join(" ", (sort @hits));
# ->{tmpls} = join(" ", @sids[(sort {$hits[$a] cmp $hits[$b]} (0 .. $#hits) )]);
# ->{topo} = join(',', sort @hits).' '.('1-1'|'1-2');
      $complexes->{summary}->{$complexes->{counter}} = {
         numunits => 2,
         score_avg => $fields[$f->{z}],
         score_max => $fields[$f->{z}],
         topo_hits => $curtopostring_hits,
         topo_sids => $curtopostring_sids,
      } ;

# Specify complex details
      {
         my $cur_complex_summary =
            $complexes->{summary}->{$complexes->{counter}} ;
         push @{$cur_complex_summary->{ocid}},
            $complexes->{counter} ;

         my $hash_these = [
            ['hits', $hit1],
            ['hits', $hit2],
            ['tmpls', $sid1],
            ['tmpls', $sid2],
         ] ;

         foreach my $item (@{$hash_these}) {
            my $type = $item->[0] ;
            my $val = $item->[1] ;

            if (!exists $complexes->{details}->{$type}->{hash}->{$val}) {
               push @{$complexes->{details}->{$type}->{list}}, $val;
                      $complexes->{details}->{$type}->{hash}->{$val} =
                     $#{$complexes->{details}->{$type}->{list}} ; }
         }

         push @{$cur_complex_summary->{hits}},
            $complexes->{details}->{hits}->{hash}->{$hit1} ;
         push @{$cur_complex_summary->{hits}},
            $complexes->{details}->{hits}->{hash}->{$hit2} ;
         push @{$cur_complex_summary->{tmpls}},
            $complexes->{details}->{tmpls}->{hash}->{$sid1} ;
         push @{$cur_complex_summary->{tmpls}},
            $complexes->{details}->{tmpls}->{hash}->{$sid2} ;
      }

      foreach my $j ( 0 .. 1) {
         my ($seq_id, $model_id, $resrange) = split(/\^/, $hits[$j]) ;
         my $hit_ind = $complexes->{details}->{hits}->{hash}->{$hits[$j]} ;
         if (!exists $complexes->{details}->{subunits}->{$hit_ind}) {
            $complexes->{details}->{subunits}->{$hit_ind} = {
               seq_id => $seq_id,
               model_id => $model_id,
               resrange => $resrange,
            } ;
         }

         push @{$complexes->{summary}->{$complexes->{counter}}->{subunits}},
            $hit_ind ;
      }
   }
   close(BINTERACTIONS) ;
   my $num_binary = $complexes->{counter} ;
   print STDERR " $num_binary\n" ;

# Build higher-order complexes out of the binary interactions using the
#  observed topologies of template complexes
   print STDERR "Building higher-order complexes: " ;
   foreach my $bdp (keys %{$binarypred}) {
      my $numints_hit = keys %{$binarypred->{$bdp}->{sid12}} ;
      if ($numints_hit < 2) {next;}

      my $t_composition ;

      my $init_content = [] ;
      my @int_list = keys %{$binarypred->{$bdp}->{sid12}} ;
      my $int_list = \@int_list;

      my $hit_list ;
      foreach my $j (0 .. $#int_list) {
        foreach my $hit (keys %{$binarypred->{$bdp}->{sid12}->{$int_list[$j]}}){
            push @{$hit_list->{$int_list[$j]}}, $hit ; }}

# Recurse over all possible interactions to build an initial list of complexes
      recurse_composition({
         complexes => $complexes,
         tmpltopo => {edges => $bdp2sid12->{$bdp},
                      chains => $sid12chains},
         bdp_id => $bdp,
         content => $init_content,
         intlist => $int_list,
         zscores => $z_scores,
         hitlist => $hit_list
      }) ;
   }
   print STDERR " ".($complexes->{counter} - $num_binary)."\n" ;

# Merger Logic:
# * all combos have already been generated, just iterate over with 2
#   loops - inner and outer - and join if a single protein bridges two of them.
# * then iterate over this until no more mergings are done, similar to a
#   hierarchical clustering algorithm
# * just make sure of the covalent checks are already being performed

# Print out raw complexes
   {
      my ($cur_fh, $cur_fn) = File::Temp::tempfile(
         "complexlist.round0.raw.$$.XXXXX", SUFFIX => '.modtie');
      display_complexes({ complexes => $complexes, out_fh => $cur_fh });
      close($cur_fh) ;
   }

# Remove redundant complexes
   print STDERR "Removing redundant complexes: " ;
   removecomplexredundancy({ complexes => $complexes }) ;
   {
      my ($cur_fh, $cur_fn) = File::Temp::tempfile(
         "complexlist.round0.nr.$$.XXXXX", SUFFIX => '.modtie') ;
      display_complexes({ complexes => $complexes, out_fh => $cur_fh });
      close($cur_fh) ;
   }
   my $num_nr = keys %{$complexes->{summary}} ;
   print STDERR " $num_nr complexes remaining\n" ;

# If specified, one round of raw complex mergings by bridging sequences
   if (exists $params->{covalent_merge_fl} ||
       (!exists $params->{covalent_merge_fl} &&
        $specs->{covalent_merge_fl} == 1)) {
      covalent_complex_merge({complexes => $complexes,
                              seq2dom   => $seq2dom});
   }

}


=head2 removecomplexredundancy()

   Title:       removecomplexredundancy()
   Function:    (in-place) removes topologically redundant complexes
   Args:        $_->{complexes} = complexes data structure

=cut

sub removecomplexredundancy {

   my $in = shift;
   my $complexes = $in->{complexes} ;
   my @deletecid ;

# Iterate over complexes, keeping track of what toplogies have been seen.
   foreach my $cid (keys %{$complexes->{summary}}) {
      my $keep_it = 0 ;
      my $fulltopo = join(',', sort @{$complexes->{details}->{hits}->{list}}[@{$complexes->{summary}->{$cid}->{hits}}]).' '.$complexes->{summary}->{$cid}->{topo_hits};

      my $md5topo = md5($fulltopo) ;

# keep current complex if 
# 1. topology not seen before OR
# 2. score_max is lower than previously seen complex OR
# 3. score_max is same, but score_avg is lower than previously seen complex

      if (exists $complexes->{md5topo}->{$md5topo} &&
          $complexes->{md5topo}->{$md5topo} != $cid ) {
         my $oldcid = $complexes->{md5topo}->{$md5topo} ;
         if ( ($complexes->{summary}->{$oldcid}->{score_max} >
               $complexes->{summary}->{$cid}->{score_max}) ||
              (($complexes->{summary}->{$oldcid}->{score_max} ==
                $complexes->{summary}->{$cid}->{score_max}) &&
               ($complexes->{summary}->{$oldcid}->{score_avg} >
                $complexes->{summary}->{$cid}->{score_avg}))) {
            $keep_it = 1 ;
         } else {
            push @deletecid, $cid ;
         }
      } else {
         $keep_it = 1 ;
      }

      if ($keep_it) {
         if (exists $complexes->{md5topo}->{$md5topo}) {
            push @deletecid, $complexes->{md5topo}->{$md5topo};}

         $complexes->{md5topo}->{$md5topo} = $cid ;
      }
   }

# Delete redundant complexes
   delete @{$complexes->{summary}}{@deletecid} ;

}


=head2 display_complexes()

   Title:       display_complexes()
   Function:    Display complexes
   Args:        $_->{complexes} = complexes data structure
                $_->{out_fh} = output file handle
                $_->{cid} = optional; if specified, only displays that complex

=cut

sub display_complexes {

   my $in = shift ;
   my $complexes = $in->{complexes} ;

   my @cids ;
   if (exists $in->{cid}) { @cids = $in->{cid} ; }
   else { @cids = sort {$a <=> $b} keys %{$complexes->{summary}}; }

   print {$in->{out_fh}} '#'.join("\t", 'compl', 'cid', 'num_subunits',
      'score_avg', 'score_max', 'template_domains', 'original_cids')."\n" ;

   print {$in->{out_fh}} '#'.join("\t", 'cid', 'subunit_num', 'seq_id',
      'model_id', 'res_range', 'template_domain', 'domain_class',
      'score_avg', 'score_max')."\n" ;
   foreach my $cid ( @cids) {
      my $complex_summary = $complexes->{summary}->{$cid} ;
      if (exists $in->{start} && ($cid < $in->{start})) {next;}
      my @outvals = ('compl', $cid,
            $complex_summary->{numunits},
            $complex_summary->{score_avg},
            $complex_summary->{score_max},
            join(',', @{$complexes->{details}->{tmpls}->{list}}[@{$complex_summary->{tmpls}}]),
            join(',', @{$complex_summary->{ocid}}) ) ;

      print {$in->{out_fh}} '#'.join("\t", @outvals)."\n" ;

      my $subno = 1 ;
      foreach my $j  (0 .. $#{$complex_summary->{subunits}}) {
         my $hit_ind = $complex_summary->{subunits}->[$j] ;
         my $subunit = $complexes->{details}->{subunits}->{$hit_ind} ;
         my $cur_sid = $complexes->{details}->{tmpls}->{list}->[$complex_summary->{tmpls}->[$j]] ;
         my @outvals = ( $cid, ($j + 1), $subunit->{seq_id},
                         $subunit->{model_id}, $subunit->{resrange},
                         $cur_sid, $complexes->{pb_pl}->{sid2class}->{$cur_sid},
                         $complex_summary->{score_avg},
                         $complex_summary->{score_max} ) ;
         print {$in->{out_fh}} join("\t", @outvals)."\n" ;
         $subno++ ;
      }
   }

}


=head2 covalent_complex_merge()

   Title:       covalent_complex_merge()
   Function:    Merge candidate complexes by identifying complexes that
                contain different parts of same protein.

   Args:        $_->{complexes}
                $_->{seq2dom}->{$seqid}->[i]->{def} - ith domain in $seqid

   Returns:     Nothing

=cut

sub covalent_complex_merge {

   my $in = shift;
   my $complexes = $in->{complexes} ,
   my $seq2dom   = $in->{seq2dom} ;

   my $maxcid_perround = { '0' => $complexes->{counter} };

   print STDERR "Preparing to merge complexes: define LEGO connectors\n" ;
   my $round_cid ;
   my $ocid_cable ; #list of connectable complexes

# Identify those complexes that are 'connectable' - containing sequences
#  with domains that do not participate in the complex - these domains
#  can be used as bridges to other predicted complexes
   my $complexcounter = $complexes->{counter} ;
   foreach my $cid ( sort {$a <=> $b} keys %{$complexes->{summary}}) {
      my $thiscom = {};
      foreach my $hit_ind (@{$complexes->{summary}->{$cid}->{subunits}}) {
         my $subunit = $complexes->{details}->{subunits}->{$hit_ind} ;
         if ($#{$seq2dom->{$subunit->{seq_id}}} == 1) {next;}
         $thiscom->{$subunit->{seq_id}}->{$subunit->{resrange}}++ ;
      }

      my $connectable = 0 ;
      foreach my $seqid (keys %{$thiscom}) {
         my @has =  keys %{$thiscom->{$seqid}} ;
         my $numhas = $#has + 1 ; # how many $seqid domains involved in complex
         my $numfull = $#{$seq2dom->{$seqid}} ; # dont add 1; starts from [1]

# If this sequence has domains that are not involved in this complex,
# it is "connectable" to other complexes involving those domains
         if ($numhas < $numfull) {
            $connectable++ ;

            my $hasstring = $seqid.'.'.join(',',sort @has) ;
            $complexes->{summary}->{$cid}->{hasarch}->{$hasstring}++ ;
         }
      }

# If connectable, define 'hasarch' as what domains are actually present.
      if ($connectable) {
         if (exists $complexes->{summary}->{$cid}->{hasarch}) {
            my $hasarchstring = join(',', sort
               keys %{$complexes->{summary}->{$cid}->{hasarch}});
         }
         push @{$ocid_cable}, $cid ;
      }
   }

   my $hasprot2cid ; #hash pointing from seq_id to connectable complexes with it
   foreach my $tcid ( @{$ocid_cable} ) {
      if (exists $complexes->{summary}->{$tcid}) {
         push @{$round_cid->[0]}, $tcid ;
         foreach my $tocid (@{$complexes->{summary}->{$tcid}->{ocid}}) {
           foreach my $t_hit (@{$complexes->{summary}->{$tocid}->{subunits}}){
               my $subunit = $complexes->{details}->{subunits}->{$t_hit} ;
               $hasprot2cid->{$subunit->{seq_id}}->{$tcid}++; }}}}


   my $round = 1 ; my $maxiterations = 3 ; my $mergechange = 1 ;
#   while ($mergechange && ($round <= $maxiterations)) {
      print STDERR "Merger Round $round\n" ;

      my $knotround = 0;
      if ($round > 1) { $knotround = $round - 1 ; }

      print STDERR "* Will attempt " ;
      my $mergethese ;
      foreach my $j (0 .. ($#{$round_cid->[$knotround]} - 1)) {
         my $cid = $round_cid->[$knotround]->[$j] ;
         if (!exists $complexes->{summary}->{$cid}) {next;}

         my $prot1;
         foreach my $hit_ind (@{$complexes->{summary}->{$cid}->{subunits}}) {
            my $subunit = $complexes->{details}->{subunits}->{$hit_ind} ;
            $prot1->{$subunit->{seq_id}}++ ; }

# find other complexes that share proteins
         my $checkthese ;
         foreach my $prot ( keys %{$prot1} ) {
            map {$checkthese->{$_}++;} (keys %{$hasprot2cid->{$prot}})} ;

         foreach my $cid2 (keys %{$checkthese}) {
            if ($cid == $cid2) {next;}
            if (!exists $complexes->{summary}->{$cid2}) {next;}

# do the complexes share any domains; if so, out
            my $incompatible_fl = 0;
            foreach my $arch (keys %{$complexes->{summary}->{$cid}->{hasarch}}){
               if (exists $complexes->{summary}->{$cid2}->{hasarch}->{$arch}) {
                  $incompatible_fl++ ; last; }}

            if ($incompatible_fl) {next;}
            my $cidstring = join(',', sort {$a <=> $b} ($cid, $cid2)) ;
            $mergethese->{$cidstring}++ ;
         }
      }
      my $num_merge = keys %{$mergethese} ;
      print STDERR "$num_merge mergers\n" ;

      my ($fullcur_fh, $fullcur_fn) = File::Temp::tempfile(
         "complexlist.merger$round.XXXXX", SUFFIX => '.modtie') ;

      my @thisroundcid ;
      $mergechange = 0 ;
      print STDERR "* Now on: -1" ;
      foreach my $cidstring (keys %{$mergethese}) {
         print STDERR "\b"x(length($mergechange - 1)).$mergechange ;
         my ($cid1, $cid2) = split(/\,/, $cidstring) ;
         $complexcounter++ ;
         mergecomplexes({
            complexes => $complexes,
            cid1 => $cid1,
            cid2 => $cid2,
            hasprot2cid => $hasprot2cid,
            complexcounter => $complexcounter,
            out_fh => $fullcur_fh
         }) ;

         push @thisroundcid, $complexcounter ;
         $mergechange++ ;
      }
      print STDERR "\n" ;

# USED TO HAVE MULTIPLE ROUND OF MERGING -- memory-intensive
#
#      $maxcid_perround->{$round} = $complexcounter ;
#      print STDERR "Newly forged complexes (round $round): $mergechange\n" ;
#
#      my ($fullcur_fh, $fullcur_fn) = File::Temp::tempfile(
#         "complexlist.covexpand$round.XXXXX", SUFFIX => '.modtie') ;
#
#      display_complexes({ complexes => $complexes,
#                          out_fh => $fullcur_fh,
#                          start=> ($maxcid_perround->{($round -1)} + 1)});
#      close($fullcur_fh) ;
#
#
#
#      print STDERR "Removing redundancy (post-round $round merging): ";
#      removecomplexredundancy({ complexes => $complexes }) ;
#      print STDERR "X\n" ;
#
#      foreach my $tcid (@thisroundcid) {
#         if (exists $complexes->{$tcid}) {
#            push @{$round_cid->[$round]}, $tcid; }}
#
#      my ($nrcur_fh, $nrcur_fn) = File::Temp::tempfile(
#         "complexlist.covexpand$round.nr.XXXXX", SUFFIX => '.modtie') ;
#      display_complexes({ complexes => $complexes,
#                          out_fh => $nrcur_fh,
#                          start=> ($maxcid_perround->{($round -1)} + 1)});
#      close($nrcur_fh) ;
#
#      $round++ ;
#   }

}


=head2 display_complexes()

   Title:       mergecomplexes()
   Function:    Merge two complexes that contain different parts of same protein
   Args:        $_->{complexes} = complexes data structure
                $_->{cid1} = complex 1
                $_->{cid2} = complex 2
                $_->{complexcounter} = cid of new complex
                $_->{hasprot2cid}->{$seq_id}->{$cid} = hash of cid content
                $_->{out_fh} = output file for merged complex display

=cut

sub mergecomplexes {

   my $in  = shift;
   my $complexes = $in->{complexes};
   my $newcid = $in->{complexcounter};
   my $cid1 = $in->{cid1} ;
   my $cid2 = $in->{cid2} ;
   my $hasprot2cid = $in->{hasprot2cid} ;
   my $out_fh = $in->{out_fh} ;

   $complexes->{summary}->{$newcid} = {};
   my $newcomplex_summ = $complexes->{summary}->{$newcid} ;
   $newcomplex_summ->{hasarch} = {} ;
   $newcomplex_summ->{numunits} = 0 ;
   foreach my $tcid ($cid1, $cid2) {
      foreach my $tocid (@{$complexes->{summary}->{$tcid}->{ocid}}) {
         foreach my $hit_ind (@{$complexes->{summary}->{$tocid}->{subunits}}) {
            my $subunit = $complexes->{details}->{subunits}->{$hit_ind} ;
            $hasprot2cid->{$subunit->{seq_id}}->{$newcid}++ ; }}

      foreach my $hasarch (keys %{$complexes->{summary}->{$tcid}->{hasarch}}) {
         $newcomplex_summ->{hasarch}->{$hasarch}++ ; }

      push @{$complexes->{summary}->{$newcid}->{ocid}},
         @{$complexes->{summary}->{$tcid}->{ocid}} ;

      $newcomplex_summ->{numunits} +=
         $complexes->{summary}->{$tcid}->{numunits} ;
   }

   foreach my $key_type ('tmpls', 'subunits', 'hits') {
      $newcomplex_summ->{$key_type} = [];
      push @{$newcomplex_summ->{$key_type}},
         @{$complexes->{summary}->{$cid1}->{$key_type}} ;
      push @{$newcomplex_summ->{$key_type}},
         @{$complexes->{summary}->{$cid2}->{$key_type}} ;
   }

   $complexes->{summary}->{$newcid}->{score_max} =
      modtie::potentials::maxval_arrref([
                     $complexes->{summary}->{$cid1}->{score_max},
                     $complexes->{summary}->{$cid2}->{score_max}]) ;

   $complexes->{summary}->{$newcid}->{score_avg} = sprintf("%.3f",
                   (($complexes->{summary}->{$cid1}->{score_avg} *
                     $complexes->{summary}->{$cid1}->{numunits} /
                     $newcomplex_summ->{numunits}) +

                    ($complexes->{summary}->{$cid2}->{score_avg} *
                     $complexes->{summary}->{$cid2}->{numunits} /
                     $newcomplex_summ->{numunits}))) ;

# Need to make 2 topostrings: topostring_sids, topostring_hits
# * use old topostring_sids to build actual connectivity - unambiguous edges
# * once connectivity specified, then relabel new graph with sorted hits index
   my $newtopostring_hits ;
   my $newtopostring_sids ;
   {
      my $com1 = $complexes->{summary}->{$cid1} ;
      my $com2 = $complexes->{summary}->{$cid2} ;
      my $com_new = $complexes->{summary}->{$newcid} ;

# make list of sids (old1,old2,new)
      my ($sid2hit, $sid2ind, $sid2sidsorted_ind, $sid12, $sids_sorted) ;
      my @coms = ($com1, $com2, $com_new) ;
      foreach my $j ( 0 .. $#coms) {
         my $com = $coms[$j] ;
         my @sids= @{$complexes->{details}->{tmpls}->{list}}[@{$com->{tmpls}}];
         my @hits= @{$complexes->{details}->{hits}->{list}}[@{$com->{hits}}];
         map {$sid2hit->[$j]->{$sids[$_]} = $hits[$_];} (0 .. $#sids) ;

         my @sids_sorted_by_hit = sort {$sid2hit->[$j]->{$a} cmp
                                        $sid2hit->[$j]->{$b} }
                                       keys %{$sid2hit->[$j]};
         map {$sid2ind->[$j]->{$sids_sorted_by_hit[$_]} = $_ ; }
            ( 0 .. $#sids_sorted_by_hit) ;

         $sids_sorted->[$j] = [sort keys %{$sid2hit->[$j]}] ;
         map {$sid2sidsorted_ind->[$j]->{$sids_sorted->[$j]->[$_]} = $_ ; }
            ( 0 .. $#{$sids_sorted->[$j]}) ;

# make list of sid connections from sid2ind (old1,old2)
         if ($j <= 1) {
         foreach my $sidind_12 (split(',', $com->{topo_sids})) {
            my ($sidind_1, $sidind_2) = split('-', $sidind_12) ;
            push @{$sid12->[$j]}, [$sids_sorted->[$j]->[$sidind_1],
                                   $sids_sorted->[$j]->[$sidind_2]] ;
         }
         }
      }

# make new topolist 
      my $new_edges ;
      foreach my $j (0 .. 1) {
         foreach my $t_sid12 (@{$sid12->[$j]}) {
            my ($t_sid1, $t_sid2) = @{$t_sid12} ;

            my $new_sid_ind1_toposid = $sid2sidsorted_ind->[2]->{$t_sid1} ;
            my $new_sid_ind1_topohit = $sid2ind->[2]->{$t_sid1} ;

            my $new_sid_ind2_toposid = $sid2sidsorted_ind->[2]->{$t_sid2} ;
            my $new_sid_ind2_topohit = $sid2ind->[2]->{$t_sid2} ;

            push @{$new_edges->{topostring_hit}}, join('-',
               sort ($new_sid_ind1_topohit, $new_sid_ind2_topohit)) ;

            push @{$new_edges->{topostring_sid}}, join('-',
               sort ($new_sid_ind1_toposid, $new_sid_ind2_toposid)) ;
         }
      }

      $complexes->{summary}->{$newcid}->{topo_hits} = 
         join(',', sort @{$new_edges->{topostring_hit}});

      $complexes->{summary}->{$newcid}->{topo_sids} = 
         join(',', sort @{$new_edges->{topostring_sid}});

#In-place redundancy removal: if this topo already observed, skip it.
# Otherwise, print it out. then delete it, to prevent memory overload.

      my $md5topo= md5(join(',', sort @{$complexes->{details}->{hits}->{list}}[@{$complexes->{summary}->{$newcid}->{hits}}]).' '.$complexes->{summary}->{$newcid}->{topo_hits});
      if (!exists $complexes->{md5topo}->{$md5topo}) {
         display_complexes({
            out_fh => $out_fh,
            cid => $newcid,
            complexes => $complexes,
         }) ;
      }
      delete $complexes->{summary}->{$newcid} ;
   }
}


sub recurse_composition {

   my $in = shift ;
   my $complexes = $in->{complexes} ;
   my $content = $in->{content} ;
   my $int_list = $in->{intlist} ;
   my $hit_list = $in->{hitlist} ;
   my $tmpltopo = $in->{tmpltopo} ;
   my $z_scores = $in->{zscores} ;

   my $filled = $#{$content} ;
   my $numints = $#{$int_list} ;

   if ($filled == $numints) {
      assess_selfcons({
         complexes => $complexes,
         bdp_id => $in->{bdp_id},
         zscores => $z_scores,
         tmpltopo => $tmpltopo,
         content => $content,
         intlist => $int_list,
         hitlist => $hit_list}) ;
      return;
   }

   my $tofill = $#{$content} + 1 ;

   if ($tofill != 0) {
      my $validity = assess_selfcons({
         complexes => $complexes,
         bdp_id => $in->{bdp_id},
         zscores => $z_scores,
         tmpltopo => $tmpltopo,
         content => $content,
         intlist => $int_list,
         hitlist => $hit_list}) ;
      #dead end elimination
      if ($validity == 0) {return ;}
   }

# preemptive DEE kill on those that are not self-consistent is fine,
# but not on those that are not connected, as addition of componetns may connect
   foreach my $j ( 0 .. $#{$hit_list->{$int_list->[$tofill]}}) {
      $content->[$tofill] = $j ;
      recurse_composition({
         complexes => $complexes,
         bdp_id => $in->{bdp_id},
         zscores => $z_scores,
         tmpltopo => $tmpltopo,
         content => $content,
         intlist => $int_list,
         hitlist => $hit_list}) ;
   }

}


sub assess_selfcons {
#purpose: Check a putative complex to make sure all adjacent arcs are
#        compatible in the target domains they claim
# Iterate over putative target complex; if not self-consistent, then kill off.
# If self-consistent and topology is connected, add to complex list

   my $in = shift ;
   my $complexes = $in->{complexes} ;
   my $content = $in->{content} ;
   my $int_list = $in->{intlist} ;
   my $hit_list = $in->{hitlist} ;
   my $tmpltopo = $in->{tmpltopo} ;
   my $z_scores = $in->{zscores} ;

   my $maxscore ;
   my $score = 0; my $numedges = 0;
   my $assn ; #points from template domain to hit
# Iterate over each interface in the template BDP complex
   foreach my $j ( 0 .. $#{$content}) { #j = interface number
# If 0, the interface is absent in the proposed target complex
      if ($content->[$j] == 0) {next;}

      my $sid12 = $int_list->[$j] ;
      $numedges++ ;

      my $t_mod12 = $hit_list->{$sid12}->[($content->[$j] - 1)];
      $score += $z_scores->{$sid12}->{$t_mod12} ;
      if (!defined $maxscore || $maxscore < $z_scores->{$sid12}->{$t_mod12}){
         $maxscore = $z_scores->{$sid12}->{$t_mod12}; }

      my ($t_mod1, $t_mod2) = split(/\n/, $t_mod12) ;
      my ($t_sid1, $t_sid2) = split(/\n/, $sid12) ;

# If a hit has already been specified for either of the template domains
# (in a previous interface in this template complex) make sure it is the same
#  as the one proposed for that domain in the current interface
      if (exists $assn->{$t_sid1}) {
         if ($assn->{$t_sid1} ne $t_mod1) {
            return 0;}
      } else { $assn->{$t_sid1} = $t_mod1;}

      if (exists $assn->{$t_sid2}) {
         if ($assn->{$t_sid2} ne $t_mod2) {
            return 0;}
      } else { $assn->{$t_sid2} = $t_mod2;}
   }

# Make lists of the template domains and their hits in the target complex
   my @sids = sort keys %{$assn} ;
   my @assns = @{$assn}{@sids} ;

# If binary complex, no need to check if it is connected
   if ($#assns < 2) {return 1;}

# If higher-order complex, make sure it is connected and figure out the topology
   my $connects = is_connected({
      hitlist => $hit_list,
      tmpltopo => $tmpltopo,
      content => $content,
      intlist => $int_list}) ;

   if ($connects->{is_connected}) {
      $score /= $numedges ;
      $score = sprintf("%.3f", $score) ;

      $complexes->{counter}++ ;

      $complexes->{summary}->{$complexes->{counter}} = {
         numunits => ($#assns + 1),
         ocid => [$complexes->{counter}],
         score_avg => $score,
         score_max => $maxscore,
         topo_hits => $connects->{topostring_hits},
         topo_sids => $connects->{topostring_sids},
      } ;

      foreach my $j (0 .. $#sids) {
         push @{$complexes->{summary}->{$complexes->{counter}}->{hits}},
            $complexes->{details}->{hits}->{hash}->{$assns[$j]} ;
         push @{$complexes->{summary}->{$complexes->{counter}}->{subunits}},
            $complexes->{details}->{hits}->{hash}->{$assns[$j]} ;

         if (!exists $complexes->{details}->{tmpls}->{hash}->{$sids[$j]}) {
            push @{$complexes->{details}->{tmpls}->{list}}, $sids[$j];
            $complexes->{details}->{tmpls}->{hash}->{$sids[$j]} =
               $#{$complexes->{details}->{tmpls}->{list}} ; }

         push @{$complexes->{summary}->{$complexes->{counter}}->{tmpls}},
            $complexes->{details}->{tmpls}->{hash}->{$sids[$j]} ;
      }
   }

   return 1;
}



=head2 is_connected()

   Title:       is_connected()
   Function:    Checks a putative complex to ensure it has a connected topology

   Args:        $_->{intlist}
                $_->{hitlist}
                $_->{content}
                $_->{tmpltopo}

   Returns:     ->{numcomp} = number of connected components in target complex
                ->{is_connected} = 0 (if numcomp > 1); otherwise 1
                ->{topostring} = nodestring (sorted target domain list).' '
                                 edgestring (sorted edgelist)

=cut

sub is_connected {

# now based on edges from binaryset - use template strx instead
   my $in = shift ;
   my $int_list = $in->{intlist} ;      #->[i] = ith $sid12 interface
   my $content  = $in->{content} ;      #->[i] = which hit to use for the ith
                                        #        interface in $int_list
                                         
   my $hit_list = $in->{hitlist} ;      #->{$sid12}->[i] = $hit12

   my $tmpltopo = $in->{tmpltopo} ;     #->{edges} = list of $sid12 in complex
                                        #->{chains}->{$sid12} = same|diff|both

   my $indeptopostring = {        #NODE = "hit" (seqid/modelid/resrange)
      nodecounter => 0   } ;      #->{nodecounter} = number of hits
                                  #->{nodelabel}->{$hit} = index
                                  #->{edgelist} = list of $hit12 interfaces
                                  #->{nodelist} = unsorted list of $hits
                                  #->{nodestring} = join(', ') sorted nodelist

# assign node numbers to domains from int_list; then fill graph using $content
   my $graph ;          # target complex graph connectivity matrix;
                        #    NODE=template domain

   my $sid2ind ;            #->{$sid} = index
   my $sid2sidsorted_ind ;  #->{$sid} = index
   my @nodes ;              #->[i] = ith $sid in target complex
   my $domspres  ;          #->{$sid} = { seq => $seq_id
                            #             mod => $model_id
                            #             res => $residue_range
                            #             hit => $hit }

#Number graph nodes based on alphanumeric sort of hits
# (Although could be ambiguous - because homodimers)
# currently sorted on $sid, but equally ambiguous for comparison
# to other topologies
# know that comparing topology strings is not perfect; will let some perfect
#       matches through
# BUT: err towards being too loose on filtering, rather than too strict
# 1: in case of multiple copies, add onto hit a 'number' just for disambiguation
#    purposes.


# Number nodes based on alphanumeric comparison of hits.
#  not perfect; partly ambiguous; same graph can have multiple topo strings
#                                 in case of multiple copies of same hit
   {
      my $sid2hit ;
      foreach my $j ( 0 .. $#{$int_list}) {
         my $sid12 = $int_list->[$j] ;
         if (!defined $content->[$j] || $content->[$j] <= 0) {next;}
         my ($sid1, $sid2) = split(/\n/, $sid12) ;
         my ($hit1, $hit2) =
            split(/\n/, $hit_list->{$sid12}->[($content->[$j] - 1)]) ;
         $sid2hit->{$sid1} = $hit1 ;
         $sid2hit->{$sid2} = $hit2 ;
      }
      my @sids_sorted_by_hit = sort {$sid2hit->{$a} cmp $sid2hit->{$b}}
                                 keys %{$sid2hit} ;
      map {$sid2ind->{$sids_sorted_by_hit[$_]} = $_ ; }
            ( 0 .. $#sids_sorted_by_hit) ;

      my @sids_sorted = sort keys %{$sid2hit} ;
      map {$sid2sidsorted_ind->{$sids_sorted[$_]} = $_ ; }
            ( 0 .. $#sids_sorted) ;
   }


# Iterate over interfaces in target complex to build topology string
   foreach my $j ( 0 .. $#{$int_list}) {
      my $sid12 = $int_list->[$j] ;
      if (!defined $content->[$j] || $content->[$j] <= 0) {next;}
      my ($sid1, $sid2) = split(/\n/, $sid12) ;

      my ($hit1, $hit2) =
            split(/\n/, $hit_list->{$sid12}->[($content->[$j] - 1)]) ;

      my ($seq1, $mod1, $res1) = split(/\^/, $hit1) ;
      my ($seq2, $mod2, $res2) = split(/\^/, $hit2) ;

      $domspres->{$sid1} = {seq => $seq1, mod => $mod1,
                            res => $res1, hit => $hit1} ;

      $domspres->{$sid2} = {seq => $seq2, mod => $mod2,
                            res => $res2, hit => $hit2} ;

      $indeptopostring->{nodelabel}->{$hit1} = $sid2ind->{$sid1} ;
      $indeptopostring->{nodelabel}->{$hit2} = $sid2ind->{$sid2} ;

      my $hitlab12sig = join('-', sort {$a <=> $b}
                             ($indeptopostring->{nodelabel}->{$hit1},
                              $indeptopostring->{nodelabel}->{$hit2})) ;
      push @{$indeptopostring->{edgelist}}, $hitlab12sig ;

      my $sidlab12sig = join('-', sort {$a <=> $b}
                             ($sid2sidsorted_ind->{$sid1},
                              $sid2sidsorted_ind->{$sid2})) ;
      push @{$indeptopostring->{sid_edgelist}}, $sidlab12sig;

      push @{$graph->[$sid2ind->{$sid1}]}, $sid2ind->{$sid2} ;
      push @{$graph->[$sid2ind->{$sid2}]}, $sid2ind->{$sid1} ;
   }

   my @domspres = keys %{$domspres} ;
   foreach my $sid (@domspres) {
      push @{$indeptopostring->{nodelist}},
         $domspres->{$sid}->{hit} ; }

   $indeptopostring->{nodestring} = join(', ', sort 
      @{$indeptopostring->{nodelist}}) ;

# Check all-vs-all target complex domains for intra interface; populate $graph
   foreach my $j (0 .. ($#domspres - 1)) {
      my $sid1 = $domspres[$j] ;

      foreach my $k (($j + 1) .. $#domspres) {
         my $sid2 = $domspres[$k] ;

# are the two target domains different ranges of the same sequence 
# and do these two domains  have an intra-chain interface in the template
         if (($domspres->{$sid1}->{seq} eq $domspres->{$sid2}->{seq}) &&
             ($domspres->{$sid1}->{res} ne $domspres->{$sid2}->{res}) &&

            ((exists $tmpltopo->{chains}->{$sid1."\n".$sid2} &&
                     $tmpltopo->{chains}->{$sid1."\n".$sid2} ne 'diff') ||

             (exists $tmpltopo->{chains}->{$sid2."\n".$sid1} &&
                     $tmpltopo->{chains}->{$sid2."\n".$sid1} ne 'diff'))) {

# now check if both target domains are from the same sequence - but not the same domain - eg if A-A make sure target(A1) is not same stretch as target(A2)

#               print "#intra FOUND: $sid1 ($domspres->{$sid1}->{seq} $domspres->{$sid1}->{res} ) $sid2 ($domspres->{$sid2}->{seq} $domspres->{$sid2}->{res} )\n" ;
            push @{$graph->[$sid2ind->{$sid1}]}, $sid2ind->{$sid2} ;
            push @{$graph->[$sid2ind->{$sid2}]}, $sid2ind->{$sid1} ;

            my $hit1 = $domspres->{$sid1}->{hit} ;
            my $hit2 = $domspres->{$sid2}->{hit} ;
            my $hitlab12sig= join('-', sort {$a <=> $b}
                                    ($indeptopostring->{nodelabel}->{$hit1},
                                     $indeptopostring->{nodelabel}->{$hit2})) ;
            push @{$indeptopostring->{edgelist}}, $hitlab12sig ;

            my $sidlab12sig = join('-', sort {$a <=> $b}
                              ($sid2sidsorted_ind->{$sid1},
                               $sid2sidsorted_ind->{$sid2})) ;
            push @{$indeptopostring->{sid_edgelist}}, $sidlab12sig;
         }
      }
   }

# Edge string includes both inter and intra interfaces
   $indeptopostring->{edgestring} = join(',',
      sort @{$indeptopostring->{edgelist}});

# Determine target complex connectivity
   my $conncomp = dfs_conncomp($graph) ;
   my $numcomp = $#{$conncomp} + 1 ;

# don't actually need nodestring on the hits topo
#  - is just sorted @tmpls sorted by @hits.
   my $backinfo = { numcomp => $numcomp,
                    topostring_hits => $indeptopostring->{edgestring},
                    topostring_sids => 
                     join(',', sort @{$indeptopostring->{sid_edgelist}}) } ;

   if ($numcomp > 1) { $backinfo->{is_connected} = 0; }
   else              { $backinfo->{is_connected} = 1; }

   return $backinfo;

}


sub dfs_conncomp {
   my $graph = shift ;

# input: graph (adj list)
# output: component number

            
   my $color  = [ ] ;
   my $parent = [ ] ;
   my $compcont = [ ] ; #component contents

   foreach my $vertex (0 .. $#{$graph}) {
      $color->[$vertex] = 0 ;           #0 = white, 1 = grey, 2 = black ;
      $parent->[$vertex] = $vertex ;    #0 = white, 1 = grey, 2 = black ;
   }
      
   my $counter = -1 ;
      
   foreach my $vertex (0 .. $#{$graph}) {
      if ($color->[$vertex] == 0) {
         $counter++ ;
         $compcont->[$counter] = [ ] ;
         dfs_visit($vertex, $graph, $counter, $parent, $color, $compcont) ;
      }
   }


   my @components;
   foreach my $k (0 .. $#{$compcont}) {
      push @components, join(',', sort {$a <=> $b} @{$compcont->[$k]}) ;
   }

   return \@components ;
}


sub dfs_visit {
   my ($vertex, $graph, $counter, $parent, $color, $compcont) = @_ ;

   no warnings 'recursion';

   $color->[$vertex] = 1 ;
   push @{$compcont->[$counter]}, $vertex ;
   if ($graph->[$vertex][0] ne 'nil') {
      foreach my $adj (@{$graph->[$vertex]}) {
         if ($color->[$adj] == 0) {
            $parent->[$adj] = $vertex ;
            dfs_visit($adj, $graph, $counter, $parent, $color, $compcont) ;
         }
      }
   }
   $color->[$vertex] = 2 ;

}


sub load_interface_clusters {

   my $in = shift ;
   my $specs ;
   if (exists $in->{specs}) { $specs = $in->{specs} ; }
   else { $specs = modtie::set_modtie_specs() ; }

   open (ICLUST, $specs->{interface_clusters_fn}) ;
   my $sidpair2clus ;
   my $fampairclus ;
   while (my $line = <ICLUST>) {
      chomp $line;
      if ($line =~ /ERROR/) {next;}
      if ($line =~ / trial /) {last;}
      my ($fam1, $fam2, $clusno, $sid1, $sid2, $intraclusno)=split(/\t/, $line);
#      print STDERR "fam1 =$fam1, fam2 = $fam2, sid1 = $sid1, sid2 = $sid2  clusno = $clusno, intraclusno = $intraclusno\n" ;
      $sidpair2clus->{$sid1."\n".$sid2} = $fam1."\n".$fam2."\n".$clusno ;
      $fampairclus->{$fam1."\n".$fam2."\n".$clusno}->{$sid1."\n".$sid2}++ ;
   }
   close(ICLUST) ;

   return {
      sidpair2clus => $sidpair2clus,
      fampairclus => $fampairclus,
   } ;

}

sub list_complexes_preload {

   my ($sid2bdp, $sid2class) ;
   {
      my @tod_res = modtie::pibase::rawselect_tod(
         "SELECT bdp_id, subset_id, class FROM subsets") ;
      foreach my $j ( 0 .. $#{$tod_res[0]}) {
         if ($tod_res[0]->[$j] eq "" ||
             $tod_res[1]->[$j] !~ /SCOP/) {next;}
         $sid2bdp->{$tod_res[1]->[$j]} = $tod_res[0]->[$j] ;
         $sid2class->{$tod_res[1]->[$j]} = $tod_res[2]->[$j] ;
      }
   }

   my ($bdp2sid12, $sid12chains) ;
   {
      my @tod_res = modtie::pibase::rawselect_tod(
         "SELECT bdp_id, subset_id_1, subset_id_2, chains ".
         "FROM intersubset_contacts") ;
      foreach my $j ( 0 .. $#{$tod_res[0]}) {
         if ($tod_res[1]->[$j] !~ /SCOP/) {next;}
         my $sid12 = $tod_res[1]->[$j]."\n".$tod_res[2]->[$j] ;
         push @{$bdp2sid12->{$tod_res[0]->[$j]}}, $sid12 ;
         $sid12chains->{$sid12} = $tod_res[3]->[$j] ;
      }
   }


   my $pb = {
      sid2bdp => $sid2bdp,
      bdp2sid12 => $bdp2sid12,
      sid12chains => $sid12chains,
      sid2class => $sid2class,
   } ;
   return $pb ;

}


=head2 memusage()

   Title:       memusage()
   Function:    Returns the VmSize of the process
      works on linux /proc/$$/status

=cut

sub memusage {
   my $memusage = `grep ^VmSize /proc/$$/status` ;
   print STDERR "memusage is $memusage\n" ;
   chop $memusage ; $memusage =~ s/^VmSize:\s*//g ;
   return $memusage ;
}


1 ;
