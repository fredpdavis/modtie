=head1 NAME

yeast.pm - Yeast-specific routines for prediction and assessment

=head1 DESCRIPTION

The yeast.pm perl library contains subroutines to assess the predicted
yeast interactions using sub-cellular localization, function annotation,
and known interactions from MIPS, BIND and CELLZOME (2006) data.

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


package modtie::yeast ;
use strict;
use warnings;
use Exporter ;
use Cwd qw/getcwd/;
our @ISA = qw/Exporter/ ;
our @EXPORT_OK = qw/load_cellzome assess_yeast_results/ ;
push @EXPORT_OK, qw/calc_p_b_prior/ ;

use modtie ;


sub assess_yeast_results {

   require DBI ;

   my $in = shift ;

   print STDERR "Loading YEASTgfp data: " ;
   my $yeastgfp = readin_yeastgfp() ;
   print STDERR "X\n" ;

   print STDERR "Loading SGD GO annotation: " ;
   my $sgdgo = readin_sgd_go() ;
   print STDERR "X\n" ;

# venn overlap of modtie binaries with BIND yeast binary interactions
# put all into SGD naming

   print STDERR "Loading modbase seqid - SGD YORF translation table: " ;
   my $seqidsgd = load_seqid_2_sgd() ;
   my $seqid_2_sgd = $seqidsgd->{seqid_2_sgd} ;
   my $sgd_2_seqid = $seqidsgd->{sgd_2_seqid} ;
   print STDERR "X\n" ;

   print STDERR "Loading pibase interfaces: ";
   my $pb;
   {
      my ($dbh) = pibaseconnect_pibase() ;
      $pb->{sid12chains} = pibase::mysql_hashload($dbh,
      "SELECT concat(subset_id_1,'\\n',subset_id_2), chains ".
      "FROM intersubset_contacts WHERE subset_id_1 LIKE '\%SCOP\%'") ;

      $pb->{sid2class} = pibase::mysql_hashload($dbh,
      "SELECT subset_id, class FROM subsets where bdp_id IS NOT NULL AND ".
      "subset_id LIKE '\%SCOP\%'") ;

      $pb->{bdp2sid} = pibase::mysql_hasharrload($dbh,
      "SELECT bdp_id, subset_id FROM subsets where bdp_id IS NOT NULL ".
      "AND subset_id LIKe '\%SCOP\%'") ;

      $pb->{fam12} = pibase::mysql_hashload($dbh,
      "SELECT concat(least(class_1,class_2), '\\n', ".
      "greatest(class_1, class_2)) ".
      "FROM intersubset_contacts WHERE num_contacts >= 1000 AND ".
      "subset_id_1 LIKE '\%SCOP\%'") ;
   }
   print STDERR "X\n" ;

   print STDERR "Loading yeast domain architectures:" ;
   my $domdata ; $domdata->{seqinfo} = {} ;

   my $seqid_2_domains = readin_seqid_aaseq_domarch({
      seqinfo => $domdata->{seqinfo},
      fn => $in->{seqid_aaseq_arch_fn}}) ;

   my $mbinfo = _parse_domarchs($domdata) ;

   my ( $scop2seq, $seq2dom, $seq2domfrag, $scop2num ) =
      ($mbinfo->{scop2seq},
       $mbinfo->{seq2dom},
       $mbinfo->{seq2domfrag},
       $mbinfo->{scop2num}) ;
   print STDERR "X\n" ;


   if (!exists $in->{binary_fn}) {
      die "assess_yeast_results(): need the binary_fn filename" ; }
   my $modtie_binint;
   my $modtie_binint_domains;

   my $head;
   @{$head->{mt_bin}->{f}} = qw/INTERFACE subset_id_1 subset_id_2 seq1 model_id_1 t_def1 class_1 seq2 model_id_2 t_def2 class_2 tmpl_contacts aln_contacts pot_file pot_source pot_type pot_dist raw z z_p1 z_p2 avg min min_tn stdev z_min z_min_tn false_pos rmsd_dom1 equivpos_dom1 numres_tmpl_dom1 numres_targ_dom1 numres_ident_dom1 numres_tmpl_bs_dom1 numres_targ_bs_dom1 numres_ident_bs_dom1 rmsd_dom2 equivpos_dom2 numres_tmpl_dom2 numres_targ_dom2 numres_ident_dom2 numres_tmpl_bs_dom2 numres_targ_bs_dom2 numres_ident_bs_dom2/ ;
   $head->{mt_bin}->{f2i} = pibase::array2hash($head->{mt_bin}->{f}) ;

   my $filters;
   my $readerrs = {
      tabtab => 0,
      zundef => 0,
      samechainsamedom => 0,
   };

   my $input_thresh = {
      aln_tmpl_contacts => 0.5,
   } ;

   my $rawpred ;
   my $mtint_maxseqid ;
   my $mtint_bestzentry ;
   print STDERR "Loading binary predictions: " ;
   open(MTINT, $in->{binary_fn}) ;
   while (my $line = <MTINT>) {
      if ($line !~ /\n$/) {next;}
      chomp $line;

      my @t  = split(/\t/, $line) ;
      if (!defined  $t[$head->{mt_bin}->{f2i}->{z}] ||
          $t[$head->{mt_bin}->{f2i}->{z}] eq '') {$readerrs->{zundef}++ ; next;}

      if (($t[$head->{mt_bin}->{f2i}->{aln_contacts}] /
           $t[$head->{mt_bin}->{f2i}->{tmpl_contacts}]) <
           $input_thresh->{aln_tmpl_contacts} ) {next;}

      my $sid1 = $t[$head->{mt_bin}->{f2i}->{subset_id_1}] ;
      my $sid2 = $t[$head->{mt_bin}->{f2i}->{subset_id_2}] ;
      my $sid12 = $sid1."\n".$sid2 ;
      if ($sid2 lt $sid1) {$sid12 = $sid2."\n".$sid1;}


      my $seqid1 = $t[$head->{mt_bin}->{f2i}->{seq1}] ;
      my $seqid2 = $t[$head->{mt_bin}->{f2i}->{seq2}] ;

      my $chains = $pb->{sid12chains}->{$sid12} ;
      if ($chains eq 'same' &&
          !($seqid1 eq $seqid2 &&
            $t[$head->{mt_bin}->{f2i}->{t_def1}] eq
            $t[$head->{mt_bin}->{f2i}->{t_def2}])) {
            
         $readerrs->{samechainsamedom}++; next;}

      my $sgd1 = $seqid_2_sgd->{$seqid1} ;
      my $sgd2 = $seqid_2_sgd->{$seqid2} ;

      my $sig1dom = $sgd1."\t".$t[$head->{mt_bin}->{f2i}->{t_def1}] ;
      my $sig2dom = $sgd2."\t".$t[$head->{mt_bin}->{f2i}->{t_def2}] ;
      my $sgdsigdom = $sig1dom."\t".$sig2dom ;
      if ($sig2dom lt $sig1dom) {
         $sgdsigdom = $sig2dom."\t".$sig1dom ; }

      my $sgdsig = $sgd1."\t".$sgd2 ;
      my $sig1 = $sgd1; my $sig2 = $sgd2 ;
      if ($sgd2 lt $sgd1) {
         $sig2 =$sgd1; $sig1 = $sgd2 ;
         $sgdsig = $sgd2."\t".$sgd1 ; }

      my $avg_seqid = sprintf("%.2f", (($t[$head->{mt_bin}->{f2i}->{numres_ident_dom1}] / $t[$head->{mt_bin}->{f2i}->{numres_targ_dom1}]) + ($t[$head->{mt_bin}->{f2i}->{numres_ident_dom2}] / $t[$head->{mt_bin}->{f2i}->{numres_targ_dom2}])) * 50) ;

      my $recnum = $#{$rawpred} + 1 ;

      if (!exists $modtie_binint->{$sgdsig}) {
         $mtint_maxseqid->{$sgdsig} = $avg_seqid ;
         $mtint_bestzentry->{$sgdsig} = $recnum ;
         $modtie_binint->{$sgdsig} = $t[$head->{mt_bin}->{f2i}->{z}] ; }


      if (!exists $modtie_binint_domains->{$sgdsig}->{$sgdsigdom}) {
         $modtie_binint_domains->{$sgdsig}->{$sgdsigdom} = $t[$head->{mt_bin}->{f2i}->{z}] ; }

      if ($t[$head->{mt_bin}->{f2i}->{z}] < $modtie_binint->{$sgdsig}) {
            $mtint_bestzentry->{$sgdsig} = $recnum ;
            $modtie_binint->{$sgdsig} = $t[$head->{mt_bin}->{f2i}->{z}] ; }

      if ($avg_seqid > $mtint_maxseqid->{$sgdsig}) {
         $mtint_maxseqid->{$sgdsig} = $avg_seqid ; }

      if ($t[$head->{mt_bin}->{f2i}->{z}] < $modtie_binint_domains->{$sgdsig}->{$sgdsigdom}) {
            $modtie_binint_domains->{$sgdsig}->{$sgdsigdom} = $t[$head->{mt_bin}->{f2i}->{z}] ; }

      $filters->{none}->{protprot}->{$sgdsig}++  ;
      $filters->{none}->{prot}->{$sgd1}++  ;
      $filters->{none}->{prot}->{$sgd2}++  ;
      $filters->{none}->{domdom}->{$sgdsigdom}++  ;
      $filters->{none}->{dom}->{$sig1dom}++  ;
      $filters->{none}->{dom}->{$sig2dom}++  ;
      push @{$rawpred}, \@t ;
   }
   print STDERR " errors (tabtab=".$readerrs->{tabtab}.
                  "; zundef=".$readerrs->{zundef}.
                  "; samechainsamedom=".$readerrs->{samechainsamedom} ;
   print STDERR "X\n" ;

   close(MTINT) ;



#   my $z_thresh = -2 ;
   my $z_thresh = -1.7 ;


   @{$head->{mt_compl}->{f}} = qw/cid subunitid seq_id model_id def bdp_id subset_id z maxz/ ;
   $head->{mt_compl}->{f2i} = pibase::array2hash($head->{mt_compl}->{f}) ;


   my $cid2comments ;
   my $bin2comments ;

   open(MTCOMPL, $in->{complexes_fn});
   print STDERR "Loading complex predictions: " ;
   my $mtcompl ;
   while (my $line = <MTCOMPL>) {
      chomp $line;
      if ($line =~ /^\#/) {next;}
      my @t = split(/\t/, $line) ;
      my $sgd = $seqid_2_sgd->{$t[$head->{mt_compl}->{f2i}->{seq_id}]} ;

      if ($t[$head->{mt_compl}->{f2i}->{maxz}] > $z_thresh) {next;}
      push @{$mtcompl->{cid2sgd}->{$t[$head->{mt_compl}->{f2i}->{cid}]}}, {
         sgd => $sgd,
         def=> $t[$head->{mt_compl}->{f2i}->{def}],
         model_id => $t[$head->{mt_compl}->{f2i}->{model_id}],
         subset_id => $t[$head->{mt_compl}->{f2i}->{subset_id}],
         seq_id => $t[$head->{mt_compl}->{f2i}->{seq_id}],
         z => $t[$head->{mt_compl}->{f2i}->{z}],
         maxz => $t[$head->{mt_compl}->{f2i}->{maxz}]
      } ;
      $mtcompl->{sgd2cid}->{$sgd}->{$t[$head->{mt_compl}->{f2i}->{cid}]}++ ;
   }
   print STDERR "X\n" ;
   close(MTCOMPL) ;


   print STDERR "   Binary filtering Z: " ;
   my $filterset_z ;
   foreach my $sgdsig ( keys %{$modtie_binint_domains} ) {
      foreach my $sgdsigdom ( keys %{$modtie_binint_domains->{$sgdsig}}) {
         my ($sgd1, $dom1, $sgd2, $dom2) = split(/\t/, $sgdsigdom) ;
         if ($modtie_binint_domains->{$sgdsig}->{$sgdsigdom} <= $z_thresh) {
            $filters->{z}->{dom}->{$sgd1."\t".$dom1}++ ;
            $filters->{z}->{dom}->{$sgd2."\t".$dom2}++ ;
            $filters->{z}->{domdom}->{$sgdsigdom}++ ;
            $filters->{z}->{prot}->{$sgd1}++ ;
            $filters->{z}->{prot}->{$sgd2}++ ;
            $filters->{z}->{protprot}->{$sgd1."\t".$sgd2}++ ;
            $filterset_z->{$sgdsig}->{$sgdsigdom}++ ;
         }
      }
   }
   print STDERR "X\n" ;


#make SCOP graph

   print STDERR "SCOP Graph: ";
   my $scopedges ;
   my ($scopedges_fa, $scopedges_sf) ;
   foreach my $fam12 (keys %{$pb->{fam12}}) {
      my ($fam1, $fam2) = split(/\n/, $fam12) ;
      my ($rfam1, $rfam2) = sort ($fam1, $fam2) ;
      print "PDBDIRECT\t$rfam1\t$rfam2\n" ;
      if ($rfam1 eq $rfam2) {next;}
      $scopedges->{$rfam1}->{$rfam2} = 'black' ;
      $scopedges_fa->{$rfam1}->{$rfam2}->{'PDBDIRECT'}++ ;

      my ($sf1) = ($fam1 =~ /([a-z]\.[0-9]+\.[0-9]+)/) ;
      my ($sf2) = ($fam2 =~ /([a-z]\.[0-9]+\.[0-9]+)/) ;
      if ($sf1 eq $sf2) {next;}

      my @ordsf = sort ($sf1, $sf2) ;
      $scopedges_sf->{$ordsf[0]}->{$ordsf[1]}->{'PDBDIRECT'}++ ;
   }

   foreach my $bdp (keys %{$pb->{bdp2sid}}) {
      my $domarr = $pb->{bdp2sid}->{$bdp} ;
      foreach my $j ( 0 .. ($#{$domarr} - 1)) {
         my $dom1 = $domarr->[$j] ;
         my $fam1 = $pb->{sid2class}->{$domarr->[$j]} ;
         foreach my $k ( ($j + 1) .. $#{$domarr}) {

            my $dom2 = $domarr->[$k] ;
            my $fam2 = $pb->{sid2class}->{$domarr->[$k]} ;

#            if (exists $pb->{sid12chains}->{$dom1."\n".$dom2} ||
#                exists $pb->{sid12chains}->{$dom2."\n".$dom1}) {next;}

            my ($rfam1, $rfam2) = sort ($fam1, $fam2) ;
            if ($rfam1 eq $rfam2) {next;}
            if (!exists $scopedges->{$rfam1}->{$rfam2}) {
               print "PDBCOCOMP\t$rfam1\t$rfam2\n" ;
               $scopedges->{$rfam1}->{$rfam2} = 'blue' ;
               $scopedges_fa->{$rfam1}->{$rfam2}->{'PDBCOCOMP'}++ ;

               my ($sf1) = ($fam1 =~ /([a-z]\.[0-9]+\.[0-9]+)/) ;
               my ($sf2) = ($fam2 =~ /([a-z]\.[0-9]+\.[0-9]+)/) ;
               my @ordsf = sort ($sf1, $sf2) ;

               if (($sf1 ne $sf2) && (!exists
                  $scopedges_sf->{$ordsf[0]}->{$ordsf[1]}->{'PDBDIRECT'})) {
                  $scopedges_sf->{$ordsf[0]}->{$ordsf[1]}->{'PDBCOCOMP'}++ ;
               }
            }
         }
      }
   }

   foreach my $cid (keys %{$mtcompl->{cid2sgd}}) {
      my @sids ;
      foreach my $subunit (@{$mtcompl->{cid2sgd}->{$cid}}) {
         push @sids, $subunit->{subset_id} ; }

      foreach my $j ( 0 .. ($#sids - 1)) {
         my $fam1 = $pb->{sid2class}->{$sids[$j]} ;
         foreach my $k ( ($j + 1) .. $#sids) {
            my $fam2 = $pb->{sid2class}->{$sids[$k]} ;
#            if (exists $pb->{sid12chains}->{$sids[$j]."\n".$sids[$k]} ||
#                exists $pb->{sid12chains}->{$sids[$k]."\n".$sids[$j]}) {next;}
            my ($rfam1, $rfam2) = sort ($fam1, $fam2) ;
            if ($rfam1 eq $rfam2) {next;}

            if (!exists $scopedges->{$rfam1}->{$rfam2}) {
               print "PREDCOCOMP\t$rfam1\t$rfam2\n" ;
               $scopedges->{$rfam1}->{$rfam2} = 'red' ;
            }

            if (!exists $scopedges_fa->{$rfam1}->{$rfam2}) {
               print "PREDCOCOMP\t$rfam1\t$rfam2\n" ;
               $scopedges_fa->{$rfam1}->{$rfam2}->{'PREDCOCOMP'}++;

               my ($sf1) = ($fam1 =~ /([a-z]\.[0-9]+\.[0-9]+)/) ;
               my ($sf2) = ($fam2 =~ /([a-z]\.[0-9]+\.[0-9]+)/) ;
               if ($sf1 eq $sf2)  {next;}
               my @ordsf = sort ($sf1, $sf2) ;

               if (!exists $scopedges_sf->{$ordsf[0]}->{$ordsf[1]}) {
                  $scopedges_sf->{$ordsf[0]}->{$ordsf[1]}->{'PREDCOCOMP'}++ ;
               }
            }
         }
      }
   }
   print STDERR "X\n" ;
#   build_scop_graph({edges => $scopedges}) ;


   print STDERR "Binary checks: \n" ;
   foreach my $key (sort {$b cmp $a} keys %{$filters->{z}}) {
      my $tnum = keys %{$filters->{z}->{$key}} ;
      print STDERR "   PASS z_thresh <= $z_thresh\t$key\t$tnum\n" ;
   }
   print STDERR "X\n" ;

#   summarize_yeast_bin_pred() ;


   print STDERR "Loading BIND yeast data: " ;
   my $bind_yeast = load_bind_yeast() ;
   my $bind_div = load_bind_division() ;
   my $bind_yeast_int = $bind_yeast->{nrints} ;
   my $bind_yeast_complex = $bind_yeast->{complex} ;
   print STDERR "X\n" ;

   print STDERR "Loading CELLZOME yeast data: " ;
   my $cellzome = load_cellzome() ;
   my $cz_binint = $cellzome->{binary};
   my $cz_complex = $cellzome->{complex};
   print STDERR "X\n" ;

#   print STDERR "Loading MIPS complexes: " ;
#   my $mips_complexes = load_mips_complexes() ;

   my $bind_binint ;
   my $exp_binint_filters ;
   foreach my $rgid (keys %{$bind_yeast_int->{tax2rgid}->{"4932\n4932"}}) {
      if (exists $bind_yeast_int->{nrints}->{$rgid}->{a_sgd} &&
          exists $bind_yeast_int->{nrints}->{$rgid}->{b_sgd}) {
         my $sgd1 = $bind_yeast_int->{nrints}->{$rgid}->{a_sgd};
         my $sgd2 = $bind_yeast_int->{nrints}->{$rgid}->{b_sgd} ;
         my $sgdsig = $sgd1."\t".$sgd2 ;
         if ($sgd2 lt $sgd1) {
            $sgdsig = $sgd2."\t".$sgd1 ; }
         $bind_binint->{$sgdsig} = $rgid ;

         my $divs = {};
         $divs->{'BIND'}++ ;
         if (exists $bind_div->{rgid2div}->{$rgid}->{'SGD'}){
            $divs->{'SGD'}++; }
         if (exists $bind_div->{rgid2div}->{$rgid}->{'MIPS'}){
            $divs->{'MIPS'}++; }


         my $flags ;

         if (exists $yeastgfp->{orf2lochash}->{$sgd1} &&
             exists $yeastgfp->{orf2lochash}->{$sgd2}) {
            my $tcoloc = 0 ;
            foreach my $l (keys %{$yeastgfp->{orf2lochash}->{$sgd1}}) {
               if (exists $yeastgfp->{orf2lochash}->{$sgd2}->{$l}) {
                     $tcoloc++;}
            }
            if ($tcoloc > 0) {$flags->{'l'}++;}
         }


         if (exists $sgdgo->{orf2go}->{$sgd1} &&
            exists $sgdgo->{orf2go}->{$sgd2}) {
            my $tcofunc  ;
            $tcofunc->{'P'} = 0 ;
            $tcofunc->{'C'} = 0 ;
            $tcofunc->{'F'} = 0 ;
            $tcofunc->{'all'} = 0 ;
            foreach my $goterm (keys %{$sgdgo->{orf2go}->{$sgd1}}) {
               if (exists $sgdgo->{orf2go}->{$sgd2}->{$goterm}) {
                  $tcofunc->{$sgdgo->{go2type}->{$goterm}}++ ;
                  $tcofunc->{all}++ ;
               }
            }
            if ($tcofunc->{'all'} > 0) {$flags->{'f'}++;}
         }

#fpd060415 fixed homo pass
         if ($sgd1 eq $sgd2) {
            $flags->{'l'}++ ;
            $flags->{'f'}++ ;
         }
         
         foreach my $flag (keys %{$flags}) {
            foreach my $div (keys %{$divs}) {
               $exp_binint_filters->{$div}->{$flag}->{protprot}->{$sgdsig}++ ;
               $exp_binint_filters->{$div}->{$flag}->{prot}->{$sgd1}++ ;
               $exp_binint_filters->{$div}->{$flag}->{prot}->{$sgd2}++ ;
            }
         }

         foreach my $div (keys %{$divs}) {
            $exp_binint_filters->{$div}->{none}->{protprot}->{$sgdsig}++ ;
            $exp_binint_filters->{$div}->{none}->{prot}->{$sgd1}++ ;
            $exp_binint_filters->{$div}->{none}->{prot}->{$sgd2}++ ;

            if (exists $flags->{'f'} && exists $flags->{'l'}) {
               $exp_binint_filters->{$div}->{'lf'}->{protprot}->{$sgdsig}++ ;
               $exp_binint_filters->{$div}->{'lf'}->{prot}->{$sgd1}++ ;
               $exp_binint_filters->{$div}->{'lf'}->{prot}->{$sgd2}++ ;
            }
         }
      }
   }

   print STDERR "   filtering BIND complexes:\n" ;
   my $bind_complex_filters ;
   open (OUTF, ">complex_filter.bind")  ;
   foreach my $bcid (keys %{$bind_yeast_int->{tax2rgid}->{"4932\n4932"}}) {
      my @sgds ;
      my $domfams ;
      foreach my $subunit (@{$bind_yeast_complex->{complex2subunits}->{$bcid}}){
         if (!exists $subunit->{'sgd'}) {next; }
         else {push @sgds, $subunit->{sgd};}

         my $t_seqid = $sgd_2_seqid->{$subunit->{sgd}} ;
#         print STDERR "bcid $bcid: ".$subunit->{sgd}." == $t_seqid\n" ;
         if (!exists $seq2dom->{$t_seqid}) {next;}

         foreach my $alldomind (1 .. $#{$seq2dom->{$t_seqid}}) {
            my $alldom = $seq2dom->{$t_seqid}->[$alldomind] ;
#            print STDERR "   domain $alldomind: ".$alldom->{scop}."\n" ;
            if ($alldom->{scop} !~ /\.[0-9]+\./) {next;}
            $domfams->{$alldom->{scop}}++ ;
         }
      }
      my $compostring = join(',', sort @sgds) ;

      my @domfams = keys %{$domfams} ;
      foreach my $j (0 .. ($#domfams - 1)) {
         my $t_fam1 = $domfams[$j] ;
         foreach my $k (($j + 1) .. $#domfams ) {
            my $t_fam2 = $domfams[$k] ;
            if ($t_fam1 eq $t_fam2) {next;}

            if (($t_fam1 =~ /\..*\..*\./) && ($t_fam2 =~ /\..*\..*\./)) {
               my @ordfam = sort ($t_fam1, $t_fam2) ;
               $scopedges_fa->{$ordfam[0]}->{$ordfam[1]}->{'BIND'}++ ;
               $scopedges_fa->{$ordfam[0]}->{$ordfam[1]}->{'ALLEXP'}++ ;
            }
            
            my ($sf1) = ($t_fam1 =~ /([a-z]\.[0-9]+\.[0-9]+)/) ;
            my ($sf2) = ($t_fam2 =~ /([a-z]\.[0-9]+\.[0-9]+)/) ;
            if ($sf1 ne $sf2) {
                  my @ordsf = sort ($sf1, $sf2) ;
                  $scopedges_sf->{$ordsf[0]}->{$ordsf[1]}->{'BIND'}++ ;
                  $scopedges_sf->{$ordsf[0]}->{$ordsf[1]}->{'ALLEXP'}++ ;
            }
         }
      }


#OHFUCK - dont remember why this cutoff placed: fpd060306_1753
# ONLY HIGHER ORDER COMLEXES (>= 3 subunits)
      if ($#sgds < 1) {next;}

      my $totalpairs = ($#sgds * ($#sgds + 1)) / 2 ;

      my $samefunc = 0 ;
      my $sameloc = 0 ; 
      my $samelf = 0 ;
      foreach my $j ( 0 .. ($#sgds - 1)) {
         my $sgd1 = $sgds[$j] ;
         foreach my $k ( ($j + 1) .. $#sgds) {
            my $sgd2 = $sgds[$k] ;
            my $thisone ;
            if (exists $yeastgfp->{orf2lochash}->{$sgd1} &&
                exists $yeastgfp->{orf2lochash}->{$sgd2}) {
               my $tcoloc = 0 ;
               foreach my $l (keys %{$yeastgfp->{orf2lochash}->{$sgd1}}) {
                  if (exists $yeastgfp->{orf2lochash}->{$sgd2}->{$l}) {
                        $tcoloc++;}
               }
               if ($tcoloc > 0) {$sameloc++; $thisone->{l}++ ;}
            }

            if (exists $sgdgo->{orf2go}->{$sgd1} &&
               exists $sgdgo->{orf2go}->{$sgd2}) {
               my $tcofunc  ;
               $tcofunc->{'P'} = 0 ;
               $tcofunc->{'C'} = 0 ;
               $tcofunc->{'F'} = 0 ;
               $tcofunc->{'all'} = 0 ;
               foreach my $goterm (keys %{$sgdgo->{orf2go}->{$sgd1}}) {
                  if (exists $sgdgo->{orf2go}->{$sgd2}->{$goterm}) {
                     $tcofunc->{$sgdgo->{go2type}->{$goterm}}++ ;
                     $tcofunc->{all}++ ;
                  }
               }
               if ($tcofunc->{'all'} > 0) {$samefunc++; $thisone->{f}++ ;}
            }

            if (exists $thisone->{f} && exists $thisone->{l}) { $samelf++ ; }
         }
      }

      $bind_complex_filters->{sgdunits}->{$bcid} = $#sgds + 1 ;
      $bind_complex_filters->{f}->{$bcid} = sprintf("%.3f", ($samefunc / $totalpairs)) ;
      $bind_complex_filters->{l}->{$bcid} = sprintf("%.3f", ($sameloc / $totalpairs)) ;
      $bind_complex_filters->{lf}->{$bcid} = sprintf("%.3f", ($samelf / $totalpairs)) ;

      my @outvals = ( $bcid,
                      $compostring,
                      $bind_complex_filters->{sgdunits}->{$bcid},
                      $bind_complex_filters->{sgdunits}->{$bcid},
                      $bind_complex_filters->{f}->{$bcid},
                      $bind_complex_filters->{l}->{$bcid},
                      $bind_complex_filters->{lf}->{$bcid} ) ;
      print OUTF join(" ", @outvals)."\n" ;
   }
   close(OUTF) ;

   my ($set_pdbdir, $set_pdbcocomp, $set_predcocomp, $set_bind) ;
   print STDERR "   filtering CELLZOME interactions:\n" ;
   foreach my $sgdsig (keys %{$cz_binint}) {
      my ($sgd1, $sgd2) = split(/\t/, $sgdsig)  ;
         my $flags ;

         if (exists $yeastgfp->{orf2lochash}->{$sgd1} &&
             exists $yeastgfp->{orf2lochash}->{$sgd2}) {
            my $tcoloc = 0 ;
            foreach my $l (keys %{$yeastgfp->{orf2lochash}->{$sgd1}}) {
               if (exists $yeastgfp->{orf2lochash}->{$sgd2}->{$l}) {
                     $tcoloc++;}
            }
            if ($tcoloc > 0) {$flags->{'l'}++;}
         }

         if (exists $sgdgo->{orf2go}->{$sgd1} &&
            exists $sgdgo->{orf2go}->{$sgd2}) {
            my $tcofunc  ;
            $tcofunc->{'P'} = 0 ;
            $tcofunc->{'C'} = 0 ;
            $tcofunc->{'F'} = 0 ;
            $tcofunc->{'all'} = 0 ;
            foreach my $goterm (keys %{$sgdgo->{orf2go}->{$sgd1}}) {
               if (exists $sgdgo->{orf2go}->{$sgd2}->{$goterm}) {
                  $tcofunc->{$sgdgo->{go2type}->{$goterm}}++ ;
                  $tcofunc->{all}++ ;
               }
            }
            if ($tcofunc->{'all'} > 0) {$flags->{'f'}++;}
         }

#fpd060415 fixed homo pass
         if ($sgd1 eq $sgd2) {
            $flags->{'l'}++ ;
            $flags->{'f'}++ ;
         }

         my $divs = {} ;
         $divs->{'cz'}++ ;
         
         if ($cz_binint->{$sgdsig} >= 5) { $divs->{'czreliable'}++;}
         
         foreach my $flag (keys %{$flags}) {
            foreach my $div (keys %{$divs}) {
               $exp_binint_filters->{$div}->{$flag}->{protprot}->{$sgdsig}++ ;
               $exp_binint_filters->{$div}->{$flag}->{prot}->{$sgd1}++ ;
               $exp_binint_filters->{$div}->{$flag}->{prot}->{$sgd2}++ ;
            }
         }

         foreach my $div (keys %{$divs}) {
            $exp_binint_filters->{$div}->{none}->{protprot}->{$sgdsig}++ ;
            $exp_binint_filters->{$div}->{none}->{prot}->{$sgd1}++ ;
            $exp_binint_filters->{$div}->{none}->{prot}->{$sgd2}++ ;

            if (exists $flags->{'f'} && exists $flags->{'l'}) {
               $exp_binint_filters->{$div}->{'lf'}->{protprot}->{$sgdsig}++ ;
               $exp_binint_filters->{$div}->{'lf'}->{prot}->{$sgd1}++ ;
               $exp_binint_filters->{$div}->{'lf'}->{prot}->{$sgd2}++ ;
            }
         }
   }

   print STDERR "   filtering CELLZOME complexes:\n" ;
   my $cz_complex_filters ;
   open (OUTF, ">complex_filter.cellzome")  ;
   foreach my $bcid (keys %{$cz_complex->{complex2subunits}}) {
      my @sgds ;
      my $domfams ;
      foreach my $subunit (@{$cz_complex->{complex2subunits}->{$bcid}}){
         push @sgds, $subunit->{sgd};

         my $t_seqid = $sgd_2_seqid->{$subunit->{sgd}} ;
#         print STDERR "bcid $bcid: ".$subunit->{sgd}." == $t_seqid\n" ;
         if (!exists $seq2dom->{$t_seqid}) {next;}

         foreach my $alldomind (1 .. $#{$seq2dom->{$t_seqid}}) {
            my $alldom = $seq2dom->{$t_seqid}->[$alldomind] ;
#            print STDERR "   domain $alldomind: ".$alldom->{scop}."\n" ;
            if ($alldom->{scop} !~ /\.[0-9]+\./) {next;}
            $domfams->{$alldom->{scop}}++ ;
         }
      }
      my $compostring = join(',', sort @sgds) ;

      my @domfams = keys %{$domfams} ;
      foreach my $j (0 .. ($#domfams - 1)) {
         my $t_fam1 = $domfams[$j] ;
         foreach my $k (($j + 1) .. $#domfams ) {
            my $t_fam2 = $domfams[$k] ;
            if ($t_fam1 eq $t_fam2) {next;}

            if (($t_fam1 =~ /\..*\..*\./) && ($t_fam2 =~ /\..*\..*\./)) {
               my @ordfam = sort ($t_fam1, $t_fam2) ;
               $scopedges_fa->{$ordfam[0]}->{$ordfam[1]}->{'CZ'}++ ;
               $scopedges_fa->{$ordfam[0]}->{$ordfam[1]}->{'ALLEXP'}++ ;
            }
            
            my ($sf1) = ($t_fam1 =~ /([a-z]\.[0-9]+\.[0-9]+)/) ;
            my ($sf2) = ($t_fam2 =~ /([a-z]\.[0-9]+\.[0-9]+)/) ;
            if ($sf1 ne $sf2) {
                  my @ordsf = sort ($sf1, $sf2) ;
                  $scopedges_sf->{$ordsf[0]}->{$ordsf[1]}->{'CZ'}++ ;
                  $scopedges_sf->{$ordsf[0]}->{$ordsf[1]}->{'ALLEXP'}++ ;
            }
         }
      }


#OHFUCK - dont remember why this cutoff placed: fpd060306_1753
# ONLY HIGHER ORDER COMLEXES (>= 3 subunits)
      if ($#sgds < 1) {next;}

      my $totalpairs = ($#sgds * ($#sgds + 1)) / 2 ;

      my $samefunc = 0 ;
      my $sameloc = 0 ; 
      my $samelf = 0 ;
      foreach my $j ( 0 .. ($#sgds - 1)) {
         my $sgd1 = $sgds[$j] ;
         foreach my $k ( ($j + 1) .. $#sgds) {
            my $sgd2 = $sgds[$k] ;
            my $thisone ;
            if (exists $yeastgfp->{orf2lochash}->{$sgd1} &&
                exists $yeastgfp->{orf2lochash}->{$sgd2}) {
               my $tcoloc = 0 ;
               foreach my $l (keys %{$yeastgfp->{orf2lochash}->{$sgd1}}) {
                  if (exists $yeastgfp->{orf2lochash}->{$sgd2}->{$l}) {
                        $tcoloc++;}
               }
               if ($tcoloc > 0) {$sameloc++; $thisone->{l}++ ;}
            }

            if (exists $sgdgo->{orf2go}->{$sgd1} &&
               exists $sgdgo->{orf2go}->{$sgd2}) {
               my $tcofunc  ;
               $tcofunc->{'P'} = 0 ;
               $tcofunc->{'C'} = 0 ;
               $tcofunc->{'F'} = 0 ;
               $tcofunc->{'all'} = 0 ;
               foreach my $goterm (keys %{$sgdgo->{orf2go}->{$sgd1}}) {
                  if (exists $sgdgo->{orf2go}->{$sgd2}->{$goterm}) {
                     $tcofunc->{$sgdgo->{go2type}->{$goterm}}++ ;
                     $tcofunc->{all}++ ;
                  }
               }
               if ($tcofunc->{'all'} > 0) {$samefunc++; $thisone->{f}++ ;}
            }

            if (exists $thisone->{f} && exists $thisone->{l}) { $samelf++ ; }
         }
      }

      $cz_complex_filters->{sgdunits}->{$bcid} = $#sgds + 1 ;
      $cz_complex_filters->{f}->{$bcid} = sprintf("%.3f", ($samefunc / $totalpairs)) ;
      $cz_complex_filters->{l}->{$bcid} = sprintf("%.3f", ($sameloc / $totalpairs)) ;
      $cz_complex_filters->{lf}->{$bcid} = sprintf("%.3f", ($samelf / $totalpairs)) ;

      my @outvals = ( $bcid,
                      $compostring,
                      $cz_complex_filters->{sgdunits}->{$bcid},
                      $cz_complex_filters->{sgdunits}->{$bcid},
                      $cz_complex_filters->{f}->{$bcid},
                      $cz_complex_filters->{l}->{$bcid},
                      $cz_complex_filters->{lf}->{$bcid} ) ;
      print OUTF join(" ", @outvals)."\n" ;
   }
   close(OUTF) ;

#   my ($set_pdbdir, $set_pdbcocomp, $set_predcocomp, $set_bind) ;
   my $fampair_nums = {
      'BIND' => 0,
      'CZ' => 0,
      'BIND-CZ' => 0,
      'PDBDIRECT' => 0,
      'PDBCOCOMP' => 0,
      'PREDCOCOMP' => 0,
      'BIND-PDBDIRECT' => 0,
      'BIND-PDBCOCOMP' => 0,
      'BIND-PREDCOCOMP' => 0,
      'CZ-PDBDIRECT' => 0,
      'CZ-PDBCOCOMP' => 0,
      'CZ-PREDCOCOMP' => 0,
   } ;

   open(DOMPAIRSF, ">domainpairs_venn_list.$$.out") ;

   my $sfpair_nums = {
      'BIND' => 0,
      'CZ' => 0,
      'BIND-CZ' => 0,
      'PDBDIRECT' => 0,
      'PDBCOCOMP' => 0,
      'PREDCOCOMP' => 0,
      'BIND-PDBDIRECT' => 0,
      'BIND-PDBCOCOMP' => 0,
      'BIND-PREDCOCOMP' => 0,
      'CZ-PDBDIRECT' => 0,
      'CZ-PDBCOCOMP' => 0,
      'CZ-PREDCOCOMP' => 0,
   } ;

   {
      foreach my $fam1 (keys %{$scopedges_fa}) {
         foreach my $fam2 (keys %{$scopedges_fa->{$fam1}}) {
            foreach my $type (keys %{$scopedges_fa->{$fam1}->{$fam2}}) {
               print DOMPAIRSF "FAMPAIR $type $fam1 $fam2\n";
               $fampair_nums->{$type}++ ;
            }
            my @types = keys %{$scopedges_fa->{$fam1}->{$fam2}} ;
            if (exists $scopedges_fa->{$fam1}->{$fam2}->{'BIND'} &&
               $#types > 0) {
               foreach my $type (@types) {
                  if ($type eq 'BIND') {next;}
                  $fampair_nums->{'BIND-'.$type}++ ;
               }
            }
            if (exists $scopedges_fa->{$fam1}->{$fam2}->{'CZ'} &&
               $#types > 0) {
               foreach my $type (@types) {
                  if ($type eq 'BIND' || $type eq 'CZ') {next;}
                  $fampair_nums->{'CZ-'.$type}++ ;
               }
            }
            if (exists $scopedges_fa->{$fam1}->{$fam2}->{'ALLEXP'} &&
               $#types > 0) {
               foreach my $type (@types) {
                  if ($type eq 'ALLEXP') {next;}
                  $fampair_nums->{'ALLEXP-'.$type}++ ;
               }
            }
         }
      }

      foreach my $sf1 (keys %{$scopedges_sf}) {
         foreach my $sf2 (keys %{$scopedges_sf->{$sf1}}) {
            foreach my $type (keys %{$scopedges_sf->{$sf1}->{$sf2}}) {
               print DOMPAIRSF "SFPAIR $type $sf1 $sf2\n";
               $sfpair_nums->{$type}++ ;
            }
            my @types = keys %{$scopedges_sf->{$sf1}->{$sf2}} ;
            if (exists $scopedges_sf->{$sf1}->{$sf2}->{'BIND'} &&
               $#types > 0) {
               foreach my $type (@types) {
                  if ($type eq 'BIND') {next;}
                  $sfpair_nums->{'BIND-'.$type}++ ;
               }
            }

            if (exists $scopedges_sf->{$sf1}->{$sf2}->{'CZ'} &&
               $#types > 0) {
               foreach my $type (@types) {
                  if ($type eq 'BIND' || $type eq 'CZ') {next;}
                  $sfpair_nums->{'CZ-'.$type}++ ;
               }
            }

            if (exists $scopedges_sf->{$sf1}->{$sf2}->{'ALLEXP'} &&
               $#types > 0) {
               foreach my $type (@types) {
                  if ($type eq 'ALLEXP') {next;}
                  $sfpair_nums->{'ALLEXP-'.$type}++ ;
               }
            }

         }
      }
   }

   foreach my $key (keys %{$fampair_nums}) {
      print "FAMPAIRNUMS\t$key\t$fampair_nums->{$key}\n" ; }

   foreach my $key (keys %{$sfpair_nums}) {
      print "SFPAIRNUMS\t$key\t$sfpair_nums->{$key}\n" ; }


   print STDERR "   filtering MODTIE complexes:\n" ;
   my $complex_filters ;
   open (OUTF, ">complex_filter.modtie")  ;
   foreach my $cid (keys %{$mtcompl->{cid2sgd}}) {
      my @sgds ;
      my $uniqsgd ;
      foreach my $subunit (@{$mtcompl->{cid2sgd}->{$cid}}) {
         $uniqsgd->{$subunit->{sgd}}++ ;
         push @sgds, $subunit->{sgd} ; }

# ONLY HIGHER ORDER COMLEXES (>= 3 subunits)

      my $numunits = $#sgds + 1;
      if ($numunits < 3) {next;}

      @sgds = keys %{$uniqsgd} ;
      my $compostring = join(',', sort @sgds) ;
      my $numtypes = $#sgds + 1 ;

# at least 2 types of proteins AND more than or equal to 3 proteins
      if ($numtypes < 2) {next;}

      my $t_expcompl ;
      foreach my $sgd ( @sgds) {
         if (exists $bind_yeast_complex->{sgd2bcid}->{$sgd}) {
            foreach my $bcid (keys %{$bind_yeast_complex->{sgd2bcid}->{$sgd}}) {
               $t_expcompl->{'BIND'.$bcid}++;}
         }
         if (exists $cz_complex->{sgd2bcid}->{$sgd}) {
            foreach my $bcid (keys %{$cz_complex->{sgd2bcid}->{$sgd}}) {
               $t_expcompl->{'CZ'.$bcid}++;}
         }
      }
      my @expcomplexes ;
      foreach my $bcid ( keys %{$t_expcompl}) {
         if ($t_expcompl->{$bcid} == $numtypes ) {
            push @expcomplexes, $bcid ; } }

      my $totalpairs = ($#sgds * ($#sgds + 1)) / 2 ;
      my $samefunc = 0 ;
      my $sameloc = 0 ; 
      my $samelf = 0 ;
      foreach my $j ( 0 .. ($#sgds - 1)) {
         my $sgd1 = $sgds[$j] ;
         foreach my $k ( ($j + 1) .. $#sgds) {
            my $sgd2 = $sgds[$k] ;
            my $thisone ;
            if (exists $yeastgfp->{orf2lochash}->{$sgd1} &&
                exists $yeastgfp->{orf2lochash}->{$sgd2}) {
               my $tcoloc = 0 ;
               foreach my $l (keys %{$yeastgfp->{orf2lochash}->{$sgd1}}) {
                  if (exists $yeastgfp->{orf2lochash}->{$sgd2}->{$l}) {
                        $tcoloc++;}
               }
               if ($tcoloc > 0) {$sameloc++; $thisone->{l}++ ;}
            }

            if (exists $sgdgo->{orf2go}->{$sgd1} &&
               exists $sgdgo->{orf2go}->{$sgd2}) {
               my $tcofunc  ;
               $tcofunc->{'P'} = 0 ;
               $tcofunc->{'C'} = 0 ;
               $tcofunc->{'F'} = 0 ;
               $tcofunc->{'all'} = 0 ;
               foreach my $goterm (keys %{$sgdgo->{orf2go}->{$sgd1}}) {
                  if (exists $sgdgo->{orf2go}->{$sgd2}->{$goterm}) {
                     $tcofunc->{$sgdgo->{go2type}->{$goterm}}++ ;
                     $tcofunc->{all}++ ;
                  }
               }
               if ($tcofunc->{'all'} > 0) {$samefunc++; $thisone->{f}++ ;}
            }

            if (exists $thisone->{f} && exists $thisone->{l}) { $samelf++ ; }
         }
      }

      if ($#expcomplexes >= 0 ) {
         $complex_filters->{num_exper_match}->{$cid} = $#expcomplexes + 1;
         my @newexpcompl ;
         foreach my $cur_bcid (@expcomplexes) {
            my $t_numsubunits ;
            if ($cur_bcid =~ /^BIND/) {
               my $t_bcid = $cur_bcid ; $t_bcid =~ s/^BIND// ;
               $t_numsubunits = $#{$bind_yeast_complex->{complex2subunits}->{$t_bcid}} + 1;
            } elsif ($cur_bcid =~ /^CZ/) {
               my $t_bcid = $cur_bcid ; $t_bcid =~ s/^CZ// ;
               $t_numsubunits = $#{$cz_complex->{complex2subunits}->{$t_bcid}} + 1;
            }
            push @newexpcompl, "$cur_bcid($t_numsubunits)" ;
         }
         $complex_filters->{exper_match}->{$cid} = join(",", @newexpcompl) ;
         if ($complex_filters->{exper_match}->{$cid} =~ /BIND/) {
            $cid2comments->{$cid}->{'BIND'}++; }
         if ($complex_filters->{exper_match}->{$cid} =~ /CZ/) {
            $cid2comments->{$cid}->{'CZ'}++; }
      } else {
         $complex_filters->{num_exper_match}->{$cid} = 0 ;
         $complex_filters->{exper_match}->{$cid} = '' ;
      }

      $complex_filters->{numunits}->{$cid} = $numunits ;
      $complex_filters->{numtypes}->{$cid} = $numtypes ;
      $complex_filters->{f}->{$cid} = sprintf("%.3f", ($samefunc / $totalpairs)) ;
      $complex_filters->{l}->{$cid} = sprintf("%.3f", ($sameloc / $totalpairs)) ;
      $complex_filters->{lf}->{$cid} = sprintf("%.3f", ($samelf / $totalpairs)) ;
      if ($complex_filters->{f}->{$cid} == 1) {
         $cid2comments->{$cid}->{'cofunc'}++; }

      if ($complex_filters->{l}->{$cid} == 1) {
         $cid2comments->{$cid}->{'coloc'}++; }

      my @outvals = ( $cid,
                      $compostring,
                      $complex_filters->{numunits}->{$cid},
                      $complex_filters->{numtypes}->{$cid},
                      $complex_filters->{f}->{$cid},
                      $complex_filters->{l}->{$cid},
                      $complex_filters->{lf}->{$cid},
                      $complex_filters->{num_exper_match}->{$cid},
                      $complex_filters->{exper_match}->{$cid} ) ;
      print OUTF join(" ", @outvals)."\n" ;
   }
   close(OUTF) ;


   print STDERR "Biological Filters:\n\n" ;
   print STDERR "Fig 1d: Frequency of Co=Localization and Co-FunctionEnrichment of Filters\n\n";
   print STDERR 'postscript("filter_histo.eps",horizontal=FALSE,onefile=FALSE,pointsize=16, height=8,width=8)'."\n";
   print STDERR 'par(mar=c(4.1, 4.1, 1, 0.5))'."\n";

   foreach my $div (qw/BIND SGD MIPS cz czreliable/) {
      my $totalnum = keys %{$exp_binint_filters->{$div}->{none}->{protprot}} ;
      my $filtnum ;
      my $p ;
      foreach my $flag (qw/f l lf/) {
         $filtnum->{$flag} = keys %{$exp_binint_filters->{$div}->{$flag}->{protprot}} ;
         $p->{$flag} = sprintf("%.2f", $filtnum->{$flag} / $totalnum) ;
#         print STDERR "$div\t$flag\t$p->{$flag}\t$filtnum->{$flag}\n" ;
      }
      print STDERR "$div<-c($p->{f}, $p->{l}, $p->{lf})\n" ;
   }


   my ($yeast_coloc, $yeast_cofunc) ;
   foreach my $sgdsig (keys %{$filterset_z}) {
      my ($sgd1, $sgd2) = split(/\t/, $sgdsig) ;
      my $seqsig = join("\n",
         sort ($sgd_2_seqid->{$sgd1},$sgd_2_seqid->{$sgd2})) ;

#fpd060415_1806 FIX homo filter ; even if not available, homo interactions
#  pass cofunc/coloc filter by default - thank you referee 2

      if ($sgd1 eq $sgd2) {
         $filters->{zl}->{protprot}->{$sgd1."\t".$sgd2}++ ;
         $filters->{zl}->{prot}->{$sgd1}++ ;
         $filters->{zl}->{prot}->{$sgd2}++ ;

         $filters->{zf}->{protprot}->{$sgd1."\t".$sgd2}++ ;
         $filters->{zf}->{prot}->{$sgd1}++ ;
         $filters->{zf}->{prot}->{$sgd2}++ ;

         $filters->{zlf}->{protprot}->{$sgd1."\t".$sgd2}++ ;
         $filters->{zlf}->{prot}->{$sgd1}++ ;
         $filters->{zlf}->{prot}->{$sgd2}++ ;

         foreach my $sgdsigdom (keys %{$filterset_z->{$sgdsig}}) {
            my ($sgd1, $dom1, $sgd2, $dom2) = split(/\t/, $sgdsigdom) ;
            $filters->{zl}->{domdom}->{$sgdsigdom}++ ;
            $filters->{zl}->{dom}->{$sgd1."\t".$dom1}++ ;
            $filters->{zl}->{dom}->{$sgd2."\t".$dom2}++ ;

            $filters->{zf}->{domdom}->{$sgdsigdom}++ ;
            $filters->{zf}->{dom}->{$sgd1."\t".$dom1}++ ;
            $filters->{zf}->{dom}->{$sgd2."\t".$dom2}++ ;

            $filters->{zlf}->{domdom}->{$sgdsigdom}++ ;
            $filters->{zlf}->{dom}->{$sgd1."\t".$dom1}++ ;
            $filters->{zlf}->{dom}->{$sgd2."\t".$dom2}++ ;
         }

         next;
      }


      if (exists $yeastgfp->{orf2lochash}->{$sgd1} &&
          exists $yeastgfp->{orf2lochash}->{$sgd2}) {
         my $tcoloc = 0 ;
         foreach my $l (keys %{$yeastgfp->{orf2lochash}->{$sgd1}}) {
            if (exists $yeastgfp->{orf2lochash}->{$sgd2}->{$l}) {
                  $tcoloc++;}
         }
         $yeast_coloc->{$sgd1}->{$sgd2} = $tcoloc ;
         $yeast_coloc->{$sgd2}->{$sgd1} = $tcoloc ;

         if ($tcoloc) {
            $bin2comments->{$seqsig}->{'coloc'}++ ;
#            print STDERR "SET COLOC FUCK: ".join(" ", sort ($sgd_2_seqid->{$sgd1},
#               $sgd_2_seqid->{$sgd2}))." has coloc\n" ;
         }
      }

      if (exists $sgdgo->{orf2go}->{$sgd1} &&
          exists $sgdgo->{orf2go}->{$sgd2}) {
         my $tcofunc  ;
         $tcofunc->{'P'} = 0 ;
         $tcofunc->{'C'} = 0 ;
         $tcofunc->{'F'} = 0 ;
         $tcofunc->{'all'} = 0 ;
         foreach my $goterm (keys %{$sgdgo->{orf2go}->{$sgd1}}) {
            if (exists $sgdgo->{orf2go}->{$sgd2}->{$goterm}) {
               $tcofunc->{$sgdgo->{go2type}->{$goterm}}++ ;
               $tcofunc->{all}++ ;
            }
            $yeast_cofunc->{$sgd1}->{$sgd2} = $tcofunc ;
            $yeast_cofunc->{$sgd2}->{$sgd1} = $tcofunc ;
         }

         if ($tcofunc->{'all'} > 0 ) {
#            print STDERR "SET COFUNC FUCK : ".join(" ", sort ($sgd_2_seqid->{$sgd1},
#               $sgd_2_seqid->{$sgd2}))." has cofunc\n" ;
            $bin2comments->{$seqsig}->{'cofunc'}++ ; }
      }

      if (exists $yeast_coloc->{$sgd1}->{$sgd2}) {
         if ($yeast_coloc->{$sgd1}->{$sgd2} > 0 ) {
            $filters->{zl}->{protprot}->{$sgd1."\t".$sgd2}++ ;
            $filters->{zl}->{prot}->{$sgd1}++ ;
            $filters->{zl}->{prot}->{$sgd2}++ ;

            foreach my $sgdsigdom (keys %{$filterset_z->{$sgdsig}}) {
               my ($sgd1, $dom1, $sgd2, $dom2) = split(/\t/, $sgdsigdom) ;
               $filters->{zl}->{domdom}->{$sgdsigdom}++ ;
               $filters->{zl}->{dom}->{$sgd1."\t".$dom1}++ ;
               $filters->{zl}->{dom}->{$sgd2."\t".$dom2}++ ;
            }
         }
      }

      if (exists $yeast_cofunc->{$sgd1}->{$sgd2}) {
         if ($yeast_cofunc->{$sgd1}->{$sgd2}->{all} > 0 ) {
            $filters->{zf}->{protprot}->{$sgdsig}++ ;
            $filters->{zf}->{prot}->{$sgd1}++ ;
            $filters->{zf}->{prot}->{$sgd2}++ ;

            if (exists $filters->{zl}->{protprot}->{$sgdsig}) {
               $filters->{zlf}->{protprot}->{$sgdsig}++ ;
               $filters->{zlf}->{prot}->{$sgd1}++ ;
               $filters->{zlf}->{prot}->{$sgd2}++ ;
            }

            foreach my $sgdsigdom (keys %{$filterset_z->{$sgdsig}}) {
               my ($sgd1, $dom1, $sgd2, $dom2) = split(/\t/, $sgdsigdom) ;
               $filters->{zf}->{domdom}->{$sgdsigdom}++ ;
               $filters->{zf}->{dom}->{$sgd1."\t".$dom1}++ ;
               $filters->{zf}->{dom}->{$sgd2."\t".$dom2}++ ;

               if (exists $filters->{zl}->{protprot}->{$sgdsig}) {
                  $filters->{zlf}->{domdom}->{$sgdsigdom}++ ;
                  $filters->{zlf}->{dom}->{$sgd1."\t".$dom1}++ ;
                  $filters->{zlf}->{dom}->{$sgd2."\t".$dom2}++ ;
               }
            }
         }
      }
   }

   {
      open(ZORFP, ">orfpairs.modtie.z.$$.out") ;
      foreach my $sgdsig (keys %{$filters->{z}->{protprot}}) {
         print ZORFP $sgdsig."\n" ; } close(ZORFP) ;

      open(ZLORFP, ">orfpairs.modtie.zl.$$.out") ;
      foreach my $sgdsig (keys %{$filters->{zl}->{protprot}}) {
         print ZLORFP $sgdsig."\n" ; } close(ZLORFP) ;

      open(ZFORFP, ">orfpairs.modtie.zf.$$.out") ;
      foreach my $sgdsig (keys %{$filters->{zf}->{protprot}}) {
         print ZFORFP $sgdsig."\n" ; } close(ZFORFP) ;

      open(ZLFORFP, ">orfpairs.modtie.zlf.$$.out") ;
      foreach my $sgdsig (keys %{$filters->{zlf}->{protprot}}) {
         print ZLFORFP $sgdsig."\n" ; } close(ZLFORFP) ;



      my $totalnum   = keys %{$filters->{z}->{protprot}};
      my $filtnum ;
      my $p ;
      if ($totalnum == 0) {print STDERR "OH FUCK : no overlaps in the filter z set\n" ;}
      else {
      foreach my $flag (qw/zf zl zlf/) {
         $filtnum->{$flag} = keys %{$filters->{$flag}->{protprot}} ;
         $p->{$flag} = sprintf("%.2f", $filtnum->{$flag} / $totalnum) ;
#         print STDERR "pred\t$flag\t$p->{$flag}\t$filtnum->{$flag}\n" ;
      }
      print STDERR "pred<-c($p->{zf}, $p->{zl}, $p->{zlf})\n" ;
      }
   }

   print STDERR "barplot(rbind(pred,BIND,cz,czreliable),beside=T,space=c(0,1.5),legend.text=c('Predicted', 'BIND','CZ', 'CZ Reliable'),names.arg=c('Co-Function','Co-Localization','Co-Loc + Co-Func'), density=c(200,0,15,15),angle=c(90,0,0,135),col=c('black','black','black','black'), ylim=c(0,1), ylab='Relative Frequency', xlab='Biological Filter', cex=0.5, cex.lab=1.25, cex.names=1, cex.axis=1.25)\n\n" ;
   print STDERR "dev.off()\n";


   {
      my $allexp_sgdsigs;
      my $numbind = 0;
      foreach my $tsig (keys %{$bind_binint}) {
         $numbind++ ;
         $allexp_sgdsigs->{$tsig}++;}

      my $numczrel = 0;
      foreach my $tsig (keys %{$cz_binint}) {
         if ($cz_binint->{$tsig} >= 5) {
            $numczrel++ ;
            $allexp_sgdsigs->{$tsig}++;}
      }
      my $numkeys = keys %{$allexp_sgdsigs} ;
      open(EXPORFSIG, ">orfpairs.bind_czrel.$$.list") ;
      foreach my $exporf12 (keys %{$allexp_sgdsigs}) {
         print EXPORFSIG "$exporf12\n" ;
      }
      close(EXPORFSIG) ;
      print STDERR "ALLEXP: TOTAL NUMBER OF BINARY INTERACTIONS (BIND V Cellzome-rel): $numkeys\n";

      print STDERR "   BIND: $numbind\n";

      print STDERR "   Cellzome-rel: $numczrel\n";
   }

   my $filters_allexpoverlap ;
   my $filters_bindoverlap ;
   my $filters_czoverlap ;
   my $filters_czreloverlap ;
   foreach my $ftype (keys %{$filters}) {
      foreach my $sgdsig (keys %{$filters->{$ftype}->{protprot}}) {
         my ($sgd1, $sgd2) = split(/\t/, $sgdsig) ;
         my $seqsig = join("\n",
            sort ($sgd_2_seqid->{$sgd1},$sgd_2_seqid->{$sgd2})) ;

         if (exists $bind_binint->{$sgdsig}  ||
             exists $cz_binint->{$sgdsig} &&
             $cz_binint->{$sgdsig} >= 5) {
            $filters_allexpoverlap->{$ftype}->{protprot}->{$sgdsig}++ ;
            $filters_allexpoverlap->{$ftype}->{prot}->{$sgd1}++ ;
            $filters_allexpoverlap->{$ftype}->{prot}->{$sgd2}++ ;
            foreach my $sgdsigdom (keys %{$filterset_z->{$sgdsig}}) {
               my ($sgd1, $dom1, $sgd2, $dom2) = split(/\t/, $sgdsigdom) ;
               $filters_allexpoverlap->{$ftype}->{domdom}->{$sgdsigdom}++ ;
               $filters_allexpoverlap->{$ftype}->{dom}->{$sgd1."\t".$dom1}++ ;
               $filters_allexpoverlap->{$ftype}->{dom}->{$sgd2."\t".$dom2}++ ;
            }
         }

         if (exists $bind_binint->{$sgdsig}) {
            $bin2comments->{$seqsig}->{'BIND'}++;

            $filters_bindoverlap->{$ftype}->{protprot}->{$sgdsig}++ ;
            $filters_bindoverlap->{$ftype}->{prot}->{$sgd1}++ ;
            $filters_bindoverlap->{$ftype}->{prot}->{$sgd2}++ ;
            foreach my $sgdsigdom (keys %{$filterset_z->{$sgdsig}}) {
               my ($sgd1, $dom1, $sgd2, $dom2) = split(/\t/, $sgdsigdom) ;
               $filters_bindoverlap->{$ftype}->{domdom}->{$sgdsigdom}++ ;
               $filters_bindoverlap->{$ftype}->{dom}->{$sgd1."\t".$dom1}++ ;
               $filters_bindoverlap->{$ftype}->{dom}->{$sgd2."\t".$dom2}++ ;
            }
         }

         if (exists $cz_binint->{$sgdsig}) {

            $filters_czoverlap->{$ftype}->{protprot}->{$sgdsig}++ ;
            $filters_czoverlap->{$ftype}->{prot}->{$sgd1}++ ;
            $filters_czoverlap->{$ftype}->{prot}->{$sgd2}++ ;
            foreach my $sgdsigdom (keys %{$filterset_z->{$sgdsig}}) {
               my ($sgd1, $dom1, $sgd2, $dom2) = split(/\t/, $sgdsigdom) ;
               $filters_czoverlap->{$ftype}->{domdom}->{$sgdsigdom}++ ;
               $filters_czoverlap->{$ftype}->{dom}->{$sgd1."\t".$dom1}++ ;
               $filters_czoverlap->{$ftype}->{dom}->{$sgd2."\t".$dom2}++ ;
            }
         }

         if (exists $cz_binint->{$sgdsig}&& $cz_binint->{$sgdsig} >= 5) {
            $bin2comments->{$seqsig}->{'CZ'}++;

            $filters_czreloverlap->{$ftype}->{protprot}->{$sgdsig}++ ;
            $filters_czreloverlap->{$ftype}->{prot}->{$sgd1}++ ;
            $filters_czreloverlap->{$ftype}->{prot}->{$sgd2}++ ;
            foreach my $sgdsigdom (keys %{$filterset_z->{$sgdsig}}) {
               my ($sgd1, $dom1, $sgd2, $dom2) = split(/\t/, $sgdsigdom) ;
               $filters_czreloverlap->{$ftype}->{domdom}->{$sgdsigdom}++ ;
               $filters_czreloverlap->{$ftype}->{dom}->{$sgd1."\t".$dom1}++ ;
               $filters_czreloverlap->{$ftype}->{dom}->{$sgd2."\t".$dom2}++ ;
            }
         }
      }
   }

   format_modbaseoutput({
      comments => {
         complexes => $cid2comments,
         binary => $bin2comments
      },
      complexes_fn => $in->{complexes_fn},
      binary_fn => $in->{binary_fn},
      project => 'yeastSGD1',
      out_fn => "modtie_complexes.$$.modtie",
   }) ;


   my $latexheader = {
      none => 'none',
      z => 'Z-score $\le$ '.$z_thresh,
      zl => 'Z + Co-Localization',
      zf => 'Z + Co-Function',
      zlf => 'Z + Co-Loc + Co-Func',
   } ;

   my $latexout ;
   foreach my $ftype (keys %{$filters}) {
      my $outstring = "$latexheader->{$ftype} " ;
      foreach my $key (sort {$b cmp $a} keys %{$filters->{$ftype}}) {
         my $tnum = keys %{$filters->{$ftype}->{$key}} ;
         $outstring .= " & ".commathrees($tnum) ;
         print STDERR "filter $ftype\t$key\t$tnum\n" ;
         if ($key eq 'protprot') {
            my $numppo_allexp = 0;
               $numppo_allexp += keys %{$filters_allexpoverlap->{$ftype}->{protprot}} ;
            my $numppo_bind = 0;
               $numppo_bind += keys %{$filters_bindoverlap->{$ftype}->{protprot}} ;
            my $numppo_cz = 0;
               $numppo_cz += keys %{$filters_czoverlap->{$ftype}->{protprot}} ;

            my $numppo_czrel = 0 ;
               $numppo_czrel += keys %{$filters_czreloverlap->{$ftype}->{protprot}} ;

            my $numpp = keys %{$filters->{$ftype}->{protprot}} ;

            my $ratioo_allexp = sprintf("%.2f", ($numppo_allexp / $numpp)) ;
            my $ratioo_bind = sprintf("%.2f", ($numppo_bind / $numpp)) ;
            my $ratioo_cz = sprintf("%.2f", ($numppo_cz / $numpp)) ;
            my $ratioo_czrel = sprintf("%.2f", ($numppo_czrel / $numpp)) ;

            print STDERR "ALLEXP protprot_overlap / protprot $ftype = $numppo_allexp ($ratioo_allexp)\n" ;
            print STDERR "BIND protprot_overlap / protprot $ftype = $numppo_bind ($ratioo_bind)\n" ;
            print STDERR "CZ protprot_overlap / protprot $ftype = $numppo_cz ($ratioo_cz)\n" ;
            print STDERR "CZrel protprot_overlap / protprot $ftype = $numppo_czrel ($ratioo_czrel)\n" ;
            print STDERR "\n" ;
            $outstring .= " ( $numppo_allexp, $numppo_bind, $numppo_cz, $numppo_czrel )"  ;
         }
      }
      $outstring .= ' \\\\' ;
      $latexout->{$ftype} = $outstring ;

      foreach my $key (sort {$b cmp $a} keys %{$filters_allexpoverlap->{$ftype}}){
         my $tnum = keys %{$filters_allexpoverlap->{$ftype}->{$key}} ;
         print STDERR "ALLEXP overlap filter $ftype\t$key\t$tnum\n" ;
      }
      foreach my $key (sort {$b cmp $a} keys %{$filters_bindoverlap->{$ftype}}){
         my $tnum = keys %{$filters_bindoverlap->{$ftype}->{$key}} ;
         print STDERR "BIND overlap filter $ftype\t$key\t$tnum\n" ;
      }
      foreach my $key (sort {$b cmp $a} keys %{$filters_czoverlap->{$ftype}}){
         my $tnum = keys %{$filters_czoverlap->{$ftype}->{$key}} ;
         print STDERR "CZ overlap filter $ftype\t$key\t$tnum\n" ;
      }
      foreach my $key (sort {$b cmp $a} keys %{$filters_czreloverlap->{$ftype}}){
         my $tnum = keys %{$filters_czreloverlap->{$ftype}->{$key}} ;
         print STDERR "CZrel overlap filter $ftype\t$key\t$tnum\n" ;
      }
      print STDERR "\n\n" ;
   }



   print STDERR "Fig 1a. filter lines:\n";
   print STDERR join("\n", ($latexout->{z}, $latexout->{zf},
                     $latexout->{zl}, $latexout->{zlf}))."\n" ;
   print STDERR "\n\n";

   my $sgdprop = readin_sgd_protein_properties() ;
#   build_yeast_graph({
#      filters => $filters,
#      bind_binint => $bind_binint,
#      genenames => $sgdprop->{orf_2_name}
#   }) ;

   print STDERR "Starting Plots\n" ;
   print STDERR "Plot 2. Total Target-Template Sequence Identity vs z-score\n" ;

   my $seqid_vs_z ;
   foreach my $j ( 0 .. $#{$rawpred}) {
      if ($rawpred->[$j]->[$head->{mt_bin}->{f2i}->{z}] > -1) {next;}

      my $seqid1= $rawpred->[$j]->[$head->{mt_bin}->{f2i}->{seq1}];
      my $seqid2= $rawpred->[$j]->[$head->{mt_bin}->{f2i}->{seq2}];
      my $sgd1 = $seqid_2_sgd->{$seqid1} ;
      my $sgd2 = $seqid_2_sgd->{$seqid2} ;
      my $sgdsig = $sgd1."\t".$sgd2 ;
      if ($sgd2 lt $sgd1) {
         $sgdsig = $sgd2."\t".$sgd1 ; }

      my $inttype = 0 ;
      if (exists $bind_binint->{$sgdsig}) {$inttype = 1 ;}
      elsif (exists $filters->{zlf}->{protprot}->{$sgdsig}) {$inttype = 2 ;}
      elsif (exists $filters->{zl}->{protprot}->{$sgdsig}) {$inttype = 3 ;}
      elsif (exists $filters->{zf}->{protprot}->{$sgdsig}) {$inttype = 4 ;}

      my $avg_seqid = sprintf("%.2f", (($rawpred->[$j]->[$head->{mt_bin}->{f2i}->{numres_ident_dom1}] / $rawpred->[$j]->[$head->{mt_bin}->{f2i}->{numres_targ_dom1}]) + ($rawpred->[$j]->[$head->{mt_bin}->{f2i}->{numres_ident_dom2}] / $rawpred->[$j]->[$head->{mt_bin}->{f2i}->{numres_targ_dom2}])) * 50) ;

      if ($inttype == 1 &&
          (($rawpred->[$j]->[$head->{mt_bin}->{f2i}->{z}] <= -3 &&
            $avg_seqid <= 30)  ||
           ($rawpred->[$j]->[$head->{mt_bin}->{f2i}->{z}] >= -2 &&
            $avg_seqid >= 80)) ) {
          my @outvals = ($seqid1,$sgd1,$seqid2,$sgd2,
                         $rawpred->[$j]->[$head->{mt_bin}->{f2i}->{z}],
                         $avg_seqid, 'OVERLAP', $mtint_maxseqid->{$sgdsig},
                         $modtie_binint->{$sgdsig}) ;

         my $tsid1 = $rawpred->[$j]->[$head->{mt_bin}->{f2i}->{subset_id_1}] ;
         my $tsid2 = $rawpred->[$j]->[$head->{mt_bin}->{f2i}->{subset_id_2}] ;
         my $tsid12 = $tsid1."\n".$tsid2 ;
         if ($tsid2 lt $tsid1) {$tsid12 = $tsid2."\n".$tsid1;}
         push @outvals, $pb->{sid12chains}->{$tsid12} ;
         print join("\t", @outvals)."\n" ;
      }


      push @{$seqid_vs_z->[$inttype]->{x}}, $rawpred->[$j]->[$head->{mt_bin}->{f2i}->{z}] ;
      push @{$seqid_vs_z->[$inttype]->{y}}, $avg_seqid ;
   }

   my $seqid_vs_bz ;
#   foreach my $sgdsig ( keys %{$filters->{none}->{protprot}})
   open (ALLMTFH, ">seqid_vs_bestz.zlf.txt") ;
   open (PREDONLYFH, ">seqid_vs_bestz.predonly.zlf.txt") ;
   open (OVERLAPFH, ">seqid_vs_bestz.overlap.zlf.txt") ;
   foreach my $sgdsig ( keys %{$filters->{zlf}->{protprot}}) {
      if ($modtie_binint->{$sgdsig} > -1) {next;}

      my $j = $mtint_bestzentry->{$sgdsig} ;

      my ($sgd1, $sgd2) = split(/\t/, $sgdsig) ;

      my $seqid1 = $sgd_2_seqid->{$sgd1} ;
      my $seqid2 = $sgd_2_seqid->{$sgd2} ;

      my $inttype = 0 ;
      if (exists $bind_binint->{$sgdsig} ||
          exists $cz_binint->{$sgdsig} &&
          $cz_binint->{$sgdsig} >= 5 ) {$inttype = 1 ;}

      my $avg_seqid = sprintf("%.2f", (($rawpred->[$j]->[$head->{mt_bin}->{f2i}->{numres_ident_dom1}] / $rawpred->[$j]->[$head->{mt_bin}->{f2i}->{numres_targ_dom1}]) + ($rawpred->[$j]->[$head->{mt_bin}->{f2i}->{numres_ident_dom2}] / $rawpred->[$j]->[$head->{mt_bin}->{f2i}->{numres_targ_dom2}])) * 50) ;

      if ($inttype == 1 &&
          (($rawpred->[$j]->[$head->{mt_bin}->{f2i}->{z}] <= -2 &&
            $avg_seqid <= 40)  ||
           ($rawpred->[$j]->[$head->{mt_bin}->{f2i}->{z}] > -2 &&
            $avg_seqid >= 80)) ) {

          my @outvals = ($seqid1,$sgd1,$seqid2,$sgd2,
                         $rawpred->[$j]->[$head->{mt_bin}->{f2i}->{z}],
                         $avg_seqid, 'OVERLAP',
                         $mtint_maxseqid->{$sgdsig},
                         $modtie_binint->{$sgdsig},
                         $rawpred->[$j]->[$head->{mt_bin}->{f2i}->{raw}],
                         $sgdprop->{orf_2_name}->{$sgd1},
                         $sgdprop->{orf_2_name}->{$sgd2},
         ) ;
         my $tsid1 = $rawpred->[$j]->[$head->{mt_bin}->{f2i}->{subset_id_1}] ;
         my $tsid2 = $rawpred->[$j]->[$head->{mt_bin}->{f2i}->{subset_id_2}] ;
         my $tsid12 = $tsid1."\n".$tsid2 ;
         if ($tsid2 lt $tsid1) {
            $tsid12 = $tsid2."\n".$tsid1;}

         push @outvals, $pb->{sid12chains}->{$tsid12} ;
         print join("\t", @outvals)."\n" ;
      }

      print ALLMTFH "$rawpred->[$j]->[$head->{mt_bin}->{f2i}->{z}] $avg_seqid".
                    "\n";

      if ($inttype == 0) {
         print PREDONLYFH "$rawpred->[$j]->[$head->{mt_bin}->{f2i}->{z}] ".
                          "$avg_seqid\n" ;
      } elsif ($inttype == 1) {
         print OVERLAPFH "$rawpred->[$j]->[$head->{mt_bin}->{f2i}->{z}] ".
                         "$avg_seqid\n" ; }
   }
   close(ALLMTFH) ;
   close(OVERLAPFH) ;
   close(PREDONLYFH) ;

   print STDERR "Fig 1b: Average seqid vs (best) z-score (for the interaction)\n\n";
   print STDERR "seqidbestz<-read.table('seqid_vs_bestz.zlf.txt', sep=' ')\n";
   print STDERR "postscript('seqid_vs_bestz.zlf.eps',horizontal=FALSE,onefile=FALSE,pointsize=16, height=8,width=8)\n";
   print STDERR "par(mar=c(4.1, 4.1, 0.5, 0.5))\n";
   print STDERR "plot(seqidbestz\$V1,seqidbestz\$V2, 'p', xlab='Z-score', ylab='Target - Template Sequence Identity', cex=0.5, cex.lab=1.25, cex.axis=1.25)\n";
   print STDERR "rug(jitter(seqidbestz\$V1), side=1,ticksize=0.015)\n";
   print STDERR "rug(jitter(seqidbestz\$V2), side=2,ticksize=0.015)\n";
   print STDERR "dev.off()\n\n\n";

   print STDERR "Fig 1c: P(overlap | z <= thresh) for z, zl, zf, and zlf\n";
   foreach my $filter (qw/z zl zf zlf/) {
      print STDERR "  $filter: " ;
      my $overlap = 0 ;
      my $all = 0 ;
      my $points ;
      foreach my $sgdsig (sort {$modtie_binint->{$a} <=> $modtie_binint->{$b}}
                           keys %{$filters->{$filter}->{protprot}}) {
         if (exists $bind_binint->{$sgdsig}) {
            $overlap++ ;
         } else {
            $all++ ;
         }
         $points->{$filter}->{$modtie_binint->{$sgdsig}} =
            sprintf("%.3f", ($overlap / $all)) ;
      }

      open(OUTFH, ">enrichment_filter_bind.$filter.txt") ;
      foreach my $thresh (sort {$a <=> $b} keys %{$points->{$filter}}) {
         print OUTFH "$thresh $points->{$filter}->{$thresh}\n" ;
      }
      close(OUTFH) ;
      print STDERR "X\n" ;
   }

   print STDERR "Fig 1c: Enrichment of Filters\n\n";
   print STDERR 'postscript("enrichment_filter_bind.eps",horizontal=FALSE,onefile=FALSE,pointsize=16, height=8,width=8)'."\n";
   print STDERR 'par(mar=c(4.1, 4.1, 0.5, 0.5))'."\n";
   print STDERR 'zlf<-read.table("enrichment_filter_bind.zlf.txt", sep=" ")'."\n";
   print STDERR 'zf<-read.table("enrichment_filter_bind.zf.txt", sep=" ")'."\n";
   print STDERR 'zl<-read.table("enrichment_filter_bind.zl.txt", sep=" ")'."\n";
   print STDERR 'z<-read.table("enrichment_filter_bind.z.txt", sep=" ")'."\n";
   print STDERR 'plot(zlf$V1,zlf$V2, "l", lty=1, xlab="Z-score Threshold",  ylab=expression(paste("p(Overlap | ", Z-score <= Threshold, ")")), cex=0.5, cex.lab=1.25, cex.axis=1.25, xlim=c(-5,-2))'."\n";
   print STDERR 'lines(zl$V1,zl$V2, lty=2)'."\n";
   print STDERR 'lines(zf$V1,zf$V2, lty=3)'."\n";
   print STDERR 'lines(z$V1,z$V2, lty=4)'."\n";
   print STDERR 'legend("bottomleft",legend=c("Z-score", "Z + Co-Function", "Z + Co-Localization", "Z + Co-Loc + Co-Func"), lty= c(4,3,2,1), bty="n")'."\n";
   print STDERR "dev.off()\n\n\n";



   print STDERR "Fig 1d: P(overlap | z <= thresh) for z, zl, zf, and zlf\n";
   foreach my $filter (qw/z zl zf zlf/) {
      print STDERR "  $filter: " ;
      my $overlap = 0 ;
      my $all = 0 ;
      my $points ;
      foreach my $sgdsig (sort {$modtie_binint->{$a} <=> $modtie_binint->{$b}}
                           keys %{$filters->{$filter}->{protprot}}) {
         if (exists $bind_binint->{$sgdsig} ||
             exists $cz_binint->{$sgdsig} &&
             $cz_binint->{$sgdsig} >= 5) {
            $overlap++ ;
         } else {
            $all++ ;
         }
         $points->{$filter}->{$modtie_binint->{$sgdsig}} =
            sprintf("%.3f", ($overlap / $all)) ;
      }

      open(OUTFH, ">enrichment_filter_bindorczrel.$filter.txt") ;
      foreach my $thresh (sort {$a <=> $b} keys %{$points->{$filter}}) {
         print OUTFH "$thresh $points->{$filter}->{$thresh}\n" ;
      }
      close(OUTFH) ;
      print STDERR "X\n" ;
   }


   print STDERR "Fig 1c: Enrichment of Filters\n\n";
   print STDERR 'postscript("enrichment_filter_bindorczrel.eps",horizontal=FALSE,onefile=FALSE,pointsize=16, height=8,width=8)'."\n";
   print STDERR 'par(mar=c(4.1, 4.1, 0.5, 0.5))'."\n";
   print STDERR 'zlf<-read.table("enrichment_filter_bindorczrel.zlf.txt", sep=" ")'."\n";
   print STDERR 'zf<-read.table("enrichment_filter_bindorczrel.zf.txt", sep=" ")'."\n";
   print STDERR 'zl<-read.table("enrichment_filter_bindorczrel.zl.txt", sep=" ")'."\n";
   print STDERR 'z<-read.table("enrichment_filter_bindorczrel.z.txt", sep=" ")'."\n";
   print STDERR 'plot(zlf$V1,zlf$V2, "l", lty=1, xlab="Z-score Threshold",  ylab=expression(paste("p(Overlap | ", Z-score <= Threshold, ")")), cex=0.5, cex.lab=1.25, cex.axis=1.25, xlim=c(-5,-2))'."\n";
   print STDERR 'lines(zl$V1,zl$V2, lty=2)'."\n";
   print STDERR 'lines(zf$V1,zf$V2, lty=3)'."\n";
   print STDERR 'lines(z$V1,z$V2, lty=4)'."\n";
   print STDERR 'legend("bottomleft",legend=c("Z-score", "Z + Co-Function", "Z + Co-Localization", "Z + Co-Loc + Co-Func"), lty= c(4,3,2,1), bty="n")'."\n";
   print STDERR "dev.off()\n\n\n";


   print STDERR "Fig 1c: P(overlap | z <= thresh) for z, zl, zf, and zlf\n";
   foreach my $filter (qw/z zl zf zlf/) {
      print STDERR "  $filter: " ;
      my $overlap = 0 ;
      my $all = 0 ;
      my $points ;
      foreach my $sgdsig (sort {$modtie_binint->{$a} <=> $modtie_binint->{$b}}
                           keys %{$filters->{$filter}->{protprot}}) {
         if (exists $cz_binint->{$sgdsig}) {
            $overlap++ ;
         } else {
            $all++ ;
         }
         $points->{$filter}->{$modtie_binint->{$sgdsig}} =
            sprintf("%.3f", ($overlap / $all)) ;
      }

      open(OUTFH, ">enrichment_filter_cz.$filter.txt") ;
      foreach my $thresh (sort {$a <=> $b} keys %{$points->{$filter}}) {
         print OUTFH "$thresh $points->{$filter}->{$thresh}\n" ;
      }
      close(OUTFH) ;
      print STDERR "X\n" ;
   }

   print STDERR "Fig 1c: Enrichment of Filters\n\n";
   print STDERR 'postscript("enrichment_filter_cz.eps",horizontal=FALSE,onefile=FALSE,pointsize=16, height=8,width=8)'."\n";
   print STDERR 'par(mar=c(4.1, 4.1, 0.5, 0.5))'."\n";
   print STDERR 'zlf<-read.table("enrichment_filter_cz.zlf.txt", sep=" ")'."\n";
   print STDERR 'zf<-read.table("enrichment_filter_cz.zf.txt", sep=" ")'."\n";
   print STDERR 'zl<-read.table("enrichment_filter_cz.zl.txt", sep=" ")'."\n";
   print STDERR 'z<-read.table("enrichment_filter_cz.z.txt", sep=" ")'."\n";
   print STDERR 'plot(zlf$V1,zlf$V2, "l", lty=1, xlab="Z-score Threshold",  ylab=expression(paste("p(Overlap | ", Z-score <= Threshold, ")")), cex=0.5, cex.lab=1.25, cex.axis=1.25, xlim=c(-5,-2))'."\n";
   print STDERR 'lines(zl$V1,zl$V2, lty=2)'."\n";
   print STDERR 'lines(zf$V1,zf$V2, lty=3)'."\n";
   print STDERR 'lines(z$V1,z$V2, lty=4)'."\n";
   print STDERR 'legend("bottomleft",legend=c("Z-score", "Z + Co-Function", "Z + Co-Localization", "Z + Co-Loc + Co-Func"), lty= c(4,3,2,1), bty="n")'."\n";
   print STDERR "dev.off()\n\n\n";


   print STDERR "Fig 1c: P(overlap | z <= thresh) for z, zl, zf, and zlf\n";
   foreach my $filter (qw/z zl zf zlf/) {
      print STDERR "  $filter: " ;
      my $overlap = 0 ;
      my $all = 0 ;
      my $points ;
      foreach my $sgdsig (sort {$modtie_binint->{$a} <=> $modtie_binint->{$b}}
                           keys %{$filters->{$filter}->{protprot}}) {
         if (exists $cz_binint->{$sgdsig} && ($cz_binint->{$sgdsig} >= 5)) {
            $overlap++ ;
         } else {
            $all++ ;
         }
         $points->{$filter}->{$modtie_binint->{$sgdsig}} =
            sprintf("%.3f", ($overlap / $all)) ;
      }

      open(OUTFH, ">enrichment_filter_czrel.$filter.txt") ;
      foreach my $thresh (sort {$a <=> $b} keys %{$points->{$filter}}) {
         print OUTFH "$thresh $points->{$filter}->{$thresh}\n" ;
      }
      close(OUTFH) ;
      print STDERR "X\n" ;
   }

   print STDERR "Fig 1c: Enrichment of Filters\n\n";
   print STDERR 'postscript("enrichment_filter_czrel.eps",horizontal=FALSE,onefile=FALSE,pointsize=16, height=8,width=8)'."\n";
   print STDERR 'par(mar=c(4.1, 4.1, 0.5, 0.5))'."\n";
   print STDERR 'zlf<-read.table("enrichment_filter_czrel.zlf.txt", sep=" ")'."\n";
   print STDERR 'zf<-read.table("enrichment_filter_czrel.zf.txt", sep=" ")'."\n";
   print STDERR 'zl<-read.table("enrichment_filter_czrel.zl.txt", sep=" ")'."\n";
   print STDERR 'z<-read.table("enrichment_filter_czrel.z.txt", sep=" ")'."\n";
   print STDERR 'plot(zlf$V1,zlf$V2, "l", lty=1, xlab="Z-score Threshold",  ylab=expression(paste("p(Overlap | ", Z-score <= Threshold, ")")), cex=0.5, cex.lab=1.25, cex.axis=1.25, xlim=c(-5,-2))'."\n";
   print STDERR 'lines(zl$V1,zl$V2, lty=2)'."\n";
   print STDERR 'lines(zf$V1,zf$V2, lty=3)'."\n";
   print STDERR 'lines(z$V1,z$V2, lty=4)'."\n";
   print STDERR 'legend("bottomleft",legend=c("Z-score", "Z + Co-Function", "Z + Co-Localization", "Z + Co-Loc + Co-Func"), lty= c(4,3,2,1), bty="n")'."\n";
   print STDERR "dev.off()\n\n\n";


}



sub readin_sgd_go {

   my $in = shift ;
   my $specs ;
   if (exists $in->{specs}) { $specs = $in->{specs} ; }
   else { $specs = modtie::set_modtie_specs() ; }

   my $sgdprop = readin_sgd_protein_properties() ;

   my $orf2go ;
   my $go2orf ;
   my $goterm2type ;
   open(FILE, $specs->{sgd_go_fn}) ;
   while (my $line = <FILE>) {
      if ($line =~ /\!/) {next;}
      chomp $line;
      my @t = split(/\t/, $line) ;
      my $sgdid = $t[1] ;
      my $goterm = $t[4] ;
      my $type = $t[8] ;


      my $orf = $sgdprop->{sgdid_2_orf}->{$sgdid} ;
      if (!defined $orf) {
#         print STDERR "WARNING: no orf id for $sgdid\n";
         next;}

      $goterm2type->{$goterm} =$type ;
      $orf2go->{$orf}->{$goterm}++;
      $go2orf->{$goterm}->{$orf}++;
   }

   return {
      orf2go => $orf2go,
      go2orf => $go2orf,
      go2type => $goterm2type
   } ;

}



sub build_yeast_graph {

   require LGL ;
#   use GD ;

   my $in = shift ;
   my $filters = $in->{filters} ;
   my $bind_binint= $in->{bind_binint} ;
   my $genenames = $in->{genenames} ;

   my $sgdurlprefix = 'http://db.yeastgenome.org/cgi-bin/locus.pl?locus=';
#   my $overlap_binint ;

   my ($urls, $edges, $edges_col, $nodes_col) ;
   foreach my $sgdsig (keys %{$filters->{z}->{protprot}}) {
      my ($sgd1, $sgd2) = split(/\t/, $sgdsig) ;
      if ($sgd1 eq $sgd2) {next;}

      my $sig1 =$sgd1; my $sig2 = $sgd2 ;
      if ($sgd2 lt $sgd1) { $sig2 =$sgd1; $sig1 = $sgd2 ;}

      if ((!exists $filters->{zlf}->{protprot}->{$sgdsig} &&
          !exists $filters->{zf}->{protprot}->{$sgdsig} &&
          !exists $filters->{zl}->{protprot}->{$sgdsig}) &&
          !exists $bind_binint->{$sgdsig}) {next;}

#      my $seq1 = $seqidsgd->{sgd_2_seqid}->{$sgd1} ;
#      my $seq2 = $seqidsgd->{sgd_2_seqid}->{$sgd2} ;

      $edges->{$sgd1}->{$sgd2}++ ;
      $nodes_col->{$sgd1} = 'LightGray' ;
      $nodes_col->{$sgd2} = 'LightGray' ;
      $urls->{$sgd1} = $sgdurlprefix.$sgd1 ;
      $urls->{$sgd2} = $sgdurlprefix.$sgd2 ;

      if (exists $bind_binint->{$sgdsig}) {
         $edges_col->{$sgd1}->{$sgd2} = 'green' ;
      } elsif (exists $filters->{zlf}->{protprot}->{$sgdsig}) {
         $edges_col->{$sgd1}->{$sgd2} = 'red' ;
      } elsif (exists $filters->{zl}->{protprot}->{$sgdsig}) {
         $edges_col->{$sgd1}->{$sgd2} = 'purple' ;
      } elsif (exists $filters->{zf}->{protprot}->{$sgdsig}) {
         $edges_col->{$sgd1}->{$sgd2} = 'orange' ;
      }
   }


#   foreach my $bind_int (keys %{$bind_binint}) {
#      my ($sgd1, $sgd2) = split(/\t/, $bind_int) ;
#      $yeasturls->{$sgd1} = $sgdurlprefix.$sgd1 ;
#      $yeasturls->{$sgd2} = $sgdurlprefix.$sgd2 ;
#      my $seq1 = $seqidsgd->{sgd_2_seqid}->{$sgd1} ;
#      my $seq2 = $seqidsgd->{sgd_2_seqid}->{$sgd2} ;
#      $mt_lglnodes_col->{$sgd1} = 'LightGray' ;
#      $mt_lglnodes_col->{$sgd2} = 'LightGray' ;
#
#      my $coloc = 'U';
#      if (exists  $yeastgfp->{orf2lochash}->{$sgd1} &&
#          exists  $yeastgfp->{orf2lochash}->{$sgd2}) {
#         $coloc = 0 ;
#         foreach my $l (keys %{$yeastgfp->{orf2lochash}->{$sgd1}}) {
#            if (exists $yeastgfp->{orf2lochash}->{$sgd2}->{$l}) {
#               $coloc++; next;}}}
#
#      my $cofunc ;
#      $cofunc->{'P'} = 'U' ;
#      $cofunc->{'C'} = 'U' ;
#      $cofunc->{'F'} = 'U' ;
#      $cofunc->{'all'} = 'U' ;
#      if (exists $sgdgo->{orf2go}->{$sgd1} &&
#          exists $sgdgo->{orf2go}->{$sgd2}) {
#         $cofunc->{'P'} = 0 ;
#         $cofunc->{'C'} = 0 ;
#         $cofunc->{'F'} = 0 ;
#         $cofunc->{'all'} = 0 ;
#         foreach my $goterm (keys %{$sgdgo->{orf2go}->{$sgd1}}) {
#            if (exists $sgdgo->{orf2go}->{$sgd2}->{$goterm}) {
#               $cofunc->{'all'}++ ;
#               $cofunc->{$sgdgo->{go2type}->{$goterm}}++ ;
#               next;
#            }
#         }
#      }
#
#      my $div= join(', ', sort keys %{$bind_div->{rgid2div}->{$bind_binint->{$bind_int}}});
#
#      my @outvals ;
#
#      if (exists $modtie_binint->{$bind_int}) {
#         @outvals = ('OVERLAP', $bind_int, $seq1, $seq2,
#                     $coloc, $cofunc->{'all'}, $cofunc->{'C'},
#                     $cofunc->{'F'}, $cofunc->{'P'}, $div) ;
#         my $sig1 =$sgd1; my $sig2 = $sgd2 ;
#         if ($sgd2 lt $sgd1) {
#            $sig2 =$sgd1; $sig1 = $sgd2 ;}
#         $mt_lgledges_col->{$sig1}->{$sig2} = 'green' ;
#      } else {
#         @outvals = ('BINDONLY', $bind_int, $seq1, $seq2,
#                     $coloc, $cofunc->{'all'}, $cofunc->{'C'},
#                     $cofunc->{'F'}, $cofunc->{'P'}, $div) ;
#      }
#      print join("\t", @outvals)."\n" ;
#   }
#

   print STDERR "Laying out Binary predictions: ";
   my $coords = LGL::edges2coords({edges => $edges}) ;

##   my ($graph_fh, $graph_fn) = tempfile("mt_yeast.XXXXX", SUFFIX => ".lgl.ps") ;
##   LGL::coords2ps({
##      fh => $graph_fh,
##      coords => $coords,
##      edges => $mt_lgledges,
##      ecols => $mt_lgledges_col,
##   }) ;

   my ($graphpng_fh, $graphpng_fn) =
      tempfile("mt_yeast_graph.XXXXX", SUFFIX => ".png") ;
   my ($graphmap_fh, $graphmap_fn) =
      tempfile("mt_yeast_graph.XXXXX", SUFFIX => ".html") ;
   my $im = LGL::coords2png({
      mag => 75,
      nrad => 3,
      ethick => 1.5,
      coords => $coords,
      edges => $edges,
      ecols => $edges_col,
      ncols => $nodes_col,
   }) ;

   binmode($graphpng_fh) ;
   print {$graphpng_fh} $im->png() ;
   close($graphpng_fh) ;

   LGL::coords2pngmap({
      mag => 75,
      nodeurls => $urls,
      nodenames => $genenames, 
      nrad => 5,
      ethick => 2,
      img_fn => $graphpng_fn,
      map_fh => $graphmap_fh,
      map_name => 'modtieyeast',
      coords => $coords,
      edges => $edges,
      ecols => $edges_col,
   }) ;
   close($graphmap_fh) ;
   my $dir = getcwd ;
   print STDERR "file://".$dir."/".$graphmap_fn."\n" ;

   print STDERR " X\n" ;
}

sub readin_yeastgfp {

   my $in = shift ; my $specs ;
   if (exists $in->{specs}) { $specs = $in->{specs} ; }
   else { $specs = modtie::set_modtie_specs() ; }

   my $fn ; 
   if (exists $in->{fn}) {
      $fn = $in->{fn} ;
   } else {
      $fn = $specs->{yeastgfp_fn} ;
   }

   open(YEASTGFP, $fn) ;
   my $field2no ;
   my $fields;
   my $orf2loc ;
   my $orf2lochash ;
   my $linenum = 0 ;
   while (my $line = <YEASTGFP>) {
      chomp $line ;
      my @f = split(/\t/, $line) ;
      if ($linenum == 0) {
         foreach my $j ( 0 .. $#f) {
            $fields->[$j] = $f[$j] ;
            $field2no->{$f[$j]} = $j ;
         }
      } else {
         if (defined $f[$field2no->{'localization summary'}] &&
            $f[$field2no->{'localization summary'}] ne '') {

            $orf2loc->{$f[$field2no->{'yORF'}]} =
               $f[$field2no->{'localization summary'}] ;
#            print "set ".$f[$field2no->{'yORF'}]." to ".
#                     $f[$field2no->{'localization summary'}]."\n" ;

            my @locs = split(/\,/, $f[$field2no->{'localization summary'}]) ;

            my $hloc ;
            foreach my $tloc (@locs) {
               if ($tloc eq 'ambiguous') {next;}
               $orf2lochash->{$f[$field2no->{'yORF'}]}->{$tloc}++ ;
            }
         }
      }
      $linenum++ ;
   }
   close(YEASTGFP) ;

   my $data = {
      orf2loc => $orf2loc,
      orf2lochash => $orf2lochash
   } ;

   return $data ;

}


sub get_yeastnames_seqid_sgd {

   my $dbh = connect_modbase() ;

   my ($seqid2sgd) = pibase::mysql_hashload($dbh,
      "SELECT seq_id, database_id FROM sequences where run = \"yeast-SGD1\"") ;

   my ($sgd2seqid) = pibase::mysql_hashload($dbh,
      "SELECT database_id, seq_id FROM sequences where run = \"yeast-SGD1\"") ;

   return {
      seqid_2_sgd => $seqid2sgd,
      sgd_2_seqid => $sgd2seqid
   } ;

}


sub load_bind_yeast {

   my $bind_complex = load_bind_complex2subunits() ;
   my $bind_nrints = load_bind_nrints() ;

   my $seqidsgd = load_seqid_2_sgd() ;
   my $seqid_2_sgd = $seqidsgd->{seqid_2_sgd} ;

   my $genbanks ;
   foreach my $rgid ( keys %{$bind_nrints->{tax2rgid}->{"4932\n4932"}}) {
      $genbanks->{$bind_nrints->{nrints}->{$rgid}->{a_id}}++ ;
      $genbanks->{$bind_nrints->{nrints}->{$rgid}->{b_id}}++ ;
   }

   foreach my $bcid ( keys %{$bind_complex->{tax2bcid}->{4932}}) {
      foreach my $subunit ( @{$bind_complex->{complex2subunits}->{$bcid}}) {
         $genbanks->{$subunit->{id}}++ ;
      }
   }

   my ($dbh_mb) = connect_modbase()  ;

   my $gb_2_sgd ;
   foreach my $gid (keys %{$genbanks}) {
      my ($seqid) = pibase::mysql_singleval($dbh_mb,
      "SELECT seq_id FROM synonyms where database_id = \"$gid\"") ;

      if (defined $seqid && exists $seqid_2_sgd->{$seqid}) {
         $gb_2_sgd->{$gid} = $seqid_2_sgd->{$seqid};}
#      else {
#         print STDERR "WARNING: couldnt translate $gid\n" ;}
   }

   foreach my $rgid ( keys %{$bind_nrints->{tax2rgid}->{"4932\n4932"}}) {
      if (exists $gb_2_sgd->{$bind_nrints->{nrints}->{$rgid}->{a_id}}) {
         $bind_nrints->{nrints}->{$rgid}->{a_sgd} =
            $gb_2_sgd->{$bind_nrints->{nrints}->{$rgid}->{a_id}} ; }

      if (exists $gb_2_sgd->{$bind_nrints->{nrints}->{$rgid}->{b_id}}) {
         $bind_nrints->{nrints}->{$rgid}->{b_sgd} =
            $gb_2_sgd->{$bind_nrints->{nrints}->{$rgid}->{b_id}} ; }

      if (exists $bind_nrints->{nrints}->{$rgid}->{a_sgd} &&
          exists $bind_nrints->{nrints}->{$rgid}->{b_sgd}) {
         $bind_nrints->{p2p}->{$bind_nrints->{nrints}->{$rgid}->{a_sgd}}->{$bind_nrints->{nrints}->{$rgid}->{b_sgd}}++ ;
         $bind_nrints->{p2p}->{$bind_nrints->{nrints}->{$rgid}->{b_sgd}}->{$bind_nrints->{nrints}->{$rgid}->{a_sgd}}++ ;
      }
   }

   foreach my $bcid ( keys %{$bind_complex->{tax2bcid}->{4932}}) {
      foreach my $subunit ( @{$bind_complex->{complex2subunits}->{$bcid}}) {
         if (exists $gb_2_sgd->{$subunit->{id}}) {
            $subunit->{sgd} = $gb_2_sgd->{$subunit->{id}} ;
            $bind_complex->{sgd2bcid}->{$gb_2_sgd->{$subunit->{id}}}->{$bcid}++ ;
         }
      }
   }

   return {
      complex => $bind_complex,
      nrints =>  $bind_nrints,
   } ;

}

sub load_cellzome {

   my $binary;
   my $complex ;

   my $in = shift ; my $specs ;
   if (exists $in->{specs}) { $specs = $in->{specs} ; }
   else { $specs = modtie::set_modtie_specs() ; }

   open(BINF, $specs->{cellzome_binary_fn}) ;
   my $name2orf ;
   while (my $line = <BINF>) {
      chomp $line;

      if ($line =~ /^\#/) {next;}
      my (undef, $orf1, $name1, $orf2, $name2, $score) = split(' ', $line) ;
      my $orfsig = $orf1."\t".$orf2 ;
      if ($orf2 lt $orf1) { $orfsig = $orf2."\t".$orf1 ; }
      $binary->{$orfsig} = $score ; $binary->{$orfsig} = $score ;
      $name2orf->{$name1} = $orf1 ;
      $name2orf->{$name2} = $orf2 ;
   }
   close(BINF) ;

   open(COMPLF, $specs->{cellzome_complexes_fn}) ;
   while (my $line = <COMPLF>) {
      chomp $line;

      if ($line =~ /^\#/) {next;}
      $line =~ s/\"//g ;
      my ($cid, $name, $coreunits, $modules, $attach,$loc) = split(/\t/, $line);
      my @subunits = split(' ', $coreunits) ;

      $modules =~ s/\&\&//g ; $modules =~ s/[0-9]+\: //g ;
      push @subunits, split(' ', $modules) ;
      push @subunits, split(' ', $attach) ;
      my $subunits;
      foreach my $j ( 0 .. $#subunits) { $subunits->{$subunits[$j]}++; }
      foreach my $subunit ( keys %{$subunits}) {
         my $ucname = uc($subunit) ;
         my $tname = $ucname ;
         if (exists $name2orf->{$ucname}) { $tname = $name2orf->{$ucname};}
         else {
            print STDERR "  warning: cellzome complex $cid, no name: $ucname\n";
         }
         push @{$complex->{complex2subunits}->{$cid}},
            { sgd => $tname } ;
         $complex->{sgd2bcid}->{$tname}->{$cid}++ ;
      }
   }
   close(COMPLF);

   return ({
      binary => $binary,
      complex => $complex,
   }) ;

}


sub load_mips_complexes {

   my $in = shift ; my $specs ;
   if (exists $in->{specs}) { $specs = $in->{specs} ; }
   else { $specs = modtie::set_modtie_specs() ; }

   my $mipscompl ;
   open(FILE, $specs->{mips_complexcat_fn}) ;
   while (my $line = <FILE>) {
      if ($line =~ /\!/) {next;}
      chomp $line;
      my ($id, $cat, $evid,$ref) = split(/\|/, $line) ;
      $mipscompl->{id2cat}->{$id}->{$cat}++ ;
      $mipscompl->{cat2id}->{$cat}->{$id}++ ;
   }

   my $filtcompl ;
   foreach my $cat (keys %{$mipscompl->{cat2id}}) {
      my $numunits = keys %{$mipscompl->{cat2id}->{$cat}} ;
      if ($numunits > 1) {
         $filtcompl->{cat2id}->{$cat} = $mipscompl->{cat2id}->{$cat} ;
         foreach my $id (keys %{$filtcompl->{cat2id}->{$cat}}) {
            $filtcompl->{id2cat}->{$id}->{$cat}++ ; } } }

   return $filtcompl ;

}


sub calc_p_b_prior {

   my $in = shift ;

   my $domdata ;
   $domdata->{seqinfo} = {} ;


   my ($dbh) = pibase::connect_pibase() ;
   my ($idom1, $idom2) = pibase::mysql_fetchcols($dbh,
      "SELECT class_1, class_2 FROM intersubset_contacts WHERE subset_id_1 LIKE \"\%SCOP\%\" and num_contacts >= 1000 and chains = \"diff\" and class_1 != class_2") ;

   my $intfa12 ;

   my $intsf12 ;
   foreach my $j ( 0 .. $#{$idom1}) {
      $intfa12->{$idom1->[$j]."\n".$idom2->[$j]}++ ;
      $intfa12->{$idom2->[$j]."\n".$idom1->[$j]}++ ;

      my @p = split(/\./, $idom1->[$j]) ;
      my $sf1 = $p[0].'.'.$p[1].'.'.$p[2] ;

      @p = split(/\./, $idom2->[$j]) ;
      my $sf2 = $p[0].'.'.$p[1].'.'.$p[2] ;
      $intsf12->{$sf1."\n".$sf2}++ ;
      $intsf12->{$sf2."\n".$sf1}++ ;
   }


   my $seqid_2_domains = readin_seqid_aaseq_domarch({
      seqinfo => $domdata->{seqinfo},
      fn => $in->{seqid_aaseq_arch}}) ;

   my $mbinfo = _parse_domarchs($domdata) ;

   my ( $scop2seq, $seq2dom, $seq2domfrag, $scop2num ) =
      ($mbinfo->{scop2seq},
       $mbinfo->{seq2dom},
       $mbinfo->{seq2domfrag},
       $mbinfo->{scop2num}) ;

   my $aaseq = $domdata->{seqinfo}->{aaseq} ;


   my $yeastgfp = readin_yeastgfp() ;

   my $names = get_yeastnames_seqid_sgd() ;

   my $counts ;
   $counts->{coifam} = 0 ;
   $counts->{coloc} = 0 ;
   $counts->{coifam_coloc} = 0 ;

   my $seq2hloc ;
   my $seq2hdom_sf ;
   my $seq2hdom_fa ;

   foreach my $seq1 (sort keys %{$mbinfo->{seq2dom}}) {
      if (!exists $names->{seqid_2_sgd}->{$seq1}) {next;}
      my $orf1 = $names->{seqid_2_sgd}->{$seq1};

      if (!exists $mbinfo->{seq2dom}->{$seq1}) {next;}

      if (!defined $yeastgfp->{orf2loc}->{$orf1}) {next;}

      my $loc1 = $yeastgfp->{orf2loc}->{$orf1} ;
      my @loc1 = split(/\,/, $loc1) ;

      my $hloc1 ;
      foreach my $tloc (@loc1) {
         if ($tloc eq 'ambiguous') {next;}
         $hloc1->{$tloc}++ ;
      }
      my $numloc1 = keys %{$hloc1} ; if ($numloc1 == 0) {next;}
      $seq2hloc->{$seq1} = $hloc1 ;

      foreach my $j (1 .. $#{$mbinfo->{seq2dom}->{$seq1}}) {
         my $tdom = $mbinfo->{seq2dom}->{$seq1}->[$j] ;
         my $scoptype = $tdom->{scop} ;
         my @p = split(/\./, $scoptype) ;

         if ($#p < 2) {next;}
         my $sf = $p[0].'.'.$p[1].'.'.$p[2] ;
         $seq2hdom_sf->{$seq1}->{$sf} = 1 ;

         if ($#p < 3) {next;}
         $seq2hdom_fa->{$seq1}->{$scoptype} = 1 ;
      }
   }

   foreach my $seq1 (sort keys %{$seq2hloc}) {
      foreach my $seq2 (sort keys %{$seq2hloc}) {
         if ($seq1 eq $seq2) {next;}

         my $coifam = 0 ;
         foreach my $d (keys %{$seq2hdom_fa->{$seq1}}) {
            if ($coifam == 1) {last;}
            foreach my $d2 (keys %{$seq2hdom_fa->{$seq2}}) {
               if (exists $intfa12->{$d."\n".$d2}) {$coifam=1; last;} } }

         my $coisf = 0 ;
         foreach my $d (keys %{$seq2hdom_sf->{$seq1}}) {
            if ($coisf == 1) {last;}
            foreach my $d2 (keys %{$seq2hdom_sf->{$seq2}}) {
               if (exists $intsf12->{$d."\n".$d2}) {$coisf=1; last;} } }

         my $coloc = 0 ;
         foreach my $l (keys %{$seq2hloc->{$seq1}}) {
            if (exists $seq2hloc->{$seq2}->{$l}) {$coloc= 1; next;}}

         if ($coifam && $coloc) {$counts->{coifam_coloc}++;}
         if ($coisf && $coloc) {$counts->{coisf_coloc}++;}
         if ($coisf) {$counts->{coisf}++;}
         if ($coifam) {$counts->{coifam}++;}
         if ($coloc) {$counts->{coloc}++;}
      }
   }

   my $numprot = keys %{$seq2hloc} ;
   print "num proteins with a modoel and loc data: ".$numprot."\n" ;
   print "coifam: $counts->{coifam}\n" ;
   print "coisf: $counts->{coisf}\n" ;
   print "coloc: $counts->{coloc}\n\n\n" ;
   print "coifam and coloc: $counts->{coifam_coloc}\n" ;
   print "p(coloc | coifam): ".sprintf("%.3f",
      ($counts->{coifam_coloc} / $counts->{coifam}))."\n\n\n" ;
   print "coisf and coloc: $counts->{coisf_coloc}\n" ;
   print "p(coloc | coisf): ".sprintf("%.3f",
      ($counts->{coisf_coloc} / $counts->{coisf}))."\n" ;

}


sub load_seqid_2_sgd {

   my $in = shift ; my $specs ;
   if (exists $in->{specs}) { $specs = $in->{specs} ; }
   else { $specs = modtie::set_modtie_specs() ; }

   open(FILE, $specs->{seqid_2_sgd_fn}) ;
   my ($seqid_2_sgd, $sgd_2_seqid) ;
   while (my $line = <FILE>) {
      chomp $line;
      my ($seqid, $sgd) = split(/\t/, $line) ;
      $seqid_2_sgd->{$seqid} = $sgd ;
      $sgd_2_seqid->{$sgd} = $seqid ;
   }


   return {
      seqid_2_sgd => $seqid_2_sgd,
      sgd_2_seqid => $sgd_2_seqid,
   } ;

}


sub load_bind_nrints {

   my $in = shift ; my $specs ;
   if (exists $in->{specs}) { $specs = $in->{specs} ; }
   else { $specs = modtie::set_modtie_specs() ; }

   open(FILE, $specs->{bind_nrints}) ;
   my $int2info ;
   my $inthash;
   my $tax2nrints;
   while (my $line = <FILE>) {
      chomp $line;
      my ($rgid, $a_type, $a_db, $a_acc, $a_id, $a_tax,
                 $b_type, $b_db, $b_acc, $b_id, $b_tax) = split(/\t/, $line) ;
      
      $int2info->{$rgid} = {
         a_type => $a_type,
         a_db => $a_db,
         a_acc => $a_acc,
         a_id => $a_id,
         a_tax => $a_tax,
         b_type => $b_type,
         b_db => $b_db,
         b_acc => $b_acc,
         b_id => $b_id,
         b_tax => $b_tax
      } ;

      $inthash->{$a_id}->{$b_id}++ ;
      $inthash->{$b_id}->{$a_id}++ ;
      $tax2nrints->{$a_tax."\n".$b_tax}->{$rgid}++ ;
   }
   close(FILE) ;

   return {
      p2p => $inthash,
      nrints => $int2info,
      tax2rgid => $tax2nrints
   } ;

}


sub load_bind_division {

   my $in = shift ; my $specs ;
   if (exists $in->{specs}) { $specs = $in->{specs} ; }
   else { $specs = modtie::set_modtie_specs() ; }

   open(FILE, $specs->{bind_division}) ;
   my $int2info ;
   my $inthash ;
   my $rgid2div ;
   my $div2rgid ;
   while (my $line = <FILE>) {
      chomp $line;
      my ($div, $rgid) = split(/\t/, $line) ;
      $rgid2div->{$rgid}->{$div}++ ;
      $div2rgid->{$div}->{$rgid}++ ;
   }
   close(FILE) ;

   return {
      div2rgid => $div2rgid,
      rgid2div => $rgid2div,
   } ;

}


sub load_bind_complex2subunits {

   my $in = shift ; my $specs ;
   if (exists $in->{specs}) { $specs = $in->{specs} ; }
   else { $specs = modtie::set_modtie_specs() ; }

   open(FILE, $specs->{bind_complex2subunits}) ;
   my $complex2subunits;
   my $subunit2bcid;
   my $tax2complex ;
   while (my $line = <FILE>) {
      chomp $line;
      my ($bcid, $type, $db, $acc, $id, $tax, $short_label, $other_names) =
                 split(/\t/, $line) ;
      push @{$complex2subunits->{$bcid}}, {
         type => $type,
         db => $db,
         acc => $acc,
         id => $id,
         tax => $tax,
         short_label => $short_label,
         other_names => $other_names
      } ;
      $subunit2bcid->{$id}->{$bcid}++ ;
      $tax2complex->{$tax}->{$bcid}++ ;
   }
   close(FILE) ;

   return {
      complex2subunits => $complex2subunits,
      subunit2bcid => $subunit2bcid,
      tax2bcid => $tax2complex,
   } ;
}


sub readin_sgd_protein_properties {

   my $in = shift ; my $specs ;
   if (exists $in->{specs}) { $specs = $in->{specs} ; }
   else { $specs = modtie::set_modtie_specs() ; }

   open(FILE, $specs->{sgd_protein_properties_fn}) ;
   my $props;
   my $sgdid2orf ;
   my $orf2sgdid ;
   while (my $line = <FILE>) {
      chomp $line;
      my ($yorf, $sgdid, $mw, $pi, $cai, $protlen, $ntermseq, $ctermseq,
          $codbias, $ala, $arg, $asn, $asp, $cys, $gln, $glu, $gly, $his,
          $ile, $leu, $lys, $met, $phe, $pro, $ser, $thr, $trp, $tyr, $val,
          $fop, $gravy, $aromaticity, $type) = split(/\t/, $line) ;
      push @{$props->{$yorf}}, {
         yorf => $yorf,
         sgdid => $sgdid,
         mw => $mw,
         pi => $pi,
         len => $protlen,
         nterm => $ntermseq,
         cterm => $ctermseq,
         type => $type
      } ;
      $orf2sgdid->{$yorf} = $sgdid ;
      $sgdid2orf->{$sgdid} = $yorf ;
   }

   open(FILE, $specs->{sgd_features_fn}) ;
   my $orf2name; my $name2orf;
   while (my $line = <FILE>) {
      chomp $line;
      my ($sgdid, $type, $qual, $name, $genename, $alias, $parentname, $sec_sgdid, $chromosome, $start_coord, $stop_coord, $strand, $genpos, $coord_ver, $seq_ver, $descr) = split(/\t/, $line) ;
      if ($type eq 'ORF') {
#         print STDERR" $name is called $genename\n" ;
         if (defined $genename && $genename ne '') {
            $name2orf->{$genename} = $name;
            $orf2name->{$name} = $genename ; } } }
   close(FILE) ;

   return {
      properties => $props,
      orf_2_sgdid => $orf2sgdid,
      orf_2_name => $orf2name,
      name_2_orf => $name2orf,
      sgdid_2_orf => $sgdid2orf,
   } ;
}


sub commathrees {
   my $num = shift ;

   my @a = ();
   while($num =~ /\d\d\d\d/)
   {
      $num =~ s/(\d\d\d)$//;
      unshift @a,$1;
   }
   unshift @a,$num;
   return join(',',@a);
}


sub build_scop_graph {

   require LGL ;
   my $in = shift;

   print STDERR "Laying out SCOP family map: ";
   my $coords = LGL::edges2coords({edges => $in->{edges}}) ;

   my ($graph_fh, $graph_fn) = tempfile("scopnet_graph.XXXXX", SUFFIX => ".lgl.ps") ;
   LGL::coords2ps({
      fh => $graph_fh,
      coords => $coords,
      edges => $in->{edges},
      ecols => $in->{edges}
   }) ;
   close($graph_fh) ;
   print STDERR "X\n" ;

}

1 ;
