=head1 NAME

SGE.pm - MODTIE routines to interact with an SGE cluster

=head1 DESCRIPTION

The SGE.pm perl library contains subroutines to submit jobs and retrieve
output from an SGE cluster

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


package modtie::SGE ;
use strict;
use Exporter ;
our @ISA = qw/Exporter/ ;
our @EXPORT_OK = qw// ;

use modtie ;
use Cwd ;


sub _clust_split_ins {

   use File::Basename qw/basename/ ;

   my $in = shift ;

   if (!exists $in->{fn}) {
      die "_clust_split_ins: input file not specified\n" ; }

   if (!exists $in->{dir}) {
      $in->{dir} = "./" ; }

   my $inbase = basename($in->{fn}) ;

   my $header_lines = '';
   my $num_lines = 0 ;
   open(INF, $in->{fn}) ;
   while (my $line = <INF>) {
      if ($line =~ /^#SET /) {
         $header_lines .= $line ;
      } else {
         $num_lines++ ;
      }
   }
   close(INF);

   my $splitlines = POSIX::ceil($num_lines / $in->{numjobs});

   if (exists $in->{minlines} && ($splitlines < $in->{minlines})) {
      $splitlines = $in->{minlines} ;
   }

   my @splitfiles ;

   $num_lines = 0 ;

   my $tasklist = '';
   my $cur_splitnum = 1 ;

   my $cur_fn = "split.$inbase.$cur_splitnum" ;
   push @splitfiles, $cur_fn ;
   $tasklist .= "'$cur_fn' " ;

   open(INF, $in->{fn}) ;
   open(OUTF, ">$in->{dir}/$cur_fn") ;
   print OUTF $header_lines ;
   while (my $line = <INF>) {
      if ($line =~ /^#SET/) {next;}
      if ($num_lines > $splitlines) {
         close(OUTF) ;
         $cur_splitnum++ ;
         $cur_fn = "split.$inbase.$cur_splitnum" ;
         push @splitfiles, $cur_fn ;
         $tasklist .= "'$cur_fn' " ;
         open(OUTF, ">$in->{dir}/$cur_fn") ;
         print OUTF $header_lines ;
         $num_lines = 0 ;
      }
      print OUTF $line ;
      $num_lines++ ;
   }
   close(OUTF) ;
   close(INF);
   $tasklist =~ s/ $// ;


   return {
      numjobs => $cur_splitnum,
      splitfiles => \@splitfiles,
      tasklist => $tasklist
   } ;
}


sub _clust_merge_outs {

   my $in = shift ;
   open(OUTFN, ">$in->{out_fn}") ;

   if (exists $in->{headers}) {
      print OUTFN '#'.join("\t", @{$in->{headers}})."\n"; }

   open(ERRFN, ">$in->{err_fn}") ;
   foreach my $j ( 1 .. $in->{numjobs}) {
      if (-s "$in->{outdir}/$in->{script_fn}.o$in->{job_id}.$j") {
         open(SGEOUT, "$in->{outdir}/$in->{script_fn}.o$in->{job_id}.$j") ;
         while (my $line = <SGEOUT>) {
            if ($line !~ /^Warning: no access to tty/ &&
                $line !~ /^Thus no job control/ &&
                $line !~ /^#sgejob run/) {
               print OUTFN $line ;
            }
         }
         close(SGEOUT) ;
      }

      if (-s "$in->{outdir}/$in->{script_fn}.e$in->{job_id}.$j") {
         open(SGEERR, "$in->{outdir}/$in->{script_fn}.e$in->{job_id}.$j") ;
         while (my $line = <SGEERR>) {
            print ERRFN $line ; }
         close(SGEERR) ;
      }
   }
   close(OUTFN) ;
   close(ERRFN) ;

   return ;
}


sub _clust_qsub {
   my $in = shift;

   my $qsub_command = "qsub ".$in->{sgescript_fn} ;
   if (exists $in->{host} && defined $in->{host} && $in->{host} ne '') {
      my $cwd = getcwd() ;
      $qsub_command = "ssh ".$in->{host}." \"cd $cwd;".$qsub_command."\""; }

#   my $status = `ssh login-eddy \"cd $cwd; qsub $in->{sgescript_fn}\"` ;
   my $status = `$qsub_command` ;
   chomp $status ;
   if ($status !~ /has been submitted/) {
      die "failed to qsub $in->{sgescript_fn}";}
#   my ($job_id) = ($status =~ /Your job ([0-9]+)\./) ;
# sometimes its Your job 12345 ("job.SGE.sh") has been submitted.
#others it is Your job 30614.1-31:1 ("job.SGE.sh") has been submitted.

   my ($job_id) = ($status =~ /Your job.* ([0-9]+)/) ;
   print STDERR "submitted $in->{sgescript_fn} job $job_id\n" ;

   return $job_id ;
}


sub _clust_qstat {

   my $in = shift;
   my $qstat_command = "qstat -j $in->{job_id} > /tmp/qstat.$$.out 2>&1" ;
   if (exists $in->{host} && defined $in->{host} && $in->{host} ne '') {
      $qstat_command = "ssh ".$in->{host}." ".$qstat_command ; }

#   system("ssh login-eddy qstat -j $in->{job_id} > /tmp/qstat.$$.out 2>&1") ;
   system($qstat_command) ;
   my $status = `cat /tmp/qstat.$$.out` ;
   system("rm /tmp/qstat.$$.out") ;
   if ($status =~ /Following jobs do not exist/) {
      return 1;
   } else {
      return 0;
   }

}




1 ;
