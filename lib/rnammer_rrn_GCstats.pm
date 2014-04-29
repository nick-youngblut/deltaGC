package rnammer_rrn_GCstats;

use 5.006;
use strict;
use warnings;

=head1 NAME

rnammer_rrn_GCstats - scripts for running rnammer_rrn_GCstats.pl

=head1 VERSION

Version 0.01

=cut

our $VERSION = '0.01';


=head1 SYNOPSIS

=head1 EXPORT


=head1 SUBROUTINES/METHODS

=cut

use base 'Exporter';
our @EXPORT_OK = '';
use Carp qw/ confess carp /;
use Data::Dumper;

=head2 parse_rnammer_gff

Parsing rnammer output in gff format.

genome => scaffold => [gene => start|end|strand]

=cut

push @EXPORT_OK, 'parse_rnammer_gff';

sub parse_rnammer_gff{
  my ($gff_file, $score_r) = @_;
  
  open IN, $gff_file or die $!;
  my %gff;
  my %filt_summary;
  while(<IN>){
    chomp;
    next if /^\s*$/;
    next if /^\s*#/;
    
    my @l = split /\t/;
    die "ERROR: line $. of gff file does not have 10 columns\n"
      unless scalar @l == 10;
    
    # filtering by score
    die "ERROR: do not recognize gene '$l[9]'\n"
      unless exists $score_r->{$l[9]};
    $filt_summary{total}{$l[9]}++;
    unless($l[6] >= $score_r->{$l[9]}){
      $filt_summary{filtered}{$l[9]}++;
      next;
    }

    push @{$gff{$l[0]}{$l[1]}{$l[9]}}, 
      { start => $l[4],
	end => $l[5],
	strand => $l[7] };
  }

  # status on filtering
  print STDERR "#--- filtering summary ---#\n";
  foreach my $feat (keys %{$filt_summary{total}}){
    $filt_summary{filtered}{$feat} = 0 unless
      exists $filt_summary{filtered}{$feat};

    printf STDERR "$feat: %i of %i filtered\n",
      $filt_summary{filtered}{$feat},
	$filt_summary{total}{$feat};
  }

  #print Dumper %gff; exit;
  
  return \%gff;
}

=head2 pairwise_dist_ssu_lsu

pairwise distances between between ssu & lsu features
foreach genome->scaffold

ITS regions added to gff hash 

=cut

push @EXPORT_OK, 'pairwise_dist_ssu_lsu';

sub pairwise_dist_ssu_lsu{
  my ($gff_r) = @_;

  my %pdists;
  foreach my $genome (keys %$gff_r){
    foreach my $scaf (keys %{$gff_r->{$genome}}){
      for my $i (0..$#{$gff_r->{$genome}{$scaf}{'16s_rRNA'}}){
	for my $ii (0..$#{$gff_r->{$genome}{$scaf}{'23s_rRNA'}}){
	  my $ssu_feat = $gff_r->{$genome}{$scaf}{'16s_rRNA'}->[$i];
	  my $lsu_feat = $gff_r->{$genome}{$scaf}{'23s_rRNA'}->[$ii];

	  # features must be on same strand
	  next unless $ssu_feat->{strand} eq $lsu_feat->{strand};

	  # determining which orientation is min
	  my $ssu_lsu =  $lsu_feat->{start} - $ssu_feat->{end};
	  $ssu_lsu = 1000000000 if $ssu_lsu < 0;
	  my $lsu_ssu =  $ssu_feat->{start} - $lsu_feat->{end};
	  $lsu_ssu = 1000000000 if $lsu_ssu < 0;
	  
	  if( $ssu_lsu < $lsu_ssu ){
	    $pdists{$genome}{$scaf}{$i}{$ii} = [$ssu_lsu,
						$ssu_feat->{end}, 
						$lsu_feat->{start},
						$ssu_feat->{strand}];
	  }
	  elsif( $ssu_lsu > $lsu_ssu ){
	    $pdists{$genome}{$scaf}{$i}{$ii} = [$lsu_ssu,
						$lsu_feat->{end},
						$ssu_feat->{start}, 
						$ssu_feat->{strand}];	    
	  }
	  else{ next; }
	}
      }
    }
  }
  
  #print Dumper %pdists; exit;
  return \%pdists;
}


=head2 determine_ITS

Determine the ITS region between ssu & lsu features.

ITS selected as min distance between ssu and lsu features,
with distance must be <= interval cutoff

=cut

push @EXPORT_OK, 'determine_ITS';

sub determine_ITS{
  use List::Util qw/min/;
  my ($gff_r, $pdists_r, $i_cutoff) = @_;

  foreach my $genome (keys %$pdists_r){
    foreach my $scaf (keys %{$pdists_r->{$genome}}){
      my @dists;
      foreach my $ssu_i (keys %{$pdists_r->{$genome}{$scaf}}){
	my @min_dist;
	foreach my $lsu_i (keys %{$pdists_r->{$genome}{$scaf}{$ssu_i}}){
	  if(! @min_dist){
	    @min_dist = @{$pdists_r->{$genome}{$scaf}{$ssu_i}{$lsu_i}};
	  }
	  else{
	    @min_dist = @{$pdists_r->{$genome}{$scaf}{$ssu_i}{$lsu_i}}
	      if $pdists_r->{$genome}{$scaf}{$ssu_i}{$lsu_i}->[0] < 
		$min_dist[0];
	  }
	}
	#print Dumper @min_dist; exit;

	# loading ITS into gff hash
	push @{$gff_r->{$genome}{$scaf}{ITS}},
	  { start => $min_dist[1],
	    end => $min_dist[2],
	    length => $min_dist[2] - $min_dist[1],
	    strand => $min_dist[3]
	    };
      }
    }
  }   

  #print Dumper %$gff_r; exit;
}


=head2 load_fasta

loading fasta file

returning hash: name=>seq

=cut

push @EXPORT_OK, 'load_fasta';

sub load_fasta{
  use File::Spec;
  my ($genome, $fasta_dir) = @_;

  my @parts = File::Spec->splitpath($genome);
  $parts[2] =~ s/\.[^.]+$|$//;
  
  # trying fasta extensions
  my @ext = qw/.fasta .fna .fa .fsa/;
  my @hit;
  foreach my $ext (@ext){
    my $file = "$fasta_dir$parts[2]$ext";
    push @hit, $file if -e $file;
  }
  die "ERROR: cannot find fasta file for genome '$genome'\n"
    unless @hit;
  die "ERROR: multiple fasta files hit genome '$genome'\n"
    if scalar @hit > 1;
  
  open IN, $hit[0] or die $!;
  
  my (%fasta, $tmpkey);
  while(<IN>){
    chomp;
    s/#.+//;
    next if  /^\s*$/;
    if(/>.+/){
      ($tmpkey = $_) =~ s/>| .+//g;  # changing key
      $fasta{$tmpkey} = "";
      }
    else{$fasta{$tmpkey} .= $_; }
    }
  close IN;
  
  #print Dumper %fasta; exit;
  return \%fasta;

}


=head2 substr_frags

Parsing out sequence fragments of interest from the genome fasta

Return: %frags: {gene} => [fragment]

=cut

push @EXPORT_OK, 'substr_frags';

sub substr_frags{
  my ($gff_r, $genome, $fasta_r) = @_;

  my %frags;
  foreach my $scaf ( keys %{$gff_r->{$genome}} ){
    die "ERROR: cannot find scaffold '$scaf' in fasta for genome '$genome'\n"
      unless exists $fasta_r->{$scaf};

    foreach my $gene (keys %{$gff_r->{$genome}{$scaf}}){
      foreach my $feat_r ( @{$gff_r->{$genome}{$scaf}{$gene}} ){
	push @{$frags{$gene}}, substr($fasta_r->{$scaf}, $feat_r->{start}-1,
				      $feat_r->{end} - $feat_r->{start} + 1);
      }
    }
  }

  #print Dumper %frags; exit;
  return \%frags;
}





=head1 AUTHOR

Nick Youngblut, C<< <ndy2 at cornell.edu> >>

=head1 BUGS

Please report any bugs or feature requests to C<bug-rnammer_rrn_GCstats at rt.cpan.org>, or through
the web interface at L<http://rt.cpan.org/NoAuth/ReportBug.html?Queue=rnammer_rrn_GCstats>.  I will be notified, and then you'll
automatically be notified of progress on your bug as I make changes.


=head1 SUPPORT

You can find documentation for this module with the perldoc command.

    perldoc rnammer_rrn_GCstats


You can also look for information at:

=over 4

=item * RT: CPAN's request tracker (report bugs here)

L<http://rt.cpan.org/NoAuth/Bugs.html?Dist=rnammer_rrn_GCstats>

=item * AnnoCPAN: Annotated CPAN documentation

L<http://annocpan.org/dist/rnammer_rrn_GCstats>

=item * CPAN Ratings

L<http://cpanratings.perl.org/d/rnammer_rrn_GCstats>

=item * Search CPAN

L<http://search.cpan.org/dist/rnammer_rrn_GCstats/>

=back


=head1 ACKNOWLEDGEMENTS


=head1 LICENSE AND COPYRIGHT

Copyright 2014 Nick Youngblut.

This program is free software; you can redistribute it and/or modify it
under the terms of either: the GNU General Public License as published
by the Free Software Foundation; or the Artistic License.

See http://dev.perl.org/licenses/ for more information.


=cut

1; # End of rnammer_rrn_GCstats
