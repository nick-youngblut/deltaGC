package deltaMass;

use 5.006;
use strict;
use warnings;

=head1 NAME

deltaMass - scripts for running deltaMass.pl

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
use Carp qw/ confess /;
use Data::Dumper;
use Text::ParseWords;
use Math::Random qw/
random_uniform_integer
random_uniform
random_normal
random_exponential/;

=head2 calc_GC

=cut

push @EXPORT_OK, 'calc_GC';

sub calc_GC {
# calculting GC of a sequence
  my ($seq) = @_;

  confess "ERROR: no sequences provided for calculating GC!\n" 
    unless defined $seq;

  # removing gaps #
  $seq =~ s/[-.]//g;

  # length
  my $len = length $seq;

  # GC
  my %GC = ( G => 0, C => 0 );
  while( $seq =~ /G/gi ){
    $GC{G}++;
  }
  while( $seq =~ /C/gi ){
    $GC{C}++;
  }
  my $GC_sum = $GC{G} + $GC{C};
  my $GC_content = $GC_sum / $len * 100;

  return $GC_content, $len;
}

=head2 load_genome_fasta

=cut

push @EXPORT_OK, 'load_genome_fasta';

sub load_genome_fasta {
  my $fasta_file = shift;
  open IN, $fasta_file or die $!;

  # loading fasta 
  my (%fasta, $tmpkey);
  while(<IN>){
    chomp;
    s/#.+//;
    next if  /^\s*$/;
    if(/>.+/){
      $tmpkey = $_;# changing key
      $tmpkey =~ s/>| .+//g; # nothing after space kept
      $fasta{$tmpkey}{seq} = "";
    }
    else{$fasta{$tmpkey}{seq} .= $_; }
  }
  close IN;

  # getting genome lengths
  map{ $fasta{$_}{len} = length $fasta{$_}{seq} } keys %fasta;

  #print Dumper %fasta; exit;
  return \%fasta;
}

=head2 load_reads

Loading reads produced by grinder.
Fasta format, but headers need to be parsed

=cut

push @EXPORT_OK, 'load_reads';

sub load_reads{
  my $reads_file = shift;
  
  open IN, $reads_file or die $!;

  my (%fasta, $ref, $Uid);
  while(<IN>){
    chomp;

    if(eof(IN)){
      $fasta{$ref}{$Uid}{seq} .= $_;
      ($fasta{$ref}{$Uid}{amp_GC}, $fasta{$ref}{$Uid}{amp_len}) 
	= calc_GC($fasta{$ref}{$Uid}{seq});
      delete $fasta{$ref}{$Uid}{seq} or die $!;
    }
    elsif(/^>/){
      # calculating GC on last sequence
      if( defined $ref ){
	($fasta{$ref}{$Uid}{amp_GC}, $fasta{$ref}{$Uid}{amp_len}) 
	  = calc_GC($fasta{$ref}{$Uid}{seq});
	delete $fasta{$ref}{$Uid}{seq} or die $!;
      }

      # loading new sequence
      my @l = split />| /, $_, 3;
      my %vals = quotewords("=| ", 0,  $l[2]);
      
      $Uid = $l[1];
      $ref = $vals{reference};
      
      ## amplicon start-end
      my ($amp_start, $amp_end, $amp_strand);
      if( $vals{amplicon} =~ /complement\((.+)\)/){ 
	($amp_end, $amp_start) = split /\.\./, $1;
	$amp_strand = -1;
      }
      else{
	($amp_start, $amp_end) = split /\.\./, $vals{amplicon};
	$amp_strand = 1;
      }
      
      $fasta{$ref}{$Uid}{seq} = "";
      $fasta{$ref}{$Uid}{amp_start} = $amp_start;
      $fasta{$ref}{$Uid}{amp_end} = $amp_end;
      $fasta{$ref}{$Uid}{amp_strand} = $amp_strand;
      $fasta{$ref}{$Uid}{desc} = $vals{description};
    }
    else{ # just part of the sequence
      next if /^\s*$/;
      $fasta{$ref}{$Uid}{seq} .= $_;
    }
  }

  close IN or die $!;

  #print Dumper %fasta; exit;
  return \%fasta;
}

=head2 get_frag_GC

creating fragments for each read & calculating GC on each

=cut

push @EXPORT_OK, 'get_frag_GC';

sub get_frag_GC{
# foreach read:
## determine random fragment size
## pull fragment from genome
### start-end around amplicon determine from a uniform distribution
## calculate fragment GC
## record GC & length
  my ($genome_seqs_r, $reads_r, $size_dist, $frag_min, $frag_max) = @_;

  map{ confess "ERROR: argument missing: $!\n" unless defined $_ } 
    ($size_dist, $frag_min, $frag_max);

  # options that should be user-provided!!!
  my $mean = 4100;
  my $stdev = 100;
  my $primer_buffer = 70;

  # interating through fragments
  foreach my $genome (keys %$reads_r){
    die "ERROR: cannot find '$genome' in genomes fasta\n"
      unless exists $genome_seqs_r->{$genome};

    foreach my $Uid (keys %{$reads_r->{$genome}}){
      my $amp_start = $reads_r->{$genome}{$Uid}{amp_start};
      my $amp_end = $reads_r->{$genome}{$Uid}{amp_end};
      my $amp_len = $reads_r->{$genome}{$Uid}{amp_len};

      # fragment size
      my $frag_size;
      if($size_dist eq 'uniform'){
       $frag_size = random_uniform_integer(1, $frag_min, $frag_max);
      }
      elsif($size_dist eq 'normal'){
	while(1){
	  $frag_size = int random_normal(1, $mean, $stdev);
	  last if $frag_size >= $frag_min && $frag_size <= $frag_max;
	}
      }
      elsif($size_dist eq 'exponential'){
	while(1){
	  $frag_size = int random_exponential(1, $mean);
	  last if $frag_size >= $frag_min && $frag_size <= $frag_max;
	}
      }
      else{ confess "ERROR: do not recognize size distribution\n"; }
      
      #print Dumper $frag_size; exit;
      
      # determine fragment start-end based on amplicon start-end
      ## amplicon center postion
      my $amp_center = int($amp_start + abs($amp_len * 0.5));
      ## fragment_start = amp_center - (frag_size * x); x = random draw from unifrom distribution 0:1
      my $x = random_uniform();
      ### start_floor = (primer_buffer + 0.5*amp_len) / frag_size
      my $floor  = ($primer_buffer + 0.5 * $amp_len) / $frag_size; 
      $x = $floor if $x < $floor;
      my $ceiling = 1 - $floor;
      $x = $ceiling if $x > $ceiling;
      my $frag_start = int( $amp_center - ($frag_size * $x) );

      ### sanity check 
      confess "ERROR: frag_start is too far from amplicon\n"
	if $frag_start - ($amp_end + $primer_buffer) > $frag_size;

      # getting fragment
      ## wrapping around genome sequence string if neede
      my $frag_seq;
      if($frag_start + $frag_size + 1 > $genome_seqs_r->{$genome}{len}){ # wrap
	my $max_len = $genome_seqs_r->{$genome}{len} - $frag_start + 1;
	$frag_seq = substr($genome_seqs_r->{$genome}{seq}, $frag_start, $max_len);
	$frag_seq .= substr($genome_seqs_r->{$genome}{seq}, 0, $frag_size - $max_len);
      }
      else{
	$frag_seq = substr($genome_seqs_r->{$genome}{seq}, $frag_start, $frag_size );
      }

      # getting fragment GC
      ($reads_r->{$genome}{$Uid}{frag_GC},$reads_r->{$genome}{$Uid}{frag_len}) =
	calc_GC( $frag_seq );
    }
  }
  
  #print Dumper $reads_r; exit;
}

=head2 write_summary

Summary table of GC 

=cut


=head1 AUTHOR

Nick Youngblut, C<< <ndy2 at cornell.edu> >>

=head1 BUGS

Please report any bugs or feature requests to C<bug-deltaMass at rt.cpan.org>, or through
the web interface at L<http://rt.cpan.org/NoAuth/ReportBug.html?Queue=deltaMass>.  I will be notified, and then you'll
automatically be notified of progress on your bug as I make changes.




=head1 SUPPORT

You can find documentation for this module with the perldoc command.

    perldoc deltaMass


You can also look for information at:

=over 4

=item * RT: CPAN's request tracker (report bugs here)

L<http://rt.cpan.org/NoAuth/Bugs.html?Dist=deltaMass>

=item * AnnoCPAN: Annotated CPAN documentation

L<http://annocpan.org/dist/deltaMass>

=item * CPAN Ratings

L<http://cpanratings.perl.org/d/deltaMass>

=item * Search CPAN

L<http://search.cpan.org/dist/deltaMass/>

=back


=head1 ACKNOWLEDGEMENTS


=head1 LICENSE AND COPYRIGHT

Copyright 2014 Nick Youngblut.

This program is free software; you can redistribute it and/or modify it
under the terms of either: the GNU General Public License as published
by the Free Software Foundation; or the Artistic License.

See http://dev.perl.org/licenses/ for more information.


=cut

1; # End of deltaMass
