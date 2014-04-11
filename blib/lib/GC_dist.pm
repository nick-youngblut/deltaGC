package GC_dist;

use 5.006;
use strict;
use warnings;

=head1 NAME

GC_dist - scripts for running GC_dist.pl

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
use Text::ParseWords;
use Bio::SeqIO;
use Math::Random qw/
random_uniform_integer
random_uniform
random_normal
random_exponential/;


=head2 parse_desc

parsing sequence description of reads from grinder

=head3 input

* description <string>
* amplicon? [false] <bool>

=head3 output

* @(reference, start, end, strand)

=cut

push @EXPORT_OK, 'parse_desc';

sub parse_desc{
  my ($desc, $amp_b) = @_;

  $desc =~ s/^>//;
  my %vals = quotewords("=| ", 0, $desc);

  # checking read descriptions
  my @chk = ('reference');
  if($amp_b){ push @chk, 'amplicon'; }
  else{ push @chk, 'position'; }
  map{ die "ERROR: cannot find '$_' in read!\n"
         unless exists $vals{$_} } @chk;


  # strand
  my $strand;
  if($vals{$chk[1]} =~ /complement/){
    $strand = -1;
  }
  else{
    $strand = 1;
  }

  # postion
  $vals{$chk[1]} =~ /(\d+)\.\.(\d+)/;

  return [$vals{reference}, $1, $2, $strand];
}


=head2 correct_fasta

Making sure sequence lines for each entry are the same length.

=cut

push @EXPORT_OK, 'correct_fasta';

sub correct_fasta{
  my $genome_file = shift;

  (my $genome_out = $genome_file) =~ s/\.[^.]+$/_fix.fasta/;

  my $seqin = Bio::SeqIO->new(-file => $genome_file, -format => 'fasta');
  my $seqout = Bio::SeqIO->new(-file => ">$genome_out", -format => 'fasta');

  while(my $seq = $seqin->next_seq){
    $seqout->write_seq($seq);
  }

  print STDERR "Corrected genome file written: '$genome_out'\n";
  return $genome_out;
}

=head2 calc_GC

=cut

push @EXPORT_OK, 'calc_GC';

sub calc_GC {
# calculting GC of a sequence
  my ($seq, $gap_frac) = @_;

  confess "ERROR: no sequences provided for calculating GC!\n"
    unless defined $seq;

  # length
  my $raw_len = length $seq;

  # removing gaps
  $seq =~ s/-+//g;
  my $len = length $seq;
  return ('NA',$len) if 1 - $len / $raw_len > $gap_frac; # if too many gaps in DNA

  # upcase
  $seq =~ tr/a-z/A-Z/;

  # scoring table (included ambiguous nucleotides)
  my %score = (G => 1,
               C => 1,
               R => 0.5,
               Y => 0.5,
               S => 1,
               K => 0.5,
               M => 0.5,
               B => 0.66,
               D => 0.33,
               H => 0.33,
               V => 0.66,
               N => 0.5
               );

  # GC
  my $q = join("", "[", keys %score, "]");
  $q = qr/$q/;
  my $GC_sum = 0;
  $GC_sum += $score{$1} while $seq =~ /($q)/g;
  my $GC_content = $GC_sum / $len * 100;

  return $GC_content;
}

=head2 calc_frag_GC_window

calculating GC over a window for the read

=cut 

push @EXPORT_OK, 'calc_frag_GC_window';

sub calc_frag_GC_window{
  my ($genome, $reads_r, $genome_seq, 
      $frag_size, $window, $jump,
      $gap_frac) = @_;

  # getting genome seq
  my $genome_len = length $genome_seq;

  my %gc;
  foreach my $read (keys %$reads_r) {
    my ($read_start, $read_end, $read_strand) = @{$reads_r->{$read}};

    #finding center of read
    my $center = $read_start + (abs($read_end - $read_start) / 2);
    my $frag_start = int($center - $frag_size * 0.5);

    # assuming linear chromosome
    $frag_start = 1 if $frag_start < 1;
    my $frag_end = $frag_start + $frag_size - 1;
    $frag_end = $genome_len if $frag_end > $genome_len;

    # getting seq
    my $frag_seq = substr($genome_seq, $frag_start - 1, abs($frag_end - $frag_start) + 1);

    # reversing frag if - strand
    $frag_seq = reverse $frag_seq if $read_strand == -1;
    
    # sliding window GC analysis 
    my $frag_len = length $frag_seq;
    for (my $i=0; $i<=($frag_len - 1); $i+=$jump) {	
      my $mid = int($i + $window/2);
      # GC
      my $frag_sub_seq = substr($frag_seq, $i, $window);
      $gc{$genome}{$read}{$mid}{GC} = length($frag_sub_seq) >= $window ?
	calc_GC( $frag_sub_seq, $gap_frac ) : next;   # GC_content may be 'NA'
      
      # position & buoyant_density
      $gc{$genome}{$read}{$mid}{pos} = $i + $frag_start + $mid;
      $gc{$genome}{$read}{$mid}{density} = $gc{$genome}{$read}{$mid}{GC} =~ /^[\d.]+$/ ?
	($gc{$genome}{$read}{$mid}{GC} * 0.098 / 100) + 1.660 : 'NA';
    }
  }
 # print Dumper %gc; exit;
  return \%gc;
}


=head2 write_output

Writing output table for hash produced by calc_frag_GC_window

Input: $gc hash; genome_description

=cut

push @EXPORT_OK, 'write_output';

sub write_output{
  my ($gc_r, $desc) = @_;

  foreach my $genome (keys %$gc_r){
    # writing out windows for each read
    my $read_num = 0;
    foreach my $read (keys %{$gc_r->{$genome}}){
      $read_num++;
      foreach my $start (sort {$a<=>$b} keys %{$gc_r->{$genome}{$read}}){
	print join("\t",$desc, $genome, $read, $read_num, $start,
		   $gc_r->{$genome}{$read}{$start}{pos},
		   $gc_r->{$genome}{$read}{$start}{GC},
		   $gc_r->{$genome}{$read}{$start}{density}), "\n";
      }
    }
  }
}


=head1 AUTHOR

Nick Youngblut, C<< <ndy2 at cornell.edu> >>

=head1 BUGS

Please report any bugs or feature requests to C<bug-GC_dist at rt.cpan.org>, or through
the web interface at L<http://rt.cpan.org/NoAuth/ReportBug.html?Queue=GC_dist>.  I will be notified, and then you'll
automatically be notified of progress on your bug as I make changes.




=head1 SUPPORT

You can find documentation for this module with the perldoc command.

    perldoc GC_dist


You can also look for information at:

=over 4

=item * RT: CPAN's request tracker (report bugs here)

L<http://rt.cpan.org/NoAuth/Bugs.html?Dist=GC_dist>

=item * AnnoCPAN: Annotated CPAN documentation

L<http://annocpan.org/dist/GC_dist>

=item * CPAN Ratings

L<http://cpanratings.perl.org/d/GC_dist>

=item * Search CPAN

L<http://search.cpan.org/dist/GC_dist/>

=back


=head1 ACKNOWLEDGEMENTS


=head1 LICENSE AND COPYRIGHT

Copyright 2014 Nick Youngblut.

This program is free software; you can redistribute it and/or modify it
under the terms of either: the GNU General Public License as published
by the Free Software Foundation; or the Artistic License.

See http://dev.perl.org/licenses/ for more information.


=cut

1; # End of GC_dist
