package deltaMass;

use 5.006;
use strict;
use warnings;

=head1 NAME

deltaMass - The great new deltaMass!

=head1 VERSION

Version 0.01

=cut

our $VERSION = '0.01';


=head1 SYNOPSIS

Quick summary of what the module does.

Perhaps a little code snippet.

    use deltaMass;

    my $foo = deltaMass->new();
    ...

=head1 EXPORT

A list of functions that can be exported.  You can delete this section
if you don't export anything, such as for a purely object-oriented module.

=head1 SUBROUTINES/METHODS

=cut

use base 'Exporter';
our @EXPORT_OK = qw/
calc_GC
/;
use Carp qw/ confess /;
use Data::Dumper;

=head2 calc_GC


=cut

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
  my %GC = ( G => 0, C => 0);
  while( $seq =~ /G/gi ){
    $GC{G}++;
  }
  while( $seq =~ /C/gi ){
    $GC{C}++;
  }
  my $GC_sum = $GC{G} + $GC{C};

  my $GC_content = $GC_sum / $len * 100;

  return $GC_content, $GC_sum, $len;
}

=head2 load_fasta

=cut

push @EXPORT_OK, 'load_fasta';

sub load_fasta {
  my $fasta_in = shift;
  open IN, $fasta_in or die $!;
  my (%fasta, $tmpkey);
  while(<IN>){
    chomp;
    s/#.+//;
    next if  /^\s*$/;
    if(/>.+/){
      $tmpkey = $_;# changing key
      $tmpkey =~ s/>| .+//g;
      $fasta{$tmpkey} = "";
    }
    else{$fasta{$tmpkey} .= $_; }
  }
  close IN;
  #print Dumper %fasta; exit;
  return \%fasta;
}

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
