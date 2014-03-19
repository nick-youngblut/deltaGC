package readFilter;

use 5.006;
use strict;
use warnings;

=head1 NAME

readFilter - scripts for running readFilter.pl

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
use Set::IntervalTree;
use List::MoreUtils qw/ each_arrayref /;

=head2 load_read_itrees

Loading read positions into itrees for detecting overlap.

Input: all read ids = @ids; all read positions @@[scaffold, start, end]

Output: hash of interval trees

=cut

push @EXPORT_OK, 'load_read_itrees';

sub load_read_itrees{
  my ($ids_r, $pos_r) = @_;
  print STDERR "Loading interval trees...\n";

  my %itrees;
  my $ea = each_arrayref($ids_r, $pos_r);
  while(my ($x, $y) = $ea->()){
    $itrees{$$y[0]} = Set::IntervalTree->new()
      unless exists $itrees{$$y[0]};

    $itrees{$$y[0]}->insert($x, $$y[1], $$y[2]);
  }
 
  print STDERR " Number of interval trees created: ",
    scalar keys %itrees, "\n";
  return \%itrees;
}

=head2 screen_reads

Screening reads based on potential overlap

=cut

push @EXPORT_OK, 'screen_reads';

sub screen_reads{
  my ($ids_r, $pos_r, $itrees_r, $window) = @_;
  print STDERR "Screening reads for overlap +/- window...\n";

  my %keep;
  my %rm;
  my $ea = each_arrayref($ids_r, $pos_r);
  while( my($x,$y) = $ea->()){
    die "ERROR: genome $$y[0] not found in iterval tree!\n"
      unless exists $itrees_r->{$$y[0]};

    next if exists $rm{$x}; # on remove list

    my $res = $itrees_r->{$$y[0]}->fetch($$y[1] - $window, 
				       $$y[2] + $window);

    if(scalar @$res > 1){ # overlapping (or nearly) reads
      for my $i ( 1..$#$res){ 
	$rm{$$res[$i]} = 1;
      }
      $keep{$$res[0]}= 1;
    }
    else{
      $keep{$$res[0]} = 1;
    }
  }
  
  # status
  print STDERR "Number of reads pre-filtering: ", scalar @$ids_r, "\n";
  print STDERR "Number remaining post-filtering: ", 
    scalar keys %keep, "\n";

  #print Dumper %keep; exit;
  return \%keep;
}



=head1 AUTHOR

Nick Youngblut, C<< <ndy2 at cornell.edu> >>

=head1 BUGS

Please report any bugs or feature requests to C<bug-readFilter at rt.cpan.org>, or through
the web interface at L<http://rt.cpan.org/NoAuth/ReportBug.html?Queue=readFilter>.  I will be notified, and then you'll
automatically be notified of progress on your bug as I make changes.




=head1 SUPPORT

You can find documentation for this module with the perldoc command.

    perldoc readFilter


You can also look for information at:

=over 4

=item * RT: CPAN's request tracker (report bugs here)

L<http://rt.cpan.org/NoAuth/Bugs.html?Dist=readFilter>

=item * AnnoCPAN: Annotated CPAN documentation

L<http://annocpan.org/dist/readFilter>

=item * CPAN Ratings

L<http://cpanratings.perl.org/d/readFilter>

=item * Search CPAN

L<http://search.cpan.org/dist/readFilter/>

=back


=head1 ACKNOWLEDGEMENTS


=head1 LICENSE AND COPYRIGHT

Copyright 2014 Nick Youngblut.

This program is free software; you can redistribute it and/or modify it
under the terms of either: the GNU General Public License as published
by the Free Software Foundation; or the Artistic License.

See http://dev.perl.org/licenses/ for more information.


=cut

1; # End of readFilter
