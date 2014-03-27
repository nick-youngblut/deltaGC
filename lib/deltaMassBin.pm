package deltaMassBin;

use 5.006;
use strict;
use warnings;

=head1 NAME

deltaMassBin - scripts for running deltaMassBin.pl

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

=head2 load_deltaGC_table

Loading deltaGC table produced by deltaMass.
Loading by genome [column 1]

=cut

push @EXPORT_OK, 'load_deltaGC_table';

sub load_deltaGC_table{
  use List::MoreUtils qw/minmax/;

  my ($infile) = @_;

  my $fh;
  $ARGV{'<filename>'} eq '-' ? $fh = \*STDIN :
    open $fh, $ARGV{'<filename>'} or die $!;

  my %tbl;
  my @densities;
  while(<>){
    if($.==1){ 
      chomp;
      print join("\t", $_, qw/bin_min bin_max count median median_rank/), "\n";
      next; }
    
    chomp;
    next if /^\s$/;
    my @l = split /\t/;

    $tbl{$l[0]}{$.} = [@l[1..$#l]];

    push @densities, $l[9];
  }

  return \%tbl, minmax(@densities);
}


=head2 binByDensity

making a hash of bins foreach buoyant density

=cut

push @EXPORT_OK, 'binByDensity';

sub binByDensity{
  my ($tbl_r, $binRanges_r) = @_;

  my %bins;
  foreach my $genome (keys %$tbl_r){
    foreach  my $Uid (keys %{$tbl_r->{$genome}}){
      my $density = ${$tbl_r->{$genome}{$Uid}}[8];
      foreach my $bin (@$binRanges_r){
	# intializing bin
	unless(exists $bins{$genome}{$bin->[2]}){
	  $bins{$genome}{$bin->[2]}{count} = 0;
	  $bins{$genome}{$bin->[2]}{row} = [('NA') x 9];
	}

	# adding to bin
	if( $density >= $bin->[0] && 
	    $density < $bin->[1]){
	  $bins{$genome}{$bin->[2]}{count}++;
	  $bins{$genome}{$bin->[2]}{row} = $tbl_r->{$genome}{$Uid};
	}
      }
    }
    delete $tbl_r->{$genome};
  }

  #print Dumper %bins; exit;
  return \%bins;
}


=head2 calcMedianRank

calculating median density & ranking genomes by median

=cut

push @EXPORT_OK, 'calcMedianRank';

sub calcMedianRank{
  use Statistics::Descriptive;
  my ($tbl_r) = @_;

  # calculating medians
  my %medians;
  foreach my $genome (keys %$tbl_r){
    my $stat = Statistics::Descriptive::Full->new();

    foreach my $Uid(keys %{$tbl_r->{$genome}}){
      $stat->add_data(${$tbl_r->{$genome}{$Uid}}[8]); # adding density values
    }
    $medians{$genome}{median} = $stat->median();
  }
  # ranking by median
  my $rank=0;
  foreach my $genome (sort{$medians{$a}{median}<=>$medians{$b}{median}} 
		      keys %medians){
    $rank++;
    $medians{$genome}{rank} = $rank;
  }
  
  #print Dumper %medians; exit;
  return \%medians;
}




=head1 AUTHOR

Nick Youngblut, C<< <ndy2 at cornell.edu> >>

=head1 BUGS

Please report any bugs or feature requests to C<bug-deltaMassBin at rt.cpan.org>, or through
the web interface at L<http://rt.cpan.org/NoAuth/ReportBug.html?Queue=deltaMassBin>.  I will be notified, and then you'll
automatically be notified of progress on your bug as I make changes.




=head1 SUPPORT

You can find documentation for this module with the perldoc command.

    perldoc deltaMassBin


You can also look for information at:

=over 4

=item * RT: CPAN's request tracker (report bugs here)

L<http://rt.cpan.org/NoAuth/Bugs.html?Dist=deltaMassBin>

=item * AnnoCPAN: Annotated CPAN documentation

L<http://annocpan.org/dist/deltaMassBin>

=item * CPAN Ratings

L<http://cpanratings.perl.org/d/deltaMassBin>

=item * Search CPAN

L<http://search.cpan.org/dist/deltaMassBin/>

=back


=head1 ACKNOWLEDGEMENTS


=head1 LICENSE AND COPYRIGHT

Copyright 2014 Nick Youngblut.

This program is free software; you can redistribute it and/or modify it
under the terms of either: the GNU General Public License as published
by the Free Software Foundation; or the Artistic License.

See http://dev.perl.org/licenses/ for more information.


=cut

1; # End of deltaMassBin
