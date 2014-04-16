package rnammer_tandem_rrn;

use 5.006;
use strict;
use warnings;

=head1 NAME

rnammer_tandem_rrn - scripts for running rnammer_tandem_rrn.pl

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
use Set::IntervalTree;
use List::MoreUtils qw/any none/;

=head2 load_gff_file_list

Converting file list to array of files.
Checking that files exist.

=cut

push @EXPORT_OK, 'load_gff_file_list';

sub load_gff_file_list{
  my $file_in = shift;
  
  my $fh;
  $file_in eq '-' ? $fh = \*STDIN :
    open $fh, $file_in or die $!;

  chomp(my @file_list = <$fh>);
  map{ die "ERROR: cannot find '$_'\n" unless -e $_ } @file_list;

  return \@file_list;
}

=head gff2itree

Parsing each gff and loading features into interval tree

=cut

push @EXPORT_OK, 'gff2itree';

sub gff2itree{
  my ($gff_list_r) = @_;

  my %itrees;
  my %ssu_loc;

  foreach my $file (@$gff_list_r){
    open IN, $file or die $!;
    while(<IN>){
      chomp;
      next if /^#/;
      next if /^\s*$/;

      my @l = split /\t/;
      die "ERROR: line $. of file '$file' does not have 9 columns!\n"
	unless @l == 9;
      
      # loading start-end by strand
      my ($start, $end, $strand) = @l[(3,4,6)];
      ($start, $end) = ($end, $start) if $start > $end; # start <= end
      my ($seqname, $feature) = @l[(0,8)];

      $itrees{$file}{$seqname}{$strand} = Set::IntervalTree->new
	unless exists $itrees{$file}{$seqname}{$strand};      
      $itrees{$file}{$seqname}{$strand}->insert($feature, $start, $end);      

      # keeping track of 16S rRNA start-end
      push @{$ssu_loc{$file}}, {
				seqname => $seqname,
				start => $start,
				end => $end,
				strand => $strand
			       } if $feature eq '16s_rRNA';
      
    }
    close IN;
  }

  #print Dumper %ssu_loc; exit;
  return \%itrees, \%ssu_loc;
}

=head2 check_tandem

Checking to see if 16S rRNA gene has 2 adjacent 23S rRNA genes.
Adjacency determined by '$dist'

=cut

push @EXPORT_OK, 'check_tandem';

sub check_tandem{
  my ($itrees_r, $ssu_loc_r, $dist) = @_;
  
  my %summary;

  foreach my $file (keys $ssu_loc_r){
    
    # sanity check 
    die "ERROR: cannot find '$file' in itrees\n"
      unless exists $itrees_r->{$file};
    # loading total copy number
    $summary{$file}{copy} = scalar @{$ssu_loc_r->{$file}};
    # counting tandems
    foreach my $r ( @{$ssu_loc_r->{$file}} ){
      
      die "ERROR: cannot find scaffold->strand for '$file'\n"
	unless exists $itrees_r->{$file}{ $r->{seqname} }{ $r->{strand} };

      # checking upstream
      my $ret_up = $itrees_r->{$file}{ $r->{seqname} }{ $r->{strand} }->fetch( $r->{start} - $dist, $r->{start} );

      # checking downstream
      my $ret_down = $itrees_r->{$file}{ $r->{seqname} }{ $r->{strand} }->fetch( $r->{end}, $r->{end} + $dist );

      # making sure genes don't overlap 
      my $ret_over = $itrees_r->{$file}{ $r->{seqname} }{ $r->{strand} }->fetch( $r->{start}, $r->{end}); 
     
      # adding to count if both downstream & upstream 23S rRNA gene
      my $q23 = '23s_rRNA';
      $summary{$file}{'tandem'}++ if any{ $_ eq $q23 } @$ret_up 
	and any{ $_ eq $q23} @$ret_down 
	  and none{ $_ eq $q23 } @$ret_over;
    }								  
  }
  
 # print Dumper %summary; exit;
  return \%summary;
}



=head1 AUTHOR

Nick Youngblut, C<< <ndy2 at cornell.edu> >>

=head1 BUGS

Please report any bugs or feature requests to C<bug-rnammer_tandem_rrn at rt.cpan.org>, or through
the web interface at L<http://rt.cpan.org/NoAuth/ReportBug.html?Queue=rnammer_tandem_rrn>.  I will be notified, and then you'll
automatically be notified of progress on your bug as I make changes.


=head1 SUPPORT

You can find documentation for this module with the perldoc command.

    perldoc rnammer_tandem_rrn


You can also look for information at:

=over 4

=item * RT: CPAN's request tracker (report bugs here)

L<http://rt.cpan.org/NoAuth/Bugs.html?Dist=rnammer_tandem_rrn>

=item * AnnoCPAN: Annotated CPAN documentation

L<http://annocpan.org/dist/rnammer_tandem_rrn>

=item * CPAN Ratings

L<http://cpanratings.perl.org/d/rnammer_tandem_rrn>

=item * Search CPAN

L<http://search.cpan.org/dist/rnammer_tandem_rrn/>

=back


=head1 ACKNOWLEDGEMENTS


=head1 LICENSE AND COPYRIGHT

Copyright 2014 Nick Youngblut.

This program is free software; you can redistribute it and/or modify it
under the terms of either: the GNU General Public License as published
by the Free Software Foundation; or the Artistic License.

See http://dev.perl.org/licenses/ for more information.


=cut

1; # End of rnammer_tandem_rrn
