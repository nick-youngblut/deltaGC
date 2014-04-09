#!/usr/bin/env perl
use strict;
use warnings;

=pod

=head1 NAME

template.pl -- basic description here

=head1 VERSION

This is version 0.0.1

=head1 USAGE

    template.pl [options]

=head1 REQUIRED ARGUMENTS

=over

=item <filename>

tab-delimited file produced by deltaGC.pl. ('-' if from STDIN)

=back

=head1 OPTIONS

=over

=item -b[inwidth] <binwidth>

fragment_buoyant_density binwidth.

Default: binwidth.default

=for Euclid:
binwidth.type: num > 0
binwidth.default: 0.001

=item -range <min>-<max>

Minimum & maximum fragment_buoyant_density

Default: min.default-max.default

=for Euclid:
min.type: num >= 0
min.default: 1.660
max.type: num >= 0
max.default: 1.770

=item --debug [<log_level>]

Set the log level. Default is log_level.default but if you provide --debug,
then it is log_level.opt_default.

=for Euclid:
    log_level.type:        int
    log_level.default:     0
    log_level.opt_default: 1

=item --version

=item --usage

=item --help

=item --man

Print the usual program information

=back

=head1 DESCRIPTION

Bin data produced by deltaGC.pl in terms of
read(amplicon)_buoyant_density & fragment_buoyant_density. 
This makes it easier to produce a heatmap
of buoyant densities with ggplot.

Also, the median fragment_buoyant_density
is calculated for each genome and each genome
is ranked by median fragment_buoyant_density.

=head1 EXAMPLES

=head1 AUTHOR

Nick Youngblut (ndy2@cornell.edu)

=head1 BUGS

There are undoubtedly serious bugs lurking somewhere in this code.
Bug reports and other feedback are most welcome.

=head1 COPYRIGHT

Copyright 2010, 2011
This software is licensed under the terms of the GPLv3

=cut

#--- modules ---#
use Data::Dumper;
use Getopt::Euclid;
use deltaGCBin qw/
load_deltaGC_table
binByDensity
calcMedianRank/;

#--- I/O error ---#
die "ERROR: max in range is < min!\n"
  unless $ARGV{-range}{min} <= $ARGV{-range}{max};

#--- MAIN ---#
my ($tbl_r, $min, $max) = load_deltaGC_table($ARGV{'<filename>'});

# calculating median density by genome & ranking
my $stats_r = calcMedianRank($tbl_r);

# binwidth=0.001
my $binRanges_r = makeBinRanges($ARGV{-binwidth}, 
				$ARGV{-range}{min}, 
				$ARGV{-range}{max});
my $bins_r = binByDensity($tbl_r, $binRanges_r);

# writing output
foreach my $genome (keys %$bins_r){
  # sanity check 
  die "ERROR: cannot find median & rank for genome: $genome\n"
    unless exists $stats_r->{$genome};

  # output by bin
  foreach my $bin (keys %{$bins_r->{$genome}}){
    print join("\t", 
	       $genome, # genome_ID
	       @{$bins_r->{$genome}{$bin}{row}}, # row values (random draw from bin)
	       split (/-/, $bin), # bin range
	       $bins_r->{$genome}{$bin}{amp_count}, # count in bin for amplicon
	       $bins_r->{$genome}{$bin}{frag_count}, # count in bin for fragment
	       $stats_r->{$genome}{median},     # median 
	       $stats_r->{$genome}{rank}        # rank by median
	       ), "\n";
  }
}


#--- Subroutines ---#
sub makeBinRanges{
  my ($binwidth, $min, $max) = @_;

  my @binRanges;
  for(my $i=$min; $i<=$max; $i+=$binwidth){
    push @binRanges, [$i, $i + $binwidth, 
		      join('-', $i, $i+$binwidth)];
  }

  return \@binRanges;
}
