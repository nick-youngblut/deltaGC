#!/usr/bin/env perl
use strict;
use warnings;

=pod

=head1 NAME

GC_dist.pl -- fast calculation of GC content for all sequence in a fasta file

=head1 VERSION

This is version 0.0.1

=head1 USAGE

    GC_dist.pl [options] 

=head1 REQUIRED ARGUMENTS

=over

=item <fasta_file>

Fasta file ('-' if provided via STDIN).

=back

=head1 OPTIONS

=item -N[_remove]

Remove 'N' characters from G+C calculation. [TRUE]

=item -w[orkers] <w>

Number of workers used by MCE.

Default: 'auto'

=for Euclid:
w.type: int > 0

=item -c[hunk] <c>

Chunk size used by MCE.

Default: 'auto'

=for Euclid:
c.type: int > 0

=item --debug [<log_level>]

Set the log level. Default is log_level.default but if you provide --debug,
then it is log_level.opt_default.

=for Euclid:
    log_level.type:        int
    log_level.default:     0
    log_level.opt_default: 1

=item --quiet

=item --version

=item --usage

=item --help

=item --man

Print the usual program information

=back

=head1 DESCRIPTION

Determine G+C content stats for a genome fasta or
collection of other sequences.

MCE is used to quickly determine G+C content by
using multiple processors (workers). 

Ambiguous IUPAC nucleotide codes will be included
in GC calculations based the the fraction of 'G' & 'C'
in the ambiguous character. See 'Warnings' about 'N' characters.

=head2 Warnings

If 'N' characters signify gaps, remove them from the
genome/read fasta files (use '-N' flag)!

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
use MCE::Map;


#--- I/O error ---#
my $fh;
$ARGV{'<fasta_file>'} eq '-' ? $fh = \*STDIN :
  open $fh, $ARGV{'<fasta_file>'} or die $!;

$ARGV{'-w'} = 'auto' unless defined $ARGV{'-w'};
$ARGV{'-c'} = 'auto' unless defined $ARGV{'-c'};


MCE::Map::init {
  chunk_size => $ARGV{'-c'}, 
    max_workers => $ARGV{'-w'}
};

#--- MAIN ---#
sub calc_gc{
  return [0,0] if /^>/;
  
  s/N//g unless defined $ARGV{'-N'};

  # upcase
  tr/a-z/A-Z/;

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
  
  my $q = join("", "[", keys %score, "]");
  $q = qr/$q/;
  my $GC_sum = 0;
  $GC_sum += $score{$1} while /($q)/g;

  return [$GC_sum, length $_];
}

my @vals = mce_map_f{ calc_gc($_) } $fh;

my($gc_sum, $len_sum) = (0,0);
foreach (@vals){
  $gc_sum += $_->[0];
  $len_sum += $_->[1];
}

print join("\t", $ARGV{'<fasta_file>'}, 
	   $gc_sum / $len_sum), "\n";
