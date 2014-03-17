#!/usr/bin/env perl
use strict;
use warnings;

=pod

=head1 NAME

deltaMass.pl -- basic description here

=head1 VERSION

This is version 0.0.1

=head1 USAGE

deltaMass.pl [options] -reads -genomes

=head1 REQUIRED ARGUMENTS

=over

=item -r[eads] <read_file>

Reads output by grinder

=for Euclid:
read_file.type: readable

=item -g[enomes] <genome_file>

Fasta file of genomes used as '-reference_file' in grinder

=for Euclid:
genome_file.type: readable

=back

=head1 OPTIONS

=over

=item -sd <size_dist> | -size_distribution <size_dist>

Distribution from which fragment sizes are randomly selected.
Options: 'uniform','normal','exponential'

Default: size_dist.default

=for Euclid:
size_dist.type: string, size_dist eq 'uniform' || size_dist eq 'normal' || size_dist eq 'exponential'
size_dist.type.error: <size_dist> must be 'uniform', 'normal', or 'exponential'
size_dist.default: 'exponential'

=item -range <min_length>-<max_length>

Minimum & maximum fragment length.

Default: min_length.default max_length.default

=for Euclid:
min_length.type: int > 0
min_length.default: 4000
max_length.type: int > 0
max_length.default: 100000

=item -mean <mean_length>

Mean fragment length.

Default: mean_length.default

=for Euclid:
mean_length.type: int > 0
mean_length.default: 4000

=item -stdev <stdev>

Fragment length standard deviation.

Default: stdev.default

=for Euclid:
stdev.type: int >=0
stdev.default: 100

=item -p[rimer_buffer] <primer_buffer>

Fragment must contain this many bp before & after the 'amplicon' (specified by grinder).

Default: primer_buffer.default

=for Euclid:
primer_buffer.type: int >=0
primer_buffer.default: 70


=item -c[orrect]

Write out new version of the genome fasta file where all sequence
lines except the last are the same length for each entry? [FALSE]

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

The rest of the documentation

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

use Data::Dumper;
use Getopt::Euclid;
use FindBin qw/$Bin/;
use lib "$Bin";
use Bio::DB::Fasta;
use deltaMass qw/
load_genome_fasta
load_reads
get_frag_GC
get_genome_GC_stats
get_total_GC_stats
write_read_info
write_stats_summary
/;


#--- I/O error ---#
die "ERROR: min > max\n"
  if $ARGV{-range}{min_length} > $ARGV{-range}{max_length};


# correcting genome fasta if neede
if( $ARGV{-correct} ){
  print STDERR "Correcting any line length inconsistencies in the genome fasta\n";
  $ARGV{-genomes} = correct_fasta($ARGV{-genomes});
}

# make genome database
print STDERR "Make genome database...\n" unless $ARGV{-quiet};

#my $genome_seqs_r = load_genome_fasta($ARGV{-genomes});
my $db = Bio::DB::Fasta->new($ARGV{-genomes});

# parsing reads file
print STDERR "Loading reads...\n" unless $ARGV{-quiet};
my $reads_r = load_reads($ARGV{-reads});

# creating fragments for each read & calculating GC
print STDERR "Creating random fragments & calculating GC...\n" unless $ARGV{-quiet};
 get_frag_GC($db, $reads_r,
             $ARGV{-sd},
             $ARGV{-range}{min_length},
             $ARGV{-range}{max_length},
             $ARGV{-mean},
             $ARGV{-stdev},
             $ARGV{-primer_buffer});


# writing all read info to file if debug
write_read_info($reads_r, $db);



#--- using R or python for stats calculation ---#

# calculate variance & CI
## each genome
# print STDERR "Calculating delta-GC statistics for each genome...\n" unless $ARGV{-quiet};
# MCE->new(
# 	 chunk_size => 1,
# 	 max_workers => $ARGV{-c} ,
# 	 user_func => \&get_genome_GC_stats,
# 	 user_args => { 'reads_r' => $reads_r }
# 	);
# my %stats;
# MCE->process( [keys %$reads_r], { gather =>\%stats } );

## all genomes combined
#print STDERR "Calculating delta-GC statistics for all genomes combined...\n" unless $ARGV{-quiet};
#get_total_GC_stats($reads_r, \%stats);

# writing out summary
#write_stats_summary(\%stats);
