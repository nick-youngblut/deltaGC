#!/usr/bin/env perl
use strict;
use warnings;

=pod

=head1 NAME

deltaGC.pl -- Simulate template fragments and calculated GC & buoyant density

=head1 VERSION

This is version 0.0.1

=head1 USAGE

deltaGC.pl [options] -reads -genomes

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

=item -a[mplicon]

Amplicon reads instead of shotgun metagenome reads? [FALSE]

=item -sd <size_dist> | -size_distribution <size_dist>

Distribution from which fragment sizes are randomly selected.
Options: 'uniform','normal','exponential', 'poisson'

Default: size_dist.default

=for Euclid:
size_dist.type: string, size_dist eq 'uniform' || size_dist eq 'normal' || size_dist eq 'exponential' || size_dist eq 'poisson' || size_dist eq 'f' || size_dist eq 'skewed-normal'
size_dist.type.error: <size_dist> must be 'exponential', 'uniform', 'normal', 'skewed-normal', 'poisson', or 'f'
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

=item -skew[ness] <skewness>

Skewness of the normal distribution (negative values = left-skew)

Default: skewness.default

=for Euclid:
skewness.type: num
skewness.default: 0.75

=item -mu <mu>

mu parameter for poisson distribution.

Default: mu.default

=for Euclid:
mu.type: num > 0
mu.default: 1

=item -df <DFn> <DFd>

Numerator and denominator for the degrees of freedom in the F-distribution.

Default: DFn.default DFd.default

=for Euclid:
DFn.type: int > 0
DFn.default: 10
DFd.type: int > 0
DFd.default: 100

=item -p[rimer_buffer] <primer_buffer>

Fragment must contain this many bp before & after the 'amplicon' (specified by grinder).

Default: primer_buffer.default

=for Euclid:
primer_buffer.type: int >=0
primer_buffer.default: 70

=item -gap[_fraction] <gap_frac>

Fraction of DNA segment that can be composed of gaps (any character besides [ATGCatgc]). 
'NA' if cutoff not met.

Default: gap_frac.default

=for Euclid:
gap_frac.type: num >= 0
gap_frac.default: 0.05

=item -c_g[enomes]

Write out new version of the genome fasta file where all sequence
lines except the last are the same length for each entry? [FALSE]

=item -c_r[eads]

Write out new version of the read fasta file where all sequence
lines except the last are the same length for each entry? [FALSE]

=item -t[hreads] <threads>

Number of genomes to process in parallel.

Default: threads.default

=for Euclid:
threads.type: int >=0
threads.default: 0

=item -i[ndex]

Force the genome and read database files (*.index) to be rebuilt 
(saves time if not rebuilt)? [TRUE]

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

Simulate the DNA fragments containing the read 
(amplicon or shotgun) of interest and calculate
the GC content of the fragment & also the predicted
buoyant density of the fragment (assuming equilibrium
has been reached). 

=head2 Notes on size distributions

The f-distribution is scaled by 5 to be 0-1, which is
inverted (1-x) and multiplied by difference in range values
(max-min fragment length).

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
use Bio::DB::Fasta;
use Parallel::ForkManager;
use Term::ProgressBar;
use GC_dist qw/correct_fasta/;
use deltaGC qw/
get_frag_GC
get_genome_GC_stats
get_total_GC_stats
write_read_info
write_stats_summary
/;


#--- I/O error ---#
die "ERROR: min > max\n"
  if $ARGV{-range}{min_length} > $ARGV{-range}{max_length};

if(defined $ARGV{'-index'}){ $ARGV{'-index'} = 0; }
else{ $ARGV{'-index'} = 1; }

print STDERR "'-amplicon' not specified. Assuming shotgun reads, which affects how read position info is parsed\n"
  unless $ARGV{'-amplicon'};


# genome fasta loading
## correcting genome fasta if neede
if( $ARGV{-c_genomes} ){
  print STDERR "Correcting any line length inconsistencies in the genome fasta\n";
  $ARGV{-genomes} = correct_fasta($ARGV{-genomes});
}

## make genome database
if(! $ARGV{-index}){ print STDERR "Loading genome database index...\n" unless $ARGV{-quiet}; }
else{ print STDERR "Making genome database...\n" unless $ARGV{-quiet}; }

my $genome_db = Bio::DB::Fasta->new($ARGV{-genomes}, 
				    -reindex=>$ARGV{'-index'});

# parsing reads file
##correcting reads if needed
if( $ARGV{-c_reads} ){
  print STDERR "Correcting any line length inconsistencies in the reads fasta\n";
  $ARGV{-reads} = correct_fasta($ARGV{-reads});
}

## make read database
if(! $ARGV{-index}){ print STDERR "Loading read database index...\n" unless $ARGV{'--quiet'}; }
else{ print STDERR "Making read database...\n" unless $ARGV{'--quiet'}; }

sub make_my_id {
  my $desc = shift;
  $desc =~ /^>(\S+).*reference=(\S+)/;
  return join("__", $2, $1);
}
my $read_db = Bio::DB::Fasta->new($ARGV{-reads}, 
				  -makeid=>\&make_my_id, 
				  -reindex=>$ARGV{'-index'});
my @read_ids = $read_db->ids();

# creating fragments for each read & calculating GC
## header 
print join("\t", qw/genome scaffold read read_GC read_buoyant_density 
		    read_start read_length
		    fragment_GC fragment_buoyant_density 
		    fragment_start fragment_length/), "\n";

## forking
my $pm = Parallel::ForkManager->new( $ARGV{'-threads'} );
$pm->run_on_finish(
    sub { 
      my ($pid, $exit_code, $ident, $exit_signal, 
	  $core_dump, $ret_r) = @_;
      map{ print join("\t", @$_), "\n"} @$ret_r;
    }
  );


my @genomes = $genome_db->ids;
print STDERR "Creating fragments & calculating GC...\n" unless $ARGV{-quiet};
my $prog = Term::ProgressBar->new(scalar @genomes) unless $ARGV{-quiet};
#foreach my $genome (@genomes){
for my $i (0..$#genomes){
  # status
  $prog->update($i) unless $ARGV{-quiet};
  my $genome = $genomes[$i]; 

  # read_IDs for that genome
  my @spec_ids = grep /$genome\__/, @read_ids;

  $pm->start and next;
  # making random fragments & calculating GC
  my $ret_r = get_frag_GC($genome, 
			  \@spec_ids,
			  $genome_db, 
			  $read_db,	      
			  $ARGV{-amplicon},
			  $ARGV{-sd},
			  $ARGV{-range}{min_length},
			  $ARGV{-range}{max_length},
			  $ARGV{-mean},
			  $ARGV{-stdev},
			  $ARGV{-skewness},
			  $ARGV{-mu},
			  $ARGV{-df}{DFn},
			  $ARGV{-df}{DFd},
			  $ARGV{-primer_buffer},
			  $ARGV{-gap_fraction}
			 );

  $pm->finish(0, $ret_r);
}
$pm->wait_all_children;


