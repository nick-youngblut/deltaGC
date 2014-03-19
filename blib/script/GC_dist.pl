#!/usr/bin/env perl
use strict;
use warnings;

=pod

=head1 NAME

GC_dist.pl -- sliding window of GC around each grinder-produced read

=head1 VERSION

This is version 0.0.1

=head1 USAGE

    GC_dist.pl [options] > GC_byWindow.txt

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

=item -s[ize] <frag_size>

Fragment size (bp) around the read.

=for Euclid:
frag_size.type: int >0
frag_size.default: 10000

=item -w[indow] <window_size>

Size (bp) of sliding window around read.

=for Euclid:
window_size.type: int > 0
window_size.default: 100

=item -j[ump] <jump_size>

Size (bp) for sliding window jump around read.

=for Euclid:
jump_size.type: int >0
jump_size.default: 50 

=item -c_genomes

Write out new version of the genome fasta file where all sequence
lines except the last are the same length for each entry? [FALSE]

=item -c_reads

Write out new version of the read fasta file where all sequence
lines except the last are the same length for each entry? [FALSE]

=item -t[hreads] <threads>

Number of threads to use for GC calculations.

Default: threads.default

=for Euclid:
threads.type: +int
threads.default: 1

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

Determine how GC is distributed around
each read produced by grinder.

This is an exploratory program to 
help understand values produced by deltaMass.

The output table can be easily plotted in R
using ggplot.

=head2 WARNINGS

Delete the genome and read DB file (*index files)
if the script is killed before the DB construction
is completed.

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
#use FindBin;
#use lib "$FindBin::Bin/lib";
use local::lib;
use Bio::DB::Fasta;
use GC_dist qw/calc_GC/;
use List::MoreUtils qw/each_array/;
use GC_dist qw/ calc_frag_GC_window /;
#use MCE::Map;
use Parallel::ForkManager;


#--- I/O error ---#


#--- MAIN ---#
# genome fasta loading
## correcting genome fasta if needed
if( $ARGV{-c_genomes} ){
  print STDERR "Correcting any line length inconsistencies in the genome fasta\n";
  $ARGV{-genomes} = correct_fasta($ARGV{-genomes});
}

## make genome database
print STDERR "Making genome database...\n" unless $ARGV{'--quiet'};
my $genome_db = Bio::DB::Fasta->new($ARGV{-genomes});


# load read info
##correcting reads if needed
if( $ARGV{-c_reads} ){
  print STDERR "Correcting any line length inconsistencies in the reads fasta\n";
  $ARGV{-reads} = correct_fasta($ARGV{-reads});
}

## make read database
print STDERR "Making read database...\n" unless $ARGV{'--quiet'};
my $read_db = Bio::DB::Fasta->new($ARGV{-reads});
my @ids = $read_db->ids();

## getting read info
my @seq_obs = map{ $read_db->get_Seq_by_id($_) } @ids;
my @read_pos = map { parse_desc() } map { $_->desc } @seq_obs; # read genome-start-end


## combining read info by genome
# genome => read_id => [start|end] => value
my $ea = each_array( @ids, @read_pos);
my %reads;
while( my ($x,$y) = $ea->()){
  $reads{$$y[0]}{$x} = [$$y[1], $$y[2]];
}

## parsing genome fragments & calculating sliding window GC values
my $pm = Parallel::ForkManager->new($ARGV{-threads});

### finish statement
my %GC;
$pm->run_on_finish ( # called BEFORE the first call to start()
  sub {
    my ($pid, $exit_code, $ident, $exit_signal, $core_dump, $ref) = @_;
 
    # retrieve data structure from child
    if ( defined $ref  ) {  # children are not forced to send anything
      $GC{$pid} = $ref;
    } 
    else {  # problems occuring during storage or retrieval will throw a warning
      print STDERR "WARNING: No message received from child process $pid!\n";
    }
  }
);

### forking 
print STDERR "Calculating GC content by window...\n";
foreach my $genome (keys %reads){
  $pm->start and next;
  print STDERR "Processing genome: $genome\n" if $ARGV{'--debug'};
  my $gc_r = calc_frag_GC_window($genome, $reads{$genome}, $genome_db, $ARGV{'-size'}, $ARGV{'-window'}, $ARGV{'-jump'} );
  $pm->finish(0, $gc_r);
}
$pm->wait_all_children;


### writing output 
print join("\t", qw/genome scaf readID read_num pos_local pos_global GC/), "\n";
foreach my $pid (keys %GC){
  foreach my $genome (keys %{$GC{$pid}}){
    # getting description
    my $desc = $genome_db->header($genome);

    # writing out windows for each read
    my $read_num = 0;
    foreach my $read (keys %{$GC{$pid}{$genome}}){
      $read_num++;
      foreach my $start (sort {$a<=>$b} keys %{$GC{$pid}{$genome}{$read}}){
	print join("\t",$desc, $genome, $read, $read_num, $start, 
		   $GC{$pid}{$genome}{$read}{$start}{pos},
		   $GC{$pid}{$genome}{$read}{$start}{GC}), "\n";
      }
    }
  }
}

#--- subroutines ---#
use Text::ParseWords;
sub parse_desc{
  s/>//;
  my %vals = quotewords("=| ", 0, $_);
  
  map{ die "ERROR: cannot find '$_' in read!\n" 
	 unless exists $vals{$_} } qw/amplicon reference/;

  $vals{amplicon} =~ /(\d+)\.\.(\d+)/;

  return [$vals{reference}, $1, $2];
}

