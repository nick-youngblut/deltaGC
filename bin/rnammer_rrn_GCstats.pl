#!/usr/bin/env perl
use strict;
use warnings;

=pod

=head1 NAME

rnammer_rrn_GCstats.pl -- calculating GC content stats for ssu, lsu, tsu & ITS features

=head1 VERSION

This is version 0.0.1

=head1 USAGE

    rnammer_rrn_GCstats.pl [options]

=head1 REQUIRED ARGUMENTS

=over

=item -g[ff] <gff_file>

gff file of rnammer output.

=for Euclid:
gff_file.type: in

=item -f[asta] <fasta_dir>

Directory containing fasta files for each genome.

=for Euclid:
fasta_dir.type: str

=back

=head1 OPTIONS

=over

=item -s[core] <16s_rRNA> <23s_rRNA> <5s_rRNA>

rnammer score cutoffs for each feature (16s_rRNA, 23s_rRNA, 5s_rRNA).

Defaults: 16s_rRNA.default 23s_rRNA.default 5s_rRNA.default

=for Euclid:
16s_rRNA.type: int > 0
23s_rRNA.type: int > 0
5s_rRNA.type: int > 0
16s_rRNA.default: 1250
23s_rRNA.default: 2500
5s_rRNA.default: 65

=item -I[[TS][_length]] <i>

The max ITS length allowed

Default: i.default

=for Euclid:
i.type: int > 0
i.default: 2000

=item -gap[_fraction] <gap_frac>

Fraction of DNA segment that can be composed of gaps ('-').
GC_content = 'NA' if cutoff not met.

Default: gap_frac.default

=for Euclid:
gap_frac.type: num >= 0
gap_frac.default: 0.05

=item -r[ound] <r>

Decimal to round output values to.

Default: r.default

=for Euclid:
r.type: int >= 0
r.default: 3

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
use Statistics::Descriptive; 
use rnammer_rrn_GCstats qw/
			    parse_rnammer_gff
			    pairwise_dist_ssu_lsu
			    determine_ITS
			    load_fasta
			    substr_frags
			  /;
use deltaGC qw/ calc_GC /;

#--- I/O error ---#
die "ERROR: cannot find $ARGV{'-fasta'}\n"
  unless -d $ARGV{'-fasta'};

my $round = $ARGV{'-round'};
$round = join("", '%.', $round, 'f');


#--- MAIN ---#
# parsing rnammer hits
my $gff_r = parse_rnammer_gff($ARGV{'-gff'}, $ARGV{'-score'});

# determining ITS regions
my $pdists_r = pairwise_dist_ssu_lsu($gff_r);
determine_ITS($gff_r, $pdists_r, $ARGV{'-I'});

# selecting sequences from fasta & calculating GC
my %GC; # {gene} => [GC]
foreach my $genome (keys %$gff_r){
  print STDERR "processing genome: $genome\n" if $ARGV{'--debug'};
  # getting fragments
  my $fasta_r = load_fasta($genome, $ARGV{'-fasta'});
  my $frags_r = substr_frags($gff_r, $genome, $fasta_r);   

  # calc GC for each fragment 
  my %GC_loc;
  foreach my $gene (keys %$frags_r){
    $GC{$gene} = Statistics::Descriptive::Full->new()
      unless exists $GC{$gene};
    $GC_loc{$gene} =  Statistics::Descriptive::Full->new()
      unless exists $GC_loc{$gene};
    foreach my $seq (@{$frags_r->{$gene}}){
      my ($gc, $len) = calc_GC($seq, $ARGV{'-gap'});
      $GC{$gene}->add_data( $gc );
      $GC_loc{$gene}->add_data( $gc );
    }
  }
  
  # genome GC stats
  if( $ARGV{'--debug'} > 1 ){
    foreach my $gene (keys %GC_loc){
      print join("\t", $genome, $gene,
		 sprintf($round, $GC_loc{$gene}->min()),
		 sprintf($round, $GC_loc{$gene}->mean()),
		 sprintf($round, $GC_loc{$gene}->median()),
		 sprintf($round, $GC_loc{$gene}->max()),
		 sprintf($round, $GC_loc{$gene}->standard_deviation())
		), "\n";
    }
  }
}

# output
print join("\t", qw/gene min mean median max stdev/),"\n";
foreach my $gene (keys %GC){
  print join("\t", $gene,
	     sprintf($round, $GC{$gene}->min()),
	     sprintf($round, $GC{$gene}->mean()),
    	     sprintf($round, $GC{$gene}->median()),
	     sprintf($round, $GC{$gene}->max()),
	     sprintf($round, $GC{$gene}->standard_deviation())
	     ), "\n";
}
