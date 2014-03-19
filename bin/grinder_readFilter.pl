#!/usr/bin/env perl
use strict;
use warnings;

=pod

=head1 NAME

grinder_readFiler.pl -- filter out (nearly) overlapping reads produced by grinder

=head1 VERSION

This is version 0.0.1

=head1 USAGE

    grinder_readFiler.pl [options] 

=head1 REQUIRED ARGUMENTS

=over

=item -r[eads] <read_file>

Reads output by grinder

=for Euclid:
read_file.type: readable

=back

=head1 OPTIONS

=over

=item -w[indow] <window>

Amplicons located within '-window' of each
other and will be filtered out with one chosen
randomly.

=for Euclid:
window.type: +int
window.default: 100

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

Dereplicate reads (primarily amplicons) produced by grinder 
before downstream analyses.

Reads are filtered out if they are overlapping or 
nearly overlapping, with one of the overlapping reads
randomly chosen as the representative sequence.

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
use FindBin;
use lib "$FindBin::Bin/lib";
use Bio::SeqIO;
use Bio::DB::Fasta;
use GC_dist qw/calc_GC/;
use MCE::Map;
use GC_dist qw/correct_fasta/;
use readFilter qw/
load_read_itrees
screen_reads/;


#--- I/O error ---#

#--- MAIN ---#
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
my @read_pos = mce_map { parse_desc() } map { $_->desc } @seq_obs; # read genome-start-end

## loading itrees
my $itrees_r = load_read_itrees(\@ids, \@read_pos);

## screening reads
my $keep_ids_r = screen_reads(\@ids, \@read_pos, $itrees_r, $ARGV{-window}); 



## writing filtered reads
my $outseq = Bio::SeqIO->new( -fh => \*STDOUT,
			      -format => 'fasta');
foreach my $seqo ( map { $read_db->get_Seq_by_id($_) } sort{$a<=>$b} keys %$keep_ids_r){
  $outseq->write_seq($seqo);
}


#--- subroutines ---#
use Text::ParseWords;
sub parse_desc{
  s/>(\d+) //;
  my %vals = quotewords("=| ", 0, $_);
  
  map{ die "ERROR: cannot find '$_' in read!\n" 
	 unless exists $vals{$_} } qw/amplicon reference/;

  $vals{amplicon} =~ /(\d+)\.\.(\d+)/;

  return [$vals{reference}, $1, $2];
}

