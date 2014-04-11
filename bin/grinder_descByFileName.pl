#!/usr/bin/env perl
use strict;
use warnings;

=pod

=head1 NAME

grinder_descByFileName.pl -- grinder read descriptions changed from chromosome to genome fasta file name

=head1 VERSION

This is version 0.0.1

=head1 USAGE

    grinder_descByFileName.pl [options] 

=head1 REQUIRED ARGUMENTS

=item <dir> 

Directory containing all of the read fasta files produced by grinder ('*-reads.fa')

=head1 OPTIONS

=over

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

The 'description="taxon"' portion of sequence headers in
the grinder read fasta files as based on the sequence name
in the genome fasta. This is a problem if genomes have
>1 chromosome or scaffold. This script gives all reads
for the same taxon the same taxon name (based on
read fasta file name). The main requirement is that
each taxon be processed separately by grinder 
(easy to do with a for loop or xargs).

The reads from each file will be written to STDOUT.

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

#--- I/O error ---#
die "ERROR: cannot find $ARGV{'<dir>'}\n"
  unless -d $ARGV{'<dir>'};

#--- MAIN ---#
opendir IN, $ARGV{'<dir>'} or die $!;
my @files = grep(/^[^.].+-reads.fa/, readdir IN);
closedir IN;

foreach my $file (@files){
  open IN, "$ARGV{'<dir>'}/$file" or die $!;
  (my $new_name = $file) =~ s/-reads.fa//;
  $new_name =~ s/(\.fasta|\.fna)$//;
  while(<IN>){
    if(/^>/){
      die "ERROR: cannot find 'description=' in $file -> line $.\n"
	unless /description=/;
      s/(description)="([^"]+)"/$1="$new_name"/;
    }
    print;
  }
}
