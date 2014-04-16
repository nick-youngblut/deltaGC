#!/usr/bin/env perl
use strict;
use warnings;

=pod

=head1 NAME

rnammer_tandem_rrn.pl -- find adjacent rrn operons using gff files from rnammer

=head1 VERSION

This is version 0.0.1

=head1 USAGE

    rnammer_tandem_rrn.pl [options]

=head1 REQUIRED ARGUMENTS

=over

=item <gff_file_list>

List of gff files produced by rnammer ('-' if STDIN).

=back

=head1 OPTIONS

=over

=item -m[[ax][_distance]] <d>

Max distance between 16S rRNA and adjacent 23S rRNA genes (bp).

Default: d.default

=for Euclid:
d.type: int > 0
d.default: 5000

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
use rnammer_tandem_rrn qw/
load_gff_file_list
gff2itree
check_tandem/;

#--- I/O error ---#

#--- MAIN ---#
# getting list of gff files from STDIN; making array
my $gff_list_r = load_gff_file_list($ARGV{'<gff_file_list>'});

# loading gff files and placing features into interval trees (by genome)
my ($itrees_r, $ssu_loc_r) = gff2itree($gff_list_r);

# foreach genome, determine if 23S rRNA gene is upstream of 16S rRNA gene
## load table of counts 
my $summary_r = check_tandem($itrees_r,$ssu_loc_r, $ARGV{'-m'});

# write summary table
print join("\t", qw/file n_copy n_tandem frac_tandem/),"\n";
foreach my $file (keys %$summary_r){
  $summary_r->{$file}{tandem} = 0 unless exists $summary_r->{$file}{tandem};
  print join("\t", $file, 
	     @{$summary_r->{$file}}{qw/copy tandem/},
	     $summary_r->{$file}{tandem} / $summary_r->{$file}{copy}
	    ), "\n";
}
