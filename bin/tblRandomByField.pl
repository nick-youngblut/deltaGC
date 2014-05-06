#!/usr/bin/env perl
use strict;
use warnings;

=pod

=head1 NAME

tblRandomByField.pl -- random sampling of rows from aggregated rows in a table

=head1 VERSION

This is version 0.0.1

=head1 USAGE

    tblRandomByField.pl [options]

=head1 REQUIRED ARGUMENTS

=over

=item <table.txt>

tab-delimited table file. ('-' if from STDIN)

=back

=head1 OPTIONS

=over

=item -c[olumn] <column>...

Column number(s) of the field(s) to select a random row from.

Default: column.default

=for Euclid:
column.type: int >= 1
column.default: 1

=item -header

Header in file? [FALSE]

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

This script simply selects a random row from
a table and writes it to STDOUT. Rows
are grouped by the specified column(s).

Make sure to sepecify a header if it exists in 
the table!

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
use Math::Random::MT::Perl;

#--- I/O error ---#
warn "Assuming no header in table file!\n"
  unless $ARGV{'-header'};

#--- MAIN ---#
  # filtering table
my $tbl_r = load_table($ARGV{'<table.txt>'}, 
		       $ARGV{'-column'},
		       $ARGV{'-header'});
select_rand($tbl_r);


#--- subroutines ---#
sub select_rand{
  my ($tbl_r) = @_;
  
  my $gen = Math::Random::MT::Perl->new(localtime);
   
  foreach my $cat (keys %$tbl_r){
    my $r = $gen->rand( $#{$tbl_r->{$cat}} );
    print ${$tbl_r->{$cat}}[$r], "\n";
  }
}

sub load_table{
# loading prokaryotes.txt table
  my ($file, $cols_r, $header) = @_;

  # 0-indexing columns
  map{ $_-- } @$cols_r;

  # file open
  my $fh;
  $file eq '-' ? $fh = \*STDIN : 
    open $fh, $file or die $!;  
  
  # table parse
  my %tbl;
  while(<$fh>){
    chomp;
    next if /^\s*$/;
  
    if($header && $.==1){
      print "$_\n";
      next;
    }
  

    my @l = split /\t/;
  
    # checking that cols exist
    map{ die "ERROR: cannot find column, ", $_ +1, 
	   " in line $.\n" unless defined $l[$_] } @$cols_r;

    push @{$tbl{join("|", @l[@$cols_r])}}, $_;
  }
  
  close $fh or die $!;

  #print Dumper %tbl; exit;
  return \%tbl;
}

