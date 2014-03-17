#!/usr/bin/env perl

use strict;
use warnings;
use Data::Dumper;
use Time::HiRes qw(time);
use MCE::Grep;

my $seq = "G" x 5000000 . "A" x 5000000;

my $time1 = time();
my $GC1 = calc_GC($seq);
my $time2 = time();
my $GC2 = calc_GC_MCE($seq);
my $time3 = time();

print $time2 - $time1, "\n";
print $time3 - $time2, "\n";
print "$GC1\t$GC2\n";


#--- subroutines ---#
sub calc_GC {
# calculting GC of a sequence
  my ($seq) = @_;

  # removing gaps #
  $seq =~ s/[-.]//g;

  # length
  my $len = length $seq;

  # GC
  my %GC = ( G => 0, C => 0 );
  while( $seq =~ /G/gi ){
    $GC{G}++
  }
  while( $seq =~ /C/gi ){
    $GC{C}++;
  }
  my $GC_sum = $GC{G} + $GC{C};
  my $GC_content = $GC_sum / $len * 100;

  return $GC_content, $len;
  }


sub calc_GC_MCE{

  my ($seq) = @_;

  # removing gaps #
  $seq =~ s/[-.]//g;

  # exploding
  my @seq = split //, $seq;

  my $GC_sum = scalar mce_grep{ $_ =~ /[GC]/i } @seq;
  
  my $len = length $seq;
  my $GC_content = $GC_sum / $len * 100;

  return $GC_content, $len;
}
