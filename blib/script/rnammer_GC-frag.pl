#!/usr/bin/env perl

=pod

=head1 NAME

rnammer_GC-frag.pl -- get GC of sub-fragment of rRNA gene & fragment encompassing rRNA gene 

=head1 SYNOPSIS

rnammer_GC-frag.pl [flags] > out.txt

=head2 Required flags

=over

=item -genome

Fasta file of genome used for rnammer analysis. 

=item -seq

Fasta file output by rnammer.

=back

=head2 Optional flags

=over

=item -start  <int>

Start of the sub-fragment of the rRNA gene. [515]

=item -end  <int>

End of the sub-fragment of the rRNA gene. [927]

=item -frag_size  <int>

Size of the fragment encompassing the rRNA gene. [4000]

=item -verbose  <bool>

Verbose output. [FALSE]

=item -h  <bool>

Print this help message & exit. [FALSE]

=back

=head2 For more information:

perldoc rnammer_GC-frag.pl

=head1 DESCRIPTION

Quick and dirty calculating of how
diagnositic the GC of an rRNA amplicon 
is to the genomic fragment it was amplified
from (the fragment that ran through the CsCl
gradient).

=head2 output

=over

=item sub-frag_GC-content

=item large-frag_GC-content

=item scaffold

=item start

=item end

=item amplicon_length

=item rnammer_score

=item molecule

=item ssu_start

=item ssu_end

=item frag_size

=back

=head1 EXAMPLES

=head2 Basic usage:

rnammer_GC-frag.pl -genome ecoli.fasta -seq rnammer_out.fasta  > output

=head1 AUTHOR

Nick Youngblut <ndy2@cornell.edu>

=head1 AVAILABILITY

email me

=head1 COPYRIGHT

Copyright 2010, 2011
This software is licensed under the terms of the GPLv3

=cut


#--- modules ---#
use strict;
use warnings;
use Pod::Usage;
use Data::Dumper;
use Getopt::Long;
use FindBin qw/$Bin/;
use lib "$Bin";


#--- args/flags ---#
pod2usage("$0: No files given.") if ((@ARGV == 0) && (-t STDIN));

my ($verbose_b, $genome_in, $seq_in, $gff_in);
my $ssu_start = 515;
my $ssu_end = 927;
my $frag_size = 4000; # 4 kb frag size
GetOptions(
#	   "gff=s" => \$gff_in,
	   "genome=s" => \$genome_in,
	   "seq=s" => \$seq_in,
	   "start=i" => \$ssu_start,
	   "end=i" => \$ssu_end,
	   "frag_size=i" => \$frag_size,
	   "verbose" => \$verbose_b,
	   "help|?" => \&pod2usage # Help
	  );

#--- I/O error ---#
#check_file($gff_in, '-gff');
check_file($genome_in, '-genome');
check_file($seq_in, '-seq');

#--- setting defaults ---#
#use lib "$FindBin::RealBin/../lib";
#use lib "$FindBin::RealBin/../lib/perl5";

#--- MAIN ---#
# parse gff
#my $gff_r = parse_gff($gff_in);

# get GC of 16S gene 
my %res;
get_ssu_GC($seq_in, \%res, $ssu_start, $ssu_end);

# get GC of fragment
my $genome = deltaGC::load_fasta($genome_in);
get_fragment_GC($genome, \%res, $frag_size);

# write output table
write_output(\%res, $ssu_start, $ssu_end, $frag_size);

#--- Subroutines ---#
sub write_output{
  my ($res_r, $ssu_start, $ssu_end, $frag_size) = @_;

  foreach my $ssu (keys %$res_r){
    print join("\t", 
	       $res_r->{$ssu}{'gene_frag'}{gc_content},
	       $res_r->{$ssu}{'large_frag'}{gc_content},
	       $res_r->{$ssu}{'gene_frag'}{scaffold},
	       $res_r->{$ssu}{'gene_frag'}{start},
	       $res_r->{$ssu}{'gene_frag'}{end},	       
	       $res_r->{$ssu}{'gene_frag'}{length},
	       $res_r->{$ssu}{'gene_frag'}{score},
	       $res_r->{$ssu}{'gene_frag'}{molecule},
	       $ssu_start,
	       $ssu_end,
	       $frag_size), "\n";
  }

}

sub get_fragment_GC{
# getting GC of the fragment
  my ($genome, $res_r, $frag_size) = @_;

 # print Dumper $genome; exit;
  foreach my $ssuFrag (keys %$res_r){
    # sanity check 
    die "ERROR: cannot find sequence: $res_r->{$ssuFrag}{'gene_frag'}{'scaffold'} in genome\n"
      unless exists $genome->{ $res_r->{$ssuFrag}{gene_frag}{scaffold} };
    # determining fragment size to parse
    my $center = int(abs( $res_r->{$ssuFrag}{gene_frag}{end} - $res_r->{$ssuFrag}{gene_frag}{start} )) / 2;
    my $trim_start = $res_r->{$ssuFrag}{gene_frag}{end} - ($center + int($frag_size * 0.5)); # geneFrag_end - (center_bpFromStart - $frag_size * 1/2)
    my $largeFrag = substr($genome->{ $res_r->{$ssuFrag}{gene_frag}{scaffold} }, $trim_start - 1, $frag_size);
       
    # calc gc
    my @gc = calc_GC($largeFrag);
    load_gc(\@gc, $ssuFrag, $res_r, 'large_frag');
   
    #print Dumper $largeFrag; exit;
  }

#print Dumper %$res_r; exit;
}

sub get_ssu_GC{
# getting the GC of each identified SSU sequences
  my ($seq_in, $res_r, $ssu_start, $ssu_end) = @_;

  use rnammer_GCfrag qw/ calc_GC /;
  open IN, $seq_in or die $!;

  my (%fasta, $tmpkey);

  my %GC;
  while(<IN>){
    chomp;
    next if /^\s*$/;
    
    if(eof(IN)){
      $fasta{$tmpkey} .= $_;
      my $seq_trim = trim_seq($fasta{$tmpkey}, $ssu_start, $ssu_end);
      my @gc = calc_GC($seq_trim);
      load_gc(\@gc, $tmpkey, $res_r, 'gene_frag');
    }
    elsif(/^>/){
      if(defined $tmpkey){
	my $seq_trim = trim_seq($fasta{$tmpkey}, $ssu_start, $ssu_end);
	my @gc = calc_GC($seq_trim);
	load_gc(\@gc, $tmpkey, $res_r, 'gene_frag');
	%fasta = ();
      }
      
      # parsing rnammer seq name
      ($tmpkey = $_) =~ s/^>//;
      my @name = parse_rnammer_seq_name($tmpkey, $res_r, 'gene_frag');
      
      $fasta{$tmpkey} = "";
    }
    else{
      $fasta{$tmpkey} .= $_;
    }
  }
  close IN;

  #print Dumper %$res_r; exit;
}

sub trim_seq{
  my ($seq, $start, $end) = @_;
  
  die "ERROR: start ($start) > end ($end)\n"
    if $start > $end;

  my $seq_len = length $seq;

  my $frag_len = $end - $start + 1;
  $frag_len = $seq_len if $frag_len > $seq_len; # boundry issues

  my $seq_trim = substr($seq, $start - 1, $frag_len);

  #print Dumper length $seq_trim;
  #print Dumper $seq_trim; exit;
  return $seq_trim;
}

sub load_gc{
# loading gc into hash of all results
  my ($gc_r, $name, $res_r, $cat) = @_;

  die "ERROR: calc_gc did not return values\n"
    unless scalar @$gc_r == 3;

  $res_r->{$name}{$cat}{gc_content} = $$gc_r[0];
  $res_r->{$name}{$cat}{gc_sum} = $$gc_r[1];
  $res_r->{$name}{$cat}{length} = $$gc_r[2];
}

sub parse_rnammer_seq_name{
  my ($name, $res_r, $cat) = @_;

  my @j = split / \//, $name;
  $j[1] =~ s/.*?=//;
  $j[2] =~ s/.*=//;
  my @jj = split /_/, $j[0];
  
  my @start_end = split /-/, $jj[2];

  $res_r->{$name}{$cat}{scaffold} = $jj[1];
  $res_r->{$name}{$cat}{start} = $start_end[0];
  $res_r->{$name}{$cat}{end} = $start_end[1];
  $res_r->{$name}{$cat}{molecule} = $j[1];
  $res_r->{$name}{$cat}{score} = $j[2];
  $res_r->{$name}{$cat}{strand} = $jj[3];

  #print Dumper %$res_r; exit;
}

sub parse_gff{
# parsing the input gff file
  my ($gff_in) = @_;

  open IN, $gff_in or die $!;

  my @gff;
  while(<IN>){
    chomp;
    next if /^\s*$/;
    next if /^#/; # no comments

    my @l = split /\t/;
    die "ERROR: line $. in $gff_in is not formatted correctly\n"
      unless scalar @l >= 9;

    push @gff, \@l;
  }

  #print Dumper @gff; exit;
  return \@gff;
}

sub check_file{
  my ($file, $names) = @_;
  die "ERROR: provide '$names'\n"
    unless defined $file;
  die "ERROR: cannot find $file\n"
    unless -e $file;
}
	   
