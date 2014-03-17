package deltaMass;

use 5.006;
use strict;
use warnings;

=head1 NAME

deltaMass - scripts for running deltaMass.pl

=head1 VERSION

Version 0.01

=cut

our $VERSION = '0.01';


=head1 SYNOPSIS

=head1 EXPORT


=head1 SUBROUTINES/METHODS

=cut

use base 'Exporter';
our @EXPORT_OK = '';
use Carp qw/ confess carp /;
use Data::Dumper;
use Text::ParseWords;
use Bio::SeqIO;
#use Statistics::Basic;
#use Statistics::PointEstimation;
use Math::Random qw/
random_uniform_integer
random_uniform
random_normal
random_exponential/;


=head2 correct_fasta

Making sure sequence lines for each entry are the same length.

=cut

sub correct_fasta{
  my $genome_file = shift;

  (my $genome_out = $genome_file) =~ s/\.[^.]+$/_fix.fasta/;

  my $seqin = Bio::SeqIO->new(-file => $genome_file, -format => 'fasta');
  my $seqout = Bio::SeqIO->new(-file => ">$genome_out", -format => 'fasta');

  while(my $seq = $seqin->next_seq){
    $seqout->write_seq($seq);
  }

  print STDERR "Corrected genome file written: '$genome_out'\n";
  return $genome_out;
}

=head2 calc_GC

=cut

push @EXPORT_OK, 'calc_GC';

sub calc_GC {
# calculting GC of a sequence
  my ($seq) = @_;

  confess "ERROR: no sequences provided for calculating GC!\n" 
    unless defined $seq;

  # removing gaps #
  $seq =~ s/[-.]//g;

  # length
  my $len = length $seq;

  # GC
  my %GC = ( G => 0, C => 0 );
  while( $seq =~ /G/gi ){
    $GC{G}++;
  }
  while( $seq =~ /C/gi ){
    $GC{C}++;
  }
  my $GC_sum = $GC{G} + $GC{C};
  my $GC_content = $GC_sum / $len * 100;

  return $GC_content, $len;
}

=head2 load_genome_fasta

=cut

push @EXPORT_OK, 'load_genome_fasta';

sub load_genome_fasta {
  my $fasta_file = shift;
  open IN, $fasta_file or die $!;

  # loading fasta 
  my (%fasta, $tmpkey);
  while(<IN>){
    chomp;
    s/#.+//;
    next if  /^\s*$/;
    if(/>.+/){
      $tmpkey = $_;# changing key
      $tmpkey =~ s/>| .+//g; # nothing after space kept
      $fasta{$tmpkey}{seq} = "";
    }
    else{$fasta{$tmpkey}{seq} .= $_; }
  }
  close IN;

  # getting genome lengths
  map{ $fasta{$_}{len} = length $fasta{$_}{seq} } keys %fasta;

  #print Dumper %fasta; exit;
  return \%fasta;
}

=head2 load_reads

Loading reads produced by grinder.
Fasta format, but headers need to be parsed

=cut

push @EXPORT_OK, 'load_reads';

sub load_reads{
  my $reads_file = shift;
  
  open IN, $reads_file or die $!;

  my (%fasta, $ref, $Uid);
  while(<IN>){
    chomp;

    if(eof(IN)){
      $fasta{$ref}{$Uid}{seq} .= $_;
      ($fasta{$ref}{$Uid}{amp_GC}, $fasta{$ref}{$Uid}{amp_len}) 
	= calc_GC($fasta{$ref}{$Uid}{seq});
      delete $fasta{$ref}{$Uid}{seq} or die $!;
    }
    elsif(/^>/){
      # calculating GC on last sequence
      if( defined $ref ){
	($fasta{$ref}{$Uid}{amp_GC}, $fasta{$ref}{$Uid}{amp_len}) 
	  = calc_GC($fasta{$ref}{$Uid}{seq});
	delete $fasta{$ref}{$Uid}{seq} or die $!;
      }

      # loading new sequence
      my @l = split />| /, $_, 3;
      my %vals = quotewords("=| ", 0,  $l[2]);
      
      $Uid = $l[1];
      $ref = $vals{reference};
      
      ## amplicon start-end
      my ($amp_start, $amp_end, $amp_strand);
      if( $vals{amplicon} =~ /complement\((.+)\)/){ 
	($amp_end, $amp_start) = split /\.\./, $1;
	$amp_strand = -1;
      }
      else{
	($amp_start, $amp_end) = split /\.\./, $vals{amplicon};
	$amp_strand = 1;
      }

      ## sanity check
      die "ERROR: read_start at 0. 1-indexing!\n"
	if $amp_start == 0;
      
      # loading postion info
      $fasta{$ref}{$Uid}{seq} = "";
      $fasta{$ref}{$Uid}{amp_start} = $amp_start;
      $fasta{$ref}{$Uid}{amp_end} = $amp_end;
     #$fasta{$ref}{$Uid}{amp_strand} = $amp_strand;
     #$fasta{$ref}{$Uid}{desc} = $vals{description};
    }
    else{ # just part of the sequence
      next if /^\s*$/;
      $fasta{$ref}{$Uid}{seq} .= $_;
    }
  }

  close IN or die $!;

  #print Dumper %fasta; exit;
  return \%fasta;
}

=head2 get_frag_GC

Creating fragments for each read & calculating GC on each.
Called with MCE.
Return (gather): a hash that replicates input hash along with frag_GC, frag_len, & frag_start

=cut

push @EXPORT_OK, 'get_frag_GC';

sub get_frag_GC_MCE{
  my $genome = $_;
  my ($self) = @_;

  my $reads_r = $self->user_args->{reads_r};
  my $db = $self->user_args->{db};
  my $ARGV = $self->user_args->{argv};

  #print Dumper $argv; exit;

  # inializing args
  my ($size_dist, $frag_min, $frag_max, $mean, $stdev, $primer_buffer) =
     ($ARGV{-sd},
      $ARGV{-range}{min_length},
      $ARGV{-range}{max_length},
      $ARGV{-mean},
      $ARGV{-stdev},
      $ARGV{-primer_buffer},
      $ARGV{-c});
  
  # sanity checks
  ## values provided
  map{ confess "ERROR: argument missing: $!\n" unless defined $_ } 
    ($size_dist, $frag_min, $frag_max, $mean, $stdev, $primer_buffer);

  # checking for genome
  my $genome_len = $db->length($genome);
  die "ERROR: cannot find '$genome' in genomes fasta\n"
      unless defined $genome_len;

  # interating through fragments
  my %res;
  foreach my $Uid (keys %{$reads_r->{$genome}}){
    my $amp_start = $reads_r->{$genome}{$Uid}{amp_start};
    my $amp_end = $reads_r->{$genome}{$Uid}{amp_end};
    my $amp_len = $reads_r->{$genome}{$Uid}{amp_len};
    
    # fragment size
    my $frag_size;
    if($size_dist eq 'uniform'){
      while(1){
	$frag_size = random_uniform_integer(1, $frag_min, $frag_max);
	last if $frag_size >= $frag_min && $frag_size <= $frag_max;
      }
    }
    elsif($size_dist eq 'normal'){
      while(1){
	$frag_size = int random_normal(1, $mean, $stdev);
	last if $frag_size >= $frag_min && $frag_size <= $frag_max;
      }
    }
    elsif($size_dist eq 'exponential'){
      while(1){
	$frag_size = int random_exponential(1, $mean);
	last if $frag_size >= $frag_min && $frag_size <= $frag_max;
      }
    }
    else{ confess "ERROR: do not recognize size distribution\n"; }
            
    # determine fragment start-end based on amplicon start-end
    ## amplicon center postion
    my $amp_center = int($amp_start + abs($amp_len * 0.5));
    ## fragment_start = amp_center - (frag_size * x); x = random draw from unifrom distribution 0:1
    my $x = random_uniform();
    ### start_floor = (primer_buffer + 0.5*amp_len) / frag_size
    my $floor  = ($primer_buffer + 0.5 * $amp_len) / $frag_size; 
    $x = $floor if $x < $floor;
    my $ceiling = 1 - $floor;                               # ceiling = 1 - floor
    $x = $ceiling if $x > $ceiling;
    my $frag_start = int( $amp_center - ($frag_size * $x) );

    ### sanity check 
    carp "WARNING: frag_start is too far from amplicon!\n\tfrag_start: $frag_start, amp_end: $amp_end, primer_buffer: $primer_buffer, frag_size: $frag_size, amp_center: $amp_center, x: $x\n"
      if $frag_start - ($amp_end + $primer_buffer) > $frag_size;

    print Dumper $db->seq($genome, $frag_start, $frag_start + $frag_size);
      
    # storing amp_start & fragment_start
  #  $res{$Uid}{frag_start} = $frag_start;
    
    # getting fragment GC & length
   # ($res{$Uid}{frag_GC}, $res{$Uid}{frag_len}) =
    #  calc_GC( $frag_seq );
    
    # copying all other values from reads
   # map{ $res{$Uid}{$_} = $reads_r->{$genome}{$Uid}{$_} } 
   #   keys %{$reads_r->{$genome}{$Uid}};

  }

  # gathering
  MCE->gather( $genome, \%res );  
}

sub get_frag_GC{
# foreach read:
## determine random fragment size
## pull fragment from genome
### start-end around amplicon determine from a uniform distribution
## calculate fragment GC
## record GC & length
  my ($db, $reads_r, $size_dist, 
      $frag_min, $frag_max, $mean, $stdev,
     $primer_buffer) = @_;

  # sanity checks
  ## values provided
  map{ confess "ERROR: argument missing: $!\n" unless defined $_ } 
    ($size_dist, $frag_min, $frag_max, $mean, $stdev, $primer_buffer);

  # interating through fragments
  foreach my $genome (keys %$reads_r){
    my $genome_len = $db->length($genome);
    die "ERROR: cannot find '$genome' in genomes fasta\n"
      unless defined $genome_len;

    foreach my $Uid (keys %{$reads_r->{$genome}}){
      my $amp_start = $reads_r->{$genome}{$Uid}{amp_start};
      my $amp_end = $reads_r->{$genome}{$Uid}{amp_end};
      my $amp_len = $reads_r->{$genome}{$Uid}{amp_len};

      # fragment size
      my $frag_size;
      my $time_out = 0;
      if($size_dist eq 'uniform'){
	while(1){
	  $time_out++;
	  die "ERROR for genome $genome: could not get sequence in length range within $time_out attempts\n"
	    if $time_out > 9999;
	  $frag_size = random_uniform_integer(1, $frag_min, $frag_max); 
	  last if $frag_size >= $frag_min && $frag_size <= $frag_max;
	}
      }
      elsif($size_dist eq 'normal'){
	while(1){
	  $time_out++;
	  die "ERROR for genome $genome: could not get sequence in length range within $time_out attempts\n"
	    if $time_out > 9999;
	  $frag_size = int random_normal(1, $mean, $stdev);
	  last if $frag_size >= $frag_min && $frag_size <= $frag_max;
	}
      }
      elsif($size_dist eq 'exponential'){
	while(1){
	  $time_out++;
	  die "ERROR for genome $genome: could not get sequence in length range within $time_out attempts\n"
	    if $time_out > 9999;
	  $frag_size = int random_exponential(1, $mean);
	  last if $frag_size >= $frag_min && $frag_size <= $frag_max;
	}
      }
      else{ confess "ERROR: do not recognize size distribution\n"; }
            
      # determine fragment start-end based on amplicon start-end
      ## amplicon center postion
      my $amp_center = int($amp_start + abs($amp_len * 0.5));
      ## fragment_start = amp_center - (frag_size * x); x = random draw from unifrom distribution 0:1
      my $x = random_uniform();
      ### start_floor = (primer_buffer + 0.5*amp_len) / frag_size
      my $floor  = ($primer_buffer + 0.5 * $amp_len) / $frag_size; 
      $x = $floor if $x < $floor;
      my $ceiling = 1 - $floor;                               # ceiling = 1 - floor
      $x = $ceiling if $x > $ceiling;
      my $frag_start = int( $amp_center - ($frag_size * $x) );

      ### sanity check 
      carp "WARNING for genome '$genome': frag_start is too far from amplicon!\n\tfrag_start: $frag_start, amp_end: $amp_end, primer_buffer: $primer_buffer, frag_size: $frag_size, amp_center: $amp_center, x: $x\n"
	if $frag_start - ($amp_end + $primer_buffer) > $frag_size;


      # getting fragment
      ## fragment start-end can span ends of genome (assuming circular genome)
      ## need to account for this; making an %@ of 1-3 start-end positions
      my %pos;
      ### left wrap
      if($frag_start < 0){
	$pos{left_wrap} = [$genome_len + $frag_start, $genome_len];
	$pos{middle} = [1, $frag_size + $frag_start];
      }
      if(! exists $pos{middle}){
	$pos{middle} = [$frag_start, $frag_start + $frag_size -1 ];
      }
      if($frag_start + $frag_size -1 > $genome_len){
	$pos{right_wrap} = [1, $frag_start + $frag_size -1 - $genome_len];
	${$pos{middle}}[1] = [$genome_len];
      }

      ### getting sequence
      my $frag_seq = "";
      map{ $frag_seq .= $db->seq($genome, ${$pos{$_}}[0], ${$pos{$_}}[1]) } keys %pos; 

      ### sanity check 
      #my $actual_frag_len = length $frag_seq;
      #carp "WARNING for genome '$genome': frag_length -> '$actual_frag_len' is not within range! Start_frag_length: $frag_size; frag_start: $frag_start; frag_start2: $frag_start2; genome_length: ", $genome_seqs_r->{$genome}{len}, "\n"
      # 	if $actual_frag_len < $frag_min or $actual_frag_len > $frag_max;
      
      # storing amp_start & fragment_start
      $reads_r->{$genome}{$Uid}{frag_start} = $frag_start;

      # getting fragment GC & length
      ($reads_r->{$genome}{$Uid}{frag_GC},$reads_r->{$genome}{$Uid}{frag_len}) =
       	calc_GC( $frag_seq );

      # sanity check
      carp "WARNING: for genome '$genome': frag_length -> ", $reads_r->{$genome}{$Uid}{frag_len},
	"is not within length range!\tStart_frag_len: $frag_size, frag_start: $frag_start\n"
	  if $reads_r->{$genome}{$Uid}{frag_len} < $frag_min or 
	    $reads_r->{$genome}{$Uid}{frag_len} > $frag_max;
    }
  }

  #print Dumper $reads_r; exit;
}

=head write_read_info

Writing out all read info if wanted.
For debugging.

=cut

push @EXPORT_OK, 'write_read_info';

sub write_read_info{
  my ($reads_r, $db) = @_;

#  print Dumper $reads_r; exit;
  
  # output file
#  open OUT, ">$outfile" or die $!;

  # header
  print join("\t", qw/genome scaffold read_ID read_GC 
			 read_length read_start 
			 frag_GC frag_length frag_start/), "\n";

  # body
  foreach my $genome (keys %$reads_r){
    (my $desc = $db->header($genome)) =~ s/^[^ ]+ //;
    foreach my $Uid (keys %{$reads_r->{$genome}}){
      print join("\t",
		     #$reads_r->{$genome}{$Uid}{desc},
		     $desc,
		     $genome, 
		     $Uid, 
		     $reads_r->{$genome}{$Uid}{amp_GC},
		     $reads_r->{$genome}{$Uid}{amp_len},
		     $reads_r->{$genome}{$Uid}{amp_start},		     
		     $reads_r->{$genome}{$Uid}{frag_GC},
		     $reads_r->{$genome}{$Uid}{frag_len},
		     $reads_r->{$genome}{$Uid}{frag_start}
		    ), "\n";
    }
  }

  #close OUT or die $!;
  #print STDERR "Read info file written to: '$outfile'\n";
}


=head1 get_genome_GC_stats

Calculating GC stats per genome. 
Foreach read: 
   calculate GC diff beteween amplicon & fragment
   load into stat object
calculate stats
load stats into return hash

=cut

push @EXPORT_OK, 'get_genome_GC_stats';

sub get_genome_GC_stats{

  my $genome = $_;
  my ($self) = @_;
  my $reads_r = $self->user_args->{reads_r};

  my %stats;
  my $genome_id;
  my $genome_stat = new Statistics::PointEstimation;
  foreach my $Uid (keys %{$reads_r->{$genome}}){

    $genome_id = join("|", $reads_r->{$genome}{$Uid}{desc}, $genome)
      unless defined $genome_id;
    
    # loading deltaGC of frag & read(amplicon)
    my $deltaGC = abs($reads_r->{$genome}{$Uid}{amp_GC} - 
			$reads_r->{$genome}{$Uid}{frag_GC});
    $genome_stat->add_data($deltaGC);
  }

  # calculting genome stats
  $stats{mean} = $genome_stat->mean();
  $stats{variance} = $genome_stat->variance();
  $stats{df} = $genome_stat->df();
  $genome_stat->set_significance(95);
  $stats{upper_clm_95} = $genome_stat->upper_clm();
  $genome_stat->set_significance(99);
  $stats{upper_clm_99} = $genome_stat->upper_clm();

  MCE->gather($genome_id, \%stats);
}


=head1 get_total_GC_stats

Calculating GC stats for all genomes combined.
Same as get_genome_GC_stats, but reads from all
genomes are combining. 

=cut

push @EXPORT_OK, 'get_total_GC_stats';

sub get_total_GC_stats{
  my ($reads_r, $stats_r) = @_;

  my $total_stat = new Statistics::PointEstimation;
  foreach my $genome (keys %$reads_r){
    foreach my $Uid (keys %{$reads_r->{$genome}}){

      # loading deltaGC of frag & read(amplicon)
      my $deltaGC = abs($reads_r->{$genome}{$Uid}{amp_GC} - 
			$reads_r->{$genome}{$Uid}{frag_GC});
      $total_stat->add_data($deltaGC);
    }
  }
  
  # total 
  $stats_r->{TOTAL}{mean} = $total_stat->mean();
  $stats_r->{TOTAL}{variance} = $total_stat->variance();
  $stats_r->{TOTAL}{df} = $total_stat->df();
  $total_stat->set_significance(95);
  $stats_r->{TOTAL}{upper_clm_95} = $total_stat->upper_clm();
  $total_stat->set_significance(99);
  $stats_r->{TOTAL}{upper_clm_99} = $total_stat->upper_clm();
  
  #print Dumper $stats_r;  exit;
}

=head2 write_stats_summary

Writing a summary table to STDOUT

Using %stats produced by sub get_GC_stats

=cut

push @EXPORT_OK, 'write_stats_summary';

sub write_stats_summary{
  my ($stats_r) = @_;

  foreach my $genome ( sort{ 
    if($a eq 'TOTAL'){ return -1; }
    elsif($b eq 'TOTAL'){ return 1; }
    else{ return $a cmp $b }
    } keys %$stats_r ){
    
    print join("\t", $genome,
	       $stats_r->{$genome}{mean},
	       $stats_r->{$genome}{variance},
	       $stats_r->{$genome}{df},
	       $stats_r->{$genome}{upper_clm_95},
	       $stats_r->{$genome}{upper_clm_99}
	       ), "\n";
  }

}





=head1 AUTHOR

Nick Youngblut, C<< <ndy2 at cornell.edu> >>

=head1 BUGS

Please report any bugs or feature requests to C<bug-deltaMass at rt.cpan.org>, or through
the web interface at L<http://rt.cpan.org/NoAuth/ReportBug.html?Queue=deltaMass>.  I will be notified, and then you'll
automatically be notified of progress on your bug as I make changes.




=head1 SUPPORT

You can find documentation for this module with the perldoc command.

    perldoc deltaMass


You can also look for information at:

=over 4

=item * RT: CPAN's request tracker (report bugs here)

L<http://rt.cpan.org/NoAuth/Bugs.html?Dist=deltaMass>

=item * AnnoCPAN: Annotated CPAN documentation

L<http://annocpan.org/dist/deltaMass>

=item * CPAN Ratings

L<http://cpanratings.perl.org/d/deltaMass>

=item * Search CPAN

L<http://search.cpan.org/dist/deltaMass/>

=back


=head1 ACKNOWLEDGEMENTS


=head1 LICENSE AND COPYRIGHT

Copyright 2014 Nick Youngblut.

This program is free software; you can redistribute it and/or modify it
under the terms of either: the GNU General Public License as published
by the Free Software Foundation; or the Artistic License.

See http://dev.perl.org/licenses/ for more information.


=cut

1; # End of deltaMass
