#!/usr/bin/perl
# (c) 2010 Magnus Bjursell

# standard module usage
#use lib '/bubo/home/h3/magnusb/perl/ActivePerl-5.10.1.1008-x86_64-linux-glibc-2.3.5-294165/perl/lib/';
use strict;
#use IO::Handle;
use feature ":5.10";
#use File::stat;


# custom module usage
#use lib '/home/magnus.bjursell/script/modules';
#use myFileIO qw(:basic);
#use myMathStat qw(max min round);
#use mySeqAnalysis qw(returnFaFHandIndex returnSequenceFromFasta);


# Read input and setup global variables
my $progName; if ( $0 =~ /([^\/]+)$/ ) { $progName = $1; }
my $usage = "\nUsage: $progName [ref fasta] [bam file]\n\n";

my $refFN   = shift; die $usage unless -e $refFN;
my $bamFN   = shift; die $usage unless -e $bamFN;

my $dataHR = {};
my $infoHR = {};
my $geneHR = {};
#my $faiHR = returnFaFHandIndex($refFN);

$infoHR->{'minConQ'}   = 20;
$infoHR->{'minVarQ'}   = 20;
$infoHR->{'minMapQ'}   = 20;
$infoHR->{'minNoRds'}  = 8;
$infoHR->{'maxDelFrc'} = 0.75;

$infoHR->{'minIndelNoRds'} = 12;
$infoHR->{'minIndelRcnt'}  = 3;
$infoHR->{'minIndelFrac'}  = 0.25;
$infoHR->{'maxOtherFrac'}  = 0.125;

{
  print_parameters();
  local $| = 1;
  my $awkCode = "{if (\$5 >= $infoHR->{'minConQ'} && \$6 >= $infoHR->{'minVarQ'} && \$7 >= $infoHR->{'minMapQ'} && \$8 >= $infoHR->{'minNoRds'}) {print};}";
  my $fullCmdLine = "samtools pileup -cAv -f $refFN $bamFN | awk '$awkCode'";
  my $pipe = myOpen("$fullCmdLine |");
  while ( my $singleEntry = <$pipe> ) {
    my @cols = split(/\t/, $singleEntry);
    print parse_pileup_line(\@cols, $infoHR);
    print STDERR "Processing $cols[0]...\n" unless exists($dataHR->{$cols[0]});
    $dataHR->{$cols[0]} = 1;
  }
  close($pipe);
  print STDERR "Done.\n\n";
  print STDERR "No SNVs: $infoHR->{'noSNVs'}\n";
  print STDERR "No Indels: $infoHR->{'noIndels'}\n\n";
}


exit;



# # #  S U B S  # # #

sub parse_pileup_line {
  my $colsAR = shift;
  my $infoHR = shift;
  die "\nIncorrect pileup format, expected >= 10 columns, got this:\n" . join("\t", @{$colsAR}) . "\n\n" unless scalar(@{$colsAR}) >= 10;
  if ( $colsAR->[2] eq "*" ) {
    return parse_pileup_indel_line($colsAR, $infoHR);
  } else {
    return parse_pileup_SNP_line($colsAR, $infoHR);
  }
}


sub parse_pileup_SNP_line {
  my $colsAR = shift;
  my $infoHR = shift;
  die "\nIncorrect pileup variant format, expected 10 columns, got this:\n" . join("\t", @{$colsAR}) . "\n\n" unless scalar(@{$colsAR}) == 10;

  my $noSeqs   = 0;
  my $noDels   = 0;
  my $noIndels = 0;
  my $cntHR    = {};
  while ( $colsAR->[8] =~ /(([\+\-])([0-9]+)((??{"[ACGTNacgtn]{$3}"}))|((\^)(.))?+([\*\.\,ACGTNacgtn])(\$)?+)/g ) {
    my $matchAR = [$1, $2, $3, $4, $5, $6, $7, $8, $9];
    $noIndels++, next if $matchAR->[2];
    my $base = "";
#    given($matchAR->[7]) {
#      when (/[\,\.]/) {$base = uc($colsAR->[2]);}
#      when (/\*/) {$base = "D";}
#      when (/[ACGTNacgtn]/) {$base = uc($_);}
#    }
    if ( $matchAR->[7] =~ /[\,\.]/ ) {
      $base = uc($colsAR->[2]);
    } elsif ( $matchAR->[7] =~ /\*/ ) {
      $base = "D";
    } elsif ( $matchAR->[7] =~ /[ACGTNacgtn]/ ) {
      $base = uc($matchAR->[7]);
    }
    $cntHR->{$base}++;
    $noSeqs++;
    $noDels++ if $matchAR->[7] =~ /\*/;
  }
  warn "Error parsing, $colsAR->[7] vs $noSeqs ($noSeqs bp, $noIndels indel starts (+/-[0-9][ACGTNacgtn]+) ) for line: " . join("\t", @{$colsAR}) . "\n" unless $colsAR->[7] == $noSeqs;

  my $cntStr = "";
  foreach my $base ( sort keys %{$cntHR} ) { $cntStr .= sprintf("%s:%d,", uc($base), $cntHR->{$base}); }
  $cntStr =~ s/,$//;

#  warn "Reference bases differ between input ref and bam file for " . join("\t", @{$colsAR}) . "...\n" unless uc(returnSequenceFromFasta($faiHR, $colsAR->[0], ($colsAR->[1] - 1), 1)) eq uc($colsAR->[2]);

  if ( $noDels / $noSeqs > $infoHR->{'maxDelFrac'} ) {
    return "";
  } else {
    $infoHR->{'noSNVs'}++;
    return sprintf ("%s\t%d\t%d\t%s=>%s(%d;%s)\t%d\t+\n", $colsAR->[0], ($colsAR->[1] - 1), $colsAR->[1], uc($colsAR->[2]), uc($colsAR->[3]), $colsAR->[7], $cntStr, $colsAR->[5] );
  }
}


sub parse_pileup_indel_line {
  my $colsAR = shift;
  my $infoHR = shift;
  die "\nIncorrect pileup indel format, expected 15 columns, got this:\n" . join("\t", @{$colsAR}) . "\n\n" unless scalar(@{$colsAR}) == 15;

  return "" if $colsAR->[12] / $colsAR->[7] > $infoHR->{'maxOtherFrac'};
  return "" if $colsAR->[7] < $infoHR->{'minIndelNoRds'};

  my $alleleHR = {$colsAR->[8] => $colsAR->[10], $colsAR->[9] => $colsAR->[11], 'O' => $colsAR->[12], };
  delete($alleleHR->{'O'}) unless $alleleHR->{'O'} > 0;
  my $pass = 0; my $length = 0;
  my $cntStr  = "";
  foreach my $keys ( sort { $alleleHR->{$b} <=> $alleleHR->{$a} } keys %{$alleleHR} ) {
    $cntStr .= sprintf("%s:%d,", uc($keys), $alleleHR->{$keys});
    next if $keys eq "*" or $keys eq "O";
    if ( $keys =~ /([\+\-])(\w+)/ ) {
      $length = max($length, ($1 eq "-" ? length($2) : 0));
      $pass = 1 if $alleleHR->{$keys} >= $infoHR->{'minIndelRcnt'} and $alleleHR->{$keys} / $colsAR->[7] >= $infoHR->{'minIndelFrac'};
    }
  }
  $cntStr =~ s/,$//;

  if ( $pass ) {
    my $start = $length == 0 ? $colsAR->[1] - 1 : $colsAR->[1];
    my $end   = $length == 0 ? $colsAR->[1] - 1 : $colsAR->[1] + $length;
    $infoHR->{'noIndels'}++;
#    return sprintf ("%s\t%d\t%d\t%s=>%s(%d;%s)\t%d\t+\n", $colsAR->[0], $start, $end, uc(returnSequenceFromFasta($faiHR, $colsAR->[0], $start, 1)), uc($colsAR->[3]), $colsAR->[7], $cntStr, $colsAR->[5] );
  } else {
    return "";
  }
}


sub print_parameters {
  print STDERR "\n";
  print STDERR "Prameters:\n";
  print STDERR "   - Min mapping quality: $infoHR->{'minMapQ'}\n";
  print STDERR "   - Min consensus quality: $infoHR->{'minConQ'}\n";
  print STDERR "   - Min variant quality: $infoHR->{'minVarQ'}\n";
  print STDERR "   - Min number of reads: $infoHR->{'minNoRds'}\n";
  print STDERR "   - Max deletion fraction: $infoHR->{'maxDelFrc'}\n";
  print STDERR "\n";

  print STDERR "   - Min number of reads for indels: $infoHR->{'minIndelNoRds'}\n";
  print STDERR "   - Min indel reads: $infoHR->{'minIndelRcnt'}\n";
  print STDERR "   - Min indel fraction: $infoHR->{'minIndelFrac'}\n";
  print STDERR "   - Max 3rd+ alleles fraction: $infoHR->{'maxOtherFrac'}\n";
  print STDERR "\n";

}


sub myOpen {
  my $fn = shift;
  if ( $fn eq "-" or $fn eq "STDIN" ) {
    return *STDIN;
  } else {
    open(my $fh, $fn) or die "\nCannot open file $fn for reading ($!)\n\n";
    return $fh;
  }
}


sub myOpenRW {
  my $fn = shift;
  open(my $fh, ">$fn") or die "\nCannot open file $fn for writing ($!)\n\n";
  return $fh;
}

sub min {
  my @num = @_;
  my $min = $num[0];
  foreach (@num) {$min = $_ if $_ < $min}
  return $min;
}

sub max {
  my @num = @_;
  my $max = $num[0];
  foreach (@num) {$max = $_ if $_ > $max}
  return $max;
}



