#!/usr/bin/perl
# (c) 2009 Magnus Bjursell

# standard module usage
use strict;
#use IO::Handle;
#use feature ":5.10";
#use File::stat;


# custom module usage
#use lib '/home/magnus.bjursell/script/modules';
#use myFileIO qw(:basic);
#use myMathStat qw(max min round);
#use mySeqAnalysis qw(isCanonicalChr);


# Read input and setup global variables
my $progName; if ( $0 =~ /([^\/]+)$/ ) { $progName = $1; }
my $usage = "\nUsage: $progName [reference fasta seq] [bed region file] [BAM file] [output dir] [covr bed output - optional]\n\n";

my $refFN    = shift; die "\nCannot find the reference file: $refFN\n" . $usage unless -f $refFN;
my $bedFN    = shift; die "\nCannot find the regions bed file: $bedFN\n" . $usage unless -f $bedFN;
my $bamFN    = shift; die "\nCannot find the input bam file: $bamFN\n" . $usage unless -f $bamFN;
my $outDN    = shift;
my $outBedFN = shift; die "\nPlease supply a filename for the output coverage filename\n" . $usage if $outBedFN =~ /\/$/;
checkSamTools();

my $dataHR = {};
my $infoHR = {'covCO' => 8, 'conQualCO' => 20};

{ my $idx = 1; $infoHR->{'chrSortOrder'} = { map { 'chr' . $_ => $idx++ } (1 .. 22, 'X', 'Y', 'M') }; }
{ $infoHR->{'chrTranslate'} = { map { $_ => 'chr' . $_ } (1 .. 22, 'X', 'Y', 'M') }; $infoHR->{'chrTranslate'}{'MT'} = "chrM"; }


# Read bed file
{
  local $| = 1;
  my ($noLines) = ( `wc -l $bedFN` =~ /(\d+)/ );
  print STDERR "\n";
  print STDERR "Reading bed file ($noLines lines; " . int($noLines/10000) . " dots; printing one dot per 10,000 lines)";
  my $bedFH = myOpen($bedFN);
  while ( my $line = <$bedFH> ) {
    print STDERR "." if $. % 10000 == 0;
    chomp($line);
    my @cols = split(/\t/, $line);
    $cols[0] = $infoHR->{'chrTranslate'}{$cols[0]} if defined($infoHR->{'chrTranslate'}{$cols[0]});
    next unless isCanonicalChr($infoHR, $cols[0]);
    for my $pos ( $cols[1] .. $cols[2] - 1 ) {
      next if $dataHR->{'bed'}{$cols[0]}{$pos};
      $dataHR->{'bed'}{$cols[0]}{$pos} = 1;
      $infoHR->{'no_bed_positions'}++;
      $infoHR->{'coverage'}{'0'}++;
    }
  }
  print STDERR "Done\n\n";
}


# Read data from pileup files
{
  local $| = 1;
  my $outBedFH = myOpenRW($outBedFN) if $outBedFN;
  print STDERR "Reading pileup file; printing one dot per 10,000,000 lines)";

  my $noDone = 1;

  my $fullCmdLine = "samtools pileup -cA -f $refFN $bamFN | awk '{print \$1\"\t\"\$2\"\t\"\$3\"\t\"\$5\"\t\"\$8}'";
  my $pipe = myOpen("$fullCmdLine |");
  while ( my $line = <$pipe> ) {
    print STDERR "." if $noDone++ % 10000000 == 0;
    chomp($line);
    my @cols = split(/\t/, $line);
    next if $cols[2] =~ /\*/;
    $cols[0] = $infoHR->{'chrTranslate'}{$cols[0]} if defined($infoHR->{'chrTranslate'}{$cols[0]});
    next unless isCanonicalChr($infoHR, $cols[0]);
    printf $outBedFH ("%s\t%d\t%d\t%s\t%s\t%d\n", $cols[0], $cols[1], $cols[1] + 1, join(":", @cols[0..4]), "+",$cols[4]) if $outBedFN;
    $infoHR->{'no_pileup_positions'}++;
    if ( $dataHR->{'bed'}{$cols[0]}{$cols[1]} ) {
      $infoHR->{'coverage'}{'all'}{$cols[4]}++;
      $infoHR->{'coverage'}{'chr'}{$cols[0]}{$cols[4]}++;
      $infoHR->{'coverage'}{'bed'}{$cols[4]}++;

      $infoHR->{'total_coverage'}{'all'} += $cols[4];
      $infoHR->{'total_coverage'}{'bed'} += $cols[4];
      $infoHR->{'total_coverage'}{'chr'}{$cols[0]} += $cols[4];
      if ( $cols[4] >= $infoHR->{'covCO'} and $cols[3] >= $infoHR->{'conQualCO'} ) {
        $infoHR->{'highqual_coverage'}{'all'}++; 
        $infoHR->{'highqual_coverage'}{'bed'}++; 
        $infoHR->{'highqual_coverage'}{'chr'}{$cols[0]}++;
      }
    } else {
      $infoHR->{'coverage'}{'all'}{$cols[4]}++;
      $infoHR->{'coverage'}{'chr'}{$cols[0]}{$cols[4]}++;

      $infoHR->{'total_coverage'}{'all'} += $cols[4];
      $infoHR->{'total_coverage'}{'chr'}{$cols[0]} += $cols[4];
      if ( $cols[4] >= $infoHR->{'covCO'} and $cols[3] >= $infoHR->{'conQualCO'} ) {
        $infoHR->{'highqual_coverage'}{'all'}++;
        $infoHR->{'highqual_coverage'}{'chr'}{$cols[0]}++;
      }
    }
  }
  close($pipe);
  close($outBedFH) if $outBedFN;
  print STDERR "Done\n\n";

}


print STDERR "Analyzing data\n\n";
# Calculate coverage histogram
{
  for ( my $x = max(keys %{$infoHR->{'coverage'}{'all'}}) - 1; $x >= 1; $x-- ) {
    $infoHR->{'coverage'}{'all'}{$x} += $infoHR->{'coverage'}{'all'}{$x + 1};
    foreach my $chr ( keys %{$infoHR->{'chrSortOrder'}} ) { $infoHR->{'coverage'}{'chr'}{$chr}{$x} += $infoHR->{'coverage'}{'chr'}{$chr}{$x + 1}; }
    $infoHR->{'coverage'}{'bed'}{$x} += $infoHR->{'coverage'}{'bed'}{$x + 1};
  }
}

print STDERR "Done\n\n";


# Print results
{
  print STDERR "Calculating genome size\n";
  my $genomeSizeHR = returnGenomeSize($refFN); #2861343702; #2858674662
  print STDERR "Done\n\n";

  print "Input bed file:\t$bedFN\n";
  print "Input BAM file:\t$bamFN\n";

  {
    my ($bamName) = ($bamFN =~ /([^\/]+)$/);
    foreach my $file ( "coverage_per_chromosome", "coverage_histogram_per_bed_and_genome", "coverage_histogram_per_chromosome" ) {
      my $name = ($outDN ? "$outDN/" : "") . "$bamName.$file.txt"; $name =~ s/\/+/\//g;
      push(@{$infoHR->{'of'}}, {'name' => $name, 'fh' => myOpenRW($name)});
      print "Output $file file:\t$name\n";
    }
    print "\n";
  }

  printf("Number of positions in bed file:\t%d\n", $infoHR->{'no_bed_positions'});
  printf("Number of positions in pileup file:\t%d\n", $infoHR->{'no_pileup_positions'});
  print "\n";

  printf("Total coverage:\t%d\n", $infoHR->{'total_coverage'}{'all'});
  printf("Average coverage over genome (genome size: %d):\t%0.3f\tX\n", $genomeSizeHR->{'genome'}, $infoHR->{'total_coverage'}{'all'} / $genomeSizeHR->{'genome'});
  printf("Average coverage over regions with coverage:\t%0.3f\tX\n", $infoHR->{'total_coverage'}{'all'} / $infoHR->{'no_pileup_positions'});
  print "\n";

  printf("Total number of high qual positions (>=8X covr, >=consensus_Q20; frac of genome and pileup pos):\t%d\t%0.3f\t%0.3f\n", $infoHR->{'highqual_coverage'}{'all'},
      $infoHR->{'highqual_coverage'}{'all'} / $genomeSizeHR->{'genome'},$infoHR->{'highqual_coverage'}{'all'} / $infoHR->{'no_pileup_positions'});
  printf("Number of high qual positions within bed regions (>=8X covr, >=consensus_Q20; frac bed pos):\t%d\t%0.3f\n", $infoHR->{'highqual_coverage'}{'bed'}, $infoHR->{'highqual_coverage'}{'bed'} / $infoHR->{'no_bed_positions'});
  print "\n";

  printf("Total coverage in bed region:\t%d\n", $infoHR->{'total_coverage'}{'bed'});


  unless ( $infoHR->{'total_coverage'}{'bed'} ) {
    print "No coverage in bed regions...\n";

  } else {
# Calculate base80 penalty

    for my $x ( 0 .. max(keys %{$infoHR->{'coverage'}{'all'}}) ) {
      if ( $infoHR->{'coverage'}{'bed'}{$x+1} / $infoHR->{'no_bed_positions'} < 0.8 ) {
        $infoHR->{'base80_coverage'} = $x + ( ($infoHR->{'coverage'}{'bed'}{$x} - ( 0.8 * $infoHR->{'no_bed_positions'})) / ($infoHR->{'coverage'}{'bed'}{$x} - $infoHR->{'coverage'}{'bed'}{$x+1}) );
        last;
      }
    }

    printf("Average coverage in bed regions:\t%0.3f\tX\n", ( $infoHR->{'total_coverage'}{'bed'} / $infoHR->{'no_bed_positions'} ) );
    printf("Base 80 coverage:\t%0.3f\tX\n", $infoHR->{'base80_coverage'} );
    printf("Base 80 penalty:\t%0.3f\n", ( ( $infoHR->{'total_coverage'}{'bed'} / $infoHR->{'no_bed_positions'} ) / $infoHR->{'base80_coverage'} ));
    print "\n";
  }


  {
# Export coverage per chromosome
    my $idx = 0;
    print "Exporting coverage per chromosome to file $infoHR->{'of'}[$idx]{'name'}...\n";
    print {$infoHR->{'of'}[$idx]{'fh'}} "Chromosome\tSize\tHigh qual pos\tTotal coverage\tAvg coverage\n";
    foreach my $chr ( sort { $infoHR->{'chrSortOrder'}{$a} <=> $infoHR->{'chrSortOrder'}{$b} } keys %{$infoHR->{'chrSortOrder'}} ) {
      printf {$infoHR->{'of'}[$idx]{'fh'}} ("%s\t%d\t%d\t%d\t%0.2f\tX\n", $chr, $genomeSizeHR->{'chr'}{$chr}, $infoHR->{'highqual_coverage'}{'chr'}{$chr}, $infoHR->{'total_coverage'}{'chr'}{$chr}, $infoHR->{'total_coverage'}{'chr'}{$chr} / $genomeSizeHR->{'chr'}{$chr});
    }
  }


  {
# Export coverage histogram per bed, region with any coverage and genome
    my $idx = 1;
    print "Exporting coverage histogram for positions in bed file, positions with coverage and genome to file $infoHR->{'of'}[$idx]{'name'}...\n";
    print {$infoHR->{'of'}[$idx]{'fh'}} "Coverage\tNo positions (bed)\tFraction of positions (bed)\tNo positions (covr)\tFraction of positions (covr)\tFraction of positions (genome)\n";
    printf {$infoHR->{'of'}[$idx]{'fh'}} ("%d\t%d\t%0.3f\t%d\t%0.3f\t%0.3f\n", 0, $infoHR->{'no_bed_positions'}, 1, $infoHR->{'no_pileup_positions'}, 1, 1);
    for my $x ( 1 .. min(max(keys %{$infoHR->{'coverage'}{'all'}}), 500) ) {
      printf {$infoHR->{'of'}[$idx]{'fh'}} ("%d\t%d\t%0.3f\t", $x, $infoHR->{'coverage'}{'bed'}{$x}, $infoHR->{'coverage'}{'bed'}{$x} / $infoHR->{'no_bed_positions'});
      printf {$infoHR->{'of'}[$idx]{'fh'}} ("%d\t%0.3f\t", $infoHR->{'coverage'}{'all'}{$x}, $infoHR->{'coverage'}{'all'}{$x} / $infoHR->{'no_pileup_positions'});
      printf {$infoHR->{'of'}[$idx]{'fh'}} ("%0.3f\n", $infoHR->{'coverage'}{'all'}{$x} / $genomeSizeHR->{'genome'});
    }
  }


  {
# Export coverage histogram per chromosome
    my $idx = 2;
    print "Exporting coverage histogram by chromosome to file $infoHR->{'of'}[$idx]{'name'}...\n";
#    print {$infoHR->{'of'}[$idx]{'fh'}} "Coverage\tFraction of positions\n";
    print {$infoHR->{'of'}[$idx]{'fh'}} "cvr"; foreach my $chr ( sort { $infoHR->{'chrSortOrder'}{$a} <=> $infoHR->{'chrSortOrder'}{$b} } keys %{$infoHR->{'chrSortOrder'}} ) { printf {$infoHR->{'of'}[$idx]{'fh'}} ("\t%s", $chr); } print {$infoHR->{'of'}[$idx]{'fh'}} "\n";
    print {$infoHR->{'of'}[$idx]{'fh'}} "0"; foreach my $chr ( sort { $infoHR->{'chrSortOrder'}{$a} <=> $infoHR->{'chrSortOrder'}{$b} } keys %{$infoHR->{'chrSortOrder'}} ) { printf {$infoHR->{'of'}[$idx]{'fh'}} ("\t%0.3f", 1); } print {$infoHR->{'of'}[$idx]{'fh'}} "\n";

    for my $x ( 1 .. min(max(keys %{$infoHR->{'coverage'}{'all'}}), 500) ) {
      print {$infoHR->{'of'}[$idx]{'fh'}} "$x";
      foreach my $chr ( sort { $infoHR->{'chrSortOrder'}{$a} <=> $infoHR->{'chrSortOrder'}{$b} } keys %{$infoHR->{'chrSortOrder'}} ) { printf {$infoHR->{'of'}[$idx]{'fh'}} ("\t%0.3f", $infoHR->{'coverage'}{'chr'}{$chr}{$x} / $genomeSizeHR->{'chr'}{$chr} ); }
      print {$infoHR->{'of'}[$idx]{'fh'}} "\n";
    }
  }

  foreach my $outFileHR ( @{$infoHR->{'of'}} ) { close($outFileHR->{'fh'}); }

}

exit;


# # #  S U B S  # # #

sub checkSamTools {
  my $ret = `samtools pileup 2>&1`;

  my $score = 0;
  foreach ( split(/\n/, $ret) ) {
    $score++ if /\-A\s+use the MAQ model for SNP calling/;
    $score++ if /\-c\s+output the SOAPsnp consensus sequence/;
  }
  die "Seems to be a problem with samtools. I need the pileup function and I am looking for the consensus and MAQ model options\n\n" unless $score == 2;
}

sub returnGenomeSize {
  my $refFN = shift;
  my $genomeSizeHR = {};
  my $refFH = myOpen($refFN);
  my $currentChr = "";
  while ( my $line = <$refFH> ) {
    chomp($line);
    if ( $line =~ /^>(\S+)/ ) {
      $currentChr = $1;
      $currentChr = $infoHR->{'chrTranslate'}{$currentChr} if defined($infoHR->{'chrTranslate'}{$currentChr});
      die "Error: Unrecognized chromosome name: $currentChr\n" unless $currentChr =~ /^chr[1-9XYM][0-9]?$/;
    } else {
      my $noBases = $line =~ tr/ACGTacgt//;
      $genomeSizeHR->{'genome'} += $noBases;
      $genomeSizeHR->{'chr'}{$currentChr} += $noBases;
    }
  }
  close($refFH);
  return $genomeSizeHR;
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

sub isCanonicalChr {
  my $infoHR = shift;
  my $chr    = shift;
  return exists($infoHR->{'chrSortOrder'}{$chr}) ? $infoHR->{'chrSortOrder'}{$chr} : undef;
}

#sub estimate_filesize {
#  my $inFN = shift;
#  my $noLines = 0;
#  my $lineLen = 0; my $sampleSize = 0; my $fileSize = stat($inFN)->size;
#  my $inFH = myOpen($inFN); while ( my $line = <$inFH> ) { $lineLen += length($line); last if $sampleSize++ > 100000; } close($inFH);
#  $noLines = int($fileSize / ( $lineLen / $sampleSize ));
#  return $noLines;
#}


