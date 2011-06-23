#!/usr/bin/perl
# (c) 2010 Magnus Bjursell

# standard module usage
use strict;
use Getopt::Long;
use Pod::Usage;
use Cwd 'abs_path';
use FindBin;
use Data::Dumper;

# custom module usage
use lib $FindBin::Bin; #use lib '/home/magnus.bjursell/script/modules';
use mySharedFunctions qw(:basic myBasename fixPath isCanonicalChr chrTranslate returnGenomeSize);


my $dataHR = {};
my $infoHR = {};

#$infoHR->{'configHR'} = readConfigFile();
$infoHR->{'optHR'}    = {};


# Read input and setup global variables
GetOptions ($infoHR->{'optHR'}, 'ref=s', 'bed=s', 'bam=s', 'outdir=s', 'outPrefix=s', 'outCoverage=s', 'minCoverage=s', 'minConcQual=s', 'minMappingQual=s', 'help', 'verbose') or pod2usage(2);
pod2usage(1) if $infoHR->{'optHR'}{'h'} or $infoHR->{'optHR'}{'help'};
loadDefaults($infoHR);

pod2usage("Cannot find the reference file: $infoHR->{'optHR'}{'ref'}.") unless -f $infoHR->{'optHR'}{'ref'};
pod2usage("Cannot find the regions bed file: $infoHR->{'optHR'}{'bed'}.") unless -f $infoHR->{'optHR'}{'bed'};
pod2usage("Cannot find the input bam file: $infoHR->{'optHR'}{'bam'}.") unless -f $infoHR->{'optHR'}{'bam'};
pod2usage("Cannot find output directory: $infoHR->{'optHR'}{'outdir'}.") unless -d $infoHR->{'optHR'}{'outdir'};

$infoHR->{'optHR'}{'outPrefix'} = myBasename($infoHR->{'optHR'}{'bam'}) unless $infoHR->{'optHR'}{'outPrefix'};

#my $outBedFN = shift; die "\nPlease supply a filename for the output coverage filename\n" . $usage if $outBedFN =~ /\/$/;
checkSamTools();
getGenomeSize($infoHR);


{ my $idx = 1; $infoHR->{'chrSortOrder'} = { map { 'chr' . $_ => $idx++ } (1 .. 22, 'X', 'Y', 'M') }; }
#{ $infoHR->{'chrTranslate'} = { map { $_ => 'chr' . $_ } (1 .. 22, 'X', 'Y', 'M') }; $infoHR->{'chrTranslate'}{'MT'} = "chrM"; }
#print Dumper($infoHR);


# Read bed file
{
  local $| = 1;
  my ($noLines) = ( `wc -l $infoHR->{'optHR'}{'bed'}` =~ /(\d+)/ );
  print "\n";
  print "Reading bed file ($noLines lines; " . int($noLines/10000) . " dots; printing one dot per 10,000 lines)";
  my $bedFH = myOpen($infoHR->{'optHR'}{'bed'});
  while ( my $line = <$bedFH> ) {
    print "." if $. % 10000 == 0;
    chomp($line);
    my @cols = split(/\t/, $line);
    $cols[0] = chrTranslate($cols[0]) if defined(chrTranslate($cols[0]));
    next unless isCanonicalChr($cols[0]);
    for my $pos ( $cols[1] .. $cols[2] - 1 ) {
      next if $dataHR->{'bed'}{$cols[0]}{$pos};
      $dataHR->{'bed'}{$cols[0]}{$pos} = 1;
      $infoHR->{'no_bed_positions'}++;
      $infoHR->{'coverage'}{'0'}++;
    }
  }
  print "Done\n\n";
}


# Read data from pileup files
{
  local $| = 1;
  my $outBedFH = myOpenRW($infoHR->{'optHR'}{'outCoverage'}) if $infoHR->{'optHR'}{'outCoverage'};
  print "Reading pileup file; printing one dot per 10,000,000 lines)";

  my $noDone = 1;

  my $fullCmdLine = "samtools pileup -cA -f $infoHR->{'optHR'}{'ref'} $infoHR->{'optHR'}{'bam'} | awk '{print \$1\"\t\"\$2\"\t\"\$3\"\t\"\$5\"\t\"\$7\"\t\"\$8}'";
  my $pipe = myOpen("$fullCmdLine |");
  while ( my $line = <$pipe> ) {
    print "." if $noDone++ % 10000000 == 0;
    chomp($line);
    my @cols = split(/\t/, $line);
    next if $cols[2] =~ /\*/;
    $cols[0] = chrTranslate($cols[0]) if defined(chrTranslate($cols[0]));
    next unless isCanonicalChr($cols[0]);
    printf $outBedFH ("%s\t%d\t%d\t%s\t%s\t%d\n", $cols[0], $cols[1], $cols[1] + 1, join(":", @cols[0..5]), "+",$cols[5]) if $infoHR->{'optHR'}{'outCoverage'};
    $infoHR->{'no_pileup_positions'}++;
    if ( $dataHR->{'bed'}{$cols[0]}{$cols[1]} ) {
      $infoHR->{'coverage'}{'all'}{$cols[5]}++;
      $infoHR->{'coverage'}{'chr'}{$cols[0]}{$cols[5]}++;
      $infoHR->{'coverage'}{'bed'}{$cols[5]}++;

      $infoHR->{'total_coverage'}{'all'} += $cols[5];
      $infoHR->{'total_coverage'}{'bed'} += $cols[5];
      $infoHR->{'total_coverage'}{'chr'}{$cols[0]} += $cols[5];
      if ( $cols[5] >= $infoHR->{'optHR'}{'minCoverage'} and $cols[3] >= $infoHR->{'optHR'}{'minConcQual'} and $cols[4] >= $infoHR->{'optHR'}{'minMappingQual'} ) {
        $infoHR->{'highqual_coverage'}{'all'}++; 
        $infoHR->{'highqual_coverage'}{'bed'}++; 
        $infoHR->{'highqual_coverage'}{'chr'}{$cols[0]}++;
      }
    } else {
      $infoHR->{'coverage'}{'all'}{$cols[5]}++;
      $infoHR->{'coverage'}{'chr'}{$cols[0]}{$cols[5]}++;

      $infoHR->{'total_coverage'}{'all'} += $cols[5];
      $infoHR->{'total_coverage'}{'chr'}{$cols[0]} += $cols[5];
      if ( $cols[5] >= $infoHR->{'optHR'}{'minCoverage'} and $cols[3] >= $infoHR->{'optHR'}{'minConcQual'} and $cols[4] >= $infoHR->{'optHR'}{'minMappingQual'} ) {
        $infoHR->{'highqual_coverage'}{'all'}++;
        $infoHR->{'highqual_coverage'}{'chr'}{$cols[0]}++;
      }
    }
  }
  close($pipe);
  close($outBedFH) if $infoHR->{'optHR'}{'outCoverage'};
  print "Done\n\n";

}


print "Analyzing data\n\n";
# Calculate coverage histogram
{
  for ( my $x = max(keys %{$infoHR->{'coverage'}{'all'}}) - 1; $x >= 1; $x-- ) {
    $infoHR->{'coverage'}{'all'}{$x} += $infoHR->{'coverage'}{'all'}{$x + 1};
    foreach my $chr ( keys %{$infoHR->{'chrSortOrder'}} ) { $infoHR->{'coverage'}{'chr'}{$chr}{$x} += $infoHR->{'coverage'}{'chr'}{$chr}{$x + 1}; }
    $infoHR->{'coverage'}{'bed'}{$x} += $infoHR->{'coverage'}{'bed'}{$x + 1};
  }
}

print "Done\n\n";


# Print results
{

# Open outputfiles
  {
    my ($bamName) = ($infoHR->{'optHR'}{'bam'} =~ /([^\/]+)$/);
    foreach my $file ( "coverage_per_chromosome", "coverage_histogram_per_bed_and_genome", "coverage_histogram_per_chromosome", "summary", "EZsummary" ) {
      my $name = fixPath(($infoHR->{'optHR'}{'outdir'} ? "$infoHR->{'optHR'}{'outdir'}/" : "") . "$infoHR->{'optHR'}{'outPrefix'}.$file.txt");
      push(@{$infoHR->{'of'}}, {'name' => $name, 'fh' => myOpenRW($name)});
      print "Output $file file:\t$name\n";
    }
    print "\n";
  }

  {
# print summary and EZsummary files
    my $idx = 3;

    print {$infoHR->{'of'}[$idx]{'fh'}} "Input reference file:\t" . abs_path($infoHR->{'optHR'}{'ref'}) . "\n";
    print {$infoHR->{'of'}[$idx]{'fh'}} "Input regions bed file:\t" . abs_path($infoHR->{'optHR'}{'bed'}) . "\n";
    print {$infoHR->{'of'}[$idx]{'fh'}} "Input bam file:\t" . abs_path($infoHR->{'optHR'}{'bam'}) . "\n";
    print {$infoHR->{'of'}[$idx]{'fh'}} "Output directory:\t" . abs_path($infoHR->{'optHR'}{'outdir'}) . "\n";
    print {$infoHR->{'of'}[$idx]{'fh'}} "\n";

    printf {$infoHR->{'of'}[$idx]{'fh'}} ("Number of positions in bed file:\t%d\n", $infoHR->{'no_bed_positions'});
    printf {$infoHR->{'of'}[$idx]{'fh'}} ("Number of positions in pileup file:\t%d\n", $infoHR->{'no_pileup_positions'});
    print {$infoHR->{'of'}[$idx]{'fh'}} "\n";

    printf {$infoHR->{'of'}[$idx]{'fh'}} ("Total coverage:\t%d\n", $infoHR->{'total_coverage'}{'all'});
    printf {$infoHR->{'of'}[$idx]{'fh'}} ("Average coverage over genome (genome size: %d):\t%0.3f\tX\n", $infoHR->{'genomeSizeHR'}{'genome'}, $infoHR->{'total_coverage'}{'all'} / $infoHR->{'genomeSizeHR'}{'genome'});
    printf {$infoHR->{'of'}[$idx]{'fh'}} ("Average coverage over regions with coverage:\t%0.3f\tX\n", $infoHR->{'total_coverage'}{'all'} / $infoHR->{'no_pileup_positions'});
    print {$infoHR->{'of'}[$idx]{'fh'}} "\n";

    printf {$infoHR->{'of'}[$idx]{'fh'}} ("Total number of high qual positions (>= %d X covr, >= consensus_Q%d, >= RMS_mapping_Q%d;), fractions of genome and fractions of pileup pos:\t%d\t%0.3f\t%0.3f\n", @{$infoHR->{'optHR'}}{'minCoverage', 'minConcQual', 'minMappingQual'}, 
        $infoHR->{'highqual_coverage'}{'all'}, $infoHR->{'highqual_coverage'}{'all'} / $infoHR->{'genomeSizeHR'}{'genome'},$infoHR->{'highqual_coverage'}{'all'} / $infoHR->{'no_pileup_positions'});
    printf {$infoHR->{'of'}[$idx]{'fh'}} ("Number of high qual positions within bed regions (>= %d X covr, >= consensus_Q%d, >= RMS_mapping_Q%d;), fractions of region bed:\t%d\t%0.3f\n", @{$infoHR->{'optHR'}}{'minCoverage', 'minConcQual', 'minMappingQual'}, 
        $infoHR->{'highqual_coverage'}{'bed'}, $infoHR->{'highqual_coverage'}{'bed'} / $infoHR->{'no_bed_positions'});
    print {$infoHR->{'of'}[$idx]{'fh'}} "\n";

    printf {$infoHR->{'of'}[$idx]{'fh'}} ("Total coverage in bed region:\t%d\n", $infoHR->{'total_coverage'}{'bed'});
    printf {$infoHR->{'of'}[$idx]{'fh'}} ("Ratio coverage in bed region / total coverage:\t%0.4f\n", $infoHR->{'total_coverage'}{'bed'} / $infoHR->{'total_coverage'}{'all'});
    printf {$infoHR->{'of'}[$idx]{'fh'}} ("Ratio high quality bed positions / total high quality positions:\t%0.4f\n", $infoHR->{'highqual_coverage'}{'bed'} / $infoHR->{'highqual_coverage'}{'all'});



    printf {$infoHR->{'of'}[$idx+1]{'fh'}} ("NoBedFilePositions:\t%d\tpos\n", $infoHR->{'no_bed_positions'});
    printf {$infoHR->{'of'}[$idx+1]{'fh'}} ("NoBamFilePositions:\t%d\tpos\n", $infoHR->{'no_pileup_positions'});
    printf {$infoHR->{'of'}[$idx+1]{'fh'}} ("TotalCoverage:\t%d\tbases\n", $infoHR->{'total_coverage'}{'all'});

    printf {$infoHR->{'of'}[$idx+1]{'fh'}} ("GenomeSize:\t%d\tbases\n", $infoHR->{'genomeSizeHR'}{'genome'});
    printf {$infoHR->{'of'}[$idx+1]{'fh'}} ("AverageCoveragePosInGenome\t%0.3f\tX\n", $infoHR->{'total_coverage'}{'all'} / $infoHR->{'genomeSizeHR'}{'genome'});
    printf {$infoHR->{'of'}[$idx+1]{'fh'}} ("AverageCoveragePosInBam\t%0.3f\tX\n", $infoHR->{'total_coverage'}{'all'} / $infoHR->{'no_pileup_positions'});

    printf {$infoHR->{'of'}[$idx+1]{'fh'}} ("HighQualMinCoverage\t%d\tX\n", $infoHR->{'optHR'}{'minCoverage'});
    printf {$infoHR->{'of'}[$idx+1]{'fh'}} ("HighQualMinConsensusQual\t%d\tQual\n", $infoHR->{'optHR'}{'minConcQual'});
    printf {$infoHR->{'of'}[$idx+1]{'fh'}} ("HighQualMinMappingQual\t%d\tQual\n", $infoHR->{'optHR'}{'minMappingQual'});

    printf {$infoHR->{'of'}[$idx+1]{'fh'}} ("TotalHighCoveragePos\t%d\tpos\n", $infoHR->{'highqual_coverage'}{'all'});
    printf {$infoHR->{'of'}[$idx+1]{'fh'}} ("FractionHighCoveragePosInGenome\t%0.3f\tfrac\n", $infoHR->{'highqual_coverage'}{'all'} / $infoHR->{'genomeSizeHR'}{'genome'});
    printf {$infoHR->{'of'}[$idx+1]{'fh'}} ("FractionHighCoveragePosInBam\t%0.3f\tfrac\n", $infoHR->{'highqual_coverage'}{'all'} / $infoHR->{'no_pileup_positions'});

    printf {$infoHR->{'of'}[$idx+1]{'fh'}} ("TotalHighCoveragePosInBed\t%d\tpos\n", $infoHR->{'highqual_coverage'}{'bed'});
    printf {$infoHR->{'of'}[$idx+1]{'fh'}} ("FractionHighCoveragePosInBed\t%0.3f\tfrac\n", $infoHR->{'highqual_coverage'}{'bed'} / $infoHR->{'no_bed_positions'});
    printf {$infoHR->{'of'}[$idx+1]{'fh'}} ("TotalCoverageInBed:\t%d\tbases\n", $infoHR->{'total_coverage'}{'bed'});

    printf {$infoHR->{'of'}[$idx+1]{'fh'}} ("RatioCoverageInBedToAllCoverage:\t%d\tbases\n", $infoHR->{'total_coverage'}{'bed'} / $infoHR->{'total_coverage'}{'all'});
    printf {$infoHR->{'of'}[$idx+1]{'fh'}} ("RatioHighQualPosInBedToAllHighQualPos:\t%d\tbases\n", $infoHR->{'highqual_coverage'}{'bed'} / $infoHR->{'highqual_coverage'}{'all'});


    unless ( $infoHR->{'total_coverage'}{'bed'} ) {
      print {$infoHR->{'of'}[$idx]{'fh'}} "No coverage in bed regions...\n";
    } else {
# Calculate base80 penalty

      for my $x ( 0 .. max(keys %{$infoHR->{'coverage'}{'all'}}) ) {
        if ( $infoHR->{'coverage'}{'bed'}{$x+1} / $infoHR->{'no_bed_positions'} < 0.8 ) {
          $infoHR->{'base80_coverage'} = $x + ( ($infoHR->{'coverage'}{'bed'}{$x} - ( 0.8 * $infoHR->{'no_bed_positions'})) / ($infoHR->{'coverage'}{'bed'}{$x} - $infoHR->{'coverage'}{'bed'}{$x+1}) );
          last;
        }
      }

      printf {$infoHR->{'of'}[$idx]{'fh'}} ("Average coverage in bed regions:\t%0.3f\tX\n", ( $infoHR->{'total_coverage'}{'bed'} / $infoHR->{'no_bed_positions'} ) );
      printf {$infoHR->{'of'}[$idx]{'fh'}} ("Base 80 coverage:\t%0.3f\tX\n", $infoHR->{'base80_coverage'} );
      printf {$infoHR->{'of'}[$idx]{'fh'}} ("Base 80 penalty:\t%0.3f\n", ( ( $infoHR->{'total_coverage'}{'bed'} / $infoHR->{'no_bed_positions'} ) / $infoHR->{'base80_coverage'} ));
      print {$infoHR->{'of'}[$idx]{'fh'}} "\n";

      printf {$infoHR->{'of'}[$idx+1]{'fh'}} ("AverageCoveragePosInBed\t%0.3f\tfrac\n", $infoHR->{'total_coverage'}{'bed'} / $infoHR->{'no_bed_positions'});
      printf {$infoHR->{'of'}[$idx+1]{'fh'}} ("Base80Coverage\t%0.3f\tX\n", $infoHR->{'base80_coverage'});
      printf {$infoHR->{'of'}[$idx+1]{'fh'}} ("Base80Penalty\t%0.3f\tpen\n", ( $infoHR->{'total_coverage'}{'bed'} / $infoHR->{'no_bed_positions'} ) / $infoHR->{'base80_coverage'});

    }
  }

  {
# Export coverage per chromosome
    my $idx = 0;
    print "Exporting coverage per chromosome to file $infoHR->{'of'}[$idx]{'name'}...\n";
    print {$infoHR->{'of'}[$idx]{'fh'}} "Chromosome\tSize\tHigh qual pos\tTotal coverage\tAvg coverage\n";
    foreach my $chr ( sort { $infoHR->{'chrSortOrder'}{$a} <=> $infoHR->{'chrSortOrder'}{$b} } keys %{$infoHR->{'chrSortOrder'}} ) {
      printf {$infoHR->{'of'}[$idx]{'fh'}} ("%s\t%d\t%d\t%d\t%0.2f\tX\n", $chr, $infoHR->{'genomeSizeHR'}{'chr'}{$chr}, $infoHR->{'highqual_coverage'}{'chr'}{$chr}, $infoHR->{'total_coverage'}{'chr'}{$chr}, $infoHR->{'total_coverage'}{'chr'}{$chr} / $infoHR->{'genomeSizeHR'}{'chr'}{$chr});
    }
    print "Done\n\n";
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
      printf {$infoHR->{'of'}[$idx]{'fh'}} ("%0.3f\n", $infoHR->{'coverage'}{'all'}{$x} / $infoHR->{'genomeSizeHR'}{'genome'});
    }
    print "Done\n\n";
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
      foreach my $chr ( sort { $infoHR->{'chrSortOrder'}{$a} <=> $infoHR->{'chrSortOrder'}{$b} } keys %{$infoHR->{'chrSortOrder'}} ) { printf {$infoHR->{'of'}[$idx]{'fh'}} ("\t%0.3f", $infoHR->{'coverage'}{'chr'}{$chr}{$x} / $infoHR->{'genomeSizeHR'}{'chr'}{$chr} ); }
      print {$infoHR->{'of'}[$idx]{'fh'}} "\n";
    }
    print "Done\n\n";
  }

  foreach my $outFileHR ( @{$infoHR->{'of'}} ) { close($outFileHR->{'fh'}); }

}

exit;


# # #  S U B S  # # #

sub loadDefaults {
  my $infoHR = shift;

  $infoHR->{'optHR'}{'minCoverage'}    = 8  unless $infoHR->{'optHR'}{'minCoverage'};
  $infoHR->{'optHR'}{'minConcQual'}    = 30 unless $infoHR->{'optHR'}{'minConcQual'};
  $infoHR->{'optHR'}{'minMappingQual'} = 30 unless $infoHR->{'optHR'}{'minMappingQual'};
}


sub getGenomeSize {
  my $infoHR = shift;

  print STDERR "Calculating genome size\n";
  $infoHR->{'genomeSizeHR'} = returnGenomeSize($infoHR->{'optHR'}{'ref'}); #2861343702; #2858674662
  print STDERR "Done\n\n";
}


sub checkSamTools {
  my $ret = `samtools pileup 2>&1`;

  my $score = 0;
  foreach ( split(/\n/, $ret) ) {
    $score++ if /\-A\s+use the MAQ model for SNP calling/;
    $score++ if /\-c\s+output the SOAPsnp consensus sequence/;
  }
  die "Seems to be a problem with samtools. I need the pileup function and I am looking for the consensus and MAQ model options\n\n" unless $score == 2;
}



#__END__


=head1 NAME

calculate_coverage_statistics.pl - Analysis of sequencing coverage across specified regions

=head1 SYNOPSIS

calculate_coverage_statistics.pl [options]

  Options:
   --help            Brief help message (not yet implemented)
   --verbose         Write some additional output (not yet implemented)
   --outdir          Output directory [required]
   --outPrefix       Output file prefix, bam file name will be used if not supplied
   --outCoverage     Output bed file for coverage (will be large; contains
                     coverage information for each position with coverage [optional]
   --ref             Reference sequence fasta file [required]
   --bed             Bed file describing the regions [required]
   --bam             Bam file containing the alignments [required]
   --shm             Use shared memory (/dev/shm)
   --minCoverage     Minimum coverage (default: 8)
   --minConcQual     Minimum concensus quality (default: 30)
   --minMappingQual  Minimum RMS mapping quality (default: 30)


Reference (fasta), regions (bed) and alignment (bam) files are required.

Note that chromosome names should match (e.g. chr1 vs 1, chrM vs MT etc.)

This program requires Samtools v 0.1.8 (other versions may work).

=cut

#=head1 OPTIONS

#=over 8

#=item B<-help>

#Print a brief help message and exits.

#=back

#=head1 DESCRIPTION




