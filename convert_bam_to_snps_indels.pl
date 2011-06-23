#!/usr/bin/perl
#/sw/comp/perl/5.12.3_kalkyl/bin/perl
# (c) 2011 Magnus Bjursell

#BEGIN { push @INC,"/sw/comp/perl/5.12.3_kalkyl/lib/"; }


# standard module usage
use strict;
use Getopt::Long;
use Pod::Usage;
use Cwd 'abs_path';
use FindBin;
use Data::Dumper;

# custom module usage
use lib $FindBin::Bin; #use lib '/home/magnus.bjursell/script/modules';
use mySharedFunctions qw(:basic returnFaFHandIndex returnSequenceFromFasta myBasename);

my $dataHR = {};
my $infoHR = {};
#my $geneHR = {};

$infoHR->{'optHR'} = {};

# Read input and setup global variables
GetOptions ($infoHR->{'optHR'}, 'ref=s', 'bed=s', 'bam=s', 'outdir=s', 'outPrefix=s', 'outCoverage=s', 'minCoverage=s', 'minVarQual=s', 'minConcQual=s', 'minMappingQual=s',
                                'minIndelNoRds=s', 'maxDelFrc=s', 'minIndelRcnt=s', 'minIndelFrac=s', 'maxOtherFrac=s', 'help', 'verbose') or pod2usage(2);
pod2usage(1) if $infoHR->{'optHR'}{'h'} or $infoHR->{'optHR'}{'help'};
loadDefaults($infoHR);

pod2usage("Cannot find the reference file: $infoHR->{'optHR'}{'ref'}.") unless -f $infoHR->{'optHR'}{'ref'};
pod2usage("Cannot find the input bam file: $infoHR->{'optHR'}{'bam'}.") unless -f $infoHR->{'optHR'}{'bam'};
pod2usage("Cannot find output directory: $infoHR->{'optHR'}{'outdir'}.") unless -d $infoHR->{'optHR'}{'outdir'};

my $faiHR = returnFaFHandIndex($infoHR->{'optHR'}{'ref'});
if ( myBasename($infoHR->{'optHR'}{'bam'}) =~ /(.+)\.bam$/ ) { $infoHR->{'optHR'}{'bamBaseName'} = $1; } else { $infoHR->{'optHR'}{'bamBaseName'} = myBasename($infoHR->{'optHR'}{'bam'}); }
$infoHR->{'optHR'}{'outPrefix'} = $infoHR->{'optHR'}{'bamBaseName'} unless $infoHR->{'optHR'}{'outPrefix'};

$infoHR->{'optHR'}{'outdir_results'} = fixPath($infoHR->{'optHR'}{'outdir'} . "/results/");
my $ret = `mkdir -p $infoHR->{'optHR'}{'outdir_results'}` unless -d $infoHR->{'optHR'}{'outdir_results'};

{
  my $bedFH = myOpenRW(fixPath(($infoHR->{'optHR'}{'outdir'} ? "$infoHR->{'optHR'}{'outdir'}/" : "") . "$infoHR->{'optHR'}{'outPrefix'}.variants.vbed"));
#  my $annovarFH = myOpen(fixPath($infoHR->{'optHR'}{'outdir'} ? "$infoHR->{'optHR'}{'outdir'}/" : "" . "$infoHR->{'optHR'}{'outPrefix'}.variants.annovar"));
  my $sumFH = myOpenRW(fixPath(($infoHR->{'optHR'}{'outdir_results'} ? "$infoHR->{'optHR'}{'outdir_results'}/" : "") . "$infoHR->{'optHR'}{'outPrefix'}.summary.txt"));
  my $aFqFH = myOpenRW(fixPath(($infoHR->{'optHR'}{'outdir_results'} ? "$infoHR->{'optHR'}{'outdir_results'}/" : "") . "$infoHR->{'optHR'}{'outPrefix'}.alleleReadFrequecies.txt"));

  print_parameters($sumFH);
  local $| = 1;
  my $awkCode = "{if (\$5 >= $infoHR->{'optHR'}{'minConcQual'} && \$6 >= $infoHR->{'optHR'}{'minVarQual'} && \$7 >= $infoHR->{'optHR'}{'minMappingQual'} && \$8 >= $infoHR->{'optHR'}{'minCoverage'}) {print};}";
  my $fullCmdLine = "samtools pileup -cAv -f $infoHR->{'optHR'}{'ref'} $infoHR->{'optHR'}{'bam'} | awk '$awkCode'";
  my $pipe = myOpen("$fullCmdLine |");
  while ( my $singleEntry = <$pipe> ) {
    my @cols = split(/\t/, $singleEntry);
    print $bedFH parse_pileup_line(\@cols, $infoHR);
    print "Processing $cols[0]...\n" unless exists($dataHR->{$cols[0]});
    $dataHR->{$cols[0]} = 1;
  }
  close($pipe);
  print "Done.\n\n";
  print $sumFH "Number of SNVs:\t$infoHR->{'stats'}{'noSNVs'}\n";
  print $sumFH "Number of homozygous SNVs:\t$infoHR->{'stats'}{'noHomSNVs'}\n";
  print $sumFH "Number of heterozygous SNVs:\t$infoHR->{'stats'}{'noHetSNVs'}\n";
  print $sumFH "Number of transitions:\t$infoHR->{'stats'}{'noTsTv'}{'transition'}\n";
  print $sumFH "Number of transversions:\t$infoHR->{'stats'}{'noTsTv'}{'transversion'}\n";
  print $sumFH "\n";

  foreach my $toBase ( sort keys %{$infoHR->{'stats'}{'toBaseCount'}} ) { print $sumFH "\t$toBase"; } print $sumFH "\n";
  foreach my $fromBase ( sort keys %{$infoHR->{'stats'}{'fromBaseCount'}} ) {
    print $sumFH "$fromBase";
    foreach my $toBase ( sort keys %{$infoHR->{'stats'}{'toBaseCount'}} ) {
      print $sumFH "\t$infoHR->{'stats'}{'noSpecChngs'}{$fromBase}{$toBase}";
    }
    print $sumFH "\n";
  }
  print $sumFH "\n\n";

  print $sumFH "Num of alleles\tCount\n";
  foreach my $numAll ( sort {$a <=> $b} keys %{$infoHR->{'stats'}{'SNValleleHisto'}} ) {
    print $sumFH "$numAll\t$infoHR->{'stats'}{'SNValleleHisto'}{$numAll}\t";
    print $sumFH join("\t", map(sprintf("%0.2f", $_ / $infoHR->{'stats'}{'SNValleleHisto'}{$numAll}), @{$infoHR->{'stats'}{'SNValleleHistoCounts'}{$numAll}})) . "\n";
  }
  print $sumFH "\n\n";

  print $aFqFH "Pct reads supporting alternative allele\tSNVs (counts)\n";
  my $maxKey = max(keys %{$infoHR->{'stats'}{'SNValleleFreqHisto'}});
  foreach my $allFrq ( 0 .. $maxKey ) { print $aFqFH "$allFrq\t$infoHR->{'stats'}{'SNValleleFreqHisto'}{$allFrq}\n"; }
  print $aFqFH "\n\n";

  print $sumFH "Number of Indels:\t$infoHR->{'stats'}{'noIndels'}\n";
  print $sumFH "Number of insertions:\t$infoHR->{'stats'}{'noIndelType'}{'insertion'}\n";

  print $sumFH "Length\tCount\n";
  foreach my $len ( sort {$a <=> $b} keys %{$infoHR->{'stats'}{'noIndelType'}{'insertionHisto'}} ) {
    print $sumFH "$len\t$infoHR->{'stats'}{'noIndelType'}{'insertionHisto'}{$len}\n";
  }
  print $sumFH "\n";

  print $sumFH "Number of deletions:\t$infoHR->{'stats'}{'noIndelType'}{'deletion'}\n";
  print $sumFH "Length\tCount\n";
  foreach my $len ( sort {$a <=> $b} keys %{$infoHR->{'stats'}{'noIndelType'}{'deletionHisto'}} ) {
    print $sumFH "$len\t$infoHR->{'stats'}{'noIndelType'}{'deletionHisto'}{$len}\n";
  }
  print $sumFH "\n\n";

  close($bedFH);
  close($sumFH);
  close($aFqFH);
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

  $infoHR->{'stats'}{'SNValleleHisto'}{scalar(keys %{$cntHR})}++;
  $infoHR->{'stats'}{'SNValleleFreqHisto'}{int( 100 * ( ($colsAR->[7] - $cntHR->{uc($colsAR->[2])}) / ($colsAR->[7] + 0.000000001) ) )}++;

  {
    my $idx = 0;
    foreach my $base ( sort {$cntHR->{$b} <=> $cntHR->{$a}} keys %{$cntHR} ) {
      $infoHR->{'stats'}{'SNValleleHistoCounts'}{scalar(keys %{$cntHR})}[$idx++] += round( 100 * ( $cntHR->{$base} / ($colsAR->[7] + 0.000000001) ), 2 ); #$cntHR->{$base};
    }
  }

#print "No alleles: " . scalar(keys %{$cntHR}) . "; cntStr: $cntStr; " . join(":", keys %{$cntHR}) . "; AlStr: $colsAR->[8];\n" if scalar(keys %{$cntHR}) >= 4;
#print "No same as ref (" . uc($colsAR->[2]) . "): " . $cntHR->{uc($colsAR->[2])} . "; Total covr: $colsAR->[7]; No alt alleles: " . int( 100 * ( ($colsAR->[7] - $cntHR->{uc($colsAR->[2])}) / ($colsAR->[7] + 0.000000001) ) ) . "\n\n" if uc($colsAR->[3]) !~ /^[ACGT]$/;

#  warn "Reference bases differ between input ref and bam file for " . join("\t", @{$colsAR}) . "...\n" unless uc(returnSequenceFromFasta($faiHR, $colsAR->[0], ($colsAR->[1] - 1), 1)) eq uc($colsAR->[2]);

  if ( $noDels / $noSeqs > $infoHR->{'optHR'}{'maxDelFrc'} ) {
    return "";
  } else {
    my $fNt = uc($colsAR->[2]);
    my $tNt = uc($colsAR->[3]);
    
    $infoHR->{'stats'}{'noSNVs'}++;
    if ( $tNt =~ /^[ACGT]$/ ) { $infoHR->{'stats'}{'noHomSNVs'}++; } else { $infoHR->{'stats'}{'noHetSNVs'}++; }
    $infoHR->{'stats'}{'noSpecChngs'}{$fNt}{$tNt}++;
    $infoHR->{'stats'}{'noTsTv'}{'transition'}++   if ($fNt eq "A" and $tNt eq "G") or ($fNt eq "G" and $tNt eq "A") or ($fNt eq "C" and $tNt eq "T") or ($fNt eq "T" and $tNt eq "C");
    $infoHR->{'stats'}{'noTsTv'}{'transversion'}++ if ($fNt eq "A" and $tNt eq "C") or ($fNt eq "C" and $tNt eq "A") or ($fNt eq "G" and $tNt eq "T") or ($fNt eq "T" and $tNt eq "G") or
                                                      ($fNt eq "A" and $tNt eq "T") or ($fNt eq "T" and $tNt eq "A") or ($fNt eq "C" and $tNt eq "G") or ($fNt eq "G" and $tNt eq "C");
    $infoHR->{'stats'}{'fromBaseCount'}{$fNt}++; 
    $infoHR->{'stats'}{'toBaseCount'}{$tNt}++;
    return sprintf ("%s\t%d\t%d\t%s=>%s(%d;%s)\t%d\t+\n", $colsAR->[0], ($colsAR->[1] - 1), $colsAR->[1], $fNt, $tNt, $colsAR->[7], $cntStr, $colsAR->[5] );
  }
}


sub parse_pileup_indel_line {
  my $colsAR = shift;
  my $infoHR = shift;
  die "\nIncorrect pileup indel format, expected 15 columns, got this:\n" . join("\t", @{$colsAR}) . "\n\n" unless scalar(@{$colsAR}) == 15;

  return "" if $colsAR->[12] / $colsAR->[7] > $infoHR->{'optHR'}{'maxOtherFrac'};
  return "" if $colsAR->[7] < $infoHR->{'optHR'}{'minIndelNoRds'};

  my $alleleHR = {$colsAR->[8] => $colsAR->[10], $colsAR->[9] => $colsAR->[11], 'O' => $colsAR->[12], };
  delete($alleleHR->{'O'}) unless $alleleHR->{'O'} > 0;
  my $pass = 0; my $length = 0;
  my $cntStr  = "";
  foreach my $keys ( sort { $alleleHR->{$b} <=> $alleleHR->{$a} } keys %{$alleleHR} ) {
    $cntStr .= sprintf("%s:%d,", uc($keys), $alleleHR->{$keys});
    next if $keys eq "*" or $keys eq "O";
    if ( $keys =~ /([\+\-])(\w+)/ ) {
      $length = max($length, ($1 eq "-" ? length($2) : 0));
      $pass = 1 if $alleleHR->{$keys} >= $infoHR->{'optHR'}{'minIndelRcnt'} and $alleleHR->{$keys} / $colsAR->[7] >= $infoHR->{'optHR'}{'minIndelFrac'};
    }
  }
  $cntStr =~ s/,$//;

  if ( $pass ) {
    my $start = $length == 0 ? $colsAR->[1] - 1 : $colsAR->[1];
    my $end   = $length == 0 ? $colsAR->[1] - 1 : $colsAR->[1] + $length;
    $infoHR->{'stats'}{'noIndels'}++;
    if ( $colsAR->[3] =~ /\+(\w+)/ ) {
      $infoHR->{'stats'}{'noIndelType'}{'insertionHisto'}{length($1)}++;
      $infoHR->{'stats'}{'noIndelType'}{'insertion'}++;
    } elsif ( $colsAR->[3] =~ /\-(\w+)/ ) {
      $infoHR->{'stats'}{'noIndelType'}{'deletionHisto'}{length($1)}++;
      $infoHR->{'stats'}{'noIndelType'}{'deletion'}++;
    }
    return sprintf ("%s\t%d\t%d\t%s=>%s(%d;%s)\t%d\t+\n", $colsAR->[0], $start, $end, uc(returnSequenceFromFasta($faiHR, $colsAR->[0], $start, 1)), uc($colsAR->[3]), $colsAR->[7], $cntStr, $colsAR->[5] );
  } else {
    return "";
  }
}


sub print_parameters {
  my $fh = shift; $fh = *STDOUT unless $fh;

  print $fh "\n";
  print $fh "Prameters:\n";
  print $fh "   - Min mapping quality: $infoHR->{'optHR'}{'minMappingQual'}\n";
  print $fh "   - Min consensus quality: $infoHR->{'optHR'}{'minConcQual'}\n";
  print $fh "   - Min variant quality: $infoHR->{'optHR'}{'minVarQual'}\n";
  print $fh "   - Min number of reads: $infoHR->{'optHR'}{'minCoverage'}\n";
  print $fh "   - Max deletion fraction: $infoHR->{'optHR'}{'maxDelFrc'}\n";
  print $fh "\n";

  print $fh "   - Min number of reads for indels: $infoHR->{'optHR'}{'minIndelNoRds'}\n";
  print $fh "   - Min number of reads supporting indel: $infoHR->{'optHR'}{'minIndelRcnt'}\n";
  print $fh "   - Min indel fraction: $infoHR->{'optHR'}{'minIndelFrac'}\n";
  print $fh "   - Max 3rd+ alleles fraction: $infoHR->{'optHR'}{'maxOtherFrac'}\n";
  print $fh "\n";

}

sub loadDefaults {
  my $infoHR = shift;

  $infoHR->{'optHR'}{'minCoverage'}    = 8     unless $infoHR->{'optHR'}{'minCoverage'};
  $infoHR->{'optHR'}{'minConcQual'}    = 30    unless $infoHR->{'optHR'}{'minConcQual'};
  $infoHR->{'optHR'}{'minVarQual'}     = 30    unless $infoHR->{'optHR'}{'minVarQual'};
  $infoHR->{'optHR'}{'minMappingQual'} = 30    unless $infoHR->{'optHR'}{'minMappingQual'};
  $infoHR->{'optHR'}{'maxDelFrc'}      = 0.75  unless $infoHR->{'optHR'}{'maxDelFrc'};

  $infoHR->{'optHR'}{'minIndelNoRds'}  = 12    unless $infoHR->{'optHR'}{'minIndelNoRds'};
  $infoHR->{'optHR'}{'minIndelRcnt'}   = 3     unless $infoHR->{'optHR'}{'minIndelRcnt'};
  $infoHR->{'optHR'}{'minIndelFrac'}   = 0.25  unless $infoHR->{'optHR'}{'minIndelFrac'};
  $infoHR->{'optHR'}{'maxOtherFrac'}   = 0.125 unless $infoHR->{'optHR'}{'maxOtherFrac'};
}




#__END__


=head1 NAME

calculate_coverage_statistics.pl - Analysis of sequencing coverage across specified regions

=head1 SYNOPSIS

convert_bam_to_snvs_indels.pl [options]

  Options:
   --help            Brief help message (not yet implemented)
   --verbose         Write some additional output (not yet implemented)
   --outdir          Output directory [required]
   --outPrefix       Output file prefix, bam file name will be used if not supplied
   --ref             Reference sequence fasta file [required]
   --bam             Bam file containing the alignments [required]
   --minCoverage     Minimum coverage (default: 8)
   --minVarQual      Minimum variant quality (default: 30)
   --minConcQual     Minimum concensus quality (default: 30)
   --minMappingQual  Minimum RMS mapping quality (default: 30)
   --minIndelNoRds   Minimum coverage from indel calling (default: 12)
   --minIndelRcnt    Minimum coverage supporting indel (default: 3)
   --minIndelFrac    Minimum fraction supporting indel (defaullt: 0.25)
   --maxOtherFrac    Maximum fraction supporting a third+ allele (default 0.125)


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






