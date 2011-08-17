#!/usr/bin/perl
# (c) 2011 Magnus Bjursell

# standard module usage
use strict;
use Getopt::Long;
use Pod::Usage;
use Cwd 'abs_path';
use FindBin;
use Data::Dumper;
use IO::Handle;
use feature ":5.10";
use File::stat;
use Term::ReadKey;

# custom module usage
use lib $FindBin::Bin;
use mySharedFunctions qw(:basic myBasename fixPath isCanonicalChr chrTranslate returnGenomeSize returnOverlap read_bed_variant_file returnFaFHandIndex);
#use lib '/home/magnus.bjursell/script/modules';
#use myFileIO qw(:basic);
#use myMathStat qw(max min round);
#use mySeqAnalysis qw(translate IUPAC2baseA reverseComplement returnOverlap readKgXrefHR returnFaFHandIndex returnSequenceFromFasta complementA returnChrSortIndex getFileHeader read_bed_variant_file readGtpInfo readFilterInfo);

my $dataHR = {};
my $infoHR = {};

$infoHR->{'configHR'} = readConfigFile();
$infoHR->{'optHR'}    = {};


# Read input and setup global variables
GetOptions ($infoHR->{'optHR'}, 'varBed=s', 'outdir=s', 'outdb=s', 'help', 'verbose') or pod2usage(2);
pod2usage(1) if $infoHR->{'optHR'}{'h'} or $infoHR->{'optHR'}{'help'};
#loadDefaults($infoHR);

pod2usage("Cannot find the variant bed file: $infoHR->{'optHR'}{'varBed'}.") unless -f $infoHR->{'optHR'}{'varBed'};
pod2usage("Cannot find output directory: $infoHR->{'optHR'}{'outdir'}.") unless -d $infoHR->{'optHR'}{'outdir'};
#pod2usage("Cannot find output database: $infoHR->{'optHR'}{'outdb'}.") unless -d $infoHR->{'optHR'}{'outdb'}; # Should work for sqlite3 db check

$dataHR->{'paramHR'}       = {'ssLen' => 2, 'upsLen' => 1000, 'dnsLen' => 1000,};
$dataHR->{'faiHR'}         = returnFaFHandIndex($infoHR->{'configHR'}{'PATH'}{'REFERENCE_GENOME'});
$dataHR->{'kgXrefHR'}      = readKgXrefHR($infoHR->{'configHR'}{'PATH'}{'UCSC_KGXREF'});
$dataHR->{'geneFilterHR'}  = readFilterInfo($infoHR->{'configHR'}{'PATH'}{'UCSC_KNOWNGENE_CANONICAL'}, "id");


# Setup output stuff
setupOutput_parseSampleName($infoHR, $infoHR->{'optHR'}{'outdir'}, $infoHR->{'optHR'}{'varBed'});


# Load variants from file
loadVariationData($dataHR, $infoHR->{'optHR'}{'varBed'});


{ # Comparing to SNP file
  local $| = 1;
  print STDERR "Comparing to $infoHR->{'configHR'}{'PATH'}{'dbSNP_131'} known variations file\n";
  my $kVarFH = myOpen($infoHR->{'configHR'}{'PATH'}{'dbSNP_131'});
  while ( my $line = <$kVarFH> ) {
    print STDERR "." if $. % 1000000 == 0;
    my $kVarHR = parseKVar($line, $infoHR->{'configHR'}{'PATH'}{'dbSNP_131'});
#    next unless $kVarHR->{'weight'} == 1;
    next if $kVarHR->{'obs'} eq "lengthTooLong";
    $infoHR->{'no_valid_var_in_file'}++;
    next unless exists($dataHR->{'varChrPos'}{$kVarHR->{'chr'}});

    my $overlapVarHR = {};
    for my $pos ( $kVarHR->{'sta'} .. $kVarHR->{'end'} - 1 ) {
      next unless exists($dataHR->{'varChrPos'}{$kVarHR->{'chr'}}{$pos});
      foreach my $varID ( keys %{$dataHR->{'varChrPos'}{$kVarHR->{'chr'}}{$pos}} ) {
        $overlapVarHR->{$varID} = $dataHR->{'varChrPos'}{$kVarHR->{'chr'}}{$pos}{$varID} unless exists($overlapVarHR->{$varID});
      }
    }
    foreach my $varHR ( values %{$overlapVarHR} ) {
      detectKnownVarOverlap($dataHR, $varHR, $kVarHR);
    }
  }
  close($kVarFH);
  print STDERR "\n   - $infoHR->{'no_valid_var_in_file'} known variants in file\nDone\n\n";
}


{ # Comparing to gene file
  local $| = 1;
  print STDERR "Comparing to $infoHR->{'configHR'}{'PATH'}{'UCSC_KNOWNGENE_FULL'} gene file\n";
  my $okStatCodes = {'cmpl' => 1, 'none' => 1,};
  my $transFH = myOpen($infoHR->{'configHR'}{'PATH'}{'UCSC_KNOWNGENE_FULL'});
  while ( my $line = <$transFH> ) {
    print STDERR "." if $. % 10000 == 0;
    my $transHR = parseGene($line, $infoHR->{'configHR'}{'PATH'}{'UCSC_KNOWNGENE_FULL'});
    next if defined($dataHR->{'geneFilterHR'}) and not defined($dataHR->{'geneFilterHR'}{$transHR->{'id'}});
    next if ( defined($transHR->{'cStaSt'}) and $okStatCodes->{$transHR->{'cStaSt'}} != 1 ) or ( defined($transHR->{'cEndSt'}) and $okStatCodes->{$transHR->{'cEndSt'}} != 1 );
    next unless returnCDSlen($transHR) % 3 == 0 or $transHR->{'cSta'} == $transHR->{'cEnd'};
    $infoHR->{'no_valid_trans_in_file'}++;


    next unless exists($dataHR->{'varChr'}{$transHR->{'chr'}});

    foreach my $varHR ( sort { $a->{'sta'} <=> $b->{'sta'} } @{$dataHR->{'varChr'}{$transHR->{'chr'}}} ) {
      last if $varHR->{'sta'} > (($transHR->{'end'} + $dataHR->{'paramHR'}{'dnsLen'}) - 1);
      next if $varHR->{'end'} - 1 < ($transHR->{'sta'} - $dataHR->{'paramHR'}{'upsLen'});

      detectTransOverlap($dataHR, $varHR, $transHR);
    }
  }
  close($transFH);
  print STDERR "\n   - $infoHR->{'no_valid_trans_in_file'} known transcripts in file (after filterring)\nDone\n\n";
}


print "Input summary:\n";
printf("   - %d input variations in $infoHR->{'optHR'}{'varBed'} file\n",  scalar(keys %{$dataHR->{'varHR'}}));
printf("   - %d known variations in file\n",  $infoHR->{'no_valid_var_in_file'});
printf("   - %d transcripts in file (after filterring)\n",  $infoHR->{'no_valid_trans_in_file'});
print "\n";


# Print out data

{ # Write variant data
  print STDERR "Summarizing and printing variant data...\n";
  $dataHR->{'pphInfoHR'}{'pphFile'}{'name'} = $infoHR->{'outDirs'}{'pph'} . "input.pph";
  $dataHR->{'pphInfoHR'}{'pphFile'}{'fh'}   = myOpenRW($dataHR->{'pphInfoHR'}{'pphFile'}{'name'});
  my $outFH_HR = {};
  my $allKeysHR = {'all' => 0, 'known' => 998, 'novel' => 999, 'lastID' => 2,};
  foreach my $varID ( sort { $a <=> $b } keys %{$dataHR->{'varHR'}} ) {
    my $varHR = $dataHR->{'varHR'}{$varID};
    my $snpAllelesHR = { map { $_ => 1 } ($varHR->{'info'}{'ref'}, IUPAC2baseA($varHR->{'info'}{'var'})) } if $varHR->{'info'}{'varType'} eq "snp";
    if ( ref($varHR->{'kVar'}) =~ /ARRAY/ and scalar(@{$varHR->{'kVar'}}) > 0 ) { $infoHR->{'no_known_vars'}++; } else { $infoHR->{'no_novel_vars'}++; }
    $infoHR->{'no_varTypesHR'}{$varHR->{'info'}{'varType'}}++;
    my $annotComb = join("_", sort keys %{$varHR->{'annot'}}); $annotComb = "no_annotations" unless length($annotComb) > 0;
    $infoHR->{'varAnnotComb'}{$annotComb}++;
    foreach my $annot ( keys %{$varHR->{'annot'}} ) {
      writeQuickOutputFileVars($outFH_HR, $infoHR, "all", $annot, $varID, "vars");
      $infoHR->{'var'}{'annot'}{$annot}{'all'}++;
      if ( ref($varHR->{'kVar'}) =~ /ARRAY/ and scalar(@{$varHR->{'kVar'}}) > 0 ) {
        foreach my $kVarHR ( @{$varHR->{'kVar'}} ) {
          next unless length($kVarHR->{'varFile'}) > 0;
          $infoHR->{'var'}{'annot'}{$annot}{$kVarHR->{'varFile'}}++;
          if ( $allKeysHR->{$kVarHR->{'varFile'}} < 1 ) { $allKeysHR->{$kVarHR->{'varFile'}} = $allKeysHR->{'lastID'}++; }
        }
        $infoHR->{'var'}{'annot'}{$annot}{'known'}++;
        writeQuickOutputFileVars($outFH_HR, $infoHR, "known", $annot, $varID, "vars");
      } else {
        $infoHR->{'var'}{'annot'}{$annot}{'novel'}++;
        writeQuickOutputFileVars($outFH_HR, $infoHR, "novel", $annot, $varID, "vars");
        printf {$dataHR->{'pphInfoHR'}{'pphFile'}{'fh'}} ("%s:%d %s #%s;%s;%s\n", $varHR->{'chr'}, $varHR->{'sta'} + 1, join("/", keys %{$snpAllelesHR}), $varHR->{'id'}, $annot, "novel")
            if $varHR->{'info'}{'varType'} eq "snp" and $annot eq "candidateGenes";
      }
    }
  }
  foreach my $outFN ( keys %{$outFH_HR} ) { close($outFH_HR->{$outFN}); }
  close($dataHR->{'pphInfoHR'}{'pphFile'}{'fh'});

  writeAnnotVarsCombinationsToFile($infoHR);
  printNoVarsTypes($infoHR);

  print "Variants\n";
  foreach my $status ( sort { $allKeysHR->{$a} <=> $allKeysHR->{$b} } keys %{$allKeysHR} ) { next if $status eq "lastID"; print "\t$status"; }
  print "\n";
  foreach my $annot ( sort keys %{$infoHR->{'var'}{'annot'}} ) {
    print $annot;
    foreach my $status ( sort { $allKeysHR->{$a} <=> $allKeysHR->{$b} } keys %{$allKeysHR} ) {
      next if $status eq "lastID";
      printf("\t%d", $infoHR->{'var'}{'annot'}{$annot}{$status});
    }
    print "\n";
  }
  print "\n\n";
  print STDERR "Done\n\n";
}


# Run polyphen-2 for selected data
if ( $dataHR->{'analysisDatabase'} =~ /^ucscKnownGene/ ) {
  runPolyPhen2_parseRet($dataHR, $infoHR);
} else {
  warn "Polyphen-2 analysis is currently only available when using the \"ucscKnownGene\" database for annotation, skipping...\n\n";
}


my $allKeysHR = {'all' => 0, 'all_recessive' => 1, 'all_homozygous' => 2,
                 'known' => 10, 'known_recessive' => 11, 'known_homozygous' => 12,
                 'novel' => 20, 'novel_recessive' => 21, 'novel_homozygous' => 22,
                 'novel_polyphen-2' => 30, 'novel_polyphen-2_recessive' => 31, 'novel_polyphen-2_homozygous' => 32,};

{ # Write transcript data and group into genes
  print STDERR "Summarizing and printing transcript data...\n";
  my $outSkippedIDsFH = myOpenRW($infoHR->{'outDirs'}{'base'} . "/no_gene_ID_found.txt");
  foreach my $transID ( sort { $a <=> $b } keys %{$dataHR->{'transHR'}} ) {
    my $transHR = $dataHR->{'transHR'}{$transID};
    $transHR->{'transcriptDescriptionHR'}{$dataHR->{'kgXrefHR'}{$transID}{'description'}} = 1;
    my $ns_ss_cInd = {};
    my $geneID = $dataHR->{'gtpHR'}{$transID};
    $dataHR->{'geneDescriptionByGeneIDHR'}{$geneID}{$dataHR->{'kgXrefHR'}{$transID}{'description'}} = 1;
    print $outSkippedIDsFH "Cannot find gene id for $transID, skipping\n" and $infoHR->{'skipped_no_gene_IDs'}++ and next unless length($geneID) > 0;
    foreach my $annot ( keys %{$transHR->{'annot'}} ) {
      foreach my $varID ( sort { $a <=> $b } keys %{$transHR->{'annot'}{$annot}} ) {
        my $varHR = $transHR->{'annot'}{$annot}{$varID};
        foreach my $allele ( sort keys %{$varHR->{'info'}{'ns'}{$transID}{'alt'}} ) {
          my $mutStr = join("", @{$varHR->{'info'}{'ns'}{$transID}{'ref'}}{'refAA', 'relStaAA_1'}, $varHR->{'info'}{'ns'}{$transID}{'alt'}{$allele}{'snpAA'});
          $transHR->{'info'}{'ns'}{'mutStrHR'}{$mutStr} = {'pph_prediction' => $varHR->{'info'}{'ns'}{$transID}{'alt'}{$allele}{'pph_prediction'}, 'homHet' => $varHR->{'info'}{'homHet'}, 'pos' => $varHR->{'info'}{'ns'}{$transID}{'ref'}{'relStaAA_1'}, 'varType' => $varHR->{'info'}{'varType'}};
        }
      }
      foreach my $transGrp ( keys %{$transHR->{'transGrp'}} ) {
        if ( $transHR->{'transGrp'}{$transGrp}{$annot}{'any'} > 0 ) {
          $infoHR->{'trans'}{'annot'}{$annot}{$transGrp}++;
          push(@{$dataHR->{'genes'}{$transGrp}{$annot}{$geneID}}, $transHR);
        }
        if ( $transHR->{'transGrp'}{$transGrp}{$annot}{'sum'} >= 1 ) {
          $infoHR->{'trans'}{'annot'}{$annot}{$transGrp . '_recessive'}++;
          push(@{$dataHR->{'genes'}{$transGrp . '_recessive'}{$annot}{$geneID}}, $transHR);
        }
        if ( $transHR->{'transGrp'}{$transGrp}{$annot}{'hom'} >= 1 ) {
          $infoHR->{'trans'}{'annot'}{$annot}{$transGrp . '_homozygous'}++;
          push(@{$dataHR->{'genes'}{$transGrp . '_homozygous'}{$annot}{$geneID}}, $transHR);
        }
      }
    }
  }
  close($outSkippedIDsFH);
  print "Transcripts\n";
  foreach my $transGrp ( sort { $allKeysHR->{$a} <=> $allKeysHR->{$b} } keys %{$allKeysHR} ) { print "\t$transGrp"; }
  print "\n";
  foreach my $annot ( sort keys %{$infoHR->{'trans'}{'annot'}} ) {
    print $annot;
    foreach my $transGrp ( sort { $allKeysHR->{$a} <=> $allKeysHR->{$b} } keys %{$allKeysHR} ) {
      printf("\t%d",$infoHR->{'trans'}{'annot'}{$annot}{$transGrp});
    }
    print "\n";
  }
  print "\n\n";
  print STDERR "Done\n\n";

}


{ # Write gene and full data outputs
  print STDERR "Printing gene data...\n";
  print "Genes\n";
  foreach my $geneGrp ( sort { $allKeysHR->{$a} <=> $allKeysHR->{$b} } keys %{$allKeysHR} ) { print "\t$geneGrp"; }
  print "\n";
  foreach my $annot ( sort keys %{$dataHR->{'genes'}{'all'}} ) {
    print $annot;
    foreach my $geneGrp ( sort { $allKeysHR->{$a} <=> $allKeysHR->{$b} } keys %{$allKeysHR} ) { 
      printf("\t%d", scalar(keys %{$dataHR->{'genes'}{$geneGrp}{$annot}}));
      writeQuickOutputFile($dataHR, $infoHR, $geneGrp, $annot, "genes");
      writeQuickOutputFile($dataHR, $infoHR, $geneGrp, $annot, "trans");
      writeFullOutputFileGenes($dataHR, $infoHR, $geneGrp, $annot) if $geneGrp eq "novel_recessive" and $annot eq "candidateGenes";
      writeFullOutputFileGenes($dataHR, $infoHR, $geneGrp, $annot) if $geneGrp eq "novel_homozygous" and $annot eq "candidateGenes";
      writeFullOutputFileGenesXML($dataHR, $infoHR, $geneGrp, $annot) if $geneGrp eq "novel_homozygous" and $annot eq "candidateGenes";
    }
    print "\n";
  }
  print "\n";
  print STDERR "Done\n\n";

}

printf("%d (%0.3f) known variants in $infoHR->{'optHR'}{'varBed'} file\n", $infoHR->{'no_known_vars'}, $infoHR->{'no_known_vars'} / scalar(keys %{$dataHR->{'varHR'}}));
printf("%d (%0.3f) novel variants in $infoHR->{'optHR'}{'varBed'} file\n", $infoHR->{'no_novel_vars'}, $infoHR->{'no_novel_vars'} / scalar(keys %{$dataHR->{'varHR'}}));
print "\n\n";

printf("Number of transcripts where no gene id was found (skipped in the gene output): %d\n", $infoHR->{'skipped_no_gene_IDs'});
print "\n\n";

printf("Number of candidate genes under the recessive model: %d\n\n", scalar(keys %{$dataHR->{'genes'}{'novel_recessive'}{'candidateGenes'}}));
print "\n\n\n";


exit;






########################
########################
#  #  #  S U B S #  #  #
########################
########################

########################
# Load default values
########################

sub loadDefaults {
  my $infoHR = shift;

#  $infoHR->{'optHR'}{'minCoverage'}    = 8  unless $infoHR->{'optHR'}{'minCoverage'};
}


########################
# Write data to files
########################

sub writeQuickOutputFileVars {
  my $outFH_HR = shift;
  my $infoHR   = shift;
  my $varGrp   = shift;
  my $annot    = shift;
  my $varID    = shift;
  my $outFileType = shift;

  my $outFNtmp = $infoHR->{'outDirs'}{'vars'} . join(".", @{$infoHR->{'sampleInfo'}}{'name', 'run'}, $outFileType, $varGrp, $annot, "txt");
  $outFH_HR->{$outFNtmp} = myOpenRW($outFNtmp) unless defined($outFH_HR->{$outFNtmp});
  print {$outFH_HR->{$outFNtmp}} "$varID\n";
}

sub writeQuickOutputFile {
  my $dataHR  = shift;
  my $infoHR  = shift;
  my $geneGrp = shift;
  my $annot   = shift;
  my $outFileType = shift;
  my $outFH = myOpenRW($infoHR->{'outDirs'}{$outFileType} . join(".", @{$infoHR->{'sampleInfo'}}{'name', 'run'}, $outFileType, $geneGrp, $annot, "txt"));
  foreach my $geneID ( keys %{$dataHR->{'genes'}{$geneGrp}{$annot}} ) {
    if ( $outFileType =~ /genes/ ) {
      my $tmpDescrHR = {};
      foreach my $transHR ( @{$dataHR->{'genes'}{$geneGrp}{$annot}{$geneID}} ) { $tmpDescrHR = { map { $_ => 1 } sort keys %{$transHR->{'transcriptDescriptionHR'}} }; }
      my $concatDescription = join("\t", sort keys %{$tmpDescrHR});
      print $outFH "$geneID\t$concatDescription\n";
    } elsif ( $outFileType =~ /trans/ ) {
      foreach my $transHR ( @{$dataHR->{'genes'}{$geneGrp}{$annot}{$geneID}} ) {
        print $outFH "$transHR->{'id'}\t" . join("\t", sort keys %{$transHR->{'transcriptDescriptionHR'}}) . "\n";
      }
    }
  }
  close($outFH);
}


sub writeFullOutputFileGenes {
  my $dataHR  = shift;
  my $infoHR  = shift;
  my $geneGrp = shift;
  my $annot   = shift;

  my $fullOutFH = myOpenRW($infoHR->{'outDirs'}{'base'} . "/" . join(".", @{$infoHR->{'sampleInfo'}}{'name', 'run'}, $geneGrp, $annot, "fullOut", "txt"));
  foreach my $geneID ( keys %{$dataHR->{'genes'}{$geneGrp}{$annot}} ) {
    print $fullOutFH "$geneID\t{\n";
    foreach my $transHR ( @{$dataHR->{'genes'}{$geneGrp}{$annot}{$geneID}} ) {
      my $transID = $transHR->{'id'};
      print $fullOutFH "\t - $transID\t$transHR->{'line'}\n"; #. join("\t", @{$transHR}{'sta', 'end', 'str',}) . "\n";
      print $fullOutFH "\t - " . join("\t", values %{$dataHR->{'kgXrefHR'}{$transID}}) . "\n" if ref($dataHR->{'kgXrefHR'}{$transID}) =~ /HASH/;
      print $fullOutFH "\t - SeqNT: " . uc($transHR->{'seqHR'}{'seq'}) . "\n";
      print $fullOutFH "\t - SeqAA: " . translate(uc($transHR->{'seqHR'}{'seq'})) . "\n";
      foreach my $varID ( keys %{$transHR->{'annot'}{$annot}} ) {
        my $varHR = $transHR->{'annot'}{$annot}{$varID};
        print $fullOutFH "\t\t - $varID\t" . join("\t", $varHR->{'line'}, @{$varHR->{'info'}}{'varType', 'homHet'});
        my $kVarStr = join(",", map { "$_->{'name'}($_->{'obs'})" } @{$varHR->{'kVar'}});
        if ( ref($varHR->{'info'}{'ns'}{$transID}) =~ /HASH/ ) {
          print $fullOutFH "\t" . join("\t", @{$varHR->{'info'}{'ns'}{$transID}{'ref'}}{'exonNo', 'relSta', 'snpFrame', 'genBase', 'refCdn', 'refAA'}) . ($kVarStr ? "\t$kVarStr\n" : "\n");
          foreach my $allele ( sort keys %{$varHR->{'info'}{'ns'}{$transID}{'alt'}} ) {
            print $fullOutFH "\t\t\t - $allele\t" . join("\t", @{$varHR->{'info'}{'ns'}{$transID}{'alt'}{$allele}}{'snpCdn', 'snpAA', 'pph_prediction'}) . "\tpphInput: $varHR->{'chr'}:" . ($varHR->{'sta'} + 1) . "\t$varHR->{'info'}{'ref'}/$allele\n";
            print $fullOutFH "\t\t\t - SeqNT: " . uc($varHR->{'info'}{'ns'}{$transID}{'alt'}{$allele}{'mutSeq'}) . "\n";
            print $fullOutFH "\t\t\t - SeqAA: " . translate(uc($varHR->{'info'}{'ns'}{$transID}{'alt'}{$allele}{'mutSeq'})) . "\n";
          }
        } else {
          print $fullOutFH ($kVarStr ? "\t$kVarStr\n" : "\n");
        }
      }
      print $fullOutFH "\n";
    }
    print $fullOutFH "}\n\n";
  }
  close($fullOutFH);
}


sub readKgXref {
  my $kgXrefFN = shift || "/home/magnus.bjursell/hg/knownGene_GRCh37_hg19/kgXref.txt";
  my $kgXrefHR = {};
  my $kgXrefFH = myOpen($kgXrefFN);
  while ( my $line = <$kgXrefFH> ) {
    chomp($line);
    my ($key, @rest) = split(/\t/, $line);
    $kgXrefHR->{$key} = [@rest];
  }
  close($kgXrefFH);
  return $kgXrefHR;
}


sub writeAnnotVarsCombinationsToFile {
  my $infoHR  = shift;

  my $outCombFH = myOpenRW($infoHR->{'outDirs'}{'base'} . "/combinations_of_annotation_tags.txt");
  print $outCombFH "Variant annotation combinations:\n";
  foreach my $annotComb (sort { $infoHR->{'varAnnotComb'}{$b} <=> $infoHR->{'varAnnotComb'}{$a} } keys %{$infoHR->{'varAnnotComb'}} ) { print $outCombFH "$annotComb\t$infoHR->{'varAnnotComb'}{$annotComb}\n"; }
  print $outCombFH "\n\n\n";
  close($outCombFH);
}


sub printNoVarsTypes {
  my $infoHR  = shift;

  print "Number of each variant types\n";
  foreach my $varType ( sort keys %{$infoHR->{'no_varTypesHR'}} ) { print "$varType\t$infoHR->{'no_varTypesHR'}{$varType}\n"; }
  print "\n\n";
}


###############################################
# Load variation data (SNPs, Indels) from file
###############################################

sub loadVariationData { 
  my $dataHR = shift;
  my $bedFN  = shift;

  print STDERR "Loading $bedFN variation file\n";
  read_bed_variant_file($dataHR, $bedFN, 1);
  print STDERR "   - Loaded " . scalar(keys %{$dataHR->{'varHR'}}) . " variants, $dataHR->{'varInfo'}{'cntVarType'}{'snp'} snps and $dataHR->{'varInfo'}{'cntVarType'}{'indel'} indels\n";
  print STDERR "Done\n\n";
}


###############################################
# Parse gene line from gene input file
###############################################

sub parseGene {
  my $line = shift;
  my $ft   = shift;
  my $hdr  = getFileHeader($ft);

  chomp($line);
  my @cols = split(/\t/, $line);
  die "\nNumber of columns do not match the header specifications for this file format for line:\n" . join("\t", sort { $hdr->{$a} <=> $hdr->{$b} } keys %{$hdr}) . "\n$line\n\n" unless scalar(@cols) == scalar(keys %{$hdr});

  my $transHR = { map { $_ => $cols[$hdr->{$_}] } sort { $hdr->{$a} <=> $hdr->{$b} } keys %{$hdr} };
  $transHR->{'eStaA'} = [split(/\,/, $transHR->{'eSta'})];
  $transHR->{'eEndA'} = [split(/\,/, $transHR->{'eEnd'})];
  $transHR->{'line'}  = $line;

  return $transHR;
}


###############################################
# Quick CGS length calculator
###############################################

sub returnCDSlen {
  my $transHR = shift;
  my $len = 0;
  for my $eNo (0 .. scalar(@{$transHR->{'eStaA'}}) - 1 ) {
    my $start = max($transHR->{'eStaA'}[$eNo], $transHR->{'cSta'});
    my $end   = min($transHR->{'eEndA'}[$eNo], $transHR->{'cEnd'});
    next if $start >= $end;
    $len += $end - $start;
  }
  return $len;
}


###############################################
# Parse DNA sequence
###############################################

sub parseSeq {
  my $dataHR  = shift;
  my $transHR = shift;

  my $seqHR = {};
  $seqHR->{'sumIntrLen'}  = [(0) x scalar(@{$transHR->{'eStaA'}})];

  for my $eNo (0 .. scalar(@{$transHR->{'eStaA'}}) - 1 ) {

    my $start = max($transHR->{'eStaA'}[$eNo], $transHR->{'cSta'});
    my $end   = min($transHR->{'eEndA'}[$eNo], $transHR->{'cEnd'});
    next if $start >= $end;
    my $len = $end - $start;

    $seqHR->{'seq'} .= uc( returnSequenceFromFasta($dataHR->{'faiHR'}, $transHR->{'chr'}, $start, $len) );
    $seqHR->{'sumIntrLen'}[$eNo + 1]  = ($seqHR->{'sumIntrLen'}[$eNo] + ($transHR->{'eStaA'}[$eNo + 1] - $transHR->{'eEndA'}[$eNo])) if $eNo < scalar(@{$transHR->{'eStaA'}}) - 1;
  }
  reverseComplement(\$seqHR->{'seq'}) if $transHR->{'str'} eq "-";
  warn "Sequence length (" . length($seqHR->{'seq'}) . ") *** IS NOT *** evenly diviceble by 3 for:\n$transHR->{'line'}\n\n" unless length($seqHR->{'seq'}) % 3 == 0;
  $transHR->{'lenNT'} = length($seqHR->{'seq'});
  $transHR->{'lenAA'} = ( length($seqHR->{'seq'}) - 3 ) / 3;
  return $seqHR;
}


###############################################
# Parse line from varians file (dbSNP)
###############################################

sub parseKVar {
  my $line = shift;
  my $ft   = shift;
  my $hdr  = getFileHeader($ft);

  chomp($line);
  my @cols = split(/\t/, $line);
  die "\nNumber of columns do not match the header specifications for this file format for line:\n" . join("\t", sort { $hdr->{$a} <=> $hdr->{$b} } keys %{$hdr}) . "\n$line\n\n" unless scalar(@cols) == scalar(keys %{$hdr});

  my $kVarHR = { map { $_ => $cols[$hdr->{$_}] } keys %{$hdr} };
  my $currLen = $kVarHR->{'end'} - $kVarHR->{'sta'};
  my $obsMaxLen = $currLen;
  $kVarHR->{'obsHR'}{$kVarHR->{'refNCBI'}} = 1 if $kVarHR->{'refNCBI'} =~ /^[ACGT]$/i;
  $kVarHR->{'obsHR'}{$kVarHR->{'refUCSC'}} = 1 if $kVarHR->{'refUCSC'} =~ /^[ACGT]$/i;
  foreach my $obs ( split(/\//, $kVarHR->{'obs'}) ) {
    $obs = "*" if $obs eq "-";
    reverseComplement(\$obs) if $kVarHR->{'str'} eq "-";
    $kVarHR->{'obsHR'}{$obs} = 1;
    $obsMaxLen = max($obsMaxLen, length($obs));
  }
  $kVarHR->{'line'}    = $line;
  $kVarHR->{'varFile'} = $ft;
  if ( $kVarHR->{'class'} ne "single" and $obsMaxLen > $currLen ) {
    $kVarHR->{'sta'} -= ($obsMaxLen - $currLen);
    $kVarHR->{'end'} += ($obsMaxLen - $currLen);
  }

  return $kVarHR;
}


###############################################
# Connect variant to gene
###############################################

sub tieVarToGene {
  my $dataHR  = shift;
  my $varHR   = shift;
  my $transHR = shift;
  my $annot   = shift;

  $varHR->{'annot'}{$annot}{$transHR->{'id'}} = $transHR unless exists($varHR->{'annot'}{$annot}{$transHR->{'id'}});
  $transHR->{'annot'}{$annot}{$varHR->{'id'}} = $varHR unless exists($transHR->{'annot'}{$annot}{$varHR->{'id'}});

  addTransGrp($varHR, $transHR, $annot, "all");

  if ( ref($varHR->{'kVar'}) =~ /ARRAY/ and scalar(@{$varHR->{'kVar'}}) > 0 ) {
    addTransGrp($varHR, $transHR, $annot, "known");
  } else {
    addTransGrp($varHR, $transHR, $annot, "novel");
  }

  $dataHR->{'transHR'}{$transHR->{'id'}} = $transHR unless exists($dataHR->{'transHR'}{$transHR->{'id'}});
  $dataHR->{'debug'}{'annotCount'}{$annot}{$varHR} = 1;
}

sub addTransGrp {
  my $varHR    = shift;
  my $transHR  = shift;
  my $annot    = shift;
  my $transGrp = shift;

  $transHR->{'transGrp'}{$transGrp}{$annot}{'any'}++;
  $transHR->{'transGrp'}{$transGrp}{$annot}{'sum'} += ($varHR->{'info'}{'homHet'} eq "hom" ? 1 : 0.5);
  $transHR->{'transGrp'}{$transGrp}{$annot}{$varHR->{'info'}{'homHet'}}++;
}


###############################################
# Detect overlap between gene and variant
###############################################

sub detectTransOverlap {
  my $dataHR  = shift;
  my $varHR   = shift;
  my $transHR = shift;

  if ( $transHR->{'cSta'} == $transHR->{'cEnd'} ) {
    tieVarToGene($dataHR, $varHR, $transHR, "ncGene") if returnOverlap([$transHR->{'sta'}, $transHR->{'end'} - 1, $varHR->{'sta'}, $varHR->{'end'} - 1]);
    return;
  }

  # Upstream
  if ( ( $transHR->{'str'} eq "+" and returnOverlap([$transHR->{'sta'} - $dataHR->{'paramHR'}{'upsLen'}, $transHR->{'sta'} - 1, $varHR->{'sta'}, $varHR->{'end'} - 1]) ) or
       ( $transHR->{'str'} eq "-" and returnOverlap([$transHR->{'end'}, $transHR->{'end'} + $dataHR->{'paramHR'}{'upsLen'} - 1, $varHR->{'sta'}, $varHR->{'end'} - 1]) ) ) {
    tieVarToGene($dataHR, $varHR, $transHR, "upStr");
  }

  # Downstream
  if ( ( $transHR->{'str'} eq "-" and returnOverlap([$transHR->{'sta'} - $dataHR->{'paramHR'}{'dnsLen'}, $transHR->{'sta'} - 1, $varHR->{'sta'}, $varHR->{'end'} - 1]) ) or
       ( $transHR->{'str'} eq "+" and returnOverlap([$transHR->{'end'}, $transHR->{'end'} + $dataHR->{'paramHR'}{'dnsLen'} - 1, $varHR->{'sta'}, $varHR->{'end'} - 1]) ) ) {
    tieVarToGene($dataHR, $varHR, $transHR, "dnStr");
  }


  return unless returnOverlap([$transHR->{'sta'}, $transHR->{'end'} - 1, $varHR->{'sta'}, $varHR->{'end'} - 1]);
  # Coding RNA
  tieVarToGene($dataHR, $varHR, $transHR, "cGene");

  # 5' UTRs
  if ( returnOverlap([$transHR->{'sta'}, $transHR->{'cSta'} - 1, $varHR->{'sta'}, $varHR->{'end'} - 1]) ) {
    tieVarToGene($dataHR, $varHR, $transHR, "5pUTR");
  }

  # 3' UTRs
  if ( returnOverlap([$transHR->{'cEnd'}, $transHR->{'end'} - 1, $varHR->{'sta'}, $varHR->{'end'} - 1]) ) {
    tieVarToGene($dataHR, $varHR, $transHR, "3pUTR");
  }

  for my $eNo (0 .. scalar( @{$transHR->{'eStaA'}} ) - 1 ) {

    # Introns
    if ( $eNo > 0 and returnOverlap([ $transHR->{'eEndA'}[$eNo - 1], $transHR->{'eStaA'}[$eNo] - 1, $varHR->{'sta'}, $varHR->{'end'} - 1]) ) {
      tieVarToGene($dataHR, $varHR, $transHR, "intron");
    }

    # Exons
    if ( returnOverlap([ $transHR->{'eStaA'}[$eNo], $transHR->{'eEndA'}[$eNo] - 1, $varHR->{'sta'}, $varHR->{'end'} - 1]) ) {
      tieVarToGene($dataHR, $varHR, $transHR, "exon");
    }

    # Donor splice sites (plus dir)
    if ( $eNo < scalar(@{$transHR->{'eStaA'}}) - 1 and returnOverlap([ $transHR->{'eEndA'}[$eNo], $transHR->{'eEndA'}[$eNo] + $dataHR->{'paramHR'}{'ssLen'} - 1, $varHR->{'sta'}, $varHR->{'end'} - 1]) ) {
      tieVarToGene($dataHR, $varHR, $transHR, "ss");
      tieVarToGene($dataHR, $varHR, $transHR, "candidateGenes");
    }

    # Acceptor splice sites (plus dir)
    if ( $eNo > 0 and returnOverlap([ $transHR->{'eStaA'}[$eNo] - $dataHR->{'paramHR'}{'ssLen'}, $transHR->{'eStaA'}[$eNo] - 1, $varHR->{'sta'}, $varHR->{'end'} - 1]) ) {
      tieVarToGene($dataHR, $varHR, $transHR, "ss");
      tieVarToGene($dataHR, $varHR, $transHR, "candidateGenes");
    }

    # Coding sequence including non-synonymous changes
    next if $varHR->{'end'} - 1 < $transHR->{'cSta'} or $varHR->{'sta'} > $transHR->{'cEnd'} - 1;
    my $start = max($transHR->{'eStaA'}[$eNo], $transHR->{'cSta'});
    my $end   = min($transHR->{'eEndA'}[$eNo], $transHR->{'cEnd'});
    next if $start >= $end;
  
    # Coding exons
    if ( returnOverlap([ $start, $end - 1, $varHR->{'sta'}, $varHR->{'end'} - 1]) ) {
      tieVarToGene($dataHR, $varHR, $transHR, "cds");
      if ( $varHR->{'info'}{'varType'} eq "snp" ) {
        detectNsTrans($dataHR, $varHR, $transHR, $eNo);
      } elsif ( $varHR->{'info'}{'varType'} eq "indel" ) {
        tieVarToGene($dataHR, $varHR, $transHR, "cIndel");
        tieVarToGene($dataHR, $varHR, $transHR, "candidateGenes");
      }

    }

  }
}


sub detectNsTrans {
  my $dataHR  = shift;
  my $varHR   = shift;
  my $transHR = shift;
  my $eNo     = shift;

  $transHR->{'seqHR'} = parseSeq($dataHR, $transHR);

  my $snvInfoHR = {};
  my $snvInfoHR->{'genomeBase'} = uc(returnSequenceFromFasta($dataHR->{'faiHR'}, $transHR->{'chr'}, $varHR->{'sta'}, 1));
  die "\n$varHR->{'line'} is not a single nucleotide substitution\n\n" unless $varHR->{'info'}{'varType'} eq "snp";
  die "\nSNP reference and genomic base do not match:\nGenomic: $snvInfoHR->{'genomeBase'}; Variation: $varHR->{'line'};\n\n" unless uc($varHR->{'info'}{'ref'}) eq $snvInfoHR->{'genomeBase'};
  
  if ( $transHR->{'str'} eq "+" ) {
    $snvInfoHR->{'relStart'}    = $varHR->{'sta'} - ($transHR->{'cSta'} + $transHR->{'seqHR'}{'sumIntrLen'}[$eNo]);
    $snvInfoHR->{'snpFrame'}    = $snvInfoHR->{'relStart'} % 3;
    $snvInfoHR->{'relCdnStart'} = $snvInfoHR->{'relStart'} - $snvInfoHR->{'snpFrame'};
    $snvInfoHR->{'refCdn'}      = uc( substr($transHR->{'seqHR'}{'seq'}, $snvInfoHR->{'relCdnStart'}, 3) );
    $snvInfoHR->{'refAA'}       = translate($snvInfoHR->{'refCdn'});
    $snvInfoHR->{'eNo'}         = $eNo;

  } else {
    $snvInfoHR->{'relStart'}    = ( length($transHR->{'seqHR'}{'seq'}) - 1 ) - ( $varHR->{'sta'} - ($transHR->{'cSta'} + $transHR->{'seqHR'}{'sumIntrLen'}[$eNo]) );
    $snvInfoHR->{'snpFrame'}    = $snvInfoHR->{'relStart'} % 3;
    $snvInfoHR->{'relCdnStart'} = $snvInfoHR->{'relStart'} - $snvInfoHR->{'snpFrame'};
    $snvInfoHR->{'refCdn'}      = uc( substr($transHR->{'seqHR'}{'seq'}, $snvInfoHR->{'relCdnStart'}, 3) );
    $snvInfoHR->{'refAA'}       = translate($snvInfoHR->{'refCdn'});
    $snvInfoHR->{'eNo'}         = (scalar(@{$transHR->{'eStaA'}}) - 1) - $eNo;

  }
  $varHR->{'info'}{'ns'}{$transHR->{'id'}}{'ref'} = {'exonNo' => $snvInfoHR->{'eNo'}, 'relSta' => $snvInfoHR->{'relStart'}, 'relStaAA_1' => int($snvInfoHR->{'relStart'} / 3) + 1, 'snpFrame' => $snvInfoHR->{'snpFrame'},
                                                     'genBase' => $snvInfoHR->{'genomeBase'}, 'refCdn' =>  $snvInfoHR->{'refCdn'}, 'refAA' =>  $snvInfoHR->{'refAA'}};

  foreach my $allele ( @{$varHR->{'allelesAR'}} ) {
    die "\n$allele from $varHR->{'line'} is not a valid single nucleotide substitution\n\n" unless $allele =~ /^[ACGT]$/i;
    next if $allele eq $varHR->{'info'}{'ref'};
    reverseComplement(\$allele) if $transHR->{'str'} eq "-";

    $snvInfoHR->{'snpCdn'} = $snvInfoHR->{'refCdn'};
    $snvInfoHR->{'snpCdn'} =~ s/(\w{$snvInfoHR->{'snpFrame'}})\w/$1$allele/;
    $snvInfoHR->{'snpAA'}  = translate($snvInfoHR->{'snpCdn'});

    warn "\nRef or SNP AA missing (RefAA: $snvInfoHR->{'refAA'}; SnpAA: $snvInfoHR->{'snpAA'};) GeneLen: " . length($transHR->{'seqHR'}{'seq'}) . "\n" . join("; ", map { $_ . ":" . $snvInfoHR->{$_} } keys %{$snvInfoHR}) . "\n   - $transHR->{'line'}\n   - $varHR->{'line'}\n\n" unless $snvInfoHR->{'refAA'} and $snvInfoHR->{'snpAA'};

    if ( $snvInfoHR->{'snpAA'} ne $snvInfoHR->{'refAA'} ) {
      my $mutSeq = uc($transHR->{'seqHR'}{'seq'}); my $trash = substr($mutSeq, $snvInfoHR->{'relStart'}, 1, $allele);

      tieVarToGene($dataHR, $varHR, $transHR, "ns");
      tieVarToGene($dataHR, $varHR, $transHR, "candidateGenes");
      tieVarToGene($dataHR, $varHR, $transHR, "nsense") if $snvInfoHR->{'snpAA'} eq "X";
      tieVarToGene($dataHR, $varHR, $transHR, "mutStop") if $snvInfoHR->{'refAA'} eq "X" and $snvInfoHR->{'snpAA'} ne "X";

      warn "\nDuplicate gene with identical ($transHR->{'id'}) id overlap same SNP ($varHR->{'line'})\n\n" if exists($varHR->{'info'}{'ns'}{$transHR->{'id'}}{'alt'}{$allele});
      $varHR->{'info'}{'ns'}{$transHR->{'id'}}{'alt'}{$allele} = {'allele' =>  $allele, 'snpCdn' =>  $snvInfoHR->{'snpCdn'}, 'snpAA' =>  $snvInfoHR->{'snpAA'}, 'mutSeq' => $mutSeq};

      my $refHR = $varHR->{'info'}{'ns'}{$transHR->{'id'}}{'ref'};
      my $altHR = $varHR->{'info'}{'ns'}{$transHR->{'id'}}{'alt'};
      $varHR->{'info'}{'ns'}{$transHR->{'id'}}{'alt'}{$allele}{'ntMutStr'} = sprintf("%s%d%s\t", ($transHR->{'str'} eq "-" ? complement($refHR->{'genBase'}) : $refHR->{'genBase'}), $refHR->{'relSta'} + 1, $allele);
      $varHR->{'info'}{'ns'}{$transHR->{'id'}}{'alt'}{$allele}{'aaMutStr'} = sprintf("%s%d%s", $refHR->{'refAA'}, int($refHR->{'relSta'} / 3) + 1, $altHR->{$allele}{'snpAA'});


#      my $mutStrAA = join("", @{$varHR->{'info'}{'ns'}{$transHR->{'id'}}{'ref'}}{'refAA', 'relStaAA_1'}) . $varHR->{'info'}{'ns'}{$transHR->{'id'}}{'alt'}{$allele}{'snpAA'};
#      $transHR->{'info'}{'ns'}{'mutStrHR'}{$mutStrAA}{'homHet'} = $varHR->{'info'}{'homHet'};
    }
  }

}


# Search against known variants

sub tieVarToKVar {
  my $dataHR  = shift;
  my $varHR   = shift;
  my $kVarHR = shift;

  if ( ref($varHR->{'kVar'}) =~ /ARRAY/ ) { foreach my $addedKVarHR ( @{$varHR->{'kVar'}} ) { return if $addedKVarHR->{'varFile'} eq $kVarHR->{'varFile'}; } }
  push(@{$varHR->{'kVar'}}, $kVarHR);
}


sub detectKnownVarOverlap {
  my $dataHR = shift;
  my $varHR  = shift;
  my $kVarHR = shift;


  if ( $varHR->{'info'}{'varType'} eq "snp" ) { # and $kVarHR->{'class'} eq "single" 
    if ( returnOverlap([$kVarHR->{'sta'}, $kVarHR->{'end'} - 1, $varHR->{'sta'}, $varHR->{'end'} - 1]) ) {
      foreach my $allele ( @{$varHR->{'allelesAR'}} ) {
        next if $allele eq $varHR->{'info'}{'ref'};
        return unless $kVarHR->{'obsHR'}{$allele} > 0;
      }
      tieVarToKVar($dataHR, $varHR, $kVarHR);
    }

  } elsif ( $varHR->{'info'}{'varType'} eq "indel" ) { # and $kVarHR->{'class'} ne "single"
    if ( returnOverlap([$kVarHR->{'sta'}, $kVarHR->{'end'} - 1, $varHR->{'sta'}, $varHR->{'end'} - 1]) ) {
      foreach my $allele ( @{$varHR->{'allelesAR'}} ) {
        next if $allele eq $varHR->{'info'}{'ref'} or $allele eq "*";
        $allele =~ s/[\+\-]//;
        return unless $kVarHR->{'obsHR'}{$allele} > 0;
      }
      tieVarToKVar($dataHR, $varHR, $kVarHR);
    }
  }
}


# Parse sample name from input bed file

sub parseBedFileName {
  my $bedFN = shift;
  my $smplInfoHR = {};
  if ( $bedFN =~ /\/?([^\.\/]+)\.([^\.]+)\./ ) {
    $smplInfoHR = {'name' => $1, 'run' => $2,};
  } else {
    $smplInfoHR = {'name' => 'sample', 'run' => 'run',};
    warn "Cannot parse sample name and run name, setting to \'sample\' and \'run\', respectively\n";
  }
  return $smplInfoHR;
}


##############################################
# Setup output structure and parse sample name
##############################################

sub setupOutput_parseSampleName { 
  my $infoHR = shift;
  my $outDN = shift;
  my $bedFN = shift;

  print STDERR "\n";
  print STDERR "Seting up directory structure: \n";
  $infoHR->{'outDirs'} = setupOutputDirectoryStructure($outDN);
  print STDERR "   - " . join("\n   - ", sort values %{$infoHR->{'outDirs'}}) . "\n\n";
  print STDERR "Done\n\n";

  print STDERR "Parsing input file name: \n";
  $infoHR->{'sampleInfo'} = parseBedFileName($bedFN);
  print STDERR "   - Sample name: $infoHR->{'sampleInfo'}{'name'}\n   - Run name: $infoHR->{'sampleInfo'}{'run'}\n";
  print STDERR "Done\n\n";
}


##############################################
# Setup output directory structure
##############################################

sub setupOutputDirectoryStructure {
  my $outDN = shift;
  my $tmpHR = {};
  $outDN =~ s/\/+$//;
  foreach my $cDir ( split(/\//, $outDN) ) {
    if ( -d "$tmpHR->{'cDirLevel'}$cDir" ) {
      removeFiles("$tmpHR->{'cDirLevel'}$cDir");
    } else {
      mkdir("$tmpHR->{'cDirLevel'}$cDir") or die "\nCannot create output dir $tmpHR->{'cDirLevel'}$cDir ($!)\n\n"; 
    }
    $tmpHR->{'cDirLevel'} .= "$cDir/";
  }
  die "\nOutput directory ($outDN) is not a valid directory\n\n" unless -d $outDN;
  my $dirsHR = {'vars' => $outDN . '/variantBased/', 'trans' => $outDN . '/transcriptBased/', 'genes' => $outDN . '/geneBased/', 'pph' => $outDN . '/polyPhenData/',};
  foreach my $dir ( keys %{$dirsHR} ) {
    if ( -d $dirsHR->{$dir} ) { removeFiles($dirsHR->{$dir}); } else { mkdir($dirsHR->{$dir}) or die "\nCannot create directory $dirsHR->{$dir} ($!)\n\n"; }
  }
  $dirsHR->{'base'} = $outDN . "/"; 
  return $dirsHR;
}

sub removeFiles {
  my $dn = shift; $dn =~ s/\/+$//;
return;
  opendir(my $dh, $dn) or die "\nCannot open directory $dn ($!)\n\n";
  while ( my $file = readdir($dh) ) {
    #unlink("$dn/$file") or die "\nCannot remove file $file located in dir $dn ($!)\n\n";
  }
  closedir($dh);
}


##############################################
# Debug print variant data
##############################################

sub debugPrintVarData {
  my $dataHR = shift;

  print STDERR "Debug print annot cnts\n";
  foreach my $annot ( sort keys %{$dataHR->{'debug'}{'annotCount'}} ) { print STDERR "$annot\t" . scalar(keys %{$dataHR->{'debug'}{'annotCount'}{$annot}}) . "\n"; }
  print STDERR "\n\n";
}


##############################################
# Polyphen
##############################################

sub runPolyPhen2_parseRet {
  my $dataHR  = shift;
  my $infoHR  = shift;

  my $startTime = time;
  my $maxWait = 100000;
#  local $SIG{INT} = sub { print STDERR "Skipping polyphen-2 analysis (ctrl-c pressed), please wait...\n"; $maxWait = 0; };
  local $| = 1;
  print STDERR "Running polyphen-2 using $dataHR->{'pphInfoHR'}{'pphFile'}{'name'} as input, press \"q\" (and wait) to skip (max wait time: " . round($maxWait / 60) . " min)\n";
  my $pphSubmitCmd = "curl -o $infoHR->{'outDirs'}{'pph'}submit.html -c $infoHR->{'outDirs'}{'pph'}submit.cks -f " .
                     "--form _ggi_project=PPHWeb2 " .
                     "--form _ggi_origin=query " .
                     "--form _ggi_batch= " .
                     "--form _ggi_batch_file=\@$dataHR->{'pphInfoHR'}{'pphFile'}{'name'} " .
                     "--form description= " .
                     "--form NOTIFYME= " .
                     "--form uploaded_sequences_1= " .
                     "--form description_of_uploaded_sequences_1= " .
                     "--form MODELNAME=HumDiv " .
                     "--form UCSCDB=hg19 " .
                     "--form SNPFILTER=1 " .
                     "--form SNPFUNC=c " .
                     "--form _ggi_target_pipeline=Submit%20Batch " .
                     "http://genetics.bwh.harvard.edu/cgi-bin/ggi/ggi2.cgi " .
                     "2> $infoHR->{'outDirs'}{'pph'}curl.submit.stderr ";

  my $pphSubmitRet = `$pphSubmitCmd`;
  my $exitCode = $?>> 8; `echo "\ncurl exit code: $exitCode\n" >> $infoHR->{'outDirs'}{'pph'}curl.submit.stderr`;
  sleep 60;
  warn "Error running polyphen, cannot submit job please see $infoHR->{'outDirs'}{'pph'}submit.html or $infoHR->{'outDirs'}{'pph'}curl.submit.stderr for details" and return unless ($exitCode) == 0;

  my @pphSubmitCks = split(/\s+/, `tail -n 1 $infoHR->{'outDirs'}{'pph'}submit.cks`);
  warn "Error in polyphen processing, please see $infoHR->{'outDirs'}{'pph'}submit.html for details\n\n" and return if $pphSubmitCks[5] ne "polyphenweb2";
#.genetics.bwh.harvard.edu	TRUE	/cgi-bin/ggi	FALSE	1295964557	polyphenweb2	93b87e0382430cd36b19995065c467e134f17d51
  my $pphSubmitSID = pop(@pphSubmitCks);
  print STDERR "Polyphen-2 submit ID: $pphSubmitSID\n";
  my $pphCheckCmd  = "curl -o $infoHR->{'outDirs'}{'pph'}manage.html -c $infoHR->{'outDirs'}{'pph'}manage.cks -b $infoHR->{'outDirs'}{'pph'}submit.cks -f " .
                     "--form sidreset=1 " .
                     "--form _ggi_project=PPHWeb2 " .
                     "--form _ggi_origin=manage " .
                     "--form _ggi_target_manage=Refresh " .
                     "http://genetics.bwh.harvard.edu/cgi-bin/ggi/ggi2.cgi " .
                     "2> $infoHR->{'outDirs'}{'pph'}curl.manage.stderr "; 
  my $pphCheckFix  = qq{sed -i 's/=\\"\\//=\\"http:\\/\\/genetics\\.bwh\\.harvard\\.edu\\//g'  $infoHR->{'outDirs'}{'pph'}manage.html };

  my $pphGetLogCmd = "curl -o $infoHR->{'outDirs'}{'pph'}log.txt -c $infoHR->{'outDirs'}{'pph'}log.cks -f http://genetics.bwh.harvard.edu/ggi/pph2/$pphSubmitSID/1/pph2-log.txt 2> $infoHR->{'outDirs'}{'pph'}curl.log.stderr";
  my $pphGetRetCmd = "curl -o $infoHR->{'outDirs'}{'pph'}results.txt -c $infoHR->{'outDirs'}{'pph'}results.cks -f http://genetics.bwh.harvard.edu/ggi/pph2/$pphSubmitSID/1/pph2-full.txt 2> $infoHR->{'outDirs'}{'pph'}curl.results.stderr";

  my $curlCmdFH = myOpenRW("$infoHR->{'outDirs'}{'pph'}curl.allCmds.txt"); print $curlCmdFH join("\n\n", $pphSubmitCmd, $pphCheckCmd, $pphCheckFix, $pphGetLogCmd, $pphGetRetCmd); close($curlCmdFH);

  my $pphCheckRet = ""; my $pphGetRetRet = ""; my $pphGetLogRet = ""; my $pphCheckFixRet = "";
  ReadMode(3);
  while ( 1 ) {
    $pphGetRetRet = `$pphGetRetCmd`;
    my $exitCode = $?>> 8; `echo "\ncurl exit code: $exitCode\n" >> $infoHR->{'outDirs'}{'pph'}curl.results.stderr`;
    last if ($exitCode) == 0;

    $pphCheckRet = `$pphCheckCmd`;
    my $exitCode = $?>> 8; `echo "\ncurl exit code: $exitCode\n" >> $infoHR->{'outDirs'}{'pph'}curl.manage.stderr`;
    $pphCheckFixRet = `$pphCheckFix`;

    $pphGetLogRet = `$pphGetLogCmd`;
    my $exitCode = $?>> 8; `echo "\ncurl exit code: $exitCode\n" >> $infoHR->{'outDirs'}{'pph'}curl.log.stderr`;

    my $soFar = 0;
    while ( 1 ) {
      my $inKey = ReadKey(-1);
      if ( $inKey =~ /q/ ) {
        ReadMode(0);
        warn "\nSkipping polyphen-2 analysis (\"q\" pressed), please wait...\n\n";
        return undef;
      }
      sleep 2;
      $soFar += 2;
      last if $soFar >= 60;
    }

    warn "\nWaited for " . round($maxWait / 60) . " minutes, skipping polyphen-2 analysis\n\n" and return if time - $startTime >= $maxWait;

    print STDERR ".";
  }
  ReadMode(0);
  print STDERR "\n\n";

  sleep 30; # Make sure all is done at the server side
  $pphCheckRet  = `$pphCheckCmd`;
  my $exitCode = $?>> 8; `echo "\ncurl exit code: $exitCode\n" >> $infoHR->{'outDirs'}{'pph'}curl.manage.stderr`;
  $pphCheckFixRet = `$pphCheckFix`;
  $pphGetLogRet = `$pphGetLogCmd`;
  my $exitCode = $?>> 8; `echo "\ncurl exit code: $exitCode\n" >> $infoHR->{'outDirs'}{'pph'}curl.log.stderr`;
  $pphGetRetRet = `$pphGetRetCmd`;
  my $exitCode = $?>> 8; `echo "\ncurl exit code: $exitCode\n" >> $infoHR->{'outDirs'}{'pph'}curl.results.stderr`;
  warn "Error running polyphen, cannot download the results file, please see $infoHR->{'outDirs'}{'pph'}log.txt for details" and return unless ($?>> 8) == 0;

  my @pphRetLines  = split(/\n/, `cat $infoHR->{'outDirs'}{'pph'}results.txt`);
  my @header = split(/\s*\t\s*/, shift(@pphRetLines)); $header[0] =~ s/^#//;

  my $pphParseHR = {};
  foreach my $pphLine ( @pphRetLines ) {
    next if $pphLine =~ /^#/;
    my @cols = split(/\s*\t\s*/, $pphLine);
    my $pphRetHR = { map { $_ => shift(@cols) } @header };
    my ($pphVarID, $pphAnnot, $pphTransGrp) = split(/;/, $pphRetHR->{'Comments'});
    my ($pphChr, $pphPos, $pphBases, $knownGeneID) = ( $pphRetHR->{'o_snp_id'} =~ /(chr[^:]+):(\d+)\.(\w+)\.(uc\w+\.\d+)/ ); #chr12:51696879.GC.uc001ryg.2

    warn "Polyphen-2 error: Cannot find transcript $knownGeneID\n" and next unless defined($dataHR->{'transHR'}{$knownGeneID});
    warn "Polyphen-2 error: Cannot find variant $pphVarID\n" and next unless defined($dataHR->{'varHR'}{$pphVarID});
    warn "Polyphen-2 error: Multiple entries for the $knownGeneID (transcript), $pphRetHR->{'o_snp_id'} <=> $pphParseHR->{$pphRetHR->{'o_snp_id'}}\n" and next if defined($pphParseHR->{$pphRetHR->{'o_snp_id'}});

    my $transHR = $dataHR->{'transHR'}{$knownGeneID};
    my $varHR   = $dataHR->{'varHR'}{$pphVarID};

    $pphParseHR->{$pphRetHR->{'o_snp_id'}} = 1;
    addTransGrp($varHR, $transHR, $pphAnnot, $pphTransGrp . '_polyphen-2') unless $pphRetHR->{'prediction'} eq "benign";
    for my $allele ( @{$pphRetHR}{'nt1', 'nt2'} ) {
      next unless defined($varHR->{'info'}{'ns'}{$knownGeneID}{'alt'}{$allele});
      $varHR->{'info'}{'ns'}{$knownGeneID}{'alt'}{$allele}{'pph_prediction'} = $pphRetHR->{'prediction'};
#      my $mutStrAA = join("", @{$varHR->{'info'}{'ns'}{$transHR->{'id'}}{'ref'}}{'refAA', 'relStaAA_1'}) . $varHR->{'info'}{'ns'}{$transHR->{'id'}}{'alt'}{$allele}{'snpAA'};
#      $transHR->{'info'}{'ns'}{'mutStrHR'}{$mutStrAA}{'pph_prediction'} = $pphRetHR->{'prediction'} if ref($transHR->{'info'}{'ns'}{'mutStrHR'}{$mutStrAA}) =~ /HASH/;
    }
  }

  print STDERR "\nDone\n\n";
}




#__END__


=head1 NAME

annotate_variants.pl - Annotate variants

=head1 SYNOPSIS

annotate_variants.pl [options]

  Options:
   --help            Brief help message (not yet implemented)
   --verbose         Write some additional output (not yet implemented)
   --varBed          Input variant bed file [required]
   --outdir          Output directory [required]
   --outdb           Output database

  The input variant bed file should use this specific format:
  chr1    866318  866319  G=>A(62;A:62)   214     +
  chr1    866510  866510  C=>*/+CCCT(99;*:68,+CCCT:31)    2021    +
  chr1    870902  870903  T=>C(27;A:1,C:26)       105     +
  chr1    871041  871041  C=>*/+A(133;*:70,+A:63) 2063    +
  chr1    876498  876499  A=>G(36;D:4,G:32)       123     +
  chr1    884100  884101  A=>M(53;A:37,C:14,D:2)  134     +

Hopefully I will document this format better in the future...

=cut

#=head1 OPTIONS

#=over 8

#=item B<-help>

#Print a brief help message and exits.

#=back

#=head1 DESCRIPTION




