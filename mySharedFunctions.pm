package mySharedFunctions;

# (c) 2008 - Magnus Bjursell

use strict;
use warnings;
no warnings qw(uninitialized);
use feature ":5.10";
use FindBin;

use lib $FindBin::Bin;
use mySharedVariables qw(:basic);


BEGIN {
  use Exporter ();
  our ($VERSION, @ISA, @EXPORT, @EXPORT_OK, %EXPORT_TAGS);
  @ISA = qw(Exporter);
  $VERSION     = 0.01;
  @EXPORT      = qw();
  %EXPORT_TAGS = ( basic => [qw(myOpen myOpenRW fixPath createDirs myBasename min max round readConfigFile returnHumanReadableTime)], html => [qw(writeHTMLheader writeHTMLend)] );
  @EXPORT_OK   = qw(myOpen myOpenRW fixPath createDirs myBasename
                    min max round insertThousandSeparator readConfigFile returnHumanReadableTime
                    reverseComplement complement complementA translate IUPAC2baseA returnIUPACcode
                    readFastaIndex returnFaFHandIndex returnSequenceFromFasta
                    read_bed_variant_file readKgXref readKgXrefHR getFileHeader retFileType readGtpInfo readFilterInfo
                    returnChrSortIndex isCanonicalChr chrTranslate
                    returnGenomeSize returnOverlap returnGCfraction
                    writeHTMLheader writeHTMLend);
}

#our @EXPORT_OK;


# File input output and path realted functions

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

sub fixPath {
  my $path = shift;
  $path =~ s/\/{2,}/\//g;
  return $path;
}

sub createDirs {
  my $dir = shift;
  $dir = fixPath($dir);
  my $ret = `mkdir -p $dir` unless -d $dir;
  return $dir;
}

sub myBasename {
  my $fullname = shift;
  if ( $fullname =~ /\/([^\/]+)$/ ) { return $1; } else { return $fullname; }
}

# Simple math and stat functions

sub min {
  return (sort {$a <=> $b} @_)[0];
}

sub max {
  return (sort {$b <=> $a} @_)[0];
}

sub round {
  return sprintf("%0.*f", @_[1,0]);
}

sub insertThousandSeparator {
  my $number = shift;
  $number = round($number, $_[0]) if $_[0];
  $number =~ s/\d{1,3}(?=(\d{3})+(?!\d))/$&,/g;
  return $number;
}

# Useful stuff

sub returnHumanReadableTime {
  my $sec = shift;
  return sprintf("%02d-%02d:%02d:%02d", int($sec/(24*60*60)), ($sec/(60*60))%24, ($sec/60)%60, $sec%60);
}


# Read config file

sub readConfigFile {
  my $configHR = {};
  my $fh = myOpen(fixPath($FindBin::Bin . "/resequencing_pipeline.conf"));
  while ( my $line = <$fh> ) {
    chomp($line);
    $line =~ s/\s*\#.+$//;
    next unless $line =~ /\S+/;
    if ( $line =~ /(\w+)\\([\w\-]+)\s*=\s*\"([^\"]*)\"/ ) {
      my ($type, $key, $value) = (uc($1), $2, $3);
      while ( $value =~ /\[([^\]]+)\]/g ) { $value =~ s/\[([^\]]+)\]/$configHR->{$type}{$1}\// if defined($configHR->{$type}{$1}) };
      $value = fixPath($value) if $type =~ /^PATH/;
      $configHR->{$type}{$key} = $value;
    }
  }
  return $configHR;
}



# Sequence transformation

sub reverseComplement {
  my $seqIn = shift;
  my $pattern = "[" . join( "", keys( %{$myBioDataHR->{'complementHR'}} ) ) . "]";
  if ( ref($seqIn) ) {
    if ( defined( wantarray ) ) {
      my $retSeq = reverse($$seqIn);
      $retSeq =~ s/($pattern)/$myBioDataHR->{'complementHR'}{$1}/g;
      return $retSeq;
    } else {
      $$seqIn = reverse($$seqIn);
      $$seqIn =~ s/($pattern)/$myBioDataHR->{'complementHR'}{$1}/g;
    }
  } else {
    $seqIn = reverse($seqIn);
    $seqIn =~ s/($pattern)/$myBioDataHR->{'complementHR'}{$1}/g;
    return $seqIn;
  }
}

sub complement {
  my $seqIn = shift;
  my $pattern = "[" . join( "", keys( %{$myBioDataHR->{'complementHR'}} ) ) . "]";
  if ( ref($seqIn) ) {
    if ( defined( wantarray ) ) {
      my $retSeq = $$seqIn;
      $retSeq =~ s/($pattern)/$myBioDataHR->{'complementHR'}{$1}/g;
      return $retSeq;
    } else {
      $$seqIn =~ s/($pattern)/$myBioDataHR->{'complementHR'}{$1}/g;
    }
  } else {
    $seqIn =~ s/($pattern)/$myBioDataHR->{'complementHR'}{$1}/g;
    return $seqIn;
  }
}

sub complementA {
  my @seqIn = @_;
  for my $x (0 .. scalar(@seqIn) - 1) {
    die "\nBase $seqIn[$x] not found in complementArray" unless exists($myBioDataHR->{'complementHR'}{$seqIn[$x]});
    complement(\$seqIn[$x]);
  }
  return @seqIn;
}

sub translate {
  my $seq = shift;
  $seq =~ s/\s//g;
  my $aaSeq;
  while ($seq =~ /([\w]{3})/g) { $aaSeq .= $myBioDataHR->{'translationHR'}{uc($1)}; }
  return $aaSeq;
}

sub IUPAC2baseA {
  my $base = substr(shift, 0, 1);
  return @{$myBioDataHR->{'IUPAC2NucleotideAR_HR'}{uc($base)}};
}

sub returnIUPACcode {
  my $baseStr = join("", sort split(//, uc(shift)));
  $baseStr =~ s/(\w)\1+/$1/g;
  return $myBioDataHR->{'myIUPACCodes'}{$baseStr};
}


# Get seq from large fasta files

sub readFastaIndex {
  my $faiFN = shift() . ".fai";
#  indexFastaFile($faiFN) unless -f $faiFN;
  die "\nCannot find the .fai index file, please index fasta file using samtools\n" unless -f $faiFN;
  my $faiFH = myOpen($faiFN);
  my $idxHR = {};  
  while ( my $line = <$faiFH> ) {
    chomp($line);
    my @cols = split(/\t/, $line);
    $idxHR->{$cols[0]} = {'len' => $cols[1], 'start' => $cols[2], 'lineLen' => $cols[3], 'lineLenN' => $cols[4], 'Nlen' => $cols[4] - $cols[3],};
  }
  return $idxHR;
}

sub returnFaFHandIndex {
  my $faFN  = shift;
  my $faFH  = myOpen($faFN);
  my $idxHR = readFastaIndex($faFN);
  return {'faFH' => $faFH, 'idxHR' => $idxHR};
}


sub returnSequenceFromFasta {
  my $faIdxHR = shift;
  my $faFH    = $faIdxHR->{'faFH'};
  die "\nNo filehandle, please initiate using the returnFaFHandIndex() subroutine\n" unless ref($faFH) eq "GLOB";
  my $chr     = shift;
  die "\nSequence name ($chr) does not exist or hash has not been initiated ( returnFaFHandIndex() )\n" unless defined($faIdxHR->{'idxHR'}{$chr});
  my $start   = shift;
  my $len     = shift;
  die "\nStart plus len is longer than the sequence\n" if $start + $len > $faIdxHR->{'idxHR'}{$chr}{'len'};
  my $seq     = "";
  my $startInFile = $faIdxHR->{'idxHR'}{$chr}{'start'} + $start + ( int($start / $faIdxHR->{'idxHR'}{$chr}{'lineLen'}) * $faIdxHR->{'idxHR'}{$chr}{'Nlen'} );
  seek($faFH, $startInFile, 0);
  while ( length($seq) < $len ) {
    $seq .= <$faFH>;
    chomp($seq);
  }
  return substr($seq, 0, $len);
}


# Read and parse variants file (my bed-like format)

# Read bed file
sub read_bed_variant_file {
  my $dataHR = shift;
  my $bedFH = myOpen(shift);
  my $modifyCoords = shift;
  my $idx = 0;
  while ( my $line = <$bedFH> ) {
    $idx++;
    chomp($line);
    my @cols = split(/\t/, $line);
    my $varHR  = {'line' => $line, 'chr' => $cols[0], 'sta' => $cols[1], 'end' => $cols[2], 'descr' => $cols[3], 'qual' => $cols[4], };
    $cols[3] =~ s/ //g;
    if ( $cols[3] =~ /([^=]+)=>([^\(]+)\((\d+);([^\)]+),?\)/ ) {
      $varHR->{'info'} = {'ref' => $1, 'var' => $2, 'depth' => $3, 'cntStr' => $4};
      if ( $varHR->{'info'}{'var'} =~ /([\*\+\-\w]+)\/([\*\+\-\w]+)/ ) {
        push(@{$varHR->{'allelesAR'}}, ($1 eq $2 ? $1 : ($1, $2)));
        $varHR->{'info'}{'homHet'} = $1 eq $2 ? "hom" : "het";
        $varHR->{'info'}{'varType'} = "indel";
        $dataHR->{'varInfo'}{'cntVarType'}{'indel'}++;
      } elsif ( $varHR->{'info'}{'var'} =~ /^(\w)$/ ) {
        push(@{$varHR->{'allelesAR'}}, IUPAC2baseA($1));
        $varHR->{'info'}{'homHet'} = $varHR->{'info'}{'var'} =~ /^[ACGT]$/ ? "hom" : "het";
        $varHR->{'info'}{'varType'} = "snp";
        $dataHR->{'varInfo'}{'cntVarType'}{'snp'}++;
      } else {
        warn "ERROR: $line does not appear to be a valid variation, skipping\n";
        next;
      }
      foreach my $baseCnt ( split(/,/, $varHR->{'info'}{'cntStr'}) ) {
        if ( $baseCnt =~ /([\*\+\-\w]+):(\d+)/ ) { $varHR->{'observed'}{uc($1)} = $2; }
      }
    }
    my $id = "$varHR->{'chr'}:$varHR->{'sta'}-$varHR->{'end'}_$varHR->{'info'}{'var'}";
    if ( exists($dataHR->{'varHR'}{$id}) ) { warn "Duplicate variants in input file, skipping $varHR->{'line'} (Line $. in file)\n)"; next; }
    if ( $modifyCoords and $varHR->{'sta'} == $varHR->{'end'} ) { $varHR->{'end'} += 2; }
    $varHR->{'id'} = $id;
    $varHR->{'idx'} = $idx;
    $dataHR->{'varHR'}{$id} = $varHR;
    push(@{$dataHR->{'varChr'}{$varHR->{'chr'}}}, $varHR);
    for my $pos ( $varHR->{'sta'} .. $varHR->{'end'} - 1 ) { $dataHR->{'varChrPos'}{$varHR->{'chr'}}{$pos}{$id} = $varHR; }
    if ( $varHR->{'sta'} == $varHR->{'end'} ) { $dataHR->{'varChrPos'}{$varHR->{'chr'}}{$varHR->{'sta'}}{$id} = $varHR; }
  }
  close($bedFH);
}


# Subs related to annotation

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

sub readKgXrefHR {
  my $kgXrefFN = shift || "/home/magnus.bjursell/hg/knownGene_GRCh37_hg19/kgXref.txt";
  my $hdrHR = getFileHeader("ucsc_kgXref");

  my $kgXrefHR = {};
  my $kgXrefFH = myOpen($kgXrefFN);
  while ( my $line = <$kgXrefFH> ) {
    chomp($line);
    my @cols = split(/\t/, $line);
    my $key  = $cols[0];
    foreach my $colHdr ( sort { $hdrHR->{$a} <=> $hdrHR->{$b} } keys %{$hdrHR} ) {
      $kgXrefHR->{$key}{$colHdr} = shift @cols;
    }
  }
  close($kgXrefFH);
  return $kgXrefHR;
}


# Return file headers
sub getFileHeader {
  my $format = shift;
  $format = exists($myFileHeadersHR->{$format}) ? $format : retFileType($format);
  return exists($myFileHeadersHR->{$format}) ? $myFileHeadersHR->{$format} : die "\n$format is not (and cannot be translated into) a recognized file format\n\n";
}


# File type identification

sub retFileType {
  my $fn = shift;
  given(" " . $fn) {
    when ( /[ \/]kgXref\.[^\/]*$/ )         { return "ucsc_kgXref"; }
    when ( /[ \/]ensGtp\.[^\/]*$/ )         { return "ucsc_ensGtp"; }
    when ( /[ \/]knownToRefSeq\.[^\/]*$/ )  { return "ucsc_knownToRefSeq"; }
    when ( /[ \/]knownGene\.[^\/]*$/ )      { return "ucsc_knownGene"; }
    when ( /[ \/]knownCanonical\.[^\/]*$/ ) { return "ucsc_knownCanonical"; }
    when ( /[ \/]ensGene\.[^\/]*$/ )        { return "ucsc_ensGene"; }
    when ( /[ \/]refGene\.[^\/]*$/ )        { return "ucsc_refGene"; }
    when ( /[ \/]knownIsoforms\.[^\/]*$/ )  { return "ucsc_knownIsoforms"; }
    when ( /[ \/]snp131\.[^\/]*$/ )         { return "ucsc_snp131"; }
    default                                 { die "\nCannot determine type for $fn\n\n"; }
  }
}


# Read Gtp files
sub readGtpInfo {
  my $fn  = shift;
  my $transCol = shift;
  my $geneCol  = shift;

  my $hdrHR = getFileHeader($fn);
  die "\nUnknown columns for transcript ($transCol) or gene ($geneCol) in file $fn\nAvailable columns: " . join("; ", sort { $hdrHR->{$a} <=> $hdrHR->{$b} } keys %{$hdrHR}) . "\n\n" unless exists($hdrHR->{$transCol}) and exists($hdrHR->{$geneCol});

  my $gtpHR = {};
  my $fh = myOpen($fn);
  while ( my $line = <$fh> ) {
    chomp($line);
    my @cols = split(/\t/, $line);
    next if length($cols[$hdrHR->{$geneCol}]) < 1;
    $gtpHR->{$cols[$hdrHR->{$transCol}]} = $cols[$hdrHR->{$geneCol}];
  }
  close($fh);
  return $gtpHR;
}


# Read Filter files
sub readFilterInfo {
  my $fn  = shift;
  my $col = shift;

  my $hdrHR = getFileHeader($fn);
  die "\nUnknown columns for filterring ($col) in file $fn\nAvailable columns: " . join("; ", sort { $hdrHR->{$a} <=> $hdrHR->{$b} } keys %{$hdrHR}) . "\n\n" unless exists($hdrHR->{$col});

  my $filterHR = {};
  my $fh = myOpen($fn);
  while ( my $line = <$fh> ) {
    chomp($line);
    my @cols = split(/\t/, $line);
    $filterHR->{$cols[$hdrHR->{$col}]} = 1;
#print STDERR "$fn; $col; $hdrHR->{$col}; $line; $cols[$hdrHR->{$col}]\n"; sleep 1;
  }
  close($fh);
  return $filterHR;
}



# Sortorder for chrs

sub returnChrSortIndex {
  my $chr = shift;
  return exists($myChrInformationHR->{'myChrSortOrder'}{$chr}) ? $myChrInformationHR->{'myChrSortOrder'}{$chr} : 999;
}


sub isCanonicalChr {
  my $chr = shift;
  return exists($myChrInformationHR->{'myCanonicalChr'}{$chr}) ? $myChrInformationHR->{'myCanonicalChr'}{$chr} : undef;
}

sub chrTranslate {
  my $chr = shift;
  return exists($myChrInformationHR->{'myChrTranslate'}{$chr}) ? $myChrInformationHR->{'myChrTranslate'}{$chr} : undef;
}

# Sequence analysis

sub returnGenomeSize {
  my $refFN = shift;
  my $genomeSizeHR = {};
  my $mySizeFileSuffix = "mySizeFile";

  if ( -f "$refFN.$mySizeFileSuffix" ) {
    $genomeSizeHR = readMySizeFile("$refFN.$mySizeFileSuffix");
  } else {
    $genomeSizeHR = calculateGenomeSize($refFN);
    createFastaSizeFile("$refFN.$mySizeFileSuffix", $genomeSizeHR);
  }
  return $genomeSizeHR;
}


sub readMySizeFile {
  my $sizeFN = shift;
  my $genomeSizeHR = {};
  my $sizeFH = myOpen($sizeFN);
  while ( my $line = <$sizeFH> ) {
    chomp($line);
    if ( $line =~ /Total_size\t(\d+)/ ) {
      $genomeSizeHR->{'genome'} = $1;
    } else {
      my @cols = split(/\t/, $line);
      $genomeSizeHR->{'chr'}{$cols[0]} = $cols[1];
    }
  }
  close($sizeFH);
  return $genomeSizeHR;
}


sub createFastaSizeFile {
  my $sizeFN       = shift;
  my $genomeSizeHR = shift;

  warn "File $sizeFN already exists.\n" and return undef if -f $sizeFN;
  my $sizeFH = myOpenRW($sizeFN);

  printf $sizeFH ("Total_size\t%d\n", $genomeSizeHR->{'genome'});
  foreach my $chr ( @{$genomeSizeHR->{'chrSortOrder'}} ) {
    printf $sizeFH ("%s\t%d\n", $chr, $genomeSizeHR->{'chr'}{$chr});
  }
  close($sizeFH);
}


sub calculateGenomeSize {
  my $refFN = shift;
  my $refFH = myOpen($refFN);
  my $genomeSizeHR = {};
  my $currentChr = "";
  while ( my $line = <$refFH> ) {
    chomp($line);
    if ( $line =~ /^>(\S+)/ ) {
      $currentChr = $1;
      push(@{$genomeSizeHR->{'chrSortOrder'}}, $currentChr);
    } else {
      my $noBases = $line =~ tr/ACGTacgt//;
      $genomeSizeHR->{'genome'} += $noBases;
      $genomeSizeHR->{'chr'}{$currentChr} += $noBases;
    }
  }
  close($refFH);
  return $genomeSizeHR;
}



sub returnOverlap {
  my $coordAR = shift;
  my $overlap;
  if ($coordAR->[0] > $coordAR->[1]) {my $tmp = $coordAR->[0]; $coordAR->[0] = $coordAR->[1]; $coordAR->[1] = $tmp;}
  if ($coordAR->[2] > $coordAR->[3]) {my $tmp = $coordAR->[2]; $coordAR->[2] = $coordAR->[3]; $coordAR->[3] = $tmp;}
  if ( $coordAR->[2] >= $coordAR->[0] and $coordAR->[2] <= $coordAR->[1] ) {
    $overlap = min($coordAR->[1], $coordAR->[3]) - $coordAR->[2] + 1;
  } elsif ( $coordAR->[3] >= $coordAR->[0] and $coordAR->[3] <= $coordAR->[1] ) {
    $overlap = $coordAR->[3] - max($coordAR->[0], $coordAR->[2]) + 1;
  } elsif ( $coordAR->[2] < $coordAR->[0] and $coordAR->[3] > $coordAR->[1] ) {
    $overlap = min($coordAR->[1] - $coordAR->[0] + 1, $coordAR->[3] - $coordAR->[2] + 1);
  } else {
    $overlap = 0;
  }
  return $overlap;
}

sub returnGCfraction {
  my $seqIn = shift;
  my $noGC;
  return undef if length($seqIn) < 1;
  if ( ref($seqIn) =~ /SCALAR/ ) {
    $noGC++ while $$seqIn =~ /[CGcg]/g;
    return ($noGC / length($$seqIn) );
  } elsif ( ref($seqIn) =~ /HASH/ ) {
    $noGC++ while $seqIn->{'seq'} =~ /[CGcg]/g;
    return ($noGC / length($seqIn->{'seq'}) );
  } else {
    $noGC++ while $seqIn =~ /[CGcg]/g;
    return ($noGC / length($seqIn) );
  }
}


# html related funtions


sub writeHTMLheader {
  my $fh    = shift;
  my $title = shift;

#print $fh "Content-type:text/html\r\n\r\n";
print $fh <<END;

<html>
  <head>
    <meta http-equiv="Content-Type" content="text/html; charset=iso-8859-1">
    <title>$title</title>
  </head>
  <body>
END
}

sub writeHTMLend {
  my $fh    = shift;

print $fh <<END;

  </body>
</html>
END
}



1;




