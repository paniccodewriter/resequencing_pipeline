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
use mySharedFunctions qw(:basic myBasename fixPath isCanonicalChr chrTranslate returnGenomeSize returnOverlap);


my $dataHR = {};
my $infoHR = {};

#$infoHR->{'configHR'} = readConfigFile();
$infoHR->{'optHR'}    = {};


# Read input and setup global variables
GetOptions ($infoHR->{'optHR'}, 'ref=s', 'gene=s', 'bam=s', 'outdir=s', 'outPrefix=s', 'outCoverage=s', 'minCoverage=s', 'minConcQual=s', 'minMappingQual=s', 'help', 'verbose') or pod2usage(2);
pod2usage(1) if $infoHR->{'optHR'}{'h'} or $infoHR->{'optHR'}{'help'};
loadDefaults($infoHR);

pod2usage("Cannot find the reference file: $infoHR->{'optHR'}{'ref'}.") unless -f $infoHR->{'optHR'}{'ref'};
pod2usage("Cannot find the regions gene file: $infoHR->{'optHR'}{'gene'}.") unless -f $infoHR->{'optHR'}{'gene'};
pod2usage("Cannot find the input bam file: $infoHR->{'optHR'}{'bam'}.") unless -f $infoHR->{'optHR'}{'bam'};
pod2usage("Cannot find output directory: $infoHR->{'optHR'}{'outdir'}.") unless -d $infoHR->{'optHR'}{'outdir'};

$infoHR->{'optHR'}{'outPrefix'} = myBasename($infoHR->{'optHR'}{'bam'}) unless $infoHR->{'optHR'}{'outPrefix'};

getGenomeSize($infoHR);

{ my $idx = 1; $infoHR->{'chrSortOrder'} = { map { 'chr' . $_ => $idx++ } (1 .. 22, 'X', 'Y', 'M') }; }

{



  for my $eNo (0 .. scalar( @{$transHR->{'eStaA'}} ) - 1 ) {

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





exit;


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
   --ref             Reference sequence fasta file [required]
   --gene            UCSC gene file describing the genes [required]
   --bam             Bam file containing the alignments [required]
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


