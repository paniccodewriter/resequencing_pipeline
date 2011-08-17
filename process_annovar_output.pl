#!/usr/bin/perl
# (c) 2011 Magnus Bjursell

# standard module usage
use strict;
use Getopt::Long;
use Pod::Usage;
use Cwd 'abs_path';
use FindBin;
use Data::Dumper;
use feature ":5.10";

# custom module usage
use lib $FindBin::Bin; #use lib '/home/magnus.bjursell/script/modules';
use mySharedFunctions qw(:basic returnFaFHandIndex returnSequenceFromFasta myBasename IUPAC2baseA);

my $dataHR = {};
my $infoHR = {};
#my $geneHR = {};

$infoHR->{'optHR'} = {};

# Read input and setup global variables
GetOptions ($infoHR->{'optHR'}, 'sampleid=s', 'instub=s', 'annovardb=s', 'outdb=s', 'outdir=s', 'outPrefix=s', 'help', 'verbose') or pod2usage(2);
pod2usage(1) if $infoHR->{'optHR'}{'h'} or $infoHR->{'optHR'}{'help'};
loadDefaults($infoHR);

pod2usage("Sample ID is required.") unless defined($infoHR->{'optHR'}{'sampleid'});
pod2usage("Cannot find annovar database directory: $infoHR->{'optHR'}{'annovardb'}.") unless -d $infoHR->{'optHR'}{'annovardb'};
#pod2usage("Cannot find the database file: $infoHR->{'optHR'}{'outdb'}.") unless -f $infoHR->{'optHR'}{'outdb'};
#pod2usage("Cannot find output directory: $infoHR->{'optHR'}{'outdir'}.") unless -d $infoHR->{'optHR'}{'outdir'};

{
  readAnnovarData($dataHR, $infoHR);
#  getFilesInDir($infoHR->{'optHR'}{'instub'});
}


#####################
# # #  S U B S  # # #
#####################

sub readAnnovarData {
  my $dataHR = shift;
  my $infoHR = shift;

  my $inFileHR = getFilesInDir($infoHR->{'optHR'}{'instub'});

  foreach my $key ( keys %{$inFileHR} ) {
    print "$key \t-> \t$inFileHR->{$key}\n";
  }

}

sub getFilesInDir {
  my $stub = shift;
  my $inFileHR = {};

  pod2usage("Cannot find input files based on stub: $infoHR->{'optHR'}{'instub'}.") unless -f "$infoHR->{'optHR'}{'instub'}.variant_function" and -f "$infoHR->{'optHR'}{'instub'}.exonic_variant_function";

  while ( my $fn = <$stub*> ) {
    next if $fn =~ /\.log$/;
    if ( $fn =~ /$stub\.(.+)$/ ) { $inFileHR->{$1} = $fn; }
  }
  return $inFileHR;
}


#sub returnAnnovarFileTypes {
#  return {
#    '' => "",
#  };

#Filename: 143-10F/testAnnotation/143-10F.merged.recal.realn.reSrt.exonic_variant_function
#Filename: 143-10F/testAnnotation/143-10F.merged.recal.realn.reSrt.hg19_ALL.sites.2010_11_dropped
#Filename: 143-10F/testAnnotation/143-10F.merged.recal.realn.reSrt.hg19_ALL.sites.2010_11_filtered
#Filename: 143-10F/testAnnotation/143-10F.merged.recal.realn.reSrt.hg19_avsift_dropped
#Filename: 143-10F/testAnnotation/143-10F.merged.recal.realn.reSrt.hg19_avsift_filtered
#Filename: 143-10F/testAnnotation/143-10F.merged.recal.realn.reSrt.hg19_ljb_pp2_dropped
#Filename: 143-10F/testAnnotation/143-10F.merged.recal.realn.reSrt.hg19_ljb_pp2_filtered
#Filename: 143-10F/testAnnotation/143-10F.merged.recal.realn.reSrt.hg19_snp132_dropped
#Filename: 143-10F/testAnnotation/143-10F.merged.recal.realn.reSrt.hg19_snp132_filtered
#Filename: 143-10F/testAnnotation/143-10F.merged.recal.realn.reSrt.log
#Filename: 143-10F/testAnnotation/143-10F.merged.recal.realn.reSrt.variant_function
#}

sub loadDefaults {
  my $infoHR = shift;

  $infoHR->{'optHR'}{'instub'} =~ s/\.+$//;
#  $infoHR->{'optHR'}{'minCoverage'}    = 8     unless $infoHR->{'optHR'}{'minCoverage'};
}


#__END__


=head1 NAME

process_annovar_output.pl - Analysis of sequencing coverage across specified regions

=head1 SYNOPSIS

process_annovar_output.pl [options]

  Options:
   --help            Brief help message (not yet implemented)
   --verbose         Write some additional output (not yet implemented)
   --sampleid        Sample id [required]
   --instub          Input filename base (including directory) for the annovar output files [required]
   --annovardb       Directory containing annovar database files [required]
   --outdb           Output directory [required]
   --outdir          Output directory
   --outPrefix       Output file prefix, bam file name will be used if not supplied




=cut

#=head1 OPTIONS

#=over 8

#=item B<-help>

#Print a brief help message and exits.

#=back

#=head1 DESCRIPTION



