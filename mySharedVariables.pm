package mySharedVariables;

# (c) 2008- Magnus Bjursell

use strict;
use warnings;

BEGIN {
  use Exporter ();
  our ($VERSION, @ISA, @EXPORT, @EXPORT_OK, %EXPORT_TAGS);
  @ISA = qw(Exporter);
  $VERSION     = 0.01;
  @EXPORT      = qw();
  %EXPORT_TAGS = (basic => [qw($myFileLocHR $myBioDataHR $myChrInformationHR $myFileHeadersHR)]);
  @EXPORT_OK   = qw($myFileLocHR $myBioDataHR $myChrInformationHR $myFileHeadersHR);
}

our @EXPORT_OK;

# exported package globals go here
our ($myFileLocHR, $myBioDataHR, $myChrInformationHR, $myFileHeadersHR);

# file and directory information

$myFileLocHR = {
    'tmpDir'               => 'tmp/',
    'blastDir'             => '/home/magnus.bjursell/tools/blast/blast-2.2.17/bin/',
    'bl2seq'               => 'bl2seq',
    'blastall'             => 'blastall',
    'myBlastPgmAR'         => [['blastn', 'blastx'], ['tblastn', 'blastp']],
    'myBlastPgmHR'         => {'nt' => {'nt' => 'blastn', 'aa' => 'blastx'}, 'aa' => {'nt' => 'tblastn', 'aa' => 'blastp'}, 'NT' => {'NT' => 'blastn', 'AA' => 'blastx'}, 'AA' => {'NT' => 'tblastn', 'AA' => 'blastp'}},
    'UCSCseqDir'           => '/home/magnus.bjursell/bioDB/UCSCgenomeBrowser/chromFa/',
    'UCSCallSeq'           => 'allChrom.fa',
    'UCSCallSeqBlastDBnt'  => '/home/magnus.bjursell/bioDB/UCSCgenomeBrowser/chromFa/blastDB/UCSC_chromDB',
    'ENScDNABlastDBnt'     => '/home/magnus.bjursell/bioDB/Ensembl/blastDB/ENS_NCBI36_49_cDNA',
    'ENSpeptBlastDBaa'     => '/home/magnus.bjursell/bioDB/Ensembl/blastDB/ENS_NCBI36_49_protein',
    'TCDBBlastDBaa'        => '/home/magnus.bjursell/bioDB/TCDB/blastDB/TCDB_protein',
    'TransportDBBlastDBnt' => '/home/magnus.bjursell/bioDB/TransportDB/blastDB/TransportDB_nucleotide',
    'TransportDBBlastDBaa' => '/home/magnus.bjursell/bioDB/TransportDB/blastDB/TransportDB_protein',
};

# assign values to globals

$myBioDataHR = {
    'replaceBaseHR' => {
        'sameCase'  => {'A' => ['C', 'G', 'T'], 'C' => ['A', 'G', 'T'], 'G' => ['A', 'C', 'T'], 'T' => ['A', 'C', 'G'], 'ALL' => ['A', 'C', 'G', 'T'],
                        'a' => ['c', 'g', 't'], 'c' => ['a', 'g', 't'], 'g' => ['a', 'c', 't'], 't' => ['a', 'c', 'g'], 'all' => ['a', 'c', 'g', 't']},
        'shiftCase' => {'A' => ['c', 'g', 't'], 'C' => ['a', 'g', 't'], 'G' => ['a', 'c', 't'], 'T' => ['a', 'c', 'g'], 'ALL' => ['a', 'c', 'g', 't'],
                        'a' => ['C', 'G', 'T'], 'c' => ['A', 'G', 'T'], 'g' => ['A', 'C', 'T'], 't' => ['A', 'C', 'G'], 'all' => ['A', 'C', 'G', 'T']},
    },

    'alphabetsHR' => {
        'dnaAR'          => ['A', 'C', 'G', 'T'],
        'dnaMatch'       => "[ATCGatcg]",
        'dnaXtMatch'     => "[ATCGXNatcgxn]",
        'dnaNotMatch'    => "[^ATCGatcg]",
        'dnaNotXtMatch'  => "[^ATCGXNatcgxn]",
        'dnaXtStarMatch' => "[ATCGXNatcgxn\*]",
    },

    'alphabetAR' => ['A', 'C', 'G', 'T'],

    'complementHR' => {
        'A' => 'T',
        'C' => 'G',
        'G' => 'C',
        'T' => 'A',
        'a' => 'T',
        'c' => 'G',
        'g' => 'C',
        't' => 'A',
    },

    'enzymeRecognisionSeqHR' => {
        'FauI'   => 'CCCGC',
        'EarI'   => 'CTCTTC',
        'AluI'   => 'AGCT',
        'HaeIII' => 'GGCC',
        'SmaI'   => 'CCCGGG',
        'StuI'   => 'AGGCCT',
    },

    'translationHR' => {'AAA' => 'K','AAC' => 'N','AAG' => 'K','AAT' => 'N',
                        'ACA' => 'T','ACC' => 'T','ACG' => 'T','ACT' => 'T',
                        'AGA' => 'R','AGC' => 'S','AGG' => 'R','AGT' => 'S',
                        'ATA' => 'I','ATC' => 'I','ATG' => 'M','ATT' => 'I',
                        'CAA' => 'Q','CAC' => 'H','CAG' => 'Q','CAT' => 'H',
                        'CCA' => 'P','CCC' => 'P','CCG' => 'P','CCT' => 'P',
                        'CGA' => 'R','CGC' => 'R','CGG' => 'R','CGT' => 'R',
                        'CTA' => 'L','CTC' => 'L','CTG' => 'L','CTT' => 'L',
                        'GAA' => 'E','GAC' => 'D','GAG' => 'E','GAT' => 'D',
                        'GCA' => 'A','GCC' => 'A','GCG' => 'A','GCT' => 'A',
                        'GGA' => 'G','GGC' => 'G','GGG' => 'G','GGT' => 'G',
                        'GTA' => 'V','GTC' => 'V','GTG' => 'V','GTT' => 'V',
                        'TAA' => 'X','TAC' => 'Y','TAG' => 'X','TAT' => 'Y',
                        'TCA' => 'S','TCC' => 'S','TCG' => 'S','TCT' => 'S',
                        'TGA' => 'X','TGC' => 'C','TGG' => 'W','TGT' => 'C',
                        'TTA' => 'L','TTC' => 'F','TTG' => 'L','TTT' => 'F',
    },

    'IUPAC2NucleotideAR_HR' => {
        'A' => ['A'],
        'G' => ['G'],
        'C' => ['C'],
        'T' => ['T'],
        'U' => ['U'],
        'R' => ['A', 'G'],
        'Y' => ['C', 'T'],
        'N' => ['A', 'C', 'G', 'T'],
        'W' => ['A', 'T'],
        'S' => ['G', 'C'],
        'M' => ['A', 'C'],
        'K' => ['G', 'T'],
        'B' => ['G', 'C', 'T'],
        'H' => ['A', 'C', 'T'],
        'D' => ['A', 'G', 'T'],
        'V' => ['A', 'G', 'C'],
    },

    'nucleotideStr2IUPAC_HR' => {
      'A' => 'A',
      'G' => 'G',
      'C' => 'C',
      'T' => 'T',
      'U' => 'U',
      'AC' => 'M',
      'AG' => 'R',
      'AT' => 'W',
      'CG' => 'S',
      'CT' => 'Y',
      'GT' => 'K',
      'ACG' => 'V',
      'ACT' => 'H',
      'AGT' => 'D',
      'CGT' => 'B',
      'ACGT' => 'N',
    },

};


$myChrInformationHR = {
  
  'myChrSortOrder' => {

    'chr1' => 0, 'chr1_gl000191_random' => 1, 'chr1_gl000192_random' => 2, 'chr2' => 3, 'chr3' => 4, 'chr4' => 5, 'chr4_ctg9_hap1' => 6, 'chr4_gl000193_random' => 7, 
    'chr4_gl000194_random' => 8, 'chr5' => 9, 'chr6' => 10, 'chr6_apd_hap1' => 11, 'chr6_cox_hap2' => 12, 'chr6_dbb_hap3' => 13, 'chr6_mann_hap4' => 14, 'chr6_mcf_hap5' => 15, 
    'chr6_qbl_hap6' => 16, 'chr6_ssto_hap7' => 17, 'chr7' => 18, 'chr7_gl000195_random' => 19, 'chr8' => 20, 'chr8_gl000196_random' => 21, 'chr8_gl000197_random' => 22, 
    'chr9' => 23, 'chr9_gl000198_random' => 24, 'chr9_gl000199_random' => 25, 'chr9_gl000200_random' => 26, 'chr9_gl000201_random' => 27, 'chr10' => 28, 'chr11' => 29, 
    'chr11_gl000202_random' => 30, 'chr12' => 31, 'chr13' => 32, 'chr14' => 33, 'chr15' => 34, 'chr16' => 35, 'chr17' => 36, 'chr17_ctg5_hap1' => 37, 'chr17_gl000203_random' => 38, 
    'chr17_gl000204_random' => 39, 'chr17_gl000205_random' => 40, 'chr17_gl000206_random' => 41, 'chr18' => 42, 'chr18_gl000207_random' => 43, 'chr19' => 44, 
    'chr19_gl000208_random' => 45, 'chr19_gl000209_random' => 46, 'chr20' => 47, 'chr21' => 48, 'chr21_gl000210_random' => 49, 'chr22' => 50, 'chrX' => 51, 
    'chrY' => 52, 'chrM' => 53, 'chrUn_gl000211' => 54, 'chrUn_gl000212' => 55, 'chrUn_gl000213' => 56, 'chrUn_gl000214' => 57, 'chrUn_gl000215' => 58, 'chrUn_gl000216' => 59, 
    'chrUn_gl000217' => 60, 'chrUn_gl000218' => 61, 'chrUn_gl000219' => 62, 'chrUn_gl000220' => 63, 'chrUn_gl000221' => 64, 'chrUn_gl000222' => 65, 'chrUn_gl000223' => 66, 
    'chrUn_gl000224' => 67, 'chrUn_gl000225' => 68, 'chrUn_gl000226' => 69, 'chrUn_gl000227' => 70, 'chrUn_gl000228' => 71, 'chrUn_gl000229' => 72, 'chrUn_gl000230' => 73, 
    'chrUn_gl000231' => 74, 'chrUn_gl000232' => 75, 'chrUn_gl000233' => 76, 'chrUn_gl000234' => 77, 'chrUn_gl000235' => 78, 'chrUn_gl000236' => 79, 'chrUn_gl000237' => 80, 
    'chrUn_gl000238' => 81, 'chrUn_gl000239' => 82, 'chrUn_gl000240' => 83, 'chrUn_gl000241' => 84, 'chrUn_gl000242' => 85, 'chrUn_gl000243' => 86, 'chrUn_gl000244' => 87, 
    'chrUn_gl000245' => 88, 'chrUn_gl000246' => 89, 'chrUn_gl000247' => 90, 'chrUn_gl000248' => 91, 'chrUn_gl000249' => 92, 
  },

  'myCanonicalChr' => {
    'chr1'  => 1,
    'chr2'  => 2,
    'chr3'  => 3,
    'chr4'  => 4,
    'chr5'  => 5,
    'chr6'  => 6,
    'chr7'  => 7,
    'chr8'  => 8,
    'chr9'  => 9,
    'chr10' => 10,
    'chr11' => 11,
    'chr12' => 12,
    'chr13' => 13,
    'chr14' => 14,
    'chr15' => 15,
    'chr16' => 16,
    'chr17' => 17,
    'chr18' => 18,
    'chr19' => 19,
    'chr20' => 20,
    'chr21' => 21,
    'chr22' => 22,
    'chrX'  => 23,
    'chrY'  => 24,
    'chrM'  => 25,
  },

  'myChrTranslate' => {
    '1'  => 'chr1',
    '2'  => 'chr2',
    '3'  => 'chr3',
    '4'  => 'chr4',
    '5'  => 'chr5',
    '6'  => 'chr6',
    '7'  => 'chr7',
    '8'  => 'chr8',
    '9'  => 'chr9',
    '10' => 'chr10',
    '11' => 'chr11',
    '12' => 'chr12',
    '13' => 'chr13',
    '14' => 'chr14',
    '15' => 'chr15',
    '16' => 'chr16',
    '17' => 'chr17',
    '18' => 'chr18',
    '19' => 'chr19',
    '20' => 'chr20',
    '21' => 'chr21',
    '22' => 'chr22',
    'X'  => 'chrX',
    'x'  => 'chrX',
    'Y'  => 'chrY',
    'y'  => 'chrY',
    'MT' => 'chrM',
    'mt' => 'chrM',
    'M'  => 'chrM',
    'm'  => 'chrM',
  },

};



$myFileHeadersHR = {
  'ucsc_knownGene'      => {'id' => 0, 'chr' => 1, 'str' => 2, 'sta' => 3, 'end' => 4, 'cSta' => 5, 'cEnd' => 6, 'eCnt' => 7, 'eSta' => 8, 'eEnd' => 9, 'prId' => 10, 'alId' => 11,},
  'ucsc_knownCanonical' => {'chr' => 0, 'sta' => 1, 'end' => 2, 'clusterId' => 3, 'id' => 4, 'protein' => 5,},
  'ucsc_ensGene'        => {'bin' => 0, 'id' => 1, 'chr' => 2, 'str' => 3, 'sta' => 4, 'end' => 5, 'cSta' => 6, 'cEnd' => 7, 'eCnt' => 8, 'eSta' => 9, 'eEnd' => 10, 'scr' => 11, 'id2' => 12, 'cStaSt' => 13, 'cEndSt' => 14, 'eFrm' => 15,},
  'ucsc_refGene'        => {'bin' => 0, 'id' => 1, 'chr' => 2, 'str' => 3, 'sta' => 4, 'end' => 5, 'cSta' => 6, 'cEnd' => 7, 'eCnt' => 8, 'eSta' => 9, 'eEnd' => 10, 'scr' => 11, 'id2' => 12, 'cStaSt' => 13, 'cEndSt' => 14, 'eFrm' => 15,},
  'ucsc_kgXref'         => {'kgID' => 0, 'mRNA' => 1, 'spID' => 2, 'spDisplayID' => 3, 'geneSymbol' => 4, 'refseq' => 5, 'protAcc' => 6, 'description' => 7,},
  'ucsc_ensGtp'         => {'gene' => 0, 'transcript' => 1, 'protein' => 2,},
  'ucsc_knownToRefSeq'  => {'name' => 0, 'value' => 1},
  'ucsc_snp131'         => {'bin' => 0, 'chr' => 1, 'sta' => 2, 'end' => 3, 'name' => 4, 'scr' => 5, 'str' => 6, 'refNCBI' => 7, 'refUCSC' => 8, 'obs' => 9, 'molType' => 10, 'class' => 11, 'valid' => 12, 'avHet' => 13, 'avHetSE' => 14, 'func' => 15, 'locType' => 16, 'weight' => 17,},
  'ucsc_knownIsoforms'  => {'cluster' => 0, 'id' => 1},
};


1;



