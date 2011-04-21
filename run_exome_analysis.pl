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
use mySharedFunctions qw(:basic);


my $dataHR = {};
my $infoHR = {};

$infoHR->{'configHR'} = readConfigFile();
$infoHR->{'optHR'}    = {};


# Read input and setup global variables
GetOptions ($infoHR->{'optHR'}, 'fq1=s', 'fq2=s', 'fqS=s', 'outdir=s', 'sample=s', 'run=s', 'shm', 'help', 'verbose') or pod2usage(2);
pod2usage(1) if $infoHR->{'optHR'}{'h'} or $infoHR->{'optHR'}{'help'};

#print Dumper($infoHR);

checkInParameters($infoHR);
#checkExternalSoftware();


setupDirStructure($infoHR);
writeSbatchScripts($infoHR);




exit;



############################################################################
###############################  S  U  B  S  ###############################
############################################################################


### Check for software ###

sub checkExternalSoftware {
  checkSamTools();
}

sub checkMosaik {
  # Coming soon
}

sub checkGATK {
  # Coming soon
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


### Setup parameters and directories

sub checkInParameters {
  my $infoHR  = shift;

  unless ( -f $infoHR->{'optHR'}{'fq1'} or -f $infoHR->{'optHR'}{'fq2'} or -f $infoHR->{'optHR'}{'fqS'} ) { pod2usage("Error: At lease one sequence data file must be submitted.") }
  if ( -f $infoHR->{'optHR'}{'fq1'} xor -f $infoHR->{'optHR'}{'fq2'} ) { pod2usage("Error: Both fq1 and fq2 must be present and match for paired end data.") }
  if ( -f $infoHR->{'optHR'}{'fqS'} ) { pod2usage("Error: Single end read processing not implremented yet.") }

  unless ( -d $infoHR->{'optHR'}{'outdir'} ) { pod2usage("Error: An existing output directory must be provided.") }
  unless ( $infoHR->{'optHR'}{'sample'} ) { pod2usage("Error: Please provide a sample name for the output name (output filenames will be based on the prefix).") }
  unless ( $infoHR->{'optHR'}{'run'} ) { pod2usage("Error: Please provide a run name for the output name (output filenames will be based on the prefix).") }

  processInParameters($infoHR);
}

sub processInParameters {
  my $infoHR  = shift;

  $infoHR->{'optHR'}{'prefix'} = "$infoHR->{'optHR'}{'sample'}.$infoHR->{'optHR'}{'run'}";
  my @okSuffix = qw(fastq fq); my $matchStr = join("|", @okSuffix);
  for my $p ( 'fq1', 'fq2', 'fqS', 'outdir' ) {
    next unless -e $infoHR->{'optHR'}{$p};
    $infoHR->{'optHR'}{$p} = abs_path($infoHR->{'optHR'}{$p});
    if ( $p =~ /^fq\w$/ ) {
      if ( $infoHR->{'optHR'}{$p} =~ /(([^\/]+)\.($matchStr))$/ ) {
        $infoHR->{'optHR'}{$p . "n"} = $1;
        $infoHR->{'optHR'}{$p . "b"} = $2;
      } else {
        pod2usage("Error: Input sequence files must be fastq and end with " . join(" or ", @okSuffix) . ".")
      }
    }
  }
}

sub setupDirStructure {
  my $infoHR  = shift;

  $infoHR->{'analysis_dirs'}{'sbatch'} = createDirs("$infoHR->{'optHR'}{'outdir'}/sbatch");
  $infoHR->{'analysis_dirs'}{'sbatch-info'} = createDirs("$infoHR->{'optHR'}{'outdir'}/sbatch/info");
  $infoHR->{'analysis_dirs'}{'sbatch-cur'} = createDirs("$infoHR->{'optHR'}{'outdir'}/sbatch/current");
  for my $basedir ('fastq', 'fastq_filter', 'mosaik', 'GATK', 'variants', 'annotation') {
    $infoHR->{'analysis_dirs'}{$basedir} = createDirs("$infoHR->{'optHR'}{'outdir'}/$basedir");
    for my $subdir ('info', 'results') {
      $infoHR->{'analysis_dirs'}{$basedir . "-$subdir"} = createDirs("$infoHR->{'optHR'}{'outdir'}/$basedir/$subdir");
    }
  }
  $infoHR->{'analysis_dirs'}{'GATK-intm'} = createDirs("$infoHR->{'optHR'}{'outdir'}/GATK/intermediary");
  $infoHR->{'analysis_dirs'}{'mosaik-mosaikDS_PE'} = createDirs("$infoHR->{'optHR'}{'outdir'}/mosaik/MosaikDupSnoop_PE");
  $infoHR->{'analysis_dirs'}{'mosaik-mosaikDS_SE'} = createDirs("$infoHR->{'optHR'}{'outdir'}/mosaik/MosaikDupSnoop_SE");
#print Dumper($infoHR->{'analysis_dirs'});
}

sub myBasename {
  my $fullname = shift;
  if ( $fullname =~ /\/([^\/]+)$/ ) { return $1; } else { return $fullname; }
}

### Sbatch functions

sub writeCurrent {
  my $infoHR  = shift;
  my $cmdAR  = shift;
  my $current = shift;

  $current =~ s/\"/\\\"/g;
  push(@{$cmdAR}, {'cmd' => qq{echo "$current" > $infoHR->{'optHR'}{'sbatch-cur-file'}}},);
}

sub writeSbatchScripts {
  my $infoHR  = shift;

  $infoHR->{'optHR'}{'sbatch-cur-SE'} = "$infoHR->{'analysis_dirs'}{'sbatch-cur'}/$infoHR->{'optHR'}{'prefix'}.SE.running";
  $infoHR->{'optHR'}{'sbatch-cur-PE'} = "$infoHR->{'analysis_dirs'}{'sbatch-cur'}/$infoHR->{'optHR'}{'prefix'}.PE.running";

  { # fastQC x 2, filter and MosaikBuild
    if ( -f $infoHR->{'optHR'}{'fq1'} and -f $infoHR->{'optHR'}{'fq2'} ) {
      my $ret = `echo "1" > $infoHR->{'optHR'}{'sbatch-cur-SE'}`;
      my $ret = `echo "1" > $infoHR->{'optHR'}{'sbatch-cur-PE'}`;
      my $scriptPgmNameAR = ['fastQC_filter', 'mosaik_PE', 'mosaik_SE', 'GATK', 'variants'];
      my $scriptNameAR = []; for my $idx ( 0 .. scalar(@{$scriptPgmNameAR}) - 1) { $scriptNameAR->[$idx] = $infoHR->{'optHR'}{'prefix'} . "." . $idx . "." . $scriptPgmNameAR->[$idx] . ".sbatch"; }


      { # FastQC and filter
        my $cpHR = {};
        my $cIdx = 0;
        my $cmdAR = [];

        my $paramAR = returnParamAR($infoHR, {'scriptNamePgm' => $scriptPgmNameAR->[$cIdx], 'C' => "thin", 't' => "50:00:00"});
        returnInitAndModuleLoad($infoHR, $cmdAR, 'Bioinfo-tools', 'FastQC', 'Mosaik');


        push(@{$cmdAR},
#          {'cmd' => "echo \"running\" > $infoHR->{'optHR'}{'sbatch-cur-prefix'}", 'xn' => 2},
          {'cmd' => "# FastQC"},
          {'cmd' => "fastqc -o $infoHR->{'analysis_dirs'}{'fastq-results'}/ $infoHR->{'optHR'}{'fq1'}", 'echo' => 1, 'cmdFN' => "$infoHR->{'analysis_dirs'}{'fastq-info'}/$infoHR->{'optHR'}{'fq1b'}.fastQC", 'cont' => 1, 'xn' => 2},
          {'cmd' => "fastqc -o $infoHR->{'analysis_dirs'}{'fastq-results'}/ $infoHR->{'optHR'}{'fq2'}", 'echo' => 1, 'cmdFN' => "$infoHR->{'analysis_dirs'}{'fastq-info'}/$infoHR->{'optHR'}{'fq2b'}.fastQC", 'cont' => 1, 'xn' => 2},

          {'cmd' => "# Filter reads"},
          {'cmd' => "$infoHR->{'configHR'}{'BIN'}{'filter_fastq'} --outdir $infoHR->{'analysis_dirs'}{'fastq_filter'}/ --prefix $infoHR->{'optHR'}{'prefix'} --suffix filter --fq1 $infoHR->{'optHR'}{'fq1'} --fq2 $infoHR->{'optHR'}{'fq2'}",
           'echo' => 1, 'cmdFN' => "$infoHR->{'analysis_dirs'}{'fastq_filter-info'}/$infoHR->{'optHR'}{'prefix'}.filter_fastq", 'cont' => 0, 'xn' => 2},

          {'cmd' => "# Submit next sbatch script(s)"},
          {'cmd' => "sbatch $infoHR->{'analysis_dirs'}{'sbatch'}/$scriptNameAR->[$cIdx + 1]"},
          {'cmd' => "sbatch $infoHR->{'analysis_dirs'}{'sbatch'}/$scriptNameAR->[$cIdx + 2]", 'xn' => 2},
        );

        $cpHR->{'inpPfx'} = "$infoHR->{'analysis_dirs'}{'fastq_filter'}/$infoHR->{'optHR'}{'prefix'}";
        $cpHR->{'infPfx'} = "$infoHR->{'analysis_dirs'}{'fastq_filter-info'}/$infoHR->{'optHR'}{'prefix'}";
        $cpHR->{'resPfx'} = "$infoHR->{'analysis_dirs'}{'fastq_filter-results'}/$infoHR->{'optHR'}{'prefix'}";
        push(@{$cmdAR},
          {'cmd' => "# FastQC after filter"},
          {'cmd' => "fastqc -o $infoHR->{'analysis_dirs'}{'fastq_filter-results'}/ $cpHR->{'inpPfx'}.1.filter.fastq", 'echo' => 1, 'cmdFN' => "$cpHR->{'infPfx'}.PE.1.filter_fastq_fastQC", 'cont' => 1, 'xn' => 2},
          {'cmd' => "fastqc -o $infoHR->{'analysis_dirs'}{'fastq_filter-results'}/ $cpHR->{'inpPfx'}.2.filter.fastq", 'echo' => 1, 'cmdFN' => "$cpHR->{'infPfx'}.PE.2.filter_fastq_fastQC", 'cont' => 1, 'xn' => 2},
          {'cmd' => "fastqc -o $infoHR->{'analysis_dirs'}{'fastq_filter-results'}/ $cpHR->{'inpPfx'}.S.filter.fastq", 'echo' => 1, 'cmdFN' => "$cpHR->{'infPfx'}.SE.S.filter_fastq_fastQC", 'cont' => 1, 'xn' => 2},
        );

        my $sbatchFH = myOpenRW("$infoHR->{'analysis_dirs'}{'sbatch'}/$scriptNameAR->[$cIdx]");
        print $sbatchFH returnSbatchScript($infoHR, $paramAR, $cmdAR);
        close($sbatchFH);
      }


      { # Mosaik PE
        my $cpHR = {};
        my $cIdx = 1;
        my $cmdAR = [];

        my $paramAR = returnParamAR($infoHR, {'scriptNamePgm' => $scriptPgmNameAR->[$cIdx], 'C' => "fat", 't' => "50:00:00"});
        returnInitAndModuleLoad($infoHR, $cmdAR, 'Bioinfo-tools', 'Mosaik');

        $cpHR->{'inpPfx'} = "$infoHR->{'analysis_dirs'}{'fastq_filter'}/$infoHR->{'optHR'}{'prefix'}";
        $cpHR->{'outPfx'} = "$infoHR->{'analysis_dirs'}{'mosaik'}/$infoHR->{'optHR'}{'prefix'}";
        $cpHR->{'infPfx'} = "$infoHR->{'analysis_dirs'}{'mosaik-info'}/$infoHR->{'optHR'}{'prefix'}";
        $cpHR->{'resPfx'} = "$infoHR->{'analysis_dirs'}{'mosaik-results'}/$infoHR->{'optHR'}{'prefix'}";

        my $shmRef = "/dev/shm/MosaikRef/"; my $shmMDS = "/dev/shm/MosaikDupSnoop/";
#        my $mRefB = myBasename($infoHR->{'configHR'}{'PATH'}{'REFERENCE_MOSAIK'});
        $cpHR->{'nMosPfx'} = "$infoHR->{'configHR'}{'PATH'}{'MOSAIK_NODE'}/$infoHR->{'optHR'}{'prefix'}";
        $cpHR->{'nInfPfx'} = "$infoHR->{'configHR'}{'PATH'}{'MOSAIK_NODE'}/info/$infoHR->{'optHR'}{'prefix'}";

        push(@{$cmdAR},
          {'cmd' => "# Setup Mosaik TMP dir"},
          {'cmd' => "mkdir -p $infoHR->{'configHR'}{'PATH'}{'MOSAIK_TEMP'}"},
          {'cmd' => "export MOSAIK_TMP=$infoHR->{'configHR'}{'PATH'}{'MOSAIK_TEMP'}"},
          {'cmd' => "mkdir -p $infoHR->{'configHR'}{'PATH'}{'MOSAIK_NODE'}/info/"},
          {'cmd' => "mkdir -p $shmRef"},
          {'cmd' => "mkdir -p $shmMDS", 'xn' => 2},
#          {'cmd' => "cp -f $infoHR->{'configHR'}{'PATH'}{'REFERENCE_MOSAIK'} $shmRef", 'xn' => 2},

          {'cmd' => "# MosaikBuild"},
          {'cmd' => "MosaikBuild $infoHR->{'configHR'}{'PAR'}{'MosaikBuild_PE'} -q $cpHR->{'inpPfx'}.1.filter.fastq -q2 $cpHR->{'inpPfx'}.2.filter.fastq -out $cpHR->{'nMosPfx'}.PE.dat",
           'echo' => 1, 'cmdFN' => "$cpHR->{'nInfPfx'}.PE.MosaikBuild", 'cont' => 0, 'xn' => 2},
          {'cmd' => "rsync -ca $infoHR->{'configHR'}{'PATH'}{'MOSAIK_NODE'} $infoHR->{'analysis_dirs'}{'mosaik'}/", 'xn' => 2},
        );

        $cpHR->{'inpPfx'} = "$infoHR->{'analysis_dirs'}{'mosaik'}/$infoHR->{'optHR'}{'prefix'}";
        push(@{$cmdAR},
          {'cmd' => "# MosaikAligner"},
          {'cmd' => "MosaikAligner $infoHR->{'configHR'}{'PAR'}{'MosaikAligner_PE'} -in $cpHR->{'nMosPfx'}.PE.dat -out $cpHR->{'nMosPfx'}.PE.aln.dat -rur $cpHR->{'nMosPfx'}.PE.unAln.fastq " .
           "-ia $infoHR->{'configHR'}{'PATH'}{'REFERENCE_MOSAIK'} -j $infoHR->{'configHR'}{'PATH'}{'REFERENCE_MOSAIK_JMP'}",
           'echo' => 1, 'cmdFN' => "$cpHR->{'nInfPfx'}.PE.MosaikAligner", 'cont' => 0, 'xn' => 2},
          {'cmd' => "rsync -ca $infoHR->{'configHR'}{'PATH'}{'MOSAIK_NODE'} $infoHR->{'analysis_dirs'}{'mosaik'}/", 'xn' => 2},
        );

        push(@{$cmdAR},
          {'cmd' => "# MosaikDupSnoop"},
          {'cmd' => "MosaikDupSnoop $infoHR->{'configHR'}{'PAR'}{'MosaikDupSnoop_PE'} -in $cpHR->{'nMosPfx'}.PE.aln.dat -od $shmMDS",
           'echo' => 1, 'cmdFN' => "$cpHR->{'nInfPfx'}.PE.MosaikDupSnoop", 'cont' => 0, 'xn' => 1},
#          {'cmd' => "cp -f $shmMDS.db $infoHR->{'analysis_dirs'}{'mosaik-mosaikDS_PE'}/", 'xn' => 2},
          {'cmd' => "rsync -ca $infoHR->{'configHR'}{'PATH'}{'MOSAIK_NODE'} $infoHR->{'analysis_dirs'}{'mosaik'}/", 'xn' => 2},
        );

        push(@{$cmdAR},
          {'cmd' => "# MosaikSort"},
          {'cmd' => "MosaikSort $infoHR->{'configHR'}{'PAR'}{'MosaikSort_PE'} -in $cpHR->{'nMosPfx'}.PE.aln.dat -out $cpHR->{'nMosPfx'}.PE.aln.mSrt.dat -dup $shmMDS",
           'echo' => 1, 'cmdFN' => "$cpHR->{'nInfPfx'}.PE.MosaikSort", 'cont' => 0, 'xn' => 2},
          {'cmd' => "rsync -ca $infoHR->{'configHR'}{'PATH'}{'MOSAIK_NODE'} $infoHR->{'analysis_dirs'}{'mosaik'}/", 'xn' => 2},
        );

        push(@{$cmdAR},
          {'cmd' => "# MosaikMerge & MosaikText"},
          {'cmd' => "rm -f $infoHR->{'optHR'}{'sbatch-cur-PE'}"},
          {'cmd' => "if [ ! \"\$(ls -A $infoHR->{'analysis_dirs'}{'sbatch-cur'})\" ]; then"}, # Alternative: [ "$(ls -l 12-10F/mosaik/*.aln.mSrt.dat|wc -l)" == 2 ] && echo "2 files" || echo "not 2 files"
          {'cmd' => "MosaikMerge $infoHR->{'configHR'}{'PAR'}{'MosaikMerge'} \$(find $infoHR->{'analysis_dirs'}{'mosaik'}/*.aln.mSrt.dat|perl -nle 'printf(\"-in %s \", \$_)') -out $infoHR->{'analysis_dirs'}{'mosaik'}/$infoHR->{'optHR'}{'sample'}.merged.dat",
           'echo' => 1, 'cmdFN' => "$cpHR->{'infPfx'}.PE.MosaikMerge", 'cont' => 0, 'xi' => 2, 'xn' => 2},
          {'cmd' => "MosaikText -in $infoHR->{'analysis_dirs'}{'mosaik'}/$infoHR->{'optHR'}{'sample'}.merged.dat -bam $infoHR->{'analysis_dirs'}{'mosaik'}/$infoHR->{'optHR'}{'sample'}.merged.bam",
           'echo' => 1, 'cmdFN' => "$cpHR->{'infPfx'}.PE.MosaikText", 'cont' => 0, 'xi' => 2, 'xn' => 2},

          {'cmd' => "# Submit next sbatch script(s)", 'xi' => 2},
          {'cmd' => "sbatch $infoHR->{'analysis_dirs'}{'sbatch'}/$scriptNameAR->[$cIdx + 2]", 'xi' => 2, 'xn' => 2},

          {'cmd' => "fi"},
          {'cmd' => "rsync -ca $infoHR->{'configHR'}{'PATH'}{'MOSAIK_NODE'} $infoHR->{'analysis_dirs'}{'mosaik'}/", 'xn' => 2},
        );

        my $sbatchFH = myOpenRW("$infoHR->{'analysis_dirs'}{'sbatch'}/$scriptNameAR->[$cIdx]");
        print $sbatchFH returnSbatchScript($infoHR, $paramAR, $cmdAR);
        close($sbatchFH);
      }


      { # Mosaik SE
        my $cpHR = {};
        my $cIdx = 2;
        my $cmdAR = [];

        my $paramAR = returnParamAR($infoHR, {'scriptNamePgm' => $scriptPgmNameAR->[$cIdx], 'C' => "fat", 't' => "50:00:00"});
        returnInitAndModuleLoad($infoHR, $cmdAR, 'Bioinfo-tools', 'Mosaik');

        $cpHR->{'inpPfx'} = "$infoHR->{'analysis_dirs'}{'fastq_filter'}/$infoHR->{'optHR'}{'prefix'}";
        $cpHR->{'outPfx'} = "$infoHR->{'analysis_dirs'}{'mosaik'}/$infoHR->{'optHR'}{'prefix'}";
        $cpHR->{'infPfx'} = "$infoHR->{'analysis_dirs'}{'mosaik-info'}/$infoHR->{'optHR'}{'prefix'}";
        $cpHR->{'resPfx'} = "$infoHR->{'analysis_dirs'}{'mosaik-results'}/$infoHR->{'optHR'}{'prefix'}";

        my $shmRef = "/dev/shm/MosaikRef/"; my $shmMDS = "/dev/shm/MosaikDupSnoop/";
#        my $mRefB = myBasename($infoHR->{'configHR'}{'PATH'}{'REFERENCE_MOSAIK'});
        $cpHR->{'nMosPfx'} = "$infoHR->{'configHR'}{'PATH'}{'MOSAIK_NODE'}/$infoHR->{'optHR'}{'prefix'}";
        $cpHR->{'nInfPfx'} = "$infoHR->{'configHR'}{'PATH'}{'MOSAIK_NODE'}/info/$infoHR->{'optHR'}{'prefix'}";

        push(@{$cmdAR},
          {'cmd' => "# Setup Mosaik TMP dir"},
          {'cmd' => "mkdir -p $infoHR->{'configHR'}{'PATH'}{'MOSAIK_TEMP'}"},
          {'cmd' => "export MOSAIK_TMP=$infoHR->{'configHR'}{'PATH'}{'MOSAIK_TEMP'}"},
          {'cmd' => "mkdir -p $infoHR->{'configHR'}{'PATH'}{'MOSAIK_NODE'}/info/"},
          {'cmd' => "mkdir -p $shmRef"},
          {'cmd' => "mkdir -p $shmMDS", 'xn' => 2},
#          {'cmd' => "cp -f $infoHR->{'configHR'}{'PATH'}{'REFERENCE_MOSAIK'} $shmRef", 'xn' => 2},

          {'cmd' => "# MosaikBuild"},
          {'cmd' => "MosaikBuild $infoHR->{'configHR'}{'PAR'}{'MosaikBuild_SE'} -q $cpHR->{'inpPfx'}.S.filter.fastq -out $cpHR->{'nMosPfx'}.SE.dat",
           'echo' => 1, 'cmdFN' => "$cpHR->{'nInfPfx'}.SE.MosaikBuild", 'cont' => 0, 'xn' => 2},
          {'cmd' => "rsync -ca $infoHR->{'configHR'}{'PATH'}{'MOSAIK_NODE'} $infoHR->{'analysis_dirs'}{'mosaik'}/", 'xn' => 2},
        );

        $cpHR->{'inpPfx'} = "$infoHR->{'analysis_dirs'}{'mosaik'}/$infoHR->{'optHR'}{'prefix'}";
        push(@{$cmdAR},
          {'cmd' => "# MosaikAligner"},
          {'cmd' => "MosaikAligner $infoHR->{'configHR'}{'PAR'}{'MosaikAligner_SE'} -in $cpHR->{'nMosPfx'}.SE.dat -out $cpHR->{'nMosPfx'}.SE.aln.dat -rur $cpHR->{'nMosPfx'}.SE.unAln.fastq " .
           "-ia $infoHR->{'configHR'}{'PATH'}{'REFERENCE_MOSAIK'} -j $infoHR->{'configHR'}{'PATH'}{'REFERENCE_MOSAIK_JMP'}",
           'echo' => 1, 'cmdFN' => "$cpHR->{'nInfPfx'}.SE.MosaikAligner", 'cont' => 0, 'xn' => 2},
          {'cmd' => "rsync -ca $infoHR->{'configHR'}{'PATH'}{'MOSAIK_NODE'} $infoHR->{'analysis_dirs'}{'mosaik'}/", 'xn' => 2},
        );

        push(@{$cmdAR},
          {'cmd' => "# MosaikDupSnoop"},
          {'cmd' => "MosaikDupSnoop $infoHR->{'configHR'}{'PAR'}{'MosaikDupSnoop_SE'} -in $cpHR->{'nMosPfx'}.SE.aln.dat -od $shmMDS",
           'echo' => 1, 'cmdFN' => "$cpHR->{'nInfPfx'}.SE.MosaikDupSnoop", 'cont' => 0, 'xn' => 1},
          {'cmd' => "cp -f $shmMDS.db $infoHR->{'analysis_dirs'}{'mosaik-mosaikDS_SE'}/", 'xn' => 2},
          {'cmd' => "rsync -ca $infoHR->{'configHR'}{'PATH'}{'MOSAIK_NODE'} $infoHR->{'analysis_dirs'}{'mosaik'}/", 'xn' => 2},
        );

        push(@{$cmdAR},
          {'cmd' => "# MosaikSort"},
          {'cmd' => "MosaikSort $infoHR->{'configHR'}{'PAR'}{'MosaikSort_SE'} -in $cpHR->{'nMosPfx'}.SE.aln.dat -out $cpHR->{'nMosPfx'}.SE.aln.mSrt.dat -dup $shmMDS",
           'echo' => 1, 'cmdFN' => "$cpHR->{'nInfPfx'}.SE.MosaikSort", 'cont' => 0, 'xn' => 2},
          {'cmd' => "rsync -ca $infoHR->{'configHR'}{'PATH'}{'MOSAIK_NODE'} $infoHR->{'analysis_dirs'}{'mosaik'}/", 'xn' => 2},
        );

        push(@{$cmdAR},
          {'cmd' => "# MosaikMerge & MosaikText"},
          {'cmd' => "rm -f $infoHR->{'optHR'}{'sbatch-cur-SE'}"},
          {'cmd' => "if [ ! \"\$(ls -A $infoHR->{'analysis_dirs'}{'sbatch-cur'})\" ]; then"},
          {'cmd' => "MosaikMerge $infoHR->{'configHR'}{'PAR'}{'MosaikMerge'} \$(find $infoHR->{'analysis_dirs'}{'mosaik'}/*.aln.mSrt.dat|perl -nle 'printf(\"-in %s \", \$_)') -out $infoHR->{'analysis_dirs'}{'mosaik'}/$infoHR->{'optHR'}{'sample'}.merged.dat",
           'echo' => 1, 'cmdFN' => "$cpHR->{'infPfx'}.SE.MosaikMerge", 'cont' => 0, 'xi' => 2, 'xn' => 2},
          {'cmd' => "MosaikText -in $infoHR->{'analysis_dirs'}{'mosaik'}/$infoHR->{'optHR'}{'sample'}.merged.dat -bam $infoHR->{'analysis_dirs'}{'mosaik'}/$infoHR->{'optHR'}{'sample'}.merged.bam",
           'echo' => 1, 'cmdFN' => "$cpHR->{'infPfx'}.PE.MosaikText", 'cont' => 0, 'xi' => 2, 'xn' => 2},

          {'cmd' => "# Submit next sbatch script(s)", 'xi' => 2},
          {'cmd' => "sbatch $infoHR->{'analysis_dirs'}{'sbatch'}/$scriptNameAR->[$cIdx + 1]", 'xi' => 2, 'xn' => 2},

          {'cmd' => "fi"},
          {'cmd' => "rsync -ca $infoHR->{'configHR'}{'PATH'}{'MOSAIK_NODE'} $infoHR->{'analysis_dirs'}{'mosaik'}/", 'xn' => 2},
        );

        my $sbatchFH = myOpenRW("$infoHR->{'analysis_dirs'}{'sbatch'}/$scriptNameAR->[$cIdx]");
        print $sbatchFH returnSbatchScript($infoHR, $paramAR, $cmdAR);
        close($sbatchFH);
      }


      { # GATK
        my $cpHR = {};
        my $cIdx = 3;
        my $cmdAR = [];

        my $paramAR = returnParamAR($infoHR, {'scriptNamePgm' => $scriptPgmNameAR->[$cIdx], 'C' => "fat", 't' => "50:00:00"});
        returnInitAndModuleLoad($infoHR, $cmdAR, 'Bioinfo-tools', 'GATK', 'Samtools');

        $cpHR->{'inpPfx'} = "$infoHR->{'analysis_dirs'}{'mosaik'}/$infoHR->{'optHR'}{'sample'}";
        $cpHR->{'infPfx'} = "$infoHR->{'analysis_dirs'}{'mosaik-info'}/$infoHR->{'optHR'}{'sample'}";
        $cpHR->{'outPfx'} = "$infoHR->{'analysis_dirs'}{'GATK'}/$infoHR->{'optHR'}{'sample'}";
        $cpHR->{'resPfx'} = "$infoHR->{'analysis_dirs'}{'GATK-results'}/$infoHR->{'optHR'}{'sample'}";
        $cpHR->{'intPfx'} = "$infoHR->{'analysis_dirs'}{'GATK-intm'}/$infoHR->{'optHR'}{'sample'}";
#        my $shmRef = "/dev/shm/MosaikRef/"; my $shmMDS = "/dev/shm/MosaikDupSnoop/"; my $mRefB = myBasename($infoHR->{'configHR'}{'PATH'}{'REFERENCE_MOSAIK'});
        push(@{$cmdAR},
          {'cmd' => "# Samtools Index"},
          {'cmd' => "samtools index $cpHR->{'inpPfx'}.merged.bam", 'echo' => 1, 'cmdFN' => "$cpHR->{'infPfx'}.merged.samtoolsIndex_1_mosaik", 'cont' => 0, 'xn' => 2},
        );

        $cpHR->{'infPfx'} = "$infoHR->{'analysis_dirs'}{'GATK-info'}/$infoHR->{'optHR'}{'sample'}";
        push(@{$cmdAR},
          {'cmd' => "# GATK recalibration step 1"},
          {'cmd' => "java -Xmx$infoHR->{'configHR'}{'PAR'}{'GATK_mem'} -jar $infoHR->{'configHR'}{'PATH'}{'GATK_jar'} $infoHR->{'configHR'}{'PAR'}{'GATK_recal1'} -R $infoHR->{'configHR'}{'PATH'}{'REFERENCE_GENOME'} --DBSNP $infoHR->{'configHR'}{'PATH'}{'dbSNP_131'} ".
           "-I $cpHR->{'inpPfx'}.merged.bam -recalFile $cpHR->{'intPfx'}.recal.data.csv",
           'echo' => 1, 'cmdFN' => "$cpHR->{'infPfx'}.merged.GATK_recal_1", 'cont' => 0, 'xn' => 2},
        );

        push(@{$cmdAR},
          {'cmd' => "# GATK recalibration step 2"},
          {'cmd' => "java -Xmx$infoHR->{'configHR'}{'PAR'}{'GATK_mem'} -jar $infoHR->{'configHR'}{'PATH'}{'GATK_jar'} $infoHR->{'configHR'}{'PAR'}{'GATK_recal2'} -R $infoHR->{'configHR'}{'PATH'}{'REFERENCE_GENOME'} -I $cpHR->{'inpPfx'}.merged.bam ".
           "-o $cpHR->{'outPfx'}.merged.recal.bam -recalFile $cpHR->{'intPfx'}.recal.data.csv",
           'echo' => 1, 'cmdFN' => "$cpHR->{'infPfx'}.merged.GATK_recal_1", 'cont' => 0, 'xn' => 2},
        );

        $cpHR->{'inpPfx'} = "$infoHR->{'analysis_dirs'}{'GATK'}/$infoHR->{'optHR'}{'sample'}";
        push(@{$cmdAR},
          {'cmd' => "# Samtools Index"},
          {'cmd' => "samtools index $cpHR->{'inpPfx'}.merged.recal.bam", 'echo' => 1, 'cmdFN' => "$cpHR->{'infPfx'}.merged.samtoolsIndex_2_recal", 'cont' => 0, 'xn' => 2},
        );

        push(@{$cmdAR},
          {'cmd' => "# GATK realignment step 1"},
          {'cmd' => "java -Xmx$infoHR->{'configHR'}{'PAR'}{'GATK_mem'} -jar $infoHR->{'configHR'}{'PATH'}{'GATK_jar'} $infoHR->{'configHR'}{'PAR'}{'GATK_realign1'} -R $infoHR->{'configHR'}{'PATH'}{'REFERENCE_GENOME'} -I $cpHR->{'inpPfx'}.merged.recal.bam ".
           "-o $cpHR->{'intPfx'}.merged.intervals",
           'echo' => 1, 'cmdFN' => "$cpHR->{'infPfx'}.merged.GATK_realign_1", 'cont' => 0, 'xn' => 2},
        );

        push(@{$cmdAR},
          {'cmd' => "# GATK realignment step 2"},
          {'cmd' => "java -Xmx$infoHR->{'configHR'}{'PAR'}{'GATK_mem'} -jar $infoHR->{'configHR'}{'PATH'}{'GATK_jar'} $infoHR->{'configHR'}{'PAR'}{'GATK_realign2'} -R $infoHR->{'configHR'}{'PATH'}{'REFERENCE_GENOME'} -I $cpHR->{'inpPfx'}.merged.recal.bam ".
           "-o $cpHR->{'inpPfx'}.merged.recal.realn.bam -targetIntervals $cpHR->{'intPfx'}.merged.intervals",
           'echo' => 1, 'cmdFN' => "$cpHR->{'infPfx'}.merged.GATK_realign_2", 'cont' => 0, 'xn' => 2},
        );

        push(@{$cmdAR},
          {'cmd' => "# Samtools Sort"},
          {'cmd' => "samtools sort $cpHR->{'inpPfx'}.merged.recal.realn.bam $cpHR->{'inpPfx'}.merged.recal.realn.reSrt", 'echo' => 1, 'cmdFN' => "$cpHR->{'infPfx'}.merged.samtoolsSort_1_realn", 'cont' => 0, 'xn' => 2},
        );

        push(@{$cmdAR},
          {'cmd' => "# Samtools Index"},
          {'cmd' => "samtools index $cpHR->{'inpPfx'}.merged.recal.realn.reSrt.bam", 'echo' => 1, 'cmdFN' => "$cpHR->{'infPfx'}.merged.samtoolsIndex_3_realn", 'cont' => 0, 'xn' => 2},
        );

        my $sbatchFH = myOpenRW("$infoHR->{'analysis_dirs'}{'sbatch'}/$scriptNameAR->[$cIdx]");
        print $sbatchFH returnSbatchScript($infoHR, $paramAR, $cmdAR);
        close($sbatchFH);

      }

    } elsif ( -f $infoHR->{'optHR'}{'fqS'} ) {
      die "Error: Single end fastq data processing not implemented yet\n";
    } else {
      die "Error: Cannot find input fastq files\n";
    }

  }


}


sub returnParamAR {
  my $infoHR = shift;
  my $optHR  = shift;
  my $paramAR = [
    "-A $infoHR->{'configHR'}{'VAR'}{'PROJECT_NAME'}",
    "-p node -n 8",
    "-C $optHR->{'C'}",
    "-t $optHR->{'t'}",
    "-J $optHR->{'scriptNamePgm'}_$infoHR->{'optHR'}{'prefix'}",
    "-o $infoHR->{'analysis_dirs'}{'sbatch-info'}/$infoHR->{'optHR'}{'prefix'}_$optHR->{'scriptNamePgm'}.sbatch.stdout",
    "-e $infoHR->{'analysis_dirs'}{'sbatch-info'}/$infoHR->{'optHR'}{'prefix'}_$optHR->{'scriptNamePgm'}.sbatch.stderr",
    "--mail-type=All",
    "--mail-user=$infoHR->{'configHR'}{'VAR'}{'E-MAIL'}"
  ];
  return $paramAR;
}


sub returnInitAndModuleLoad {
  my $infoHR = shift;
  my $cmdAR  = shift;
  push(@{$cmdAR}, {'cmd' => "echo \"Running on: \$(hostname)\"", 'xn' => 1}, {'cmd' => "# Loading modules"},);
  foreach my $mod ( @_ ) {
    die "Module $mod not recognized, please define (including version) in conf file\n" unless defined($infoHR->{'configHR'}{'MOD'}{$mod});
    push(@{$cmdAR}, {'cmd' => "$infoHR->{'configHR'}{'MOD'}{'LoadCmd'} $infoHR->{'configHR'}{'MOD'}{$mod}"},);
  }
  $cmdAR->[scalar(@{$cmdAR}) - 1]{'xn'} = 1;
  push(@{$cmdAR}, {'cmd' => "# Running analysis", 'xn' => 1},);
}


sub returnDateString {
  my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time);
  return sprintf("%04d-%02d-%02d_%02d-%02d-%02d", $year + 1900, $mon + 1, $mday, $hour, $min, $sec);
}


sub returnCmd {
  my $cmdHR   = shift;
  my $cmd     = fixPath($cmdHR->{'cmd'});
  my $cmdFile = fixPath($cmdHR->{'cmdFN'});

  $cmd =~ s/\n+/ /g; $cmd =~ s/ {2,}/ /g;
  $cmd =~ s/\/{2,}/\//g;
  $cmd =~ s/\&//g;
  return qq{$cmd 1> $cmdFile.stdout 2> $cmdFile.stderr} . ($cmdHR->{'cont'} ? " &" : "");
}


sub returnEchoCmd {
  my $cmdHR   = shift;
  my $cmd     = fixPath($cmdHR->{'cmd'});
  my $cmdFile = fixPath($cmdHR->{'cmdFN'});

  $cmd = returnCmd($cmdHR);
#  $cmd =~ s/\"/\\\"/g;
  return qq{echo "$cmd" > $cmdFile.cmdLine};
}


sub returnSbatchScript {
  my $infoHR  = shift;
  my $paramAR = shift;
  my $cmdAR   = shift;

  my $sbatchString = "#! /bin/bash -l\n";

  foreach my $param ( @{$paramAR} ) {
    $sbatchString .= "#SBATCH $param\n";
  }
  $sbatchString .= "\n";

  foreach my $cmdHR ( @{$cmdAR} ) {
    if ( $cmdHR->{'cmdFN'} ) {
      if ( $cmdHR->{'echo'} ) { $sbatchString .= (" " x $cmdHR->{'xi'}) . returnEchoCmd($cmdHR) . "\n\n"; }
      $sbatchString .= (" " x $cmdHR->{'xi'}) . returnCmd($cmdHR) . ("\n" x (1 + $cmdHR->{'xn'}));
    } else {
      $sbatchString .= (" " x $cmdHR->{'xi'}) . fixPath($cmdHR->{'cmd'}) . ("\n" x (1 + $cmdHR->{'xn'}));
    }
  }

  $sbatchString .= "wait\n\n";

  return $sbatchString;
}



#__END__


=head1 NAME

run_exome_analysis.pl - Automated analysis of sequence data using Mosaik and GATK

=head1 SYNOPSIS

run_exome_analysis.pl [options]

  Options:
   --help            Brief help message (not yet implemented)
   --verbose         Write some additional output (not yet implemented)
   --outdir          Output directory [required]
   --sample          Output file prefix [required]
   --run             Output file prefix [required]
   --fq1             PE fastq file, direction 1
   --fq2             PE fastq file, direction 2
   --fqS             SE fastq file (not yet implemented)
   --shm             Use shared memory (/dev/shm)

At least one fastq file is required. If paired end (PE) data is used, both fq1 and fq2 are required and must match internally.

=cut

#=head1 OPTIONS

#=over 8

#=item B<-help>

#Print a brief help message and exits.

#=back

#=head1 DESCRIPTION






