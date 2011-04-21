# FROM run_exome_analysis.pl

sub setParameters {
  my $infoHR = shift;

  $infoHR->{'run'}{'MosaikBuild'}{'par'}    = "-mfl 400 -st illumina";
  $infoHR->{'run'}{'MosaikAligner'}{'par'}  = "-p 8 -m all -mhp 100 -act 20 -bw 35 -mm 4";
  $infoHR->{'run'}{'MosaikDupSnoop'}{'par'} = "-afl -rmm";
  $infoHR->{'run'}{'MosaikSort'}{'par'}     = "-afl -rmm -sa -mem 12000000";
  $infoHR->{'run'}{'MosaikMerge'}{'par'}    = "-mem 12000000";
  $infoHR->{'run'}{'MosaikText'}{'par'}     = "";

  $infoHR->{'run'}{'GATK_recal1'}{'par'}    = "-l INFO -et STDOUT -T CountCovariates -cov ReadGroupCovariate -cov QualityScoreCovariate -cov CycleCovariate -cov DinucCovariate";
  $infoHR->{'run'}{'GATK_recal2'}{'par'}    = "-l INFO -et STDOUT -T TableRecalibration";
  $infoHR->{'run'}{'GATK_realign1'}{'par'}  = "-l INFO -et STDOUT -T RealignerTargetCreator";
  $infoHR->{'run'}{'GATK_realign2'}{'par'}  = "-l INFO -et STDOUT -T IndelRealigner";
}

sub setupBiodbPaths {
  my $infoHR = shift;

  $infoHR->{'biodb'}{'refFA'}           = "~/hg/GRCh37_hg19/GRCh37_hg19_allChr_nl.fasta";
  $infoHR->{'biodb'}{'mosaikRef'}       = "~/hg/GRCh37_hg19/mosaik/GRCh37_hg19_allChr_nl.dat";
  $infoHR->{'biodb'}{'mosaikJmp'}       = "~/hg/GRCh37_hg19/mosaik/GRCh37_hg19_allChr_nl_15";
  $infoHR->{'biodb'}{'dbSNP'}           = "~/hg/snp_GRCh37_hg19/snp131.noAlt.sorted.txt";
}

sub setupSampleNameAndPaths {
  my $infoHR = shift;

  $infoHR->{'dirs'}{'fastq'}   = "fastq";

  $infoHR->{'dirs'}{'mosaik'}  = "mosaik";
  $infoHR->{'dirs'}{'info'}    = "info";
  $infoHR->{'dirs'}{'ds'}      = "MosaikDupSnoop";
  $infoHR->{'files'}{'dsDB'}   = ".db";

  $infoHR->{'dirs'}{'GATK'}    = "GATK";

  $infoHR->{'dirs'}{'vars'}    = "variants";

  $infoHR->{'dirs'}{'shm'}      = "/dev/shm";
  $infoHR->{'dirs'}{'pwd'}      = `pwd`; chomp($infoHR->{'dirs'}{'pwd'});
  $infoHR->{'dirs'}{'basename'} = `basename $infoHR->{'dirs'}{'pwd'}`; chomp($infoHR->{'dirs'}{'basename'});

  $infoHR->{'filenameSyntax'} = ['sampleName', 'runName', 'PE', 'dir', 'suffix']; # 09-11147.hiSeq20100624.PE.1.fastq
}

sub setupLocalMosaikDirs {
  my $infoHR   = shift;
  my $setupShm = shift;
  my $debug    = shift;
  my $localDirsHR = {};

  createDirs($infoHR->{'dirs'}{'mosaik'}) unless $debug;
  $infoHR->{'dirs'}{'mosaik_info'} = "$infoHR->{'dirs'}{'mosaik'}/$infoHR->{'dirs'}{'info'}"; createDirs($infoHR->{'dirs'}{'mosaik_info'}) unless $debug;
  $infoHR->{'dirs'}{'mosaik_MDS'}  = "$infoHR->{'dirs'}{'mosaik'}/$infoHR->{'dirs'}{'MDS'}"; createDirs($infoHR->{'dirs'}{'mosaik_MDS'}) unless $debug;

  my $timeStamp = localtime(); $timeStamp =~ s/\s+/\-/g; $timeStamp =~ s/:/\./g;
  $infoHR->{'files'}{'logfile'}    = "$infoHR->{'dirs'}{'mosaik_info'}/$infoHR->{'smplInfo'}{'sampleName'}_" . $timeStamp . "_" . $$ . "_MosaikPipeline.log";
  $infoHR->{'misc'}{'majSpacer'} = "=" x 40;
  $infoHR->{'misc'}{'minSpacer'} = "-" x 40;

  if ( $setupShm ) {
    $infoHR->{'dirs'}{'shm_sample'} = "$infoHR->{'dirs'}{'shm'}/$infoHR->{'smplInfo'}{'sampleName'}"; createDirs($infoHR->{'dirs'}{'shm_sample'}) unless $debug;
    $infoHR->{'dirs'}{'shm_mosaik'} = "$infoHR->{'dirs'}{'shm_sample'}/$infoHR->{'dirs'}{'mosaik'}"; createDirs($infoHR->{'dirs'}{'shm_mosaik'}) unless $debug;
    $infoHR->{'dirs'}{'shm_info'}   = "$infoHR->{'dirs'}{'shm_mosaik'}/info"; createDirs($infoHR->{'dirs'}{'shm_info'}) unless $debug;
    $infoHR->{'dirs'}{'shm_MDS'}    = "$infoHR->{'dirs'}{'shm_mosaik'}/MosaikDupSnoop"; createDirs($infoHR->{'dirs'}{'shm_MDS'}) unless $debug;
    $infoHR->{'dirs'}{'shm_ia'}     = "$infoHR->{'dirs'}{'shm'}/mosaikRef"; createDirs($infoHR->{'dirs'}{'shm_ia'}) unless $debug;
  }
}

sub creatSbatchScript {
  my $infoHR  = shift;
  my $paramHR = shift;



  my $echoCmd = returnMosaikEchoCmd($infoHR);
  my $cmd     = returnMosaikRunCmd($infoHR);

  my $shmCpCmd = "cp -f $infoHR->{'biodb_R'}{'mosaikRef'} $infoHR->{'dirs'}{'shm_ia'}/" if $infoHR->{'optHR'}{'shm'};


  return <<END;
#! /bin/bash -l
#SBATCH -A $infoHR->{'configHR'}{'VAR'}{'PROJECT_NAME'};
#SBATCH -p node -n 8
#SBATCH -C fat
#SBATCH -t 50:00:00
#SBATCH -J $infoHR->{'remote'}{$infoHR->{'MosaikPgmAR'}[0]}{'jobname'}
#SBATCH -o $infoHR->{'remote'}{$infoHR->{'MosaikPgmAR'}[0]}{'jobStdout'}
#SBATCH -e $infoHR->{'remote'}{$infoHR->{'MosaikPgmAR'}[0]}{'jobStderr'}
#SBATCH --mail-type=All
#SBATCH --mail-user=$infoHR->{'optHR'}{'mail'}


echo "Running on: \\\$(hostname)"

# Load modules
module load bioinfo-tools
module load mosaik-aligner/1.0.1388

# cd to sample dir
cd $infoHR->{'dirs'}{'remSample'}/

# Setup shm dirs
mkdir -p $infoHR->{'dirs'}{'shm_info'}
mkdir -p $infoHR->{'dirs'}{'shm_MDS'}
mkdir -p $infoHR->{'dirs'}{'shm_ia'}
$shmCpCmd

# Log command
$echoCmd

# Run command
$cmd

# Exit status
exit_status=\\\${1:-\\\$?}
echo "$infoHR->{'remote'}{$infoHR->{'MosaikPgmAR'}[0]}{'jobname'}:\\\${exit_status}"
exit \\\${exit_status}

wait

END
}

#        {'cmd' => "echo \"Running on: \$(hostname)\"", 'xn' => 1},
#        {'cmd' => "# Loading modules"},
#        {'cmd' => "$infoHR->{'configHR'}{'MOD'}{'LoadCmd'} $infoHR->{'configHR'}{'MOD'}{'Bioinfo-tools'}"},
#        {'cmd' => "$infoHR->{'configHR'}{'MOD'}{'LoadCmd'} $infoHR->{'configHR'}{'MOD'}{'FastQC'}"},
#        {'cmd' => "$infoHR->{'configHR'}{'MOD'}{'LoadCmd'} $infoHR->{'configHR'}{'MOD'}{'Mosaik'}", 'xn' => 1},
#        {'cmd' => "# Running analysis", 'xn' => 1},


## Exit status
#exit_status=${1:-$?}
#echo "971-08D_MosaikAligner:${exit_status}"
#exit ${exit_status}







