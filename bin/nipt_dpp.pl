#!/usr/bin/env perl
################################################################################
# File Name: nipt_dpp.pl
# Author: samson-xu
# mail: xy_xu@foxmail.com
# Created Time: Mon 22 Jun 2020 02:57:53 PM CST
################################################################################

use strict;
use warnings;
use File::Basename;
use Getopt::Long;
use Data::Dumper;
use POSIX;
use FindBin qw($Bin);
use lib "$Bin/../lib";
use ConfigParse;
use SampleStat;
use WriteShell;
use ReadsAlign;

# File or tool path check
my $config = path_check("$Bin/config.txt");

# Global variable
my ($help, $stat, $fastqc_help, $fastp_help, $backtrack_help, $mem_help, $fusion_help, %sample_shell, $main_shell);
my $project = strftime("%Y%m%d-%H%M%S",localtime());
my $workDir = $ENV{'PWD'};
my $prefix = "";
my $sample_match_filter = "";
my $platform = "ILLUMINA";
my $ref = $config->{'hg19'};
my $step = 'cfb';
my $thread = '35';
my $run = 'no';
my $rm = 'yes';
my $fastqc_arg = '';
my $fastp_arg = "-q 5 -u 50 -n 4 -e 10 -l 20 -w 4";
my $align_way = 'mem';
my $align_arg = '';

# Guide
my $guide_separator = "=" x 150;
my $indent = " " x 15;
my $parameter_separator = "*" x 70;
my $guide=<<INFO;
$guide_separator
=
=$indent$indent$indent NAME: NIPT Data Preprocessing Pipeline(nipt_dpp) 
=$indent$indent$indent AUTHOR: xuxiangyang(xy_xu\@foxmail.com) 
=$indent$indent$indent VERSION: v1.0	2020/11/18
=                                          
$guide_separator


FUNCTIONS
$indent c. Quality control and check(FastQC) of input data(FASTQ).
$indent f. Adapter cut and low quality sequence filter of fastq(Fastp).
$indent b. Fastq Alignment and quality control of sequence alignment results.

PARAMETER
$indent $0 [options] bcl_dir|fq.lst|ubam.lst

$parameter_separator Basic $parameter_separator 
$indent --help                        Print this guide information 
$indent --project <str>               Project name, default "$project"
$indent --workDir <str>               Work directory, default "$workDir"
$indent --prefix <str>                Prefix for this run, default "$prefix"
$indent --smf <str>                   Sample id match this string will be filtered, default "$sample_match_filter"
$indent --platform <str>              Platform/technology used to produce the read, such as ILLUMINA, SOLID, IONTORRENT, HELICOS and PACBIO, default "$platform"
$indent --ref <str>                   Reference genome absolute path, default "$ref"
$indent --step <str>                  Set step for run, default "$step"
$indent --run <str>                   whether run pipeline, yes or no, default "$run"
$indent --rm <str>                    whether rm tmp files, yes or no, default "$rm"
$indent --thread <i>                  Set the number of threads for the program to run, default "$thread"
$indent --stat                        Wether stat sample information, default not stat
$parameter_separator Filter $parameter_separator 
$indent --fastqc_help                 Print fastqc help information
$indent --fastqc_arg <str>            Fastqc argument setting, default "$fastqc_arg"
$indent --fastp_help                  Print fastp help information
$indent --fastp_arg <str>             Fastp argument setting, default "$fastp_arg"
$parameter_separator Align $parameter_separator 
$indent --align_way <str>             Select align algorithm, 'backtrack', 'mem', default "$align_way"
$indent --backtrack_help              Print BWA-backtrack help information
$indent --mem_help                    Print BWA-mem help information
$indent --align_arg <str>             Align argument setting, this has to correspond to the align_way, default "$align_arg"

NOTE
$indent 1. Input must be bcl dir, fastq list or unmapped bam list
$indent 2. If input is bcl dir, there must be SampleSheet.csv in the directory
$indent 3. If input is fastq list, format like: SampleId    fq1    fq2 or SampleId    fq1
$indent 4. If input is unmapped bam list, format like: SampleId    bam
$indent 5. Fastq quality system should be phred 33

EXAMPLE
$indent Example1: $0 --step cfb --run n bcl_dir 
$indent Example2: $0 --step cfb --run n fq.lst 
$indent Example3: $0 --step cfb --run n ubam.lst 
$guide_separator

INFO

# Parameter
GetOptions(
	"h|help" => \$help,
	"project=s" => \$project,
	"workDir=s" => \$workDir,
	"prefix=s" => \$prefix,
	"smf=s" => \$sample_match_filter,
	"platform" => \$platform,
	"ref=s" => \$ref,
	"step=s" => \$step,
	"thread=i" => \$thread,
	"run=s" => \$run,
	"rm=s" => \$rm,
	"stat" => \$stat,
	"fastqc_help" => \$fastqc_help,
	"fastqc_arg=s" => \$fastqc_arg,
	"fastp_help" => \$fastp_help,
	"fastp_arg=s" => \$fastp_arg,
	"align_way=s" => \$align_way, 
	"backtrack_help" => \$backtrack_help,
	"mem_help" => \$mem_help,
	"align_arg=s" => \$align_arg,
);

if (@ARGV == 0) {
	die `$config->{fastqc} -h` . "\n" if (defined $fastqc_help);
	die `$config->{fastp}` . "\n" if (defined $fastp_help);
	die `$config->{bwa} aln` . "\n" if (defined $backtrack_help);
	die `$config->{bwa} mem` . "\n" if (defined $mem_help);
}
die $guide if (@ARGV == 0 || defined $help);

# Main
$project = ".$project" if ($project !~ m/^\./);
my $projectDir = "$workDir/$project";
my $input = shift;
unless (-e $input) {
	print STDERR "$input isn't exist! Plesse check it!\n";
	exit;
}
system("mkdir -p $projectDir") == 0 || die $!;
$main_shell = "# Run nipt_dpp pipeline for all samples\n";
# load fastq info
my %sampleInfo;
if ($input eq "fq.lst" ) {
	open FQ, $input or die $!;
	while (<FQ>) {
		next if (/#/);
		chomp;
		my @arr = split /\s+/;
		next if ($arr[0] =~ m/$sample_match_filter/);
		@{$sampleInfo{$arr[0]}{'fastq'}} = ($arr[1]) if (@arr == 2);	
		@{$sampleInfo{$arr[0]}{'fastq'}} = ($arr[1], $arr[2]) if (@arr == 3);	
	}
	close FQ;
} elsif ($input eq "ubam.lst") {
	my $bam2fastq_shell = "";	
	system("mkdir -p $projectDir/fastq") == 0 || die $!;
	open UBAM, $input or die $!;
	while (<UBAM>) {
		next if (/#/);
		chomp;
		my @arr = split /\s+/;
		next if ($arr[0] =~ m/$sample_match_filter/);
		$bam2fastq_shell .= "$config->{'samtools'} fastq --threads 4 -0 $projectDir/fastq/$arr[0].1.fq.gz $arr[1]\n";
		@{$sampleInfo{$arr[0]}{'fastq'}} = ("$projectDir/fastq/$arr[0].1.fq.gz");
	}
	close UBAM;
	parallel_shell($bam2fastq_shell, "$projectDir/bam2fastq.sh", $thread, 4);
	$main_shell .= "sh $projectDir/bam2fastq.sh >$projectDir/bam2fastq.sh.o 2>$projectDir/bam2fastq.sh.e\n";
} else {
	system("mkdir -p $projectDir/fastq") == 0 || die $!;
	my $bcl2fastq_shell=<<BCL;
# convert bcl to fastq
$config->{'bcl2fastq'} -R $input -o $projectDir/tmp -r $thread -p $thread -w $thread --no-lane-splitting
rm $projectDir/tmp/Undetermined*
mv $projectDir/tmp/*gz $projectDir/tmp/*/*gz $projectDir/fastq
rm -rf $projectDir/tmp
BCL
	open SS, "$input/SampleSheet.csv" or die $!;
	my $data_label = 0;
	while (<SS>) {
		if (/Data/) {
			$data_label = 1;
			next;
		}
		if ($data_label) {
			next if (/#/);
			next if (/^Sample/i);
			next if (/^\s+$/);
			chomp;
			my @arr = split ",";
			next if ($arr[1] =~ m/$sample_match_filter/);
			@{$sampleInfo{$arr[1]}{'fastq'}} = ("$projectDir/fastq/$arr[1]*R1*gz", "$projectDir/fastq/$arr[1]*R2*gz");
		}
	}
	close SS;
	write_shell($bcl2fastq_shell, "$projectDir/bcl2fastq.sh");
	$main_shell .= "sh $projectDir/bcl2fastq.sh >$projectDir/bcl2fastq.sh.o 2>$projectDir/bcl2fastq.sh.e\n";
}
#print Dumper \%sampleInfo;
my $sample_total = keys %sampleInfo;

my ($batch_fastq_shell, $batch_align_shell);

foreach my $sampleId (sort {$a cmp $b} keys %sampleInfo) {
	# Fastq quality control
	my $fastq = $sampleInfo{$sampleId}{'fastq'};
	if ($step =~ /c/) {
		my $fastqcDir = "$projectDir/$sampleId/00.fastqc";
		my $fastqc_shell=<<FASTQC;
if [ -e $fastq->[1] ]
then
	$config->{fastqc} -t 2 -o $fastqcDir $fastqc_arg $fastq->[0] $fastq->[1]	
else
	$config->{fastqc} -t 2 -o $fastqcDir $fastqc_arg $fastq->[0]
fi
rm $fastqcDir/*.zip
FASTQC
		write_shell($fastqc_shell, "$fastqcDir/$sampleId.fastqc.sh");
		$batch_fastq_shell .= "sh $fastqcDir/$sampleId.fastqc.sh >$fastqcDir/$sampleId.fastqc.sh.o 2>$fastqcDir/$sampleId.fastqc.sh.e\n";
	}
	# Fastq filter
	if ($step =~ /f/) {
		my $filterDir = "$projectDir/$sampleId/01.filter";
		my $filter_shell=<<FILTER;
if [ -e $fastq->[1] ]
then
	$config->{fastp} -i $fastq->[0] -o $filterDir/$sampleId.clean.1.fq.gz -I $fastq->[1] -O $filterDir/$sampleId.clean.2.fq.gz --detect_adapter_for_pe $fastp_arg -j $filterDir/$sampleId.fastq.json -h $filterDir/$sampleId.fastq.html -R '$sampleId fastq report'
else
	$config->{fastp} -i $fastq->[0] -o $filterDir/$sampleId.clean.1.fq.gz $fastp_arg -j $filterDir/$sampleId.fastq.json -h $filterDir/$sampleId.fastq.html -R '$sampleId fastq report'
fi
FILTER
		#$filter_shell .= "perl -I '$Bin/../lib' -MReadsStat -e \"reads_stat('$filterDir/$sampleId.fastq.json')\"\n";
		write_shell($filter_shell, "$filterDir/$sampleId.filter.sh");
		$batch_fastq_shell .= "sh $filterDir/$sampleId.filter.sh >$filterDir/$sampleId.filter.sh.o 2>$filterDir/$sampleId.filter.sh.e\n";
		@{$sampleInfo{$sampleId}{'clean'}} = ("$filterDir/$sampleId.clean.1.fq.gz", "$filterDir/$sampleId.clean.2.fq.gz");
	}
	# Alignment
	if ($step =~ /b/) {
		my $alignDir = "$projectDir/$sampleId/02.align";
		#@{$sampleInfo{$sampleId}{'clean'}} = ($fastq->[0], $fastq->[1]) unless ($sampleInfo{$sampleId}{'clean'});
		my $align_program = $config->{'bwa'};
		reads_align($align_program, $config->{'samtools'}, $config->{'gatk3'}, $thread, $sampleId, $sampleInfo{$sampleId}{'clean'},
                    $ref, $platform, $config->{'dbsnp'}, $config->{'mills'}, $config->{'tindels'}, $align_way, $align_arg, $alignDir, $rm);
		#$sample_shell{$sampleId} .= "sh $alignDir/$sampleId.align.sh >$alignDir/$sampleId.align.sh.o 2>$alignDir/$sampleId.align.sh.e\n";
		$batch_align_shell .= "sh $alignDir/$sampleId.align.sh >$alignDir/$sampleId.align.sh.o 2>$alignDir/$sampleId.align.sh.e\n";
		$sampleInfo{$sampleId}{'align'} = "$alignDir/$sampleId.final.bam"; 
	}
}
parallel_shell($batch_fastq_shell, "$projectDir/fastq.deal.sh", $thread, 3) if ($step =~ /c/ or $step =~ /f/);
parallel_shell($batch_align_shell, "$projectDir/align.sh", $thread, 8) if ($step =~ /b/);
$main_shell .= "sh $projectDir/fastq.deal.sh >$projectDir/fastq.deal.sh.o 2>$projectDir/fastq.deal.sh.e\n" if ($step =~ /c/ or $step =~ /f/);
$main_shell .= "sh $projectDir/align.sh >$projectDir/align.sh.o 2>$projectDir/align.sh.e\n" if ($step =~ /b/);

#foreach my $sampleId (sort {$a cmp $b} keys %sampleInfo) {
#	# SNP/InDel detection
#	if ($step =~ /v/) {
#	}
#}
#
#if ($step =~ /c|f|b|v/) {
#	foreach my $sampleId (sort {$a cmp $b} keys %sampleInfo) {
#		write_shell($sample_shell{$sampleId}, "$projectDir/$sampleId/$sampleId.sh");
#		$main_shell .= "sh $projectDir/$sampleId/$sampleId.sh >$projectDir/$sampleId/$sampleId.sh.o 2>$projectDir/$sampleId/$sampleId.sh.e\n";
#	}
#}
#
#if ($step =~ /f/) {
#	$main_shell.=<<PSFQ;
#paste $projectDir/*/01.filter/*.fq.stat.txt | awk '{for(i=3; i<=NF; i+=2){\$i=""}; print \$0}' | sed "s/\\s\\+/\\t/g" > $projectDir/sample.fq.stat.xls
#PSFQ
#}
#
#if ($step =~ /b/) { 
#	$main_shell.=<<PSBM;
#paste $projectDir/*/02.align/*.bam.stat.txt | awk '{for(i=3; i<=NF; i+=2){\$i=""}; print \$0}' | sed "s/\\s\\+/\\t/g" > $projectDir/sample.bam.stat.xls
#PSBM
#}
#
if ($step =~ /f/ and $step =~ /b/) {
	#$main_shell .= "cat $projectDir/sample.fq.stat.xls $projectDir/sample.bam.stat.xls | sed '7d' | sed '8,10d'> $projectDir/sample.stat.xls\n";
	$main_shell.=<<MV;
#mkdir -p $workDir/${prefix}Result 
#mv $projectDir/*/01.filter/{*.fq.gz,*.html,*.json} $workDir/${prefix}Result
#mv $projectDir/*/02.align/*bam* $workDir/${prefix}Result
MV
}

#$main_shell .= "rm $projectDir/*/01.filter/*.gz\n" if ($step =~ /f/ and $rm =~ m/y/i);

write_shell($main_shell, "$projectDir/main.sh");

stat_log($sample_total, $Bin) if (defined $stat);

system("nohup sh $projectDir/main.sh >$projectDir/main.sh.o 2>$projectDir/main.sh.e &") == 0 || die $! if ($run =~ m/y/i);
