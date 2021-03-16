#!/usr/bin/env perl
################################################################################
# File Name: region_call.pl
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
use WriteShell;

# Global variable
my ($help, $main_shell, $region);
my $project = strftime("%Y%m%d-%H%M%S",localtime());
my $workDir = $ENV{'PWD'};
my $thread = 36;
my $run = 'no';
my $ref = '';
my $cp = '';
my $samtools = '';
my $call = '';
my $gtz = '';

# Guide
my $guide_separator = "=" x 150;
my $indent = " " x 15;
my $parameter_separator = "*" x 70;
my $guide=<<INFO;
$guide_separator
=
=$indent$indent$indent NAME: NIPT Region Call Pipeline(region_call) 
=$indent$indent$indent AUTHOR: xuxiangyang(xy_xu\@foxmail.com) 
=$indent$indent$indent VERSION: v1.0	2020/12/22
=                                          
$guide_separator


FUNCTIONS
$indent 1. Target region select.
$indent 2. Region SNP call.

PARAMETER
$indent $0 [options] bam.lst

$parameter_separator Basic $parameter_separator 
$indent --help                        Print this guide information 
$indent --project <str>               Project name, default "$project"
$indent --workDir <str>               Work directory, default "$workDir"
$indent --run <str>                   whether run pipeline, yes or no, default "$run"
$indent --thread <i>                  Set the number of threads for the program to run, default "$thread"
$indent --ref <str>                   Reference genome absolute path, default "$ref"
$indent --cp <str>                    Data cp soft absolute path, default "$ref"
$indent --samtools <str>              Samtools absolute path, default "$samtools"
$indent --call <str>                  SNP call soft absolute path, default "$call"
$indent --gtz <str>                   GTZ soft absolute path, default "$gtz"
$indent --region <str>                Target region or bed file for variant call, samtools-like region, chr:start-end, such as chr1:6000000-8000000, chr1 and so on 

NOTE
$indent 1. Input must be bam path list, one sample per line

EXAMPLE
$indent Example1: $0 --run y --thread 36 --ref hg19.fa --cp cp --samtools /path/for/samtools --call /path/for/call --gtz /path/for/gtz --region chr1:6000000-8000000 bam.lst
$guide_separator

INFO

# Parameter
GetOptions(
	"h|help" => \$help,
	"project=s" => \$project,
	"workDir=s" => \$workDir,
	"run=s" => \$run,
	"thread=i" => \$thread,
	"ref=s" => \$ref,
	"cp=s" => \$cp,
	"samtools=s" => \$samtools,
	"call=s" => \$call,
	"gtz=s" => \$gtz,
	"region=s" => \$region,
);

die $guide if (@ARGV == 0 || defined $help);

# Main
#$project = ".$project" if ($project !~ m/^\./);
my $projectDir = "$workDir/$project";
my $input = shift;
unless (-e $input) {
	print STDERR "$input isn't exist! Plesse check it!\n";
	exit;
}
system("mkdir -p $projectDir") == 0 || die $!;
# load region
my %rg;
if (-e $region) {
	open RG, $region or die $!;
	my $count = 0;
	while (<RG>) {
		next if (/^#/);
		next if (/^\s*$/);
		chomp;
		$count++;
		my @arr = split /\s/;
		my $value = "$arr[0]:$arr[1]-$arr[2]";
		$rg{$count} = $value;
	}
	close RG;
} else {
	if ($region) {
		$rg{1} = $region;
	} else {
		$rg{1} = "";
	}
}
# Target region select
my $batch_select_shell = "";
my %select_bam_lst = ();
open LIST, $input or die $!;
while (<LIST>) {
	next if (/#/);
	next if (/^\s*$/);
	chomp;
	my @arr = split "/";
	my $ossDir = dirname($_);
	my $bam = $_;
	my $sampleId = $arr[-2];
	my $libDir = "$projectDir/$arr[-3]";
	my $sample_shell=<<SS; 
mkdir -p $libDir
$cp cp $bam $libDir -u 
$samtools index $libDir/$sampleId.bam
SS
	foreach my $key (sort {$a<=>$b} keys %rg) {
		$sample_shell .= "$samtools view -b $libDir/$sampleId.bam $rg{$key} > $libDir/$sampleId.$rg{$key}.bam\n";
		$sample_shell .= "$samtools index $libDir/$sampleId.$rg{$key}.bam\n";
		$select_bam_lst{$key} .= "$libDir/$sampleId.$rg{$key}.bam\n";
	}		
	#$sample_shell .= "$gtz --ref $ref --donot-pack-ref -o $libDir/$sampleId.bam.gtz -e $libDir/$sampleId.bam\n";
	#$sample_shell .= "$cp cp $libDir/$sampleId.bam.gtz $ossDir/ -u\n";
	$sample_shell .= "rm $libDir/$sampleId.bam*\n";
	write_shell($sample_shell, "$libDir/$sampleId.select.sh");
	$batch_select_shell .= "sh $libDir/$sampleId.select.sh >$libDir/$sampleId.select.sh.o 2>$libDir/$sampleId.select.sh.e\n";
}
close LIST;

parallel_shell($batch_select_shell, "$projectDir/select.sh", $thread, 1);

foreach my $key (sort {$a<=>$b} keys %select_bam_lst) {
	open NB, ">$projectDir/call_bam.$key.lst" or die $!;
	print NB $select_bam_lst{$key};
	close NB;
}

# Region SNP call
my $call_shell = "";
foreach my $key (sort {$a<=>$b} keys %rg) {
	$call_shell .= "$call basetype --input $projectDir/call_bam.$key.lst --output $workDir/nipt.$rg{$key} --reference $ref --region $rg{$key} --mapq 30 --thread 4 --batch 10 --rerun\n";
}
write_shell($call_shell, "$projectDir/call.sh");

$main_shell = "# Run region_call pipeline for all samples\n";
$main_shell .= "sh $projectDir/select.sh >$projectDir/select.sh.o 2>$projectDir/select.sh.e\n";
$main_shell .= "sh $projectDir/call.sh >$projectDir/call.sh.o 2>$projectDir/call.sh.e\n";
#$main_shell.=<<MV;
#MV

write_shell($main_shell, "$projectDir/main.sh");

system("sh $projectDir/main.sh >$projectDir/main.sh.o 2>$projectDir/main.sh.e") == 0 || die $! if ($run =~ m/y/i);
