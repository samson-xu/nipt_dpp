#!/usr/bin/env perl
################################################################################
# File Name: bam_split.pl
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
use Parallel::ForkManager;
use FindBin qw($Bin);


# Global variable
my ($help, $region, $check, $rm);
my $project = strftime("%Y%m%d-%H%M%S",localtime());
my $workDir = $ENV{'PWD'};
my $ref = '';
my $thread = 36;
my $gtz = '';
my $samtools = '';
my $call = '';

# Guide
my $guide_separator = "=" x 150;
my $indent = " " x 15;
my $parameter_separator = "*" x 70;
my $guide=<<INFO;
$guide_separator
=
=$indent$indent$indent NAME: Bam Split Pipeline(bam_split) 
=$indent$indent$indent AUTHOR: xuxiangyang(xy_xu\@foxmail.com) 
=$indent$indent$indent VERSION: v1.0	2021/02/05
=                                          
$guide_separator


FUNCTIONS
$indent 1. Bam split base on region bed.
$indent 2. Split Bam check.

PARAMETER
$indent $0 [options] *.bam

$parameter_separator Basic $parameter_separator 
$indent --help                        Print this guide information 
$indent --project <str>               Project name, default "$project"
$indent --workDir <str>               Work directory, default "$workDir"
$indent --ref <str>                   Reference genome absolute path, default "$ref"
$indent --region <str>                Bed file for bam split, samtools-like region, chr:start-end, such as chr1:6000000-8000000, chr1 and so on 
$indent --thread <i>                  Set the number of threads for the program to run, default "$thread"
$indent --gtz <str>                   GTZ soft absolute path, default "$gtz"
$indent --samtools <str>              Samtools absolute path, default "$samtools"
$indent --call <str>                  Call absolute path, default "$call"
$indent --check                       whether check split bam, IONTORRENT should check 
$indent --rm                          rm all tmp files 

NOTE
$indent 1. Input file must be bam format.
$indent 2. Must be set path for gtz, samtools, call

EXAMPLE
$indent Example1: $0 --project result --workDir /path/for/workDir --ref hg19.fa  --region hg19.1m.final.bed --thread 36 --gtz /path/for/gtz --samtools /path/for/samtools --call /path/for/call *.bam 
$guide_separator

INFO

# Parameter
GetOptions(
	"h|help" => \$help,
	"project=s" => \$project,
	"workDir=s" => \$workDir,
	"ref=s" => \$ref,
	"region=s" => \$region,
	"thread=i" => \$thread,
	"gtz=s" => \$gtz,
	"samtools=s" => \$samtools,
	"call=s" => \$call,
	"check" => \$check,
	"rm" => \$rm,
);

die $guide if (@ARGV == 0 || defined $help);

# Main
my @regions;
open RG, $region or die $!;
while (<RG>) {
	next if (/#/);
	next if (/^\s+$/);
	chomp;
	my @arr = split /\s+/;
	push @regions, "$arr[0]:$arr[1]-$arr[2]";
}
close RG;

my $time;
$time = strftime("%Y/%m/%d %H:%M:%S",localtime());
print "\n\nStart Bam Split at $time!\n\n";
my $projectDir = "$workDir/$project";
system("mkdir -p $projectDir/") == 0 || die $!;
system("mkdir -p $projectDir/log/") == 0 || die $! if (defined $check);
foreach my $input (@ARGV) {
	my $file = "";
	my $lib = "";
	if ($input =~ /gtz$/) {
		system("$gtz -d -f --ref $ref -O $workDir -p $thread $input") == 0 || die $!;
		my @arr = split /\//, $input;
		$file = $arr[-1]; 
		$file =~ s/.gtz//;
		$lib = $arr[-3];
		system("$samtools index -@ $thread $workDir/$file") == 0 || die $!;
	} else {
		print STDERR "$input input format error! Plesse check it!\n";
		exit;
	}
	my $pm = new Parallel::ForkManager($thread);
	foreach my $part (@regions) {
			$pm->start and next;
			bsplit("$workDir/$file", $lib, $part, $samtools);
			$pm->finish;
	}
	$pm->wait_all_children;
	system("rm -rf $workDir/*.bam $workDir/*.bai") == 0 || die $! if (defined $rm);
}
$time = strftime("%Y/%m/%d %H:%M:%S",localtime());
print "\n\nEnd Bam Split at $time!\n\n";

sub bsplit {
	my $bam = shift;
	my $lib = shift;
	my $rg = shift;
	my $samtools = shift;
	my @arr = split /\//, $bam;
	my $pre = $arr[-1];
	$pre =~ s/.bam//;
	system("mkdir -p $projectDir/$rg/$lib/") == 0 || die $!;
	system("$samtools view -b $bam $rg > $projectDir/$rg/$lib/$pre.$rg.bam") == 0 || die $!;
	system("$samtools index $projectDir/$rg/$lib/$pre.$rg.bam") == 0 || die $!;
	#split bam check
	if (defined $check) {
		system("echo $projectDir/$rg/$lib/$pre.$rg.bam > $workDir/bam.lst ") == 0 || die $!;
		my $exitcode = system("$call basetype --input $workDir/bam.lst --output $workDir/$rg --reference $ref --region $rg --mapq 30 --thread 1 --batch 1 --maf 0.00001 --rerun");
		system("rm -rf $workDir/$rg*") == 0 || die $!;
		$exitcode = $exitcode >> 8;
		if ($exitcode != 0) {
			system("echo ERROR:$projectDir/$rg/$lib/$pre.$rg.bam $rg $exitcode >> $projectDir/log/$lib.log") == 0 || die $!;
			system("$samtools view -H -b $bam $rg > $projectDir/$rg/$lib/$pre.$rg.bam") == 0 || die $!;
			system("$samtools index $projectDir/$rg/$lib/$pre.$rg.bam") == 0 || die $!;
		}
	}	
}
