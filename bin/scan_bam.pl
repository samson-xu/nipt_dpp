#!/usr/bin/env perl
################################################################################
# File Name: scan_bam
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
my $help;
my $project = strftime("%Y%m%d-%H%M%S",localtime());
my $workDir = $ENV{'PWD'};
my $ref = '';
my $thread = 36;
my $samtools = '';
my $call = '';

# Guide
my $guide_separator = "=" x 150;
my $indent = " " x 15;
my $parameter_separator = "*" x 70;
my $guide=<<INFO;
$guide_separator
=
=$indent$indent$indent NAME: Scan Bam Pipeline(scan_bam) 
=$indent$indent$indent AUTHOR: xuxiangyang(xy_xu\@foxmail.com) 
=$indent$indent$indent VERSION: v1.0	2021/04/07
=                                          
$guide_separator


FUNCTIONS
$indent 1. Split Bam check.

PARAMETER
$indent $0 [options] bam.lst

$parameter_separator Basic $parameter_separator 
$indent --help                        Print this guide information 
$indent --project <str>               Project name, default "$project"
$indent --workDir <str>               Work directory, default "$workDir"
$indent --ref <str>                   Reference genome absolute path, default "$ref"
$indent --thread <i>                  Set the number of threads for the program to run, default "$thread"
$indent --samtools <str>              Samtools absolute path, default "$samtools"
$indent --call <str>                  Call absolute path, default "$call"

NOTE
$indent 1. Input file must be bam list.
$indent 2. Must be set path for samtools, call

EXAMPLE
$indent Example1: $0 --project result --workDir /path/for/workDir --ref hg19.fa --thread 36 --samtools /path/for/samtools --call /path/for/call bam.lst 
$guide_separator

INFO

# Parameter
GetOptions(
	"h|help" => \$help,
	"project=s" => \$project,
	"workDir=s" => \$workDir,
	"ref=s" => \$ref,
	"thread=i" => \$thread,
	"samtools=s" => \$samtools,
	"call=s" => \$call,
);

die $guide if (@ARGV == 0 || defined $help);

# Main
my $time;
my @bam_list;
$time = strftime("%Y/%m/%d %H:%M:%S",localtime());
print "\n\nStart Scan Bam at $time!\n\n";

my $projectDir = "$workDir/$project";
system("mkdir -p $projectDir/") == 0 || die $!;

open LIST, $ARGV[0] or die $!;
while (<LIST>) {
	chomp;
	push @bam_list, $_;
}
close LIST;

system("mkdir -p $projectDir/log/") == 0 || die $!;
system("mkdir -p $workDir/run/") == 0 || die $!;
my $region = `head -1 $ARGV[0] | awk -F '/' '{print \$4}'`;
chomp($region);
system("> $projectDir/log/$region.log") == 0 || die $!;
system("> $projectDir/log/$region.err") == 0 || die $!;

my $count = 0;
my $pm = new Parallel::ForkManager($thread);
foreach my $file (@bam_list) {
		$count++;
		$pm->start and next;
		scan($file);
		$pm->finish;
}
$pm->wait_all_children;

$time = strftime("%Y/%m/%d %H:%M:%S",localtime());
print "\n\nEnd Scan Bam at $time!\n\n";

sub scan {
	my $bam = shift;
	my @arr = split /\//, $bam;
	my $lib = $arr[-2];
	my $rg = $arr[-3];
	my $pre = $arr[-1];
	$pre =~ s/.bam//;
	#split bam check
	system("echo $bam > $workDir/run/$pre.bam.lst ") == 0 || die $!;
	system("echo $pre $count >> $projectDir/log/$rg.log") == 0 || die $!;
	my $exitcode = system("$call basetype --input $workDir/run/$pre.bam.lst --output $workDir/run/$pre --reference $ref --region $rg --mapq 30 --thread 1 --batch 1 --load");
	system("rm -rf $workDir/run/$pre*") == 0 || die $!;
	$exitcode = $exitcode >> 8;
	if ($exitcode != 0) {
		system("mkdir -p $projectDir/$rg/$lib/") == 0 || die $!;
		system("echo ERROR:$bam $exitcode >> $projectDir/log/$rg.err") == 0 || die $!;
		system("$samtools view -H -b $bam > $projectDir/$rg/$lib/$pre.bam") == 0 || die $!;
		system("$samtools index $projectDir/$rg/$lib/$pre.bam") == 0 || die $!;
	}
}
