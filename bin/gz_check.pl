#!/usr/bin/env perl
################################################################################
# File Name: gz_check.pl
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
my $thread = 16;

# Guide
my $guide_separator = "=" x 150;
my $indent = " " x 15;
my $parameter_separator = "*" x 70;
my $guide=<<INFO;
$guide_separator
=
=$indent$indent$indent NAME: GZ Check Pipeline(gz_check) 
=$indent$indent$indent AUTHOR: xuxiangyang(xy_xu\@foxmail.com) 
=$indent$indent$indent VERSION: v1.0	2021/04/08
=                                          
$guide_separator


FUNCTIONS
$indent 1. cvg or vcf gz file check.

PARAMETER
$indent $0 [options] *.gz

$parameter_separator Basic $parameter_separator 
$indent --help                        Print this guide information 
$indent --project <str>               Project name, default "$project"
$indent --workDir <str>               Work directory, default "$workDir"
$indent --thread <i>                  Set the number of threads for the program to run, default "$thread"

NOTE
$indent 1. Input file must be cvg or vcf gz file.

EXAMPLE
$indent Example1: $0 --project result --workDir /path/for/workDir --thread 16 *gz 
$guide_separator

INFO

# Parameter
GetOptions(
	"h|help" => \$help,
	"project=s" => \$project,
	"workDir=s" => \$workDir,
	"thread=i" => \$thread,
);

die $guide if (@ARGV == 0 || defined $help);

# Main
my $time;
my @bam_list;
$time = strftime("%Y/%m/%d %H:%M:%S",localtime());
print "\n\nStart GZ check at $time!\n\n";

my $projectDir = "$workDir/$project";
system("mkdir -p $projectDir/") == 0 || die $!;
system("mkdir -p $workDir/run/") == 0 || die $!;
system("> $projectDir/gz.check.txt") == 0 || die $!;
my $pm = new Parallel::ForkManager($thread);
foreach my $input (@ARGV) {
		$pm->start and next;
		GZcheck($input);
		$pm->finish;
}
$pm->wait_all_children;

$time = strftime("%Y/%m/%d %H:%M:%S",localtime());
print "\n\nEnd GZ check at $time!\n\n";

sub GZcheck {
	my $file = shift;
	my $pre = basename($file);
	$pre =~ s/.gz$//;
	system("zgrep '^#' $file > $workDir/run/$pre.head") == 0 || die $!;
	my $md5 = `md5sum $workDir/run/$pre.head | awk '{print \$1}'`; 
	chomp($md5);
	my $stat = `zgrep -v '^#' $file | awk '{print NF}' | uniq -c | awk '{print \$2,\$1}'`;
	chomp($stat);
	system("echo $pre $md5 $stat >> $projectDir/gz.check.txt") == 0 || die $!;
	system("rm -rf $workDir/run/$pre*") == 0 || die $!;
}
