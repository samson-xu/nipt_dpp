package ReadsAlign;

use File::Basename;
use Data::Dumper;
use WriteShell;

require Exporter;
our @ISA = qw(Exporter);
our @EXPORT = qw(reads_align);

sub reads_align {
	my $bwa = shift;
	my $samtools = shift;
	my $gatk3 = shift; 
	my $thread = shift;
	my $prefix = shift;
	my $fastq = shift;
	my $ref = shift;
	my $platform = shift;
	my $dbsnp = shift;
	my $mills = shift;
	my $tindels = shift;
	my $way = shift;
	my $arg = shift;
	my $outDir = shift;
	my $rm = shift;
	my $fq1 = $fastq->[0];
	my $fq2 = $fastq->[1];
	my $shell = "# bwa for alignment\n";
	if ($way eq 'backtrack') {
		$arg ||= "-t $thread";
		$shell .= "$bwa aln $arg $ref $fq1 > $outDir/$prefix.aln_sa1.sai &\n";
		if ($fq2) {
			$shell.=<<ALN;
if [ -e $fq2 ]
then
	$bwa aln $arg $ref $fq2 > $outDir/$prefix.aln_sa2.sai &
	wait
	$bwa sampe -r '\@RG\\tID:$prefix\\tSM:$prefix\\tPL:$platform' $ref $outDir/$prefix.aln_sa1.sai $outDir/$prefix.aln_sa2.sai $fq1 $fq2 > $outDir/$prefix.sam
else
	wait
	$bwa samse -r '\@RG\\tID:$prefix\\tSM:$prefix\\tPL:$platform' $ref $outDir/$prefix.aln_sa1.sai $fq1 > $outDir/$prefix.sam
fi
rm $outDir/$prefix.aln_sa*
ALN
		} else {
			$shell .= "wait\n";
			$shell .= "$bwa samse -r '\@RG\\tID:$prefix\\tSM:$prefix\\tPL:$platform' $ref $outDir/$prefix.aln_sa1.sai $fq1 > $outDir/$prefix.sam\n";
			$shell .= "rm $outDir/$prefix.aln_sa*\n";
		}
	} else {
		$arg = "-K 100000000 -t $thread -Y" if ($arg eq '');
		if ($fq2) {
			$shell.=<<MEM;
if [ -e $fq2 ]
then
	$bwa mem $arg -R '\@RG\\tID:$prefix\\tSM:$prefix\\tPL:$platform' $ref $fq1 $fq2 > $outDir/$prefix.sam
else
	$bwa mem $arg -R '\@RG\\tID:$prefix\\tSM:$prefix\\tPL:$platform' $ref $fq1 > $outDir/$prefix.sam
fi
MEM
		} else {
			$shell .= "$bwa mem $arg -R '\@RG\\tID:$prefix\\tSM:$prefix\\tPL:$platform' $ref $fq1 > $outDir/$prefix.sam\n";
		}
	}
	$shell.=<<DEAL;
$samtools sort -n -T $outDir --no-PG --threads $thread -o $outDir/$prefix.sortn.bam $outDir/$prefix.sam
rm $outDir/$prefix.sam
$samtools fixmate -m --no-PG --threads $thread $outDir/$prefix.sortn.bam $outDir/$prefix.fixmate.bam 
rm $outDir/$prefix.sortn.bam
$samtools sort -T $outDir --no-PG --threads $thread -o $outDir/$prefix.sort.bam $outDir/$prefix.fixmate.bam
rm $outDir/$prefix.fixmate.bam
$samtools markdup -T $outDir --no-PG --threads $thread $outDir/$prefix.sort.bam $outDir/$prefix.markdup.bam
rm $outDir/$prefix.sort.bam
$samtools index -@ $thread $outDir/$prefix.markdup.bam
# Indel附近重新比对
## 确定重新比对的区域
java -Xmx5g -jar $gatk3 \\
-T RealignerTargetCreator \\
-nt $thread \\
-R $ref \\
-I $outDir/$prefix.markdup.bam \\
-o $outDir/$prefix.realn.intervals \\
-known $mills \\
-known $tindels
## 重新比对
java -Xmx5g -jar $gatk3 \\
-T IndelRealigner \\
-R $ref \\
-targetIntervals $outDir/$prefix.realn.intervals \\
-I $outDir/$prefix.markdup.bam \\
-o $outDir/$prefix.realn.bam \\
-known $mills \\
-known $tindels
rm $outDir/$prefix.markdup.bam*
# 比对文件碱基质量值校正
## 生成质量校正所需输入文件grp
java -Xmx5g -jar $gatk3 \\
-T BaseRecalibrator \\
-nct $thread \\
-R $ref \\
-I $outDir/$prefix.realn.bam \\
-knownSites $dbsnp \\
-knownSites $mills \\
-o $outDir/$prefix.realn.grp
## 输出质量校正后的数据
java -Xmx5g -jar $gatk3 \\
-T PrintReads \\
-nct $thread \\
-R $ref \\
-I $outDir/$prefix.realn.bam \\
-BQSR $outDir/$prefix.realn.grp \\
-o $outDir/$prefix.bam
rm $outDir/$prefix.realn.*

$samtools stats $outDir/$prefix.bam > $outDir/$prefix.bam.stats 

$samtools quickcheck $outDir/$prefix.bam > $outDir/$prefix.bam.check

DEAL
	write_shell($shell, "$outDir/$prefix.align.sh");
}

1;
