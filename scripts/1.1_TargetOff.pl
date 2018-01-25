#!/usr/bin/perl
#Author: ***
#Email:***@genomics.org.cn
#File Name:
#Description:
#  
#Edit History:
#2016-01-30 17:00:03  File created.
use strict;
die "perl $0 [input bam.file][input bed][input maximum_quality] [PE TorF ]"if @ARGV<3;
my ($bam,$bed,$maxq,$PE)=@ARGV;
$maxq=60 unless defined($maxq);
my @TG;
open IN2,"$bed"or die $!;
while(<IN2>){
	chomp;
	my @t=split;
	$t[0]=~s/chr//;$t[0]=~s/X/23/;
	push(@TG,[@t]);
}
unshift(@TG,[(1,0,0)]);

my $tmps;
my $i=1;
if($bam=~/bam/){
	open IN,"samtools view $bam |"or die $!;
}elsif($bam=~/sam/){
	open IN,"$bam"or die $!;
}
while (my $LINE=<IN>){
	chomp($LINE);
	next if $LINE=~/^\s*$/;
	my @t=split /\t/,$LINE;
	$_=$t[0];
	next  if (/^@/);
	next unless(defined($t[1]));
	next if $t[1] >=1024 or $t[1]==4;
	next if ($t[5] eq "*");
	next if $t[4]<$maxq*0.5;
	my $b=$t[2];
	$b=~s/chr//;$b=~s/X/23/;
	next if $b =~/[a-zA-Z]/;
	$t[2]="chr$b";
	if($PE eq "F"){
		if($t[1]==16){
			$t[3]=$t[3]-length($t[9])+170;
		}
	}
	if($PE eq "T"){
		if($t[8]<0){
			$t[3]=$t[7];
		}
	}
	if($tmps){
		my @tm=split;
		if($b==$tm[0] and $t[3]<$tm[1]){
			while($t[3]>$TG[$i]->[1] and $b==$TG[$i]->[0] and $i-1>0 and $b==$TG[$i-1]->[0]){
				$i--;
			}
		}
	}
	$tmps="$b\t$t[3]";

	if($b<$TG[$i]->[0]){
		print $t[2],"\t",$t[3],"\n";
		next;
	}
	while($i<$#TG and ($b>$TG[$i]->[0] or ($b==$TG[$i]->[0] and $t[3]>$TG[$i]->[2]))){
		$i++;
	}
	if($t[3]<$TG[$i]->[1] or $b<$TG[$i]->[0]){
		print $t[2],"\t",$t[3],"\n";
	}
}