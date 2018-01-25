#!/usr/bin/perl -w
#Author: ***
#Email:***@genomics.cn
#File Name:
#Description:
#  
#Edit History:
#2016-01-30 17:00:03  File created.
use strict;
die "perl $0 [input bamFile][input maximum_quality] [input PE TorF]" if(@ARGV<2);
my ($bamFile,$maxq,$PE)=@ARGV;
$maxq=60 unless defined($maxq);
if($bamFile=~/.bam/){
	open IN,"samtools view $bamFile |"or die $!;
}else{
	open IN,"$bamFile"or die $!;
}
while(my $LINE=<IN>){
	chomp($LINE);
	next if $LINE=~/^\s*$/;
	my @line=split /\t/,$LINE;
	next if($line[0]=~/^@/);
	next unless(defined($line[1]));
	my $b=$line[2];
	$b=~s/chr//;
	$b=~s/X/23/;
	$b=~s/Y/24/;
	next if($b=~/[a-zA-Z]/);
	next if $line[1]>=1024 or $line[1]==4;
	next if ($line[5] eq "*");
	next if ($line[4] <($maxq*0.5));
	my $start;
	if($PE eq "F"){
		my $num=length($line[9]);
		if($line[1]==16){
			$start=$line[3]+$num-170;
		}else{
			$start=$line[3]
		}
	}elsif($PE eq "T"){
		next if $line[8]==0;
		if($line[8]>0){
			$start=$line[3];
		}elsif($line[8]<0){
			$start=$line[7];
		}
	}
	print "chr$b\t$start\n"; 
}
close IN;
