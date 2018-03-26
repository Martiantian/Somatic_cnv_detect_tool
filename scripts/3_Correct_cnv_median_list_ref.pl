#!/bin/perl -w
#Author: wangxiaofeng
#Email: wangxiaofeng@genomics.cn
#File Name:
#Description:
#  
#Edit History:
#2016-01-30 17:00:03  File created.
use strict;
die"perl $0 [input cnr_stac.list] [input mean_stac_name]"if @ARGV<2;
my $filelist=shift;
my $mean_cnv=shift;
$mean_cnv||="control_mean.cnr.stac";


my %hash;
open IN1,"$mean_cnv"or die $!;
while(<IN1>){
	chomp;
	my @array=split;
	$hash{"$array[0]_$array[1]"}=$array[-1];
}
close IN1;


open IN2,"$filelist";
while(my $file=<IN2>){
	chomp($file);
	my @array=split /\//,$file;
	my $name=(split /\./,$array[-1])[0];
	pop @array;
	my $dir=join "\/",@array;
	&correct_cnv($file,"$dir/$name");
}
close IN2;

sub correct_cnv{
	my $file=shift;
	my $name=shift;
	open IN3,"$file";
	open OUT,">$name.CNV";
	while(<IN3>){
		chomp;
		my @lines=split;
		if(exists($hash{"$lines[0]_$lines[1]"})){
			my $control=$hash{"$lines[0]_$lines[1]"};
			my $cnv=$lines[-1]/$control;
#				$lines[-1]=$cnv;
			print OUT join("\t",@lines),"\t$control\t$cnv","\n";
		}
	}
	close IN3;
}
