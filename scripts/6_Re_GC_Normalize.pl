#! /bin/perl  -w 
#Author: qianzhaoyang 
#        shichang
#Email: qianzhaoyang@genomics.cn
#       shichang@genomics.cn
#File Name:
#Description:
#  
#Edit History:
#2017-12-10 14:20:03  File created.
use strict;
die "perl $0 [input GCref file][input Filelist] [input lm-2.R] [input Rscript__path] [input outdir]"if @ARGV<5;
my $ref=shift;
my $file=shift;
my $lm_R=shift;
my $Rscript=shift;
my $dir=shift;
my @gc;
open REF,"<$ref" or die $!;
while (<REF>){
	my @array=split;
	$array[0]=~ s/X/23/;
	$array[0]=~ s/Y/24/;
	$array[0]=~ s/chr//;
	$gc[$array[0]][$array[1]]=$array[2];
	}
close REF;
open FILE,"<$file" or die $!;
while(<FILE>){
	chomp;
	open INPUT,"$_"or die $!;
	my $name=(split/\//,$_)[-1];
	$name=(split/\./,$name)[0];
	open OUTPUT,">$dir/$name.GC"or die $!;
	while(<INPUT>){
		chomp;
		my @data=split;
		print OUTPUT "$data[0]\t$data[1]\t$data[2]\t$gc[$data[0]][$data[1]]\n";
	}
	close OUTPUT;
	my $out=`$Rscript $lm_R $dir/$name.GC`;
	my @out2=split /\s+/,$out;
	open IN,"<$dir/$name.GC"or die $!;
	open OUT,">>$dir/$name.reNGC"or die $!;
	while(<IN>){
		chomp;
		my @seq=split;
		print OUT join("\t",@seq[0 .. 3]),"\t",$out2[-8]*$seq[-2]/&max(0.1,($out2[-3]+$out2[-2]*$seq[-1]+$out2[-1]*$seq[-1]**2)),"\n";
	}
	close OUT;
	sub max{return $_[0]>$_[1]?$_[0]:$_[1]}
}
close FILE;
