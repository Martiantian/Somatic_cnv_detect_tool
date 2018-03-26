#!/bin/perl -w
#Author: wangxiaofeng
#        shichang
#Email: wangxiaofeng@genomics.cn
#       shichang@genomics.cn
#File Name:
#Description:
#	
#Edit History:
#2016-01-30 17:00:03	File created.
use strict;
die "perl $0 [input y_ye.list][input bin_length(bp)][input confidence_alpha][input detect_threshold][input plot window length(kb)][input xlsdir][input min_CNV_length(bp)][input date] [input Detectplot_script] [input Rscript_path]"if @ARGV<9;
my ($filelist,$bin_length,$alpha,$threshold,$wl,$xlsdir,$cnv_length,$date,$Detectplot,$Rscript)=@ARGV;
$cnv_length=4000000 unless(defined($cnv_length));

open IN,"$filelist"or die $!;
my $i=0;
$xlsdir.="_$date";
if(-d $xlsdir){
	`rm -rf $xlsdir`;
}else{
	`mkdir $xlsdir`;
}

if(-e "detect_cnv_line_Pvalue.$date.xls"){
	`rm -rf detect_cnv_line_Pvalue.$date.xls`;
}

while(<IN>){
	chomp;
	$i++;
	my $name=(split /\//,$_)[-1];
	$name=~s/\.O.nr.*.E.ratio//;
	`$Rscript $Detectplot $name $threshold $wl $alpha $_ $bin_length $xlsdir/$name.cnv_line_p_value.xls $xlsdir/detect_cnv_line_Pvalue.$date.xls $xlsdir/merge_sd.$date.xls $i`;
	}
close IN;







