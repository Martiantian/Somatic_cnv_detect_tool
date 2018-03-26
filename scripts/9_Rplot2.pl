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
die "perl $0 [input p_value.list][input figdir][input Rplot_script] [input Rscript_path]"if @ARGV<4;
my ($filelist,$dir,$Rplot,$Rscript)=@ARGV;
open IN,"<$filelist" or die $!;
while (<IN>){
	chomp;
	my $name=(split/\//,$_)[-1];
	$name=(split/\./,$name)[0];
	`$Rscript $Rplot $_ $dir/$name $name`;
}
close IN;

