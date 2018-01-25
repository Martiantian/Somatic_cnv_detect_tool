#!/bin/perl -w
use strict;
#Author: ***
#Email:***@genomics.org.cn
#File Name:
#Description:
#  
#Edit History:
#2016-01-30 17:00:03  File created.
die "perl $0 [input Filelist][input 4.1_Normalize.R] [input Rscript_path]"if @ARGV<3;
my $filelist=shift;
my $normalize=shift;
my $Rscript=shift;
open IN,"$filelist"or die $!;
while(<IN>){
	chomp;
	my @array=split;
	my $file=$array[0];#$file=~s/\/.*//;$file.='.CNV';
	`$Rscript $normalize $file `;
	open INN,"normalized_nb"or die $!;
	my $nb=<INN>;
	chomp($nb);
	$nb=(split /\t/,$nb)[0];
	close INN;
	`rm -rf normalized_nb`;
	&normalize($file,$nb);
}
close IN;

sub normalize{
	my $file1=shift;
	my $nb=shift;
	my $out=$file1.".nr";
	open IN1,"$file1"or die $!;
	open OUT,">$out"or die $!;
	while (my $line=<IN1>){
		chomp($line);
		my @lines=split/\t/,$line;
		my $n=@lines;
		my $t=$lines[$n-1];#
		next unless $t>=0;###
		$lines[$n-1]=$t/$nb;#
		print OUT join("\t",$lines[0],$lines[1],$lines[$n-1]),"\n";
	}
	close IN1;
	close OUT;
}
