#!/bin/perl -w
#Author: ***
#Email:***@genomics.org.cn
#File Name:
#Description:
# this program after extractInfo.pl chrno start length
#Edit History:
#2016-01-30 17:00:03  File created.
use strict;
die"perl $0 [input file_CNV.list][debin_length][bin_length]"if @ARGV<3;
my ($filelist,$debin_length,$bin_length)=@ARGV;
my $nb=$bin_length/$debin_length;
open IN,"<$filelist"or die $!;
while(<IN>){
	chomp;
	my $file=$_;
	open OUT,">$file.M"or die $!;
	my (@position,@numb);
	for my $x(1..24){
		my @chr;
		push(@position,[@chr]);
		push(@numb,[@chr]);
	}	

	open IN1,"$file"or die $!;
	while(<IN1>){
		chomp;
		my @t=split;
		my $b=$t[0];
		$position[$t[0]]->[int($t[1]/$nb)]+=$t[2]; 
		if($t[2]>0){
			$numb[$t[0]]->[int($t[1])/$nb]++;
		}
	}
	close IN1;
	for my $i (1 .. 24){
		if($position[$i]){
			for my $j (0 .. (scalar(@{$position[$i]})-1)){
				if($position[$i]->[$j]){
					my $reads=$position[$i]->[$j];
					if($numb[$i]->[$j]>=4){
					print OUT "$i\t$j\t$reads\t$numb[$i]->[$j]\t",$reads/$numb[$i]->[$j],"\n";
					}
				}
			}
		}
	}
	close OUT;
}
close IN;
