#!/bin/perl -w
#Author: wangxiaofeng
#Email:wangxiaofeng@genomics.cn
#File Name:
#Description:
#  
#Edit History:
#2016-01-30 17:00:03  File created.
use strict;
use FindBin qw($Bin);
die "[input ref_hg19]  [input gc_file_type 1 or 2] [input stat_gc_length(bp)]"if @ARGV<2;
my $ref=shift;
my $type=shift;
my $len=shift;
my $path = $Bin;
my @chr=(1..24);
if ($type==1){
	$len=170 unless(defined($len));
	my %hash=%{&getSeq($ref)};
	my %stat;
	open OUT,">>$path/data/hg19.$len.1.gc"or die $!;
	foreach my $k(@chr){
		$k=~s/23/X/;
		$k=~s/24/Y/;
		$k="chr".$k;
		foreach my $kk(0..(length($hash{$k})-1)){
			if(substr($hash{$k},$kk,1) eq "N"){
				print OUT "0\n";
			}else{
				my $seq=substr($hash{$k},$kk,$len);
				my $GC=($seq=~tr/GCgc/GCgc/);
				my $gc=$GC/length($seq);
				my $g=sprintf("%.3f",$gc);
				$g=$g*1000;
				if ($k ne "chrX" && $k ne "chrY"){
					$stat{$g}+=1;
					}
				print OUT "$g\n";
			}
		}
	}
	open OUT2 ,">>$path/data/hg19.$len.no_sex.gc.stac" or die $!;
	foreach my $m (sort{$a<=>$b} keys %stat){
		print OUT2 "$m\t$stat{$m}\n";
		}
	close OUT;
	close OUT2;
	`gzip $path/data/hg19.$len.1.gc`;
}elsif($type==2){
	$len=1000000 unless(defined($len));
	my %hash=%{&getSeq($ref)};
	open OUT,">$path/data/hg19.$len.2.gc"or die $!;
	foreach my $k(@chr){
		$k=~s/23/X/;
		$k=~s/24/Y/;
		$k="chr".$k;
		my $start=0;
		while($start<=length($hash{$k})){
			my $seq=substr($hash{$k},$start,$len);
			$seq=~s/N//g;
			my $GC=($seq=~tr/GCgc/GCgc/);
			my $gc;
			my $n=$start/$len;
			if(length($seq)>0){
				$gc=$GC/length($seq);
				my $g=sprintf"%.6f",$gc;
				print OUT "$k\t$n\t$g\n";
			}else{
				print OUT "$k\t$n\t0\n";
			}
			$start+=$len;
		}
	}
}else{
	die "Please input gc_file_type 1 or 2 !";
}

sub getSeq {
	my $file = shift;
	my %HASH;
	my $i = 0;
	open IN,"<$file" or die $!;
	while (<IN>) {
		chomp;
		if (/^>/) {
				$i=$_;$i=~s/>//g;
		} else {
				s/\s+//g;
				$HASH{$i} .= $_ ;
		}
	}
close IN;
return \%HASH;
}