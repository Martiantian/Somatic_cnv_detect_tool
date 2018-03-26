#!/bin/perl -w
#Author: wangxiaofeng
#        shichang
#Email: wangxiaofeng@genomics.cn
#       shichang@genomics.cn
#File Name:
#Description:
# this program after extractInfo.pl chrno start length
#Edit History:
#2016-01-30 17:00:03  File created.
use strict;
die"perl $0 [input read_gc_file][input read_gc.stac_file][input ref_read_gc.stac_file][input chrY_interval][input bin_length][input outname]"if@ARGV<6;
my $read_gc_file=shift;
my $rgc_stac_file=shift;
my $control=shift;
my $interval=shift;
my $bin_len=shift;
$bin_len||=1000000;
my $outname=shift;

my %cnr;
my $name=(split /\//,$read_gc_file)[-1];
$name=(split /\./,$name)[0];
open OUTM,">$outname.mer"or die $!;
open OUT,">$outname.read_gc.cnr"or die $!;
open OUTC,">$outname.cnr.stac"or die $!;
open IN2,"<$control"or die $!;
my (@gcp,@reads);
while(my $LINE=<IN2>){
	chomp($LINE);
	my @array=split/\t/,$LINE;
	push @gcp,$array[0];
	push @reads,$array[1];
}
close IN2;
my $n=@gcp;

foreach my $i(0..($n-1)){
	open IN1,"<$rgc_stac_file"or die $!;
	while(my $line=<IN1>){
		chomp($line);
		my @array=split/\t/,$line;
		if($array[0]==$gcp[$i]){
			my $cnr=$array[1]/$reads[$i];
			$cnr=sprintf("%.3f",$cnr);
			print OUTM "$line\t$reads[$i]\t$cnr\n";
			$cnr{$array[0]}=$cnr;
			last;
		}
	}
	close IN1;
}
close OUTM;



open IN2,"<$read_gc_file"or die $!;
while(<IN2>){
	chomp;
	my @lines=split;

	if(exists($cnr{$lines[-1]}) and $lines[-1]>=200 and $lines[-1]<=700){
		$lines[-1]=$cnr{$lines[-1]};
		print OUT join("\t",@lines),"\n";
	}
}
close OUT;
close IN2;

#my @position;
my %position;

open IN,"<$outname.read_gc.cnr" or die $!;
while(<IN>){
	chomp;
	my @array=split;
	my $ch;
	if ($array[0] eq "chrX"){$ch=23;}
	if($array[0] eq "chrY"){$ch=24;}
	for my $k (1 .. 22){
		if ($array[0] eq "chr$k"){
			$ch=$k;
			last;
		}
	}
	if($ch==24){
		open INI,"<$interval" or die $!;
		while(<INI>){
			chomp;
			my @lines=split;
			if($array[1]>=$lines[1] and $array[1]<=$lines[2]){
				$position{$ch}{$lines[-1]}+=1/$array[-1];
			}
		}
		close INI;
	}else{
		if($ch){
			$position{$ch}{int($array[1]/$bin_len)}+=1/$array[-1];
		}
	}
}
close IN;


for my $i (1 .. 24){
	if($position{$i}){
		for my $j(sort{$a<=>$b} keys %{$position{$i}}){
				my $reads=int($position{$i}{$j});
				print OUTC "$i\t$j\t$reads\n";
		}
	}
}

