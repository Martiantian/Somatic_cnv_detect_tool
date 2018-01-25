#!/bin/perl -w
#Author: ***
#Email:***@genomics.org.cn
#File Name:
#Description:
# this program after extractInfo.pl chrno start length
#Edit History:
#2016-01-30 17:00:03  File created.
use strict;

die "perl $0 [input info_file][input ref_hg19.170reads.chr.gc.gz][input Chr_Length.txt][outdir/outname]"if @ARGV<4;
my ($file,$ref,$chrlen,$outfile)=@ARGV;
my %chr_start;
open OUTC,"<$chrlen"or die $!;
<OUTC>;
my $chr_len=0;
while(<OUTC>){
	chomp;
	my @array=split;
	$chr_start{$array[0]}=$chr_len;
	$chr_len+=$array[1];
}####
#$chr_start{"chrY"}=$chr_len;

####%chr_start={chr=>length}
close OUTC;
open OUT,">$outfile.read_gc"or die $!;
open OUT1,">$outfile.read_gc.stac"or die $!;
&spl($file,$ref);
close OUT;
close OUT1;

sub spl{
	my $info=shift;
	my $ref=shift;
	my %hash1;
	open IN,"gzip -dc $ref |" or die $!;
	my $i=0;

	open INS,"$info"or die $!;
	my $line;
	while(my $LINE=<INS>){
		chomp($LINE);
		my @array=split /\t/,$LINE;#
		$array[0]=~s/chr23/chrX/;
		$array[0]=~s/chr24/chrY/;
		next if($array[0] eq "chrM");
		while($i<($chr_start{$array[0]}+$array[1])and $i<3095677412){
			$line=<IN>;
			chomp($line);
			$i++;
		}

		print OUT "$array[0]\t$array[1]\t$line\n";
		next if($array[0] eq "chrY" or $array[0] eq "chrM" or $array[0] eq "chrX");
		if($line){
			$hash1{"$line"}++;
		}
	}
	close INS;

	foreach my $k(sort{$a<=>$b} keys %hash1){
		print OUT1 "$k\t$hash1{$k}\n";
	}
}
