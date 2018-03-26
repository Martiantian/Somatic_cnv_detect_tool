#!/bin/perl 
#Author: wangxiaofeng
#        shichang
#Email: wangxiaofeng@genomics.cn
#       shichang@genomics.cn
#File Name:
#Description:
#  
#Edit History:
#2016-01-30 17:00:03  File created.
use strict;
use warnings;
use Cwd;
my $dir=cwd;
die "perl $0 [input filelist][input alpha_confidence][input cnv_size][input standard_sequenced_data][input expected_minmum_proportion][input used_sequence_data_file][input plot_windows_length][input outfile_data_dir][input 7_Normalize.R script][input Rscript__path][input detect_chrnums[1 or 1..23]]" if @ARGV<10;
my ($filelist,$alpha,$cnvsize,$readcount,$proportion,$usdf,$wlength,$datadir,$NormalizeR,$Rscript,@chrname)=@ARGV;
$alpha=0.05 unless ($alpha);
$wlength=250 unless ($wlength);
$datadir="$dir/result_data" unless ($datadir);
@chrname=(1..23) unless @chrname;


open RF, "|R --vanilla --slave";
select RF;

open TMP,">$datadir/tmp.file"or die $!;
print TMP "$datadir/result$alpha\n";
open OUT,">$datadir/result$alpha"or die $!;
print OUT "Name\tChr_NO\tRawdata_bins\tTruedata_bins\tFit_parameter\tLes\tCV\tCV_control\n";
my (@yeeline,@em,@ey,@sumy,@sumy2,@emin);

open OUTCV,">$datadir/cv.xls"or die $!;
####################################
open INT,"$usdf"or die $!;
my %hash_n;
while(<INT>){
	chomp;
	my @array=split;
	next if($array[-1] eq 'total');
	my $name=$array[-1];
	$name=(split /\//,$name)[-1];
	$name=~s/.sort.info//;
	$hash_n{$name}=$readcount/$array[0];
}
close INT;
####################################

open IN1,"<$filelist"or die $!;
while(<IN1>){
	chomp;
	my @array=split;
	my $file=$array[0];
	my $name=(split /\//,$file)[-1];
	$name=(split /\./,$name)[0];
	@yeeline=();@em=();@ey=();@sumy=();@sumy2=();@emin=();
	&myline($file,$name,"$datadir/$name.O.E.ratio",1,@chrname);
	`$Rscript $NormalizeR $datadir/$name.O.E.ratio $alpha $datadir/$name.O.nr$alpha.ratio $datadir/$name`;
	&myline("$datadir/$name.O.nr$alpha.ratio",$name,"$datadir/$name.O.nr$alpha.E.ratio",2,@chrname);
}
close OUT;
close IN1;
close RF;

sub initialize_cv{
	my $dtime=shift;
	my $sam=shift;
	my @y=@_;
	my $cv_bin_num=10;
	$cv_bin_num=int(10*$hash_n{$sam});
	if($cv_bin_num>40){
		$cv_bin_num=40;
	}
	if($cv_bin_num<10){
		$cv_bin_num=10;
	}
	my $mean_cv_num=6;
	my @chrname=(1..22);
	my $cv0;
	if($dtime==1){
		my @icv;
		foreach my $chrnum(@chrname){
			my $nnum=@{$y[$chrnum-1]};
			if($nnum==0){
				next; 
			}
			my $i=1;
			my ($tsumy,$tsumy2);
			foreach my $k(1..$nnum){
				$tsumy+=$y[$chrnum-1][$k-1];
				$tsumy2+=$y[$chrnum-1][$k-1]**2;
				$sumy[$k-1]=$tsumy;
				$sumy2[$k-1]=$tsumy2;
				if($k == $cv_bin_num*$i){
					my ($tmp_sumy2,$tmp_sumy);
					if($i==1){
						$tmp_sumy2=$sumy2[$k-1];
						$tmp_sumy=$sumy[$k-1];
					}else{
						my $tmp_num=$k-$cv_bin_num;
						$tmp_sumy2=$sumy2[$k-1]-$sumy2[$tmp_num-1];
						$tmp_sumy=$sumy[$k-1]-$sumy[$tmp_num-1];
					}
					my $tmp_mean=$tmp_sumy/$cv_bin_num;
					my $tmp_var=($tmp_sumy2-$cv_bin_num*$tmp_mean**2)/($cv_bin_num-1);
					my $cv=$tmp_var**0.5/$tmp_mean if($tmp_mean>0);
					push @icv,$cv;
					$i++;
				}
			}
		}
		@icv=sort{$a<=>$b} @icv;
		my $sum_icv;
		$mean_cv_num=int(@icv/6);
		foreach my $k(@icv[0..($mean_cv_num-1)]){
			$sum_icv+=$k;
		}
		$cv0=$sum_icv/$mean_cv_num;
	}else{
		open IN,"<$datadir/$sam.cv0.normal.txt"or die $!;
		$cv0=<IN>;
		chomp($cv0);
		close IN;
	}
	return($cv0);
}

sub myline{
	my $file=shift;
	my $sam=shift;
	my $out_y_ye=shift;
	my $dtime=shift;
	my @chrname=@_;
	open OUT3,">$out_y_ye"or die $!;

##########stat y( uneq 0)num per chromosome
	open IN2,"<$file"or die $!;
	my (@x,@y,@chr,@count,@nnum);
	foreach my $k(@chrname){
		push @count,0;
		push @nnum,0;
	}

	while(my $line=<IN2>){
		if($line=~/^\s*$/){
			$line=~s/\s*//;
			next;
		}
		chomp($line);
		my @tmp=split /\s+/,$line;
		if(grep{$_ eq $tmp[0]} @chrname){
			$count[$tmp[0]-1]++;
			next if($tmp[-1] eq 0);
			$chr[$tmp[0]-1][$nnum[$tmp[0]-1]]=$tmp[0];
			$x[$tmp[0]-1][$nnum[$tmp[0]-1]]=$tmp[1];
			$y[$tmp[0]-1][$nnum[$tmp[0]-1]]=$tmp[-1];
			$nnum[$tmp[0]-1]++;
		}else{
			next;
		}
	}
	close IN2;

##########initialize cv########################## 
	my $icv=&initialize_cv($dtime,$sam,@y);

	foreach my $chrnum(@chrname){
		if($count[$chrnum-1]==0){
			print OUT "$sam\tchr$chrnum\traw_number:$count[$chrnum-1]\n";
			next;
		}else{
			if($nnum[$chrnum-1] == 0){
				print OUT "$sam\tchr$chrnum\tfilterd_number:$nnum[$chrnum-1]\n";
				next;
			}
		}
		@yeeline=();
		@sumy=();
		@sumy2=();
		@em=();
		@ey=();
		my ($tsumy,$tsumy2);
		foreach my $k(1..$nnum[$chrnum-1]){
			$tsumy+=$y[$chrnum-1][$k-1];
			$tsumy2+=$y[$chrnum-1][$k-1]**2;
			$sumy[$k-1]=$tsumy;
			$sumy2[$k-1]=$tsumy2;
		}

##################################################################################
		my @k_error=&emin($sam,$dtime,$icv,1,$nnum[$chrnum-1],$chrnum,@{$y[$chrnum-1]});
		print OUT "$sam\tchr$chrnum\t$count[$chrnum-1]\t$nnum[$chrnum-1]\temin($k_error[0],1,$nnum[$chrnum-1],\@y)\t$k_error[1]\t$icv\t$k_error[2]\n";

		my $tmpy=$yeeline[0]; 
		my $tmpn;
		my $tmpj=0;
		for(my $i=0;$i<$nnum[$chrnum-1];$i++){
			if($i>=$tmpj){
				$tmpn=0;
				for(my $j=$i;$j<$nnum[$chrnum-1];$j++){
					if($tmpy==$yeeline[$j]){
						$tmpn++;
						$tmpj=$j+1;
					}else{
						$tmpy=$yeeline[$j];
						$tmpj=$j;
						last;
						}
				}
			}
			print OUT3 "$chrnum\t$x[$chrnum-1][$i]\t$y[$chrnum-1][$i]\t$yeeline[$i]\t$tmpn\n";
		}
	}
	close OUT3;
}

sub emin{
	my $name=shift @_;
	my $dtime=shift @_;
	my $cv=shift @_;
	my $start=shift @_;
	my $end=shift @_;
	my $chrnum=shift @_;
	my @ym=@_;
	&tpoint(@ym);
	my @rawp=($start..$end);
	my @points;
	my $cv_before=$cv;
	my $k=1;
	my ($emin,$tindex,$cve,$meane,$sd);
	while(1){
		if($k==1){
			@points=($start,$end);
			splice @rawp,($end-1),1;
			$emin=$em[$start-1][$end-1];
			$emin[$k]=$emin;
		}else{
			my @temp;
			my $ni=0;
			foreach my $i(@rawp){
				my @tpoints=@points;
				push @tpoints,$i,$i+1;
				@tpoints=sort {$a <=> $b} @tpoints;
				my $ntpoints=@tpoints;
				my $j=0;
				while($j<($ntpoints-1)){
					$temp[$ni]+=$em[$tpoints[$j]-1][$tpoints[$j+1]-1];
					$j+=2;
				}
				$ni++;
			}
			$emin=$temp[0];
			$tindex=0;  
			foreach my $i(1..($ni-1)){
				if($emin>$temp[$i]){
					$tindex=$i;
					$emin=$temp[$i];
				}
			}
			$emin[$k]=$emin;
			my $point=splice @rawp,$tindex,1;
			push @points,$point,$point+1;
			@points=sort {$a <=> $b} @points;
#########test the prior selected breakpoints is more significant than the new selected breakpoint#################
			my $tp=0;
			my $tntest=@points;
			my $numtest=$tntest;
			$tntest-=1;
			my $testemin=$emin[$k-1]; 
			for(my $l=1;$l<$tntest;$l+=2){
				my @testp=@points;
				my @trawp=@rawp;
				push @trawp,$testp[$l];
				@trawp=sort {$a<=>$b} @trawp;
				splice @testp,$l,2;

				my $m=0;
				my $ntest=@testp;
				while($m<($ntest-1)){
					$temp[$tp]+=$em[$testp[$m]-1][$testp[$m+1]-1];##
					$m+=2;
				}
				if($temp[$tp]>=$testemin){##
					next;
				}else{
					@points=@testp;
					@rawp=@trawp;
					$tntest-=2;
					$k--;
					$emin=$temp[$tp];##
					$emin[$k]=$emin;##
					$l-=1;
				}
				$tp++;
			} 
#################################################################################
		}
		my $j=0;
		my $npoints=@points;
		while($j<($npoints-1)){
			if($points[$j]==$points[$j+1]){
				$yeeline[$points[$j]-1]=$ym[$points[$j]-1];
			}else{
				foreach($points[$j]..$points[$j+1]){
					$yeeline[$_-1]=$ey[$points[$j]-1][$points[$j+1]-1];
				}
			}
			$j+=2; 
		}
		my @err;
		foreach($start..$end){
			$err[$_-1]=$ym[$_-1]-$yeeline[$_-1]+1;#######make the error values same with the observed values 
		}
		my ($sume,$sume2); 
		foreach($start..$end){
			$sume+=$err[$_-1];
			$sume2+=($err[$_-1]**2);
		}
		my $nume=@ym;
		$meane=$sume/$nume if $nume>0;
		my $vare=($sume2-$nume*($meane**2))/($nume-1) if $nume-1>0;
		$vare=0 if $nume-1==0;
		$cve=$vare**0.5/$meane if $meane!=0;
		$sd=$vare**0.5;
		print OUTCV "$chrnum\t$k\t$cve\t$meane\t$sd\n";

	################################fitting cutoff####################################################
		my $lambda=$hash_n{$name};
		if(defined($lambda)){
			if($lambda>1.1){
				$lambda=1.1;
			}elsif($lambda<0.9){
					$lambda=0.9;
			}
		}

		if($k==1){
			if($cve<=$lambda*$cv){
				last;
			}else{
				$k++;
				$cv_before=$cve;
			}
		}else{
			if(($cv_before**2-$cve**2)*($nume-1)<($proportion/2)**2*3 or $cve<=$lambda*$cv){
				last;
			}else{
				$k++;
				$cv_before=$cve;
			}
		}
	}
################################################################################################
	my @rt=($k,$emin,$cve);
	return @rt;
}

sub tpoint{
	my @y=@_;
	my $n=@y;
	my ($yee,$emin1);
	for(my $start=1;$start<=$n;$start++){
		for(my $end=$start;$end<=$n;$end++){
			if($start==$end){
				$em[$start-1][$start-1]=0;
				$ey[$start-1][$start-1]=$y[$start-1];
			}else{
				my $nnum1=$end-$start+1;
				if($start==1){
					$yee=$sumy[$end-1]/$nnum1 if $nnum1>0;
					 $emin1=$sumy2[$end-1]-$nnum1*$yee**2;
				}else{
					$yee=($sumy[$end-1]-$sumy[$start-2])/$nnum1 if $nnum1>0;
					$emin1=($sumy2[$end-1]-$sumy2[$start-2])-$nnum1*$yee**2;
				}
				$em[$start-1][$end-1]=$emin1;
				$ey[$start-1][$end-1]=$yee;
			}
		}
	}
}