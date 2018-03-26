#!/bin/perl -w
use strict;
use Cwd;
use FindBin qw($Bin);
use Getopt::Long;

=head1 NAME:

	Program:	somatic_cnv_detect_tool.pl
	Version:	0.1
	Author:		qianzhaoyang(qianzhaoyang@genomics.cn)
				wangxiaofeng(wangxiaofeng@genomics.cn)
				shichang(shichang@genomics.cn)
	Modified Date:	2014-05-19
					2017-11-10
					2018-01-10
	Description:	cfDNA sequencing cnv analyze pipeline


=head1 Usage:

	somatic_cnv_detect_tool.pl   [options]

	Options:
		-i|--input         <STR>    (required)absolute path for sample file list(sort.markdup.bam),col[0]=name col[1]=absolute path  
		-o|--out           <STR>    (required)output directory 
		-a|--alpha         <FLOAT>  (required)alpha_confidence_level 
		-d|--debin         <INT>    (required) debin length for first GC normalization(bp) 
		-c|--control       <STR>    absolute path for control list(sort.markdup.bam),col[0]=control col[1]=absolute path, default F
		-b|--bin           <INT>    bin length(bp),default 1000000
		-p|--pair          <STR>    pair end sequencing: T or F,default F 
		-m|--maxq          <INT>    filter_bam_maxq,default 60
		-l|--length        <INT>    windows length to calculate the GC-content distribution (bp),default 170 for cfDNA data
		-t|--threshold     <FLOAT>  threshold to detect CNV,default 0.01(0~1)
		-r|--readcount     <INT>    readcount to calculate lambda for terminate the loop,default 20000000
		   --proportion    <FLOAT>  required chimeric proportion for terminate the loop,default 0.03
		   --cnvsize       <INT>    required CNV size more than n times bin length for terminate the loop,default 3
		   --off_target    <STR>    off_target_data:T or F,default F
		   --bed           <STR>    a file of target.bed
		   --merge_con     <STR>    a file merge_contrl.cnr.stac.file.list,default F
		   --remove        <STR>    remove tmp_data after detect process:T or F,default T

=head1 Example:

    somatic_cnv_detect_tool.pl -i ~/a/b/c/bam.list -a 0.01 -d 100000 -o result 

=cut



my ($bamFilelist,$outdir,$alpha,$debin_length);
my ($confilelist,$bin_length,$PE,$maxq,$gc_length,$threshold,$readcount,$proportion,$cnvsize,$off_target,$bed,$merge_con_path,$remove);
my $path = $Bin;

GetOptions(
	"i|input=s" => \$bamFilelist,
	"c|control=s" => \$confilelist,
	"o|out=s" => \$outdir,
	"a|alpha=f" => \$alpha,
	"d|debin=i" => \$debin_length,
	"b|bin=i" => \$bin_length,
	"p|pair=s"  => \$PE,
	"m|maxq=s" => \$maxq,
	"l|length=i" => \$gc_length,
	"t|threshold=f" => \$threshold,
	"r|readcount=i" => \$readcount,
	"proportion=f" => \$proportion,
	"cnvsize=i" => \$cnvsize,
	"off_target=s" => \$off_target,
	"bed=s" => \$bed,
	"merge_con=s" => \$merge_con_path,
	"remove=s" => \$remove,
);
$confilelist="F" unless defined($confilelist);
$bin_length=1000000 unless defined($bin_length);
$PE="F" unless defined($PE);
$maxq=60 unless defined($maxq);
$gc_length=170 unless defined($gc_length);
$threshold=0.01 unless defined($threshold);
$readcount=20000000 unless defined($readcount);
$proportion=0.03 unless defined($proportion);
$cnvsize=3 unless defined($cnvsize);
$off_target="F" unless defined($off_target);
$merge_con_path="F" unless defined($merge_con_path);
$remove="T" unless defined($remove);
my $date=`date +%y%m%d`;
chomp $date;


die `pod2text $0` unless($bamFilelist && $alpha && $debin_length && $outdir );
die ".bed file for off_target data must be specified or does not exist! $0" if ($off_target ne "F" && !defined($bed));

#####default data files#####
my $chr_length="$path/data/Chr_Length.txt";
my $info_chrY_male="$path/data/chrY_interval_unique_q30_100kmer_1kb_1M_nofemale.info ";
my $gc_stac_ref="$path/data/hg19.$gc_length.no_sex.gc.stac";
my $gc_ref_1="$path/data/hg19.$gc_length.1.gc.gz";
my $gc_ref_2="$path/data/hg19.$bin_length.2.gc";

if ($gc_length !=170 && ! -s ($gc_ref_1 || $gc_stac_ref)){
		die "ref gc file for -l $gc_length must be specified or does not exist! $0";
	}

die "ref gc file for -b $bin_length must be specified or does not exist! $0" if ($gc_length !=1000000 && ! -s $gc_ref_2);
#####scripts#####
my $TargetOff="$path/scripts/1.1_TargetOff.pl";
my $ExtractInfo="$path/scripts/1.2_ExtractInfo.pl";
my $Stac_read_gc_onref="$path/scripts/2.1_Stac_read_gc_onref.pl";
my $Stac_reads_gc_cnr="$path/scripts/2.2_Stac_reads_gc_cnr.pl";
my $Correct="$path/scripts/3_Correct_cnv_median_list_ref.pl";
my $Normalize="$path/scripts/4_Normalize_CNV_median_list.pl";
my $Stat="$path/scripts/5_Stat_CNV_list_M.pl";
my $Regc="$path/scripts/6_Re_GC_Normalize.pl";
my $Fit="$path/scripts/7_Segmentation_fit.pl";
my $Detect="$path/scripts/8_Detect_plot_cnv_total_chr.pl";
my $Rplot="$path/scripts/9_Rplot2.pl";
my $Rplot1="$path/scripts/2_Rplot.R";
my $Mergedata="$path/scripts/3_Merge.R";
my $Meandata="$path/scripts/3_Mean.R";
my $NormalizeR="$path/scripts/4_Normalize.R";
my $Lm2="$path/scripts/6_Lm-2.R";
my $Normalize2R="$path/scripts/7_Normalize.R";
my $Detectplot="$path/scripts/8_Detect_plot_cnv.R";
my $Rplotr="$path/scripts/9_Rplot2.R";


our %config;
my $config = "$path/scdt_cfg.pl";
require $config;


`mkdir $outdir`;
chdir ("$outdir");
my $dir=cwd;
`mkdir tmp_data`;
`mkdir Figure`;
my $datadir="$dir/tmp_data";
my $figdir="$dir/Figure";
`mkdir sh`;

#########list the control cnr stac file########
my $control_merge_cnr_stac;
my $cnr_stac_control_list;
if($merge_con_path eq 'F'){
	open OUTC,">control_ref_$debin_length.merge_cnr_stac.list"or die $!;
}else{
	open CCI,"$merge_con_path"or die $!;
	$control_merge_cnr_stac=<CCI>;
	chomp($control_merge_cnr_stac);
	close CCI;
}


open SH, "> $dir/sh/run_scdt.sh" or die $!;
print SH "#!/bin/bash\n\n";
print SH "log=$dir/scdt.log\ndate '+%y-%m-%d %H:%M:%S\tanalysis begin' > \$log\n\n";
print SH "###########################################################\n";
print SH "##########  1.get_sample_info_rgc_cnr_stac  ###############\n";
print SH "###########################################################\n";
my $outname;
my $i=0;
open IN,"$bamFilelist"or die $!;
##########1.get_info_rgc_cnr_stac_##########
while(my $bamfile=<IN>){
	chomp($bamfile);
	$i++;
	my @array=split(/\t/,$bamfile);
	my $file=$array[-1];
	$outname=$array[0];
	$cnr_stac_control_list="F";
	open OUT,">$dir/sh/1.get_sample_info_rgc_cnr_stac_$outname.$debin_length.sh"or die $!;
	print OUT "#!/bin/sh\n";

#### bwa mem and get info ####
	my $sortbam=$file;

	if($off_target eq "T"){ 
		print OUT $config{"perl"}," $TargetOff $sortbam $bed $maxq $PE >$datadir/$outname.no_bed.info\n";
	}else{
		print OUT $config{"perl"}," $ExtractInfo $sortbam $maxq $PE >$datadir/$outname.info\n";
	}

	print OUT "sort -n -k 1.4 -k 2 $datadir/$outname.*info >$datadir/$outname.sort.info\n";


	if($off_target eq "T"){
		print OUT "rm -rf $datadir/$outname.no_bed.info\n";
	}else{
		print OUT "rm -rf $datadir/$outname.info\n";
	}
	print OUT $config{"perl"}," $Stac_read_gc_onref $datadir/$outname.sort.info $gc_ref_1 $chr_length $datadir/$outname\n";
	print OUT $config{"perl"}," $Stac_reads_gc_cnr $datadir/$outname.read_gc $datadir/$outname.read_gc.stac $gc_stac_ref $info_chrY_male $debin_length $datadir/$outname.$debin_length $figdir\n";
	print OUT $config{"Rscript"}," $Rplot1 $datadir/$outname.$debin_length.mer $figdir/$outname $outname\n";
	close OUT;

	`chmod +x sh/1.get_sample_info_rgc_cnr_stac_$outname.$debin_length.sh`;
	print SH "nohup ./1.get_sample_info_rgc_cnr_stac_$outname.$debin_length.sh >1.get_sample_info_rgc_cnr_stac_$outname.$debin_length.out 2>1.get_sample_info_rgc_cnr_stac_$outname.$debin_length.err&\\\n";

}
close IN;

if($merge_con_path eq 'F' and $cnr_stac_control_list ne 'F'){
	close CC;
}
print SH "wait\ndate '+%y-%m-%d %H:%M:%S\t1.get_info_rgc_cnr_stac done' >> \$log\n\n";

if ($confilelist ne "F"){
	print SH "############################################################\n";
	print SH "#########  1.5.get_control_info_rgc_cnr_stac  ##############\n";
	print SH "############################################################\n";
	`mkdir -p $datadir/control`;
	my $condir="$datadir/control";
	my $outname;
	my $j=0;
	open CON,"$confilelist"or die $!;
		while(my $bamfile=<CON>){
		chomp($bamfile);
		$j++;
		my @array=split(/\t/,$bamfile);
		my $file=$array[-1];
		if($j==1){
			open CC,">$condir/cnr_stac.control.$debin_length.list"or die $!;
		}
		$outname=(split /\//,$file)[-1];
		$outname=(split /\./,$outname)[0];
		print CC "$condir/$outname.$debin_length.cnr.stac\n";
		$cnr_stac_control_list="$condir/cnr_stac.control.$debin_length.list";

		open OUT,">$dir/sh/1.5.get_control_info_rgc_cnr_stac_$outname.$debin_length.sh"or die $!;
		print OUT "#!/bin/sh\n";

	#### bwa mem and get info ####
		my $sortbam=$file;

		if($off_target eq "T"){ 
			print OUT $config{"perl"}," $TargetOff $sortbam $bed $maxq $PE >$condir/$outname.no_bed.info\n";
		}else{
			print OUT $config{"perl"}," $ExtractInfo $sortbam $maxq $PE >$condir/$outname.info\n";
		}

		print OUT "sort -n -k 1.4 -k 2 $condir/$outname.*info >$condir/$outname.sort.info\n";


		if($off_target eq "T"){
			print OUT "rm -rf $condir/$outname.no_bed.info\n";
		}else{
			print OUT "rm -rf $condir/$outname.info\n";
		}
		print OUT $config{"perl"}," $Stac_read_gc_onref $condir/$outname.sort.info $gc_ref_1 $chr_length $condir/$outname\n";
		print OUT $config{"perl"}," $Stac_reads_gc_cnr $condir/$outname.read_gc $condir/$outname.read_gc.stac $gc_stac_ref $info_chrY_male $debin_length $condir/$outname.$debin_length \n";
		print OUT $config{"Rscript"}," $Rplot1 $datadir/control/$outname.$debin_length.mer $figdir/$outname $outname\n";
		close OUT;

		`chmod +x sh/1.5.get_control_info_rgc_cnr_stac_$outname.$debin_length.sh`;
		print SH "nohup ./1.5.get_control_info_rgc_cnr_stac_$outname.$debin_length.sh >1.5.get_control_info_rgc_cnr_stac_$outname.$debin_length.out 2>1.5.get_control_info_rgc_cnr_stac_$outname.$debin_length.err&\\\n";

	}
	close CON;
	print SH "wait\ndate '+%y-%m-%d %H:%M:%S\t1.5.get_control_info_rgc_cnr_stac done' >> \$log\n\n";
}

	
	print SH "###############################################################\n";
	print SH "####2.Baseline construction and copy ratio calculation  #######\n";
	print SH "###############################################################\n\n";

#########2.get_ref_stac#######
open OUT,">$dir/sh/2.get_ref_stac.$debin_length.sh"or die $!;
print OUT "#!/bin/sh\n";
print OUT "ls $datadir/*$debin_length.cnr.stac >$datadir/cnr_stac.$debin_length.list\n";


##use the control sample to produce the filter bin( the reads variate big):conditon the median and sd values
if($merge_con_path eq 'F'){
	if($cnr_stac_control_list ne 'F'){  
		$control_merge_cnr_stac="$dir/$date.merge.control.$debin_length.cnr.stac";
		print OUT $config{"Rscript"}," $Mergedata $cnr_stac_control_list $control_merge_cnr_stac\n";
	}else{
		$control_merge_cnr_stac="$dir/$date.merge.all.$debin_length.cnr.stac";
		print OUT $config{"Rscript"}," $Mergedata $datadir/cnr_stac.$debin_length.list $control_merge_cnr_stac\n";
	}
}
print OUT $config{"Rscript"}," $Meandata $control_merge_cnr_stac $dir/control.cnr.stac\n";
print OUT $config{"perl"}," $Correct $datadir/cnr_stac.$debin_length.list $dir/control.cnr.stac\n ";
#print OUT $config{"perl"}," $Correct $merge_con_path $datadir/cnr_stac.$debin_length.list $cnr_stac_control_list $control_merge_cnr_stac $dir/control.cnr.stac\n ";
##get: control_mean.cnr.stac name.CNV

print OUT "ls $datadir/*.CNV >$datadir/CNV.$debin_length.list\n";
print OUT $config{"perl"}," $Normalize $datadir/CNV.$debin_length.list $NormalizeR ",$config{"Rscript"},"\n";##get: name.CNV.nr
close OUT;
`chmod +x sh/2.get_ref_stac.$debin_length.sh`;
print SH "nohup ./2.get_ref_stac.$debin_length.sh >2.get_ref_stac.$debin_length.out 2>2.get_ref_stac.$debin_length.err&\\\n";
print SH "wait\ndate '+%y-%m-%d %H:%M:%S\t2.get_ref_stac done' >> \$log\n\n";
	

	print SH "###############################################################\n";
	print SH "########################  3.Segment fit  ######################\n";
	print SH "###############################################################\n";
open  OUT,">$dir/sh/3.Rplot_CNV_$debin_length.sh"or die $!;
print OUT "#!/bin/sh\n";
# print OUT "ls $datadir/*.read_gc.stac >$datadir/read_gc.stac.list\nls $datadir/*.mer >$datadir/mer_list\n";
# my $dirFigure="$dir/Figure";
# print OUT "mkdir $dirFigure\n";
# print OUT $config{"perl"}," $path/scripts/Rplot_cnr_gc_ref.pl $datadir/mer_list $i $dirFigure/RCr_GC.tiff\n";
# print OUT $config{"perl"}," $path/scripts/Rplot_sample.pl $datadir/cnr_stac.$debin_length.list $dirFigure/all_reads_cnr.$debin_length.pdf\n";

print OUT "ls $datadir/*.CNV.nr >$datadir/CNV_nr.$debin_length.list\n";	
print OUT $config{"perl"}," $Stat $datadir/CNV_nr.$debin_length.list $debin_length $bin_length\n";#result:name.nr.correct.CNV.M
print OUT "ls $datadir/*.M >$datadir/CNV_nr.M.list\n";
print OUT $config{"perl"}," $Normalize $datadir/CNV_nr.M.list $NormalizeR ",$config{"Rscript"},"\n";
print OUT "ls $datadir/*.M.nr >$datadir/CNV_nr.M.nr.list\n";

print OUT $config{"perl"}," $Regc $gc_ref_2 $datadir/CNV_nr.M.nr.list $Lm2 ",$config{"Rscript"}," $datadir\n";
print OUT "ls $datadir/*.reNGC >$datadir/reNGC.list\n";
# print OUT $config{"perl"}," $path/scripts/Rplot.pl $datadir/reNGC.list $dirFigure/all_ratios.$bin_length.pdf\n";
print OUT "wc -l $datadir/*.sort.info >$datadir/stat_reads_number.xls\n";
print OUT $config{"perl"}," $Fit $datadir/reNGC.list $alpha $cnvsize $readcount $proportion $datadir/stat_reads_number.xls 250000000/$bin_length $datadir $Normalize2R ",$config{"Rscript"},"\n";

print OUT "ls $datadir/*.O.nr$alpha.E.ratio >$datadir/y.ye.list\n";
print OUT $config{"perl"}," $Detect $datadir/y.ye.list $bin_length $alpha $threshold 250 $dir/detect_result $cnvsize*$bin_length $date $Detectplot ",$config{"Rscript"},"\n";
print OUT "ls $dir/detect_result_$date/*.cnv_line_p_value.xls > $dir/detect_result_$date/p_value.list\n";
print OUT $config{"perl"}," $Rplot $dir/detect_result_$date/p_value.list $figdir $Rplotr ",$config{"Rscript"},"\n";

close OUT;
`chmod +x sh/3.Rplot_CNV_$debin_length.sh`;
print SH "nohup ./3.Rplot_CNV_$debin_length.sh >3.Rplot_CNV_$debin_length.out 2>3.Rplot_CNV_$debin_length.err&\\\n";
print SH "wait\ndate '+%y-%m-%d %H:%M:%S\t3.Segment fit done' >> \$log\n\n";
if ($remove eq "T"){
	print SH "if [-s $dir/detect_result_$date/detect_cnv_line_Pvalue.$date.xls ]; then\n rm -rf $datadir\n fi \n";
}
if($merge_con_path eq 'F'){
	if(defined($control_merge_cnr_stac)){
		print OUTC "$control_merge_cnr_stac\n";
	}
	close OUTC;
}
