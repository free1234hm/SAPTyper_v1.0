
#########################
#Project: SAP discoverage
#Version: 1.0
#Created on:
#Author: wusf
#########################

#V2: 兼容不同的搜索引擎

#v4: 兼容DIA的搜库结果

$par_file=$ARGV[0];
my  ($sec,$min,$hour,$mday,$mon,$year) = (localtime)[0..5];
($sec,$min,$hour,$mday,$mon,$year) = (
    sprintf("%02d", $sec),
    sprintf("%02d", $min),
    sprintf("%02d", $hour),
    sprintf("%02d", $mday),
    sprintf("%02d", $mon + 1),
    $year + 1900
);
$datetime="$year$mon$mday.$hour.$min.$sec";
#$datetime="$year$mon$mday.$hour"."H.$min"."M.$sec"."S";
$outlog="SAP.M3.$datetime.log";
open (OUT_log, ">$outlog")||die;
print OUT_log "Start time: $datetime\n\n";
#close OUT_log;

#命令目录
use File::Basename;
$program_dir=dirname(__FILE__);
#print "$program_dir\n";
#$pipeline_path='r:\work\SAPproject\pipeline';
$cmd_path="$program_dir\\";

#$cmd_path='r:\work\SAPproject\pipeline\\';
$par_file='SAP.S3.data_ana.par' if (!($par_file));

open (IN_par, "$par_file")||die;

#$pfind_i=0;
#$in_par=<IN_par>;
$raw_i=0;
while ($in_par=<IN_par>)
{
	chomp $in_par;
	#print "$in_par==\n";
	if ($in_par=~/^SNP_annotation\t(.+?)\t/)
	{
		$snp_anno=$1;
	}
	elsif ($in_par=~/^SAP_annotation\t(.+?)\t/)
	{
		$sap_anno=$1;
		#print "===$sap_anno\n";
	}
	elsif ($in_par=~/^Search_engine\t(.+?)\t/)
	{
		$search_engine=$1;
		#print "===$sap_anno\n";
	}
	elsif ($in_par=~/^Proteome results:/)
	{
		#$pfind_result_marker=1;
		#print "$in_par\n";
		$in_par=<IN_par>;
		#print "$in_par\n";
		chomp $in_par;
		if ($in_par=~/^\t(.+)/)
		{
			$proteome_search_result_file=$1;
		}
	}
	elsif ($in_par=~/^Experimental design and raw files/)
	{
		#$pfind_result_marker=1;
		while ($in_par=<IN_par>)
		{
			chomp $in_par;
			last if ($in_par!~/^\t(.+)/);
			undef @arr1;
			@arr1=split("\t", $in_par);
			$raw_file=$arr1[4];
			$group=$arr1[1];
			$sample=$arr1[2];
			last if (!($group));
			$pro_raw_arr[$raw_i]="$raw_file";
			$raw_i++;
			$raw2group{$raw_file}=$group;
			$raw2sample{$raw_file}=$sample;
			$sample2raw{$sample}.=",$raw_file";
			$group2raw{$group}.=",$raw_file";
		}
	#	next;
	}
	elsif ($in_par=~/^Wes results /)
	{
		#$pfind_result_marker=1;
		while ($in_par=<IN_par>)
		{
			chomp $in_par;
			last if ($in_par!~/^\t(.+)/);
			undef @arr1;
			@arr1=split("\t", $in_par);
			$QC_file=$arr1[2];
			$group=$arr1[1];
			last if (!($group));
			if ($QC_file)
			{
				$QCfile2group{$QC_file}=$group;
				$group2QCfile{$group}=$QC_file;
			}
			else
			{
				last;
			}
		}
	#	next;
	}
	elsif ($in_par=~/^Protein_annotation\t(.+?)\t/)
	{
		$pro_anno=$1;
	}
	elsif ($in_par=~/^SNP2rsID\t(.+?)\t/)
	{
		$SNP2rsID_file=$1;
	}
	elsif ($in_par=~/^Gene annotation\t(.+?)\t/)
	{
		$gene_anno=$1;
	}
	#$in_par=<IN_par>;
}


#拆分搜索结果文件，并且读入拆分的结果
$cmd="$cmd_path"."sep.proteome.searching.result.v2.pl $proteome_search_result_file $search_engine";
print "CMD: $cmd\n";
system ("$cmd");
#读入拆分结果文件，是利用参数文件中的raw的列表寻找对应的拆分文件，如果有找不到的，报错并退出
$pfind_i=0;
for ($i=0; $i<@pro_raw_arr; $i++)
{
	$raw_file=$pro_raw_arr[$i];
	$group=$raw2group{$raw_file};
	$sample=$raw2sample{$raw_file};
	$raw_file=~s/^.+\\//;
	$raw_file=~s/\.raw$//;
	$pro_seq_result_file='proteome.'."$raw_file".'.sep.txt';
	$pro_seq_result_arr[$pfind_i]=$pro_seq_result_file;
	#$pfind_file_arr[$pfind_i]=$pro_seq_result_file;
	$pfind_i++;

	#add the following code for the sample and group marker
	$sample_exp="sample.$sample";
	$sample_dir{$sample_exp}=1;
	$group_exp="group.$group";
	$sample_dir{$group_exp}=1;

}
print "Read raw file name complete: $pfind_i raw files...\n";
#process each separated searching results
for ($i=0; $i<$pfind_i; $i++)
{
	$pro_seq_result=$pro_seq_result_arr[$i];
	$cmd="$cmd_path"."determine.var.ref.result.v0.7.pl $snp_anno $sap_anno $search_engine $pro_seq_result";
	print OUT_log "CMD: $cmd\n";
	print "CMD: $cmd\n";
	system ("$cmd");   #对鉴定到的突变进行分组
}

#process the results produced by above step, according to the group and sample
$cmd="$cmd_path"."create.comp.table.v5.pl $par_file";
system ("$cmd"); 

#Annotation:
#M0的<out file>中含有line到蛋白的对应关系
#protein annotation
#rs ID
#*.vcf.depth.QC.txt 文件的信息
$cmd="$cmd_path"."SAP.annotation.v4.pl $par_file";
system ("$cmd");



#Generate report:
$cmd="$cmd_path"."generate.report.v4.pl $par_file";
system ("$cmd");

####

($sec,$min,$hour,$mday,$mon,$year) = (localtime)[0..5];
($sec,$min,$hour,$mday,$mon,$year) = (
    sprintf("%02d", $sec),
    sprintf("%02d", $min),
    sprintf("%02d", $hour),
    sprintf("%02d", $mday),
    sprintf("%02d", $mon + 1),
    $year + 1900
);
$datetime="$year$mon$mday.$hour.$min.$sec";
#open (OUT_log, ">>$outlog")||die;
print OUT_log "\nEnd time: $datetime\n\n";
close OUT_log;
