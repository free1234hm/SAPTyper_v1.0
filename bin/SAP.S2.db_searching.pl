
#####################
#Project: SAP Project
#Program: 
#Version: 1.0
#Created on:
#Author: wusf
#####################

#实现Maxquant和pFind的自动搜索。由于不涉及到比较的问题，因此直接做鉴定的参数

#参数文件还是由相应的软件编辑，只重置Raw文件和搜索数据库。

$search_engine=1; #1-maxquant, 2-pFind
$ori_search_engine_par='r:\work\SAPproject\pipeline\mqpar.xml';
$maxquant_CMD='p:\Software\proteomics\search_engine\MaxQuant\MaxQuant.exe';
$pfind_CMD='d:\pFindStudio\pFind3\bin\Searcher.exe';

$module_par='SAP_ana.M1.par.txt';
open (IN_par, "$module_par")||die;
####从M3中拷贝来的程序---start
#$raw_i=0;
$exp_i=0;
$in_par=<IN_par>;
while ($in_par)
{
	chomp $in_par;
	#print "$in_par==\n";
	if ($in_par=~/^Protein database:\t([^\t]+)/)
	{
		$pro_db=$1;
	}
	#raw文件的合并到experimental design里面了，防止出现不一致的情况
	#elsif ($in_par=~/^MS Raw file list:/)
	#{
	#	#$pfind_result_marker=1;
	#	while ($in_par=<IN_par>)
	#	{
	#		chomp $in_par;
	#		if ($in_par=~/^\t(.+)/)
	#		{
	#			$raw_file_arr[$raw_i]=$1;
	#			$raw_i++;
	#		}
	#		else
	#		{
	#			next;
	#		}
	#	}
	#	next;
	#}
	elsif ($in_par=~/^Experimental design /)
	{
		#$pfind_result_marker=1;
		while ($in_par=<IN_par>)
		{
			chomp $in_par;
			if ($in_par=~/^\t(.+)/)
			{
				$exp_info=$1;
				undef @arr1;
				@arr1=split("\t", $exp_info);
				$group[$exp_i]=$arr1[0];
				$sample[$exp_i]=$arr1[1];
				$fraction[$exp_i]=$arr1[2];
				$raw_file_arr[$exp_i]=$arr1[3];
				$exp_i++;
			}
			else
			{
				next;
			}
		}
		next;
	}
	$in_par=<IN_par>;
}
close IN_par;
####从M3中拷贝来的程序---end

if ($search_engine==1)
{
	open (IN_MQ_par, "$ori_search_engine_par")||die;
	open (OUT_MQ_par, ">project.mqpar.xml")||die;

	while ($in=<IN_MQ_par>)
	{
		chomp $in;
		if ($in=~s/^         <fastaFilePath>.+?</         <fastaFilePath>$pro_db</)
		{
			print OUT_MQ_par "$in\n";
			next;
		}
		elsif ($in=~/^   <filePaths>/)
		{
			print OUT_MQ_par "$in\n";
			for ($i=0; $i<$exp_i; $i++)
			{
				print OUT_MQ_par "      <string>$raw_file_arr[$i]";
				print OUT_MQ_par "</string>";
				print OUT_MQ_par "\n";
			}
			while ($in=<IN_MQ_par>)
			{
				last if ($in!~/^      <string>/);
			}
		}
		elsif ($in=~/^   <experiments>/)
		{
			print OUT_MQ_par "$in\n";
			for ($i=0; $i<$exp_i; $i++)
			{
				print OUT_MQ_par "      <string>$sample[$i]";
				print OUT_MQ_par "</string>";
				print OUT_MQ_par "\n";
			}
			while ($in=<IN_MQ_par>)
			{
				last if ($in!~/^      <string>/);
			}
		}
		elsif ($in=~/^   <fractions>/)
		{
			print OUT_MQ_par "$in\n";
			for ($i=0; $i<$exp_i; $i++)
			{
				print OUT_MQ_par "      <short>$fraction[$i]";
				print OUT_MQ_par "</short>";
				print OUT_MQ_par "\n";
			}
			while ($in=<IN_MQ_par>)
			{
				last if ($in!~/^      <short>/);
			}
		}
		print OUT_MQ_par "$in\n";
	}
	close IN_MQ_par;
	close OUT_MQ_par;
	#Database searching using maxquant
	$CMD="$maxquant_CMD project.mqpar.xml";
	system ("$CMD");
}
if ($search_engine==2)
{
	

}

