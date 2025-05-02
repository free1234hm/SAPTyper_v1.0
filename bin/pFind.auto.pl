
#####################
#Project: 
#Program: 
#Version: 1.0
#Created on:
#Author: wusf
#####################

use File::Basename;
$program_dir=dirname(__FILE__);

#命令执行方式和输入文件：
#pFind.auto.pl <search_par> <pParse_par> <seq_db> <sample_info> <add_con> <out_dir>
#<search_par>: 
#<pParse_par>:
#<seq_db>: 
#<sample_info>: 从proteome data来
#<add_con>： 
#<out_dir>： 从project_dir来，后面跟db_search子目录

#污染库的地址：.\data\contaminate.sequence.txt 
$con_db="$program_dir".'\data\contaminate.sequence.txt';
print "Con_db: $con_db\n";
$ref_search_par=$ARGV[0];
$ref_pparse_par=$ARGV[1];
$seq_db=$ARGV[2];
$sample_info=$ARGV[3];
$add_con=$ARGV[4];  #1是要加，0是不加
$out_dir=$ARGV[5];


if ($add_con==1)
{
	$fasta_file=$out_dir.'\seq.db.con.fasta';
	print "CMD=type $seq_db $con_db > $fasta_file";
	system ("type $seq_db $con_db > $fasta_file");
}
else
{
	$fasta_file=$out_dir.'\seq.db.fasta';
	system ("type $seq_db > $fasta_file");
}

$search_outdir=$out_dir.'\result\\';
$search_outdir=~s/\\\\/\\/g;
print "CMD=md $search_outdir\n";
system ("md $search_outdir");

open (IN_sample, "$sample_info")||die;
$sample_num=0;
while ($in=<IN_sample>)
{
	chomp $in;
	next if ($in=~/^Person	/);
	undef @arr1;
	@arr1=split("\t", $in);
	$samlpe_arr[$sample_num]=$arr1[3];
	$sample_num++;
}
close IN_sample;

$search_par=$out_dir.'\pFind.cfg';
open (IN_search, "$ref_search_par")||die;
open (OUT_search, ">$search_par")||die;
while ($in=<IN_search>)
{
	chomp $in;
	if ($in=~/^fastapath=/)
	{
		$in="fastapath=$fasta_file";
	}
	elsif ($in=~/^outputpath=/)
	{
		$in="outputpath=$search_outdir";
	}
	elsif ($in=~/^outputname=/)
	{
		$in="outputname=pFindTask";
	}
	elsif ($in=~/^msmsnum=/)
	{
		$in="msmsnum=$sample_num";
	}
	elsif ($in=~/^msmspath1=/)
	{
		$in="";
		for ($i=0; $i<$sample_num; $i++)
		{
			$j=$i+1;
			$msfile=$samlpe_arr[$i];
			while ($msfile=~s/ $//g) {}
			$msfile=~s/\.raw$/_HCDFT\.pf2/;
			$in.="msmspath$j=$msfile\n";
		}
		$in=~s/\n$//;
	}
	elsif ($in=~/^msmspath[0-9]+=/)
	{
		next;
	}
	print OUT_search "$in\n";
}
close IN_search;
close OUT_search;

$pparse_par=$out_dir.'\pParse.cfg';
open (IN_pparse, "$ref_pparse_par")||die;
open (OUT_pparse, ">$pparse_par")||die;
while ($in=<IN_pparse>)
{
	chomp $in;
	if ($in=~/^datanum=/)
	{
		$in="datanum=$sample_num";
	}
	elsif ($in=~/^datapath1=/)
	{
		$in="";
		for ($i=0; $i<$sample_num; $i++)
		{
			$j=$i+1;
			$in.="datapath$j=$samlpe_arr[$i]\n";
		}
		$in=~s/\n$//;
	}
	elsif ($in=~/^datapath[0-9]+=/)
	{
		next;
	}
	print OUT_pparse "$in\n";
}
close IN_pparse;
close OUT_pparse;

$dep_program_dir="$program_dir".'\dep_program.dir';
open (IN, "$dep_program_dir")||die;
while ($in=<IN>)
{
	chomp $in;
	print "==$in\n";
	undef @arr1;
	@arr1=split("=", $in);
	$dep_program_hash{$arr1[0]}=$arr1[1];
}
close IN;

$pfind_dir=$dep_program_hash{'pFind'};
chdir ($pfind_dir);
print "path=$pfind_dir\n";
#pParse.exe d:\pFindWorkspace\pFindTask2\autorun_test\pParse.cfg
$cmd="pParse.exe $pparse_par";

system ("$cmd");  #debug

#Searcher.exe d:\pFindWorkspace\pFindTask2\autorun_test\pFind.cfg
$cmd="Searcher.exe $search_par";
print "pFind dir: $pfind_dir\nCMD=$cmd\n";
system ("$cmd");  #debug

chdir ($search_outdir);
$search_result=$out_dir.'\result\pFind-Filtered.spectra';
$search_result=~s/\\\\/\\/g;
$pfind_trans="$program_dir".'\pfind.trans.result.pl ';
$cmd="$pfind_trans $search_result";
print "CMD: $cmd\n";
system ($cmd);
