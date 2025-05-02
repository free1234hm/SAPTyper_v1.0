
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
#MaxQuant.auto.pl <search_par> <pParse_par> <seq_db> <sample_info> <add_con> <out_dir>

#污染库的地址：.\data\contaminate.sequence.txt 
$con_db="$program_dir".'\data\contaminate.sequence.txt';
print "Con_db: $con_db\n";
$ref_search_par=$ARGV[0];
$seq_db=$ARGV[1];
$sample_info=$ARGV[2];
$add_con=$ARGV[3];  #1是要加，0是不加
$out_dir=$ARGV[4];

#if ($add_con==1)
if (0)  #Maxquant的污染数据不用独立处理数据库，因此此处可以不用管
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
	$sample_name[$sample_num]=$arr1[1];
	$sample_fraction[$sample_num]=$arr1[2];
	$sample_num++;
}
close IN_sample;

$search_par=$out_dir.'\mqpar.project.xml';
print "mqpar file: $search_par\n";
open (IN_search, "$ref_search_par")||die;
open (OUT_search, ">$search_par")||die;
$omit_string_line=0;
while ($in=<IN_search>)
{
	chomp $in;
	if ($in=~/^[ ]+<fastaFilePath>/)
	{
		$in=~s/^([ ]+<fastaFilePath>).*?</$1$seq_db</;
	}
	elsif ($in=~/^[ ]+<includeContaminants>/)
	{
		if ($add_con==1)
		{
			$in=~s/^([ ]+<includeContaminants>).*?</$1True</;
		}
		else
		{
			$in=~s/^([ ]+<includeContaminants>).*?</$1False</;
		}
	}
	elsif ($in=~/^[ ]+<fixedCombinedFolder>/)
	{
		$in=~s/^([ ]+<fixedCombinedFolder>).*?</$1$search_outdir</;
		#$in="outputpath=$search_outdir";
	}
#	elsif ($in=~/^outputname=/)
#	{
#		$in="outputname=pFindTask";
#	}
	elsif ($in=~/^[ ]+<filePaths>/)
	{
		for (my $i=0; $i<$sample_num; $i++)
		{
			$in.="\n      <string>$samlpe_arr[$i]".'</string>';
		}
		#$omit_string_line=1;
	}
	elsif ($in=~/^[ ]+<experiments>/)
	{
		for (my $i=0; $i<$sample_num; $i++)
		{
			$in.="\n      <string>$sample_name[$i]".'</string>';
		}
		#$omit_string_line=1;
	}
	elsif ($in=~/^[ ]+<fractions>/)
	{
		for (my $i=0; $i<$sample_num; $i++)
		{
			$in.="\n      <short>$sample_fraction[$i]".'</short>';
		}
		#$omit_string_line=1;
	}
	elsif ($in=~/^[ ]+<ptms>/)
	{
		for (my $i=0; $i<$sample_num; $i++)
		{
			$in.="\n      ".'<boolean>False</boolean>';
		}
		#$omit_string_line=1;
	}
	elsif ($in=~/^[ ]+<paramGroupIndices>/)
	{
		for (my $i=0; $i<$sample_num; $i++)
		{
			$in.="\n      ".'<int>0</int>';
		}
		#$omit_string_line=1;
	}
	elsif ($in=~/^[ ]+<referenceChannel>/)
	{
		for (my $i=0; $i<$sample_num; $i++)
		{
			$in.="\n      ".'<string></string>';
		}
		#$omit_string_line=1;
	}
#	elsif ($in=~/^msmspath1=/)
#	{
#		$in="";
#		for ($i=0; $i<$sample_num; $i++)
#		{
#			$j=$i+1;
#			$msfile=$samlpe_arr[$i];
#			$msfile=~s/\.raw$/_HCDFT\.pf2/;
#			$in.="msmspath$j=$msfile\n";
#		}
#		$in=~s/\n$//;
#	}
#	elsif ($in=~/^msmspath[0-9]+=/)
#	{
#		next;
#	}
	if (($in=~/^[ ]+<string>/) or ($in=~/^[ ]+<short>/) or ($in=~/^[ ]+<boolean>/))
	{
		if ($omit_string_line==0)
		{
			print OUT_search "$in\n";
		}
	}
	else
	{
		$omit_string_line=0;
		print OUT_search "$in\n";
		if (($in=~/^[ ]+<filePaths>/) or ($in=~/^[ ]+<experiments>/) or ($in=~/^[ ]+<fractions>/) or ($in=~/^[ ]+<ptms>/) or ($in=~/^[ ]+<paramGroupIndices>/) or ($in=~/^[ ]+<referenceChannel>/))
		{
			$omit_string_line=1;
		}
	}
	
}
close IN_search;
close OUT_search;

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

$pfind_dir=$dep_program_hash{'MaxQuant'};
chdir ($pfind_dir);
print "path=$pfind_dir\n";

#Searcher.exe d:\pFindWorkspace\pFindTask2\autorun_test\pFind.cfg
$cmd="MaxQuantCmd.exe $search_par";
print "Maxquant dir: $pfind_dir\nCMD=$cmd\n";
system ("$cmd");  #debug

#chdir ($search_outdir);
#$search_result=$out_dir.'\result\pFind-Filtered.spectra';
#$search_result=~s/\\\\/\\/g;
#$pfind_trans="$program_dir".'\pfind.trans.result.pl ';
#$cmd="$pfind_trans $search_result";
#print "CMD: $cmd\n";
#system ($cmd);
