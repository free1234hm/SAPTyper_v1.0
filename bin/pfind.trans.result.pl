
#####################
#Project: 
#Program: 
#Version: 1.0
#Created on:
#Author: wusf
#####################

#将pFind的直接搜库结果转化为导出的result.txt

$infile=$ARGV[0];
$outfile="result.txt";

open (IN, "$infile")||die;
open (OUT, ">$outfile")||die;

print OUT '#	Title	Charge	Sq	Mod_Sites	Score	Spectra Mass	Q value	Theory Sq Mass	Delta Mass	Delta Mass (PPM)	Specific Flag	Label Flag	Target_Decoy	Protein AC';
print OUT "\n";
$serial=1;
while ($in=<IN>)
{
	next if ($in=~/^File_Name	Scan_No/);
	chomp $in;
	undef @arr1;
	@arr1=split("\t", $in);
	$title=$arr1[0];
	$charge=$arr1[3];
	$sq=$arr1[5];
	$mod=$arr1[10];
	$score=$arr1[9];
	$spec_mass=$arr1[2];
	$Qvalue=$arr1[4];
	$Thero_mass=$arr1[6];
	$Delta_mass=$arr1[7];
	if ($spec_mass>0)
	{
		$Delta_mass_ppm=1000000*$arr1[7]/$spec_mass;
	}
	else
	{
		$Delta_mass_ppm='-';
	}
	$target_decoy=$arr1[15];
	if ($target_decoy =~ /target/i)
	{
		$target='TRUE';
	}
	$pro_AC=$arr1[12];
	print OUT "$serial\t$title\t$charge\t$sq\t$mod\t$score\t$spec_mass\t$Qvalue\t$Thero_mass\t$Delta_mass\t$Delta_mass_ppm\tS\t1\t$target\t$pro_AC\n";
	$serial++;
}
close IN;
close OUT;
