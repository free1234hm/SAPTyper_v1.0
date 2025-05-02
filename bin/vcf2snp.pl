
#####################
#Project: 
#Program: 
#Version: 1.0
#Created on:
#Author: wusf
#####################

#Get SNP information from vcf files. The vcf files are from the WES (whole exon sequencing) results.
#The program will produce two kind of files:
#1. *.vcf.list.txt: merge all the vcf results into one file, which will be used to create the customic protein sequence database
#2. *.vcf.depth.QC.txt: detail information of vcf result, which will used for the further analysis.


$GQ_cutoff=99;

$par_file=$ARGV[0];

$outfile_SNP='all.vcf.list.txt';
open (OUT_SNP, ">$outfile_SNP")||die;
print OUT_SNP "\#Chr	Start	End	Ref	Alt\tMarker\n";

open (IN_par, "$par_file")||die;
while ($in_par=<IN_par>)
{
	chomp $in_par;
	undef @arr1;
	@arr1=split("\t", $in_par);
	$vcf_file=$arr1[1];  #为了和后面统一，把这格式改了，marker放前面去了
	$vcf_marker=$arr1[0];
	next if ($in_par=~/^Marker	File/);



	#$infile='FL.raw.vcf';
	#$infile=$ARGV[0];
	$infile=$vcf_file;
	#$sample_marker=$infile;
	#$sample_marker=~s/\.raw\.vcf$//;
	$outfile_QC="$vcf_marker".'.vcf.depth.QC.txt';



	open (IN, "$infile")||die;
	open (OUT_QC, ">$outfile_QC")||die;
	open (OUT_unmatch, ">$vcf_marker.multi_alt.txt")||die;
	#open (OUT2, ">$outfile.unmatch")||die;

	print OUT_QC "SNP_info\tGT\tAD\tDP\tGQ\n";
	while ($in=<IN>)
	{
		chomp $in;
		next if ($in=~/^\#/);
		undef @arr1;
		@arr1=split("\t", $in);
		$chr=@arr1[0];
		$pos=@arr1[1];
		$rsid=@arr1[2];
		$ref=@arr1[3];
		$alts=@arr1[4];
		$ref_len=length($ref);
		$alt_len=length($alts);
		undef @arr1_1;
		@arr1_1=split(",", $alts);
		for ($k=0; $k<@arr1_1; $k++)
		{
			$alt=$arr1_1[$k];
			$alt_len=length($alt);

			if (($ref_len==1) and ($alt_len==1))
			{
				;
			}
			elsif (($ref_len>1) and ($alt_len==1))
			{
			#	if ($ref=~s/^$alt//)
			#	{
			#		$alt='-';
			#	}
			}
			elsif (($ref_len==1) and ($alt_len>1))
			{
			#	if ($alt=~s/^$ref//)
			#	{
			#		$ref='-';
			#	}
			}
			else
			{
				print OUT_unmatch "$in\n";
			}


			$format=@arr1[8];
			$data=@arr1[9];
			undef %format2value;
			undef @arr2;
			@arr2=split(':', $format);
			undef @arr3;
			@arr3=split(':', $data);
			for ($i=0; $i<@arr2; $i++)
			{
				$format2value{$arr2[$i]}=$arr3[$i];
			}

			$SNP_info="$chr".'_'."$pos$ref>$alt";

			if ($format2value{'GQ'}>=$GQ_cutoff)
			{
				print OUT_QC "$SNP_info\t$format2value{'GT'}\t$format2value{'AD'}\t$format2value{'DP'}\t$format2value{'GQ'}\n";
				$chr2=$chr;
				$chr2=~s/chr//;
				$chr_loc2marker{$SNP_info}.="\/$vcf_marker";
				$chr_loc2info{$SNP_info}="$chr2\t$pos\t$pos\t$ref\t$alt";
	#			print OUT2 "$chr2\t$pos\t$pos\t$ref\t$alt\n";
			}
		}
		

	}

	close IN;
	close OUT_QC;
	close OUT_unmatch;

}

while (($key, $value)=each(%chr_loc2marker))
{
	$markers=$value;
	$SNP_info=$key;
	$markers=~s/^\///;
	print OUT_SNP "$chr_loc2info{$SNP_info}\t$markers\n";
}
close OUT_SNP;
