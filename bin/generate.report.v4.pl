
#####################
#Project: 
#Program: 
#Version: 1.0
#Created on:
#Author: wusf
#####################

#v2: 20210208根据WuJL要求修订report格式
#v3: 增加判断是否存在全外的参考型解析数据，如果有，另外再生成一份结果
#v4: 把输出的结果改为每个sample和group都可以同时和各个wes的结果匹配

$par_file=$ARGV[0];

open (IN_par, "$par_file")||die;

#$pfind_i=0;
#$in_par=<IN_par>;
$raw_i=0;
while ($in_par=<IN_par>)
{
	chomp $in_par;
	#print "$in_par==\n";
	if ($in_par=~/^\SNP_annotation\t(.+?)\t/)
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
			$pro_raw_arr[$raw_i]="$raw_file";
		#	print "$raw_file\n";
			$raw_i++;
			$raw2group{$raw_file}=$group;
			$raw2sample{$raw_file}=$sample;
			$sample2group{$sample}=$group;
			$sample2raw{$sample}.=",$raw_file";
			$group2raw{$group}.=",$raw_file";
		}
		next;
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
			#$QC_file=$arr1[2];
			$group=$arr1[1];
			$QC_file="$group.vcf.depth.QC.txt";
			$ref_wes_anno=$arr1[3];  #新增的全外参考型
			if ($QC_file)
			{
				$QCfile2group{$QC_file}=$group;
				$group2QCfile{$group}=$QC_file;
				$group2wesref{$group}=$ref_wes_anno;  #新增的全外参考型
			}
			else
			{
				last;
			}
		#	$QCfile2group{$QC_file}=$group;
		#	$group2QCfile{$group}=$QC_file;
		}
		next;
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


for ($i=0; $i<@pro_raw_arr; $i++)
{
	$raw_file=$pro_raw_arr[$i];
	$group=$raw2group{$raw_file};
	$sample=$raw2sample{$raw_file};
	#add the following code for the sample and group marker
	$sample_exp="sample.$sample";
	$sample_dir{$sample_exp}=1;   #do not need the sample information, using group only.
	$group_exp="group.$group";
	$sample_dir{$group_exp}=1;
	#print "$group_exp\n";
}

undef @wes_group_arr;
$wes_i=0;
while (($key, $value)=each(%QCfile2group))
{
	$QC_file=$key;
	$group=$value;
	$wes_group_arr[$wes_i]=$group;
	$wes_i++;
	print "QC group: =$group=\n";
	print "QC file: $QC_file\n";
	open (IN_QCinfo, "$QC_file")||die; #*.vcf.depth.QC.txt 文件的信息,需要遍历各个group
	while ($in=<IN_QCinfo>)
	{
		chomp $in;
		undef @arr1;
		@arr1=split("\t", $in);
		$snp=$arr1[0];
		$QCinfo="$arr1[1]\t$arr1[2]\t$arr1[3]\t$arr1[4]";
		$group_snp2QCinfo{$group}{$snp}=$QCinfo;
	}
	close IN_QCinfo;

	$ref_wes_anno=$group2wesref{$group};
	if ($ref_wes_anno)
	{
		open (IN_refwes, "$ref_wes_anno")||die;
		while ($in=<IN_refwes>)
		{
			chomp $in;
			undef @arr1;
			@arr1=split("\t", $in);
			$chr=$arr1[1];
			$pos=$arr1[2];
			$ref=$arr1[3];
			$alt=$arr1[4];
			$filter=$arr1[11];
			if (($filter =~ /PASS/i) or ($filter =~ /LowQuality/i))
			{
				next;
			}
			$genotype=$arr1[5];
			if (($genotype ne 'AA') and ($genotype ne 'GG') and ($genotype ne 'TT') and ($genotype ne 'CC'))
			{
				next;
			}
			$detail_info=$arr1[10];
			$snp="$chr".'_'."$pos$ref".'>'."$alt";
			$snp_slim="$chr".'_'."$pos$ref".'>';
			undef @arr2;
			@arr2=split(":", $detail_info);
			$GT=$arr2[0];
			$AD=$arr2[1];
			$DP=$arr2[2];
			$GQ=$arr2[3];
			$refwes_QCinfo="$GT\t$AD\t$DP\t$GQ";
			$group_snp2refwes{$group}{$snp}=$refwes_QCinfo;
			$group_snp2refwes{$group}{$snp_slim}=$refwes_QCinfo;
			$snp_slim2snp{$snp_slim}=$snp;
		}
		close IN_refwes;
	}
}

print "Loading the parameter file complete...\n";
############

while (($key, $value)=each(%sample_dir))
{

	$sample_group_exp=$key;
	$sample=$sample_group_exp;
	$group="";
	if ($sample=~/^sample/)
	{
		$sample=~s/^sample\.//;
		$group=$sample2group{$sample};
	}
	elsif ($sample=~/^group/)
	{
		$sample=~s/^group\.//;
		$group=$sample;
	}

	#新增一个wes group的循环，每个样本都算一遍
	for ($ii=0; $ii<@wes_group_arr; $ii++)
	{
		$wes_group=$wes_group_arr[$ii];

		print "==$sample_group_exp\t$wes_group\n";

		#$group_exp=$key;
		#$group=$group_exp;
		#$group=~s/group\.//;
		$report_file="$sample_group_exp.$wes_group.report.xls";
		$report_file2="$sample_group_exp.$wes_group.wesref.report.xls";
		$anno_file="$sample_group_exp.result.detail.stat.annotation.xls";
		#$anno_file2="$sample_group_exp.$wes_group.result.detail.stat.wesref.annotation.xls";
		#print "$anno_file\n";
		undef @anno_file_arr;
		undef @report_file_arr;
		$anno_file_arr[0]=$anno_file;
		$report_file_arr[0]=$report_file;
		$kk=1;
		if ($group2wesref{$wes_group})
		{
			$anno_file_arr[1]=$anno_file;
			$report_file_arr[1]=$report_file2;
			$kk++;
		}
		for ($k=0; $k<$kk; $k++)
		{

			$anno_file_tmp=$anno_file_arr[$k];
			$report_file_tmp=$report_file_arr[$k];
			print "Anno_file: $anno_file_tmp\n";
			print "Report_file: $report_file_tmp\n";
			open (IN, "$anno_file_tmp")||die;
			open (OUT, ">$report_file_tmp")||die;
			print OUT "SNP	GT_pro	GT\tFreq\tGene\trsID	AD\tDP	SAP\tline\tref_uni	var_uni	ref_mul	var_mul	M_ref_mul	M_var_mul\tSNP_wes\n";
			$title=1;
			while ($in=<IN>)
			{
				chomp $in;
				undef @arr1;
				@arr1=split("\t", $in);
				if ($title)
				{
					$title=0;
					$id_col=0;
					for ($i=0; $i<@arr1; $i++)
					{
						$col_name=$arr1[$i];
						if ($col_name eq 'GT')
						{
							$GT_col=$i;
						}
						elsif ($col_name eq 'ref_uni')
						{
							$ref_col=$i;
						}
						elsif ($col_name eq 'var_uni')
						{
							$var_col=$i;
						}
						elsif ($col_name eq 'ref_mul')
						{
							$ref_col2=$i;
						}
						elsif ($col_name eq 'var_mul')
						{
							$var_col2=$i;
						}
						elsif ($col_name eq 'M_ref_mul')
						{
							$ref_col3=$i;
						}
						elsif ($col_name eq 'M_var_mul')
						{
							$var_col3=$i;
						}
						elsif ($col_name eq 'SNP')
						{
							$snp_col=$i;
						}
						elsif ($col_name eq 'SAP')
						{
							$sap_col=$i;
						}
		##
						elsif ($col_name eq 'Freq')
						{
							$freq_col=$i;
						}
						elsif ($col_name eq 'Gene')
						{
							$GeneName_col=$i;
						}
						elsif ($col_name eq 'rsID')
						{
							$rsID_col=$i;
						}
		##
						elsif ($col_name eq 'AD')
						{
							$AD_col=$i;
						}
		##
						elsif ($col_name eq 'DP')
						{
							$DP_col=$i;
						}
						elsif ($col_name eq 'SNP_wes')
						{
							$SNP_wes_col=$i;
						}
					}
					next;
				}
				else
				{
					#$GT=$arr1[$GT_col];
					$snp=$arr1[$snp_col];
					$snp_slim=$snp;
					$snp_slim=~s/>.+/>/;
					$QCinfo="";
					if ($k==0)
					{
						$QCinfo=$group_snp2QCinfo{$wes_group}{$snp};
					}
					else
					{
						$QCinfo=$group_snp2refwes{$wes_group}{$snp_slim};
					}
					undef @arr3;
					@arr3=split("\t", $QCinfo);
					$GT=$arr3[0];
					$AD=$arr3[1];
					$DP=$arr3[2];
					$GQ=$arr3[3];
					if ($GT)
					{
						$id=$arr1[$id_col];
						$snp=$arr1[$snp_col];
						$sap=$arr1[$sap_col];
						$ref=$arr1[$ref_col];
						$var=$arr1[$var_col];
						$ref2=$arr1[$ref_col2];
						$var2=$arr1[$var_col2];
						$ref3=$arr1[$ref_col3];
						$var3=$arr1[$var_col3];
						$GT2='-';
						if ($ref)
						{
							if ($var)
							{
								$GT2='0/1';
							}
							elsif (($var2) or ($var3))
							{
								$GT2='0/[1]';
							}
							else
							{
								$GT2='0/0';
							}
						}
						elsif (($ref2) or ($ref3))
						{
							if ($var)
							{
								$GT2='[0]/1';
							}
							elsif (($var2) or ($var3))
							{
								$GT2='[0]/[1]';
							}
							else
							{
								$GT2='[0]/[0]';
							}
						}
						else
						{
							if ($var)
							{
								$GT2='1/1';
							}
							elsif (($var2) or ($var3))
							{
								$GT2='[1]/[1]';
							}
							else
							{
								$GT2='-';
							}
						}
						#print OUT "$id\t$snp\t$sap\t$GT2\t$GT\n";
						#print OUT "	SNP	SAP	GT_exam	GT\n";
						print OUT "$snp\t $GT2\t $GT\t$arr1[$freq_col]\t$arr1[$GeneName_col]\t$arr1[$rsID_col]\t$AD\t$DP\t$sap\t$id\t$arr1[$ref_col]\t$arr1[$var_col]\t$arr1[$ref_col2]\t$arr1[$var_col2]\t$arr1[$ref_col3]\t$arr1[$var_col3]\t$arr1[$SNP_wes_col]\n";
					}
				}
			}
			close IN;
			close OUT;
		}
	}
}

