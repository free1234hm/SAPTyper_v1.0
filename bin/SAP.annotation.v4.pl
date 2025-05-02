
#####################
#Project: 
#Program: 
#Version: 1.0
#Created on:
#Author: wusf
#####################

#v2: 增加判断是否存在全外的参考型解析数据，如果有，另外再生成一份结果
#v3: 把输出的注释结果改为每个sample和group都可以同时和各个wes的结果匹配
#v4: 把和全外匹配的部分删掉，加入判断基因型的步骤

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
			$raw_i++;
			$raw2group{$raw_file}=$group;
			$raw2sample{$raw_file}=$sample;
			$sample2group{$sample}=$group;
			$sample2raw{$sample}.=",$raw_file";  #由于Raw文件是绝对路径，不能那么写了
			$group2raw{$group}.=",$raw_file";  #由于Raw文件是绝对路径，不能那么写了
			#$raw_file2=$raw_file;
			#$raw_file2=~s/^(.+\\)/$1$sample2raw{$sample}/;
			
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
			#$QC_file=~s/^.+\\//;
			#$QC_file.='.depth.QC.txt' if ($QC_file!~/\.depth\.QC\.txt$/);
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
	$sample_dir{$sample_exp}=1;
	$group_exp="group.$group";
	$sample_dir{$group_exp}=1;
}

print "Loading the parameter file complete...\n";
############

print "$snp_anno\n";
open (IN_line2pro, "$snp_anno")||die;  #M0的<out file>中含有line到蛋白的对应关系
while ($in=<IN_line2pro>)
{
	chomp $in;
	undef @arr1;
	@arr1=split("\t", $in);
	$pro=$arr1[0];
	$line=$arr1[5];
	$snp=$arr1[6];
	$freq=$arr1[4];
	$sap="$arr1[0]\/$arr1[2]$arr1[1]$arr1[3]";
	$line2pro{$line}=$pro;
	$line2snp{$line}=$snp;
	$line2snp_anno{$line}="$snp\t$sap\t$freq";
}
close IN_line2pro;

print "$gene_anno\n";
open (IN_gene_anno, "$gene_anno")||die; #protein annotation
while ($in=<IN_gene_anno>)
{
	chomp $in;
	undef @arr1;
	@arr1=split("\t", $in);
	$pro=$arr1[0];
	#$line=$arr1[5];
	$pro2anno{$pro}=$in;
}
close IN_gene_anno;

open (IN_rsID, "$SNP2rsID_file")||die; #rs ID
while ($in=<IN_rsID>)
{
	chomp $in;
	undef @arr1;
	@arr1=split("\t", $in);
	$snp=$arr1[0];
	$rsid=$arr1[1];
	#$line=$arr1[5];
	$snp2rsid{$snp}=$rsid;
}
close IN_rsID;

#新增的检测基因型的判断，同时输出结果
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

	$stat_file="$sample_group_exp".'.result.detail.stat.xls';
	print "Stat file: $stat_file\n";
	open (IN_stat, "$stat_file")||die;
	$stat_anno_file="$sample_group_exp".'.result.detail.stat.annotation.xls';

	print "Stat anno file: $stat_anno_file\n";
	open (OUT_stat_anno, ">$stat_anno_file")||die;

	$title=1;
	while ($in=<IN_stat>)
	{
		chomp $in;
		if ($title)
		{
			print OUT_stat_anno "$in\tSNP	SAP	Freq	Protein	GeneID	SwissProtID	Gene	Des	rsID	GT_Pro\n";
			$title=0;
			next;
		}
		undef @arr1;
		@arr1=split("\t", $in);
		$line=$arr1[0];
		$pro=$line2pro{$line};
		$snp=$line2snp{$line};
		$snp_slim=$snp;
		$snp_slim=~s/>.+/>/;
		$snp_wes=$snp_slim2snp{$snp_slim};
		$snp_anno2=$line2snp_anno{$line};
		$snp_anno2="\t\t" if (!($snp_anno2));
		$pro_anno=$pro2anno{$pro};
		$pro_anno="				" if (!($pro_anno));
		$rsid=$snp2rsid{$snp};
		#$rsid="" if (!($rsid));
		$ref_uni=$arr1[1];
		$var_uni=$arr1[2];
		$ref_mul=$arr1[3];
		$var_mul=$arr1[4];
		$M_ref_mul=$arr1[5];
		$M_var_mul=$arr1[6];

		$GT_exam="";
		if ($ref_uni)
		{
			if ($var_uni)
			{
				$GT_exam='0/1';
			}
			elsif (($var_mul) or ($M_var_mul))
			{
				$GT_exam='0/[1]';
			}
			else
			{
				$GT_exam='0/0';
			}
		}
		elsif (($ref_mul) or ($M_ref_mul))
		{
			if ($var_uni)
			{
				$GT_exam='[0]/1';
			}
			elsif (($var_mul) or ($M_var_mul))
			{
				$GT_exam='[0]/[1]';
			}
			else
			{
				$GT_exam='[0]/[0]';
			}
		}
		else
		{
			if ($var_uni)
			{
				$GT_exam='1/1';
			}
			elsif (($var_mul) or ($M_var_mul))
			{
				$GT_exam='[1]/[1]';
			}
		}

		$refwes_QCinfo=$group_snp2refwes{$wes_group}{$snp_slim};
		$refwes_QCinfo="\t\t\t" if (!($refwes_QCinfo));
		print OUT_stat_anno "$in\t$snp_anno2\t$pro_anno\t$rsid\t $GT_exam\n";

	}
	close IN_stat;
	close OUT_stat_anno;
}



#以下是删除的用对于样本的全外做注释的结果
=pod
#$QCfile2group{$QC_file}=$group;
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
	for ($i=0; $i<@wes_group_arr; $i++)
	{
		$wes_group=$wes_group_arr[$i];

		$stat_file="$sample_group_exp".'.result.detail.stat.xls';
		print "Stat file: $stat_file\n";
		open (IN_stat, "$stat_file")||die;
		$stat_anno_file="$sample_group_exp.$wes_group".'.result.detail.stat.annotation.xls';

		print "Stat anno file: $stat_anno_file\n";
		open (OUT_stat_anno, ">$stat_anno_file")||die;

		$ref_wes_anno=$group2wesref{$wes_group};
		if ($ref_wes_anno)
		{
			$refwes_anno_file="$sample_group_exp.$wes_group".'.result.detail.stat.wesref.annotation.xls';
			open (OUT_refwes_anno, ">$refwes_anno_file")||die;
		}

		$title=1;
		while ($in=<IN_stat>)
		{
			chomp $in;
			if ($title)
			{
				print OUT_stat_anno "$in\tSNP	SAP	Freq	Protein	GeneID	SwissProtID	Gene	Des	rsID	GT	AD	DP	GQ\n";
				print OUT_refwes_anno "$in\tSNP	SAP	Freq	Protein	GeneID	SwissProtID	Gene	Des	rsID	GT	AD	DP	GQ\tSNP_wes\n" if (ref_wes_anno);
				$title=0;
				next;
			}
			undef @arr1;
			@arr1=split("\t", $in);
			$line=$arr1[0];
			$pro=$line2pro{$line};
			$snp=$line2snp{$line};
			$snp_slim=$snp;
			$snp_slim=~s/>.+/>/;
			$snp_wes=$snp_slim2snp{$snp_slim};
			$snp_anno2=$line2snp_anno{$line};
			$snp_anno2="\t\t" if (!($snp_anno2));
			$pro_anno=$pro2anno{$pro};
			$pro_anno="				" if (!($pro_anno));
			$rsid=$snp2rsid{$snp};
			#$rsid="" if (!($rsid));
			$QCinfo=$group_snp2QCinfo{$wes_group}{$snp};
			$refwes_QCinfo=$group_snp2refwes{$wes_group}{$snp_slim};
			$refwes_QCinfo="\t\t\t" if (!($refwes_QCinfo));
			print OUT_stat_anno "$in\t$snp_anno2\t$pro_anno\t$rsid\t$QCinfo\n";
			print OUT_refwes_anno "$in\t$snp_anno2\t$pro_anno\t$rsid\t$refwes_QCinfo\t$snp_wes\n" if (ref_wes_anno);

		}
		close IN_stat;
		close OUT_stat_anno;
		close OUT_refwes_anno  if (ref_wes_anno);
	}
}
=cut
