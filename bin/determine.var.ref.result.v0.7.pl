
#####################
#Project: 
#Program: 
#Version: 1.0
#Created on:
#Author: wusf
#####################


#这程序以前的名称叫get.var.ref.only.result.pl 
#v0.2：在以前的基础上加上判断是否变异属于等质量或和修饰等质量
#v0.3：2020-03-22 07:58:06 把line中的Ref和Var分开列，否则在mix部分没法区分，待完成。
#v0.4: 2020-06-27 16:41:34 增加一个输入参数，控制搜索引擎的种类，1为pfind，2为maxquant
#v0.5: 2020-12-17 兼容DIA  spectranaut的结果，搜索引擎参数值为3
#v0.7：2021-06-27 在line的信息中添加变异位点在肽段序列中的位置

#先获取
#需要输入的文件包括：
#1. SNP的注释，包含各变异的频率：前面用了id2site.line.v3.txt文件，后面用了*.exonic_variant_function.format.txt 文件，也是由l:\work\SAPproject\bin\get.id2site.pl 程序产生
$snp_anno_file=$ARGV[0];
#$snp_anno_file="id2site.line.v3.txt";
#2. 变异的注释文件
#$sap_anno_file="all.anno.1.v3.txt";
$sap_anno_file=$ARGV[1];
#3. pFind的搜库结果
#$dbSearch_result="hair-SV-v3-con.txt";
$search_engine=$ARGV[2];
$dbSearch_result=$ARGV[3];
$pro_seq_db=$ARGV[4];   #新增一个序列数据库的链接



sub judge_MS_replaceable_aa   #判断ref和alt的氨基酸是否是在质谱上可替代的氨基酸。可替代返回1，不可替代返回-1
{
	my ($refaa, $altaa)=@_;
	#要排除掉的情况：I/L，Q/K，N/D，Q/E，M/F
	my $return=-1;
	$return=1 if (($refaa=~/^I$/i) and ($altaa=~/^L$/i));
	$return=1 if (($refaa=~/^L$/i) and ($altaa=~/^I$/i));

	$return=1 if (($refaa=~/^Q$/i) and ($altaa=~/^K$/i));
	$return=1 if (($refaa=~/^K$/i) and ($altaa=~/^Q$/i));

	$return=1 if (($refaa=~/^N$/i) and ($altaa=~/^D$/i));
	$return=1 if (($refaa=~/^D$/i) and ($altaa=~/^N$/i));

	$return=1 if (($refaa=~/^Q$/i) and ($altaa=~/^E$/i));
	$return=1 if (($refaa=~/^E$/i) and ($altaa=~/^Q$/i));

	$return=1 if (($refaa=~/^M$/i) and ($altaa=~/^F$/i));
	$return=1 if (($refaa=~/^F$/i) and ($altaa=~/^M$/i));
	
	return $return;
}

sub data_precess #针对搜索引擎进行数据预处理，包括去除不需要的反库等，另外，按特定格式整理数据
{
	my ($search_engine, $dbSearch_result, $format_result)=@_;
	my $seq_data;
	print "result file: $dbSearch_result\n";
	open (IN_tmp, "$dbSearch_result")||die;
	open (OUT_tmp, ">$format_result")||die;
	print OUT_tmp "Sample\tPepseq\tProgroup\n";
	my $title_marker=1;
	if ($search_engine==1)  #pFind的搜库结果
	{
		while ($in=<IN_tmp>)
		{
			chomp $in;
			undef @arr1;
			
			@arr1=split("\t", $in);
			$sample=$arr1[1];
			$pepseq=$arr1[3];
			$progroup=$arr1[14];
			print OUT_tmp "$sample\t$pepseq\t$progroup\n";
		}
	}
	elsif ($search_engine==2)   #maxquant的搜库结果，针对的是evidence.txt文件
	{
		while ($in=<IN_tmp>)
		{
			chomp $in;
			undef @arr1;
			@arr1=split("\t", $in);
			if ($in=~/^Sequence\t/)
			{
				for ($i=0; $i<@arr1; $i++)
				{
					$colname=$arr1[$i];
					if ($colname eq 'Sequence')
					{
						$pepseq_col=$i;
					}
					elsif ($colname eq 'Modifications')
					{
						$modi_col=$i;
					}
					elsif ($colname eq 'Charge')
					{
						$charge_col=$i;
					}
					elsif ($colname eq 'Mass error [ppm]')
					{
						$error_ppm_col=$i;
					}
					elsif ($colname eq 'Mass error [Da]')
					{
						$error_da_col=$i;
					}
					elsif ($colname eq 'Raw file')
					{
						$raw_col=$i;
					}
					elsif ($colname eq 'Proteins')
					{
						$protein_col=$i;
					}
					elsif ($colname eq 'Reverse')
					{
						$rev_col=$i;
					}
					elsif ($colname eq 'MS/MS scan number')
					{
						$scan_col=$i;
					}
					elsif ($colname eq 'Potential contaminant')
					{
						$contaminant_col=$i;
					}

				}
			}
			else
			{
				$sample="$arr1[$raw_col]"."\.$arr1[$scan_col]";
				$pepseq=$arr1[$pepseq_col];
				$progroup=$arr1[$protein_col];
				$progroup=~s/\;/\//g;
				
				$rev_marker=$arr1[$rev_col];
				$contaminant_marker=$arr1[$contaminant_col];
			#	$contaminant_marker=$arr1[$protein_col];
			#	$contaminant_marker=$arr1[$protein_col];
			#	$contaminant_marker=$arr1[$protein_col];
			#	$contaminant_marker=$arr1[$protein_col];
			#	$contaminant_marker=$arr1[$protein_col];
			#	$contaminant_marker=$arr1[$protein_col];
			#	$contaminant_marker=$arr1[$protein_col];
				if (!($rev_marker))
				{
					$seq_data="$sample\t$pepseq\t$progroup";
					if (!($seq_result{$seq_data}))
					{
						print OUT_tmp "$sample\t$pepseq\t$progroup\n";
					}
					$seq_result{$seq_data}=1;
				}

			}

		}
		
	}
	elsif ($search_engine==3)   #DIA spectranut的结果
	{
		while ($in=<IN_tmp>)
		{
			chomp $in;
			undef @arr1;
			@arr1=split("\t", $in);
			if ($title_marker==1)
			{
				for ($i=0; $i<@arr1; $i++)
				{
					$str_tmp=$arr1[$i];
					if ($str_tmp eq 'R.FileName')
					{
						$sample_col=$i;
					}
					elsif ($str_tmp eq 'PG.ProteinGroups')
					{
						$progroup_col=$i;
					}
					elsif ($str_tmp eq 'PEP.StrippedSequence')
					{
						$pepseq_col=$i;
					}
					elsif ($str_tmp eq 'EG.IsDecoy')
					{
						$decoy_col=$i;
					}
					elsif ($str_tmp eq 'EG.Qvalue')
					{
						$Qvalue_col=$i;
					}

				}
				$title_marker=0;
			}
			else
			{
				for ($i=0; $i<@arr1; $i++)
				{
					$sample=$arr1[$sample_col];
					$pepseq=$arr1[$pepseq_col];
					$progroup=$arr1[$progroup_col];
					$progroup=~s/\;/\//g;
					$progroup.='/';
					$decoy=$arr1[$decoy_col];
					$Qvalue=$arr1[$Qvalue_col];
					if (($decoy =~ /FALSE/i) and ($Qvalue<0.01))
					{
						print OUT_tmp "$sample\t$pepseq\t$progroup\n";
					}
				}
			}
		}
	}
	elsif ($search_engine==4)   #Simple format 格式
	{
		while ($in=<IN_tmp>)
		{
			chomp $in;
			print OUT_tmp "$in\n";
		}
	}

	close IN_tmp;
	close OUT_tmp;
}

#open (IN_freq, "id2site.line.v3.txt")||die;  #从此文件获取各个变异的概率
open (IN_freq, "$snp_anno_file")||die;  #从此文件获取各个变异的概率
while ($in=<IN_freq>)
{
	chomp $in;
	undef @arr1;
	@arr1=split("\t", $in);
	$id=$arr1[0];
	$arr1[3]=~s/delins//;   #有些标签含有这种情况（如delinsGLGGA），但在蛋白中直接没有delins字符了
	$id_var_name="$id\t$arr1[2]$arr1[1]$arr1[3]";
	$id_var2freq{$id_var_name}=$arr1[4];
	$id_var2line{$id_var_name}=$arr1[5];
	$ref_aa=$arr1[2];
	$alt_aa=$arr1[3];
	#判断是否是质量相似，或和其修饰的结果相似
	$MS_replaceable{$id_var_name}=judge_MS_replaceable_aa($ref_aa, $alt_aa);
}
close IN_freq;

#open (IN_anno, "all.anno.1.v3.txt")||die;  #从此文件获取数据库中变异序列的ID和变异名称
open (IN_anno, "$sap_anno_file")||die;  #从此文件获取数据库中变异序列的ID和变异名称
while ($in=<IN_anno>)
{
	chomp $in;
	undef @arr1;
	@arr1=split("\t", $in);
	$idver=$arr1[0];
	$var=$arr1[1];
	$idver2var{$idver}=$var;
	#print "M-1: $idver\t$var\n" if ($idver =~ /ENST00000394001/);
	if ($idver=~/^(.+?)\.(.+)/)
	{
		$id=$1;
		$ver=$2;
		$id2ver{$id}.=",$ver";
	}
}
close IN_anno;

#debug--begin
#$id_debug="ENST00000391415";
#print "$id=$id2ver{$id_debug}=\n";
#debug--end

#插入个搜索结果的格式化过程：
$format_result=$dbSearch_result;
$format_result=~s/txt$//i;
$format_result.='format.txt';

print "result file1: $dbSearch_result\n";
data_precess($search_engine, $dbSearch_result, $format_result);   #调用子程序进行数据预处理

#open (IN_data, "hair-SV-v3-con.txt")||die;  #从此文件中获取鉴定结果
open (IN_data, "$format_result")||die;  #从此文件中获取鉴定结果
#输出结果：
$outfile=$dbSearch_result;
$outfile=~s/txt$//i;
$outfile_ref=$outfile;
$outfile_var=$outfile;
$outfile_mix=$outfile;
$outfile_ref2=$outfile;
$outfile_var2=$outfile;
$outfile_mix2=$outfile;
$outfile_var_ex=$outfile;  #因为质量数可替换，要排除掉的
$outfile_ref_ex=$outfile;  #因为质量数可替换，要排除掉的
$outfile_var_sigle_line=$outfile;
$outfile_ref_sigle_line=$outfile;
$outfile_var_multi_line=$outfile;
$outfile_ref_multi_line=$outfile;

$outfile.='result.txt';
$outfile_ref.="result.ref.txt";
$outfile_var.="result.var.txt";
$outfile_mix.="result.mix.txt";
$outfile_ref2.="result.ref.blank.txt";
$outfile_var2.="result.var.blank.txt";
$outfile_mix2.="result.mix.blank.txt";
$outfile_var_ex.="result.var_ex.txt";
$outfile_ref_ex.="result.ref_ex.txt";
$outfile_var_sigle_line.="result.var.single_line.txt";
$outfile_ref_sigle_line.="result.ref.single_line.txt";
$outfile_var_multi_line.="result.var.multi_line.txt";
$outfile_ref_multi_line.="result.ref.multi_line.txt";

#open (OUT, ">hair-SV-con.result.txt")||die;
#open (OUT_ref, ">hair-SV-con.result.ref.txt")||die;
#open (OUT_var, ">hair-SV-con.result.var.txt")||die;
#open (OUT_mix, ">hair-SV-con.result.mix.txt")||die;
#open (OUT_ref2, ">hair-SV-con.result.ref.2.txt")||die;
#open (OUT_var2, ">hair-SV-con.result.var.2.txt")||die;
#open (OUT_mix2, ">hair-SV-con.result.mix.2.txt")||die;

open (OUT, ">$outfile")||die;
open (OUT_ref, ">$outfile_ref")||die;
open (OUT_ref_ex, ">$outfile_ref_ex")||die;
open (OUT_var, ">$outfile_var")||die;
open (OUT_var_ex, ">$outfile_var_ex")||die;
open (OUT_mix, ">$outfile_mix")||die;
open (OUT_ref2, ">$outfile_ref2")||die;
open (OUT_var2, ">$outfile_var2")||die;
open (OUT_mix2, ">$outfile_mix2")||die;

open (OUT_ref_sl, ">$outfile_ref_sigle_line")||die;
open (OUT_ref_ml, ">$outfile_ref_multi_line")||die;
open (OUT_var_sl, ">$outfile_var_sigle_line")||die;
open (OUT_var_ml, ">$outfile_var_multi_line")||die;

while ($in=<IN_data>)
{
	chomp $in;
	undef @arr1;
	$MS_replaceable_marker=0;
	@arr1=split("\t", $in);
	#$sample=$arr1[1];
	#$pepseq=$arr1[3];
	#$progroup=$arr1[14];
	#$decoy=$arr1[13];
	$sample=$arr1[0];
	$pepseq=$arr1[1];
	$progroup=$arr1[2];
	#$decoy=$arr1[3];
	print OUT "$sample\t$pepseq\t$progroup\t";
	$sample_info="$sample\t$pepseq\t$progroup\t";
	#print OUT_ref "$sample\t$pepseq\t$progroup\t";
	#print OUT_var "$sample\t$pepseq\t$progroup\t";
	#print OUT_mix "$sample\t$pepseq\t$progroup\t";
	#next if ($decoy eq 'False');
	#用progroup判断是否是Ref或Var独有的。识别ENST的ID，对每个ID做个hash表，列出Ref和各个突变的版本，然后判断
	#1. 如果是ENST打头，且该ID没有.sv的后缀，则存在Ref－－调取该Ref对应的var的结果，判断是否缺失某个sv的结果
	#2. 如果是ENST打头，且有.sv后缀
	undef @arr2;
	@arr2=split("\/", $progroup);
	undef %pro_var_hash;
	$con_marker=0;
	$rev_marker=0;
	for ($i=0; $i<@arr2; $i++)
	{
		$protein=$arr2[$i];
		if ($protein=~/^ENST/)
		{
			if ($protein=~/^(.+)\.(.+)/)
			{
				$proname=$1;
				$varname=$2;
				$pro_var_hash{$proname}.=",$varname";
			}
			else
			{
				$pro_var_hash{$protein}.=",Ref";
			}
		}
		elsif ($protein=~/^CON_/)
		{
			$con_marker=1;
		}
		elsif ($protein=~/^REV_/)
		{
			$rev_marker=1;
		}
		else
		{
			print "Error1: $protein=\t$decoy=\n" if (($protein) and ($protein!~/^Protein AC/) and ($protein!~/^Progroup/));
		}
	}
	print OUT "$con_marker\t$rev_marker\t";
	$first_col="$con_marker\t$rev_marker\t";
	$pro_var="";
	$line_all="";
	$line_ref="";  #v0.3新增的
	$line_var="";  #v0.3新增的
	while (($key, $value)=each(%pro_var_hash))
	{
		$protein=$key;
		$var_all=$value;
		$var_all=~s/^,//;
		$var_list=$id2ver{$protein};
		undef @arr3;
		@arr3=split(",", $var_all);
		$Ref_marker=0;
		for ($i=0; $i<@arr3; $i++)
		{
			$var=$arr3[$i];
			if ($var eq 'Ref')
			{
				$Ref_marker=1;
			}
			else
			{
				$var_list=~s/$var//;  #将有的变异替换掉
			}
		}
		#print "$id=2=$var_list=\n";
		while ($var_list=~s/^,//) {}
		while ($var_list=~s/,$//) {}
		while ($var_list=~s/,,/,/) {}
		while ($var_all=~s/,,/,/) {}
		while ($var_all=~s/^,//) {}
		while ($var_all=~s/,$//) {}
		#print "$id=3=$var_list=\n";
		#print OUT "$protein\.$var_list\.$Ref_marker," if ($Ref_marker==1);
		if ($Ref_marker==1)
		{
			print OUT "$protein\.Ref\.$var_list,";
			undef @arr4;
			@arr4=split(",", $var_list);
			$line="";
			for ($i=0; $i<@arr4; $i++)
			{
				if ($arr4[$i])
				{
					$idver="$protein\.$arr4[$i]";
					$var_tmp=$idver2var{$idver};
					$id_var_name="$protein\t$var_tmp";
					if ($line)
					{
						$line.='&'.$id_var2line{$id_var_name};
					}
					else
					{
						$line=$id_var2line{$id_var_name};
					}
					$MS_replaceable_marker=1 if ($MS_replaceable{$id_var_name}==1);
				}
			}
			$line_all.="$line,";
			$line_ref.="$line,";
			$pro_var.="$protein\.Ref\.$var_list,";
		}
		#print OUT "$protein\.$var_all\.$Ref_marker," if ($Ref_marker==0);
		if ($Ref_marker==0)
		{
			print OUT "$protein\.Var\.$var_all,";
			undef @arr4;
			@arr4=split(",", $var_all);
			$line="";
			for ($i=0; $i<@arr4; $i++)
			{
				if ($arr4[$i])
				{
					$idver="$protein\.$arr4[$i]";
					$var_tmp=$idver2var{$idver};
					$id_var_name="$protein\t$var_tmp";
					#print "M0: $id_var_name\t=$var_tmp\t$idver\n" if ($pepseq eq 'CDLEWQNQEYQVLLDVR');
					if ($line)
					{
						$line.='&'.$id_var2line{$id_var_name};
					}
					else
					{
						$line=$id_var2line{$id_var_name};
					}
					$MS_replaceable_marker=1 if ($MS_replaceable{$id_var_name}==1);
				}
			}
			$line_all.="$line,";
			$line_var.="$line,";
			$pro_var.="$protein\.Var\.$var_all,";
		}
		#print "M0: $line_all\t$line_var\t$pro_var\t$var_all\n" if ($pepseq eq 'CDLEWQNQEYQVLLDVR');
		if ($var_list=~/sv/i)
		{
		#	print "$var_list\t$Ref_marker\n";
		}
	}
	#对$line_all的变量去冗余
	while ($line_all=~s/,,/,/) {}
	while ($line_all=~s/,$//) {}
	undef @arr4;
	@arr4=split(",", $line_all);
	$line_all2="";
	undef %line_tmp_hash;
	for ($i=0; $i<@arr4; $i++)
	{
		$line=$arr4[$i];
		if (!($line_tmp_hash{$line}))
		{
			$line_all2.="$line,";
			$line_tmp_hash{$line}=1;
		}
	}
	$line_all2=~s/,$//;

	#对$line_ref的变量去冗余
	while ($line_ref=~s/,,/,/) {}
	while ($line_ref=~s/,$//) {}
	undef @arr4;
	@arr4=split(",", $line_ref);
	$line_ref2="";
	undef %line_tmp_hash;
	for ($i=0; $i<@arr4; $i++)
	{
		$line=$arr4[$i];
		if (!($line_tmp_hash{$line}))
		{
			$line_ref2.="$line,";
			$line_tmp_hash{$line}=1;
		}
	}
	$line_ref2=~s/,$//;

	#对$line_var的变量去冗余
	while ($line_var=~s/,,/,/) {}
	while ($line_var=~s/,$//) {}
	undef @arr4;
	@arr4=split(",", $line_var);
	$line_var2="";
	undef %line_tmp_hash;
	for ($i=0; $i<@arr4; $i++)
	{
		$line=$arr4[$i];
		if (!($line_tmp_hash{$line}))
		{
			$line_var2.="$line,";
			$line_tmp_hash{$line}=1;
		}
	}
	$line_var2=~s/,$//;

	if ($pro_var!~/Ref\.,/)  #没有匹配上‘Ref.,’，说明没有完全匹配到没有变异位置的参考型中。即所匹配的肽段包含了有变异的位点
	{
		$line_num=0;
		#print "M1: $line_num\t$pro_var\t$line_all2\n" if ($pepseq eq 'CDLEWQNQEYQVLLDVR');
		while ($line_all2=~s/\&\&/\&/g){}
		$line_all2=~s/\&$//;
		$line_all2=~s/^\&//;
		while ($line_all2=~/line/gi)
		{
			$line_num++;
		}
		while ($line_all2=~/\&/g)  #这种字符是表示不同的line之间是共同存在同一个肽段中，不算是多匹配，一般在Ref中
		{
			$line_num--;
		}
		#print "M2: $line_num\t$pro_var\n" if ($pepseq eq 'CDLEWQNQEYQVLLDVR');
		if (($pro_var=~/Ref\.sv/) and ($pro_var!~/Var\.sv/))
		{
			#print OUT_ref "$sample_info$first_col$pro_var\t$line_all2\n";
			if ($MS_replaceable_marker==1)
			{
				print OUT_ref_ex "$sample_info$first_col$pro_var\t$line_all2\t$line_ref2\t$line_var2\n";
			}
			else
			{
				print OUT_ref "$sample_info$first_col$pro_var\t$line_all2\t$line_ref2\t$line_var2\n";
				if ($line_num>0)
				{
					if ($line_num==1)
					{
						print OUT_ref_sl "$sample_info$first_col$pro_var\t$line_all2\t$line_ref2\t$line_var2\n";
					}
					else
					{
						print OUT_ref_ml "$sample_info$first_col$pro_var\t$line_all2\t$line_ref2\t$line_var2\n";
					}
				}
				else
				{
					print "Error2: line num <=0. $line_all2\t$line_ref2\t$line_var2\t==$pro_var==\n";
				}
			}
		}
		if (($pro_var=~/Var\.sv/) and ($pro_var!~/Ref\.sv/))
		{
			if ($MS_replaceable_marker==1)
			{
				print OUT_var_ex "$sample_info$first_col$pro_var\t$line_all2\t$line_ref2\t$line_var2\n";
			}
			else
			{
				print OUT_var "$sample_info$first_col$pro_var\t$line_all2\t$line_ref2\t$line_var2\n";
				print "$sample_info$first_col$pro_var\t$line_all2\t$line_ref2\t$line_var2\n" if ($pepseq eq 'CDLEWQNQEYQVLLDVR');
				if ($line_num>0)
				{
					if ($line_num==1)
					{
						print OUT_var_sl "$sample_info$first_col$pro_var\t$line_all2\t$line_ref2\t$line_var2\n";
					}
					else
					{
						print OUT_var_ml "$sample_info$first_col$pro_var\t$line_all2\t$line_ref2\t$line_var2\n";
					}
				}
				else
				{
					print "Error3: line num <=0. =$line_num=\t$line_all2\t$line_ref2\t$line_var2\t==$pro_var==\n$in\n";
				}
			}
		}
		if (($pro_var=~/Var\.sv/) and ($pro_var=~/Ref\.sv/))
		{
			print OUT_mix "$sample_info$first_col$pro_var\t$line_all2\t$line_ref2\t$line_var2\n";
		}
	}
	else   #同时匹配上了没有变异位点的区域，这部分的结果应该是不要的吧
	{
		if (($pro_var=~/Ref\.sv/) and ($pro_var!~/Var\.sv/))
		{
			print OUT_ref2 "$sample_info$first_col$pro_var\t$line_all2\t$line_ref2\t$line_var2\n";
		}
		if (($pro_var=~/Var\.sv/) and ($pro_var!~/Ref\.sv/))
		{
			print OUT_var2 "$sample_info$first_col$pro_var\t$line_all2\t$line_ref2\t$line_var2\n";
			#if ($MS_replaceable_marker==1)
			#{
			#	print OUT_var2_ex "$sample_info$first_col$pro_var\t$line_all2\n";
			#}
			#else
			#{
			#	print OUT_var2 "$sample_info$first_col$pro_var\t$line_all2\n";
			#}
		}
		if (($pro_var=~/Var\.sv/) and ($pro_var=~/Ref\.sv/))
		{
			print OUT_mix2 "$sample_info$first_col$pro_var\t$line_all2\t$line_ref2\t$line_var2\n";
		}
	}


	print OUT "\n";
}

close IN_data;
close OUT;
close OUT_ref;
close OUT_var;
close OUT_var_ex;
close OUT_mix;
close OUT_ref2;
close OUT_var2;
close OUT_var2_ex;
close OUT_mix2;

close OUT_ref_sl;
close OUT_ref_ml;
close OUT_var_sl;
close OUT_var_ml;
