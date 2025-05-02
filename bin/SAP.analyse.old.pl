
#####################
#Project: 
#Program: 
#Version: 1.0
#Created on:
#Author: wusf
#####################

#总调度程序，进行如下步骤处理
#1. 读取VB配置得到的参数文件，根据模块选择，以及各个参数进行后续的分析
#2. 按模块选择的情况，分不同的部分去操作
#这三个模块分步骤运行，只有在操作参数时需要判断前面的步骤是否有


use File::Basename;
$program_dir=dirname(__FILE__);

$par_file=$ARGV[0];
#$log_file=$ARGV[1];

=pod
#pre-check  预检查程序---不在这，应该在VB中直接执行，执行完后判断是否预检查完成
$cmd="$program_dir".'\pre.check.pl '."$par_file";
system ("$cmd");
=cut

$anno_files="$program_dir".'\anno_data_file.txt';
open (IN, "$anno_files")||die;
while ($in=<IN>)
{
	chomp $in;
	print "==$in\n";
	undef @arr1;
	@arr1=split("=", $in);
	$anno_file_hash{$arr1[0]}=$arr1[1];
}
close IN;


#open (IN_project, "project.par")||die;
open (IN_project, "$par_file")||die;
while ($in=<IN_project>)
{
	chomp $in;
	undef @arr1;
	@arr1=split("=", $in);
	$par_term2value_hash{$arr1[0]}=$arr1[1];
}
close IN_project;

$current_path=$par_term2value_hash{'result_fold'};
chdir ($current_path);

if ($par_term2value_hash{'DS_search_engine1'} eq 'True')
{
	$DS_search_engine='pFind';
	$search_engine_id=1;
}
elsif ($par_term2value_hash{'DS_search_engine2'} eq 'True')
{
	$DS_search_engine='MaxQuant';
	$search_engine_id=2;
}
elsif ($par_term2value_hash{'DS_search_engine3'} eq 'True')
{
	$DS_search_engine='Spectronaut';
	$search_engine_id=3;
}

if ($par_term2value_hash{'DA_search_engine1'} eq 'True')
{
	$DA_search_engine='pFind';
	$search_engine_id=1;
}
elsif ($par_term2value_hash{'DA_search_engine2'} eq 'True')
{
	$DA_search_engine='MaxQuant';
	$search_engine_id=2;
}
elsif ($par_term2value_hash{'DA_search_engine3'} eq 'True')
{
	$DA_search_engine='Spectronaut';
	$search_engine_id=3;
}
elsif ($par_term2value_hash{'DA_search_engine4'} eq 'True')
{
	$DA_search_engine='SimpleFormat';
	$search_engine_id=4;
}

#以下开始逐个模块运行
#system ("type nul > $log_file");  #对log文件清零

if ($par_term2value_hash{'db_creation'}==1)  #第一个模块
{
	#产生vcf2snp程序需要的参数文件(直接采用‘wes_file’的文件)
	$vcf_par=$par_term2value_hash{'wes_file'};
	$cmd='perl '.$program_dir.'\vcf2snp.pl '.$vcf_par;
	system ("$cmd");#debug

	#执行第二个命令
	#$S1_outfile=$par_term2value_hash{'outfile_marker'};  #改成固定名称，不再需要用户自己定义
	$S1_outfile='var.result.txt';
	$add_pubilc=$par_term2value_hash{'add_var'};
	if ($add_pubilc==1)
	{
		$public_var_file=$par_term2value_hash{'public_var'};
	}
	else
	{
		$public_var_file='';
	}
	$add_ref=$par_term2value_hash{'add_ref'};
	$add_ref=0 if ($add_ref!=1);
	$add_wes=$par_term2value_hash{'add_wes'};
	$add_wes=0 if ($add_wes!=1);
	$cmd='perl '.$program_dir.'\SAP.S1.wes2prodb.v2.pl all.vcf.list.txt '."$S1_outfile 1 $add_ref $public_var_file $add_wes";  #all.vcf.list.txt是vcf2snp的程序的输出结果
	print "$cmd\n";
	system ("$cmd");  #debug
}
if ($par_term2value_hash{'db_searching'}==1)  #第二个模块
{
	if ($par_term2value_hash{'db_creation'}==1)  #有第一个模块
	{
		#搜索数据库从第一个模块中获得
		$db_file="$current_path".'\result\merge.refdb.customic.fa';
		$db_file=~s/\\\\/\\/;
	}
	else
	{
		#搜索数据库从参数文件中获得
		$db_file=$par_term2value_hash{'seq_db'};
	}

	#搜库模块，等待添加
	print "search engine: $DS_search_engine\n";
	if ($DS_search_engine eq 'pFind')
	{
		$ref_par=$par_term2value_hash{'ref_parameter'};
		$pparse_par=$par_term2value_hash{'pparse_parameter'};
		$pro_sample_file=$par_term2value_hash{'pro_file'};
		$add_con=$par_term2value_hash{'add_con'};
		$search_result_dir="$par_term2value_hash{'result_fold'}\\dbsearch";
		$cmd="md $search_result_dir";
		print "$cmd\n";
		system ("$cmd");  #debug
		$cmd="perl $program_dir\\pFind.auto.pl $ref_par $pparse_par $db_file $pro_sample_file $add_con $search_result_dir";
		print "$cmd\n";
		system ("$cmd");  #debug
		#解析pFind的搜索结果，整合成export的result的格式--放在pfind.auto.pl程序里面了
		print "pFind search finished....\n";
		$search_result_tmp=$search_result_dir.'\result\result.txt';
		$search_result_tmp=~s/\\\\/\\/g;
	}
	elsif ($DS_search_engine eq 'MaxQuant')
	{
		$ref_par=$par_term2value_hash{'ref_parameter'};
		$pro_sample_file=$par_term2value_hash{'pro_file'};
		$add_con=$par_term2value_hash{'add_con'};
		$search_result_dir="$par_term2value_hash{'result_fold'}\\dbsearch";
		$cmd="md $search_result_dir";
		print "$cmd\n";
		system ("$cmd");  #debug
		$cmd="perl $program_dir\\maxquant.auto.pl $ref_par $db_file $pro_sample_file $add_con $search_result_dir";
		print "$cmd\n";
		system ("$cmd");  #debug
		#解析pFind的搜索结果，整合成export的result的格式--放在pfind.auto.pl程序里面了
		print "Maxquant search finished....\n";
		$search_result_tmp=$search_result_dir.'\combined\txt\evidence.txt';
		$search_result_tmp=~s/\\\\/\\/g;
	}
}
if ($par_term2value_hash{'data_ana'}==1)  #第三个模块
{
	if ($par_term2value_hash{'db_creation'}==1)  #有第一个模块
	{
		#SNP和SAP注释从第一个模块的结果中获得
		$SNP_anno="$current_path".'\result\\'.$S1_outfile;
		if ($par_term2value_hash{var_type1} eq 'True')
		{
			$SAP_anno="$current_path".'\result\all.anno.singleSAP.txt';
		}
		else
		{
			$SAP_anno="$current_path".'\result\all.anno.all.txt';  #这部分后面的程序不兼容，删除
		}
	}
	else
	{
		$SNP_anno=$par_term2value_hash{'SNP_Anno'};
		$SAP_anno=$par_term2value_hash{'SAP_Anno'};
		
		#产生vcf2snp程序需要的参数文件(直接采用‘wes_file’的文件)   #如果没有第一步，这不会产生，会导致后续结果出不来
		$vcf_par=$par_term2value_hash{'wes_file'};
		$cmd='perl '.$program_dir.'\vcf2snp.pl '.$vcf_par;
		system ("$cmd");#debug
	}
	if ($par_term2value_hash{'db_searching'}==1)  #有第二个模块
	{
		#搜索引擎和搜索结果从第二个模块的结果中获得
		$search_engine=$DS_search_engine;
		$search_result=$search_result_tmp;  
	}
	else
	{
		$search_engine=$DA_search_engine;
		$search_result=$par_term2value_hash{'pro_result'};
	}

	#需要创建SAP.S3.data_ana.par文件，放在project目录下
	open (OUT_par, ">SAP.S3.data_ana.par")||die;

	print OUT_par '#parameters for M3 module';
	print OUT_par "\nSNP_annotation	$SNP_anno";
	print OUT_par '	#End, 这文件可以是从M0_2的输出的<out file>即可';

	print OUT_par "\n\nSAP_annotation	$SAP_anno";
	print OUT_par '	#End，这文件是M0_2输出的all.anno.singleSAP.txt 文件';

	print OUT_par "\n\nSearch_engine	$search_engine_id	";
	print OUT_par '#1. pFind; 2. Maxquant

Wes results (Individual, Wes QC file):
';

	$vcf_par=$par_term2value_hash{'wes_file'};
	open (IN_vcf, "$vcf_par")||die;
	while ($in_vcf=<IN_vcf>)
	{
		chomp $in_vcf;
		next if ($in_vcf=~/^Marker	File/);
		print OUT_par "\t$in_vcf\n";
	}
	close IN_vcf;
	print OUT_par '
Experimental design and raw files(Individual,Sample,Fraction,Raw file):
';

	$proteome_par=$par_term2value_hash{'pro_file'};
	open (IN_pro, "$proteome_par")||die;
	while ($in_pro=<IN_pro>)
	{
		chomp $in_pro;
		next if ($in_pro=~/^Person	Sample	/);
		print OUT_par "\t$in_pro\n";
	}
	close IN_pro;
	
	print OUT_par '
Proteome results:	#这些raw文件是合并搜索的，最终给出一个搜索结果
';
	print OUT_par "\t$search_result\n";
	print OUT_par "
Protein_annotation	$anno_file_hash{'pro_anno'}	\#End

SNP2rsID	$anno_file_hash{'snp2rs'} 	\#END

Gene annotation	$anno_file_hash{'gene_anno'}	\#END
";
	close OUT_par;

	#执行程序SAP.S3.data_ana.pl
	$cmd="perl $program_dir".'\SAP.S3.data_ana.pl SAP.S3.data_ana.par';
	print "CMD=$cmd\n";
	system ("$cmd");
}

