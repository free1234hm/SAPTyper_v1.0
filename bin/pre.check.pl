
#####################
#Project: 
#Program: 
#Version: 1.0
#Created on:
#Author: wusf
#####################
use File::Basename;
my $program_dir=dirname(__FILE__);
use Cwd;
my $current_dir = getcwd;
$pipeline_path=$program_dir;

#此程序是SAPTyper的预检查程序，检查需要的输入文件是否存在/或命名是否正确，数据格式是否符合要求，数据文件间的相容性等，将结果输出到一个文件，供VB程序读取
$par_file=$ARGV[0];
#$par_file='r:\work\SAPproject\pipeline\interface\debug\all_test_20210417\project.2.par';

#读项目参数文件，确定有哪些参数和文件需要判断
open (IN_par, "$par_file")||die;
while ($in=<IN_par>)
{
	chomp $in;
	undef @arr1;
	@arr1=split("=", $in);
	$term2value{$arr1[0]}=$arr1[1];
}
close IN_par;

chdir ($term2value{'result_fold'});

$dep_program_dir="$program_dir".'\dep_program.dir';
open (IN, "$dep_program_dir")||die;
while ($in=<IN>)
{
	chomp $in;
	#print "==$in\n";
	undef @arr1;
	@arr1=split("=", $in);
	$dep_program_hash{$arr1[0]}=$arr1[1];
}
close IN;

$anno_files="$program_dir".'\anno_data_file.txt';
open (IN, "$anno_files")||die;
while ($in=<IN>)
{
	chomp $in;
	#print "==$in\n";
	undef @arr1;
	@arr1=split("=", $in);
	$anno_file_hash{$arr1[0]}=$arr1[1];
}
close IN;

#判断各个文件是否存在，以及文件格式是否合规的子程序

sub check_wes_file
{
	$file=$term2value{'wes_file'};
		$wes_file_exist=-e $file;
		if (!($wes_file_exist))
		{
			return 0;
		}

	open (IN, "$file")||die;
	$line=0;
	$marker=0;
	$file_marker=1;
	while ($in=<IN>)
	{
		chomp $in;
		while ($in=~s/[\t ]$//g)
		{}
		next if (!$in);

		if ($line==0)
		{
			if ($in =~/^Marker	File/)
			{
				$title_marker=1;
			}
			else
			{
				$title_marker=0;
				last;
			}
		}
		else
		{
			undef @arr1;
			@arr1=split("\t", $in);
			if ($arr1[0])
			{
				$vcf_exist=-e $arr1[1];
				if (($vcf_exist) and ($arr1[1]=~/\.vcf$/))
				{
					$vcf_marker=1;
				}
				else
				{
					$vcf_marker=0;
					last;
				}
			}
			else
			{
				$file_marker=0;
				last;
			}
		}
		$line++;
	}
	if (($title_marker) and ($vcf_marker) and ($file_marker))
	{
		return 1;
	}
	else
	{
		return 0;
	}
}

sub check_proteome_file  #用check_wes_file改出来的，变量没变，内容变了
{
	$file=$term2value{'pro_file'};
		$proteome_file_exist=-e $file;
		if (!($proteome_file_exist))
		{
			return 0;
		}
	open (IN, "$file")||die;
	$line=0;
	$marker=0;
	$file_marker=1;
	while ($in=<IN>)
	{
		chomp $in;
		while ($in=~s/[\t ]$//g)
		{}
		next if (!$in);

		if ($line==0)
		{
			if ($in =~/^Person	Sample	Fraction	File/)
			{
				$title_marker=1;
			}
			else
			{
				$title_marker=0;
				last;
			}
		}
		else
		{
			undef @arr1;
			@arr1=split("\t", $in);
			if (($arr1[0]) and ($arr1[1]) and ($arr1[2]))
			{
				$vcf_exist=-e $arr1[3];
				if (($vcf_exist) and ($arr1[3]=~/\.raw[ ]*$/))
				{
					$vcf_marker=1;
				}
				else
				{
					$vcf_marker=0;
					last;
				}
			}
			else
			{
				$file_marker=0;
				last;
			}
		}
		$line++;
	}
	#print "$title_marker\t=$vcf_marker\t==$file_marker\n";
	if (($title_marker) and ($vcf_marker) and ($file_marker))
	{
		return 1;
	}
	else
	{
		return 0;
	}
}

sub check_public_var
{
	$file=$term2value{'public_var'};
	$file_exist=-e $file;
	if ($file_exist)
	{
		return 1;
	}
	else
	{
		return 0;
	}
}

sub check_ref_seq  #需要检查默认的ref seq的地址
{
	$pro_db="$pipeline_path".'\data\coding.protein.clean.fa';
	$file=$pro_db;
	$file_exist=-e $file;
	if ($file_exist)
	{
		return 1;
	}
	else
	{
		return 0;
	}
}

sub check_pfind_par
{
	$file=$term2value{'ref_parameter'};
	$file_exist=-e $file;
	if (($file_exist) and ($file=~/\.cfg$/))
	{
		return 1;
	}
	else
	{
		return 0;
	}
}

sub check_pparse_par
{
	$file=$term2value{'pparse_parameter'};
	$file_exist=-e $file;
	if (($file_exist) and ($file=~/\.cfg$/))
	{
		return 1;
	}
	else
	{
		return 0;
	}
}

sub check_maxquant_par
{
	$file=$term2value{'ref_parameter'};
	$file_exist=-e $file;
	if (($file_exist) and ($file=~/\.xml$/))
	{
		return 1;
	}
	else
	{
		return 0;
	}
}

sub check_add_con_par  #需要检查默认路径下的污染库数据
{
	$con_db="$program_dir".'\data\contaminate.sequence.txt';
	#$pro_db="$pipeline_path".'\data\coding.protein.clean.fa';
	$file=$con_db;
	$file_exist=-e $file;
	if ($file_exist)
	{
		return 1;
	}
	else
	{
		return 0;
	}
}

sub check_seq_db
{
	$file=$term2value{'seq_db'};
	$file_exist=-e $file;
	if ($file_exist)
	{
		return 1;
	}
	else
	{
		return 0;
	}
}

sub check_SNP_anno
{
	$file=$term2value{'SNP_Anno'};
	$file_exist=-e $file;
	if ($file_exist)
	{
		return 1;
	}
	else
	{
		return 0;
	}
}

sub check_SAP_anno
{
	$file=$term2value{'SAP_Anno'};
	$file_exist=-e $file;
	if ($file_exist)
	{
		return 1;
	}
	else
	{
		return 0;
	}
}

sub check_proteome_result
{
	$file=$term2value{'pro_result'};
	$file_exist=-e $file;
	if ($file_exist)
	{
		return 1;
	}
	else
	{
		return 0;
	}
}




#判断
$check_result="";
if ($term2value{'db_creation'}==1)
{
	$wes_file_state=check_wes_file();
	$public_var_state=check_public_var();
	$ref_seq_state=check_ref_seq();
	if (!$wes_file_state)
	{
		$check_result.="Error in wes file\n";
	}
	if (!$public_var_state)
	{
		$check_result.="Error in public var file\n";
	}
	if (!$ref_seq_state)
	{
		$check_result.="Error in reference sequences file\n";
	}
	if ($term2value{'db_searching'}==1)
	{
		$proteome_file_state=check_proteome_file();
		if (!$proteome_file_state)
		{
			$check_result.="Error in proteome parameter file\n";
		}
		if ($term2value{DS_search_engine1} eq 'True') #pfind
		{
			$pfind_par_state=check_pfind_par();
			$pparse_par_state=check_pparse_par();
			if (!$pfind_par_state)
			{
				$check_result.="Error in pfind parameter file\n";
			}
			if (!$pparse_par_state)
			{
				$check_result.="Error in pparse parameter file\n";
			}
		}
		elsif ($term2value{DS_search_engine2} eq 'True') #maxquant
		{
			$maxquant_par_state=check_maxquant_par();
			if (!$maxquant_par_state)
			{
				$check_result.="Error in maxquant parameter file\n";
			}
		}
		if (($add_con{DS_search_engine2}==1) and ($term2value{DS_search_engine1} eq 'True'))
		{
			$add_con_state=check_add_con_par();
			if (!$add_con_state)
			{
				$check_result.="Error in contaminate file\n";
			}
		}

		if ($term2value{'data_ana'}==1)
		{
			#没有额外需要判断的
		}
		else
		{
			#没有额外需要判断的
		}
	}
	else
	{
		if ($term2value{'data_ana'}==1)
		{
			#没有这种情况
		}
		else
		{
			#没有额外需要判断的
		}
	}
}
else
{
	if ($term2value{'db_searching'}==1)
	{
		$wes_file_state=check_wes_file();
		$proteome_file_state=check_proteome_file();
		if (!$wes_file_state)
		{
			$check_result.="Error in wes file\n";
		}
		if (!$proteome_file_state)
		{
			$check_result.="Error in proteome parameter file\n";
		}

		if ($term2value{DS_search_engine1} eq 'True') #pfind
		{
			$pfind_par_state=check_pfind_par();
			$pparse_par_state=check_pparse_par();
			if (!$pfind_par_state)
			{
				$check_result.="Error in pfind parameter file\n";
			}
			if (!$pparse_par_state)
			{
				$check_result.="Error in pparse parameter file\n";
			}
		}
		elsif ($term2value{DS_search_engine2} eq 'True') #maxquant
		{
			$maxquant_par_state=check_maxquant_par();
			if (!$maxquant_par_state)
			{
				$check_result.="Error in maxquant parameter file\n";
			}
		}
		if ($add_con{DS_search_engine2}==1)
		{
			$add_con_state=check_add_con_par();
			if (!$add_con_state)
			{
				$check_result.="Error in contaminate file\n";
			}
		}
		$seq_db_check=check_seq_db();

		if ($term2value{'data_ana'}==1)
		{
			$SNP_anno_state=check_SNP_anno();
			$SAP_anno_state=check_SAP_anno();
			if (!$SNP_anno_state)
			{
				$check_result.="Error in SNP annotation file\n";
			}
			if (!$SAP_anno_state)
			{
				$check_result.="Error in SAP annotation file\n";
			}
		}
		else
		{
			#没有额外需要判断的
		}
	}
	else
	{
		if ($term2value{'data_ana'}==1)
		{
			$SNP_anno_state=check_SNP_anno();
			$SAP_anno_state=check_SAP_anno();
			if (!$SNP_anno_state)
			{
				$check_result.="Error in SNP annotation file\n";
			}
			if (!$SAP_anno_state)
			{
				$check_result.="Error in SAP annotation file\n";
			}
			$proteome_result_state=check_proteome_result();
			if (!$proteome_result_state)
			{
				$check_result.="Error in proteome result file\n";
			}
		}
		else
		{
			#没有额外需要判断的
		}
	}
}

if ($term2value{'db_searching'}==1)  #增加判断预设的
{
	if ($term2value{DS_search_engine1} eq 'True') #pfind
	{
		$pfind_dir=$dep_program_hash{'pFind'};
		$pfind_pparse="$pfind_dir\\"."Searcher.exe";
		$pfind_search="$pfind_dir\\"."pParse.exe";
		if (!((-e $pfind_pparse) and (-e $pfind_search)))
		{
			$check_result.="Error in pFind program\n";
		}
	}
	elsif ($term2value{DS_search_engine2} eq 'True') #maxquant
	{
		$maxquant_dir=$dep_program_hash{'MaxQuant'};
		$maxquant_search="$maxquant_dir\\"."MaxQuantCmd.exe";
		if (!(-e $maxquant_search))
		{
			$check_result.="Error in MaxQuant program\n";
		}
	}
}

if ($term2value{'data_ana'}==1)  #增加判断预设的
{
	$pro_anno_file=$anno_file_hash{'pro_anno'};
	$snp2rs_file=$anno_file_hash{'snp2rs'};
	$gene_anno_file=$anno_file_hash{'gene_anno'};
	if (!(-e $pro_anno_file))
	{
		$check_result.="Error in protein annotation file\n";
	}
	if (!(-e $snp2rs_file))
	{
		$check_result.="Error in SNP2rs annotation file\n";
	}
	if (!(-e $gene_anno_file))
	{
		$check_result.="Error in gene annotation file\n";
	}
}

print "$check_result";

open (OUT, ">precheck.result.txt")||die;
print OUT "$check_result";
close OUT;