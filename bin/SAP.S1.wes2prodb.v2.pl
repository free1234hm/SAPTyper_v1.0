
#把以前写的代码搬过来，把相应的需要的数据库等也搬过来到该流程的子目录下
#原代码路径：r:\work\SAPproject\software\annovar\create.ref.var.database.pl 

#########################
#Project: SAP Project
#这程序是个总的调度程序，包括以下功能：
#1. 对全外测序结果的SNP转化为SAP
#2. 将新生成的SAP和ExAC获得的SAP合并，并做标签
#3. 根据SAP列表，生成带有SAP的序列（可以有多种方案可选）
#4. 整合成Ref+Var的序列数据库
#Version: 1.0
#Created on:
#Author: wusf
#########################

use File::Basename;
$program_dir=dirname(__FILE__);
#print "$program_dir\n";
$pipeline_path='r:\work\SAPproject\pipeline';
$pipeline_path=$program_dir;

use Cwd;
$current_path=getcwd();
$current_path=~s/\//\\/g;
print "current path: $current_path\n";
#$project_path='l:\work\SAPproject\old_data_from_FL\\';  #这路径用于存放结果，建议每个项目设置一个文件夹，相关的数据都会拷贝到这文件夹中
#$project_path='l:\work\SAPproject\Data\EAS_0.001\\';  #这路径用于存放结果，建议每个项目设置一个文件夹，相关的数据都会拷贝到这文件夹中
system ("md result");
$project_path="$current_path".'\result';  #在用当前路径下建立一个result的目录作为结果路径
$project_path=~s/\\\\/\\/;
#$annovar_path='l:\work\SAPproject\software\annovar\\';
$annovar_path="$pipeline_path".'\bin\\';  #在pipeline的path中拷贝了一份，测试

$pro_db="$pipeline_path".'\data\coding.protein.clean.fa';
#$public_SNP='l:\work\SAPproject\software\annovar\ExAC_nonTCGA.r1.sites.vep.vcf.extract.EAS001.v2.txt';
#$public_SAP='l:\work\SAPproject\software\annovar\ExAC.EAS001.v2.exonic_variant_function';
#$public_SAP='l:\work\SAPproject\software\annovar\ExAC.EAS001.v2.exonic_variant_function.format.txt';
#$public_SAP='l:\work\SAPproject\software\annovar\blank.file.txt';
$public_SAP="$pipeline_path".'\data\ExAC.EAS001.new.exonic_variant_function.format.txt';  #新的从屈武斌那获取的ExAC数据，用相对路径
#$customic_SNP='l:\work\SAPproject\old_data_from_FL\rs.var.format.txt';
#$customic_SNP='l:\work\SAPproject\Data\EAS_0.001\hg19_exac03_0.001.txt';  #SNP数据
$customic_SNP="$current_path"."\\$ARGV[0]";  #输入的SNP文件是第一个命令行参数
#$merge_SNP='l:\work\SAPproject\software\annovar\ExAC_nonTCGA.r1.sites.vep.vcf.extract.EAS001.v2.FL.txt';  #这是合并后的文件名
#$merge_SAP='l:\work\SAPproject\software\annovar\ExAC.EAS001.v2.FL.exonic_variant_function.format.txt';  #这是合并后的文件名
#$merge_SAP='l:\work\SAPproject\Data\EAS_0.001\ExAC.EAS001.new.exonic_variant_function.format.txt';
$merge_SAP=$ARGV[1];  #输出文件是第二个命令行参数
$customic_SNP_format=$ARGV[2]; #
$customic_SNP_format=1 if (!($customic_SNP_format));  #1－测序结果的vcf格式；2－
$add_ref=$ARGV[3];
$public_var_file=$ARGV[4];
$public_SAP=$public_var_file;
$add_wes=$ARGV[5];
#=pod
#1. 对全外测序结果的SNP转化为SAP
#修改下，转换个思路。公共SNP都是要转化的，完全可以不用重复转化，因此直接导入转化完的公共SAP，只是把新生成的转化了即可
#open (IN_pub_SNP, "$public_SNP")||die;
#open (IN_pub_SAP, "$public_SAP")||die;
#open (IN_new_SNP, "$new_SNP")||die;
#open (OUT_merge_SNP, ">$merge_SNP")||die;
open (OUT_merge_SAP, ">$merge_SAP")||die;

chdir $annovar_path;

$out_title='customic.SNP';
$cmd="$pipeline_path".'\bin\annotate_variation.pl -out '."$out_title".' -build hg19 -hgvs '."$customic_SNP".'  humandb/ -dbtype ensGene';
print "CMD: $cmd\n";
system ("$cmd");  #运行annovar的命令，获取个性化SNP转化的SAP列表。这结果应该会得到"customic.SNP.exonic_variant_function"文件

#2. 将新生成的SAP和ExAC获得的SAP合并，并做标签
#进行SAP输出结果的整理
$process_annovar_out_cmd="$pipeline_path".'\bin\get.id2site.pl';
$cmd="$process_annovar_out_cmd $out_title".'.exonic_variant_function';
print "CMD: $cmd\n";

$customic_SAP_table="$out_title".".exonic_variant_function.format.txt";
print "SAP_table: $customic_SAP_table\n";

system ("$cmd"); #运行get.id2site.pl命令，将个性化SAP的结果整理出来

#把运行结果拷贝到结果目录下
$cmd="copy $pipeline_path\\bin\\customic.SNP.* $current_path";
print "CMD: $cmd\n";

chdir $current_path;
system ("$cmd");

if ($add_wes)
{
	open (IN_customic_SAP, "$customic_SAP_table")||die;
	$i=0;
	while ($in=<IN_customic_SAP>)
	{
		chomp $in;
		undef @arr1;
		@arr1=split("\t", $in);
		next if ($arr1[3]=~/X/i);  #如果碰到突变为终止密码子的，跳过
		$customic_pro_var[$i]="$arr1[0]\t$arr1[1]\t$arr1[2]\t$arr1[3]";
		$customic_source{$customic_pro_var[$i]}="$arr1[4]";
		$customic_lines{$customic_pro_var[$i]}="$arr1[5]";
		$customic_loc{$customic_pro_var[$i]}="$arr1[6]";
		$customic_var_hash{$customic_pro_var[$i]}=1;
		$i++;
	}
	close IN_customic_SAP;
}
if ($public_SAP)
{
	open (IN_pub_SAP, "$public_SAP")||die;
	while ($in=<IN_pub_SAP>)
	{
		chomp $in;
		undef @arr1;
		@arr1=split("\t", $in);
		$pub_pro_var="$arr1[0]\t$arr1[1]\t$arr1[2]\t$arr1[3]";
		$pub_source="$arr1[4]";
		$pub_lines="$arr1[5]";
		$pub_loc="$arr1[6]";
		if ($customic_var_hash{$pub_pro_var})
		{
			#print OUT_merge_SAP "$pub_pro_var\t$pub_source\t$pub_lines\t$pub_loc\t$customic_source{$pub_pro_var}\t$customic_loc{$pub_pro_var}\n";
			print OUT_merge_SAP "$pub_pro_var\t$pub_source\t$pub_lines\/$customic_source{$pub_pro_var}\t$pub_loc\n";  #把SNP的来源marker和line的行放一块
			$customic_var_hash{$pub_pro_var}=2;
		}
		else
		{
			print OUT_merge_SAP "$pub_pro_var\t$pub_source\t$pub_lines\t$pub_loc\n";
		}
	}
	close IN_pub_SAP;
}

for ($i=0; $i<@customic_pro_var; $i++)
{
	$customic_pro_var2=$customic_pro_var[$i];
	if ($customic_var_hash{$customic_pro_var2}==1)
	{
		print OUT_merge_SAP "$customic_pro_var2\t$customic_source{$customic_pro_var2}\t$customic_lines{$customic_pro_var2}_$customic_source{$customic_pro_var2}\t$customic_loc{$customic_pro_var2}\n";
	}
}

close OUT_merge_SAP;
#close IN_pub_SAP;
#=cut

chdir $current_path;
$cmd="copy $merge_SAP $project_path";
print "CMD: $cmd\n";
system ("$cmd");

#3. 根据SAP列表，生成带有SAP的序列（可以有多种方案可选）
chdir $project_path;
$create_var_seq_cmd="$pipeline_path".'\bin\create.var_pro.alternative.pl';
$cmd="$create_var_seq_cmd $merge_SAP $pro_db";
print "CMD: $cmd\n";
system ("$cmd");

#4. 整合成Ref+Var的序列数据库
if ($add_ref)
{
	$cmd="type $pro_db $project_path".'\all.seq.singleSAP.fa'." >$project_path".'\merge.refdb.customic.fa';
}
else
{
	$cmd="type $project_path".'\all.seq.singleSAP.fa'." >$project_path".'\merge.refdb.customic.fa';
}
print "CMD=$cmd\n";
#die;
system ("$cmd");
