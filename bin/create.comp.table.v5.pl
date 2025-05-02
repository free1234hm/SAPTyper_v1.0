
#v3: 加入peptide的信息
#v4: 不同sample间的数据有串行，可能某些地方没清零，改bug
#v5： 改为流程中可以嵌入的程序，读入参数文件，输出每个样本，每个个体的结果

#命令目录
use File::Basename;
$program_dir=dirname(__FILE__);
#print "$program_dir\n";
#$pipeline_path='r:\work\SAPproject\pipeline';
$cmd_path="$program_dir\\";

#$cmd_path='r:\work\SAPproject\pipeline\\';
$par_file='SAP.S3.data_ana.par';

#读入参数文件：
open (IN_par, "$par_file")||die;

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
			$raw_file=~s/^.+\\//;
			$group=$arr1[1];
			$sample=$arr1[2];
			$pro_raw_arr[$raw_i]="$raw_file";
		#	print "Raw: $raw_i\t$pro_raw_arr[$raw_i]\n";
			$raw_i++;
			$raw2group{$raw_file}=$group;
			$raw2sample{$raw_file}=$sample;
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
			$QC_file=$arr1[2];
			$group=$arr1[1];
			print "$group==\t$QC_file==\n";
			if ($QC_file)
			{
				$QCfile2group{$QC_file}=$group;
				$group2QCfile{$group}=$QC_file;
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

$files='.sep.result.mix.blank.txt,.sep.result.mix.txt,.sep.result.ref.blank.txt,.sep.result.ref.multi_line.txt,.sep.result.ref.single_line.txt,.sep.result.ref.txt,.sep.result.ref_ex.txt,.sep.result.txt,.sep.result.var.blank.txt,.sep.result.var.multi_line.txt,.sep.result.var.single_line.txt,.sep.result.var.txt,.sep.result.var_ex.txt';
@file_arr= split (',', $files);
#print "@file_arr\n";
#for ($i=0; $i<@file_arr; $i++)
#{
#	print "$file_arr[$i]\n";
#}

#die;
for ($i=0; $i<@pro_raw_arr; $i++)
{
	$raw_file=$pro_raw_arr[$i];
	$group=$raw2group{$raw_file};
	$sample=$raw2sample{$raw_file};
	$raw_file=~s/\.raw$//;
	#print "Raw2: $raw_file\n";
	for ($j=0; $j<@file_arr; $j++)
	{
		$pro_seq_result_file='proteome.'."$raw_file"."$file_arr[$j]";  #因为用了绝对路径，不能那么写---可以不管，前面把绝对路径去掉了
		#$pro_seq_result_file="$raw_file"."$file_arr[$j]";
		$pro_seq_result_file=~s/^(.+\\)/$1proteome\./;
		$sample_exp="sample.$sample";
		$sample_dir{$sample_exp}.="$pro_seq_result_file\n";
		$group_exp="group.$group";
		$sample_dir{$group_exp}.="$pro_seq_result_file\n";
		#print "$sample_exp\t$group_exp\n";
	}
	#$pro_seq_result_file='proteome.'."$raw_file".'.sep.txt';
}

####################
=pod
system ("dir hair-*.result.*.txt /b > result.all.dir");

open (IN_dir, "result.all.dir")||die;
while ($in_dir=<IN_dir>)
{
	chomp $in_dir;
	if ($in_dir=~/^hair-(.+?)\./)
	{
		$sample=$1;
		$sample_exp="$1".'_'."$2";
		$sample_dir{All}.="$in_dir\n";
		$sample_dir{$sample}.="$in_dir\n";
		$sample_dir{$sample_exp}.="$in_dir\n";
	}
}
=cut

while (($key1, $value1)=each(%sample_dir))
{
	$sample=$key1;
	$dir_list=$value1;
	data_process($sample, $dir_list);
}
close IN_dir;

#####

sub data_process
{
	#system ("dir hair-F202-0724-*.result.*.txt /b > result.all.dir");
	
	#open (IN_dir, "result.all.dir")||die;
	#while ($in_dir=<IN_dir>)
	my ($sample, $dir_list)=@_;
	
	open (OUT_dir, ">$sample.result.dir")||die;
	print OUT_dir "$dir_list";
	close OUT_dir;
	
	my %loc_hash;
	my %datatype_loc;
	my %datatype_loc2pepseq;
	my $in_dir;
	my $data_type;
	my $data_type2;
	my $lines;
	my @datatype_arr;
	my $pepseq_all;
	my $uni_pepseq;
	my $pepseq;
	my %pepseq2specnum;
	my %spec_hash;
	
	undef @arr0;
	@arr0=split("\n", $dir_list);
	for ($ii=0; $ii<@arr0; $ii++)
	{
		$in_dir=$arr0[$ii];
		chomp $in_dir;
		next if (!($in_dir));
		print "$in_dir\n";
		open (IN, "$in_dir")||die;
		$data_type=$in_dir;
		print "data type: $data_type\n";
		$data_type=~s/^.+result\.(.+)\.txt/$1/;
		next if ($data_type eq 'var');
		next if ($data_type eq 'ref');
		#print "data type: $data_type\n";
		if ($data_type eq 'var.single_line')
		{
			#$data_type_hash{var_uni}=1;
			$data_type2='var_uni';
		}
		elsif ($data_type eq 'ref.single_line')
		{
			#$data_type_hash{ref_uni}=1;
			$data_type2='ref_uni';
		}
		elsif ($data_type=~/ref/)
		{
			#$data_type_hash{ref_mul}=1;
			$data_type2='ref_mul';
		}
		elsif ($data_type=~/var/)
		{
			#$data_type_hash{var_mul}=1;
			$data_type2='var_mul';
		}
		elsif ($data_type=~/mix/)
		{
			#$data_type_hash{ref_mul}=1;
			#$data_type_hash{var_mul}=1;
			$data_type2='mix';
		}
		else
		{
			print "Skip: $data_type     source: $in_dir\n";
		}
		while ($in=<IN>)
		{
			chomp $in;
			undef @arr1;
			@arr1=split("\t", $in);
			$lines=$arr1[6];
			$pepseq=$arr1[1];   #v3新增
			$spec=$arr1[0];
		#	if ($pepseq eq 'PQCCQSVCCQPTCCR')
		#	{
		#		print "==$sample\t$spec\n";
		#	}
			if ($spec_hash{$pepseq}!~/\n$spec\n/)
			{
				$spec_hash{$pepseq}.="\n$spec\n";
				$pepseq2specnum{$pepseq}++;
			}
			$pro_var_id==$arr1[5];   #v3新增
			undef @arr2;
			$lines2=$lines;
			$lines2=~s/\;/\t/g;
			$lines2=~s/\,/\t/g;
			$lines2=~s/\&/\t/g;
			@arr2=split("\t", $lines2);
			if ($data_type2 eq 'mix')  #对于mix的，要额外处理去判断哪些位点是ref，哪些位点是var
			{
				#ref
				$data_type_tmp='M_ref_mul';
				$lines=$arr1[7];
				undef @arr2;
				$lines2=$lines;
				$lines2=~s/\;/\t/g;
				$lines2=~s/\,/\t/g;
				$lines2=~s/\&/\t/g;
				@arr2=split("\t", $lines2);
				for ($i=0; $i<@arr2; $i++)
				{
					$loc=$arr2[$i];
					$loc_hash{$loc}=1;
					$datatype_loc{$data_type_tmp}{$loc}=1;
					$datatype_loc2pepseq{$data_type_tmp}{$loc}.=",$pepseq";
					#$datatype_loc{$data_type2}{$loc}=1;
				}
	
				#var
				$data_type_tmp='M_var_mul';
				$lines=$arr1[8];
				undef @arr2;
				$lines2=$lines;
				$lines2=~s/\;/\t/g;
				$lines2=~s/\,/\t/g;
				$lines2=~s/\&/\t/g;
				@arr2=split("\t", $lines2);
				for ($i=0; $i<@arr2; $i++)
				{
					$loc=$arr2[$i];
					$loc_hash{$loc}=1;
					$datatype_loc{$data_type_tmp}{$loc}=1;
					$datatype_loc2pepseq{$data_type_tmp}{$loc}.=",$pepseq";
				}
				
	
=pod
				undef @arr3;
				$pro_var=$arr1[5];
				$pro_var=~s/,sv/\.sv/g;
				$pro_var=~s/\;/\t/g;
				$pro_var=~s/\,/\t/g;
				$pro_var=~s/\&/\t/g;
				@arr3=split("\t", $pro_var);
				for ($j=0; $j<@arr3; $j++)
				{
					$pro=$arr3[$j];
					$loc=$arr2[$j];
					$loc_hash{$loc}=1;
					if ($pro=~/\.Ref\./)
					{
						$data_type_tmp='M_ref_mul';
						$datatype_loc{$data_type_tmp}{$loc}=1;
					}
					elsif ($pro=~/\.Var\./)
					{
						$data_type_tmp='M_var_mul';
						$datatype_loc{$data_type_tmp}{$loc}=1;
					}
					else
					{
						print "Error for the mix data: =$pro= =$pro_var= =$j= =$arr1[0]\n";
					}
				}
=cut			
			}
			else
			{
				for ($i=0; $i<@arr2; $i++)
				{
					$loc=$arr2[$i];
					$loc_hash{$loc}=1;
					$datatype_loc{$data_type2}{$loc}=1;
					$datatype_loc2pepseq{$data_type2}{$loc}.=",$pepseq";
				}
			}
		}
	}
	
	open (OUT, ">$sample.result.stat.xls")||die;
	open (OUT2, ">$sample.result.detail.stat.xls")||die;
	$datatype_arr[0]='ref_uni';
	$datatype_arr[1]='var_uni';
	$datatype_arr[2]='ref_mul';
	$datatype_arr[3]='var_mul';
	$datatype_arr[4]='M_ref_mul';
	$datatype_arr[5]='M_var_mul';
	for ($i=0; $i<@datatype_arr; $i++)
	{
		print OUT "\t$datatype_arr[$i]";
		print OUT2 "\t$datatype_arr[$i]";
	}
	print OUT "\n";
	print OUT2 "\n";
	while (($key, $value)=each(%loc_hash))
	{
		$loc=$key;
		next if (!($loc));
		print OUT "$loc";
		print OUT2 "$loc";
		for ($i=0; $i<@datatype_arr; $i++)
		{
			$data_type=$datatype_arr[$i];
			$result=0;
			$result2=$datatype_loc{$data_type}{$loc};
			$result=$result2 if ($result2);
			print OUT "\t$result";
	
			#v3新增，肽段的统计
			$pepseq_all=$datatype_loc2pepseq{$data_type}{$loc};  
			undef @arr1;
			@arr1=split(",", $pepseq_all);
			undef $uni_pepseq;
			for ($j=0; $j<@arr1; $j++)
			{
				$pepseq=$arr1[$j];
				if ($pepseq)
				{
					if (!($uni_pepseq))
					{
						$uni_pepseq="$pepseq".'('."$pepseq2specnum{$pepseq}".')';
					}
					elsif ($uni_pepseq!~/\b$pepseq\b/)
					{
						$uni_pepseq.=",$pepseq".'('."$pepseq2specnum{$pepseq}".')';
					}
				}
			}
			$uni_pepseq=0 if (!($uni_pepseq));
			print OUT2 "\t$uni_pepseq";
		}
		print OUT "\n";
		print OUT2 "\n";
	}
	close OUT;
	close OUT2;
}
		