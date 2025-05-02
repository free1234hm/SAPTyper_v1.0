
#利用位点变异信息和蛋白质序列，构建两类变异数据库：
#1. 一个变异一条蛋白序列
#2. 一条蛋白一条变异序列

#20210610：碰到突变为终止密码子时，也跳过，因为后面处理会出错

#$var_file="id2site.line.v3.txt";
$var_file=$ARGV[0];
#$pro_db='l:\work\SAPproject\bin\coding.protein.clean.fa';
$pro_db=$ARGV[1];
open (IN_var, "$var_file")||die;
#open (IN_var, "id2site.tmp.txt")||die;
open (IN_db, "$pro_db")||die;

while ($in=<IN_db>)
{
	chomp $in;
	if ($in=~/^>([^ ]+)/)
	{
		$id=$1;
	}
	elsif ($in)
	{
		$seq=$in;
		$seq=~s/\*$/X/;
		$id2seq{$id}=$seq;
	}
}

$i=0;
while ($in=<IN_var>)
{
	chomp $in;
	$id2site[$i]=$in;
	$i++;
}
undef @sort_id2site;
@sort_id2site=sort{$a cmp $b}@id2site;

open (OUT_seq_all, ">all.seq.all.fa")||die;
open (OUT_anno_all, ">all.anno.all.txt")||die;
open (OUT_seq1, ">all.seq.singleSAP.fa")||die;
open (OUT_anno1, ">all.anno.singleSAP.txt")||die;


for ($i=0; $i<@sort_id2site; $i++)
{
	$pro_id=$sort_id2site[$i];
	$pro_id=~s/^(.+?)\t.+/$1/;
	if (($pro_id_old ne $pro_id)) #新的蛋白出现
	{
		#print "pro=$pro_id_old\n";
		if ($pro_id_old)  #需要对前一个蛋白操作
		{
			$seq=$id2seq{$pro_id_old};
			$digest_site=$id2digestsite_hash{$pro_id_old};
			@sub_results=SingleProAna($pro_id_old,$seq,$digest_site,\@varsite_singlePro,$freq_cutoff);
			@seq_fa1=@{$sub_results[0]};
			@var_anno1=@{$sub_results[1]};
			@seq_fa_all=@{$sub_results[2]};
			@var_anno_all=@{$sub_results[3]};
			for ($k=0; $k<@seq_fa_all; $k++)
			{
				print OUT_seq_all "$seq_fa_all[$k]";
			}
			for ($k=0; $k<@var_anno_all; $k++)
			{
				print OUT_anno_all "$var_anno_all[$k]";
			}
			for ($k=0; $k<@seq_fa1; $k++)
			{
				print OUT_seq1 "$seq_fa1[$k]";
			}
			for ($k=0; $k<@var_anno1; $k++)
			{
				print OUT_anno1 "$var_anno1[$k]";
			}
		}
		undef @varsite_singlePro;
		$j=0;
	}
	$varsite_singlePro[$j]=$sort_id2site[$i];
	$j++;
	$pro_id_old=$pro_id;
}
#处理最后一个蛋白的结果
#print "pro=$pro_id_old\n";
if ($pro_id_old)  #需要对前一个蛋白操作
{
	$seq=$id2seq{$pro_id_old};
	$digest_site=$id2digestsite_hash{$pro_id_old};
	@sub_results=SingleProAna($pro_id_old,$seq,$digest_site,\@varsite_singlePro,$freq_cutoff);
	@seq_fa1=@{$sub_results[0]};
	@var_anno1=@{$sub_results[1]};
	@seq_fa_all=@{$sub_results[2]};
	@var_anno_all=@{$sub_results[3]};
			for ($k=0; $k<@seq_fa_all; $k++)
			{
				print OUT_seq_all "$seq_fa_all[$k]";
			}
			for ($k=0; $k<@var_anno_all; $k++)
			{
				print OUT_anno_all "$var_anno_all[$k]";
			}
			for ($k=0; $k<@seq_fa1; $k++)
			{
				print OUT_seq1 "$seq_fa1[$k]";
			}
			for ($k=0; $k<@var_anno1; $k++)
			{
				print OUT_anno1 "$var_anno1[$k]";
			}
	#for ($k=0; $k<@seq_fa; $k++)
	#{
	#	print OUT_seqall "$seq_fa[$k]";
	#}
	#for ($k=0; $k<@var_anno; $k++)
	#{
	#	print OUT_annoall "$var_anno[$k]";
	#}
}

close OUT_seq_all;
close OUT_anno_all;
close OUT_seq1;
close OUT_anno1;

#以下是单个蛋白的分析程序，把这单个蛋白的分析程序做成一个子程序，就可以不用改多少了
sub SingleProAna
{
	my ($proid,$seq,$digest_site,$varsite,$freq_cutoff) = @_;
	my @varsite=@$varsite;
	#以ENST00000531678蛋白为例（有32个变异位点，算比较多的例子了，而且比较密集）
	#$seq='MITFLPIIFSSLVVVTFVIGNFANGFIALVNSIEWFKRQKISFADQILTALAVSRVGLLWVLLLNWYSTVLNPAFNSVEVRTTAYNIWAVINHFSNWLATTLSIFYLLKIANFSNFIFLHLKRRVKSVILVMLLGPLLFLACHLFVINMNEIVRTKEFEGNMTWKIKLKSAMYFSNMTVTMVANLVPFTLTLLSFMLLICSLCKHLKKMQLHGKGSQDPSTKVHIKALQTVISFLLLCAIYFLSIMISVWSFGSLENKPVFMFCKAIRFSYPSIHPFILIWGNKKLKQTFLSVFWQMRYWVKGEKTSSP';
	#$digest_site='37,38,40,55,81,109,122,123,124,126,154,156,165,167,169,204,207,208,214,222,226,265,268,284,285,287,298,302,305';
	#open (IN_var_site, "example.id2site.txt")||die;

	my $i;
	my $in;
	my $site_old;
	my @site_array;
	my @site_detail_array;
	my $ver;
	#定义两个变量，分别为单个变异的，和所有变异的
	my @varnamevercomb1;
	my @varnamevercomb_all;
	my %site_var_hash;
	my $i1=0;
	#while ($in=<IN_var_site>)
	for ($i=0; $i<@varsite; $i++)
	{
		$in=$varsite[$i];
		undef @tmp_array;
		@tmp_array=split("\t", $in);
		#$proid=$tmp_array[0];
		$tmp_array[4]=1 if ($tmp_array[4]!~/[0-9\.]/);  #加了一列，最后一列是频率。如果没提供频率数据，就用1，对结果不影响
		#print "freq: $tmp_array[4]\n";
		next if ($tmp_array[4]<$freq_cutoff);   #如果单个位点的频率低于卡值，直接跳过
		next if ($tmp_array[3] eq 'fs');  #如果是移码突变，跳过
		next if ($tmp_array[3] eq 'X');  #如果是突变为终止密码子，因为在蛋白质组很难识别出变异，也跳过
		$tmp_array[3]=~s/delins//;   #把这字符删除
		#print "freq2: $tmp_array[4]\n";
		$varnamevercomb1[$i1]="$tmp_array[2]$tmp_array[1]$tmp_array[3]";
		$i1++;
		$site_tmp=$tmp_array[1];
		$site_freq=$tmp_array[4];
		if ($site_freq>$site_freq_hash{$site_tmp})
		{
			$site_var_hash{$tmp_array[1]}="$tmp_array[2]$tmp_array[1]$tmp_array[3]";
			$site_freq_hash{$site_tmp}=$site_freq;
		}
	#	$i++;
	}
	$i=0;
	while (($key, $value)=each(%site_var_hash))
	{
		$num_tmp=$value;
		$num_tmp=~s/[a-z]//gi;  #把数字提出来，把字母踢掉。加入数组中，才能真正排序
		$site_array[$i]="$num_tmp\t$value";
		$i++;
	}

	my @sort_site_array;
	@sort_site_array="";
	#print "=@site_array=\n";
	@sort_site_array=sort{$a<=>$b}@site_array;
	for ($i=0; $i<@sort_site_array; $i++)
	{
		$sort_site_array[$i]=~s/^.+?\t//;
		$varnamevercomb_all[0].="$sort_site_array[$i],";
	}
	$varnamevercomb_all[0]=~s/,$//;
	#print "==@sort_site_array==\n";


	#my @varnamevercomb=GetVarComVer(\%SiteCombination,\%N0_dig_site,\%Nterm_dig_site,\@site_detail_array);  #产生每个版本序列的变异度组合---#把这改了就行，改为每条序列一个，或者每个位点一个
	my @verseq1=GetVarSeq($seq,\@varnamevercomb1);  #通过变异的组合，以及序列，获得各个版本的变异序列结果
	my @verseq_all=GetVarSeq($seq,\@varnamevercomb_all);  #通过变异的组合，以及序列，获得各个版本的变异序列结果

	my @seq_fa1;
	my @var_anno1;
	for ($i=0; $i<@verseq1; $i++)
	{
#		print OUT "seq$i: $verseq[$i]\n";
#		print OUT_seq ">$proid.v$i alt pro\n$verseq[$i]\n";
		$verseq1[$i]=~s/X.*//;   #由于stop会被定义为X，因此，X后面的删掉。然后判断是否还有序列，没有序列就跳过
		if ($verseq1[$i])
		{
			$seq_fa1[$i]=">$proid.sv$i alt pro\n$verseq1[$i]\n";
			$var_anno1[$i]="$proid.sv$i\t$varnamevercomb1[$i]\n";
	#		print OUT_anno "$proid.v$i\t$varnamevercomb[$i]\n";
		}
	}
	my @seq_fa_all;
	my @var_anno_all;
	for ($i=0; $i<@verseq_all; $i++)
	{
#		print OUT "seq$i: $verseq[$i]\n";
#		print OUT_seq ">$proid.v$i alt pro\n$verseq[$i]\n";
		$verseq_all[$i]=~s/X.*//;   #由于stop会被定义为X，因此，X后面的删掉。然后判断是否还有序列，没有序列就跳过
		#$verseq_all[$i]=~s/X.*//;
		if ($verseq_all[$i])
		{
			$seq_fa_all[$i]=">$proid.av$i alt pro\n$verseq_all[$i]\n";
			$var_anno_all[$i]="$proid.av$i\t$varnamevercomb_all[$i]\n";
	#		print OUT_anno "$proid.v$i\t$varnamevercomb[$i]\n";
		}
	}
#	close OUT_seq;
#	close OUT_anno;

	#print "==%SiteCombination\n";
#	my $site_num=@site_array;
#	my $var_site;
#	for ($i=0; $i<$site_num; $i++)  #逐个位点遍历
#	{
#		$var_site=$site_array[$i];
#	#	print OUT "$var_site:\t$SiteCombination{$var_site}\n";
#	}
	return (\@seq_fa1,\@var_anno1, \@seq_fa_all,\@var_anno_all);
}

sub GetVarSeq
{
	my ($seq,$varnamevercomb) = @_;
	my @varnamevercomb=@$varnamevercomb;
	my @var_arr;
	my $i;
	my $j;
	my $varname;
	my $seqlen=length($seq);
	my @seq_ver;
	for ($i=0; $i<@varnamevercomb; $i++)
	{
		@var_arr="";
		@var_arr=split(',',$varnamevercomb[$i]);
		$seq_tmp=$seq;  #每种组合重新取原序列
		for ($j=(@var_arr-1); $j>=0; $j--)  #倒序去做，以免因为有插入或缺失导致氨基酸位置计数出问题
		{
			$varname=$var_arr[$j];
			$varname=~s/-[0-9]+//;
			if ($varname=~/([a-z]+)([0-9]+)([a-z]+)/i)
			{
				$ref_aa=$1;
				$site=$2;
				$alt_aa=$3;
			#	if ($alt_aa eq 'fs') #如果碰到移码突变，则不产生该序列
			#	{
			#		$del_marker=1;
			#		next;
			#	}
				$alt_aa=~s/^delins//;  #碰到插入的情况，把这字符删除就行
				$reflen=length($ref_aa);
				$substr_aa=substr($seq_tmp,$site-1, $reflen);
				if ($ref_aa ne $substr_aa)
				{
					print "$pro_id_old\tError: ref-$ref_aa, substr-$substr_aa, alt-$alt_aa, site=$site, varComb=$varnamevercomb[$i]\n";
					print "$seq_tmp\n";
					print "$seq\n";
				#	next;  #碰到不匹配的情况，不处理，进行下一个。不管用！
				}
				$seq1=substr($seq_tmp,0,$site-1);
				$seq2=substr($seq_tmp,$site-1+$reflen);
			#	print "var:$varname\t$site\t$reflen\nSeq1=$seq1\nSeq2=$seq2\n";
				$seq_tmp="$seq1$alt_aa$seq2";
			}
			else
			{
				print "Error: $varname\n";
			}
		}
		$seq_ver[$i]=$seq_tmp;
	#	print "varname: $varnamevercomb[$i]\nseq=$seq_ver[$i]=--=\n";
	}
	return @seq_ver;
}
