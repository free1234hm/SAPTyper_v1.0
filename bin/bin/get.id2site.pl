#此程序只用于产生ID和位点，并进行排序
#另外加一列，包含染色体和定位，并把相同染色体和定位的结果给合并了（避免因为密码子简并性形成多个行号，导致后期识别成多个位点）

$infile=$ARGV[0];
#open (IN, "ExAC.EAS001.v2.exonic_variant_function")||die;
open (IN, "$infile")||die;
open (OUT, ">$infile.format.txt")||die;

$k=0;
while ($in=<IN>)
{
	chomp $in;
	undef @tmp_array;
	@tmp_array=split("\t", $in);
	$line=$tmp_array[0];
	$anno=$tmp_array[2];
	$freq=$tmp_array[8];
	$chr_id=$tmp_array[3];
	$start_loc=$tmp_array[4];
	$ref_nt=$tmp_array[6];
	$alt_nt=$tmp_array[7];
	#$loc_info="chr$chr_id\@$start_loc";
	$var_info="chr$chr_id"."_"."$start_loc$ref_nt".'>'."$alt_nt";
	undef @tmp_array2;
	@tmp_array2=split(",", $anno);

	$pronum=@tmp_array2;
	for ($i=0; $i<$pronum; $i++)
	{
		$proanno=$tmp_array2[$i];
		if ($proanno=~/(ENST.+?):/)
		{
			$id=$1;
		#	print "$id\n";
		}
		if ($proanno=~/:p\.([a-z]+)([0-9]+)([a-z]+)/i)
		{
			$refaa=$1;
			$site=$2;
			$altaa=$3;
			next if ($refaa eq $altaa);
			$id2site[$k]="$id	$site	$refaa	$altaa	$freq	$line	$var_info";
			$k++;
		}
	}
}

@sort_id2site=sort{$a cmp $b}@id2site;

for ($i=0; $i<$k; $i++)
{
	undef @arr1;
	@arr1=split("\t", $sort_id2site[$i]);
	$info="$arr1[0]\t$arr1[1]\t$arr1[2]\t$arr1[3]";
	if ($info eq $info_old)
	{
		$freq+=$arr1[4];
		$line.=",$arr1[5]";
		#$loc_info.=",$arr1[6]";
		$var_info.=",$arr1[6]";
	}
	else
	{
		if ($info_old)
		{
			print OUT "$info_old\t$freq\t$line\t$var_info\n";
		}
		$freq=$arr1[4];
		$line=$arr1[5];
		#$loc_info=$arr1[6];
		$var_info=$arr1[6];
		$info_old=$info;
	}
	#print OUT "$sort_id2site[$i]\n";
}
print OUT "$info_old\t$freq\t$line\t$var_info\n";
