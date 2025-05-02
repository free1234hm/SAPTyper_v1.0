
$infile=$ARGV[0];
$search_engine=$ARGV[1];
$search_engine=1 if (!($search_engine)); #如果没有输入，默认是pFind

open (IN, "$infile")||die;

$title_marker=1;
if ($search_engine == 1)  #pFind
{
	while ($in=<IN>)
	{
		chomp $in;
		if ($in=~/^\#	Title/)
		{
			$title=$in;
		}
		else
		{
			undef @arr1;
			@arr1=split("\t", $in);
			$sample_tmp=$arr1[1];
			$sample="";
			if ($sample_tmp=~/^(.+?)\.[0-9]+\.[0-9]+\.[0-9]+\.[0-9]+\.dta$/)
			{
				$sample=$1;
			}
			$sample=~s/^(.+?)\..+/$1/;
			print "Sample error: =$sample=\n" if (!($sample));
			$sample2lines{$sample}.="$in\n";
		}
	}
}
elsif ($search_engine == 2)  #Maxquant
{
	while ($in=<IN>)
	{
		chomp $in;
		undef @arr1;
		@arr1=split("\t", $in);
		if ($in=~/^Sequence\t/)
		{
			$title=$in;
			for ($i=0; $i<@arr1; $i++)
			{
				$colname=$arr1[$i];
				if ($colname eq 'Raw file')
				{
					$rawfile_col=$i;
				}
			}
		}
		else
		{
			$sample_tmp=$arr1[$rawfile_col];
			$sample=$sample_tmp;
			print "Sample error: =$sample=\n" if (!($sample));
			$sample2lines{$sample}.="$in\n";
		}
	}
}
elsif ($search_engine == 3)  #DIA spectranaut searching results
{
	while ($in=<IN>)
	{
		chomp $in;
		undef @arr1;
		@arr1=split("\t", $in);

		if ($title_marker==1)
		{
			$title=$in;
			$title_marker=0;
			for ($i=0; $i<@arr1; $i++)
			{
				$str_tmp=$arr1[$i];
				if ($str_tmp eq 'R.FileName')
				{
					$sample_col=$i;
				}
			}
		}
		else
		{
			$sample_tmp=$arr1[$sample_col];
			$sample="";
			if ($sample_tmp=~/^(.+)_Slot/)
			{
				$sample=$1;
			}
			#$sample=~s/^(.+?)\..+/$1/;
			print "Sample error: =$sample=\n" if (!($sample));
			$sample2lines{$sample}.="$in\n";
		}
	}
}
elsif ($search_engine == 4)  #simple format
{
	while ($in=<IN>)
	{
		chomp $in;
		undef @arr1;
	#	print "$in\n";
		@arr1=split("\t", $in);
		if ($in=~/^Sample\t/)
		{
			$title=$in;
		}
		else
		{
			$sample_tmp=$arr1[0];
		#	print "sample1=$arr1[0]\n";
			$sample=$sample_tmp;
			$sample=~s/^([0-9a-z\_]+).*/$1/i;
		#	print "sample=$sample\n";
			print "Sample error: =$sample=\n" if (!($sample));
			$sample2lines{$sample}.="$in\n";
		}
	}
}

close IN;

open (OUT_dir, ">sep.result.list.tmp")||die;
while (($key, $value)=each(%sample2lines))
{
	$sample=$key;
	open (OUT, ">proteome.$sample.sep.txt")||die;
	print OUT_dir "proteome.$sample.sep.txt\n";
	print OUT "$title\n";
	print OUT "$value";
	close OUT;
}
close OUT_dir;
