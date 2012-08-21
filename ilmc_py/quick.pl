
open(IN,"pydata2.m");
open(OUT,">pdmod.m");

while(<IN>)
{
	$line = $_;
	@arr = split(',',$line);
	@arr2 = ();
	$i = 0;
	foreach $el (@arr)
	{
		

		if($el =~ m/nan/)
		{
			$el2 = "nan";
			push(@arr2,$el2)
		}
		elsif ($el =~ m/(-?)(\d\.\d{3})\d+(e.+)/)
		{
			$el2 = "$1$2$3";
			push(@arr2,$el2);
		}	

		if ($i < 10)
		{
			
			print "$el\t\t$el2\n";
		}
		$i++;
	}

	$line2 = join(",",@arr2);
	print OUT "$line2\n";
}

close(IN);
close(OUT);
