use FileHandle;
use Getopt::Long;
$counter = 0;
@current = ();
@raw = ();
@residues = ();
@tmp_line=();
$x = 0;
GetOptions('query=s' => \$query, 'data=s' => \$data, 'length=i' => \$length);
open(SEQ, "$query") || die "couldn't open $query";
while ($line = <SEQ>)
{
	$tmp_seq = '';
	while ($seq_line = <SEQ>)
	{    
		    if (!($seq_line =~ />/))
		    {
			$tmp_seq = ($tmp_seq . $seq_line)
		    }
		    	else {last;}
	}	
	$query_seq = $tmp_seq;
	last;  
}    
chomp $query_seq;
@raw = split(//, $query_seq);
foreach (@raw)
{
	if ($_ =~ /[a-zA-Z]/)
	{
		$residues[$x] = $_;
		$x++;
	}
}
close SEQ;
open(SIMBAL, "$data") || die "couldn't open $data";
while ($line = <SIMBAL>)
{

    
    @split_line = split(/\t/,$line);
    for($f=0; $f<=7; $f++)
    {	
	  if ($split_line[1] == $length)
	  {	
		chomp $split_line[$f];		
		$correct[$counter][$f] = $split_line[$f];
#		print "just added $split_line[$f] to correct[$counter][$f]\n";
	  }	
    }
    if ($split_line[1] == $length) {$counter++;}
}
close SIMBAL;

for ($i=0; $i<@residues; $i++)
{
	if (@current >= $length)
	{
		$foo = pop @current;
#	print "just removed $foo from current\n";
	}
	if ($correct[$i][1])
	{
		unshift @current, $correct[$i][3];
#	print "just added $correct[$i][3] to current\n";
	}
	if (@current) 
	{

		$total = 0;
		for($w=0; $w<=@current; $w++)
		{
#			print"element $w of current is $current[$w]\n";
			$total = $total + $current[$w];
		}
		$average = $total/@current;
		$position = $i + 1;
		print "$residues[$i]\t$position\t$average\n";
	}
}

