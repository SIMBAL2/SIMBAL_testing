use Math::Round;
use FileHandle;
use Getopt::Long;

my @old = '';
my @weights = '';
my @tmp_line = '';
my @output = '';
my $score = 0;
my $h = 0.93;                                 # h for heritable - can be adjusted - is the percent that a row inherits from previous row
my @split_line;
my $debug = 0;
my $entry = 0;

GetOptions('input=s' => \$input, 'output=s' => \$output, 'debug' => \$debug);

open(INPUT, "$input") || die "couldn't open $input";
while(<INPUT>)
{
    next if (/^$/);   
    @split_line = split(/\t/,$_);
    for($f=0; $f<=8; $f++)
    {
	    $tmp_line[$f] = $split_line[$f];
    }
    for($f=0; $f<=8; $f++)
    {
	chomp($tmp_line[$f]);	    
	$old[$entry][$f] = $tmp_line[$f];
    }
#    if ($debug) {print "entry = $entry, old[entry][3] = $old[$entry][3]\n";}    
	$entry++;
}
close INPUT;
assign_scores();
open (OUTPUT, ">>$output");
for ($i=0; $i<=$entry; $i++)
{
	$row = $old[$i][7];
	$col = $old[$i][8];	
	print OUTPUT ("$old[$i][0]\t$old[$i][1]\t$old[$i][2]\t$weights[$row][$col]\t$old[$i][4]\t$old[$i][5]\t$old[$i][6]\t$old[$i][7]\t$old[$i][8]\n");
#	if ($debug) {print "($row, $col) $weights[$row][$col]\n";}
}
close OUTPUT;

sub assign_scores
{
	for ($x=$entry; $x>=$0; $x--)                       # go through entries in reverse order (longest sequences first)
	{
		if ($debug) {print ("entry = $x, $old[$x][0], $old[$x][1], $old[$x][2], $old[$x][3]\n");}		
		$row = $old[$x][7];
		$col = $old[$x][8];
		if ($weights[$row+1][$col-1] && $weights[$row+1][$col])
		{
			$parent_one = $weights[$row+1][$col-1];
			$parent_two = $weights[$row+1][$col];
			$base = $old[$x][3];
			$score = nearest(.001, (.5 * $h * ($parent_one + $parent_two)) + ((1-$h)*$base));

			$weights[$row][$col] = $score;
			if ($debug) {print "case 1 ($row, $col): p1 = $parent_one p2 = $parent_two base = $base score = $score\n";}
			if ($debug) {print "x = $x, old[x][1] = $old[$x][1], old[x][2] = $old[$x][2], old[x][3] = $old[$x][3]\n";}
		}
		elsif ($weights[$row+1][$col-1] && !$weights[$row+1][$col])
		{
			$parent_one = $weights[$row+1][$col-1];
			$parent_two = 0;
			$base = $old[$x][3];
			$score = nearest(.001, ($h * ($parent_one + $parent_two)) + ((1-$h)*$base));
			$weights[$row][$col] = $score;
			if ($debug) {print "case 2 ($row, $col): p1 = $parent_one p2 = $parent_two base = $base score = $score\n";}
		}

		elsif (!$weights[$row+1][$col-1] && $weights[$row+1][$col])
		{
			$parent_one = 0;
			$parent_two = $weights[($row + 1)][$col];
			$base = $old[$x][3];
			$score = nearest(.001, ($h * ($parent_one + $parent_two)) + ((1-$h)*$base));
			$weights[$row][$col] = $score;
			if ($debug) {print "case 3 ($row, $col): p1 = $parent_one p2 = $parent_two base = $base score = $score\n";}
		}
		else 
		{
			$score =  nearest(.001, $old[$x][3]);
			$weights[$row][$col] = $score;
			if ($debug) {print "case 4 ($row, $col): p1 = $parent_one p2 = $parent_two base = $base score = $score\n";}
		}
	}
}















