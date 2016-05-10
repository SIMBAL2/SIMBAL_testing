
## Program:       profiles_to_site_scripter.pl
## Author:        Daniel Haft,  JCVI
## Creation Date: March 19, 2008
##
## program to generate the scripts that run sliding windows against training sets

use strict;
use warnings;
## use lib '/usr/local/devel/ANNOTATION/CEPHES_math_lib/lib/5.8.8/';     # math for binomial, 32 bit.
## use lib '/usr/local/projects/TUNE/malay/perl5lib/lib';                        # math for binomial, 64 bit.
use lib '/usr/local/packages/perl-5.16.1/lib/5.16.1';
use Math::Cephes qw(:dists);         # sets up complemented binomial distribution function:    bdtrc
my $DEBUG = 0;                       # 1 to print diagnostic progress reports, 0 to suppress

use Getopt::Std;          # for getting options
my %args;                 # receive getopts arguments
my $usage = "
USAGE:

REQUIRED
 -t [  ]  T-RUE:  FILENAME for fasta_file of TRUE  portion of training set
 -f [  ]  F-ALSE: for fasta_file of FALSE portion of training set
 -s [  ]  s-ingle sequence file, fasta format

OPTIONS
 -h               print this message and exit
 -j [  ]  j-ump : step-size for progressively changing window       DEFAULT = 10
 -m [  ]  m-inimum size for sliding window                          DEFAULT = 10
 -M [  ]  M-AXIMUM size for sliding window                          DEFAULT = 10
 -w [  ]  w-alk : the step size for sliding the window along
 -n [  ]  n-ame for project
 -x [  ]  scratch (X) directory path
 -p [  ]  probability : option to force use of probability p, rather than deriving from training set
 -z       zap :  don't delete the BLAST file after using it (default behavior is to delete)
 -E [  ]  E-value figure for blastp - default of 10.0 invites noise contamination of results.
";

my ($true_name, $false_name, $jump, $walk, $mini, $maxi, $tester, $odds, $zap, $expect);

# PROCESS OPTIONS:
&getopts('ht:f:j:m:M:s:w:n:x:p:zE:',\%args);
if ($args{'h'}) {
    die $usage;              # print the help message.
}

if ($args{'t'}) {
    $true_name =$args{'t'} ;
} 
else {
    die "Supply filename of TRUE sequences (-t) in fasta format, or get help (-h))\n";
}
if (! -e $true_name)  {
    die "fasta file of true-positives (-t option) not located from path given\n";
}

if ($args{'f'}) {
    $false_name =$args{'f'} ;
} 
else {
    die "Supply filename of FALSE sequences (-f) in fasta format, or get help (-h))\n";
}
if (! -e $false_name)  {
    die "fasta file of false-positives (-f option) not located from path given\n";
}

if ($args{'s'}) {               # sequence file to test
    $tester = $args{'s'} ;
} 
else {
    die "supply filename of a Sequence (-s) , or get help (-h))\n";
}
if (! -e $tester)  {
    die "fasta file of sequence to test, $tester,  not located from path given\n";
}

if ($args{'j'}) {
    $jump = $args{'j'};  
    if ($jump < 1) { 
        die "option (-j) for jump in window size must be > 1\n";
    }
} 
else {
    $jump = 10;
    print "Using default:  ";
}
print "window size progressing by jumps of $jump\n";

if ($args{'w'}) {
    $walk = $args{'w'};  
    if ($walk < 1) { 
        die "option (-w), step size for walking along, must be > 1\n";
    }
} 
else {
    $walk = 10;
    print "Using default:  ";
}
print "sliding window along in steps of size $walk\n";

if ($args{'m'}) {
    $mini = $args{'m'};  
    if ($mini < 5) { 
        die "a minimum window smaller than 5 makes no sense. \n";
    }
} 
else {
    $mini = 10;
    print "Using default:  ";
}
print "STARTING with window size $mini\n";

if ($args{'M'}) {
    $maxi = $args{'M'};  
    if ($maxi < 5) { 
        die "a window smaller than 5 makes no sense. \n";
    }
} 
else {
    $maxi = 10;
    print "Using default:  ";
}
print "ENDING with window no larger than $maxi\n";

if ($args{'z'}) {
    $zap = 0
} 
else {
    $zap =1
}           #  z for ZAP: delete blast hits file

if ($args{'p'}) {      # probability p for use in binomial distribution.  Hand-set p can bias for greater stringency.
    $odds = $args{'p'};  
    if ($odds < 0 || $odds > 1) { 
        die "Odds outside range 0-1 makes no sense \n";
    }
} 
else {
    $odds = -1;        # flag to use derived odds instead of pre-set
}
if ($args{'E'}) {
    $expect = $args{'E'};  
    if ($expect <= 0) {
         die "QUITTING: an E (expect) value must be positive!\n";
    }
} 
else {
    $expect = 10;
    print "NOTICE: default E-value of 10 invites noise-contaminated results. We recommend 0.1 \n";
}
print "STARTING with window size $mini\n";


#  Finished processing options. Check paths and output files

my ($scratch_path, $proj_name);   #DEFAULTS  "/usr/local/scratch/PPP_dir/" and "TEST_"
if ($args{'x'}) {$scratch_path = $args{'x'}} else {$scratch_path = "/usr/local/scratch/"}
if ($args{'n'}) {$proj_name = $args{'n'}} else {$proj_name = "DELETE_ME"} 
my $scratch_blast = $scratch_path . $proj_name . "_BASE";
my $report_file = $proj_name . "_REPORT";
#if (-e $report_file) {die "Sorry, I won't let you overwrite existing report file $report_file"}
open(SUMMARY,">$report_file") || die "Can't open report file $report_file for writing: quitting";


#  NowRead in the three sequence files
my ($line, $count, $true_count, $false_count, @seq, @truelib, @falselib); 
$count = &read_fasta($tester, \@seq);
if ($count != 1) {
    die "Wrong number, $count, of sequences in test sequence file $tester\n"
}
$true_count  = &read_fasta($true_name,  \@truelib);
$false_count = &read_fasta($false_name, \@falselib);
my $p;
if ($odds > 0) {
    $p = $odds;       # value was manually set, option -p
    print "Probability is $p, overriding calculation based on $true_count positives, $false_count  negatives\n";
} 
else {      
    $p = $true_count / ($true_count + $false_count);
    print "Probability is $p, based on $true_count positives, $false_count  negatives\n";
}

if (-e $scratch_blast) {    # delete the searchable scratch file, if it exists, then write in append mode
    unlink $scratch_blast
}
&append_fasta("TRUE_", $scratch_blast, \@truelib);
&append_fasta("FALSE_",$scratch_blast, \@falselib);
print "/usr/local/packages/ncbi-blast+/bin/makeblastdb -in $scratch_blast -dbtype prot\n";
system("/usr/local/packages/ncbi-blast+/bin/makeblastdb -in $scratch_blast -dbtype prot") && die "It ain't work.";                      

#  remember: window sizes from $mini to $maxi, jumping window size by $jump. Window slides by $walk
my $seq_len = length($seq[0]{seq});
my ($win_size, $j);                      # $j is array position of start of window

my ($sample, $blastout, $log_best, $midpoint, $row, $col);
my $search_me = "temp_" . $proj_name;
my @stats;                 # array for getting back parse results
$log_best = 0;

for ($win_size=$mini; $win_size<=$maxi; $win_size+=$jump) {
    $row++;     # allow for tracking row and column with simple numbering - allows for easier optional post-processing
    $col=0;
    for ($j=0; $j<=$seq_len-1-$win_size+$walk; $j+=$walk) {
    	$col++;
    	$sample = substr($seq[0]{seq}, $j, $win_size);
    	open(SEARCHSTRING,">$search_me");
    	print SEARCHSTRING (">",$proj_name,"_",$win_size,"_",$j,"\n",$sample,"\n"); 
    	close (SEARCHSTRING);
    	$midpoint = $win_size/2 + $j + 1/2;     # We number sequence from 1, but Perl counts from 0.
    	$blastout = "OUT_" . $proj_name . "_" . $win_size . "_" .$midpoint;
    #	print  "blastp $scratch_blast  $search_me -span1 -E $expect > $blastout\n";  # span1 protects vs. duplication
    #       system "blastp $scratch_blast  $search_me -span1 -E $expect > $blastout";    # span1 protects vs. duplication
    	print 	"/usr/local/packages/ncbi-blast+/bin/blastp -query $search_me  -db $scratch_blast -outfmt 7\n";
    	system "/usr/local/packages/ncbi-blast+/bin/blastp -query $search_me -db $scratch_blast -outfmt 7 > $blastout";
    	print "$blastout\n";
    	&parse_blast ($blastout, $p, \@stats);
    	
    	if ($stats[0] != 0) { 
    	    $log_best = -log($stats[0])/2.30258509
    	}
    	$log_best = 0 unless (defined $log_best);   ## new to try to make printf error go away
    	print SUMMARY ("$proj_name\t$win_size\t$midpoint\t$log_best\t$stats[1]\t$stats[2]\t$stats[3]\t$row\t$col\n");
    	
    	if ($zap) {
    	    unlink $blastout
    	}
    } 
}

# program parses WU-BLAST  blastp  output
sub parse_blast  {

    my ($blastfile, $odds, $stats_ref) = @_;
    my ($line, $yes_ctr, $no_ctr, $current, $best, $yes_at_best, $no_at_best, $temp_yes, $temp_no, @split_line, $current_eval, $current_bit);
    $best = 1;            # before taking negative log, so smaller is better.
    $current = 1;         # required if we want our time-saving check later to not cause an error
    ($yes_ctr, $no_ctr, $yes_at_best, $no_at_best) = (0) x 4;
    ($temp_yes, $temp_no) = (0) x 2;
    $current_eval = $current_bit = 0;
    @split_line = '';

    open(READ_BLAST, $blastfile) || die "Can't read supposed blast output file\n";

    while ($line  = <READ_BLAST>) {
    	if ($line =~ /^\#/) {
    	    next;	
    	} 
    	else {
    	    @split_line = split(/\s+/, $line);  
    	    
    	    if ($line =~ /\sTRUE_/) {
        		$temp_yes++;    			    	
    	    } 
    	    elsif ($line =~ /\sFALSE_/) {
        		$temp_no++;
    	    } 
    	    else {
    		  last;          # exit while loop
	       }
    	}
        if ($split_line[10] == $current_eval && $split_line[11] == $current_bit) {
    	    next;
    	} 
    	else  {
    	    $current_eval = $split_line[10];
    	    $current_bit = $split_line[11];
    	    $yes_ctr += $temp_yes;	   
    	    $no_ctr += $temp_no;
    	    if ($temp_yes) {	                                     
    		  $current = do_binomial($yes_ctr, $yes_ctr + $no_ctr, $p);	
    	    } 
    	    $temp_no = 0;
    	    $temp_yes = 0;
    	}
    	if ($current < $best) {
    	    $best = $current; $yes_at_best = $yes_ctr; $no_at_best = $no_ctr;
    	} 
    }
    @$stats_ref = ($best, $yes_at_best, $no_at_best, $p);
    print ("best = $best, yes_at_best = $yes_at_best, no_at_best = $no_at_best\n");
	
    close(READ_BLAST);
}

sub append_fasta {
    my ($label, $blastfile, $seq_array_ref) = @_;
    open(BLASTABLE,">>$blastfile") || die "Can't write $blastfile";
    
    foreach my $guy (@$seq_array_ref) {
    	$guy->{seq} =~ s/\s*//g;        # remove all spaces from sequence
    	$guy->{seq} =~ s/(.{60})/$1\n/g;
    	$guy->{seq} =~ s/\n\Z//;        # if multiple of 60, remove extra return at end of sequence
    	print BLASTABLE (">",$label,"$guy->{acc}\n$guy->{seq}\n"); 
    }  
}

sub read_fasta {
    my $ctr = 0;
    my ($filename, $seq_array_ref) = @_;  
    open (FASTAFILE, $filename) || die "Can't open sequence file $filename\n";
    my ($line, $accession, $name);
    
    while ($line = <FASTAFILE>) {
    	chomp ($line);
    	if ($line =~ /^\>/) {
    	    $ctr++;
    	    $line =~ s/\>//;
    	    $accession = $name = $line;
    	    $accession =~ s/\s.*//;       # delete everything after first space
    	    $name =~ s/\S*\s//;           # delete everyting up through first space
    	    $seq_array_ref->[$ctr-1]{acc}  = $accession;
    	    $seq_array_ref->[$ctr-1]{name} = $name;
    	} 
    	else {
    	    $seq_array_ref->[$ctr-1]{seq} .= $line;   
    	}
    }
    close (FASTAFILE);
    $seq_array_ref->[$ctr-1]{seq} =~ s/\s//g;     #stripping whitespace: should force to uppercase and check
    $DEBUG &&  print "> ACCESSION $accession :::  NAME $name\nSEQUENCE :: $seq_array_ref->[$ctr-1]{seq}\n ";
    return $ctr;
}

sub do_binomial {
    my ($y, $n, $p) = @_;         # $n means total, not negatives
    my $answer;
    if ($y<1) {
        $answer = 1
    } 
    else {
        $answer = bdtrc($y-1, $n ,$p)
    }  # no point calculating if no yesses.
 #   print ("I did a binomial! : $y - yes. $n - no. $p - prob. $answer\n");
    return $answer;
}

##  This program makes certain assumptions:
##
##  1.   access to WU-BLAST program "blastp" in your path: see  http://blast.wustl.edu/
##     (if you use NCBI-BLAST instead, options and parsing will differ some)
##
##  2.   The CEPHES math library, necessary for the complemented binomial distribution cunction  bdtrc
##
