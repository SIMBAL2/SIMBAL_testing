# TRAINING SET BUILDER BY DAVID RENFREW HAFT - COPYRIGHT 2015-2016
# ALL RIGHTS RESERVED
# SOFTWARE PUBLISHED UNDER THE GPL (GNU GENERAL PUBLIC LISCENSE) 

# This program makes use of cURL to download genomes from NCBI via FTP
# If you are only using already downloaded files, cURL is not required
# Otherwise, please install cURL before using this software

# This program requires HMMER3 (hmmer.org)
# LINES 27, 311, AND 312 SHOULD BE MODIFIED BASED ON THE CORRECT PATHS IN THE USER ENVIRONMENT
###################################################################################################################################################################
($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime($^T);
$mon++;	 #days count from 1, months count from 0
$year -= 100; #years count from 1900
$start_time = "$mon\/$mday\/$year $hour:$min:$sec\n";
###################################################################################################################################################################
use FileHandle;
use Getopt::Long;
use Carp;
###################################################################################################################################################################
my $keep = '';
my $problematic = 0;              # Is given a non-zero value if cURL fails a download
my $input_acc = '';
my $input_list = '';
my $hmm_dir =  "./";
my $config_file = "./training_set.config";
my $refseq_location = "/usr/local/projdata/9672/projects/DB/ProkCompleteRefSeq";	# Change to your own preferred location - this will be used as the default location to check for downloaded genomes
our @rules = ();
our @answers = ();
my $uses_coords = '';             # This variable gets set if any rules the user provides require information about coordinates


my $debug_rule = '';
my %targets =();
my %attribute = ();
my $debug = '';
my $usage = "
OPTIONS:
--help -> display this message and quit
--debug -> verbose runtime output and creation of logfiles
--keep -> don't delete downloaded genomes after processing (recommended)
--local -> program will look for and download genomes to the current directory 
--location <some file path> -> program will look at and download genomes to the location given here
(if no location is given, JCVI default is /usr/local/projdata/9672/projects/DB/ProkCompleteRefSeq)

REQUIRED:
--list -> input as list of accession numbers
AND/OR
--acc ->  single accession number
OR
-complete -> use all RefSeq completed genomes. if using this flag, do not combine with other forms of input

REQUIRED:
--config -> config file describing what sequences to collect and how to sort them - see README\n(defaults to looking for training_set.config in the current directory)
 
\n";


GetOptions ('help' => \$help,
	    'acc=s' => \$input_acc,
	    'keep' => \$keep,
	    'debug' => \$debug,
	    'list=s' => \$input_list,
	    'hmm_dir=s' => \$hmm_dir,
            'config=s' => \$config_file,
	    'location=s' => \$given_location,
            'local' => \$local,
	    'complete' => \$complete	);

mkdir("_TRAINING") unless (-e "_TRAINING");   # while the position of downloaded/read genomes is adjustable, output goes locally
mkdir("_YES") unless (-e "_YES");
mkdir("_NO") unless (-e "_NO");
mkdir("_TMP") unless (-e "_TMP");
mkdir("_JUNK") unless (-e "_JUNK");

			
if (!$input_list && !$input_acc && !$complete)
{
    carp ("!!!!!!!!!!!!!!ERROR, THIS PROGRAM REQUIRES A FORM OF INPUT!!!!!!!!!!!!!\n");
    $help = 1;
}

if ($complete && ($input_list || $input_acc))
{
	die ("the -complete flag is mutually exclusive with the -list/-acc flags\n");
}

if ($help)
{
    carp ($usage);
    die;
}
my ($target_hmm, $target_score) = parse_config();                    ## Read the config file made by user, and populate 2 arrays - @yes_rules and @no_rules
if ($debug) {print "|||$rules[1][0][0][1]     $rules[0][0][0][1]|||\n"};
if ($debug) {print "|||$rules[1][1][0][1]     $rules[0][1][0][1]|||\n"};

if ($complete && $local)                                             ## If the user wants to download a second copy of this DB, they need to do so very explicitly becuase they probably shouldn't be
{
	die ("The -complete and -local flags are mutually exclusive.\nIf you really want to download all RefSeq completed genes locally, specify the local path explicitly with the -db flag");
}

if ($input_list)
{
	@desired_genomes = read_list($input_list);                   ## If @desired_genomes in non-empty, it means we are going to need to download files
}
if ($input_acc)
{
	push(@desired_genomes, $input_acc);                          ## Doing it this way allows multiple input techniques to be combined in a single run
}

if ($complete)                                                       ## Download RefSeq genome list, analyze every genome marked Complete
{
	`curl ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/assembly_summary.txt > genome_list`;
	`grep "Complete" genome_list > complete_list`;
	@desired_genomes = read_list("./complete_list"); 	     ## This file (complete_list) can be used as a record of what genomes existed when the complete flag was invoked (assuming you don't overwrite it)
}


%paths = ();
 if ($local)
{
	  $location = ".";
}
elsif ($given_location)
{
	chomp($given_location);
	$location = "$given_location";
}	 
else
{
	$location = "$refseq_location";
}

if (@desired_genomes)                                                
{                                                                    ## most of the code happens inside this subroutine - which is admittedly not ideal
    ftp_download();                                                  
}                                                
if ($problematic)
{
	carp"Some genomes failed to download - their accessions are stored in a list in whatever directory you ran this program\n If you want, you can supplement this run by using that file as the input for the -list flag\n";
}                                                                                                                             
                                                                                                           
###################################################################################################################################################################
sub ftp_download
{
                                                                     # collect accessions of user-requested genomes
    `curl ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/assembly_summary.txt > genome_list`;                            
    mkdir("$location/ftp_downloads") unless (-e "$location/ftp_downloads"); 
    open(GENOMES, "./genome_list") || die("couldn't get genome information\n");	                       # not an error that this has . not $location - the FTP lookup file doesn't belong in the DB
    while($line = <GENOMES>)
    {
		if (!($line =~ /^#/))                                                         	       # ignore the header line
		{
		    @split_line = split(/\t/, $line);                                                  # splits on tabs - some entries contain whitespace
		       
			$paths{$split_line[0]}[0] = $split_line[19];                                   # creates a hash of ftp paths keyed to accession numbers
		}
    }
  
    foreach(@desired_genomes)                                                                          # some logic is by genome, some is by target sequence
    {
	 $downloaded = 0;                                                                              # even without -keep we won't delete anything we didn't download ourselves
	 $is_sortable = 0;                                                                             # set to 1 if a genomes has at least 1 target that is a YES or a NO
	 @answers = ();                                                                                # non-distance tests don't need to be repeated for each target - save answers here
	 $acc = $_;
	 ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time);                      # to let output have a human-readable timestamp
	 $mon++;       #months count from 0
         $year -= 100; #years count from 1900	 
	 if ($debug) {print ("$acc:  $mon\/$mday\/$year $hour:$min:$sec\t\tstarted at $start_time\n");}	 
	 $target_seq = '';
	 $genome_ftp_directory = $paths{$_}[0];                                                        # retrieve the ftp path keyed to the accession number
	 chomp($genome_ftp_directory);
	 $genome_ftp_directory =~ /all\/(.*)/;
	 $fasta_path = ($genome_ftp_directory."/$1"."_protein.faa.gz");                                # build the file path for pulling the fasta from the
	                                                                                               # ftp directory
	 $gff_path = ($genome_ftp_directory."/$1"."_genomic.gff.gz");                                  # build the file path for pulling the gff (coordinates)
	                                                                                               # from the ftp directory
	 	 	  	                     
	   if (-e "$location/ftp_downloads/$_.fasta.gz") {`rm $location/ftp_downloads/$_.fasta.gz`;}                   # zipped files are only left behind if cURL failed a download
	   if (-e "$location/ftp_downloads/$_.gff.gz") {`rm $location/ftp_downloads/$_.gff.gz`;}
	   if (! -e "$location/ftp_downloads/$_.fasta")                                                                # don't download it if we have it
	   {	      
		 print "curl $fasta_path > $location/ftp_downloads/$_.fasta.gz\n"; 
		 `curl $fasta_path > $location/ftp_downloads/$_.fasta.gz`;                                             # individually named so that user can keep files if desired
		  if (! -e "$location/ftp_downloads/$_.fasta.gz") {report_problem($_); next;} 
		           	 `gunzip  $location/ftp_downloads/$_.fasta.gz`;                                                        # files are downloaded zipped
		  $downloaded = 1;
	   }
 	   if (! -e "$location/ftp_downloads/$_.gff")                                                                  # don't download it if we have it
	   {	      
		 print "curl $gff_path > $location/ftp_downloads/$_.gff.gz\n"; 
	         `curl $gff_path > $location/ftp_downloads/$_.gff.gz`;                                                 # individually named so that user can keep files if desired
		 if (! -e "$location/ftp_downloads/$_.gff.gz") {report_problem($_); next;} 
		          	 `gunzip $location/ftp_downloads/$_.gff.gz`;                                                           # files are downloaded zipped
		 $downloaded = 1;
	   }
	 $fasta_path = "$location/ftp_downloads/$_.fasta";                                                             # setting these variables so that we don't have to keep passing paths (that version of the code was a nightmare)
	 $gff_path = "$location/ftp_downloads/$_.gff";
	 %targets = ();                                                                                                # zero the hash since this is a new genome
	 find_target($target_hmm, $target_score, "$location/ftp_downloads/$_.fasta", "$location/ftp_downloads/$_.gff");
	 foreach $identifier (sort keys(%targets))                                                                     # a hack to stop perl from having something keyed to '' which it likes to do
	 { 
	     length($identifier) or next;
	     sort_targets($identifier);
	 }
	 if (!$is_sortable)                                                                                             # if a genome contained no targets that are YES or NO, no need to continue. this is also a hack that lets case==2
	 {                                                                                                              # indicate either NEITHER or FAR - if it is neither but the genome contains at least one target, then it is a YES
		print ("$acc lacks targets that meet entire YES or NO criteria - not included in training set\n");      # that failed a distance rule - thus - FAR
		next;
	 }	
	 foreach $identifier (sort keys(%targets))                                                                      
	 {                                                                                                               # again, a hack to skip the null key - might not be necessary but doesn't hurt anything
     	     length($identifier) or next;
	     $case = $targets{$identifier}[4];                                                                          # the verdict on a target is stored in the HoA that is %targets
	     if ($debug)
	     {
		 open(DEBUG, ">>_JUNK/keys.log");                                                                       # one of several strange log files produced if you run -debug. can be helpful if something goes wrong
		 print DEBUG ("$acc | $identifier\n");
		 close DEBUG;
             }		 
	     if ($case == 0)                                                                           # CASE: NO
	     {
		 $targets{$identifier}[0] =~ />/;
		 $no_seq = ">NO_ATTRIBUTE|$acc|$identifier|$'";                                        # $' means everything after the matched region (one of perl's regular expressions related special variables)
		 open(NO, ">>_NO/no_attribute.seqs");
		 print NO ("$no_seq\n\n");
		 close NO;
		 open(LIST, ">>_TRAINING/training.list");
		 print LIST ("NO_ATTRIBUTE: $acc\t$identifier\n");
		 close LIST;
		 open(TRAINING, ">>_TRAINING/training.seqs");
		 print TRAINING ("$no_seq\n\n");
		 close TRAINING;
		 if ($debug) {print "$acc was NO_ATTRIBUTE\n";}
	     }
	     if ($case == 1)                                                                           # CASE: YES
	     {
		 $targets{$identifier}[0] =~ />/;
		 $yes_seq = ">YES|$acc|$identifier|$'";                                                # $' means everything after the matched region 
		 open(YES, ">>_YES/yes.seqs");
		 print YES ("$yes_seq\n\n");
		 close YES;
		 open(LIST, ">>_TRAINING/training.list");
		 print LIST ("YES: $acc\t$identifier\n");
		 close LIST;
		 open(TRAINING, ">>_TRAINING/training.seqs");
		 print TRAINING ("$yes_seq\n\n");
		 close TRAINING;
		 if ($debug) {print "$acc was YES\n";}
	     }
	     if ($case == 2)                                                                         # CASE: NEITHER
	     {                                                                                       # why is a code that means neither leading to a FAR designation? because the only way a target flagged as neither is in a genome that
		 $targets{$identifier}[0] =~ />/;                                                    # reacjes this code is if it was declared sortable and the only way that happens is if all non-distance rules pass (thus, far)
		 $far_seq = ">FAR|$acc|$identifier|$'";                                              # $' means everything after the matched region
		 open(NO, ">>_NO/far.list");
		 print NO ("FAR: $acc\t$identifier\n");
		 close NO;
		 open(FAR, ">>_NO/far.seqs");
		 print FAR ("$far_seq\n\n");
		 close FAR;
		 open(TRAINING, ">>_TRAINING/training.seqs");
		 print TRAINING ("$far_seq\n\n");
		 close TRAINING;
		 if ($debug) {print "$acc was FAR\n";}
	     }
	     if ($case == 3)                                                                         # CASE: BOTH
	     {
		 open(JUNK, ">>_JUNK/junk.list");
		 print JUNK ("BOTH: $acc\n");
		 close JUNK;
		 if ($debug) {print "$acc was BOTH\n";}
	     }
	 }
       
	 cleanup($_) unless ($keep || !$downloaded);                                                 # Removes the GFF files and fasta files, unless told not to
     }	                                                                                             # Even if we aren't keeping stuff, we won't delete anything we didn't just download (user can't accidentally delete fasta files by mistake)
}                                                                                                    # (well, they can, but not using my software)
###################################################################################################################################################################
sub read_list
{
	my @acc_from_list = ();
	my $list = shift;
	chomp($list);
	open(IN,$list) || die "could not open list input";
	while ($line = <IN>)                                  # returns true unless it hits EOF - will read entire file line by line
        {
		if ($line =~ /^(.*)#/)                        # if line doesn't end in #, does nothing. if it ends with # it will only keep the part before it
		{                                             # this is a convention I add - you can stick comments in your fasta files and my software won't mind
			$line = $1;
			chomp($line);                         # cleanup $line
		}

		if (($line =~ /([A-Z]+_?[0-9]+\.?[0-9]*)/))   # skip lines that aren't valid accessions (NOTE - if users start using fasta files with different names this line will have to go, or only apply with the -acc or -complete flags)
		{
		    push (@acc_from_list, $1);                # doing it this way allows the -list and -acc flags to be combined in a single run (even though you shouldn't do that)
		}
		else
		{
		    next;
		}    
	}   # end while <IN>
	return (@acc_from_list);                              # returns array containing all the accession numbers from the list
}
###################################################################################################################################################################
sub call_hmm                                                                       
												        
{
    my $hmm = shift;
    my $fasta = shift;
    chomp ($hmm);        
    `/usr/local/packages/hmmer-3.1b2/bin/hmmsearch --tblout _TMP/hmm_results.tmp $hmm $fasta`; 									# change to your hmmer3 hmmsearch location
    if ($debug) {print  " /usr/local/packages/hmmer-3.1b2/bin/hmmsearch --tblout _TMP/hmm_results.tmp $hmm $fasta \n";}   		# change to your hmmer3 hmmsearch location

}

###################################################################################################################################################################
sub parse_config                      # subroutine for reading the config file. If users want to add additional functionality, be sure to consider how those rules will be written
{                                     # and try not to change existing syntax, even for private purposes - we want config files to be unambiguous descriptions of the sorting that was done

    open(CONFIG,"$config_file") || die("Couldn't open config file");
    my $is_yes_rule = 0;
    my $is_no_rule = 0;
    my $rule_part_ctr = 0;
    my $target_hmm = 0;
    my $target_score = 0;
    while ($line = <CONFIG>)
    {
	if (!($line =~ /^#/))                                       ## This allows the user to include comments in the config file
	{
	    $is_yes_rule = 0;
	    $is_no_rule = 0;
	    $rule_part_ctr = 0;
	    
	    if ($line =~ /\((.*)\)/)                               ## If line has (), remove them and split on commas
	    {
		if ($debug) {print "\n\n$1\n\n";}
		@split_line = split(/,/,$1);
		if ($line =~ /YES/)
		{
		    for ($k=0; $k<@split_line; $k++)              ## prefix each part of compound YES rule with YES:
		    {
			$split_line[$k] = "YES:" . $split_line[$k];
		    }
		}
		if ($line =~ /NO/)
		{
		    for ($w=0; $w<@split_line; $w++)             ## prefix each part of compound NO rule with NO:
		    {
			$split_line[$w] = "NO:" . $split_line[$w];
		    }
		} 
	    }

	    else 
	    {
		$split_line[0] = $line;
	    }
	   
	    foreach (@split_line)
	    {
	
		if($_ =~ /^(TARGET:\s*)(.*)/)
		{       
		    $target_hmm = $2;                                          # capture HMM to be used as TARGET rule
		    last;
		}
	    
		if($_ =~ /^(TARGET SCORE:\s*)(\d+(\.\d+)?)/)
		{
		    $target_score = $2;                                       # capture SCORE to be used as TARGET score
		    last;
		}
	
		if($_ =~ /^(YES:\s*)(\[(\d+)\])?(.*)([><])(\d+(\.\d+)?)/)             # regex for parsing user rules (remember, $n refers to everything between the nth left paren and the matching right paren)
		{
		    $is_yes_rule = 1;
                    $debug_rule = ("1: $1\n2: $2\n3: $3\n4: $4\n5: $5\n6: $6");       # Useful if you think your rules are being parsed incorrectly 
		    if ($debug) {print ("$debug_rule\n");}
		    if ($2) {$uses_coords = 1;}
		    $rules[1][$yes_rule_ctr][$rule_part_ctr][0] = $4;                  # The HMM path
		    $rules[1][$yes_rule_ctr][$rule_part_ctr][1] = $6;                  # The HMM required score
		    if ($5 eq ">")
		    {
			$rules[1][$yes_rule_ctr][$rule_part_ctr][2] = 1;               # 1 asserts > score, 0 asserts < score
		    }                                                                  # while we don't forbid combining a "<" with a distance rule, what does that even mean? 
		    else                                                               # that behavior is completely untested and will never be supported (unless you do a really good job convincing me that you need it)
		    {
			$rules[1][$yes_rule_ctr][$rule_part_ctr][2] = 0;               # a < rule will fail if any hit above the given score is detected. used in most no rules, but you sometimes might want it in a yes rule
		    }
		    $rules[1][$yes_rule_ctr][$rule_part_ctr][3] = $3;                  # captures the distance requirement if one was given. part of the reason that the the various fields are saved in @rules in a different order than
		    $rule_part_ctr++;                                                  # in the config file is so that the only-sometimes-present distance component could be at the end
		                                                                       # while it is incremented after every rule-part, a non-zero value is only every used before being reset if a rule has a logical OR in it 
		                                                                       # which is represented by multiple rules on one line, grouped in () and separated by commas. the YES: or NO: goes before the ()
		}
		elsif($_ =~ /^(NO:\s*)(\[(\d+)\])?(.*)([><])(\d+(\.\d+)?)/)            # Not typing all that again - see the parsing for YES
		{                                                                      # IMPORTANT - yes and no rules both go into @rules, but yes goes into $rules[1] and no goes into $rules[0] 
		    $is_no_rule = 1;                                                   # if you have added new functionality yourself and are getting baffling results, that is probably the problem, trust me
                    $debug_rule = ("1: $1\n2: $2\n3: $3\n4: $4\n5: $5\n6: $6");
		    if ($debug) {print ("$debug_rule\n");}
		    if ($2) {$uses_coords = 1;}
		    $rules[0][$no_rule_ctr][$rule_part_ctr][0] = $4;                                        
		    $rules[0][$no_rule_ctr][$rule_part_ctr][1] = $6;
                    if ($5 eq ">")
		    {
			$rules[0][$no_rule_ctr][$rule_part_ctr][2] = 1;
		    }
		    else
		    {
			$rules[0][$no_rule_ctr][$rule_part_ctr][2] = 0;                     
		    }
		    $rules[0][$no_rule_ctr][$rule_part_ctr][3] = $3;
		    $rule_part_ctr++;

		}
		else
		{
		    die("Unable to parse config file line: $line");                      # We chose to make this a fatal error because why waste the user's time by running a config file that isn't what they meant
		}    
	
	    }
	    $yes_rule_ctr += $is_yes_rule;                                               # These are used both for rule storage and making sure both YES and NO are represented. 
	    $no_rule_ctr += $is_no_rule;                                                 # It is nice when things work out nicely like that
	}
    }
    if (!$yes_rule_ctr) {die "No YES rules found in config file\n";}                      # We don't let the user run a config file that doesn't have both YES and NO rules
    if (!$no_rule_ctr) {die "No NO rules found in config file\n";}                        # if you honestly don't care, have YES be what you want and have NO use some random HMM and requrire a score of 1000 or higher - it will just not find
    return($target_hmm, $target_score);                                                                                                                                                                                               # anything
}
# rules are stored like this:
# [0 or 1 corresponding to No or Yes][Rule number][Rule part number (only non-zero when a rule contains an OR)][various parts of the rule go here]
# for the last index, 0 contains the HMM path, 1 contains the score, 2 contains the sign of the rule, 3 contains the distance (if specified)

######################################################################################################################################################################
sub find_target                                                      # a sub that locates all sufficiently high-scoring target sequences from a genome and sticks them in %targets along with coordinate information and sequence data
{                                                                    
    my $target_hmm = shift;
    my $target_score = shift;
    my $fasta_path = shift;
    my $tmp_start = 0;
    my $tmp_stop = 0;
    call_hmm($target_hmm, $fasta_path);                               # runs the HMM used for target detection - results are temporarily stored in ./_TMP
    target_data_from_hmm($target_score);                              # looks at the results file in ./_TMP and populates %targets with the sequences that score above $target_score
    foreach $identifier (sort keys(%targets))                         # now that we have 0 or more targets, time to give them positional information from the gff file
    {
	if ($debug) {print "$identifier\n";}
	$tmp_seq = '';                                                # $tmp_seq, $line, and $seq_line work together for a easily reproducable way to pull specific sequences from a multi-fasta file
								      # you probably don't want to mess with that unless you really have to
	$tmp_molecule = '';
	$tmp_midpoint = '';
	($tmp_molecule, $tmp_start, $tmp_stop) = gff_lookup($identifier);   # usually, $identifier is a WP - pull positional information from the gff file
	$tmp_midpoint = (($tmp_start + $tmp_stop)/2);
	$found = 0;
	open(FIND_SEQ,"$fasta_path");
	while ($line = <FIND_SEQ>)
	{
	    if ($line =~ />$identifier/ && $identifier)                # without the "&& $identifier", if a blank identifier sneaks in (and with Perl it very well might), it will match the next sequence in the file
	    {
		$found = 1;
		$tmp_seq = $line;
		while ($seq_line = <FIND_SEQ>)
		{    
		    if (!($seq_line =~ />/))
		    {
			$tmp_seq = ($tmp_seq . $seq_line)              # (seriously even I find this looping a bit confusing and I wrote it. but I've used it in several programs and it is well-tested)
		    }
		    else {last;}
		}	
	    }
	    if ($found) {$targets{$identifier}[0] = $tmp_seq; last;}
	  
        }    
	$targets{$identifier}[1] = $tmp_molecule;                    # capture molecule information - if two sequences are on different molecules, they aren't close, I don't care what their coords say
	$targets{$identifier}[2] = $tmp_midpoint;                    # this line isn't actually helpful

    }	
}
#########################################################################################################################################################################
sub target_data_from_hmm
{

    my $threshold = shift;                                                                     ## we only care about targets >= this threshold score
    my @split_line = ();                                                                       ## used with split to make parsing easy
    my $hmm_score = 0;
    my $match_identifier = '';
    %results = ();
    open(HMMR,"_TMP/hmm_results.tmp") || die "failed to open hmm results";            		## open temporary file containing hmm_results
    while ($line = <HMMR>)                                                            		## loop through file one line at a time
    {
	if (!($line =~ /^#/))                                                         		## line doesn't start with "#", so is a result line
	{
	    @split_line = split(/\s+/, $line);                                       		## splits on whitespace, now each value is a different array element
#	    if (($split_line[8] > $hmm_score) && ($split_line[5] > $threshold))                 ## commented out not deleted in case users want to re-enable consideration of domain score - you probably don't want this though
												## use the domain extractor tool instead
	    if ($split_line[5] > $threshold)
	    {
#		$hmm_score = $split_line[8];                                                    ## if you are looking for the best HMM hit, you want something with protein score over threshold, but highest domain score
		$hmm_score = $split_line[5];						        ## when you want all targets, not the best, this logic is actually harmful
		$split_line[0] =~ />?(\w+\.?\w*)/ ;
		$match_identifier = $1;
		if ($debug)
		{
		     open(DEBUG, ">>_JUNK/hmm.log");
		     print DEBUG "$split_line[0] | $match_identifier\n$split_line[1]\n$split_line[7]\n\n";
		     close DEBUG;
		     open(DEBUG2, ">>_JUNK/target.log");
		     print DEBUG2 "$match_identifier: $hmm_score > $threshold\n $split_line[21] \n $split_line[0]\n\n";
		     close DEBUG2;
		 }
		$targets{$match_identifier}[3] = $hmm_score;                                    ## I think this is only used for debugging purposes
	    }	
	}    
    }
}    
########################################################################################################################################################################
sub sort_targets                                                                   # very straightforward sub - launches the yes and no tests, and has handling for all 4 possible outcomes
{
    		my $identifier = shift;
		my $is_yes = ''; my $is_no = '';		                   # we don't use a zero because then we couldn't notice the (hopefully never observed) case where we fail to make a decision 
										   # we don't actually check for that case, but that may be soon implemented, and even now would show up in log files       
		($is_yes, $is_no) = do_tests($identifier);                         # do_tests is a sub that... does tests. Two return values - for both yes and no, a 0 or 1 corresponding to success or failure on the tests
		if ($debug) {print  ("$acc: is_yes?: $is_yes\t is_no?: $is_no\n");}
 		if (!$is_yes && !$is_no)
    		{
		$targets{$identifier}[4] = 2;                                      # 2: Genome met neither yes nor no condition
    		} elsif ($is_yes && $is_no)
    		{
		$targets{$identifier}[4] = 3;                                      # 3: Genome met both yes and no conditions
    		} elsif ($is_no)
    		{
		$targets{$identifier}[4] = 0;                                      # 0: Genome met only no condition
		$is_sortable = 1;                                                  # This flag gets set when a target is actually a YES or NO - without it, no targets from the genome are used
    		} else
    		{
		$targets{$identifier}[4] = 1;                                      # 1: Genome met only yes condition
		$is_sortable = 1;                                                  # This flag gets set when a target is actually a YES or NO - without it, no targets from the genome are used
    		}							           # If this flag is set, then a target marked as NEITHER is actually called FAR - the only way that a sequence can fail both criteria
	                                                                           # when at the same time the genome contains sortable sequences is if the sequence in question failed a distance rule
}                                                                                  # a bit hacky, but it simplifies the data structure by a fair bit
#########################################################################################################################################################################
sub do_tests                                                                       # Simple subroutine that actually launches the yes and no tests, then returns both results
{
	my $current_target = shift;
	my $is_yes = check_yes_rules($current_target);
	my $is_no = check_no_rules($current_target);
	return ($is_yes, $is_no);
}
#########################################################################################################################################################################
sub check_yes_rules                                                # My apologies - this logic is kinda hard to follow, but it is also computationally efficient (and fun to write)
{
    my $i = 0;      #rule number                                                  
    my $j = 0;      #rule part
    my $current = shift;
    while (1)                                                       # the while loop never ends, but also it never traditionally loops. Instead of ending we just return, and we loop by using next.
    {
	if (!defined($rules[1][$i][$j][0]))                         # a check to see if we've run out of rules in the given context - in the same rule (OR) or overall (AND), depending on where in the logic we hit the next;
	{
	    if ($debug) {print "about to return 1, yes, $i $j\n";}  # if you have exhausted all the rules and your most recent test was a success, then you passed
	    return 1;
        }
	if (evaluate_yes($i, $j, $current))
	{
	    if ($debug)
	    {
		open(DEBUG, ">>_JUNK/tests.debug");
		print DEBUG "ySUCCESS: $rules[1][$i][$j][0] , $rules[1][$i][$j][1], $rules[1][$i][$j][2] |i|$i|j|$j\n\n";   #y means it was a yes rule and SUCCESS means the rule was satisfied
		close DEBUG;
	    }	
	    $i++;                                                           # if you pass a test, go to the next one
	    $j=0;                                                           # this short-circuts a rule that contains a logical OR  (multiple rules on a line)
	    next;
	}
	else
	{
	    if ($debug)
	    {
		open(DEBUG, ">>_JUNK/tests.debug");
		print DEBUG "yFAIL: $rules[1][$i][$j][0] , $rules[1][$i][$j][1], $rules[1][$i][$j][2] |i|$i|j|$j\n\n";        #y means it was a yes rule and FAIL means the rule wasn't satisfied
		close DEBUG;
	    }	
	    $j++;                                                         # if you fail a test, check to see if the rule has multiple parts
	}    
	if (!defined($rules[1][$i][$j][0]))                               # a check to see if we've run out of rules in the given context - in the same rule (OR) or overall (AND), depending on where in the logic we hit the next;
	{
	    if ($debug) {print "about to return 0, yes, $i $j\n";}        # if you have exhaused all the rules and your most recent test was a failure, then you failed
	    return 0;                                                     # this short-circuits a config file that contains a logical AND (multiple lines of rules)
	}
	else
	{
	    next;
	}    
	    
    }	

}

#########################################################################################################################################################################
sub check_no_rules                                                 # like check_yes_rules() but with 0s not 1s
{
    my $i = 0;
    my $j = 0;
    my $current = shift;
    while (1)
    {
	if (!defined($rules[0][$i][$j][0]))
	{
	    if ($debug) {print "about to return 1, no, $i $j\n";}   # if you have exhausted all the rules and your most recent test was a success, then you passed
	    return 1;
        }
	if (evaluate_no($i, $j, $current))
	{
	    if ($debug)
	    {
		open(DEBUG, ">>_JUNK/tests.debug");
		print DEBUG "nSUCCESS: $rules[0][$i][$j][0] , $rules[0][$i][$j][1], $rules[0][$i][$j][2] |i|$i|j|$j\n\n";
		close DEBUG;
	    }	
	    $i++;                                                 # if you pass a test, go to the next one
	    $j=0;                                                 # this short-circuts a rule that contains a logical OR
	    next;
	}
	else
	{
	    if ($debug)
	    {
		open(DEBUG, ">>_JUNK/tests.debug");
		print DEBUG "nFAIL: $rules[0][$i][$j][0] , $rules[0][$i][$j][1], $rules[0][$i][$j][2] |i|$i|j|$j\n\n";
		close DEBUG;
	    }	
	    $j++;                                                  # if you fail a test, check to see if the rule has multiple parts
	}    
	if (!defined($rules[0][$i][$j][0]))
	{
	    if ($debug) {print "about to return 0, no, $i $j\n";}  # if you have exhaused all the rules and your most recent test was a failure, then you failed
	    return 0;
	}
	else
	{
	    next;
	}    
	    
    }	

}
#####################################################################################################################################################################
sub evaluate_yes                                                                 # sub for calculating a rule - distance rules need to be calculated once per target, but other rules can be done once per genomes
{
    my $x = shift;
    my $y = shift;
    my $target = shift;
    if ($rules[1][$x][$y][3])                                                    # check to see if a rule is a distance rule - if it is, compute it, don't look for a stored value
    {
	if ($debug) {print "I am a yes distance rule [x = $x| y = $y]\n";}
	return evaluate_attribute_distance($rules[1][$x][$y][0], $rules[1][$x][$y][1], $rules[1][$x][$y][2], $rules[1][$x][$y][3], $target);	
    }
    else
    {
	if(!defined($answers[$x][$y][1]))                                        # for rules that aren't distance rules, check if we've already stored a value from earlier in the processing of the genome
	{
		$answers[$x][$y][1] = evaluate_attribute_distance($rules[1][$x][$y][0], $rules[1][$x][$y][1], $rules[1][$x][$y][2], $rules[1][$x][$y][3], $target);  # store the value for later lookup
		if ($debug) {print "I am a yes non-distance rule of unknown value [x = $x| y = $y]\n";}
	}
	else                                                                                                                                                         # if we get here, we already know the value
	{
		if ($debug) {print "I am a yes non-distance rule of known value $answers[$x][$y][1]      [x = $x| y = $y] \n";}
	}
	return $answers[$x][$y][1];                                                                                                                                   # whether we looked it up or just computed it, return the value
    }
}
###################################################################################################################################################################
sub evaluate_no                                                                 # same as evaluate_yes, but with 0s not 1s
{
    my $x = shift;
    my $y = shift;
    my $target = shift;
    if ($rules[0][$x][$y][3])
    {
	if ($debug) {print "I am a no distance rule [x = $x| y = $y]\n";}
	return evaluate_attribute_distance($rules[0][$x][$y][0], $rules[0][$x][$y][1], $rules[0][$x][$y][2], $rules[0][$x][$y][3], $target);
    }
    else
    {
	if(!defined($answers[$x][$y][0]))
	{
		$answers[$x][$y][0] = evaluate_attribute_distance($rules[0][$x][$y][0], $rules[0][$x][$y][1], $rules[0][$x][$y][2], $rules[0][$x][$y][3], $target);
		if ($debug) {print "I am a no non-distance rule of unknown value [x = $x| y = $y]\n";}
	}
	else
	{
		if ($debug) {print "I am a no non-distance rule of known value $answers[$x][$y][0]          [x = $x| y = $y]\n";}
	}
	return $answers[$x][$y][0];
    }
}

######################################################################################################################################################################
sub evaluate_attribute_distance                       # determines if any sufficiently high-scoring instance of the given attribute is within range of target sequence
{
	%attribute = '';	
	my $hmm_path = shift;
	my $hmm_score = shift;
	my $hmm_sign = shift;                         # 1 means the rule is checking >, 0 means the rule is checking < 
	my $hmm_distance = shift;
	my $target = shift;
	my $exists = 0;                               # 1 means at least 1 attribute was detected above the score cutoff
	call_hmm($hmm_path, $fasta_path);             # call_hmm has no return value - instead writes to a file in _TMP, which is searched by the _from_hmm subs
	$exists = attribute_data_from_hmm($hmm_score);# return value corresponds to whether or not a sufficiently high scoring hmm-hit was found
	foreach $identifier (sort keys(%attribute))
    	{
		if ($debug) {print "$identifier\n";}
		$tmp_seq = '';
		$tmp_molecule = '';
		($tmp_molecule, $tmp_start, $tmp_stop) = gff_lookup($identifier);                            # collects position information from the gff file
		$found = 0;
		open(FIND_SEQ,"$fasta_path");
		while ($line = <FIND_SEQ>)                                                                   # here is that code for sequence collection again
		{
	    		if ($line =~ />$identifier/ && $identifier)                                          # && $identifier stops us from always including the first sequence
	    		{
				$found = 1;
				$tmp_seq = $line;
				while ($seq_line = <FIND_SEQ>)
					{    
		    				if (!($seq_line =~ />/))
		    				{
							$tmp_seq = ($tmp_seq . $seq_line)
		    				}
		    				else {last;}
					}	
	    		}
	    		if ($found) {$attribute{$identifier}[0] = $tmp_seq; last;}                            # we store the sequence data of the attribute - we might end up wanting that
	                                                                                                      
        	}    

		$attribute{$identifier}[1] = $tmp_molecule;
		$attribute{$identifier}[2] = $tmp_start;
		$attribute{$identifier}[3] = $tmp_stop;

    	}	
	if ($exists && $hmm_sign)                                                                              # only in the case where you want an attribute to exists and it does do you actually need to check distance
	{
	    if ($hmm_distance) {$result = check_distance($target, $hmm_distance);}                             # don't check distance if your rule isn't a distance rule
	    else {$result = 1;}                                                                                # if you weren't a distance rule, if it exists, it exists and is close
	    return $result;                                                                                    # non-distance rules return 1, distance rules return 1 or 0 depending on... distance
	}
	elsif (!$exists && $hmm_sign)                                                                          # if it doesn't exists but was supposed to, the rule failed
	{
	    return 0;
	}
	elsif (!$exists && !$hmm_sign)                                                                         # if it doesn't exist and wasn't supposed to, the rule succeeded
	{	
		return 1;
	}
	else                                                                                                   # for those keeping score, this last case is ($exists && !$hmm_sign)
	{
		return 0;
	}
	
}
######################################################################################################################################################################
sub attribute_data_from_hmm                                                                     ## sub that reads the hmm results file in _TMP, populates %attribute, and returns 1 if a high scoring match is found
{
    my $threshold = shift;
    my @split_line = ();
    my $hmm_score = 0;
    my $match_identifier = '';
    my $found = 0;
    %results = '';
    open(HMMR,"_TMP/hmm_results.tmp") || die "failed to open hmm results";            		## open temporary file containing hmm_results
    while ($line = <HMMR>)                                                            		## loop through file one line at a time
    {
	if (!($line =~ /^#/))                                                         		## line doesn't start with "#", so is a result line
	{
	    @split_line = split(/\s+/, $line);                                       		## splits on whitespace, now each value is a different array element
#	    if (($split_line[8] > $hmm_score) && ($split_line[5] > $threshold))                 ## commented out not deleted in case users want to re-enable consideration of domain score - you probably don't want this though
												## use the domain extractor tool instead
	    if ($split_line[5] > $threshold)
	    {
#		$hmm_score = $split_line[8];
		$hmm_score = $split_line[5];                                                   ## since we aren't doing domain calculations, lets just keep the overall score
		$split_line[0] =~ />?(\w+\.?\w*)/ ;                                            
		$match_identifier = $1;
		$found = 1;
		if ($debug)
		{
		     open(DEBUG, ">>_JUNK/hmm.log");
		     print DEBUG "$split_line[0] | $match_identifier\n$split_line[1]\n$split_line[7]\n\n";
		     close DEBUG;
		     open(DEBUG2, ">>_JUNK/target.log");
		     print DEBUG2 "$match_identifier: $hmm_score > $threshold\n $split_line[21] \n $split_line[0]\n\n";
		     close DEBUG2;
		 }
		$attribute{$match_identifier}[4] = $hmm_score;
	    }	
	}    
    }
    return $found;
}
########################################################################################################################################################################
sub check_distance                                                               # checks if the given target is within the given distance of ANY instance of the current attribute
{                                                                                # %attribute gets cleared after each rule, so we don't need to specify attribute
	my $gene_identifier = shift;
	my $required_distance = shift;
	my ($molecule, $start, $stop) = gff_lookup($gene_identifier);            # yes, we looked this up earlier, but shhhhhhhh
	my $target_midpoint = (($start + $stop)/2);                              # simple midpoint calculation
	my $close = 0;
	foreach $identifier (sort keys(%attribute))
	{
		$attribute_midpoint = (($attribute{$identifier}[2] + $attribute{$identifier}[3])/2);
		$sum_half_lengths = ((abs($start - $stop) + abs($attribute{$identifier}[2] - $attribute{$identifier}[3]))/2);   # we compute sum half lengths so that we can measure from endpoint to endpoint without determining orientation
		if ($debug)
		{
		print("start: $start\tstop: $stop\tstart2: $attribute{$identifier}[2]\tstop2: $attribute{$identifier}[3]\tsum_half_lengths: $sum_half_lengths\n"); 
		}
	        if (!($molecule eq $attribute{$identifier}[1]))                     # check if on the same molecule
		{
			next;                                                       # if you aren't on the same molecule, then I don't care what your coordinates say
		}
	
		if (abs($target_midpoint - $attribute_midpoint) > ($sum_half_lengths + $required_distance))
		{
			next;
		}
		else
		{
			$close = 1;
			last;
		}
	}
	return $close;
}

#######################################################################################################################################################################
sub gff_lookup
{
    my $lookup_gene = shift;
    my @split_line = ();
    open(GFF, $gff_path) || die"failed to open gff file at $gff_path\n";
    while ($line = <GFF>)
    {
	if (!($line =~ /^#/))                                                         		## line doesn't start with "#", so is a result line
	{
	    @split_line = split(/\t/,$line);
	    if ($split_line[8] =~ /$lookup_gene/)
	    {
		if ($debug)
		{
		    open(DEBUG, ">>_JUNK/debug.log");
		    print DEBUG "$split_line[0]\n$split_line[1]\n$split_line[2]\n$split_line[3]\n$split_line[4]\n";
		    close DEBUG;
		}    
		return ($split_line[0], $split_line[3], $split_line[4]);
            }	
	}    

    }	

}    
###################################################################################################################################################################
sub cleanup                                             
{
    $cleanup_target = shift;
    {
	`rm ./ftp_downloads/$cleanup_target*`                                 
    }

}    
##################################################################################################################################################################
sub report_problem
{
	$problem = shift;
	open(PROBLEM, ">>FAILED_TO_DOWNLOAD");
	print PROBELM "$problem\n";
	close PROBLEM;
	$problematic = 1;
}

