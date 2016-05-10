use FileHandle;
use Getopt::Long;
#use warnings;
my $junk = "extraction.tmp";
my $usage = "
This program requires the following flags:
-seqs <file path to the FASTA file>
-hmm <file path to the HMM>
-domain <domain score cutoff>
-protein <protein score cutoff>

Optional flag:
-tmp <override the default location for the HMM results>
";
my %domains = ();
GetOptions ( 'seqs=s' => \$seqs,
	     'hmm=s' => \$hmm,
             'domain=s' => \$domain_cutoff,
	     'protein=s' => \$protein_cutoff,
	     'tmp=s' => \$junk	);

if (!$hmm || !$domain_cutoff || !$protein_cutoff) {die "$usage\n";}
if (!$seqs) {die "$usage\n";}
`hmmsearch -T $protein_cutoff --domtblout $junk  $hmm  $seqs`;

read_domtblout($junk);
search_fasta($seqs);

#########################################################################################################################################################################
sub read_domtblout
{
	my $source = shift;
	open (TABLE, "$source") || die "cannot read domtblout file\n";
	while ($line = <TABLE>)
	{
		if (!($line =~ /^#/))                                                         	       # ignore the header line
		{

		        @split_line = split(/\s+/, $line);                                              # splits on tabs - some entries contain whitespace
			if ($split_line[13] >= $domain_cutoff)
			{
				push(@{$domains{$split_line[0]}[0]}, $split_line[19]);                    # adds "env from" coordinate to array
				push(@{$domains{$split_line[0]}[1]}, $split_line[20]);                    # adds "env to" coordinate to array
				# %domains is a hash, but $domains{identifier}[0] and $domains{$identifier}[1] are both arrays
				# this way, all domains from one sequence are stored with the same hash key, but can easily be processed iteratively
			}

		}
    
	}
	close TABLE;
}
##########################################################################################################################################################################
sub search_fasta
{
	my $fasta = shift;
	open (FASTA, "$fasta") || die "cannot read sequence file\n";
	my $tmp_seq = '';	
	while ($line = <FASTA>)
	{
		if ($line =~ /^>(.*?)\s/)
		{
			if ($domains{$1})
			{
				$identifier = $1;
				$tmp_seq = '';
				$header = $line;
				while ($seq_line = <FASTA>)
				{
					if (!($seq_line =~ />/))
		    			{
						chomp($seq_line);							
						$tmp_seq = ($tmp_seq . $seq_line)
		    			}
		    			else {last;}

				}

				for ($i = 0; $i <= $#{$domains{$identifier}[0]}; $i++)
				{
					$from = ($domains{$identifier}[0][$i] - 1);
					$to = $domains{$identifier}[1][$i];
					$length = ($to - $from + 1);
#					$tmp_seq =~ /.{$from}(.{$length})/;
					$sub_seq = substr $tmp_seq, $from, $length;
					if ($sub_seq && $length) {chomp($header); print("$header\($from\.$to\)\n$sub_seq\n\n");}
				}
			}
			else
			{
				next;
			} 
		}
		
		
	}
} 











	
