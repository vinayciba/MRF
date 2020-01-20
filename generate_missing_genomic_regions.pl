
# This script takes three column mummer output and calculates the gap postions between sucessive matchin regions


my $file = $ARGV[0];	# READ MUMMER OUTPUT FILE
my $genome_size = $ARGV[1]; # TAKE PRE CALCULATED SIZE OF THE REFERENCE FASTA
my $output_prefix = $ARGV[2]; # ASSIGNING OUTPUT PREFIX

chomp($genome_size);

# read input file from the prompt
open (IN, $file) or die("We have got a situation!! Problem with mummer output!!");
#open a filehandle to mummer output file

# go through the file line by line    

while (my $line = <IN>) 
	{
# split the current line on tabs
	my @words = split(/\t/, $line);
	push(@Ref_genome_start, $words[0]);
	push(@Query_genome_start, $words[1]);
	push(@Hit_region_length, $words[2]);
	}

my @Ref_genome_end;
my @Query_genome_end;
my @Non_hit_region;
my @Q_genome_miss_start;
my @Q_genome_miss_end;
my $curr_line;
my $fin_line;
my $tot_mr_len =0;
my $trim_cds_percent;
my $Ref_init_end;
my $start_line;
my $genome_gap;
my $end_line;
my $genome_fin_start;
my $count = scalar(@Ref_genome_start) - 1;
my $count_act = scalar(@Ref_genome_start);

# get the count of items in the array list of Ref_genome_start
# and execute the loop for each item

for(my $i=0,$j=1;$i<$count_act;$i++,$j++)
	{
	$Ref_genome_end[$i] = ($Ref_genome_start[$i] + $Hit_region_length[$i]) - 1;
	$Query_genome_end[$i] = ($Query_genome_start[$i] + $Hit_region_length[$i]) - 1;
	$Non_hit_region[$j] = ($Ref_genome_start[$j] - $Ref_genome_end[$i]) - 1;
	$Q_genome_miss_start[$j] = $Ref_genome_end[$i] + 1;
	$Q_genome_miss_end[$j] = ($Q_genome_miss_start[$j] + $Non_hit_region[$j]) - 1;
	$Non_hit_region[0] = 0;
	$Q_genome_miss_start[0] = 0;
	$Q_genome_miss_end[0] = 0;

  	if($Non_hit_region[$i] > 0)
  	 	{
   		$tot_mr_len = $tot_mr_len + $Non_hit_region[$i];
   		}
#  	if($Non_hit_region[$i] <= 2)
#   		{
#		$Q_genome_miss_start[$i] = "";
#		$Q_genome_miss_end[$i] = "";
#		$Non_hit_region[$i] = "";
#		}

	 if($Ref_genome_start[0] ne 1)
                {
                $Ref_init_end = $Ref_genome_start[0] - 1;
                $start_line = "1\t$Ref_genome_start[0]\t0\t0\t$Ref_init_end\t0\t0";
                $Non_hit_region[0] =  $Ref_init_end;
                $Q_genome_miss_start[0] = 1;
                $Q_genome_miss_end[0] =  $Ref_genome_start[0] - 1;
                }

	
#	$genome_gap = ($genome_size - $Ref_genome_end[$count]);
	$curr_line = "$Ref_genome_start[$i]\t$Ref_genome_end[$i]\t$Query_genome_start[$i]\t$Query_genome_end[$i]\t$Non_hit_region[$i]\t$Q_genome_miss_start[$i]\t$Q_genome_miss_end[$i]";
	$fin_line =  $fin_line . "\n" . $curr_line;
#	$genome_fin_start = $Ref_genome_end[$count] + 1;
#	$end_line =  "0\t0\t0\t0\t$genome_gap\t$genome_fin_start\t$genome_size";
	}
$fin_line =~ s/^(.*\n){1}//;	
$genome_gap = ($genome_size - $Ref_genome_end[$count]);
$genome_fin_start = $Ref_genome_end[$count] + 1;
 $Q_genome_miss_start[$count] = $genome_fin_start;
        $Q_genome_miss_end[$count] = $genome_size;


$end_line =  "0\t\t0\t0\t$genome_gap\t$genome_fin_start\t$genome_size";

	if($genome_gap != 0)
                {
                $fin_line = $fin_line . "\n" . $end_line;
                }



# close the filehandle and exit
close IN;


$mr_header = "Ref_start \t Ref_end \t Query_start \t Query_end \t Non_hit_region \t Query_miss_start \t Query_miss_end \n";

my $mr_final = "$mr_header $fin_line\n";

my $mr_pre = $output_prefix . "_mr.txt";


## printing output to a file

open (fh, ">", "$mr_pre");
print fh $mr_final;
close(fh) or "Couldn't close the file";


exit;
