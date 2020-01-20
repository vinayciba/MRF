use List::MoreUtils qw(uniq);

my $file = $ARGV[0];
# This script takes mummer output and a reference GFF file as input and prints the completely and partially missing proteins
# open a filehandle to data.txt
open (IN, $file) or die("We have got a situation!! Problem with mummer output!!");
# go through the file line by line    
my $genome_size = $ARGV[2];
my $output_prefix = $ARGV[3];
chomp($genome_size);
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
my $genome_gap;
my $start_line;
my $end_line;
my $genome_fin_start;
my $count = scalar(@Ref_genome_start) - 1;
my $count_act = scalar(@Ref_genome_start);
my $Ref_init_end;
my $cds_count =0;
my @prot_par_count;

# get the count of items in the array list of Ref_genome_start
# and execute the loop for each item

for(my $i=0,$j=1;$i<scalar(@Ref_genome_start);$i++,$j++)
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
	if($Non_hit_region[$i] <= 2)
 		{
		$Q_genome_miss_start[$i] = "";
		$Q_genome_miss_end[$i] = "";
		$Non_hit_region[$i] = "";
		} 
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
	$genome_gap = ($genome_size - $Ref_genome_end[$count]);
	$genome_fin_start = $Ref_genome_end[$count] + 1;
	$Ref_genome_start[$count] = end;
	$Q_genome_miss_start[$count_act] = $genome_fin_start;
	$Q_genome_miss_end[$count_act] = $genome_size;

	$end_line =  "0\t0\t0\t0\t$genome_gap\t$genome_fin_start\t$genome_size";
	$fin_line =~ s/^(.*\n)1}//;
	if($genome_gap != 0)
		{
		$fin_line = $fin_line . "\n" . $end_line;
		}

# close the filehandle and exit
close IN;
my $gfffile =  $ARGV[1];
open(FH, $gfffile);
@array = <FH>;
close FH;
open(OUT,'>','featuretable.gff3');
print OUT @array[0..$#array];
close OUT;

# open a filehandle to data.txt
open (IN, "featuretable.gff3") or die;
my @product;
# my $mrc = 0;
my $cds = 0;
my $cmr = "";
my $pmr = "";
my $cmr_tmp = 0;
my $pmr_tmp = 0;
my $tot_cmr_len =0;
my $tot_pmr_len =0;
# go through the file line by line    
while (my $line = <IN>)
	{
next if /^#/;
# split the current line on tabs
    	my @position = split(/\t/, $line);
  	push(@category, $position[2]);
	push(@CDS_start, $position[3]);
    	push(@CDS_end, $position[4]);
    	push(@trash, $position[8]);
	}
foreach my $k (@CDS_start)
	{
my $miss = 0;
if ($category[$cds] =~ m/CDS/)
{
	my ($pdt) = $trash[$cds] =~ m/product=([^=]+)\;/;
        my ($protein_id) = $trash[$cds] =~ m/protein_id=([^=]+)[\;|\s]/;
# my ($protein_id) = $trash[$cds] =~ m/protein_id=([^=]+\s+|\;)/;
# my ($protein_id) = $trash[$cds] =~ m/protein_id=([^=]+)\; | protein_id=([^=]+\s+); | protein_id=([^=]+)\s+/;
        my ($note) = $trash[$cds] =~ m/Note=([^=]+\w+)\;/;
# my ($protein_id) = $trash[$cds] =~ m/protein_id=([^=]+\s+)/;
	for(my $a=0;$a<=scalar(@Ref_genome_start);$a++)
		{
		 		if ($Q_genome_miss_start[$a] != "")
				{
				if ($CDS_start[$cds] >= $Q_genome_miss_start[$a] && $CDS_end[$cds] <= $Q_genome_miss_end[$a]) #CONDITION 1: FETCHING COMPLETELY DELETED PROTEINS
					{
					my $cmr_len =0;
					$cmr_len = $CDS_end[$cds] - $CDS_start[$cds] +1;
						if ($cds == 0)
						{
							 $tot_cmr_len = $tot_cmr_len + $cmr_len;
                                                         $cmr_tmp = $CDS_end[$cds];
						}
						elsif ($cds != 0 && $CDS_start[$cds] > $cmr_tmp)
						{
								
							$tot_cmr_len = $tot_cmr_len + $cmr_len;
							$cmr_tmp = $CDS_end[$cds];
						}
						elsif ($CDS_start[$cds] <=  $cmr_tmp && $CDS_end[$cds] >  $cmr_tmp)
						{
							 $tot_cmr_len = $tot_cmr_len + ($CDS_end[$cds] -  $cmr_tmp);
							 $cmr_tmp = $CDS_end[$cds];
						}
					$cmr = $cmr . "\n" . $Q_genome_miss_start[$a] . "\t" . $Q_genome_miss_end[$a] .  "\t" . $CDS_start[$cds] . "\t" . $CDS_end[$cds] . "\t" . $cmr_len . "\t" . $pdt . "\t" . $protein_id . "\t" . $note . "\n";
					$miss = 1;
					}
				if ($CDS_start[$cds] >= $Q_genome_miss_start[$a] && $CDS_start[$cds] <= $Q_genome_miss_end[$a] && $CDS_end[$cds] >= $Q_genome_miss_end[$a]) #CONDITION 2: FETCHING PARTIALLY DELETED PROTEINS - CASE 1
					{
					my $pmr_len =0;
					my $cds_len =0;
					my $cds_percent;
					$pmr_len = $Q_genome_miss_end[$a] - $CDS_start[$cds] + 1;
					$cds_len = $CDS_end[$cds] - $CDS_start[$cds] + 1;
					$cds_percent = ($pmr_len * 100)/$cds_len;
					$trim_cds_percent = sprintf("%.2f", $cds_percent);
					if($trim_cds_percent >= 0.00)
						{
						if ($cds == 0)
                                                {
                                                         $tot_pmr_len = $tot_pmr_len + $pmr_len;
                                                         $pmr_tmp = $Q_genome_miss_end[$a];
                                                }
                                                elsif ($cds != 0 && $CDS_start[$cds] > $pmr_tmp)
                                                {

                                                        $tot_pmr_len = $tot_pmr_len + $pmr_len;
                                                        $pmr_tmp = $Q_genome_miss_end[$a];
                                                }
                                                elsif ($CDS_start[$cds] <=  $pmr_tmp && $Q_genome_miss_end[$a] >  $pmr_tmp)
                                                {
                                                         $tot_pmr_len = $tot_pmr_len + ($Q_genome_miss_end[$a] -  $pmr_tmp);
                                                         $pmr_tmp =  $Q_genome_miss_end[$a];
                                                }

 						$pmr = $pmr . $Q_genome_miss_start[$a] . "\t" . $Q_genome_miss_end[$a] . "\t" . $pmr_len . "\t" . $cds_len . "\t" . $CDS_start[$cds] . "\t" . $CDS_end[$cds] . "\t" . $trim_cds_percent . "\t" . $pdt . "\t" . $protein_id . "\t" . $note . "\n";
						$miss = 1;
                                                push(@prot_par_count, $protein_id);
						}
					}
				if ($CDS_start[$cds] < $Q_genome_miss_start[$a] && $CDS_end[$cds] >= $Q_genome_miss_start[$a] && $CDS_end[$cds] <= $Q_genome_miss_end[$a]) #CONDITION 3: FETCHING PARTIALLY DELETED PROTEINS - CASE 3
					{
					my $pmr_len =0;
					my $cds_len =0;
					my $cds_percent;
					$pmr_len = $CDS_end[$cds] - $Q_genome_miss_start[$a] + 1;
					$cds_len = $CDS_end[$cds] - $CDS_start[$cds] + 1;
					$cds_percent = ($pmr_len * 100)/$cds_len;
					$trim_cds_percent = sprintf("%.2f", $cds_percent);
					if($trim_cds_percent >= 0.00)
						{
						if ($cds == 0)
                                                {
                                                         $tot_pmr_len = $tot_pmr_len + $pmr_len;
                                                         $pmr_tmp = $CDS_end[$cds];
                                                }
                                                elsif ($cds != 0 && $Q_genome_miss_start[$a] > $pmr_tmp)
                                                {

                                                        $tot_pmr_len = $tot_pmr_len + $pmr_len;
                                                        $pmr_tmp = $CDS_end[$cds];
                                                }
                                                elsif ($Q_genome_miss_start[$a] <=  $pmr_tmp && $CDS_end[$cds] >  $pmr_tmp)
                                                {
                                                         $tot_pmr_len = $tot_pmr_len + ($CDS_end[$cds] -  $pmr_tmp);
                                                         $pmr_tmp = $CDS_end[$cds];
                                                }
						
						$pmr = $pmr . "\n" . $Q_genome_miss_start[$a] . "\t" . $CDS_end[$cds] . "\t" . $pmr_len . "\t" . $cds_len . "\t" . $CDS_start[$cds] . "\t" . $CDS_end[$cds] . "\t" . $trim_cds_percent . "\t" . $pdt . "\t" . $protein_id . "\t" . $note. "\n";
						$miss = 1;
                                                push(@prot_par_count, $protein_id);
						}
					}
				if ($CDS_start[$cds] <= $Q_genome_miss_start[$a] && $CDS_end[$cds] >= $Q_genome_miss_end[$a]) #CONDITION 4: FETCHING PARTIALLY DELETED PROTEINS - CASE 3
					{
					my $pmr_len =0;
					my $cds_len =0;
					my $cds_percent;
					$pmr_len = $Q_genome_miss_end[$a] - $Q_genome_miss_start[$a] + 1;
					$cds_len = $CDS_end[$cds] - $CDS_start[$cds] + 1;
					$cds_percent = ($pmr_len * 100)/$cds_len;
					$trim_cds_percent = sprintf("%.2f", $cds_percent);
					if($trim_cds_percent >= 0.00)
						{
						if ($cds == 0)
                                                {
                                                         $tot_pmr_len = $tot_pmr_len + $pmr_len;
                                                         $pmr_tmp =  $Q_genome_miss_end[$a];
                                                }
                                                elsif ($cds != 0 && $Q_genome_miss_start[$a] > $pmr_tmp)
                                                {

                                                        $tot_pmr_len = $tot_pmr_len + $pmr_len;
                                                        $pmr_tmp = $Q_genome_miss_end[$a];
                                                }
                                                elsif ($Q_genome_miss_start[$a] <=  $pmr_tmp &&  $Q_genome_miss_end[$a] >  $pmr_tmp)
                                                {
                                                         $tot_pmr_len = $tot_pmr_len + ($Q_genome_miss_end[$a] -  $pmr_tmp);
                                                         $pmr_tmp = $Q_genome_miss_end[$a];
                                                }

						$pmr = $pmr . "\n" . $Q_genome_miss_start[$a] . "\t" . $Q_genome_miss_end[$a] . "\t" . $pmr_len . "\t" . $cds_len . "\t" . $CDS_start[$cds] . "\t" . $CDS_end[$cds] . "\t" . $trim_cds_percent . "\t" . $pdt . "\t" . $protein_id . "\t" . $note. "\n";
						$miss = 1;
                                                push(@prot_par_count, $protein_id);
						}
					}
				}
			
		} 
 $cds_count++;
                 if($miss == 0)
                        {
                        $cpcr =  $cpcr . "\n" . $CDS_start[$cds] . "\t" . $CDS_end[$cds] . "\t" . $pdt . "\t" . $protein_id . "\t" . $note;
                        }
	
}
	$cds++;
	}

close IN;


$pmr =~ s/(^|\n)[\n\s]*/$1/g;
$pmr =~ s/^\s+//;
$cmr =~ s/(^|\n)[\n\s]*/$1/g;
$cmr =~ s/^\s+//;
$cpcr =~ s/^\s+//;
$cmr_count = $cmr =~ tr/\n//;
my $pmr_count = scalar(uniq(@prot_par_count));
my $complete_msg = "\n**********************Completely-missing-coding-regions**********************\n\n";

my @newcomplete_header = (MR_start, MR_end, CDS_start, CDS_end, MR_length, CDS_length, Product, Protein_id, Note);
my $complete_header = "MR_start\tMR_end\tCDS_start\tCDS_end\tMR_length\tCDS_Product\tProtein_id\tNote";
my $partial_header = "MR_start\tMR_end\tMR_length\tCDS_length\tCDS_start\tCDS_end\tCDS_Proportion(%)\tCDS_Product\tProtein_id\tNote \n";
my $partial_msg = "**********************Partially-missing-coding-regions********************** \n\n";
my $line_one = "Total completely missing coding regions length : $tot_cmr_len bp\n";
my $line_two = "Total partially missing coding regions length : $tot_pmr_len bp\n";
my $mr_len = $tot_cmr_len + $tot_pmr_len;
my $line_three = "Total missing regions length : $mr_len bp\n";


my $cpcr_pre = $output_prefix . "_cpcr.txt";
my $mcr_pre = $output_prefix . "_mcr.txt";

open (fh, ">", "$cpcr_pre");
print fh $cpcr;
close(fh) or "Couldn't close the file";


$fin_line = "$line_one$line_two$line_three $complete_msg $complete_header \n$cmr \n$partial_msg $partial_header$pmr";
# my $testout = join "\t", @newcomplete_header;
#print "$fin_line \n";

open (fh, ">", "$mcr_pre");
print fh $fin_line;
close(fh) or "Couldn't close the file";


exit;
