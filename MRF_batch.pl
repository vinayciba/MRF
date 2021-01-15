# This script takes mummer output and a reference GFF file as inputs and prints the complete and partial missing proteins.
# The program is divided into two parts, the first part gets the missing regions table from the mummer output.
# The second part fetches the coding sequences in the missing regions generated from the first part of the program. 

use strict;
#use warnings;
use List::Util ();
use List::MoreUtils qw(uniq);
use Getopt::Long;

my $help = "\n Usage: \n \t perl missing_regions_finder.pl [options] \n \n Options: 
	--inputfile 	 \t mummer output file 
	--genome_size 	 \t size of the reference genome 
	--gff3	 	 \t Reference gff3 file
	--cds_threshold  \t filter partial CDS [default: 0] 
	--fmatch_len 	 \t false match length [default: 15]
	--negoff 	 \t negative offset [default: 1]
	--posoff 	 \t positive offset [default: 1] \n";


my $queryName;
my $referenceName;
my $file;
my $genome_size;
my $gfffile;
my $cds_threshold;
my $fmatch_len;
my $negoff;
my $posoff;

GetOptions (    'queryName=s'           =>      \$queryName ,
		'referenceName=s'       =>      \$referenceName ,
		'inputfile=s'           =>      \$file ,
                'genome_size=i'         =>      \$genome_size ,
                'gff3:s'	        =>      \$gfffile ,
                'cds_threshold:f'       =>      \$cds_threshold ,
                'fmatch_len:i'          =>      \$fmatch_len ,
                'negoff:i'              =>      \$negoff ,
                'posoff:i'              =>      \$posoff
          );

my @words;			#Store mummer output in as array 
my @Ref_genome_start;		#Start position of reference match
my @Ref_genome_end;		#End position of reference match
my @Query_genome_start;		#Start position of query match
my @Query_genome_end;		#End position of query match
my @Non_hit_region;		#The sequence/region in the reference that has no corresponding match in the query
my @non_hit_qregion;		#The sequence/region in the query that has no match in the query
my @Q_genome_miss_start;	#Start position of the missing region of query in the reference genome's format
my @Q_genome_miss_end;		#End position of the missing region of query in the reference genome's format
my @Hit_region_length;		#Length of match
my $curr_line;			#stores iteration of missing region
my $fin_line;			#stores cumulative missing regions	
my $end_line;			#stores the last line of the missing region
my $tot_mr_len =0;		#Total missing regions length
my $genome_gap;			#The difference between the end coordinate and genome size of the reference
my $genome_fin_start;		#the start position of the last missing region when the the last match doesn't stretch to genome size
my $Ref_init_end;		#When reference match doesnt start at 1st bp, this variable is used to store the end positon of first missing region
my $cds_count =0;
my @prot_par_count;
#my $skip;
my @false_match;
my $string = "";
my $pdt_list = "#protein_name";
my $protein_list = "#protein_id";
my $pdt_lost = $queryName;;
if (! defined($file))
{
print $help;
}
else			#the remainder of code is wrapped in this block
{
$string = $file . ".fmt";

formatter($file); 	#Call to subroutine to format the mummer output file

sub formatter
{
	my ($file) = @_;
	BEGIN { $/ = "\n"; $\ = "\n"; }
	open (IN, $file) or die("mummer output file missing!!");
	truncate "$string", 0;
	while (defined($_ = <IN>)) 
	{
        	chomp $_;
		our(@F) = split(' ', $_, 0);
		open (fh, ">>", "$string");
		print fh  join("\t", @F) if $. > 1 and @F;
		close(fh) or "Couldn't close the file";
	}
}


chomp($genome_size);

open (IN, $string) or die("mummer output file missing!!");

# go through the file line by line

while (my $line = <IN>)
{
	my @words = split(/\t/, $line);      # split the current line on tabs
        push(@Ref_genome_start, $words[0]);
        push(@Query_genome_start, $words[1]);
        push(@Hit_region_length, $words[2]);
}

unless (defined $fmatch_len) {   $fmatch_len = 10; }
unless (defined $negoff) {       $negoff = 1; }
unless (defined $posoff) {       $posoff = 1; }

# get the count of items in the array list of Ref_genome_start
# and execute the loop for each item


for(my $i=0,my $j=1;$i<scalar(@Ref_genome_start);$i++,$j++)
{
	$Ref_genome_end[$i] = ($Ref_genome_start[$i] + $Hit_region_length[$i]) - 1;
        $Query_genome_end[$i] = ($Query_genome_start[$i] + $Hit_region_length[$i]) - 1;
	if($Hit_region_length[$i] <= $fmatch_len)	#if the hit length is less than the defined false match length, they will be screened in this block
	{
		screen_false_match($i, $j);
	}
	else
        {
                generate_mr($i, $j); 		# The main subroutine which generates missing regions 
        }
}

 if($Ref_genome_start[0] == 1)			#This block deals with first element of the array
{
                $Non_hit_region[0] = 0;
                $Q_genome_miss_start[0] = 0;
                $Q_genome_miss_end[0] = 0;
}
        my $count = scalar(@Ref_genome_start) - 1;
        my $count_act = scalar(@Ref_genome_start);

        $genome_gap = ($genome_size - $Ref_genome_end[$count]); 
        $genome_fin_start = $Ref_genome_end[$count] + 1;
        $Q_genome_miss_start[$count_act] = $genome_fin_start;
        $Q_genome_miss_end[$count_act] = $genome_size;
        $end_line =  "0\t0\t0\t0\t$genome_gap\t$genome_fin_start\t$genome_size"; 	#The last element of the array
        $fin_line =~ s/^(.*\n)1}//;
        if($genome_gap != 0)
                {
                $fin_line = $fin_line . "\n" . $end_line;
                }
# close the filehandle and exit
        close IN;



############################ SUBROUTINE BEGINS ########################
sub screen_false_match
{
	my ($i, $j) = @_;
	my $x;
	my $y;
	my $true_match = 0;
        $x = $i - 1;				#Negative offset
	$y = $i + 1;				#Positive offset
	my $neg = $i - $negoff;
	my $pos = $i + $posoff;
	while(1)
	{
		if($Query_genome_start[$x] <= $Query_genome_start[$i] &&  $Query_genome_start[$i] <= $Query_genome_start[$y])
		{						#if the above condtion is true, then it is not a false match
			if($x == $neg && $y == $pos)		#while loop exits if this condtion is true, the loop exits here
                        {
	                        generate_mr($i, $j);
	                        $true_match = 1;
                                last;
                        }
			if($x > $neg)
			{
				$x--;
			}
			if($y < $pos)
                        {
				$y++;
                        }
		}				
		else
		{
			last;
		}
	}	

		if($true_match == 0)			#false matches are eliminated in this loop
		{
#print "false match $i \n";
			push(@false_match, $Query_genome_start[$i]);
			splice(@Ref_genome_start, $i, 1);
			splice(@Ref_genome_end, $i, 1);
			splice(@Query_genome_start, $i, 1);
			splice(@Query_genome_end, $i, 1);
			splice(@Hit_region_length, $i, 1);
			 splice(@Non_hit_region, $i, 1);
			 splice(@non_hit_qregion, $i, 1);
			splice(@Q_genome_miss_start, $i, 1);
			 splice(@Q_genome_miss_end, $i, 1);
			$Ref_genome_end[$i] = ($Ref_genome_start[$i] + $Hit_region_length[$i]) - 1;
		        $Query_genome_end[$i] = ($Query_genome_start[$i] + $Hit_region_length[$i]) - 1;
			my $r = $i -1;
			$Non_hit_region[$i] = ($Ref_genome_start[$i] - $Ref_genome_end[$r]) - 1;
	                $non_hit_qregion[$i] = ($Query_genome_start[$i] - $Query_genome_end[$r]) - 1;
			$Q_genome_miss_start[$i] = $Ref_genome_end[$r] + 1;
	                $Q_genome_miss_end[$i] = ($Q_genome_miss_start[$i] + $Non_hit_region[$i]) - 1;
			if($Hit_region_length[$i] <= $fmatch_len)
			{
				screen_false_match($i, $j);
               	        }
	                else
               	        {
	                        generate_mr($i, $j);
       	                }


		}
}

################# SUBROUTINE screen_false_match($i, $j) ENDS #####################################


############################ SUBROUTINE BEGINS ########################
sub generate_mr
{
	my ($i, $j) = @_;
        $Non_hit_region[$j] = ($Ref_genome_start[$j] - $Ref_genome_end[$i]) - 1;
        $non_hit_qregion[$j] = ($Query_genome_start[$j] - $Query_genome_end[$i]) - 1;
        if($Non_hit_region[$i] == 1 && $non_hit_qregion[$i] ==1)		#masking single base mismatches
        {
		$Non_hit_region[$i] = "X";
                $Q_genome_miss_start[$i] = "-";
                $Q_genome_miss_end[$i] = "-";
         }
	 if($Non_hit_region[$i] > 0)						#masking negative hit regions
         {
                $tot_mr_len = $tot_mr_len + $Non_hit_region[$i];
         }
	if($i==0)
	{
		if($Ref_genome_start[0] ne 1)					#handling case where match doesn't start at 1st bp of REference genome
                {
                        $Ref_init_end = $Ref_genome_start[0] - 1;
                        $Non_hit_region[0] =  $Ref_init_end;
                        $Q_genome_miss_start[0] = 1;
                        $Q_genome_miss_end[0] =  $Ref_init_end;
                }
	}

	$Q_genome_miss_start[$j] = $Ref_genome_end[$i] + 1;
        $Q_genome_miss_end[$j] = ($Q_genome_miss_start[$j] + $Non_hit_region[$j]) - 1;
        $curr_line = "$Ref_genome_start[$i]\t$Ref_genome_end[$i]\t$Query_genome_start[$i]\t$Query_genome_end[$i]\t$Non_hit_region[$i]\t$Q_genome_miss_start[$i]\t$Q_genome_miss_end[$i]";
         $fin_line =  $fin_line . "\n" . $curr_line; 
		
	if($Non_hit_region[$i] <= 2)
               	{
                	$Q_genome_miss_start[$i] = "";
                	$Q_genome_miss_end[$i] = "";
                	$Non_hit_region[$i] = "";
               	}
}
	
################# SUBROUTINE generate_mr($i, $j) ENDS ####################################


####################### SECOND PART OF THE PROGRAM BEGINS##############
#Reads the gff3 file and the checks for the presence of coding sequences in the missing regions

if (defined($gfffile))
	{
open (IN, $gfffile) or die;		#Read reference GFF3 file
my @category;
my @CDS_start;
my @CDS_end;
my @strand;
my @trash;
my $karyotop = "chr – chr1 1 0" . " " . $genome_size . " " . "vlpurple";
my $karyotype;
my $features;
#my $value = 1;
my $sense = "";
my $antisense = "";
my $complete_highlights = "";
my $partial_highlights = "";
#my @product;				#
my $trim_cds_percent;			#Variable holding the percentage of CDS lost value
# my $mrc = 0;		
my $cds = 0;				#To keep the cds count 
my $cmr = "";				#Completely missing regions
my $pmr = "";				#Partially missing regions
my $cpcr = "";				#Completely present coding regions
my $cmr_tmp = 0;			#Temporary variable holding CDS end position
my $pmr_tmp = 0;			#Temporary variable holding CDS end postion
my $tot_cmr_len =0;
my $tot_pmr_len =0;
my @protein_count;
my @cpcr_count;
my @cpcs_count;
my $partial_label;
my @partial_lab;
	unless (defined $cds_threshold) {       $cds_threshold = 0.00; }
	while (my $line = <IN>)
		{											#While loop begins
		next if /^#/;
    		my @position = split(/\t/, $line);							# split the current line on tabs
	  	push(@category, $position[2]);
		push(@CDS_start, $position[3]);
	    	push(@CDS_end, $position[4]);
		push(@strand, $position[6]);
    		push(@trash, $position[8]);
		}											#While loop ends

	foreach my $k (@CDS_start)									#Loop through each row in gff3 file
		{
		my $miss = 0;
		if ($category[$cds] =~ m/^CDS$/)								#Verify it is an entry for CDS
			{
			if ($strand[$cds]  =~ m/\+/)
        			{
			        $sense = $sense . "chr1" .  " " . $CDS_start[$cds] . " " . $CDS_end[$cds] . " " . "blue" . "\n";
			        }
		        elsif($strand[$cds]  =~ m/\-/)
			        {
			        $antisense = $antisense . "chr1" .  " " . $CDS_start[$cds] . " " . $CDS_end[$cds] . " " . "green" . "\n";
			        }
			my ($pdt) = $trash[$cds] =~ m/product=([^=]+)\;/;				#Capture product
        		my ($protein_id) = $trash[$cds] =~ m/protein_id=([^=]+)[\;|\s]/;		#Capture protein id
			push(@protein_count, $protein_id);
	        	my ($note) = $trash[$cds] =~ m/Note=([^=]+\w+)\;/;				#capture note
			$karyotype = $karyotype . "band" . " " . "chr1" .  " " . "band" . $cds . " " . "band" . $cds . " " . $CDS_start[$cds] . " " . $CDS_end[$cds] . " " . "vdgrey" . "\n";
			$features = $features . "chr1" . " " . $CDS_start[$cds] . " " . $CDS_end[$cds] . "\n";
			my $pdt_fmtd = $pdt;
#			$pdt_list = $pdt_list . "\t" . $pdt_fmtd;
			$pdt_fmtd =~ s/\s+/-/g;
			$pdt_list = $pdt_list . "\t" . $pdt_fmtd;
			$protein_list =  $protein_list . "\t" . $protein_id;
			my $cds_tmp_num = -1;
			my $tmp_pmr_len;
			for(my $a=0;$a<=scalar(@Ref_genome_start);$a++)
				{
		 		if ($Q_genome_miss_start[$a] != "")
					{
					if ($CDS_start[$cds] >= $Q_genome_miss_start[$a] && $CDS_end[$cds] <= $Q_genome_miss_end[$a])
						{
						my $cmr_len =0;
						my $tmp_len =0;
						$cmr_len = mr_length($CDS_end[$cds], $CDS_start[$cds]);
						if ($cds == 0)
							{
							$tmp_len = $cmr_len;
							complete_mr($cmr_len, $CDS_end[$cds], $a, $pdt, $protein_id, $note, $tmp_len);
							}
						elsif ($cds != 0 && $CDS_start[$cds] > $cmr_tmp)
							{
							$tmp_len = $cmr_len;
							complete_mr($cmr_len, $CDS_end[$cds], $a,$pdt, $protein_id, $note, $tmp_len);
							}
						elsif ($CDS_start[$cds] <=  $cmr_tmp && $CDS_end[$cds] >  $cmr_tmp)
							{
							$tmp_len = ($CDS_end[$cds] -  $cmr_tmp);
							 complete_mr($cmr_len, $CDS_end[$cds], $a,$pdt, $protein_id, $note, $tmp_len);
							}
						 elsif ($CDS_start[$cds] <= $cmr_tmp && $CDS_end[$cds] <= $cmr_tmp)
                                                        {
                                                        $tmp_len = 0;
                                                        my $tmpCDS = $cds - 1;
                                                        complete_mr($cmr_len, $CDS_end[$tmpCDS], $a,$pdt, $protein_id, $note, $tmp_len);
                                                        }
							$complete_highlights =  $complete_highlights . "chr1" . " " . $CDS_start[$cds] . " " . $CDS_end[$cds] . " " . $pdt_fmtd . "\n";
						$miss = 1;
						}
					if ($CDS_start[$cds] >= $Q_genome_miss_start[$a] && $CDS_start[$cds] <= $Q_genome_miss_end[$a] && $CDS_end[$cds] >= $Q_genome_miss_end[$a])
						{
						my $pmr_len =0;
						my $cds_len =0;
						my $cds_percent;
#						my $tmp_pmr_len;
						$pmr_len = mr_length($Q_genome_miss_end[$a], $CDS_start[$cds]);
						$cds_len = mr_length($CDS_end[$cds], $CDS_start[$cds]);
						mr_percentage($pmr_len, $cds_len);
						if($trim_cds_percent >= $cds_threshold)
							{
							if ($cds == 0)
                        		                        {
							 	partial_mr($pmr_len, $Q_genome_miss_end[$a], $cds, $a, $trim_cds_percent, $pdt, $protein_id, $note);
                                                		}
	                                                elsif ($cds != 0 && $CDS_start[$cds] > $pmr_tmp)
                                                		{
								partial_mr($pmr_len, $Q_genome_miss_end[$a], $cds, $a, $trim_cds_percent, $pdt, $protein_id, $note);
        		                                        }
                        	                        elsif ($CDS_start[$cds] <=  $pmr_tmp && $Q_genome_miss_end[$a] >  $pmr_tmp)
                                		                {
								$pmr_len = ($Q_genome_miss_end[$a] -  $pmr_tmp);
								partial_mr($pmr_len, $Q_genome_miss_end[$a], $cds, $a, $trim_cds_percent, $pdt, $protein_id, $note);
                                                		}
							$pmr = $pmr . $Q_genome_miss_start[$a] . "\t" . $Q_genome_miss_end[$a] . "\t" . $CDS_start[$cds] . "\t" . $CDS_end[$cds] . "\t" . $pmr_len . "\t" . $cds_len . "\t" . $trim_cds_percent . "\t" . $pdt . "\t" . $protein_id . "\t" . $note . "\n";
							if($cds_tmp_num != $cds)
								{
								$pdt_lost = $pdt_lost . "\t" . $pmr_len;
								$cds_tmp_num = $cds;
								$tmp_pmr_len = $pmr_len;
								}
							else
								{
								$tmp_pmr_len =  $tmp_pmr_len + $pmr_len;
								my @trim = split(/\t([^\t]+)$/, $pdt_lost);
								$pdt_lost = $trim[0];
								$pdt_lost = $pdt_lost . "\t" . $tmp_pmr_len;
								}
							if($trim_cds_percent >= $cds_threshold)
								{
							$partial_highlights = $partial_highlights . "chr1" . " " . $Q_genome_miss_start[$a] . " " . $Q_genome_miss_end[$a] . " " . $pdt_fmtd . "\n";
							$partial_label =  "chr1" . " " . $CDS_start[$cds] . " " . $CDS_end[$cds] . " " .  $pdt_fmtd . "\n";
							push(@partial_lab, $partial_label);
								}
							$miss = 1;
							push(@prot_par_count, $protein_id);
							}
						}
					if ($CDS_start[$cds] < $Q_genome_miss_start[$a] && $CDS_end[$cds] >= $Q_genome_miss_start[$a] && $CDS_end[$cds] <= $Q_genome_miss_end[$a])
						{
						my $pmr_len =0;
						my $cds_len =0;
						my $cds_percent;
#						my $tmp_pmr_len;
						$pmr_len = mr_length($CDS_end[$cds], $Q_genome_miss_start[$a]);
						$cds_len = mr_length($CDS_end[$cds], $CDS_start[$cds]);
						mr_percentage($pmr_len, $cds_len);
						if($trim_cds_percent >=  $cds_threshold)
							{
							if ($cds == 0)
                                        		        {
								partial_mr($pmr_len, $CDS_end[$cds], $cds, $a, $trim_cds_percent, $pdt, $protein_id, $note);
                		                                }
                                	                elsif ($cds != 0 && $Q_genome_miss_start[$a] > $pmr_tmp)
                                        		        {
							 	partial_mr($pmr_len, $CDS_end[$cds], $cds, $a, $trim_cds_percent, $pdt, $protein_id, $note);
		                                                }
                	                                elsif ($Q_genome_miss_start[$a] <=  $pmr_tmp && $CDS_end[$cds] >  $pmr_tmp)
                        		                        {
								$pmr_len = ($CDS_end[$cds] -  $pmr_tmp);
								partial_mr($pmr_len, $CDS_end[$cds], $cds, $a, $trim_cds_percent, $pdt, $protein_id, $note);
                                                		}
						
							$pmr = $pmr . "\n" . $Q_genome_miss_start[$a] . "\t" . $Q_genome_miss_end[$a] . "\t" . $CDS_start[$cds] . "\t" . $CDS_end[$cds] . "\t" . $pmr_len . "\t" . $cds_len . "\t" . $trim_cds_percent . "\t" . $pdt . "\t" . $protein_id . "\t" . $note. "\n";
							 if($cds_tmp_num != $cds)
                                                                {
                                                                $pdt_lost = $pdt_lost . "\t" . $pmr_len;
                                                                $cds_tmp_num = $cds;
                                                                $tmp_pmr_len = $pmr_len;
                                                                }
                                                        else
                                                                {
                                                                $tmp_pmr_len =  $tmp_pmr_len + $pmr_len;
								my @trim = split(/\t([^\t]+)$/, $pdt_lost);
                                                                $pdt_lost = $trim[0];
                                                                $pdt_lost = $pdt_lost . "\t" . $tmp_pmr_len;
                                                                }
							if($trim_cds_percent >= $cds_threshold)
								{
							$partial_highlights = $partial_highlights . "chr1" . " " .  $Q_genome_miss_start[$a] . " " .  $Q_genome_miss_end[$a] . " " . $pdt_fmtd . "\n";
							 $partial_label =  "chr1" . " " . $CDS_start[$cds] . " " . $CDS_end[$cds] . " " .  $pdt_fmtd . "\n";
							  push(@partial_lab, $partial_label);
								}
							$miss = 1;
							push(@prot_par_count, $protein_id);
							}
						}
					if ($CDS_start[$cds] <= $Q_genome_miss_start[$a] && $CDS_end[$cds] >= $Q_genome_miss_end[$a])
						{
						my $pmr_len =0;
						my $cds_len =0;
						my $cds_percent;
#						my $tmp_pmr_len;
						$pmr_len = mr_length($Q_genome_miss_end[$a], $Q_genome_miss_start[$a]);
						$cds_len = mr_length($CDS_end[$cds], $CDS_start[$cds]);
						mr_percentage($pmr_len, $cds_len);
						if($trim_cds_percent >= $cds_threshold)
							{
							if ($cds == 0)
                                        		        {
							   	partial_mr($pmr_len, $Q_genome_miss_end[$a], $cds, $a, $trim_cds_percent, $pdt, $protein_id, $note);
                                                		}
	                                                elsif ($cds != 0 && $Q_genome_miss_start[$a] > $pmr_tmp)
        		                                        {
								partial_mr($pmr_len, $Q_genome_miss_end[$a], $cds, $a, $trim_cds_percent, $pdt, $protein_id, $note);
                                		                }
                                                	elsif ($Q_genome_miss_start[$a] <=  $pmr_tmp &&  $Q_genome_miss_end[$a] >  $pmr_tmp)
		                                                {
								$pmr_len = ($Q_genome_miss_end[$a] -  $pmr_tmp);
								 partial_mr($pmr_len, $Q_genome_miss_end[$a], $cds, $a, $trim_cds_percent, $pdt, $protein_id, $note);
                                		                }
							$pmr = $pmr . "\n" . $Q_genome_miss_start[$a] . "\t" . $Q_genome_miss_end[$a] . "\t" . $CDS_start[$cds] . "\t" . $CDS_end[$cds] . "\t" . $pmr_len . "\t" . $cds_len . "\t" . $trim_cds_percent . "\t" . $pdt . "\t" . $protein_id . "\t" . $note. "\n";
							 if($cds_tmp_num != $cds)
                                                                {
                                                                $pdt_lost = $pdt_lost . "\t" . $pmr_len;
                                                                $cds_tmp_num = $cds;
                                                                $tmp_pmr_len = $pmr_len;
                                                                }
                                                        else
                                                                {
                                                                $tmp_pmr_len =  $tmp_pmr_len + $pmr_len;
								my @trim = split(/\t([^\t]+)$/, $pdt_lost);
                                                                $pdt_lost = $trim[0];
                                                                $pdt_lost = $pdt_lost . "\t" . $tmp_pmr_len;
                                                                }
							if($trim_cds_percent >= $cds_threshold)
								{
							$partial_highlights = $partial_highlights . "chr1" . " " . $Q_genome_miss_start[$a] . " " . $Q_genome_miss_end[$a] . " " . $pdt_fmtd . "\n";
							 $partial_label =  "chr1" . " " . $CDS_start[$cds] . " " . $CDS_end[$cds] . " " .  $pdt_fmtd . "\n";
							 push(@partial_lab, $partial_label);
								}
							$miss = 1;
							push(@prot_par_count, $protein_id);
							}
					}	#if mr is true block ends
				}		#for loop ends
			}			#if cds is true block ends
			$cds_count++; 
		 	if($miss == 0)
                        	{
	                        $cpcr =  $cpcr . "\n" . $CDS_start[$cds] . "\t" . $CDS_end[$cds] . "\t" . $pdt . "\t" . $protein_id . "\t" . $note;
				 push(@cpcr_count, $protein_id);
				$pdt_lost = $pdt_lost . "\t" . "0";
        	                }
		}				#foreach loop ends
		$cds++;
	}					#if gff3 is true block ends

close IN;


sub complete_mr
	{
	my ($cmr_len, $cmrtmp, $a, $pdt, $protein_id, $note, $tmp_len) = @_;
	$tot_cmr_len = $tot_cmr_len + $tmp_len;
        $cmr_tmp = $cmrtmp;
	$cmr = $cmr . "\n" . $Q_genome_miss_start[$a] . "\t" . $Q_genome_miss_end[$a] .  "\t" . $CDS_start[$cds] . "\t" . $CDS_end[$cds] . "\t" . $cmr_len . "\t" . $cmr_len . "\t" . "100" . "\t" . $pdt . "\t" . $protein_id . "\t" . $note . "\n";
	$pdt_lost = $pdt_lost . "\t" . $tmp_len;
	}

sub partial_mr
        {
	 my ($pmr_len, $pmrtmp, $cds, $a, $trim_cds_percent, $pdt, $protein_id, $note) = @_;
	 $tot_pmr_len = $tot_pmr_len + $pmr_len;
	 $pmr_tmp = $pmrtmp;
        }

sub mr_length
	{
	my ($x, $y) = @_;
	my $z = ($x - $y) + 1;
	return $z;
	}

sub mr_percentage
	{
	my $cds_percent;
	 my ($pmr_len, $cds_len) = @_;
	 $cds_percent = ($pmr_len * 100)/$cds_len;
         $trim_cds_percent = sprintf("%.2f", $cds_percent);
	}

my $karyotype_final =  $karyotop . "\n" . $karyotype;

open(OUT,'>','circos_files/karyotype.txt');
print OUT $karyotype_final;
close OUT;

open(OUT,'>','circos_files/features_high.txt');
print OUT $features;
close OUT;

open(OUT,'>','circos_files/sense.txt');
print OUT $sense;
close OUT;

open(OUT,'>','circos_files/antisense.txt');
print OUT $antisense;
close OUT;

open(OUT,'>','circos_files/cmr_high.txt');
print OUT $complete_highlights;
close OUT;

open(OUT,'>','circos_files/pmr_high.txt');
print OUT $partial_highlights;
close OUT;

#$partial_label = join "\n", List::Util::uniqstr( split /\n/, $partial_label);
my @uniq_partial_lab = uniq(@partial_lab);
$partial_label = join("", @uniq_partial_lab);
my $label_highlights =  $complete_highlights . $partial_label;

open(OUT,'>','circos_files/label_high.txt');
print OUT $label_highlights;
close OUT;

my $prot_count = scalar(uniq(@protein_count));
my $cpcs_count = scalar(uniq(@cpcr_count));
$pmr =~ s/(^|\n)[\n\s]*/$1/g;
$pmr =~ s/^\s+//;
$cmr =~ s/(^|\n)[\n\s]*/$1/g;
$cmr =~ s/^\s+//;
$cpcr =~ s/^\s+//;
my $cmr_count = $cmr =~ tr/\n//;
my $pmr_count = scalar(uniq(@prot_par_count));
#my $cpcs_count =  $cpcr =~ tr/\n//;
#$cpcs_count =  $cpcs_count + 1;
my $complete_msg = "**********************Completely-missing-coding-regions**********************\n\n";

#my @newcomplete_header = (MR_start, MR_end, CDS_start, CDS_end, MR_length, CDS_length, Product, Protein_id, Note);
my $complete_header = "MR_start\tMR_end\tCDS_start\tCDS_end\tMR_length\tCDS_Product\tProtein_id\tNote";
my $partial_header = "MR_start\tMR_end\tMR_length\tCDS_length\tCDS_start\tCDS_end\tCDS_Proportion(%)\tCDS_Product\tProtein_id\tNote \n";
my $partial_msg = "**********************Partially-missing-coding-regions********************** \n\n";
my $cds_num = "No. of coding sequences in the reference genome: $prot_count \n";
my $cmr_num = "No. of completely missing coding sequences: $cmr_count \n";
my $pmr_num = "No. of partially missing coding sequences: $pmr_count \n";
my $cpcs_num = "No. of coding sequences with no deletion: $cpcs_count \n";
my $line_one = "Length of completely deleted coding sequences : $tot_cmr_len bp\n";
my $line_two = "Length of partially deleted coding sequences : $tot_pmr_len bp\n";
my $mr_len = $tot_cmr_len + $tot_pmr_len;
my $line_three = "Total deleted coding sequences length : $mr_len bp\n\n";
my $mr_header = "Ref_start \t Ref_end \t Query_start \t Query_end \t Non_hit_region \t Query_miss_start \t Query_miss_end \n";
my $summary_hdr = "#accession \t Num CDS in Ref. \t Num completely missing CDS \t Num Partially missing CDS \t Completely missing CDS length \t Partially missing CDS length \t Total missing CDS length\n";
my $summary = "$queryName \t $prot_count \t $cmr_count \t $pmr_count \t $tot_cmr_len \t $tot_pmr_len \t $mr_len";
my $summary_out = $summary_hdr . $summary;
#my $cpcs_count =  $cpcr =~ tr/\n//;
my $cpcs = $file . ".cpcr";
open (fh, ">", "$cpcs");
print fh $cpcr;
close(fh) or "Couldn't close the file";

my $complete_mcs = "";
my $incomplete_mcs = "";
$complete_mcs = $cds_num . $cmr_num . $pmr_num . $cpcs_num . $line_one . $line_two . $line_three . $complete_msg . $complete_header . "\n" . $cmr . "\n" . $partial_msg . $partial_header . $pmr;

$incomplete_mcs = $cmr . $pmr;

my $mcr = $file . ".mcr";

#print $complete_mcs;
open (fh, ">", "$mcr");
print fh $incomplete_mcs;
close(fh) or "Couldn't close the file";

my $summaryfile = $file . ".summary";

open (fh, ">", "$summaryfile");
print fh $summary_out;
close(fh) or "Couldn't close the file";

$pdt_list = $pdt_list . "\t" . "total";
my $pdt_lost = $pdt_lost . "\t" . $mr_len;
my $pdtLostList =  $pdt_list . "\n" . $protein_list . "\n" . $pdt_lost;

my $summaryCDS = $file . ".cds";

open (fh, ">", "$summaryCDS");
print fh $pdtLostList;
close(fh) or "Couldn't close the file";

}

#print "$pdt_list \n  $pdt_lost \n";
#print "$fin_line \n";
my $mr = $file . ".mr";
my $mr_header = "Ref_start \t Ref_end \t Query_start \t Query_end \t Non_hit_region \t Query_miss_start \t Query_miss_end";
my $complete_mr = $mr_header . "\n" . $fin_line;

open (fh, ">", "$mr");
print fh $complete_mr;
close(fh) or "Couldn't close the file";



my $mj_spacing;
my $mj_color = "red";
my $mj_ls = "33p";
my $mj_format;
my $mi_spacing;
my $mi_color = "black";
my $mi_ls = "28p";
my $mi_format;
my $majortick;
my $minortick;
if ($genome_size <= 10000)
        {
         $mj_spacing = "0.05u";
         $mj_format = "%.1f";
         $mi_spacing = "0.02u";
         $mi_format = "%0.01f";
        }
if ($genome_size > 10000 && $genome_size <= 100000)
        {
         $mj_spacing = "0.1u";
         $mj_format = "%d";
         $mi_spacing = "0.05u";
         $mi_format = "%0.1f";
        }
if ($genome_size > 100000)
        {
         $mj_spacing = "1u";
         $mj_format = "%d";
         $mi_spacing = "0.5u";
         $mi_format = "%d";
        }
$majortick = tick($mj_spacing, $mj_color, $mj_ls, $mj_format);
$minortick = tick($mi_spacing, $mi_color, $mi_ls, $mi_format);

sub tick
{
my ($spc, $color, $ls, $format) = @_;


return "<tick>
spacing = $spc
color = $color
show_label = yes
label_size = $ls
label_font = bold
label_offset = 0p
format = $format
</tick>";

}
sub generateticks
{
return "
show_ticks = yes
show_tick_labels = yes

<ticks>
radius = dims(ideogram,radius_outer)
multiplier = 0.001

label_offset = 5p
thickness = 3p
size = 20p

label_separation = 5p

 $majortick \n
 $minortick

</ticks>
"
}

my $generate_ticks = generateticks();

#print $generate_ticks;

open(OUT,'>','circos_files/ticks.conf');
print OUT $generate_ticks;
close OUT;


my $mircos_filename = $file . ".MirCos.png"; 

my $printCircosConf = genCircosConf();

open(OUT,'>','circos_files/circos.conf');
print OUT $printCircosConf;
close OUT;


sub genCircosConf
{

return "
karyotype = /your/working/directory/circos_files/karyotype.txt

chromosomes_units = 10000
chromosomes_display_defaults = yes
chromosomes_radius = chr1:0.90r;

<<include /your/working/directory/circos_files/high.conf>>

<image>
dir = /your/working/directory/circos_files
file  = $mircos_filename
# radius of inscribed circle in image
radius         = 1800p
background     = white
# by default angle=0 is at 3 o'clock position
angle_offset   = -90
</image>



<ideogram>
<spacing>

default = 0u
break = 0u

</spacing>


thickness = 80p

radius = 0.80r
show_label = no
label_font = bold
label_with_tag = yes
label_radius = dims(ideogram,radius) + 0.05r
label_size = 48
labell_parallel = yes
label_case = upper

stroke_thickness = 3
stroke_color = white
fill = yes

show_bands = yes
fill_bands = yes

</ideogram>

<<include /your/working/directory/circos_files/ticks.conf>>
<plots>
<plot>
type = text
color = black
label_font = normal
label_size = 35p

file = /your/working/directory/circos_files/label_high.txt

r1 = 1.50r
r0 = 1.08r

padding = 1p
rpadding = 1p
show_links = yes
link_dims = 1p,2p,3p,2p,1p,1p
link_thickness = 3p
link_color = red
label_snuggle = yes
max_snuggle_distance = 2r
snuggle_sampling = 1
snuggle_tolerance = 0.25r
snuggle_link_overlap_test = yes
snuggle_link_overlap_tolerance = 2p
snuggle_refine = yes

</plot>

</plots>

<<include /your/circos//installation/direcotry/etc/colors_fonts_patterns.conf>>

<<include /your/circos/installation/direcotry/etc/housekeeping.conf>>
"
}

}

exit;

