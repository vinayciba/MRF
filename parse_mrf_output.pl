#!/usr/bin/perl
use List::MoreUtils qw(uniq);
use List::MoreUtils qw(first_index);
use List::Util qw(first);
use Getopt::Long;
use Data::Dumper;
my $usage = "\n Usage: \n \t perl parse_mrf_output.pl [options] \n \n Options:
        --prot_ID      \t search by protein ids
        --prot_name    \t search by protein names 
	--filter    \t filter by missing coding sequences length [default: 0]
	--number    \t number of missing coding sequences to show in heatmap [default: 15]
	--first     \t list top n proteins
	--last     \t list last n proteins
	--range	   \t list by range[m,n]
	--all	   \t list all proteins
	--top	   \t lists only top n gehomes [default: 10]
        --help	       \t show this help		\n";


my $top=10;
my $protID;
my $protNames;
my $filter=0;
my $number=15;
my $head;
my $tail;
my $range;
my $all;
my $help;



GetOptions (    'prot_ID=s'           =>      \$protID ,
                'prot_name=s'       =>      \$protNames,
		'filter=i'	=>	\$filter,
		'number=i'	=>	\$number,
		'first=i'	=>	\$head,
		'last=i'	=> 	\$tail,
		'range=s'	=>	\$range,
		'all!'		=>	\$all,
		'top=i'		=>	\$top,
		'help!'	 =>	\$help
          );

my @summaryfiles = <*.summary>;
my @cdsfiles = <*.cds>;
my @summaryAll;
my @cdslist;
my @cdsAll;
my @indexes;
my @totColumn;
#if(defined ($help))
#{
#print $help;
#}

if ($help) {
    print $usage;
    exit(1);
}


#if (scalar(@ARGV) < 2)
#{
#    print "No valid arguments!\n";
#    print $usage;
#    exit(1);
#}


foreach my $file (@summaryfiles)
	{
	open (IN, $file) or die("$file file not found!!");
	while (my $line = <IN>)
		{
		next if $line =~ m/^#/;
		push(@summaryAll, $line);
		}
	close IN;
	}


foreach my $cdsfile (@cdsfiles)
        {
        open (IN, $cdsfile) or die("$cdsfile file not found!!");
        while (my $line = <IN>)
                {
                if ($line =~ m/^#/)
			{
                	push(@cdslist, $line);
			}
		else
			{
			 push(@cdsAll, $line);
			}
                }
        close IN;
        }

my $summary_hdr = "accession \t Num_CDS_in_Ref. \t Num_completely_missing_CDS \t Num_Partially_missing_CDS \t Completely_missing_CDS_length \t Partially_missing_CDS_length \t Total_missing_CDS_length\n";
my @sortedSumAll = sort { (split(' ', $b))[6] <=> (split(' ', $a))[6] } @summaryAll;
unshift(@sortedSumAll, $summary_hdr);

my @cdsNoHdr = @cdsAll;
my @cdsOneHdr = @cdsAll;
my $cds_hdr1 = @cdslist[0];
my $cds_hdr2 = @cdslist[1];
unshift(@cdsAll, $cds_hdr2);
unshift(@cdsAll, $cds_hdr1);

my $cds_hdr1NC = $cds_hdr1;
$cds_hdr1NC =~ s/^#//;
unshift(@cdsOneHdr, $cds_hdr1NC);

my @protCount = split(/\t/, $cds_hdr2);
my $totProteins = (scalar(@protCount) - 1);
my $columnCount = scalar(@protCount);
#print "total prots: $totProteins";

open(OUT,'>','summaryAll.txt');
print OUT  @sortedSumAll;
close OUT;

open(OUT,'>','cdsAll.txt');
print OUT  @cdsAll;
close OUT;

if ($all) {
    print @cdsAll;
    exit(1);
}


 foreach my $entry (@cdsNoHdr)
        {
                my @fields = split(/\t/, $entry);
                push @totColumn,$fields[$columnCount];
        }
my $totGenomes = scalar(@totColumn);
my %genomeIndexHash;
for(my $t=0;$t<$totGenomes;$t++)
{
	$genomeIndexHash{$t} = $totColumn[$t];	
}

#print values  %genomeIndexHash;
#my @sortedTotColumn = reverse sort { $a <=> $b } @totColumn;
#print "@sortedTotColumn \n";
my @sortedGenColIndex;
foreach my $nm (sort {  $genomeIndexHash{$b} <=>  $genomeIndexHash{$a} } keys %genomeIndexHash)
{
	push @sortedGenColIndex,$nm;
#print "$nm :  $genomeIndexHash{$nm} \n";       
}

my $NUM = 1;
my @cdsTop;
foreach my $idx (@sortedGenColIndex)
{
#print "$idx \t @cdsNoHdr[$idx]\n";

if($NUM <= $top)
{
push @cdsTop,  @cdsNoHdr[$idx]
}
$NUM++;
}
my @cdsTopOneHdr = @cdsTop;
unshift(@cdsTop, $cds_hdr2);
unshift(@cdsTop, $cds_hdr1);
unshift(@cdsTopOneHdr, $cds_hdr1NC);
#print "@cdsTop \n";

##remove cds with zero deletion
my @allIndex;
my @filtered;
my @retain;
my $j=1;
my %colAvgHash;
while($j<=$totProteins)
{
	my @column;
	foreach my $entry (@cdsNoHdr)
	{
		my @fields = split(/\t/, $entry);
		push @column,$fields[$j];
	}
	my $colAvg = average(@column);
#print "$colAvg Avg\n";
	if($colAvg > 0)
	{
		$colAvgHash{$j} = $colAvg;
		 push @retain,$j;
	}

	@avgVals = values %colAvgHash;
#print "@avgVals \n";

	if($filter > 0)
	{
		foreach my $value (@column)
		{
			if ($value >= $filter)
			{
				push @filtered,$j;
#			print "$value \t $filter \t $j \n";
				last;
			} 
		}
	push @allIndex, $j;
	}
	$j++;
}



sub average
{
	my @array = @_;
	my $sum; # create a variable to hold the sum of the array's values
	foreach (@array)
	{
		$sum += $_; # add each element of the array
	} 
	return $sum/scalar(@array); # divide sum by the number of elements in the
}

#@avgVals = values %colAvgHash;
#print "@avgVals \n";
#print Dumper(\%colAvgHash);

my $x = 0;
my @tempIndexes;
foreach my $name (sort { $colAvgHash{$b} <=> $colAvgHash{$a} } keys %colAvgHash)
{
	if($x<$number)
	{
#    printf "%-8s %s\n", $name, $colAvgHash{$name};
		push @tempIndexes,$name;
#sort @indexes;
#print "@indexes \n";
	}
	$x++;
}

@indexes = sort @tempIndexes;
#print @indexes;

#print "@indexes \n";
#get index/indices of $protein_id or $protein_name
#get columns of @cdsAll based on index
#print those columns as  a single array

#my @names = split(/\t/, $cds_hdr1);
#my @ids = split(/\t/, $cds_hdr2);
#my @nameSearch = split(/\,/, $protNames);
#my @idSearch = split(/\,/, $protID);

#getIndexOfIds(@idSearch);
#getIndexOfNames(@nameSearch);


#my @indexes;

if(defined($protID))
{
@indexes = ();
my @ids = split(/\t/, $cds_hdr2);
my @idSearch = split(/\,/, $protID);
getIndexOfIds(\@idSearch, \@ids);
}

if(defined($protNames))
{
@indexes = ();
my @names = split(/\t/, $cds_hdr1);
my @nameSearch = split(/\,/, $protNames);
getIndexOfNames(\@nameSearch, \@names);
}


if($filter > 0)
{
#print "filter \n";

###This hash can be used to filter one array from another
#my %indices;
#print "@purge \n";
#@indices{ @allIndex } = ();            # All files are the keys.
#delete @indices{ @purge }; # Remove the links.
#@indexes = keys %indices;
#print "@indexes \n";
@indexes = @filtered;
}



if(defined($head))
{
	@indexes = ();
	if($head <= $totProteins)
	{
		my $i=1;
		while($i<=$head)
		{
			 push @indexes,$i;
			 $i++;
		}
	}
	else
	{
                my $i=0;
                while($i<$totProteins)
                {
                         push @indexes,$i;
                         $i++;
                }
        }
}

if(defined($tail))
{
@indexes = ();
my $invTail = -$tail;
$invTail--;
my $invTotProteins = -$totProteins;
        if($tail <= $totProteins)
        {
                my $i=-2;
                while($i>=$invTail)
                {
                         push @indexes,$i;
                         $i--;
                }
        }
        else
        {
                my $i=-i;
                while($i>$totProteins)
                {
                         push @indexes,$i;
                         $i--;
                }
        }
}

if(defined($range))
{
	@indexes = ();
	my @ranges = split(/\,/, $range);
	my $rangeStart = $ranges[0];
	my $rangeEnd = $ranges[1];
	if($rangeEnd > $totProteins)
	{
		print "invalid range!!";
	}
	else
	{
		while($rangeEnd >= $rangeStart)
                {
                         push @indexes,$rangeStart;
                         $rangeStart++;
                }

	}
}
sub getIndexOfIds
{
	my ($q,$i) = @_;
	my @queries = @$q;
	my @ids = @$i;
	foreach my $place (@queries)
	{
		push @indexes, first_index { $_ eq $place } @ids;
	}

}


sub getIndexOfNames
{
        my ($q,$i) = @_;
	my @queries = @$q;
        my @names = @$i;
        foreach my $place (@queries)
        {
#		 push @indexes, first_index { $_ eq $place } @names;
                push @indexes, first_index { $_  =~ m/\Q$place\E/ } @names;
        }

}

#my @indexes;
#foreach my $place (@idSearch) {
#    push @indexes, first_index { $_ eq $place } @ids;
#}
#print "\n";
#print @indexes;


my $sliced;
#my $filtered;
while($row=shift(@cdsTopOneHdr))
{

my @words = split(/\t/, $row);
$sliced .= $words[0];
foreach my $column (@indexes) {

$sliced .= "\t" . $words[$column];

}
$sliced = $sliced . "\n";

#$filtered .= "\t" . $words[0];
foreach my $column (@retain) {

$filtered .= "\t" . $words[$column];

}
#$filtered = $filtered . "\n";
}

#print "$sliced \n";
#print "$filtered \n";

open(OUT,'>','cdsHeatmap.txt');
print OUT  $sliced;
close OUT;

my $screenOut;
while($row=shift(@cdsTop))
{
	my @words = split(/\t/, $row);
	$screenOut .= $words[0];
	foreach my $column (@indexes)
	{
		$screenOut .= "\t" . $words[$column];
	}
	$screenOut = $screenOut . "\n";

}

print "$screenOut \n";
