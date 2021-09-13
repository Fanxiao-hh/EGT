#perl AI_time.pl -t timeline-final -r rankorder -fp dir -AI [num] -hu [num] -term

use warnings;
use strict;
use Getopt::Long;

our $timeline; 
our $order;
our $rankcoverage;
our $filepath;
our $ai=10;
our $hu=30;
our $chs=80;
our $term;


our %hashname;
our %hashtime;
our %hashrank;
our %hashfathertaxonid;
our %hashskip;
our %hashinclude;
our %hashtransfer;
our %hashtransfergene;
our %hashterm;
our %hashnodeinfor;
our %hashcoverage;

if($ARGV[0]=~m/-h/) {&Usage;exit;}

GetOptions(
"timeline|t:s"     => \$timeline,
"rankorder|r:s"  => \$order,
"rankcoverage|rc:s" => \$rankcoverage,
"filepath|fp:s"  => \$filepath,
"AI:s" => \$ai,
"hu:s" => \$hu,
"chs:s" => \$chs,
"term:s" => \$term,
);

unless (defined $timeline){&Usage;exit;}
unless (defined $order){&Usage;exit;}
unless (defined $rankcoverage){&Usage;exit;}
unless (defined $filepath){&Usage;exit;}
unless (defined $term){&Usage;exit;}
unless (defined $chs){&Usage;exit;}
open IN,"< $timeline" or die $!;#open time line file;
while(<IN>){
	chomp;
	my $line=$_;
	if(!/^\#/){
	my @timelinefile=split /\t/, $line;
	my $rank=$timelinefile[0]; my $name=$timelinefile[1]; my $taxonid=$timelinefile[2]; my $time=$timelinefile[3];
	$hashname{$taxonid}=$name;
	$hashtime{$taxonid}=$time;
	$hashrank{$taxonid}=$rank;
#print $name."\t".$taxonid."\t".$time."\n";
	}
}

close IN;

open FI,"< $order" or die $!;
while(<FI>){
	next if (/^#/);
chomp;
my $line=$_;
my @elements=split /\t/,$line;
my $level=$elements[0];
my $skiptaxonid=$elements[1];
my $skipname=$elements[2];
my $skiprank=$elements[3];
my $includetaxonid=$elements[4];
my $includename=$elements[5];
my $includerank=$elements[6];
my $speice=$elements[7];
my $rankcover=$elements[9];


$hashname{$skiptaxonid}=$skipname;
$hashname{$includetaxonid}=$includename;
$hashrank{$skiptaxonid}=$skiprank;
$hashrank{$includetaxonid}=$includerank;
$hashfathertaxonid{$skiptaxonid} = $includetaxonid;
$hashskip{$level}{$speice}=$skiptaxonid;
$hashinclude{$level}{$speice}=$includetaxonid;
$hashcoverage{$level}{$speice}=$rankcover;
#print $speice."\n";
}
close FI;

foreach my $level (sort {$a <=> $b} keys %hashinclude){
my $subhash =$hashinclude{$level};
	foreach my $speice (sort keys %{$subhash}){
		my $taxonid=();
		my $fathertaxonid=(); 
		my $time=(); 
		my $fathertime=();
		#print $level."\t".$speice."\n";
		if(exists $hashinclude{$level}{$speice}){
			
			$taxonid=$hashinclude{$level}{$speice};
			if($taxonid=~m/131567/){
				$fathertime="NA";
				$time="NA";
				}else{
					$time=$hashtime{$taxonid};
			$fathertaxonid=$hashfathertaxonid{$taxonid};
			#print $time."\t".$taxonid."\t".$fathertaxonid."\n";
			if(exists $hashtime{$fathertaxonid}){
#print $taxonid."\t".$fathertaxonid."\n";
				$fathertime=$hashtime{$fathertaxonid};
			while($fathertime=~m/NA/){
				if($fathertaxonid =~m/131567/){
					$fathertime="NA";
					}else{
					$fathertaxonid=$hashfathertaxonid{$fathertaxonid};
					$fathertime=$hashtime{$fathertaxonid};
				 }
				}
			}else{ 
				print "nohash for ".$taxonid."\t".$fathertaxonid."\thashfathertime\n";
			}
		}
			}else{
				print "no"."\t".$level."\t".$speice."\n";
			}
			
		
			open TER," <$term" or print "can not find termfile\n";
			while(<TER>){
				chomp;
				my $line=$_;
				$hashterm{$line}=1;
			}
			if($hashcoverage{$level}{$speice} > $rankcoverage){
			my $file= $filepath."/"."level-".$level."-".$speice.".HGT.txt.filter.filterout";
			if(-e $file ){
		open FHI, " <$file" or print "can not find $file\n";
		while(<FHI>){
			chomp;
			my $line=$_;
			if($line=~m/\#/){}elsif($line=~m/(\S+).*/){
				my @hgtele=split /\t/, $line;
				#my $num=@hgtele;if($num!=12){print $line."\n";}
				my $e_hu = $hgtele[3];
				my $e_ai = $hgtele[6];
				my $e_chs = $hgtele[10];
				#print $file."\t".$line."\n";
				if(defined $e_hu and defined $e_ai and defined $e_chs){
					
				if ($e_hu > $hu and $e_ai > $ai and $e_chs > $chs){
				my @line_element=split /\t+/, $line;	
				my $gene=$line_element[0];
				my $geneinfor=$gene."\t".$taxonid."\t".$fathertaxonid."\t".$time."\t".$fathertime;
				my $nodeinfor=$taxonid."\t".$fathertaxonid."\t".$time."\t".$fathertime;
	
				my $reciptor=$hashname{$taxonid};
				my $donor=$line_element[11];
				$hashnodeinfor{$reciptor}=$nodeinfor;
				foreach my $termkey (sort keys %hashterm){
					if($donor=~m/$termkey/){
						if(exists $hashtransfer{$reciptor}{$termkey}){
							$hashtransfer{$reciptor}{$termkey}++;
							}else{
							$hashtransfer{$reciptor}{$termkey}=1;
							}
					}
				}
				$hashtransfergene{$reciptor}{$donor}{$gene}=$geneinfor;
				

			}
		}elsif(! defined $e_hu){print "the undefined line is\t".$file."\n";}
}
	}
}
}
}
}
my $outgenefile="a".$ai."_h".$hu."_chs".$chs."_rc".$rankcoverage.".skip_gene";
open OFGU, ">$outgenefile";
print OFGU "RECIPTOR\tDONOR\tGENE\tRECIPTORid\tRECIPTORfatherid\tRECIPTORtime\tRECIPTORfathertime\n";
foreach my $reciptor ( sort keys %hashtransfergene){
	my $subhashtg=$hashtransfergene{$reciptor};
	foreach my $donor (sort keys %{$subhashtg}){
		my $subsubhashtg=$hashtransfergene{$reciptor}{$donor};
		foreach my $gene (sort keys %{$subsubhashtg}){
			print OFGU $reciptor."\t".$donor."\t".$hashtransfergene{$reciptor}{$donor}{$gene}."\n";
		}
	}
}

my $outfile="a".$ai."_h".$hu."_chs".$chs."_rc".$rankcoverage.".skip_stat";
open OFU, ">$outfile";

print OFU "node\ttaxonid\tfathertaxonid\ttime\tfathertime\t";
foreach my $termkey (sort keys %hashterm){
		print OFU $termkey."\t";
	}
print OFU "\n";
foreach my $reciptor ( sort keys %hashtransfer){
	print OFU $reciptor."\t".$hashnodeinfor{$reciptor}."\t";
	
	foreach my $termkey (sort keys %hashterm){
		if(exists $hashtransfer{$reciptor}{$termkey}){
			print OFU $hashtransfer{$reciptor}{$termkey}."\t";
			}else{print OFU "0\t";}
		
	}
print OFU "\n";

}
print "[info] assign time to hgt and flitered hgts according to the gene coverage within the rank levels. \n";
sub Usage {
print <<End ;
The scripts is $0;
Usage:  
     $0 -t timeline-final -r rankorder -fp dir -term terms_file -AI [num] -hu [num]
               -h/help	print the Usage;
               -t/timeline	timeline-file of each nodes/taxonomy id.
               -r/rankorder	the rank of node id from speices level to kingdom level.
               -fp/filepath	the path to hgtfiles. 
               -AI Alien Index to filter [default 10].
               -hu	hU to filter [default 30].
               -term terms file including kingdom names.
               -rc rankcoverage.
               -chs CHS.
End
}







