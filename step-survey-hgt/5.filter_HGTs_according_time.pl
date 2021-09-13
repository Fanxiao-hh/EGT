
use strict;
use warnings;
use Getopt::Long;
use List::MoreUtils ':all';

our $valuefile;
our $AIresult;
our $orthofile;
our $speciesorder;
our $term;


our %hash;
our %hashs;
our %hashterm;
our %hashvalue;
our %hashduiyingname;
our %hashtransgenes;
our %hashtransogs;

if($ARGV[0]=~m/-h/) {&Usage;exit;}

GetOptions(
"AIresult|A:s"     => \$AIresult,
"valuefile|v:s" => \$valuefile,
"orthofile|og:s" => \$orthofile,
"speciesorder|so:s" => \$speciesorder,
"term|t:s" => \$term,

);

unless (defined $AIresult){&Usage;exit;}
unless (defined $valuefile){&Usage;exit;}
unless (defined $orthofile){&Usage;exit;}
unless (defined $speciesorder){&Usage;exit;}
unless (defined $term){&Usage;exit;}

open INS, "<$speciesorder" or die $!;#read the specie order;
	while(<INS>){
	chomp;
	if(/(\S+)\s+(\S+)\s+(\S+)/){
		my $species_id=$1;
		my $duiyingname=$3;
		$hashduiyingname{$species_id}=$duiyingname;
		}
	}
print "[info_$0] finish to read the species order\n";

open TER," <$term" or print "can not find termfile\n";
	while(<TER>){
				chomp;
				my $line=$_;
				$hashterm{$line}=1;
			}
print "[info_$0] finish to read the terms file\n";

open VF,"<$valuefile" or die $!;
while(<VF>){
	chomp;
	my $value=$_;
	$hashvalue{$value}=1;
}

open IN1,"<$orthofile" or die $!; #this loop read orthlogous files;
	while(<IN1>){
	chomp;	
	if(/^([0-9]+)\s+(.*)\s+\$/){
		my $OG=$1;
		my $line=$2;
		my @ele =split / /, $line;
		foreach my $elements(@ele){
			$hash{$elements}=$OG;
			foreach my $value (keys %hashvalue){$hashs{$value}{$OG}=$line;}	
}
}
}
close IN1;

print "[info_$0] successfully read orthologfiles\n";

			foreach my $value (sort keys %hashvalue){
				print "[info_$0] filter gene within $value million years.\n";
open IN2, "<$AIresult" or die $!; #this loop read ai-result and filter the genes according to the time;

	while(<IN2>){
		next if /\#/;
		chomp;
		my $line =$_;
		my @elements=split /\t/,$line;
		my $donor=$elements[1];
		my $gene=$elements[2];
		$gene=~s/\@\@\S+//;
		my $species_id=$gene;
		$species_id=~s/_[0-9]+//;
		my $genetaxid=$hashduiyingname{$species_id};
	
		my $reciptortime=$elements[5];
		my $reciptorfathertime=$elements[6];


		if(!exists $hash{$gene}){next;}
			my $og=$hash{$gene};

				if($reciptortime=~/NA/){

					if($reciptorfathertime<$value){
	
						my $newgenetaxid=();
						my $increasenum=0;
						if($genetaxid=~m/(\S+)--(\S+)--(\S+)/){
							$newgenetaxid=$1;
						}elsif($genetaxid=~m/(Chromalveolate)—(\S+)--(\S+)/){
							$newgenetaxid=$1;
						}

						foreach my $termkey (sort keys %hashterm){
							if($donor=~m/$termkey/){
								$hashtransgenes{$value}{$newgenetaxid}{$termkey}{$gene}=1;
								$hashtransogs{$value}{$newgenetaxid}{$termkey}{$og}=1;
								$increasenum++;
							}

						}
						if($increasenum==0){
							$hashtransgenes{$value}{$newgenetaxid}{"other-Eukaryotes"}{$gene}=1;
							$hashtransogs{$value}{$newgenetaxid}{"other-Eukaryotes"}{$og}=1;
						}
						$hashs{$value}{$og}=~s/ $gene / /;	
					}
				}else{
					if($reciptortime<$value){
						my $newgenetaxid=();
						my $increasenum=0;
						if($genetaxid=~m/(\S+)--(\S+)--(\S+)/){
	 						$newgenetaxid=$1;
						}elsif($genetaxid=~m/(Chromalveolate)—(\S+)--(\S+)/){
							$newgenetaxid=$1;
						}

						foreach my $termkey (sort keys %hashterm){
							if($donor=~m/$termkey/){
								$hashtransgenes{$value}{$newgenetaxid}{$termkey}{$gene}=1;
								$hashtransogs{$value}{$newgenetaxid}{$termkey}{$og}=1;
								$increasenum++;
							}

						}
		
						if($increasenum==0){
						$hashtransgenes{$value}{$newgenetaxid}{"other-Eukaryotes"}{$gene}=1;
						$hashtransogs{$value}{$newgenetaxid}{"other-Eukaryotes"}{$og}=1;
						}
	
		
						$hashs{$value}{$og}=~s/ $gene / /;	
	
					}
				}

			


	}
close IN2;
}
print "[info_$0] successfully filter the gene according to thier time\n";

foreach my $value (sort keys %hashvalue){
	my $outfile=$AIresult."-".$value."-ortho"; #this loop make a new orthologous file;
	open OUT, ">$outfile" ; 
	foreach my $key (sort keys %{$hashs{$value}}){
		print OUT $key."\t";
		print OUT $hashs{$value}{$key};
    	print OUT "\$\n";
		my $left=$hashs{$value}{$key};
		$left=~s/\$//;
		my @remain=split /\s/,$left;
 		foreach my $species_id(@remain){
 	
			my $geneid=$species_id;

 			$species_id=~s/([0-9]+)_[0-9]+/$1/;
			my $genetaxid=$hashduiyingname{$species_id};

#print $species_id."@@@@\n";
#print $genetaxid."\ttax\n";
#print $geneid."\n";

			my $termkey="unHGT";

			if($genetaxid=~m/(\S+)--(\S+)--(\S+)/){
			my $newgenetaxid=$1;
			$hashtransgenes{$value}{$newgenetaxid}{$termkey}{$geneid}=1;
			$hashtransogs{$value}{$newgenetaxid}{$termkey}{$key}=1;
	#print $newgenetaxid."\t".$termkey."\t".$geneid."\n";
			}elsif($genetaxid=~m/(Chromalveolate)—(\S+)--(\S+)/){
				my $newgenetaxid=$1;
				$hashtransgenes{$value}{$newgenetaxid}{$termkey}{$geneid}=1;
				$hashtransogs{$value}{$newgenetaxid}{$termkey}{$key}=1;
			}
 		}	
	}

my $outstatfile=$AIresult."-".$value."-ortho.stat"; #this loop make a summary file;
my $outstatfilew=$AIresult."-".$value."-ortho.statw";
	open OUT2, ">$outstatfile" ;
	open OUT3, ">$outstatfilew" ;
	print OUT2 "kingdom\tDonor\tnum\n" ;
	print OUT3 "kingdom\tDonor\tnum\n" ;
	foreach my $genetaxid (sort keys %{$hashtransgenes{$value}}){
		foreach my $donor (sort keys %{$hashtransgenes{$value}{$genetaxid}}){
			my @genes=keys %{$hashtransgenes{$value}{$genetaxid}{$donor}};
			my @ogs=keys %{$hashtransogs{$value}{$genetaxid}{$donor}};
			my $genecount=@genes;
			my $ogcount=@ogs;

			
			if($genecount>0){
				print OUT2 $genetaxid."\t".$donor."\t".$genecount."\n";
				print OUT3 $genetaxid."\t".$donor."\t".$ogcount."\n";
			}else{
				print OUT2 $genetaxid."\t".$donor."\t0\n";
				print OUT3 $genetaxid."\t".$donor."\t0\n";	
			}
			
		}	
	}
	print "[info] finish export $outstatfile \n";
	my $log=$AIresult."_".$value.".log";
	open LOG,"> $log";
	print LOG "[info] finish export $log;finish 0.PCA-filter.pl. \n";
	print "[info] finish export $log;finish 0.PCA-filter.pl. \n";
}

my $log0=$AIresult."_timelist.log";
	open LOG0,"> $log0";
	print LOG0 "[info] finish filter hgt accoding time. \n";
	print "[info] finish filter hgt accoding time. \n";


sub Usage {
print <<End ;
The scripts is $0;
Usage:  
  perl $0 -A a30_h90gene -v 100 -og I1.2-out_orth-new-2-ok
               -h/help	print the Usage;
               -A/AIresult	a30h90 result;
               -v/valuefile 	time value below which you dicide to drop;
               -og 	orthologous file;
               -so	speciesorder file;
               -t/terms;
End
}

