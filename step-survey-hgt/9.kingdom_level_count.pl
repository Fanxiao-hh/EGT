use Getopt::Long;

if($ARGV[0]=~m/-h/) {&Usage;exit;}

our %hash;
our %hashtitle;
our %newhashtitle;
our %num;
our $orthofile;

GetOptions(
"orthofile|og:s" => \$orthofile,
);
unless (defined $orthofile){&Usage;exit;}
open IN, "<$orthofile" or die $!;
while(<IN>){
chomp;
my $line=$_;
$line=~s/\s+$//;
if($line=~m/^id\s+.*/){
my @title=split /\s+/,$line;
our $titlength=@title;
$titlength=$titlength-1;
#our $hashtitle{"titlength"}=$titlength;
my $n=0;
while($n < $titlength){
	$n++;
	$hashtitle{$n}=$title[$n];
	
}
}elsif($line=~m/^(OG.[0-9]+)\s+.*/){
my $og=$1;
my @count=split /\s+/,$line;
my $oglength=@count;
$oglength=$oglength-1;

if($oglength eq $titlength){
my $n=0;
while($n < $oglength){
	$n++;
	#$hash{$og}{$hashtitle{$n}}=$count[$n];
	if($hashtitle{$n}=~m/(\S+)--(\S+)--(\S+)/){
 		my $newtitle=$1;
 		$newhashtitle{$newtitle}=1;
 		if(exists $hash{$og}{$newtitle}){
 			$hash{$og}{$newtitle}+=$count[$n];
 			$num{$og}{$newtitle}+=1;
 			}else{
 				$hash{$og}{$newtitle}=$count[$n];
 				$num{$og}{$newtitle}=1;}
 	}elsif($hashtitle{$n}=~m/(Chromalveolate)—(\S+)--(\S+)/){
		my $newtitle=$1;
		$newhashtitle{$newtitle}=1;
 		if(exists $hash{$og}{$newtitle}){
 			$hash{$og}{$newtitle}+=$count[$n];
 			$num{$og}{$newtitle}+=1;
 		}else{$hash{$og}{$newtitle}=$count[$n];$num{$og}{$newtitle}=1;}
 	}elsif($hashtitle{$n}=~m/(Total.)/){
		my $newtitle="Z-".$1;
		$newhashtitle{$newtitle}=1;
 		if(exists $hash{$og}{$newtitle}){
 			$hash{$og}{$newtitle}+=$count[$n];$num{$og}{$newtitle}+=1;
 		}else{$hash{$og}{$newtitle}=$count[$n];$num{$og}{$newtitle}=1;}

 		}else{
 		print $n."\n";
 	}
}
	
}else{
	#print $og."\n";
	next;
}
}
}

my $outfile=$orthofile."_kingdom";
open OUT ,">$outfile" or die $!;

print OUT "id\t";
foreach my $key (sort keys %newhashtitle){
	print OUT $key."\t";
}
print OUT "\n";

foreach my $key (sort keys %hash){
	print OUT $key."\t";
	$subhash=$hash{$key};
	foreach my $subkey (sort keys %{$subhash}){

		my $average=();
		my $value=();
		my $num=();
		if(exists ${$subhash}{$subkey}){
 $value=${$subhash}{$subkey};

if($num{$key}{$subkey} >0 ){$num=$num{$key}{$subkey};
 $average=$value/$num;
}else{ $average="NA";}

			print OUT $average."\t";}else{print OUT "NA\t";}
		
	}
	print OUT "\n";
}


my $log4=$orthofile.".log4";

open LOG,"> $log4";
print LOG "finish 4.kingdom_level_count.pl \n";


sub Usage {
print <<End ;
The scripts is $0;
Usage:  
  perl $0 -ID sequenceID -og a30h90 
               -h/help	print the Usage;
               -og 	 speceis_level_count;  
End
}

