#open IN,"<$ARGV[0]" or die $!;
use Getopt::Long; 
 my ($outpath,$path);
GetOptions(
'o|outpath:s' => \$outpath,
'p|path:s' => \$path,
  );
unless (defined $outpath){&Usage;exit;}
unless (defined $path){&Usage;exit;}
	my $species=0;
	while($species < 1157){
		my $level=1;
		while($level < 37){
		my $file =$path."/"."level-".$level."-taxonaid-Species".$species.".HGT.txt.filter";
		if(-e $file){
# print $species."\t".$level."\n";		
			open IN,"< $file" or die $file."is not exists\n";
			my $outfile=$outpath."/"."level-".$level."-taxonaid-Species".$species.".HGT.txt.filter.filterout";
			open OUT, "> $outfile";
			while(<IN>){
			if(/(\S+)\s+.*/){
			#print $species."\t".$level."\n";
			my $id=$1;
			my $line=$_;
			@elements=split /\t/,$line;
			my $bitout=$elements[4];
			$outvalue{$species}{$level}{$id}=$bitout;
			$hash{$species}{$level}{$id}=$line;
			my $lastlevel=$level-1;
if($level eq 1){
print OUT $line;
#print $id."\t".$outvalue{$species}{$level}{$id}."\t".$level."\n";
}else{ 
	my $changlevel=$level;
	my $count=0;
while($changlevel >1){
	$changlevel--;
	if(exists $outvalue{$species}{$changlevel}{$id} ){
		if( $outvalue{$species}{$changlevel}{$id} == $outvalue{$species}{$level}{$id} ){
$count++;
#print $id."\t".$outvalue{$species}{$lastlevel}{$id}."\t".$outvalue{$species}{$level}{$id}."\t".$level."\n";		
	
}
}
}
if($count eq 0){
	#print $id."\t".$outvalue{$species}{$lastlevel}{$id}."\t".$outvalue{$species}{$level}{$id}."\t".$level."\n";
print OUT $line;
}
}
}
}
}
		$level++;
}
	$species++;
}

my $outfile=$0."log";
open OUT, ">$outfile";
print OUT "hgt files is filtered via removing the repeat events\n";
print "[info] hgt files is filtered via removing the repeat events\n";

sub Usage {
print <<End ;
The scripts is $0;
Usage:  
     $0 -p path -o outpath
               -h/help	print the Usage;
               -p HGT.txt files path;
               -o fliter output path;
End
}

