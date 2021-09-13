use warnings;
use strict;
use Getopt::Long;

our %hash;
our $path;

GetOptions(
"path|p:s"     => \$path,
);
unless (defined $path){&Usage;exit;}
#this script is to check if the file is completed; 
opendir (DIR, $path)|| die $!;
while(my $file=readdir(DIR)){
	next if ($file =~m/(^\.$)/);
	next if ($file =~m/(^\.\.$)/);
open(my $FILE, "<$path/$file");
while(<$FILE>){
	chomp;
	my $line=$_;
	if(/\#/){next ;}
	
	my @E = split (/\t/, $line); 
	if($E[-1]!~m/\S+\;\S+\;\S+_[0-9]+/){
$hash{$file}=$E[-1];
		
	}
}
} 
 my $log="check_completness.log";
 my $err="check_completness.err";
 my $length=%hash;
 if($length==0){
 	open OUT,">$log";
 	print OUT "all hgt candidate files are completed.\n";
 	print "[info] HGT candidates files are checked out .\n";
 }else{
 	open OUTE,">$err";
 	foreach my $key (keys %hash){
 		print OUTE $key."\n";
 	}
 }
 sub Usage {
print <<End ;
The scripts is $0;
Usage:  
     perl $0 -p path;

End
}

