
use Getopt::Long;
 
our $sequenceID;
our $orthofile;
our %hashorth;
our %hashseq;
our %hashnumber;


GetOptions(
"sequenceID|ID:s" => \$sequenceID,
"orthofile|og:s" => \$orthofile,
);

unless (defined $sequenceID){&Usage;exit;}
unless (defined $orthofile){&Usage;exit;}


open IN,"<$sequenceID" or die $!; ##read the sequenceID;
while(<IN>){
chomp;
if(/(\S+)\s+(\S+)\s+(\S+)/){
my $number=$1; 
my $orthid=$2;
my $seqid=$3;
#$hashorth{$number}=$orthid;
#$hashseq{$number}=$seqid;
$hashnumber{$orthid}=$number;
#print $orthid."\t".$number."\n";
}else{
print $_;
}

}

close IN;
print "[INFO] finish to read the Sequenceid file\n";



open IN2,"<$orthofile" or die $!; ##I1.2-out
my $orthofile_number=$orthofile."_number";
open OUT4,">$orthofile_number";
while(<IN2>){
if(/^([0-9]+)\t(.*)(\s+)?\$/){
my $OG = "OG.".$1;
my $group = $2;

 @{$OG} = split / /,$group;
print OUT4 $OG."\t";
foreach my $key (sort {$a<=>$b} @{$OG}){

print OUT4 $hashnumber{$key}." ";

}
print OUT4 "\$\n";
}
	}

my $log=$orthofile.".log";
open LOG,"> $log";
print LOG "finish 6.replacegeneid.pl. \n";
print "[INFO]	finish replacing genesid.\n";

sub Usage {
print <<End ;
The scripts is $0;
Usage:  
  perl $0 -ID sequenceID -og a30h90
               -h/help	print the Usage;
               -ID/sequenceID	;
               -og 	a30h90 orthologous file.
End
}

