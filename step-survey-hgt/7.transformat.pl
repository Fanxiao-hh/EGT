
use Getopt::Long;

if($ARGV[0]=~m/-h/) {&Usage;exit;}

our $sequenceID;
our $speciesorder;
our $orthofile;
our $min=5;
GetOptions(
"sequenceID|ID:s" => \$sequenceID,
"speciesorder|so:s" => \$speciesorder,
"orthofile|og:s" => \$orthofile,
"minimalcount|min:s" => \$min,
);

unless (defined $sequenceID){&Usage;exit;}
unless (defined $speciesorder){&Usage;exit;}
unless (defined $orthofile){&Usage;exit;}


open IN,"<$sequenceID" or die $!; ##read the sequenceID;
while(<IN>){
chomp;
if(/(\S+)\s+(\S+)\s+(\S+)/){
my $number=$1; 
my $orthid=$2;$hashorth{$number}=$orthid;
my $seqid=$3;$hashseq{$number}=$seqid;

}else{
print $_;
}

}

close IN;
print "[INFO] finish to read the sequenceID\n";

open IN1, "<$speciesorder" or die $!;#read the specie order;
$order_species==0;
while(<IN1>){
$order_species++;
chomp;
if(/(\S+)\s+(\S+)\s+(\S+)/){
$species_id=$1;
$duiyingfile=$2;
$duiyingname=$3;
$hashspeciesorder{$order_species}=$species_id;
$hashduiyingname{$species_id}=$duiyingname;
}
}
print "[INFO] finish to read the species order\n";
open IN2,"<$orthofile" or die $!;##read the Ixx-out file ;
my $orth=$orthofile."_orth";
my $seq=$orthofile."_seq";
my $gene_count=$orthofile."_Gene_count";

open OUT1,">$orth";
open OUT2,">$seq";
open OUT3,">$gene_count";
while(<IN2>){
if(/^(OG.[0-9]+)\s+(.*)\s+\$/){
my $OG=$1;
my $group=$2;
@{$OG}=split / /,$group;
my $lengthOG=@{$OG};
$OGlenth{$OG}=$lengthOG;
if($lengthOG>=$min){


print OUT1 $OG."\t";

print OUT2 $OG."\t";
foreach my $key(sort {$a<=>$b} @{$OG}){
print OUT1 $hashorth{$key}." ";
print OUT2 $hashseq{$key}." ";
if($hashorth{$key}=~m/([0-9]+)_([0-9]+)/){
$species_id2=$1;
$count{$OG}{$species_id2}{$key}=1;
}
}
print OUT1 "\n";
print OUT2 "\n";
}
}
}
 print "[INFO] complete export the orth and seq\n";	
print OUT3 "id\t";
foreach my $nkey (sort {$a<=>$b} keys %hashspeciesorder){
print OUT3 $hashduiyingname{$hashspeciesorder{$nkey}}."\t";
}
print OUT3 "Total.\n";
foreach my $ckey( sort {$count{$a}<=>$count{$b}} keys %count){
if($OGlenth{$ckey} >= $min){


print OUT3 $ckey."\t";
#print $ckey."\n";
foreach my $nkey (sort {$a<=>$b} keys %hashspeciesorder){
if(exists $count{$ckey}{$hashspeciesorder{$nkey}}){
 $length= keys(%{$count{$ckey}{$hashspeciesorder{$nkey}}});
}else{
 $length=0;
}
print OUT3 $length."\t";
}

print OUT3 $OGlenth{$ckey}."\n";
}
}

print "[INFO] finish to export count\n";

my $log2=$orthofile.".log2";
open LOG,"> $log2";
print LOG "finish 2.transformat.pl. \n";



sub Usage {
print <<End ;
The scripts is $0;
Usage:  
  perl $0 -ID sequenceID -og a30h90 -so specieorder
               -h/help	print the Usage;
               -ID/sequenceID	;
               -og 	a30h90 orthologous file;
               -so	specieorder file.
               -min [num].
End
}

