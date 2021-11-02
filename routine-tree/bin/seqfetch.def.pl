#! /usr/bin/perl 
# 
# no limit to the size of database
# 

#check command arguments
$argm=join ("\*", @ARGV);
$cmarg=$#ARGV+1;
if ($argm =~ /\-help/i||$cmarg<2) {
   die "usage: seqfetch.def.pl 1-ListOfDefinition 2-Database\n";
} 


$def  = shift or die "input is wrong \?";
$dbs  = shift or die "input is wrong \?";
$mode = shift || "head"; 

#readin sequences definitions and database from which sequences will be hooked
open (DEF, "$def") or die "cant open $def\n";
open (DBS, "$dbs") or die "cant open $dbs\n";
open (OUT, ">$def.fa.txt") or die "cant open $def.fa\n";

#put all sequences definitions in a HASH
$indef=0;
%definitions=();

if ($mode eq "all") {

   while (<DEF>) {
     chomp ($_);
     my @a = split /\s+/;    
     foreach (@a) {   
       next if /^\s*$/;
       $definitions{$_}=1;
        ++$indef;
     }
   }

} else {
   while (<DEF>) {
     if (/^\s*(\S+)/) {
        ++$indef;
        $definitions{$1}=1;
     }
  }
}
#go through the data file and check if each sequence is in the definition list
#if so, export it to outfile
my $seq=$def="";
my $seqout=0;
while (<DBS>) {
  last unless %definitions;
  if (/>(\S+)/) {
     if ($definitions{$def}) {
         ++$seqout;
         $seq=&Multiline($seq);       
         print OUT ">$def\n$seq";
         delete $definitions{$def};
     }
     $def=$1;  
     $seq="";
  } else {
     chomp ($_);
     s/\-//g;
     $seq .= $_;
  }
}
     if ($definitions{$def}) {
         ++$seqout;
         $seq=&Multiline($seq);
         print OUT ">$def\n$seq";
         delete $definitions{$def};
     }

#summary
print "$seqout\/$indef sequences have been retrieved\n";
 #list all definitions whoes sequences are not present in input database
 if (%definitions) {
    foreach (keys %definitions) {
      print "-> $_\n";
    }
    print "-> failed to be retrieved form database $db\n";
 }



#=========================================================================
#break down long sequencs into series of 100-bp fragments
sub Multiline {
    my $break=100;
    my $seq = shift;
    my $len = length ($seq);
    my $lenidx = $len-1;
    my $i;
    my $seqout=();
       for($i=0; $i<=($len-$break); $i+=$break){
           $seqout .= substr ($seq, $i, $break);
           $seqout .= "\n";
    }
    if ($i != $lenidx-$break) {
       my $lenleft= $lenidx-$i+1;
       $seqout .= substr ($seq, $i, $lenleft);
       $seqout .= "\n";
    }
    return $seqout;
}
