use strict;
use warnings;
use Getopt::Long;


my ($file,$pathn);
our (%parent,%rank,%names,%crashhash);

GetOptions(
"file|f:s"  => \$file,
'p|path:s'  => \$pathn,
);

unless (defined $pathn){&Usage;exit;}
unless (defined $file){&Usage;exit;}

## parse nodes and names;
open (my $NODES, "$pathn/nodes.dmp") or die $!;
	while(<$NODES>){
		chomp;
		next if /\#/;
		my @elements = map { s/^\s+|\s+$//gr } split (/\|/, $_); 
		my $taxid=$elements[0];
		my $parentid=$elements[1];
#while($taxid){print $taxid."\t".$parentid."\n";}
		my $rk=$elements[2];
		if($taxid != $parentid){
			$parent{$taxid}=$parentid;
			$rank{$taxid}=$rk;
	}
}
close $NODES;

open (my $NAMES, "$pathn/names.dmp") or die $!;
	while(<$NAMES>){
		chomp;
		next if /\#/;
		my @elements = map { s/^\s+|\s+$//gr } split (/\|/, $_); 
		my $taxid=$elements[0];
		 $names{$taxid} = $elements[1] if ($elements[3] eq "scientific name");
	
}
close $NAMES;

open (my $MERGED , "$pathn/merged.dmp") or die $!;
while (<$MERGED>) {
      chomp;
      next if /\#/;
      my @F = map { s/^\s+|\s+$//gr } split (/\|/, $_);
      $parent{$F[0]} = $F[1]; ## key= old taxid; value= new taxid
      ## this will behave as if old taxid is a child of the new one, which is OK I guess
    }


my $outfile=$file.".mark";
my $outlog=$outfile.".log";
open(my $FILE, "$file") if (-e $file);
open (my $OUT,">$outfile");
while(<$FILE>){
		chomp;
	my $line=$_;
	if(/\#/){print $OUT $line."\n";next ;}
	
	my @E = split (/\t/, $line); 

	if($E[1]=~m/\S+\;\S+\;\S+_([0-9]+)/){
		my $taxid=$1;
if(defined $taxid){
		my $crashrank = tax_walk($taxid);
		my $newsource=$E[1].$crashrank;
		$line=~s/$E[1]/$newsource/;
		
	}
	
	}

if($E[2]=~m/[0-9]+_[0-9]+@@([0-9]+)/){
    my $taxid=$1;
if(defined $taxid){
    my $crashrank = tax_walk($taxid);
    my $newreciptor=$E[2].$crashrank;
    $line=~s/$E[2]/$newreciptor/;
    
  }
  
  }


print $OUT $line."\n";


	
}
close $FILE;

sub tax_walk {
	my $taxidin = $_[0];
	my $taxidin_rank = $rank{$taxidin};

  $crashhash{"33630"}="_CRASH_Alveolata";
  $crashhash{"543769"}="_CRASH_Rhizaria";
  $crashhash{"33634"}="_CRASH_Stramenopiles";
  $crashhash{"3027"}="_CRASH_Cryptophyta";
  $crashhash{"2830"}="_CRASH_Haptophyta";
  $crashhash{"38254"}="_CRASH_Glaucophyta";
  $crashhash{"2157"}="_Archaea";
  $crashhash{"2"}="_Bacteria";
  $crashhash{"4751"}="_Fungi";
  $crashhash{"33208"}="_Metazoa";
  $crashhash{"2763"}="_Rhodophyta";
  $crashhash{"33090"}="_Viridiplantae";
  $crashhash{"38254"}="_Glaucophyta";
  my $crashrank="_na";

  while (1) {
  	if(exists $crashhash{$taxidin}){

    		$crashrank=$crashhash{$taxidin};
    	
    	last;
   	} elsif ($taxidin == 1) {
    	#print "last\n";
      last;
    } else {
    	#print "nonon\n";
    	if(defined $parent{$taxidin}){
    	 $taxidin = $parent{$taxidin};
      $taxidin_rank = $rank{$taxidin};	
      }else{last;}
     

    }
  }

  return $crashrank;

}



open OUTLOG,">$outlog" ;
print OUTLOG "CRASH species are markded. \n";
print "[info]	CRASH species are markded. \n";


sub Usage {
print <<End ;
The scripts is $0;
Usage:  
  perl $0 -f file -p pathn;
              -f hgt AI result file;
              -p path to taxnomy files;
End
}