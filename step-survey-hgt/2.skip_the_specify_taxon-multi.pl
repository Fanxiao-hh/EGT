#this script is to delete the hgt candidates within the same specified rank, eg. phylum; 
use strict;
use warnings;
use Getopt::Long;
use List::MoreUtils ':all';

our ($files,$pathn,$skiprank,$workpath);
our (%parent,%rank,%names);

GetOptions(
"files|f:s"  => \$files,
"skiprank|sr:s" => \$skiprank,
'p|pathn:s'  => \$pathn, 
"wp|workpath:s" => \$workpath,
);

unless (defined $pathn){&Usage;exit;}
unless (defined $files){&Usage;exit;}
unless (defined $skiprank){&Usage;exit;}
unless (defined $workpath){&Usage;exit;}

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

    print "[info] node files are successfully read by $0\n";
open (my $F,"<$files");

while(<$F>){
if(/\#/){next;}	
chomp;
my $file=$_;

my $outfile=$file.".filter";
open(my $FILE, "<$file") if -e $file;
open (my $OUT,">$outfile");
while(<$FILE>){
	chomp;
	my $line=$_;

	if(/\#/){print $OUT $line."\n";next ;}
	
	my @E = split (/\t/, $line); 
	if(@E>0){
	my $queryid=$E[0];$queryid=~s/[0-9]+_[0-9]+@@//;
	my $sourceid=$E[-1]; $sourceid=~s/.*_([0-9]+)/$1/;

	if(defined $queryid && defined $sourceid){
		
		if($skiprank=~m/phylum/){
     
		my @walkquery = tax_walk_to_phylum($queryid);
		my @walksource = tax_walk_to_phylum($sourceid);
		my %hash_a = map{$_=>1} @walkquery;
		my @common = grep {$hash_a{$_}} @walksource;
		print $OUT $line."\n" if @common==0;
		}
		
		if($skiprank=~m/class/){
		my @walkquery = tax_walk_to_class($queryid);
		my @walksource = tax_walk_to_class($sourceid);
		my %hash_a = map{$_=>1} @walkquery;
		my @common = grep {$hash_a{$_}} @walksource;
		print $OUT $line."\n" if @common==0;
		}

		if($skiprank=~m/order/){
		my @walkquery = tax_walk_to_order($queryid);
		my @walksource = tax_walk_to_order($sourceid);
		my %hash_a = map{$_=>1} @walkquery;
		my @common = grep {$hash_a{$_}} @walksource;
		print $OUT $line."\n" if @common==0;
		}

		if($skiprank=~m/family/){
		my @walkquery = tax_walk_to_family($queryid);
		my @walksource = tax_walk_to_family($sourceid);
		my %hash_a = map{$_=>1} @walkquery;
		my @common = grep {$hash_a{$_}} @walksource;
		print $OUT $line."\n" if @common==0;
		}

	}else{
		print $files."\n";
	}

}else{
	last;
}
	
}
close $FILE;
}
close $F;
my $ofile=$workpath."/".$0.".log";
open OU,">$ofile";

print OU "[info] successfully skip the hgt within the same $skiprank\n";
print "[info] successfully skip the hgt within the same $skiprank\n"; 


sub tax_walk_to_family {
	my $taxidin = $_[0];
	#print $taxidin."\n";
	my $father = $parent{$taxidin};
	
  my $father_rank = $rank{$father};
  #print $father."\t".$father_rank."\n";
  my @CRASH=qw(33630 543769 33634 3027 2830);#Alveolata,#Rhizaria,#Stramenopiles,#Cryptophyta,#Haptophyta
  my @walk_way=();
  push @walk_way,$taxidin;
  push @walk_way,$father;
  #my ($phylum,$kingdom,$superkingdom) = ("undef","undef","undef");

  while (1) {
  	#print $father."_1\n";
	if ($father_rank eq "family") {
    	#print "order\n";
      last;
    }elsif ($father_rank eq "order") {
    	#print "order\n";
      last;
    }elsif ($father_rank eq "class") {
    	#print "order\n";
      last;
    }elsif ($father_rank eq "phylum") {
    	#print "phylum\n";
      last;
    } elsif(grep /^$father$/, @CRASH){
    	#print $father."\t"."@CRASH\n";
    	last;
   	}elsif ($father_rank eq "kingdom") {
   		#print "kingdom\n";
      last;
    } elsif ($father_rank eq "superkingdom") {
    	#print "superkingdom\n";
      last;
    } elsif ($father == 1) {
    	#print "last\n";
      last;
    } else {
    	#print "nonon\n";
    	if(defined $parent{$taxidin}){
      $father = $parent{$father};
      $father_rank = $rank{$father};
      #print $father."\n";
      push @walk_way,$father;}else{last;}
    }
  }

  return @walk_way;

}


sub tax_walk_to_order {
	my $taxidin = $_[0];
	#print $taxidin."\n";
	my $father = $parent{$taxidin};
	
  my $father_rank = $rank{$father};
  #print $father."\t".$father_rank."\n";
  my @CRASH=qw(33630 543769 33634 3027 2830);#Alveolata,#Rhizaria,#Stramenopiles,#Cryptophyta,#Haptophyta
  my @walk_way=();
  push @walk_way,$taxidin;
  push @walk_way,$father;
  #my ($phylum,$kingdom,$superkingdom) = ("undef","undef","undef");

  while (1) {
  	#print $father."_1\n";
	if ($father_rank eq "order") {
    	#print "order\n";
      last;
    }elsif ($father_rank eq "class") {
    	#print "order\n";
      last;
    }elsif ($father_rank eq "phylum") {
    	#print "phylum\n";
      last;
    } elsif(grep /^$father$/, @CRASH){
    	#print $father."\t"."@CRASH\n";
    	last;
   	}elsif ($father_rank eq "kingdom") {
   		#print "kingdom\n";
      last;
    } elsif ($father_rank eq "superkingdom") {
    	#print "superkingdom\n";
      last;
    } elsif ($father == 1) {
    	#print "last\n";
      last;
    } else {
    	#print "nonon\n";
    	if(defined $parent{$taxidin}){
      $father = $parent{$father};
      $father_rank = $rank{$father};
      #print $father."\n";
      push @walk_way,$father;}else{last;}
    }
  }

  return @walk_way;

}



sub tax_walk_to_class {
	my $taxidin = $_[0];
	#print $taxidin."\n";
	my $father = $parent{$taxidin};
	
  my $father_rank = $rank{$father};
  #print $father."\t".$father_rank."\n";
  my @CRASH=qw(33630 543769 33634 3027 2830);#Alveolata,#Rhizaria,#Stramenopiles,#Cryptophyta,#Haptophyta
  my @walk_way=();
  push @walk_way,$taxidin;
  push @walk_way,$father;
  #my ($phylum,$kingdom,$superkingdom) = ("undef","undef","undef");

  while (1) {
  	#print $father."_1\n";
	if ($father_rank eq "class") {
    	#print "order\n";
      last;
    }elsif ($father_rank eq "phylum") {
    	#print "phylum\n";
      last;
    } elsif(grep /^$father$/, @CRASH){
    	#print $father."\t"."@CRASH\n";
    	last;
   	}elsif ($father_rank eq "kingdom") {
   		#print "kingdom\n";
      last;
    } elsif ($father_rank eq "superkingdom") {
    	#print "superkingdom\n";
      last;
    } elsif ($father == 1) {
    	#print "last\n";
      last;
    } else {
    	#print "nonon\n";
    	if(defined $parent{$taxidin}){
      $father = $parent{$father};
      $father_rank = $rank{$father};
      #print $father."\n";
      push @walk_way,$father;}else{last;}
    }
  }

  return @walk_way;

}


sub tax_walk_to_phylum {
	my $taxidin = $_[0];
	#print $taxidin."\n";
	my $father = $parent{$taxidin};
	
  my $father_rank = $rank{$father};
  #print $father."\t".$father_rank."\n";
  my @CRASH=qw(33630 543769 33634 3027 2830);#Alveolata,#Rhizaria,#Stramenopiles,#Cryptophyta,#Haptophyta
  my @walk_way=();
  push @walk_way,$taxidin;
  push @walk_way,$father;
  #my ($phylum,$kingdom,$superkingdom) = ("undef","undef","undef");

  while (1) {
  	#print $father."_1\n";
    if ($father_rank eq "phylum") {
    	#print "phylum\n";
      last;
    } elsif(grep /^$father$/, @CRASH){
    	#print $father."\t"."@CRASH\n";
    	last;
   	}elsif ($father_rank eq "kingdom") {
   		#print "kingdom\n";
      last;
    } elsif ($father_rank eq "superkingdom") {
    	#print "superkingdom\n";
      last;
    } elsif ($father == 1) {
    	#print "last\n";
      last;
    } else {
    	#print "nonon\n";
    	if(defined $parent{$taxidin}){
      $father = $parent{$father};
      $father_rank = $rank{$father};
      #print $father."\n";
      push @walk_way,$father;}else{last;}
    }
  }

  return @walk_way;

}

sub Usage {
print <<End ;
The scripts is $0;
Usage:   
  perl $0 -f files -p pathn;
              -f hgt candidates files;
              -p path to taxnomy files;
              -sk specify a range of rank within which the HGT should be deleted: phylum/class/order/family;
              -wp work path;
End
}