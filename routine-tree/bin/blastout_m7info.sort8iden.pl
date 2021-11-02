#! /usr/bin/perl
# take treated m7 blast output file (*.info)
# sort hits according to sequence identity among those with HSP >= $len_cut; 

$len_cut = 200;

$in = shift;
open  (IN, $in) or die "cof $in by blastout_m7info.sort8iden.pl\n";
open  (OUT,">$in.sort_iden.txt");


%all = ();
while (<IN>) { 
      my $line = $_;
      my @a = split /\t/;
      my $iden = $a[6];
      my $hspl = $a[5];
      $all{$.}=$line;
      $idn{"$line"}=$iden if $hspl>=$len_cut;
}
#print scalar keys %idn, " hit with HSP >=$len_cut aa\n";

my @sorted = sort { $idn{$b} <=> $idn{$a} } keys %idn;

foreach ( sort {$a <=> $b} keys %all ) {  #
      my $nb   = $_;
      my $line = $all{$_};
      my @a = split (/\t/, $line);
      my $iden = $a[5];
      if ($iden>=$len_cut) {
          my $put = shift (@sorted);
          $all{$nb} = $put;
          print OUT $put;
      } else {
          print OUT $line;
      }
}


