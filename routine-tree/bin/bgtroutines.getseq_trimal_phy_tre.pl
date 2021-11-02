#! /usr/bin/perl
# given a set of aligned sequences in fasta format, work it through trimal
# remove short sequence (ie, <40 AA) from alignment


my @a = @ARGV;
my $tmp_file = "routine.fst2phylip.tmp_log.txt";

my $db=shift(@a);

my $count=1;
foreach my $file (@a) {
        print "#-> start $count file\n";
        &FST2PHYLIP($file);

        $count++;
}


sub FST2PHYLIP {
          $in = shift   ;
          open (IN, $in) or die "cant open file $in\n";
print "#-> $in\n";
print "#        muscle3.8.31\n";
          my $query="";
          $in =~/(^\S+)_2refseq/;
          $query = $1;
#          system ("mafft   --auto $in >$in.mft0 2>$tmp_file");
          system ("seqfetch.def.pl $in $db");
          system ("cat $query >>$in.fa.txt");
          system ("perl -i -pe \'s\/\-\-\\S\+\-\-\/\-\/\' $in.fa.txt");
	  system ("muscle3.8.31 -in $in.fa.txt -out $in.mft0 2>$tmp_file");
          system ("u.move_column.pl $in.mft0 0.2 $in.mft");          

print "#        trimAL\n";
          system ("trimal.longnames.pl $in.mft");
          #system ("cp $in.mft $in.mft\-gb.txt");
          
          open (IN2, "$in.mft\-gb.txt") or die "can open file $in.mft\-gb.txt\n";
          open (OT2, ">$in.mft\-gb2.txt");          
          while (<IN2>) {
                s/\|\S+//g; 
                s/\.\_/\_/g;
                s/\_str|\_subsp//g;
                if ($.!=1 and /^(\>\S+)/) {
                     print OT2 "\n$1\n"; 
                } elsif (/^(\>\S+)/g) {
                     print OT2 "$1\n";                 
                } elsif ($.!=1) {
                        chomp($_);
                     print OT2 "$_";                         
                }
          } 
          print OT2 "\n";          
          close (OT2);
          close (IN2);          

          
print "#        Polish\n";
          open (IN3,  "$in.mft\-gb2.txt") or die "can open file $in.mft\-gb2.txt\n";
          open (OT3, ">$in.mft\-gb2.txt8");          
          while (<IN3>) { #next unless /16578/;
                if (/^\>/g) {
                    print OT3 "$_";
                } else {
                   my $line = $_;
                   $line = &Remove_orphan (2,  6, 2, $line);
                   
                   $line = &Remove_orphan (3, 10, 0, $line);
                   $line = &Remove_orphan (0, 10, 3, $line);
                   
                   $line = &Remove_orphan (1, 10, 3, $line);
                   $line = &Remove_orphan (3, 10, 1, $line);
                   
                   $line = &Remove_orphan (3, 15, 3, $line);
                   $line = &Remove_orphan (5, 20, 5, $line);
                   
                   print OT3 "$line";
                }
          } 
          close (OT3);
          close (IN3);
          #print OT3 "\n";
                                       
          system ("fst.filterbylen.pl $in.mft\-gb2.txt8 0.1 $in.mft\-gb3.txt");
          
          #system ("tcoffee.rmsite_cutoff5.sh $in.mft\-gb3.txt");
          system ("cp $in.mft\-gb3.txt $in.mft\-gb3.txt.tcoffee");
             
          system ("u.move_column.pl $in.mft\-gb3.txt.tcoffee 0.7 $in.mft\-gb3.txt.tcoffee2");
          system ("fst.filterbylen.absolute.pl $in.mft\-gb3.txt.tcoffee2 50 $in.mft\-gb3.txt.tcoffee3");
          system ("fasta2relaxedPhylip.pl $in.mft\-gb3.txt.tcoffee3 >$tmp_file 2>&1");
          system ("phy.rmorphan.keep_order.pl $in.mft\-gb3.txt.tcoffee3.phylip");
          system ("cp $in.mft\-gb3.txt.tcoffee3.phylip.sortbydef.txt $in.phy.txt");
          system ("rm $in.mft $in.mft\-gb*txt* $tmp_file $in.mft0");
          #system ("iqtree -s $in.phy.txt -m MFP -mset WAG,LG,JTT -bb 1150");
	  system ("FastTree -quiet -nopr -lg -spr 4 -mlacc 2 -slownni $in.phy.txt >$in.phy.tre 2>$tmp_file");
#print "#-> done\n";
          system ("rm $in.mft?gb**");
}



sub Remove_orphan {
                my  $gapsize1    = shift;
                my  $aminoacids  = shift; die if $aminoacids < 1;
                my  $gapsize2    = shift;
                my  $line        = shift;
                
                my $amin       = 'zzzzz';
                my $segm       = 'zzzzz';
                
                foreach my $i (1..$aminoacids) {
                    $aminoacids  = $i;
                      $gap1 = '-'  x $gapsize1;
                      $amin = '\w' x $aminoacids;
                      $gap2 = '-'  x $gapsize2;
                      $pattern = "$gap1$amin$gap2";
                      $pattern = '^'.$pattern if $gapsize1 <2;    #pattern occurs at end of seqs
                      $pattern = $pattern.'$' if $gapsize2 <2;
                      $segm = '-'  x ($gapsize1 + $aminoacids + $gapsize2);
                      $line =~ s/$pattern/$segm/g;
                      #$pattern =~ s/\\//g;
                      #print "$pattern\n$segm\n\n";
                }
                return $line;
}
