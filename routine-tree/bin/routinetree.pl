#! /usr/bin/perl
#  given a sequence(s), do blast search against $db and retrieve sequences
#  compared to routin* :
#  grab taxa from query sequence, eg, Plantae.Rhodophyta--Bangiophyceae--Porphyra_umbilicalis.escontig00077
#                                     ^                                         ^
#  hit beloning to the taxa will not be listed and retrieved
#  parallel evolution with blastout_m8.phylaClass2.pl
##
## cmp v3; handling multi-seq query file (using 'blastout_m8.phylaClass3.pl'). hit-quota is splitted among query sequences
## cmp v3; 1) generate output using sorted and unsorted (by identity) gene lists 
##         2) use new version of blast-output based sequence identity (blastout.identity_8alignments.pl)
## cmp v4 classic;
##         1) allow database to be specified
##         2) allow number of processor to be specificed
##
## cmp v4  1) allow up to 10 sequences for each class and each phylum
## cmp v5  1) combine identity- and bitscore-based sorting results
## cmp v6  1) allow more controlled sequence selection. The numbers of sequences for each class, phylum
##            and self-species can be specified as input arguments
##

 
use warnings;
use strict;
use Getopt::Long;
#######################definitions;
our $db = "/Volumes/EGT2/routine-tree/bin/database/simplified_db";
our $th =  3;
our $tt = 64;
our $py = 10;
our $cs = 10;
our $sf = 6;
our $step="auto";
our $path ;
our $list;

###########################obtain the parameters;

if($ARGV[0]=~m/-h/) {&Usage;exit;}


GetOptions(

  "database|db:s" => \$db,
  "threads|th:s" => \$th,
  "tt:s" => \$tt,
  'py:s'  => \$py,
  "cs:s" => \$cs,
  "sf:s" => \$sf,
  "pt:s" => \$path,
  "fl:s" => \$list,
  "step:s" => \$step,
  );

unless (defined $db){&Usage;exit;}
unless (defined $th){&Usage;exit;}
unless (defined $tt){&Usage;exit;}
unless (defined $py){&Usage;exit;}
unless (defined $cs){&Usage;exit;}
unless (defined $sf){&Usage;exit;}
#unless (defined $files){&Usage;exit;}

print "database to be searched is $db\n";
print "$th processors are to be used for blastp search\n";
print "$tt in total, $py each phylum, $cs each class, $sf self-species\n";
print "everything looks fine? (return key for cancel, other keys to contiue)\n";

##################################################

my $c = 0;
my @a = ();
if(defined $list){
open INL,"<$list" or die "cant open $list\n";
while(<INL>){
	chomp;
	my $seqname=$_;
	push (@a,$seqname);
}
}else{
	 @a = @ARGV;
}

##################   ###################      #############
#					main progress				          #
###########    ##################     #####################

our $tmp_file = $path."/routine.fst2phylip.tmp_log.txt";
our $db2=$path."/routine.qry2refseq2hitfasta.temp.fa.txt";

print "$0 -db $db -th $th -tt $tt -py $py -cs $cs -sf $sf -pt $path $a\n";
if($step eq "auto"){
	#if (-e "$file\_2refseq.info") {} else {}
}elsif($step eq "123" ){
	my $cc = &STEPONE(@a);
	my $count = &STEPTWO(@a);
	print "-> step-2 finish seq trimal for $count file\n";
	my $countbuild = &STEPTHREE(@a);
	print "-> step-3 finish tree build for $countbuild file\n";	
}elsif($step eq "23"){
	my $count = &STEPTWO(@a);
	print "--> step-2 finish seq trimal for $count file\n";
	my $countbuild = &STEPTHREE(@a);
	print "-> step-3 finish tree build for $countbuild file\n";	
}elsif($step eq "3"){
	my $countbuild = &STEPTHREE(@a);
	print "-> step-3 finish tree build for $countbuild file\n";
}elsif($step eq "1"){
	my $cc = &STEPONE(@a);
}elsif($step eq "12"){
	my $cc = &STEPONE(@a);
	my $count = &STEPTWO(@a);
	print "-> step-2 finish seq trimal for $count file\n";
}elsif($step eq "2"){
	my $count = &STEPTWO(@a);
	print "-> step-2 finish seq trimal for $count file\n";
}






sub STEPONE{##################################################step-1 blastp search
print "-> step-1 blastp search\n";
my @a =@_;
foreach my $file (@a) {
        $c++;
        open (IN,  $path."/".$file) or die "cant open $file\n";
        my $taxa_remove = "aaa";  
        while (<IN>) {
              if (/>\S+-([a-z]+\_[a-z]+)/gi) { 
                   $taxa_remove = $1; 
              }  elsif (/>\S+-([a-z]+\.)/gi) { 
                   $taxa_remove = $1; 
              }
              print "######query species is $taxa_remove\n";
              last;
        } close (IN);

        print "####-> $c blast mission $file\n";
        system ("rm -rf $path/$file\_2refseq $path/$file\_2refseq.info $path/$file\_2refseq.info.sort_iden.txt $path/$file\_2refseq.info.txt ");
        #do blast
        
#       my $blast  = "blastall -p blastp -e 1e-5 -v 0 -b 5000 -F F -a $th -i $file -d $db -o $file\_2refseq";
		my $blast  = " blastp -evalue 1e-5 -num_descriptions 10000 -num_alignments 10000 -num_threads $th -query $path/$file -db $db -out $path/$file\_2refseq";        	
		system ("$blast");
        #parse blast output
        system ("blastout.identity_8alignments_2.pl $path/$file\_2refseq");
        system ("blastout_m7info.sort8iden.pl $path/$file\_2refseq.info");
        #select sequences
        system ("blastout_m8.phylaClass4.species.pl  $path/$file\_2refseq.info               $taxa_remove $tt $py $cs $sf");
        system ("blastout_m8.phylaClass4.species.pl  $path/$file\_2refseq.info.sort_iden.txt $taxa_remove $tt $py $cs $sf");
        system ("cat $path/$file\_2refseq.info.sort_iden.txt.txt >> $path/$file\_2refseq.info.txt");
        #system ("rm $file\_2refseq.info.sort_iden.txt.txt");
    	system ("cat $path/$file\_2refseq.info.txt >> $path/routine.qry2refseq2hitfasta.temp");

}
print "##-> retrieve sequences from $db\n";
system ("seqfetch.def.pl $path/routine.qry2refseq2hitfasta.temp $db");
return $c;
}


sub STEPTWO{##############################################step-2 seq trimal
# given a set of aligned sequences in fasta format, work it through trimal
# remove short sequence (ie, <40 AA) from alignment

print "#-> step-2 seqs trimal\n";

my @a =@_;

my $count=0;
foreach my $file (@a) {
      
        &FST2PHYLIP($file);

        $count++;
}
return $count;
}


sub STEPTHREE{##############################################step-3 reconstruct tree
my @a = @_;
my $countbuild=0;
foreach my $file (@a) {
       print $file."\n";
        &BUILDTREE($file);

        $countbuild++;
}
return $countbuild;
}



##################################################sub
sub FST2PHYLIP {

         my  $in = shift   ;
          open (IN,  $path."/".$in) or die "cant open file $in\n";
print "#-> $in\n";
print "#        muscle\n";
          my $query=$in;
          $in=$in."_2refseq.info.txt";
#         system ("mafft   --auto $in >$in.mft0 2>$tmp_file");

          system ("seqfetch.def.pl $path/$in $db2");
          system ("cat $path/$query >> $path/$in.fa.txt");
          print $in."\t".$query."\n";
          system ("perl -i -pe \'s\/\-\-\\S\+\-\-\/\-\/\' $path/$in.fa.txt");
		  system ("muscle -in $path/$in.fa.txt -out $path/$in.mft0 2> $tmp_file");
          system ("u.move_column.pl $path/$in.mft0 0.2 $path/$in.mft");          

print "#        trimAL\n";
          system ("trimal.longnames.pl $path/$in.mft");
          #system ("cp $in.mft $in.mft\-gb.txt");
          
          open (IN2, " $path/$in.mft\-gb.txt") or die "can open file $in.mft\-gb.txt\n";
          open (OT2, "> $path/$in.mft\-gb2.txt");          
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
          open (IN3,  "$path/$in.mft\-gb2.txt") or die "can open file $in.mft\-gb2.txt\n";
          open (OT3, ">$path/$in.mft\-gb2.txt8");          
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
                                       
          system ("fst.filterbylen.pl $path/$in.mft\-gb2.txt8 0.1 $path/$in.mft\-gb3.txt");
          
          #system ("tcoffee.rmsite_cutoff5.sh $in.mft\-gb3.txt");
          system ("cp $path/$in.mft\-gb3.txt $path/$in.mft\-gb3.txt.tcoffee");
             
          system ("u.move_column.pl $path/$in.mft\-gb3.txt.tcoffee 0.7 $path/$in.mft\-gb3.txt.tcoffee2");
          system ("fst.filterbylen.absolute.pl $path/$in.mft\-gb3.txt.tcoffee2 50 $path/$in.mft\-gb3.txt.tcoffee3");
          system ("fasta2relaxedPhylip.pl $path/$in.mft\-gb3.txt.tcoffee3 > $tmp_file 2>&1");
          system ("phy.rmorphan.keep_order.pl $path/$in.mft\-gb3.txt.tcoffee3.phylip");
          system ("cp $path/$in.mft\-gb3.txt.tcoffee3.phylip.sortbydef.txt $path/$in.phy.txt");
          system ("rm $path/$in.mft $path/$in.mft\-gb*txt* $tmp_file $path/$in.mft0");


}


        sub BUILDTREE{
        	my $in = shift   ;
        	print "#-> $in\n";
          	open (IN, $in) or die "cant open file $in\n";
			print "#-> $in\n";
        	#system ("iqtree -s $in.phy.txt -m MFP -mset WAG,LG,JTT -bb 1150");
	    	system ("FastTree -quiet -nopr -lg -spr 4 -mlacc 2 -slownni -out $path/$in.phy.tre $path/$in\_2refseq.info.txt.phy.txt > $tmp_file");
	    	#print "#-> done\n";
       		system ("rm $path/$in*mft**");
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
                      my $gap1 = '-'  x $gapsize1;
                      $amin = '\w' x $aminoacids;
                      my $gap2 = '-'  x $gapsize2;
                      my $pattern = "$gap1$amin$gap2";
                      $pattern = '^'.$pattern if $gapsize1 <2;    #pattern occurs at end of seqs
                      $pattern = $pattern.'$' if $gapsize2 <2;
                      $segm = '-'  x ($gapsize1 + $aminoacids + $gapsize2);
                      $line =~ s/$pattern/$segm/g;
                      #$pattern =~ s/\\//g;
                      #print "$pattern\n$segm\n\n";
                }
                return $line;
}

##################################################


sub Usage {
print <<End ;
The scripts is $0;
Usage:  
	perl $0 
	-db <database|str> -th <threads|num> -tt <total_Seq|num> -py <seqs_in_each_phylum|num> 
	-cs <seqs_in_each_class|num>  -sf <seqs_in_self_species|num> -pt <work_path|str>
	-fl  <seqs_file_list|str> -step <steps to run|num>
	-h/help     
default running: perl $0 -db database/1157 -th 3 -tt 64 -py 10 -cs 10 -sf 6 <seq_file>;

End

}
