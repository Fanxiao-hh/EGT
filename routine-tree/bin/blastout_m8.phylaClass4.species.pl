#! /usr/bin/perl
# take the output of ps.blsxml.avrgIdtt.pl, and pick hit sequences under the
# constraint for number sequences from each phylum and class
# sequence identifier looks like:
# Plantae.Viridiplantae--Cucurbitales--Cucumis_sativus.GI.449440522
# cmp: first version
#      repeat the round of sequence retrieval till total sequences are 75 
# cmp: 1
#      able to receive taxa_name by which the associated sequences are not considered while sorting
# cmp: 2 
#      deal with more than one query sequences; hit quota is split among queries
# cmp: 3 
#      lines 70, 74, 121, 125: qrytaxa and hittaxa were retrieved at species
#                              level (eg, homo_sapian) instead of generus level (eg, homo)
#      lines 35, 69, 120: taxa provided by command line were not considered
#
# cmp: 0
#      put constraint on sequence number fron each taxa (=$self___s)

my $in      = shift;
my $def     = shift || 'zzzzzzzzz';  # hit match $def will be ignored

my $max_seq   = shift || 90;
my $phylum_s  = shift || 12;
my $class__s  = shift || 3;
my $self___s  = shift || 3;
my $iden_cut  = "0.3";


open  (IN, $in) or die "cof $in by blastout_m8.phylaClass4.species.pl\n";
open  (OUT,">$in.txt");

%blast=();
while (<IN>) {
       my @a = split /\t/;
       next unless $a[6] > $iden_cut;
       my $hit = $a[0];                         #next if $hit=~/$def/;
       my $qry = $a[2];
       next unless ( $a[1] =~ /^\d+$/ and $a[3] =~ /^\d+$/ );
       push ( @{ $blast{$qry} }, $_ );              
}
my $nb_qry = 1;
   $nb_qry = scalar keys %blast;

$max_seq = $max_seq/$nb_qry;

foreach (keys %blast) {
        my $qry  = $_;
        my @array = @{ $blast{$qry} };
        &Pick (\@array);
}



sub Pick  {   #------------------------------------------------------------------------------------- 

my $input = shift ;
my @ins   = @$input;
my %phyla  = ();
my %class  = ();
my %taxon  = ();
my $total  = 0 ;
my $selfs  = 0 ;
my @dump   = ();

my %passed=();
foreach (@ins) {
       my $line = $_;
       my @a = split /\t/;
       my $hit = $a[0];                         #next if $hit=~/$def/;
       my $qry = $a[2];       
       my $qry_id = "";
          if ($qry =~ /^[^\-]+\-([a-z]+\_*[a-z]+\S+)/ig) {;
              $qry_id = $1;           $qryaxa =~ s/X$//;
       }                                                                          #print "$qryaxa\t";
       my $hitaxa = "";
       my $hit_id = "";
          if ($hit =~ /^\S+--\S+--([a-z]+\_*[a-z]+)(\S+)/ig) {
              $hitaxa = $1;
              $hit_id = "$hitaxa$2";
          }                                                                       #print "$hitaxa\n";
          
       next unless $a[6] > $iden_cut;
       if ($qry_id) { next if ($hit_id eq $qry_id) };
       
       my $phylum = "";
       my $classs = "NA";
       my $taxass = "xx";
       
       
       if ($hit =~ /^([^\-]+)\-\-(.*)\-\-([a-z]+\_[a-z]+)/i) {
          $taxass = $3;
          $phylum = $1;
          $classs = $2;
          $classs = "NA" unless $classs =~ /\w/i;
          
       } elsif ($hit =~ /^([^\-]+)\-([a-z]+)/i) {
          $phylum = $1;       
          $classs = $2;          
       }   
          
    if ( $phyla{$phylum} >=$phylum_s or $class{$classs} >= $class__s or $taxon{$taxass} >=$self___s or $total >=$max_seq or $selfs >=$self___s) {
        #print "1 $hit\t$phylum $classs $taxass\n";
              push (@dump,$line); #print "a";
          } else {
               if ($def =~ /\w/ and $hitaxa eq $def ) {  
                                                            print "s $def\t$hit\n";
                   $total++;
                   $selfs++;
                   print OUT $line;
               }else {
                   next if exists $passed{$hit};
                   $phyla{$phylum}++;
                   $class{$classs}++;
                   $taxon{$taxass}++;
                   $total++;
                   print OUT $line;
               }
          }
}




    print OUT 'V' x 60, "\n";
while (1) {
    last if $total >= ($max_seq * 4)/5;
    last if $#dump <= 1;
           print "# ", $#dump, " lines recycled\n";

           
           my %phyla  = ();
           my %class  = ();
           my %taxon  = ();
           my @tempt  = ();
           foreach (@dump) {
                    my $line = $_;     
                    my @a = split /\t/;
                    
       my $hit = $a[0];                         next if $hit=~/$def/;
       my $qry = $a[2];       
       my $qryaxa = "";       
          if ($qry =~ /^[^\-]+\-([a-z]+\_*[a-z]+)/ig) {;
              $qryaxa = $1;           $qryaxa =~ s/X$//;
       }                                                                          #print "$qryaxa\t";
       my $hitaxa = "";
          if ($hit =~ /^\S+--\S+--([a-z]+\_*[a-z]+)/ig) {
              $hitaxa = $1;                
          }                                                                       #print "$hitaxa\n";
       if ($qry) { next if ($hitaxa eq $qryaxa) };
       
       next unless $a[6] > $iden_cut;
                    
                    my $phylum = "";
                    my $classs = "NA";
                    my $taxass = "xx";
               #if ($hit =~ /^([^\-]+)\-\-(.*)\-\-([a-z]+)/i) {
       if ($hit =~ /^([^\-]+)\-\-(.*)\-\-([a-z]+\_[a-z]+)/i) {
          $taxass = $3;
          $phylum = $1;
          $classs = $2;
          $classs = "NA" unless $classs =~ /\w/i;
          
       } elsif ($hit =~ /^([^\-]+)\-([a-z]+)/i) {
          $phylum = $1;       
          $classs = $2;          
       }   
                    
                    if ( $phyla{$phylum} >=$phylum_s/2 or $class{$classs} >= $class__s/2 or $taxon{$taxass} >=$self___s or $total >= ($max_seq * 4)/5 ) {
                        #print "2 $hit\t$phylum $classs $taxass\n";

                              push (@tempt,$line);
                    } else {
                              next if exists $passed{$hit};                    
                              $phyla{$phylum}++;
                              $class{$classs}++;
                              $taxon{$taxass}++;
                              $total++;
                              print OUT $line;
                    }              
           }
           @dump   = ();
           @dump   = @tempt;
    print OUT '$' x 60, "\n";
           #last;
}
print "# $total total sequences retrieved\n";
    print OUT '^' x 60, "\n";
}   #-----------------------------------------------------------------------------------------------
