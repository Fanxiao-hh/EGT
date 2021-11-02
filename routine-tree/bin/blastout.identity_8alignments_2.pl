#! /usr/bin/perl
# substract query sequence from standard blast output
# join non-redundant query fragments if more than one alignments were 
# established between query and subject

$in = shift;
open (IN, "$in");
#$in=~s/\.txt//;
#open (OT,  ">$in.qryseq.txt");
#open (OT2, ">$in.qrybed.txt");
open (OT3,  ">$in.info");

$evalue_cut = "10";


$debug = 0;
$identifierQ = 0;
$identifierH = 0;
$evalue_pass = 0;
while (<IN>) {

       if (/^\S+\s+\d+\s+\S+\d\n/) {next;}                                              #print "------ $_ 00000 $identifier\n"; sleep (1);
       
       if (/BLASTP/) {
           # output statistics when encounter blast output for the next query started
           if (%all) {
              my ($out_line, $out_bed,$idenAA,$HSP,$iden,$posi) = &Parse_fragments_identity (\%all,$frame);
              #print OT ">$def\n$out_line";
              #print OT2 "$def\t$out_bed" ;
              my $hit = shift (@hits)    ;
              my $hit_coverage = $HSP/$hit_len;
                 $hit_coverage = substr ($hit_coverage, 0, 4);
              my $qry_coverage = $HSP/$qry_len;
                 $qry_coverage = substr ($qry_coverage, 0, 4);
              print OT3 "$hit\t$hit_len\t", "$def\t$qry_len\t", "$idenAA\t$HSP\t$iden\t", "$hit_coverage\t$qry_coverage\n";
           }
           %all=();
           $frame = "";
       }
       
       if (/Query=\s+(\S+)/i) {                                                         #print "*19\n";
           $align=0;    
           $def = $1;
           $identifierQ = 1; 
           <IN> ;                                                         #print "$def\n";
          $sss = <IN>;
       if ($sss=~m/Length=(\d+)/g){
            $qry_len = $1;
           $identifierQ = 0; 
           }                                                  #print "*23\n";
                                                                      #print "$def\n$qry_len\n" ;
       # query identifier is of multiple lines    
       } elsif ( $identifierQ) {
           if (/^([^\>]\s?\S+)\s+/g) { $def .= $1; }                                       #print "*28 $def\n"; sleep (2);
       #
       } elsif (/^>\s?(\S+)/g) {                                                           #print "*29\n";
           $encounter = $1;                                                             #print "\n#-->$encounter\n" if $debug;
           $identifierH = 1;
       } elsif (/Length=(\d+)/g) {                                                    #sleep (2); print "*32\n";
           push (@hits, $encounter);       
           $identifierH = 0;  
           # output statistics when encounter "Length = XXXX"
           if (%all) {
              my ($out_line, $out_bed,$idenAA,$HSP,$iden,$posi) = &Parse_fragments_identity (\%all,$frame);
              #print OT ">$def\n$out_line";
              #print OT2 "$def\t$out_bed" ;
              my $hit = shift (@hits)    ;
              my $hit_coverage = $HSP/$hit_len;
                 $hit_coverage = substr ($hit_coverage, 0, 4);
              my $qry_coverage = $HSP/$qry_len;
                 $qry_coverage = substr ($qry_coverage, 0, 4);
              print OT3 "$hit\t$hit_len\t", "$def\t$qry_len\t", "$idenAA\t$HSP\t$iden\t", "$hit_coverage\t$qry_coverage\n";
              #print "OT3  $hit\t$hit_len\t", "$def\t$qry_len\t", "$idenAA\t$HSP\t$iden\t", "$hit_coverage\t$qry_coverage\n";              
           }
           
           %all=();
           $frame = "";

           $hit_len = $1;
       # hit identifier is of multiple lines           
       } elsif ($identifierH) {
           if (/^\s+(\S+)\s+/g) {$encounter .= $1;}
       #    
       } elsif (/Expect = (\S+)/g) {
           $expect = $1;
           $evalue_pass = 1 if $expect <= $evalue_cut;
           print "#-> grant mercy, e-value small\n" unless $expect < $evalue_cut;
       } elsif (/Identities = (\d+)\/(\d+) .+ Positives = (\d+)\/\d+ .+ Gaps = (\d+)\//) {
           $identity = $1;
           $ttlength = $2;
           $positive = $3;
           $gaps     = $4;   
           if ($evalue_pass) {
               $align++;           
               $all{$align}->{similarity} = "$1 $2 $3 $4"; 

           }
           $evalue_pass = 0;
       } elsif (/Identities = (\d+)\/(\d+) .+ Positives = (\d+)\/\d+/) {
           $identity = $1;
           $ttlength = $2;
           $positive = $3;
           if ($evalue_pass) {
               $align++;
               $all{$align}->{similarity} = "$1 $2 $3 0"; 
           }
           $evalue_pass = 0;
       } elsif (/Frame\s+\=\s+(\S)/ and !$frame)  {
           $frame = $1;                                               #print "frame: $frame\n";
       } elsif (/Query:\s+(\d+)\s+(\S+)\s+(\d+)/) {                   #print "Query: $1\n";
           push ( @{ $all{$align}->{border}  }, $1, $3 );
           $all{$align}->{seq} .= $2."\n";
       } elsif (/Sbjct:\s+(\d+)\s+(\S+)\s+(\d+)/) {
           push ( @{ $all{$align}->{border_h} }, $1, $3 );
           #$sbjseq{$align} .= $2;
       } elsif (//) {
       
       }
}
           # output the last batch of data
           # output statistics when encounter "Length = 999"
           if (%all) {
              my ($out_line, $out_bed,$idenAA,$HSP,$iden,$posi) = &Parse_fragments_identity (\%all,$frame);
              print OT ">$def\n$out_line";
              print OT2 "$def\t$out_bed" ;
              my $hit = shift (@hits)    ;
              my $hit_coverage = $HSP/$hit_len;
                 $hit_coverage = substr ($hit_coverage, 0, 4);
              my $qry_coverage = $HSP/$qry_len;
                 $qry_coverage = substr ($qry_coverage, 0, 4);
              print OT3 "$hit\t$hit_len\t", "$def\t$qry_len\t", "$idenAA\t$HSP\t$iden\t", "$hit_coverage\t$qry_coverage\n";
           }




sub Parse_fragments_identity{
     my $r_hash = shift;
     my $frame  = shift;
     my %hash   = %{$r_hash};
     my %seq    = ();
     my %exon   = ();
     my $out    = "";
     my $out2   = "";
     
     my $identities = 0;
     my $positives  = 0;
     my $gaps       = 0;
     my $ttlengths  = 0;
     
     # go through each HSP
     foreach (sort {$a<=>$b} keys %hash) {                                  #print "key: $_\n";
        # export %hash information to novel local hashes such as %exon and %seq
        next unless exists $hash{$_};
        my @borders = @{ $hash{$_}->{border} };
        my $intv_strand = $borders[0].'_'.$borders[-1];
           @borders     = sort {$a<=>$b} @borders;
        my $intv        = $borders[0].'_'.$borders[-1];                     #print "\$intv $intv\n";
           $seq{$intv_strand}  = $hash{$_}->{seq};
           $exon{$intv_strand} = $borders[0]."\t".$borders[-1];
           
        my @borders_h = @{ $hash{$_}->{border_h} };
        my $intv_strand_h = $borders_h[0].'_'.$borders_h[-1];
           @borders_h     = sort {$a<=>$b} @borders_h;
        my $intv_h        = $borders_h[0].'_'.$borders_h[-1];               #print "\$intv_h $intv_h\n";

        # calculate statistics   
        my $statistics = $hash{$_}->{similarity};
        my @numbers    = split (/ /, $statistics);                          print "--> individual\t@numbers\n" if $debug;
        $identities += $numbers[0];
        $ttlengths  += $numbers[1];
        $positives  += $numbers[2];
        $gaps       += $numbers[3];
                                                                            print "#->  $identities\t$positives\t$gaps\t$ttlengths\n" if $debug;
        
        
        delete $hash{$_};
        
        # go through the rest of HSPs
        # delete HSPs if overlap with the current HSP
        my %hash_tmp    = %hash;
        foreach (keys %hash_tmp) {                                          #print "nested key: $_\n";
           my @borders2 = @{ $hash{$_}->{border} };
              @borders2 = sort {$a<=>$b} @borders2;
           my $intv2    = $borders2[0].'_'.$borders2[-1];   
           my $overlap = &Overlap_or_not_reference($intv, $intv2, 6);  # allowing 60bp overlap
           
           my @borders2_h = @{ $hash{$_}->{border_h} };
              @borders2_h = sort {$a<=>$b} @borders2_h;
           my $intv2_h    = $borders2_h[0].'_'.$borders2_h[-1];   
           my $overlap_h = &Overlap_or_not_reference($intv_h, $intv2_h, 6);  # allowing 60bp overlap
           
           #print "----- $overlap ----- $overlap_h\n";

           if ($overlap) {
              delete $hash{$_};   #print 'delete $hash{$_}; *****', "\n";
           } elsif ($overlap_h) {
              delete $hash{$_};   #print 'delete $hash{$_}; #####', "\n";        
           }
        }
     }
     if ($frame eq '-') {
         foreach (sort {$b <=> $a} keys %seq) {
                my ($a,$b) = split /\_/;
                next if ($a-$b) < 0;
                $out .= "$seq{$_}\n";
                $out2 .= "$exon{$_}\t$frame\n";
                #print ">$_\n$seq{$_}\n";
         }
     } else {
         foreach (sort {$a <=> $b} keys %seq) {
                my ($a,$b) = split /\_/;
                next if ($a-$b) > 0;
                $out .= "$seq{$_}\n";
                $out2 .= "$exon{$_}\t$frame\n";                
                #print ">$_\n$seq{$_}\n";
         }     
     }
     $ttlengths = $ttlengths-$gaps;
     my $p_identity = $identities/($ttlengths);
     my $p_positive = $positives/($ttlengths);                      print "-----1> $p_identity\t$p_positive\n" if $debug;
     return ($out,$out2,$identities,$ttlengths,$p_identity,$p_positive);  
}


sub Overlap_or_not_reference {
# judge redundance based on overlapping percentage
# allow 30-site flexibility
    my $intv_ref = shift;
    my $intv_sub = shift;
    my $cutoff   = shift;
    my $overlap = 0;
    my ($l_ref,$r_ref) = split (/\_/, $intv_ref);
    my ($l_sub,$r_sub) = split (/\_/, $intv_sub);
    my $length_ref = abs ($r_ref-$l_ref+1);
    my $length_sub = abs ($r_sub-$l_sub+1);                         #print "#interval $l_ref\t$l_sub\t$r_sub\t$r_ref\n";
    if       ($l_ref<=$l_sub and $r_sub<=$r_ref) {
        $overlap = 1;
    } elsif  ($l_sub<=$l_ref and $r_ref<=$r_sub) {
        $overlap = 1 if $length_ref >=$cutoff;        #die "impossible, reference is covered by subject\n"; #
    }  elsif ($l_ref<=$l_sub and $l_sub<=$r_ref) {
        my $overlap_len = $r_ref-$l_sub;
        $overlap = 1 if $overlap_len >=$cutoff;
    } elsif  ($l_ref<=$r_sub and $r_sub<=$r_ref) {
        my $overlap_len = $r_sub-$l_ref;
        $overlap = 1 if $overlap_len >=$cutoff;    
    }
    return $overlap;
}

