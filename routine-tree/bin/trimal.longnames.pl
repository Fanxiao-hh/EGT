#! /usr/bin/perl

my $in = shift;
my $out= "$in\-gb.txt";
open (OUT, ">$out");

         open (Aling1, $in);
         my %align1 = &fasta_gene(\*Aling1);
         my $sequence_number = scalar keys %align1;
#         die "!!! less than four sequences in the alignment\n" if $sequence_number <= 3;
         open (Aling_digit, ">$in.gb");

         my $seq_count=0;
         foreach (keys %align1) {
            my $seq_def = $_;
               $seq_count++;
            my $hash_key ="itree".$seq_count;
            s/\:|\(|\)//g;
            $digit2name{$hash_key}=$_;
            print Aling_digit ">$hash_key\n$align1{$seq_def}\n"; 
         }
         close Aling1;
         close Aling_digit;
         
        #my $gblock_task = "trimal -in gblocks.temp.$in -automated1 -htmlout gblocks.temp.$in.html >gblocks.temp.$in\-gb";
         my $gblock_task = "trimal -in $in.gb -automated1 -htmlout $in.gb.html >$in\-gb";         
         system ("$gblock_task");
#         unlink "alignments/$query.aln_digit-gb.htm";

    
    if (-s "$in\-gb") {
         open (Gbed, "$in\-gb");                               #check if alignment is empty
         while (<Gbed>) {
                     if (/(itree\d+)/g) {
                            my $leaf =$1;
                            my $seq_def=$digit2name{$leaf};  
                            s/(itree\d+)/$seq_def/g;                     
                     }
                     print OUT $_;
         
         }
   }      

#system ("fasta2relaxedPhylip.pl $out");
#system ("phy.sortbydef.pl $out.phylip");


sub fasta_gene {
#take the GLOB of filehandle of fasta sequence file
#return a HASH with key/value-pairs to be seq_definition/sequences-pairs
    my $in=shift;
    my %hash=(); 
    my $gene;
    while (<$in>){
      next if /^\s*$/;
      chomp($_);   
      if (/\>(\S+)/){
       $gene=$1;
       $hash{$gene}="";
      } else {
      s/b|j|x|z|O|U/\-/gi;
      my $seq = $_;
      $hash{$gene} .= "$seq\n";
      }
    }
    return %hash;
}

