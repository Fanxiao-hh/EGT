#! /usr/bin/perl 
# delet colum with too much gaps

# input must be fasta format

# output format is fasta


$argm=join ("\*", @ARGV);
$cmarg=$#ARGV+1;
if ($argm =~ /\-+help/i||$cmarg<3) {
   die "usage: u.move_column.pl 1_infile(fasta) 2_cutoff(allowed percentage of gaps) 3_outfile\n";
} 


$infile = shift or die "input is wrong \?";
$cutoff = shift or die "input is wrong \?";

$outfile  = shift or die "input is wrong \?";



open (IN, "$infile") or die "cof";

open (OUT, ">$outfile")or die "cof";



%matrix=();  



while (<IN>) {

  if (/\>(\S+)/) {

push (@seq_title, $1);

$seq_title_tem = $1;

}

else { chomp;
s/x/\-/ig;                                                       #<<<---remove frameshift left by GENEWISE prediction

$matrix{$seq_title_tem} .= $_;

}

}


if ($cutoff =~ /^[\d\.]+$/) {

   @good_col = &choose_column_no (\%matrix, \@seq_title, $cutoff);
} else {
   @good_col = &choose_column_byReferenceSeq (\%matrix, \@seq_title, $cutoff);
}
$survival=$#good_col+1;

if ($cutoff =~ /^[\d\.]+$/) {
print "$survival columns contains residues(nucleotides) larger than $cutoff\n";
} else {
print "$survival columns contains residues(nucleotides) in sequence ($cutoff) have been exported\n";
}


&out_fasta (\%matrix, \@seq_title, \@good_col);





##-------------------------------------------------------------------

sub out_fasta {

  my ($r_matrix, $r_seq_title, $r_good_col) =@_;

      foreach (@$r_seq_title) {

      print OUT '>'.$_."\n" ;            #print $matrix{$_};

      foreach $j (@$r_good_col) {

      print OUT substr $$r_matrix{$_}, $j, 1; 

} 

      print OUT "\n";

}

}



##--------------------------------------------------------------------

sub choose_column_no {

  my ($r_matrix, $r_seq_title, $cutoff) = @_;

 $len_aln_block = length ($$r_matrix{$$r_seq_title[0]});



 #print '^'.$r_matrix.$r_seq_title.$cutoff.$len_aln_block.'*'."\n";

for $i (0..$len_aln_block - 1) {$data=0;

       foreach (@$r_seq_title) {

      $point = substr $$r_matrix{$_}, $i, 1;

      if ($point =~ /[a-z]/i) {++$data;}

}

     if ($data > $#$r_seq_title * $cutoff) {push (@col, $i);}

}

return @col;

}
##-------------------------------------------------------------------
sub choose_column_byReferenceSeq {

  my ($r_matrix, $r_seq_title, $refseq) = @_;
  my @col=();

 $len_aln_block = length ($$r_matrix{$$r_seq_title[0]});

 die "there is no $refseq in the input data\n" unless $$r_matrix{$refseq};



for $i (0..$len_aln_block - 1) {$data=0;

      $point = substr $$r_matrix{$refseq}, $i, 1;
      die "there is no $i",  "th residue/nucleotide for $refseq\n" unless $point;

      if ($point =~ /[a-z]/i) {push (@col, $i);}

}

return @col;

}
