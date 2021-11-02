#! /usr/bin/perl 
# fltfstbylen.pl
$usage = "fltfstbylen.relative.pl___retrieve sequence longer than specified percentage of the longest seq\n",
         "usage: fltfstbylen.relative.pl 1-input(fasta) 2-minimum percetage 3-output\n";

  my $in    = shift or die "$usage";
  my $minlen= shift or die "$usage";
  my $out   = shift or die "$usage";

open (IN, "$in"); 
open (OUT, ">$out"); 


    my %hash=(); 

    my $gene;
    my $seqno=0;
    my $pass=0;
 

    while (<IN>){

       if (/\>(\S+)/){
       ++$seqno;

       $gene=$1;
       $hash{$gene}->{'num'}=$seqno;
       $hash{$gene}->{len}=0;
       $hash{$gene}->{seq}=();

       } else {
       next if /^\s*$/;

       $seqbk=$_;
       chomp;
       s/\-|\s|\r//g;   

       my $seq=$_;

       my $len = length ($seq);
       $hash{$gene}->{seq} .= $seqbk;

       $hash{$gene}->{len} += $len;

       }

    }

foreach (keys %hash) {
  push ( @len_sum, $hash{$_}->{len} );
}
@len_sum = sort {$b <=> $a} @len_sum;
$longest = $len_sum[0];
$min     = $longest * $minlen;

 foreach (sort { $hash{$a}->{'num'} <=> $hash{$b}->{'num'} } keys %hash) {
    if ($_ !~ /plugin/i) {
       if ($hash{$_}->{len} >= $min) {
          ++$pass;
          print OUT '>'.$_."\n",
                 $hash{$_}->{seq};
       }
    }
 }

print "there are $seqno seqs in input data\n",
      "$pass seqs are longer (not including gaps) than $min ($minlen\*$longest)\n";
