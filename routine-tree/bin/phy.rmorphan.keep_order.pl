#! /usr/bin/perl

my @a = @ARGV;


foreach (@a) { ###


my $in = $_;
open (IN, $in);
open (OUT, ">$in.sortbydef.txt");

my $count = 1;
my %all=();
while (<IN>) {
       next if /^\s*$/;
       if (/^(\d+\s+\d+)/) {
           $header = $1;
       } else {
           my @a = split /\s+/;
           my $def = $a[0];
              $def =~ s/\_//g;
              $def =~ s/\s//g; 
           $all{$count}->{'seq'}=$a[1];
           $all{$count}->{'def'}=$a[0];
           $count++;
       }
}

my %all2= %all;
my @tmp = ();
@tmp = sort { length($all2{$b}->{'def'}) <=> length($all2{$a}->{'def'}) } keys %all2;
my $longest_def = "";
$longest_def = length ( $all2{$tmp[0]}->{'def'} );

#foreach (keys %all2) {
#        $longest_def = length ( $all2{$_}->{'def'} );
        #last;
#        print "$longest_def\n";
#}



print OUT "$header\n";
foreach (sort {$a <=> $b} keys %all) {
        my $def = $all{$_}->{'def'};
           $def =~ s/\_$//;
        my $add = " " x 100;
        my $def = $def.$add;
        my $def = substr ($def, 0, $longest_def);
        
        my $line = $all{$_}->{'seq'};
           $line = &Remove_orphan (2, 10,  2, $line);
                   
           $line = &Remove_orphan (3,  10, 0, $line);
           $line = &Remove_orphan (10, 15, 0, $line);
           
           $line = &Remove_orphan (0, 10,  3, $line);
           $line = &Remove_orphan (0, 15, 10, $line);
                   
           $line = &Remove_orphan (1,  15, 5, $line);
           $line = &Remove_orphan (5,  15, 1, $line);
                   
           $line = &Remove_orphan (3, 15, 3, $line);
           $line = &Remove_orphan (5, 20, 5, $line);

        print OUT "$def  $line\n";
}

}###



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