
my $str =    'a0123456789*';
my $t = '2*';
if ($t=~m/\*/){
$str =~s/\*//;
}

print $str."\n";