#! /usr/bin/perl
use strict; use warnings;
use File::Spec;
my $dir = File::Spec->rel2abs($0); $dir =~ s/\/[^\/]+$//;

my (%l, %c); 

die "Usage: perl $0 Island File\n" unless @ARGV == 1;
for (`cat $ARGV[0]`) {
 chomp;
 my @f = split "\t";
 %l = (dna => $f[0], type => $f[2], L => $f[3], R => $f[4], dir => $f[6], ori => 'F');
 while ($f[8] =~ s/^([^=;]+)=([^;]+);//g) {$l{$1} = $2};
 /ID=([^;]+);/;
 my $name = $1;
 $name =~ s/\|//g;
 $name =~ s/\'//g;
 $name =~ s/fa//;
 $name =~ s/\?/trna/;
 $name =~ s/\s+//;
 $name =~ s/\,/x/g;
 $name =~ s/\(/x/g;
 $name =~ s/\)/x/g;
($l{L}, $l{R}) = ($1, $2) if /endL=([^;]+).*endR=([^;]+)/; 
if ($l{dir} =~ /-/) {$l{ori} = 'R'}
 if ($l{compose} eq 'simple'){
  my $cmd = "$dir/collectSeq.pl -h -i genome.fa -e $l{dna} -L $l{L} -R $l{R} -d $l{ori} -n $name";
  warn "$cmd\n";
  my $fa = `$cmd`;
  print $fa;}
 else{
system "mkdir Isles/$name.x";
if (-r "Isles/$name.x/$name.fa") {next} 
else {
open OUT, ">>Isles/$name.x/$name.fa" or die "$name"; 
 my $head = substr($name, 0, 15);
   %c = (dna1 => $l{coord}, dna2 => $l{coord}, L1 => $l{coord}, L2 => $l{coord}, R1 => $l{coord}, R2 => $l{coord}, dir1 => "F", dir2 => "F");
   $c{dna1} =~ s/(Str_\d+_scaf\d+)\/\d+\-\d+\+.*/$1/;    
   $c{L1} =~ s/Str_\d+_scaf\d+\/(\d+)\-\d+\+.*/$1/; 
   $c{R1} =~ s/Str_\d+_scaf\d+\/\d+\-(\d+)\+.*/$1/;

   #$c{dna1} =~ s/(\D+\d+\.\d*)\/\d+\-\d+\+.*/$1/;    
   #$c{L1} =~ s/\D+\d+\.\d*\/(\d+)\-\d+\+.*/$1/; 
   #$c{R1} =~ s/\D+\d+\.\d*\/\d+\-(\d+)\+.*/$1/;
   my $cmd = "$dir/collectSeq.pl -h -i genome.fa -e $c{dna1} -c $c{L1}-$c{R1} -n L.${head}"; 
   warn "$cmd\n"; my $fa = `$cmd`;print OUT $fa;
   $c{dna2} =~ s/Str_\d+_scaf\d+\/\d+\-\d+\+(Str_\d+_scaf\d+)\/\d+\-\d+/$1/;
   $c{L2} =~ s/Str_\d+_scaf\d+\/\d+\-\d+\+Str_\d+_scaf\d+\/(\d+)\-\d+/$1/;
   $c{R2} =~ s/Str_\d+_scaf\d+\/\d+\-\d+\+Str_\d+_scaf\d+\/\d+\-(\d+)/$1/;

   #$c{dna2} =~ s/\D+\d+\.\d*\/\d+\-\d+\+(\D+\d+\.\d*)\/\d+\-\d+/$1/;
   #$c{L2} =~ s/\D+\d+\.\d*\/\d+\-\d+\+\D+\d+\.\d*\/(\d+)\-\d+/$1/;
   #$c{R2} =~ s/\D+\d+\.\d*\/\d+\-\d+\+\D+\d+\.\d*\/\d+\-(\d+)/$1/;
my $cmd2 = "$dir/collectSeq.pl -h -i genome.fa -e $c{dna2} -c $c{L2}-$c{R2} -n R.${head}";
   warn "$cmd2\n"; my $fa = `$cmd2`;print OUT $fa;}
close OUT;} 
}
