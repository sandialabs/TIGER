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
 #die "$l{compose}\n";
 if ($l{compose} eq 'simple'){
  my $cmd = "$dir/collectSeq.pl -h -i genome.fa -e $l{dna} -L $l{L} -R $l{R} -d $l{ori} -n $name";
  warn "$cmd\n";
  my $fa = `$cmd`;
  print $fa;
 }
 else {
  system "mkdir Isles/$name.x";
  if (-r "Isles/$name.x/$name.fa") {next}
  open OUT, ">>Isles/$name.x/$name.fa" or die "$name";
  my $head = substr($name, 0, 15);
  # coord=JDTM01000012.1/5168-1+JDTM01000010.1/40589-36807;
  die "Can't parse coord $l{coord}\n" unless $l{coord} =~ /^([^\/]+)\/(\d+)-(\d+)\+([^\/]+)\/(\d+)-(\d+)/;
  #die "$1 $2 $3 $4 $5 $6\n";
  my $cmd  = "$dir/collectSeq.pl -h -i genome.fa -e $1 -c $2-$3 -n L.${head}";
  my $cmd2 = "$dir/collectSeq.pl -h -i genome.fa -e $4 -c $5-$6 -n R.${head}";
  warn "$cmd\n$cmd2\n";
  my $fa1 = `$cmd`;
  my $fa2 = `$cmd2`;
  print OUT $fa1 . $fa2;
  close OUT;
 }
}
