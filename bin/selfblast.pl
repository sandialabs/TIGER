#! /usr/bin/perl
use warnings; use strict;

die "Usage: $0 fastafile [gff with add'l data] > outfile  # Find inverted repeat termini\n" unless $ARGV[0] and -f $ARGV[0];
my ($window, $tol) = (100, 10);

my ($head, $seq, %queries, %gff);
for (`cat $ARGV[0]`) {
 if (/^(>\S+)/) {Blast($head, $seq) if $seq; $seq = ''; $head = "$1\n"}
 else {chomp; $seq .= $_}
}
Blast($head, $seq) if $seq;

if ($ARGV[2]) {for (`cat $ARGV[2]`) {next unless /ID=([^;]+)/; chomp; $gff{$1} = $_}}
for my $q (sort keys %queries) {
 @{$queries{$q}} = sort {$$a[5] <=> $$b[5]} @{$queries{$q}};
 if ($gff{$q}) {die unless $gff{$q} =~ /intsum=([^;]*).*crossover=([^;]*).*tnps=([^;]*).*/; push @{$queries{$q}[0]}, $1, $2, $3}
 print join("\t", $q, @{$queries{$q}[0]}), "\n";
 # is607   78.431      2       51      2022    1975    2027   5
}

sub Blast {
 my ($head, $seq, $len, $max) = ($_[0], $_[1], length($_[1]));
 my $query; if ($len <= $window) {$query = $head . Revcomp($seq) . "\n"} else {$query = $head . Revcomp(substr($seq, -1 * $window)) . "\n"}
 `echo '$query' | makeblastdb -in - -out test -title test -dbtype nucl &> /dev/null`;
 if ($len <= $window) {$query = "$head$seq\n"} else {$query = $head . substr($seq, 0, $window) . "\n"}
 my @hits = `echo '$query' | blastn -task blastn -dust no -evalue 10 -word_size 4 -gapopen 3 -gapextend 5 -penalty -5 -reward 4 -query - -db test -outfmt 6 2> /dev/null`;
 # is607	is607	78.431	51	7	2	1975	2022	51	2	3.32e-05	36.8
 for (@hits) {
  #print "$_";
  my @f = split "\t";
  next if $f[9] < $f[8];
  next if $f[6] > $tol or $f[8] > $tol;
  my $max = $f[6] - 1; 
  $max = $f[8] - 1 if $f[8] > $f[6];
  for (8..9) {$f[$_] = $len - $f[$_] + 1}
  push @{$queries{$f[0]}}, [@f[2,6..9], $len, $max];
 }
}

sub Revcomp {my $ret = reverse $_[0]; $ret =~ tr/ACGTacgt/TGCAtgca/; return $ret}
