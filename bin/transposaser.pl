#! /usr/bin/perl
use strict; use warnings;
# Takes existing .gff and reindexes for a subregion
use 5.10.0; use Cwd 'abs_path';

open (OUT, "> transposase.txt");
die "usage: perl $0 IslandFile\n" unless @ARGV == 1;
my $dir = abs_path($0); $dir =~ s/\/[^\/]+$//;
my $infile = $ARGV[0];
for (`cat $infile`) {
 my @f = split "\t";
 die $_ unless /ID=(([^;\.]+)[^;]+)/; 
 my %l = (dna => $f[0], type => $f[2], L => $f[3], R => $f[4], dir => $f[6]);
 while ($f[8] =~ s/^([^=;]+)=([^;]+);//g) {$l{$1} = $2};
 my ($name, $L, $R, $dna) = ($l{ID}, $f[3], $f[4], $f[0]);
  $name =~ s/\|//; 
 my $file = "Isles/$name/$name.gff";
 my @segment;
 for (`cat $file`) {
  chomp; 
  my @f = split "\t";
  next unless /annot=(T[^;]+)/ or /pfam1=(rve);/;
  push @segment, $1;
 }
 print OUT "$name\t", join(',', @segment), "\n";
}
 close (OUT); 

