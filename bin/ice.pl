#! /usr/bin/perl
use strict; use warnings;
use 5.10.0; use Cwd 'abs_path';

my ($nick, %l, %ICE, @ICE, @ICE_list, @other, @HYP, $ICEtot, $otherTot, $HYPtot, $perICE, $perHYP, $perother, $tot, $HYP_ICE, @ICEpfam);
my ($ICEcall) = (0);

die "Usage: perl $0 FastaFile\n" unless @ARGV >= 1;
my $dir = abs_path($0); $dir =~ s/\/[^\/]+$//;
my $file = $ARGV[0]; 
die unless $file =~ /(\S+)\.[^\.]+$/;
$nick = $1; 

mkdir 'protein';
system "$dir/hmmsearch --domtbl protein/${nick}_ICE.domtbl $dir/../db/ICE.hmm $file &> /dev/null";
system "perl $dir/../tater.pl -extraDoms protein/${nick}_ICE.domtbl $file" unless -f "$nick.gff";

for (`cat $dir/../db/ICEhmm.txt`){chomp; my @g = split "\t"; push @ICE_list, $g[0];}

for (`cat $nick.gff`){
  chomp; 
  next unless /\S+\t\S+\t\S+\t\S+\t\S+\t\S+\t\S+\t\S+/;
  my @f = split "\t";  
  %l = (dna => $f[0], type => $f[2], L => $f[3], R => $f[4], dir => $f[6]);
  while ($f[8] =~ s/^([^=;]+)=([^;]+);//g) {$l{$1} = $2};
  if ($l{pfam1} ~~ \@ICE_list) {push @ICE, $l{ID}; push @ICEpfam, $l{pfam1}}
   elsif ($l{product} eq 'hypothetical protein'){push @HYP, $l{ID}}
   else {push @other, $l{ID}};  
  $ICEtot = scalar @ICE; 
  $otherTot = scalar @other; 
  $HYPtot = scalar @HYP; 
 }

my $len = 0; for (`cat $file`) {next if /^>/; chomp; $len += length $_}
$tot = $ICEtot + $HYPtot + $otherTot;
$perICE = sprintf'%.0f',100 * $ICEtot/$tot;
$perother = sprintf'%.0f',100 * $otherTot/$tot;
$perHYP	= sprintf'%.0f',100 * $HYPtot/$tot;
$HYP_ICE = $perICE + $perHYP;

if ($len <= 10000){
  if ($ICEtot >= 7 or $perICE >= 15){$ICEcall = 1}
  elsif ($ICEtot > 2 or $perICE >= 12){$ICEcall=2}
} else {
  if ($ICEtot >= 7 or $perICE >= 15){$ICEcall = 1}
  elsif ($ICEtot > 2){$ICEcall = 2}
  elsif ($perICE > 7) {$ICEcall = 2}
}

open (OUT, "> ICE.txt");
print OUT "#nick\t\%ICE\t\%HYP\t\%other\tICEtot\tHYPtot\totherTot\ttot\tICEpfams\tICEcall\n";
$nick =~ s/.*\///;
print OUT "$nick\t$perICE\t$perHYP\t$perother\t$ICEtot\t$HYPtot\t$otherTot\t$tot\t", join(',', @ICEpfam), "\t$ICEcall\n"; 
close (OUT);



