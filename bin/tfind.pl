#! /usr/bin/perl
use strict; use warnings;
use Cwd 'abs_path';

die "Usage: perl $0 tax directory\ndirectory must contain file called genome.fa\n" unless @ARGV == 2;  # Autodetect tax
my ($tax, $dir, $binpath) = (@ARGV, $0);
for ($dir, $binpath) {$_ = abs_path($_)}
$binpath =~ s/\/([^\/]+)$//;
#die "$dir $tax\n";
my (%bads, %badrfams, %seqs, $entry, %lens, $cmd, %struct);
mkdir $dir; chdir $dir; mkdir 'trna';
#`rm -r trna/rfind*`;
unless (-f "trna/rfind.gff")   {$cmd = "perl $binpath/rfind.pl tmrna genome.fa trna &> trna/rfind.log";    warn "$cmd\n"; `$cmd`}
unless (-f "trna/rfam.gff")    {$cmd = "perl $binpath/rfam.pl genome.fa trna &> trna/rfam.log";            warn "$cmd\n"; `$cmd`}
unless (-f "trna/introns.gff") {$cmd = "perl $binpath/introns.pl genome.fa $tax &> trna/introns.log"; warn "$cmd\n"; `$cmd`}
unless (-f "trna/tmrna.gff")   {$cmd = "perl $binpath/tmrna.pl genome.fa trna $tax &> trna/tmrna.log";     warn "$cmd\n"; `$cmd`}
unless (-f "trna/trna.gff")    {$cmd = "perl $binpath/trna.pl genome.fa trna $tax &> trna/trna.log";       warn "$cmd\n"; `$cmd`}
#exit if -f "trna/ttm.gff";
chdir 'trna';
mkdir 'anal';
`awk -F"\\t" '\$6 == 100 {print}' rfind.gff > anal/tmPure.gff`;
my %support; my %supptype = (qw/cmscan rfam aragorn1.2.40 aragorn/);
for (`intersectBed -nonamecheck -wo -f 0.1 -F 0.1 -e -s -a anal/tmPure.gff -b rfam.gff`, `intersectBed -nonamecheck -wo -f 0.1 -F 0.1 -e -s -a anal/tmPure.gff -b tmrna.gff`) {
 my @f = split "\t"; my $id = "$f[0].$f[3]";
 unless ($support{$id}{$f[10]} and $support{$id}{$f[10]} > $f[14]) {
  $support{$id}{$supptype{$f[10]}} = $f[14];
  $support{$id}{rfam_best} = $1 if $f[17] =~ /model=([^;]+)/ and $f[10] eq 'cmscan';
  $struct{$id} = $1 if /struct=([^;]+)/;
  #warn "$id $supptype{$f[10]} $support{$id}{$supptype{$f[10]}}\n";
 }
}
`intersectBed -nonamecheck -wa -a trna.gff -b anal/tmPure.gff | awk -F"\\t" '\$6 < 55 {print \$0 "reject=tmrna;"}' > anal/TFalseTm`;
`intersectBed -nonamecheck -wa -a trna.gff -b introns.gff | awk '{print \$0 "reject=gpi;"}' > anal/TFalseGpi`;
`cat anal/TFalseTm anal/TFalseGpi | sort -k1,1 -k4,4n -k5,5rn | intersectBed -nonamecheck -v -a trna.gff -b stdin > anal/Tremain.gff`;
`intersectBed -nonamecheck -wo -f 0.1 -F 0.1 -e -s -a anal/Tremain.gff -b anal/Tremain.gff | awk -F"\\t" '\$9 != \$18 \&\& \$6 <= \$15 {print}' | cut -f 1-9 | awk '{print \$0 "reject=better trna;"}' > anal/TFalseSelfoverlap`;
#system "pwd; ls ../";
#for (`cat ../genome.fa`) {if (/^>(\S+)/) {$entry = $1} else {chomp; $seqs{$entry} .= $_}}
for (`cat ../genome.fa`) {if (/^>(\S+)/) {$entry = $1} else {chomp; $seqs{$entry} .= $_}}
for (keys %seqs) {$lens{$_} = length $seqs{$_}}

open OUT, ">rejects.gff";
for (`cat rfind.gff tmrna.gff rfam.gff | sort -k1,1 -k4,4n -k5,5rn | intersectBed -nonamecheck -v -a stdin -b anal/tmPure.gff | cat - anal/TF* | sort -k1,1 -k4,4n -k5,5rn`) {
 my @f = split "\t";
 print OUT $_, next unless $f[1] eq 'cmscan';
 chomp;
 my ($dna, $L, $R, $ori, $Ltrunc, $Rtrunc, $seq) = ($f[0], $f[3]-301, $f[4]+300, $f[6], 0, 0, '');
 if ($L < 0)           {$Ltrunc = -1*$L; $L = 0}
 if ($R > $lens{$dna}) {$Rtrunc = $R-$lens{$dna}; $R = $lens{$dna}}
 $seq = uc (substr $seqs{$dna}, $L, $R-$L);
 substr $seq, -1 * (300-$Rtrunc), 0, ',,';
 substr $seq, 300-$Ltrunc, 0, ',,';
 $seq = Revcomp($seq) if $ori eq '-';
 print OUT "${_}seq=$seq;\n";
}
close OUT;
open OUT, ">ttm.gff";
open FA, ">ttm.fa";
for (`intersectBed -nonamecheck -v -a anal/Tremain.gff -b anal/TFalseSelfoverlap | cat - anal/tmPure.gff introns.gff | sort -k1,1 -k4,4n -k5,5rn`) {
 my %t = %{Gff($_)};
 push @{$t{ORDER}}, 'product', 'questionable';
 $t{questionable} = '';
 if (/;aa=His;/ and !$t{trunc_start}) {
  my @coord = ($t{ends}[0]-$t{ori}, $t{ends}[1]);
  if ($coord[0] >= 1 and $coord[0] <= $lens{$t{dna}}) {
   my @bp = (uc(substr($seqs{$t{dna}}, $coord[0]-1, 1)), uc(substr($seqs{$t{dna}}, $coord[1]-1, 1)));
   if ($t{ori} < 0) {for (@bp) {$_ = Revcomp($_)}}
   if ($bp[0] eq 'G') {
    if ($t{ori} < 0) {$t{R} ++} else {$t{L} --}
    for (qw/a_site j_site disc_fill disc_end gpi_L gpi_R ivs_L ivs_R/) {$t{$_} ++ if $t{$_}}
    $t{note} = join(',', $t{note}, 'G-1 included');
    $t{seq} = 'G' . $t{seq};
    if ($t{trunc_end} or $bp[1] ne 'C' or $bp[1] ne 'T') {$t{struct} = '.' . $t{struct}}
    else {$t{struct} =~ s/(.*).$/>$1</}
   } else {$t{note} = join(',', $t{note}, 'no G-1')}
  }
 }
 $t{note} =~ s/^,*(.*[^,]),*$/$1/ if $t{note};
 $t{product} = 'tmRNA';
 if ($t{type} eq 'tmRNA' and $t{form} eq 'prm') {
  if ($t{ivsL}) {$t{seq} =~ /^[^,]*,[^,]*,.{$t{ivsL}}(...)/; $t{cca} = $1}
 } elsif ($t{type} eq 'tmRNA') {
  if ($t{seq} =~ s/^([^,]*,[^,]*,[^,]*)([^,]{3}),/$1,/) {$t{cca} = lc $2}
  if ($t{ori} == 1) {$t{R} -= 3} else {$t{L} += 3}
 } else {  # trna
  /;aa=([^;]+);ac=([^;]+);/;
  $t{product} = "$t{type}-$1($2)";
  $t{questionable} = 'lowCoveScore' if $t{sco} < 50;
  $t{questionable} = $t{note} if $t{note};
  $t{questionable} = 'Undet' if $1 eq 'Undet';
  my $L = $t{ends}[1];
  if ($t{ori} < 0) {$L -= 4}
  elsif ($L+3 > $lens{$t{dna}}) {$L = -1}
  if ($L >= 0) {$t{cca} = lc (substr $seqs{$t{dna}}, $L, 3)}
  else {$t{cca} = ''}
  $t{cca} = lc (Revcomp(uc $t{cca})) if $t{ori} < 0;
 }
 $t{seq} =~ s/^[^,]*,[^,]*,([^,]+),.*/$1/;
 $t{ori} =~ s/1$//;
 print OUT join("\t", @t{qw/dna src type L R sco ori/}, '.', '');
 for (@{$t{ORDER}}) {print OUT "$_=$t{$_};"} #; die "$_ $t{dna}.$t{L}\n" unless defined $t{$_}}
 if ($t{type} eq 'tmRNA') {
  $t{j_site} = length($t{seq}) - 8; $t{j_site} = $t{ivsL} - 8 if $t{ivsL} and $t{form} eq 'prm';
  print OUT "j_site=$t{j_site};";
  my $id = "$t{dna}.$t{L}";
  print OUT "struct=$struct{$id};" if $struct{$id};
  print OUT "aragorn=$support{$id}{aragorn};" if $support{$id}{aragorn};
  print OUT "rfam=$support{$id}{rfam};rfam_best=$support{$id}{rfam_best};" if $support{$id}{rfam};
 }
 print OUT  "\n";
 if ($t{type} eq 'tmRNA' and $t{form} eq 'prm') {print FA ">$t{ID}\n$t{seq}\n"} else {print FA ">$t{ID}\n$t{seq}$t{cca}\n"}
}
close OUT;
close FA;

sub Revcomp {my $seq = reverse $_[0]; $seq =~ tr/ACGT/TGCA/; return $seq}
sub Gff {
 chomp;
 my @f = split "\t";
 my %t; @t{qw/dna src type L R sco ori/} = @f[0..6];
 $t{ori} .= 1;
 @{$t{ends}} = @f[3,4]; @{$t{ends}} = @f[4,3] if $f[6] eq '-';
 while ($f[8] =~ s/^([^=]+)=([^;]*;)//g) {my $key = $1; $t{$key} = $2; $t{$key} =~ s/;//; push @{$t{ORDER}}, $key}  # ; print "$key='$t{$key}'\t'$f[8]'\n"}
 return \%t;
}
