#! /usr/bin/perl
use strict; use warnings;

die "Usage: $0 infasta outdir [AMG taxa]\n" if @ARGV < 1;
use Cwd 'abs_path';
my ($infile, $outdir, $binpath) = ($ARGV[0], $ARGV[1], $0); for ($infile, $outdir, $binpath) {$_ = abs_path($_)}
$binpath =~ s/\/([^\/]+)$//; my $scriptname = $1;
my $lib = $binpath; $lib =~ s/[^\/]*$/lib/;
my ($tax, %gencode, %final, %t, $id, %serials, $dna, $collect) = ('B');

$tax = $ARGV[2] if $ARGV[2];
LoadCode();
chdir $outdir;
open FINAL, ">tmrna.gff" or die "Can't open outfile $outdir/tmrna.gff: $!\n";
mkdir 'tmrna';
chdir 'tmrna';

Tmrna();
#system "perl $binpath/rfind.pl tmrna $infile $outdir/tmrna" unless -f 'tmrna.rfind';
exit unless keys %final;
open OUT, ">tmrna.raw.gff";
for my $id (sort keys %final) {
 $serials{final} ++;
 print OUT join("\t", $t{$id}{dna}, 'aragorn1.2.40', @{$t{$id}}{qw/type L R score ori/}, '.', "ID=aragorn.$serials{final};");
 for (qw/ivs_L ivs_R trunc_start trunc_end tag tag_L tag_R seq cca struct/)
 {if (defined $t{$id}{$_}) {print OUT "$_=$t{$id}{$_};"} else {print OUT "$_=;"}}
 warn join("\n", @{$t{$id}{notes}}) if $t{$id}{notes};
 print OUT 'notes=', join(',', @{$t{$id}{notes}}), ';' if $t{$id}{notes};
 print OUT "\n";
}
close OUT;
my %omits;
for (`intersectBed -a tmrna.raw.gff -b tmrna.raw.gff -wo -nonamecheck | awk -F"\\t" '\$9 != \$18 && \$6 <= \$15 {print}'`) {  # Note: will fail if spaces introduced into $f[8]
 #print $_;
 chomp; my @f = split "\t"; /ID=([^;]+);.*ID=([^;]+)/;
 if ($f[5] != $f[14]) {$omits{$f[8]} ++}
 else {warn "Overlapping tmrnas $1 and $2 have equal scores: $outdir/tmrna/tmrna.raw.gff\n"; $omits{(sort $f[8], $f[17])[0]} ++}
}
for (`cat tmrna.raw.gff`) {chomp; my @f = split "\t"; next if $omits{$f[8]}; print FINAL "$_\n"}
close FINAL;

# SUBROUTINES
sub Tmrna {
 #warn "Finding aragorn tmRNA gene calls\n";
 system "$binpath/aragorn1.2.40 -w -br -seq -e -m -l -o tmrna.aragorn $infile" unless -f "tmrna.aragorn";  # aragorn1.2.40
 for (`cat tmrna.aragorn`) {
  chomp;
  if (/^>(\S+)/) {$dna = $1; $id = ''}
  if (/^\d\s+tmRNA/) {
   warn "Cannot parse tmRNA line $_\n" unless /^\d+ +tmRNA\S* +c*\[([\d\-]+),(\d+)\]\s+(\S+)\s+(\d+),(\d+)\s+(\S+)/;
   $serials{tmrna} ++;
   $id = "tmrna.$serials{tmrna}";
   %{$t{$id}} = (dna => $dna, dir => 1, ori => '+', L => $1, R => $2, score => $3, tag_L => $4, tag_R => $5, tag => $6, type => 'tmRNA');
   $t{$id}{tag} =~ s/\*+$//;
   if (/c\[/) {$t{$id}{dir} = -1; $t{$id}{ori} = '-'}
   if ($t{$id}{L} < 1) {if ($t{$id}{dir} == 1) {$t{$id}{trunc_start} = -1*$t{$id}{L}} else {$t{$id}{trunc_end} = -1*$t{$id}{L}} $t{$id}{L} = 1}
   $t{$id}{type} = 'permuted_tmRNA' if /tmRNA\*/;
  } elsif (/^\d+/) {$id = ''}
  next unless $id;
  #my $chr print "$_\n";
  if (/^([a-z\.]+)\|([a-z\.]+)/) {@{$t{$id}}{qw/codseq accseq/} = ($1, $2)}
  elsif (/^([\(\.].*)\|(\S+)/)   {@{$t{$id}}{qw/codstr accstr/} = ($1, $2)}
  elsif (/^([a-z\.]+)$/) {$t{$id}{seq} .= $_}
 }
 for my $id (sort keys %t) {
  die $t{$id}{codstr} unless $t{$id}{codseq};
  push @{$t{$id}{notes}}, "Structure based on gapped codseq $t{$id}{codseq}" if $t{$id}{codseq} =~ /\./;
  push @{$t{$id}{notes}}, "Structure based on gapped accseq $t{$id}{accseq}" if $t{$id}{accseq} =~ /\./;
  my $seq = uc $t{$id}{seq}; my $len = length $seq;
  if ($t{$id}{tag_R} > $len) {
   if ($t{$id}{type} eq 'permuted_tmRNA') {$t{$id}{tag_R} = $len-6}
   else {$t{$id}{tag_R} = $len-28}
   push @{$t{$id}{notes}}, "Unstopped tag RF truncated to reasonable length";
  }
  my $tag = Translate(substr($seq, $t{$id}{tag_L}-1, $t{$id}{tag_R}-$t{$id}{tag_L}+1));
  while ($tag =~ s/Z$//) {$t{$id}{tag_R} -= 3}
  $seq = substr($seq, 0, $t{$id}{tag_L}-1) . lc(substr($seq, $t{$id}{tag_L}-1, $t{$id}{tag_R}-$t{$id}{tag_L}+1)) . substr($seq, $t{$id}{tag_R});
  for (@{$t{$id}}{qw/codstr accstr/}) {s/[a-z]/./g; s/\(/>/g; s/\)/</g}
  #print "$t{$id}{accstr}, $t{$id}{accseq}\n";
  my ($accdot, $coddot) = (0,0); $coddot = length $1 if $t{$id}{codstr} =~ /^(\.+)/; $accdot = length $1 if $t{$id}{accstr} =~ /(\.+)$/;
  my $excess = $accdot - $coddot - 1; $excess = 0 if $excess < 0; 
  $t{$id}{accstr} =~ s/\.{$excess}$//;
  $t{$id}{cca} = $1 if $t{$id}{accseq} =~ s/(.{$excess})$//;
  #if    ($t{$id}{accstr} =~ s/(\.)\.{3}$/$1/) {$t{$id}{accseq} =~ s/(...)$//; $t{$id}{cca} = $1}
  #elsif ($t{$id}{accstr} =~ s/(\.)\.{2}$/$1/) {$t{$id}{accseq} =~ s/(..)$//; $t{$id}{cca} = $1}
  #elsif ($t{$id}{accstr} =~ s/(\.)\.{1}$/$1/) {$t{$id}{accseq} =~ s/(.)$//; $t{$id}{cca} = $1}
  my %s;
  for (qw/cod acc/) {$t{$id}{seq} =~ /$t{$id}{"${_}seq"}/; $s{"${_}pre"} = $-[0]; $s{"${_}end"} = $+[0];}  # Coordinates of trna-like segments on pre-tmrna
  if ($t{$id}{type} eq 'permuted_tmRNA') {
   @{$t{$id}}{qw/ivs_L ivs_R/} = ($s{accend}+1, $s{codpre});
   $seq = substr($seq, 0, $t{$id}{ivs_L}-1) . lc(substr($seq, $t{$id}{ivs_L}-1, $t{$id}{ivs_R}-$t{$id}{ivs_L}+1)) . substr($seq, $t{$id}{ivs_R});
   my ($L,$i,$R) = ($s{accpre}, $s{codpre} - $s{accend}, length($seq) - $s{codend});
   $t{$id}{struct} = join('', "[$L]", $t{$id}{accstr}, "[$i]", $t{$id}{codstr}, "[$R]");
  } else {
   if ($t{$id}{cca}) {
    my $lencca = length $t{$id}{cca};
    $seq =~ s/.{$lencca}$//;
    if ($t{$id}{ori} eq '+') {$t{$id}{R} -= $lencca} else {$t{$id}{L} += $lencca}
   }
   my $i = $s{accpre} - $s{codend}; 
   $t{$id}{struct} = join('', $t{$id}{codstr}, "[$i]", $t{$id}{accstr});
  }
  $seq =~ s/^\.+//;
  $seq =~ s/\.+$//;
  $seq =~ s/\./N/g;
  $t{$id}{seq} = $seq;
  $final {$id} ++;
 }
}

#__END__
sub LoadCode {
 %gencode = qw(GCT A GCC A GCA A GCG A TGT C TGC C GAT D GAC D GAA E GAG E TTT F TTC F GGT G GGC G GGA G GGG G CAT H CAC H ATT I ATC I
   ATA I AAA K AAG K CTT L CTC L CTA L CTG L TTA L TTG L ATG M AAT N AAC N CCT P CCC P CCA P CCG P CAA Q CAG Q CGT R CGC R CGA R CGG R
   AGA R AGG R TCT S TCC S TCA S TCG S AGT S AGC S ACT T ACC T ACA T ACG T GTT V GTC V GTA V GTG V TGG W TAT Y TAC Y TAA Z TAG Z TGA Z);
 if ($tax eq 'M') {($tax, $gencode{TAG}) = qw/B Y/}  # ncbi gencode=4
 if ($tax eq 'G') {($tax, $gencode{TAG}) = qw/B G/}  # ncbi gencode=25
 #%gencode = qw (gct a gcc a gca a gcg a tgt c tgc c gat d gac d gaa e gag e ttt f ttc f ggt g ggc g gga g ggg g cat h cac h att i atc i ata i aaa k aag
 # k ctt l ctc l cta l ctg l tta l ttg l atg m aat n aac n cct p ccc p cca p ccg p caa q cag q cgt r cgc r cga r cgg r aga r agg r tct s tcc s tca s tcg s
 # agt s agc s act t acc t aca t acg t gtt v gtc v gta v gtg v tgg w tat y tac y taa z tag z tga z) ;
}

sub Translate {  # Translate dna sequence into protein
 my ($dna) = @_; 
 my $pep = '';
 my @codons = unpack ("a3" x int(length($dna)/3), $dna);
 for my $codon (@codons) {
  if ($codon =~ /[^ACGT]/) {$pep .= 'X'; next}
  if ($codon !~ /.../) {next}
  $pep .= $gencode{$codon};
 }
 return $pep;
}
