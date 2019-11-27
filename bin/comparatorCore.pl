#! /usr/bin/perl
use strict; use warnings;
# ToDo: Allow treating DNA as circle, alter at least q2 size to precisely maxcross (currently 2 bp longer, see couplingSeqTest/CTAGC-island)
# ToDo: collectSeq of genome to avoid blastdbcmd, see comparatorCore2.pl
# requires blast2.6 (gives accession.version form of subject)

die "Usage: perl $0 genome.fasta entry coordinate genomeDB minisle maxisle maxoverlap queryLength1 queryLength2 mobtype\n" unless @ARGV == 10;
my ($fa, $entry, $coord, $db, $minisle, $maxisle, $maxcross, $qlen1, $qlen2, $mobtype) = @ARGV; # min: 2000 for gi, 500 for is; max: 200000 for gi, 15000 for is; maxcross: 250 for gi, 30 for is
my ($minhit, $max_seqs) = (500, 1000000);
my (%sides, %isles, %reasons, $q1line, %actualQlen1, %xos);
$db = '../' . $db unless $db =~ /^\//; # Correct relative path to db
my $outdir = $entry; $outdir =~ s/^[a-z]+\|[0-9]+\|[a-z]+\|//; $outdir =~ s/\|.*//; $outdir .= ".$coord.$mobtype";

system "makeblastdb -dbtype nucl -in $fa -out genome -parse_seqids > /dev/null" unless -f 'genome.nhr';
my $cmd	= "blastdbcmd -db genome -dbtype nucl -entry $entry -outfmt %l";
print "$cmd\n";
my $gnmlen = `$cmd`;
mkdir $outdir; chdir $outdir;
for my $side (qw/L R/) { # Prepare q1 queries; q1 blast
 $actualQlen1{$side} = $qlen1;
 $actualQlen1{$side} -= 1 - ($coord-$qlen1+1)     if $coord-$qlen1+1 < 1       and $side eq 'L'; # Adjusts q1 len when reaching genomeL
 $actualQlen1{$side} -= $coord+$qlen1-1 - $gnmlen if $coord+$qlen1-1 > $gnmlen and $side eq 'R'; # Adjusts q1 len when reaching genomeR
 unless (-f "q1$side.fa") {
  my ($dir, $range) = ('plus',  ($coord-$actualQlen1{$side}+1) . "-$coord");
  ($dir, $range)    = ('minus', "$coord-" . ($coord+$actualQlen1{$side}-1)) if $side eq 'R';
  system "echo '>q1$side query:$entry:$range,$dir' > q1$side.fa";
  system "blastdbcmd -db ../genome -dbtype nucl -entry $entry -range $range -strand $dir -outfmt %s >> q1$side.fa";
 }
 next if $actualQlen1{$side} == 0;
 my $perc_qcov = int(100*$minhit/$actualQlen1{$side}); # Aiming for percent_qcoverage yielding 500 bp (minhit), but rounding errors make inaccurate
 $perc_qcov = "100" if $perc_qcov > 100;
 system "blastn -db $db -query q1$side.fa -max_target_seqs $max_seqs -qcov_hsp_perc $perc_qcov -outfmt 6 | sort -k 2,2 -k8,8nr -k 7,7n > q1$side.blast" unless -f "q1$side.blast";
}

open OUT, ">uninterrupteds.txt";
open FATE, ">fates.txt";
for my $side (qw/L R/) { # Process q1 blast hits; run and process q2 blast
 my ($gnmdir, $gnmstrand, $gnmname) = (1, 'plus', 'gnm'); ($gnmdir, $gnmstrand) = (-1, 'minus') if $side eq 'R';
 my (%rejects, %outs, %q1hits, @batch);
 if (-z "q1$side.blast") {print OUT "# $side side of coordinate: no Blast hits\n"; next;}
 open IN, "q1$side.blast";
 while (<IN>) {
  chomp;
  s/^(\S+)\t[a-zA-Z]+\|[0-9]+\|[a-zA-Z]+\|([^\|]+)\|\S*/$1\t$2/; # Strips GI line
  s/^(\S+)\t([^\.]+)\S*/$1\t$2/; # Removes version from accession
  $q1line = $_;
  #same island:             q1R     NC_011745.1     100.00  3000    0       0       1       3000    3052799 3049800 0.0     5541 (also, additional hit, to crossover)
  #uninterrupted reference: q1R     NZ_KE136876.1   99.94   1703    1       0       1       1703    564759  563057  0.0     3140 (no additional hit)
  #different island:        q1R     NC_011993.1     97.64   1822    31      12      1       1816    2771011 2769196 0.0     3116 (with additional hit, to crossover)
  my @f = split "\t";
  next if $rejects{$f[1]};
  $rejects{$f[1]} ++; # Allow only rightmost query hit to the reference genome
  $sides{$side} ++; 
  if ($f[7] > $actualQlen1{$side} - 10) {print FATE Fate('same island'); next} # Reject genome with same island (hit reaches coordinate with 10 bp tolerance)
  if ($f[3] < $minhit) {print FATE Fate('q1 hit too short'); next} # Match length is too short
  my ($refdir, $refstrand, $refname) = (1, 'plus', $f[1]);
  if ($f[8] > $f[9]) {($refdir, $refstrand) = (-1, 'minus')}
  my $start = $f[9] - $refdir * $maxcross; my $stop = $f[9] + $refdir * ($qlen2-$maxcross-1);
  my ($L, $R) = ($start, $stop); if ($refstrand eq 'minus') {($L, $R) = ($R, $L)}
  if ($L < 1) {print FATE Fate('incomplete q2'); next}; # Sequence request exceeds L end of reference; save 'exceeds R end' for later
  push @batch, "$f[1] $L-$R $refstrand\n";
  #die "$side $refname\n";
  %{$q1hits{$refname}} = (range => "$start-$stop", q1line => $q1line, refdir => $refdir, refstrand => $refstrand, refname => $refname);
 }
 close IN;
 unless (-f "q2$side.batch") {open BATCH, ">q2$side.batch"; print BATCH @batch; close BATCH}
 system "cat q2$side.batch | blastdbcmd -db $db -dbtype nucl -entry_batch - | sed ':a;N;/^>/M!s/\\n//;ta;P;D' |" . # Sed turns multi-line FastA entries to single-line
  "awk '!/^>/{next}{getline seq}length(seq)==$qlen2 {print \$0 \"\\n\" seq}' > q2$side.fa" unless -f "q2$side.fa";  # Awk removes entries ne $qlen2 (i.e., that exceed accession R end)
 system "cat q2$side.fa |  perl -pe 's/^>([^\\s:]+)\\S*/>$1/' | perl -pe 's/^>[a-zA-Z]+\\|[0-9]+\\|[a-zA-Z]+\\|([^\\|]+)\\|\\S*/>$1/'" if $1; # Strips GI line if there is one
 for (`grep -Po '^>\\S+' q2$side.fa`) {s/^>[a-zA-Z]+\|[0-9]+\|[a-zA-Z]+\|([^\|]+)\|\S*/>$1/; /^>([^\s:\.]+)/; $q1hits{$1}{fa} ++}
 system "blastn -db ../genome -query q2$side.fa -outfmt 6 | grep $entry | awk '\$4 >= $minhit' | awk '\$7 < " . ($maxcross+3) . "' | sort -k1,1 -k12,12nr > q2$side.blast" unless -f "q2$side.blast"; # will warn "Warning: [blastn] Query is Empty!" if no q1 hits; not an error message to worry about
 open IN2, "q2$side.blast";
 while (<IN2>) { # Collect top hit for each q2
  chomp;
  s/^[a-zA-Z]+\|[0-9]+\|[a-zA-Z]+\|([^\|]+)\|\S*/$1/; # Strips GI line
  s/^([^\s:]+)\S*/$1/;
  die "hi $_ $1 $2 $q1hits{$2}{range}\n" unless s/^(([^\s:\.]+)\S*)/$1:$q1hits{$2}{range}/ and $q1hits{$2} and $q1hits{$2}{range};
  $q1hits{$2}{refname} = $1;
  $q1hits{$2}{q2hit} = $_ unless $q1hits{$2}{q2hit}; # Keep only top hit
 }
 close IN2; 
 mkdir "overlap";
 for my $q1 (sort keys %q1hits) { #print "$q1, ", join(',', keys %{$q1hits{$q1}}), "\n"; next;
  # NZ_KE136876.1:563157-560158  NC_011745.1   100.00  2945    0       0       56      3000    3016097 3013153 0.0     5439
  $q1line = $q1hits{$q1}{q1line};
  die "$side $q1 $q1line\n" unless $q1line;
  unless ($q1hits{$q1}{fa}) {print FATE Fate('incomplete q2'); next} # q2 queries exceeding R end were omitted from faFile by awk filter, and not recorded in fa hash
  unless ($q1hits{$q1}{q2hit}) {print FATE Fate('no 1st-pass q2 hits'); next}
  my @f = split "\t", $q1hits{$q1}{q1line};
  my @g = split "\t", $q1hits{$q1}{q2hit};
  my ($refdir, $refstrand, $refname) = ($q1hits{$q1}{refdir}, $q1hits{$q1}{refstrand}, $q1hits{$q1}{refname});
  my $crossover = $maxcross - $g[6] + 2;
  if ($gnmdir * ($g[8] - $coord) < 0) {print FATE Fate('q2 same side as q1'); next} # q2 hit is on same side of coord as q1
  if (($g[9]-$g[8])/abs($g[9]-$g[8]) != $gnmdir) {print FATE Fate('flanks oppositely oriented'); next}	# Two island flanks must be in same orientation
  my ($prox, $distal) = ($coord-$gnmdir*($actualQlen1{$side}-$f[7]+int($crossover/2)), $g[8]+$gnmdir*int($crossover/2));
  my $islelen = abs($distal - $prox);
  if ($islelen > $maxisle) {print FATE Fate('isle too long' ); next} # Reject if island too long
  if ($islelen < $minisle) {print FATE Fate('isle too short'); next}
  if ($crossover >= $maxcross) {
   my ($L, $R, $proxM) = ($g[8] - $gnmdir * 2000, $g[8] + $gnmdir * ($crossover - 1), $coord - $gnmdir * ($actualQlen1{$side} - $f[7] + ($crossover - 1)/2)); # Distal overlap plus 2kb internal
   #die "$q1hits{$q1}{q1line}\n$q1hits{$q1}{q2hit}\n$side\t$f[7]\t$gnmdir\t$actualQlen1{$side}\nL: $L, R: $R, proxM: $proxM\n";
   ($L, $R) = ($R, $L) if $L > $R; $L = 1 if $L<1; $R = $gnmlen if $R>$gnmlen;
   system "blastdbcmd -db ../genome -dbtype nucl -entry '$entry' -range $L-$R -strand $gnmstrand > overlap/$L-$R" unless -f "overlap/$L-$R";
   system "blastn -db ../genome -query overlap/$L-$R -outfmt 6 > overlap/$L-$R.genome" unless -f "overlap/$L-$R.genome";
   open XO, "overlap/$L-$R.genome";
    while (<XO>) {
     chomp; my @h = split "\t";
     next unless Btwn ($h[8], $h[9], $proxM); # first hit that includes midpoint of prox
     my $newXo = $R-$L+2-$h[6];
     # die "$L-$R $h[8], $h[9], $proxM $newXo $crossover\n";
     $crossover = $newXo if $newXo > $crossover;
     last;
    }
   close XO;
   #if ($f[7]-$f[6]+1-$crossover < 50 or $g[7]-$maxcross < 50) {print FATE Fate('ref match not extended at both ends'); next}
   if ($f[7]-$f[6]+1-$crossover < 50 or $g[7]-$crossover < 50) {print FATE Fate('ref match not extended at both ends'); next} # Should be using the new crossover, not maxcross
   # find biggest hit spanning proxL and that match to ref extends at least 50 bp beyond this, reset $crossover accordingly
  }
  my ($proxL, $proxR, $refL, $refR, $distL, $distR) = ($coord - $gnmdir * ($actualQlen1{$side} - $f[7] + $crossover + 50 - 1), $coord - $gnmdir * ($actualQlen1{$side} - $f[7] - 50),
   $f[9] - $refdir * ($crossover + 50 - 1), $f[9] + $refdir * (50), $g[8] - $gnmdir * 50, $g[8] + $gnmdir * ($crossover + 50 - 1));
  if ($crossover > $maxcross) {
   my $adjust = $crossover - $maxcross - 1;
   ($distL, $distR) = ($g[8] - ($gnmdir * (50 + $adjust)), $g[8] + ($gnmdir * ($crossover + 50 - 1 + $adjust)))
  }
  #print "$q1hits{$q1}{q1line}\n$q1hits{$q1}{q2hit}\n$side\t$crossover\t$gnmdir\n$prox, $distal\nproxL: $proxL, proxR: $proxR\ndistL: $distL, distR: $distR\n";
  print FATE Fate('good');
  my $islesum = "$islelen:$prox-$distal";
  $isles{$islesum} ++;
  #if ($refL < 1) {$refL = 1} # to prevent refL from being negative if too close to end
  my $out .= ">$islesum; $q1; $q1line\nqueryStrand=$gnmstrand; referenceStrand=$refstrand; proximalL=$proxL, proximalR=$proxR; referenceL=$refL, referenceR=$refR; distalL=$distL, distalR=$distR\n";
  $out .= getSeq('../genome', $entry, $proxL, $proxR, $gnmstrand, $crossover);
  #die "$f[1] $q1hits{$q1}{refname} $entry $crossover $proxL, $proxR, $gnmstrand\n";
  $out .= getSeq($db, $q1hits{$q1}{refname},  $refL,  $refR,  $refstrand, $crossover);
  $out .= getSeq('../genome', $entry, $distL, $distR, $gnmstrand, $crossover);
  $out .= "crossover=$crossover; bitsum=" . ($f[11]+$g[11]) . "; hit=$q1hits{$q1}{q2hit}\n";
  push @{$outs{$islesum}{$f[11]+$g[11]}}, $out;
 }
 close IN;
 print OUT "# $side side of coordinate: $sides{$side} accessions hit; hit filtering fates:";
 for (sort {$reasons{$b} <=> $reasons{$a}} keys %reasons) {print OUT " '$_'=$reasons{$_};"} print OUT "\n";
 for my $isle (sort {keys(%{$outs{$b}}) <=> keys(%{$outs{$a}})} keys %outs) {
  for (sort {$outs{$isle}{$b} <=> $outs{$isle}{$a}} keys %{$outs{$isle}}) {print OUT join('', @{$outs{$isle}{$_}})}
 }
 %reasons = ();
}
close OUT;
close FATE;

sub Btwn {my ($L, $R, $M) = @_; ($L,$R) = ($R,$L) if $L>$R; return 1 if $M>=$L and $M<=$R; return 0;}
sub getSeq {
 my ($db, $entry, $L, $R, $strand, $crossover) = @_;
 if ($L > $R) {($L, $R) = ($R, $L)}
 my $seq = `blastdbcmd -db $db -dbtype nucl -entry '$entry' -range $L-$R -strand $strand -outfmt %s`;
 #die "$seq\n";
 chomp $seq;
 my $markseq = uc(substr($seq, 0, 50)) . lc(substr($seq, 50, $crossover)) . uc(substr($seq, 50+$crossover, 50)) . "\n";
 return $markseq;
}

sub Fate {my $reason = $_[0]; $reasons{$reason} ++; die "$reason\n" unless $q1line; return "$q1line\t$reason\n";}
