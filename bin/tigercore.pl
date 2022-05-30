#! /usr/bin/perl
use strict; use warnings;
use Cwd 'abs_path';
# requires blast2.6 (gives accession.version form of subject)
# CROSS: use whole genome for blast2, use dnastats for length

die "Usage: perl $0 genome.fasta entry coordinate genomeDB minisle maxisle maxoverlap queryLength1 queryLength2 mobtype mode\n" unless @ARGV >= 10;
my ($fa, $entry, $coord, $db, $minisle, $maxisle, $maxcross, $qlen1, $qlen2, $mobtype, $mode) = @ARGV; # min: 2000 for gi, 500 for is; max: 200000 for gi, 15000 for is; maxcross: 250 for gi, 30 for is
die "Mode must be 'intact', 'cross' or 'circleOrigin'\n" unless $mode =~ /^(intact|cross|circleOrigin)$/;
my ($minhit, $extra, $max_seqs) = (500, '', 1000000);
$extra = "| grep $entry" unless $mode eq 'cross';
my (%sides, %isles, %reasons, $q1line, %actualQlen1, %xos, %gnmlens);
$db = abs_path($db);  # Reference genome DB
my $outdir = $entry; $outdir =~ s/^[a-z]+\|[0-9]+\|[a-z]+\|//; $outdir =~ s/\|.*//; $outdir .= ".$coord.$mobtype";

RunCommand("makeblastdb -dbtype nucl -in $fa -out genome -parse_seqids > /dev/null", 'genome.nhr');
my $cmd	= "blastdbcmd -db genome -dbtype nucl -entry all -outfmt '%a %l'";
#print STDERR "Running blastdbcmd -db genome -dbtype nucl -entry all -outfmt '%a %l'\n\n";
for (`$cmd`) {die $_ unless /^(\S+) (\d+)$/; $gnmlens{$1} = $2}
mkdir $outdir; chdir $outdir;
print "command '$cmd' yielded ", scalar(keys %gnmlens), " entry lengths\n";
for my $side (qw/L R/) { # Prepare q1 queries; q1 blast
 $actualQlen1{$side} = $qlen1;
 $actualQlen1{$side} -= 1 - ($coord-$qlen1+1)     if $coord-$qlen1+1 < 1       and $side eq 'L'; # Adjusts q1 len when reaching genomeL
 $actualQlen1{$side} -= $coord+$qlen1-1 - $gnmlens{$entry} if $coord+$qlen1-1 > $gnmlens{$entry} and $side eq 'R'; # Adjusts q1 len when reaching genomeR
 unless (-f "q1$side.fa") {  # q1 R end is within putative island
  my $Lend = $coord-$actualQlen1{$side}+1;
     $Lend = 1 if $Lend < 1;
  my $Rend = $coord+50;
     $Rend = $gnmlens{$entry} if $Rend > $gnmlens{$entry};
  my ($dir, $range) = ('plus',  ($Lend) . '-' . ($coord+50));
     $Lend = $coord-50;
     $Lend = 1 if $Lend < 1;
     ($dir, $range) = ('minus', ($Lend) . '-' . ($coord+$actualQlen1{$side}-1)) if $side eq 'R';
  RunCommand("echo '>q1$side query:$entry:$range,$dir' > q1$side.fa", "q1$side.fa");
  #print STDERR "blastdbcmd -db ../genome -dbtype nucl -entry $entry -range $range -strand $dir -outfmt %s >> q1$side.fa\n\n";
  RunCommand("blastdbcmd -db ../genome -dbtype nucl -entry $entry -range $range -strand $dir -outfmt %s >> q1$side.fa", '');
 }
 next if $actualQlen1{$side} == 0;
 my $perc_qcov = int(100*$minhit/$actualQlen1{$side}); # Aiming for percent_qcoverage yielding 500 bp (minhit), but rounding errors make inaccurate
 $perc_qcov = "100" if $perc_qcov > 100;
 RunCommand("blastn -db $db -query q1$side.fa -max_target_seqs $max_seqs -qcov_hsp_perc $perc_qcov -outfmt 6 | sort -k 2,2 -k8,8nr -k 7,7n > q1$side.blast", "q1$side.blast");
}

open OUT, ">uninterrupteds.txt";
open FATE, ">fates.txt";
for my $side (qw/L R/) { # Process q1 blast hits; run and process q2 blast
 my ($gnmdir, $gnmstrand, $gnmname) = (1, 'plus', 'gnm'); ($gnmdir, $gnmstrand) = (-1, 'minus') if $side eq 'R';
 my (%rejects, %outs, %q1hits, @batch);
 if (-z "q1$side.blast") {print OUT "# $side side of coordinate: no Blast hits\n"; next;}
 for (`cat q1$side.blast`) {
  chomp;
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
  if ($f[8] > $f[9]) {($refdir, $refstrand) = (-1, 'minus')}  # q2 L end spans attB
  my $start = $f[9] - $refdir * $maxcross; my $stop = $f[9] + $refdir * ($qlen2-$maxcross-1);
  my ($L, $R) = ($start, $stop); if ($refstrand eq 'minus') {($L, $R) = ($R, $L)}
  if ($L < 1) {print FATE Fate('incomplete q2'); next}; # Sequence request exceeds L end of reference; save 'exceeds R end' for later
  push @batch, "$refname $L-$R $refstrand\n";
  #die "$start-$stop, $q1line, $refdir, $refstrand, $refname" if $refname =~ /c524204/; #/UFOJ01000004.1/;
  %{$q1hits{$refname}} = (range => "$start-$stop", q1line => $q1line, refdir => $refdir, refstrand => $refstrand, refname => $refname);
 }
 unless (-f "q2$side.batch") {open BATCH, ">q2$side.batch"; print BATCH @batch; close BATCH}
 #print STDERR "cat q2$side.batch | blastdbcmd -db $db -dbtype nucl -entry_batch - | sed ':a;N;/^>/M!s/\\n//;ta;P;D' |" . # Sed turns multi-line FastA entries to single-line
 # "awk '!/^>/{next}{getline seq}length(seq)==$qlen2 {print \$0 \"\\n\" seq}' > q2$side.fa\n\n";
 RunCommand("cat q2$side.batch | blastdbcmd -db $db -dbtype nucl -entry_batch - | sed ':a;N;/^>/M!s/\\n//;ta;P;D' |" . # Sed turns multi-line FastA entries to single-line
  "awk '!/^>/{next}{getline seq}length(seq)==$qlen2 {print \$0 \"\\n\" seq}' > q2$side.fa", "q2$side.fa");  # Awk removes entries ne $qlen2 (i.e., that exceed accession R end)
 for (`grep '^>' q2$side.fa`) {/^>([^\s:]+)/; $q1hits{$1}{fa} ++}
 
 RunCommand("blastn -db ../genome -query q2$side.fa -outfmt 6 $extra | awk '\$4 >= $minhit' | awk '\$7 < " . ($maxcross+3) .
  "' | sort -k1,1 -k12,12nr > q2$side.blast", "q2$side.blast");
  # Above will warn "Warning: [blastn] Query is Empty!" if no q1 hits; not an error message to worry about
 #next;
 for (`cat q2$side.blast`) { # Collect top hit for each q2
  chomp;
  /^([^\s:]+)/;
  $q1hits{$1}{q2hit} = $_ unless $q1hits{$1}{q2hit}; # Keep only top hit for each query2
 }
 #for (keys %q1hits) {print "$_\n"}; exit;
 mkdir "overlap";
 for my $q1 (sort keys %q1hits) {
  # NZ_KE136876.1:563157-560158  NC_011745.1   100.00  2945    0       0       56      3000    3016097 3013153 0.0     5439
  $q1line = $q1hits{$q1}{q1line}; 
  die "$side $q1 $q1line\n" unless $q1line;
  unless ($q1hits{$q1}{fa}) {print FATE Fate('incomplete q2'); next} # q2 queries exceeding R end were omitted from faFile by awk filter, and not recorded in fa hash
  unless ($q1hits{$q1}{q2hit}) {print FATE Fate('no 1st-pass q2 hits'); next}
  my @f = split "\t", $q1hits{$q1}{q1line};
  my @g = split "\t", $q1hits{$q1}{q2hit};
  my ($refdir, $refstrand, $refname) = ($q1hits{$q1}{refdir}, $q1hits{$q1}{refstrand}, $q1hits{$q1}{refname});
  my ($hitdir, $hitstrand, $hitentry) = (($g[9]-$g[8])/abs($g[9]-$g[8]), 'plus', $g[1]); $hitstrand = 'minus' if $hitdir < 0;
  if ($mode eq 'intact' and $gnmdir * ($g[8] - $coord) < 0) {print FATE Fate('q2 same side as q1'); next} # q2 hit is on same side of coord as q1
  if ($hitentry eq $entry and $hitdir != $gnmdir) {print FATE Fate('flanks oppositely oriented'); next}  # Two same-contig island flanks must be in same orientation
  my ($hitsum1, $hitsum2) = (join('/', @f[2,6,7,8,9]), join('/', @g[2,6,7,8,9])) ;
  my $compos = 'simple';
  if ($hitentry ne $entry) {$compos = 'cross'} elsif ($gnmdir * ($g[8] - $coord) < 0) {$compos = 'circlejxn'}
  my $crossover = $maxcross - $g[6] + 2;
  my ($prox, $distal) = ($coord-$gnmdir*($actualQlen1{$side}-$f[7]+int($crossover/2)), $g[8]+$hitdir*int($crossover/2));
  my (@halves, $islelen);
  if ($compos eq 'simple') {$islelen = abs($distal - $prox)} else {
   @halves = (End($entry, $gnmdir), End($hitentry, -1*$hitdir));
   $islelen = abs($prox - $halves[0]) + abs($distal - $halves[1]);
  }
  if ($islelen > $maxisle) {print FATE Fate('isle too long' ); next} # Reject if island too long
  if ($islelen < $minisle) {print FATE Fate('isle too short'); next}
  if ($crossover >= $maxcross) {
   print FATE Fate('crossover too long'); next;
   #my ($L, $R, $proxM) = (
   # $g[8] - $hitdir * 2000,
   # $g[8] + $hitdir * ($crossover - 1),
   # $coord - $gnmdir * ($actualQlen1{$side} - $f[7] + ($crossover - 1)/2)
   #); # Distal overlap plus 2kb internal
   #($L, $R) = ($R, $L) if $L > $R; $L = 1 if $L < 1; $R = $gnmlens{$hitentry} if $R > $gnmlens{$hitentry};
   #my $base = "overlap/$hitentry.$L-$R";
   #RunCommand("blastdbcmd -db ../genome -dbtype nucl -entry '$hitentry' -range $L-$R -strand $hitstrand > $base", $base);
   #RunCommand("blastn -db ../genome -query $base -outfmt 6 > $base.genome", "$base.genome");
   #for (`cat $base.genome`) {
   # chomp; my @h = split "\t";
   ## next unless $h[1] eq $entry;  # New to crossTiger, fixes old problem?
   # next unless Btwn ($h[8], $h[9], $proxM); # first hit that includes midpoint of prox
   # my $newXo = $R-$L+2-$h[6];
    # die "$L-$R $h[8], $h[9], $proxM $newXo $crossover\n";
   # $crossover = $newXo if $newXo > $crossover;
   # last;
   #}
   #if ($f[7]-$f[6]+1-$crossover < 50 or $g[7]-$crossover < 50) {print FATE Fate('ref match not extended at both ends'); next}
   # find biggest hit spanning proxL and that match to ref extends at least 50 bp beyond this, reset $crossover accordingly
  }
  my ($proxL, $proxR, $refL, $refR, $distL, $distR) = ($coord - $gnmdir * ($actualQlen1{$side} - $f[7] + $crossover + 50 - 1),
   $coord - $gnmdir * ($actualQlen1{$side} - $f[7] - 50), $f[9] - $refdir * ($crossover + 50 - 1),
   $f[9] + $refdir * (50), $g[8] - $hitdir * 50, $g[8] + $hitdir * ($crossover + 50 - 1));
  $proxL = 1 if $proxL < 1;
  $proxR = $gnmlens{$entry} if $proxR > $gnmlens{$entry};
  $refL = 1 if $refL < 1;
  $refR = $gnmlens{$entry} if $refR > $gnmlens{$entry};
  $distL = 1 if $distL < 1;
  $distR = $gnmlens{$entry} if $distR > $gnmlens{$entry};

  #if ($crossover > $maxcross) {
  # my $adjust = $crossover - $maxcross - 1;
  # ($distL, $distR) = ($g[8] - ($hitdir * (50 + $adjust)), $g[8] + ($hitdir * ($crossover + 50 - 1 + $adjust)))
  #}
  my %coords = (proximal => [$proxL, $proxR, $gnmstrand], ref => [$refL, $refR, $refstrand], distal => [$distL, $distR, $hitstrand]); 
  for (qw/proximal ref distal/) {@{$coords{$_}} = ($coords{$_}[1], $coords{$_}[0]) if $coords{$_}[2] eq 'minus'}
  print FATE Fate('good');
  my $islesum = "$islelen:$entry/$prox-$distal";
  my $bittot = $f[11]+$g[11];
  if (@halves) {$islesum = "$islelen:$entry/$prox-$halves[0]+$hitentry/$halves[1]-$distal"}
  $isles{$islesum} ++;
  #if ($refL < 1) {$refL = 1} # to prevent refL from being negative if too close to end
  my $out .= ">$islesum; bitsum=$bittot($f[11]+$g[11]); length=$islelen; crossover=$crossover; compose=$compos\n" .
   "proximal:$entry,$hitsum1,$coords{proximal}[0]-$coords{proximal}[1] reference:$q1,$coords{ref}[0]-$coords{ref}[1] " .
   "distal:$hitentry,$hitsum2,$coords{distal}[0]-$coords{distal}[1]\n";
  $out .= getSeq('../genome', $entry,    $proxL, $proxR, $gnmstrand, $crossover);
  $out .= getSeq($db,         $refname,  $refL,  $refR,  $refstrand, $crossover);
  $out .= getSeq('../genome', $hitentry, $distL, $distR, $hitstrand, $crossover);
  $outs{$islesum}{$out} = $bittot;
 }
 print OUT "# $side side of coordinate: $sides{$side} accessions hit; hit filtering fates:";
 for (sort {$reasons{$b} <=> $reasons{$a}} keys %reasons) {print OUT " '$_'=$reasons{$_};"} print OUT "\n";
 for my $isle (sort {keys(%{$outs{$b}}) <=> keys(%{$outs{$a}})} keys %outs) {
  for (sort {$outs{$isle}{$b} <=> $outs{$isle}{$a}} keys %{$outs{$isle}}) {print OUT $_}
 }
 %reasons = ();
}
close OUT;
close FATE;

sub End {my ($entry, $dir) = @_; return 1 if $dir == -1; return $gnmlens{$entry}}
sub RunCommand {
 my ($command, $checkfile) = @_;
 if ($checkfile and -e $checkfile) {print "Skipping command: $command\n"; return}
 print "Running command: $command\n";
 my $out = system($command);
 if ($out) {print "Command '$command' failed with error message $out\n"; exit}
 else {print "Command '$command' succeeded\n"}
}

sub Btwn {my ($L, $R, $M) = @_; ($L,$R) = ($R,$L) if $L>$R; return 1 if $M>=$L and $M<=$R; return 0;}
sub getSeq {
 my ($db, $entry, $L, $R, $strand, $crossover) = @_;
 if ($L > $R) {($L, $R) = ($R, $L)}
 #print STDERR "Running blastdbcmd -db $db -dbtype nucl -entry '$entry' -range $L-$R -strand $strand -outfmt %s\n\n";
 my $seq = `blastdbcmd -db $db -dbtype nucl -entry '$entry' -range $L-$R -strand $strand -outfmt %s`;  # Always upper case
 print "No sequence obtained for $db entry $entry: $L-$R $strand\n" unless $seq;
 chomp $seq;
 $seq =~ tr/YRSWKMBDHV/NNNNNNNNNN/;
 my $begin = uc(substr($seq , 0 , length($seq)));
 my $middle = '';
 my $end = '';
 if (50 <= length($seq)) { $begin = uc(substr($seq, 0, 50)); }
 if ($crossover and 50 <= length($seq)) { $middle = lc(substr($seq, 50, length($seq) - 50)); }
 if ($crossover and 50+$crossover <= length($seq)) { $middle = lc(substr($seq, 50, $crossover)); }
 if (50+$crossover <= length($seq)) { $end = uc(substr($seq, 50+$crossover, length($seq)-(50+$crossover))); }
 if (100+$crossover <= length($seq)) { $end = uc(substr($seq, 50+$crossover, 50)) }
 my $markseq = $begin . $middle . $end . "\n";
 return $markseq;
}

sub Fate {my $reason = $_[0]; $reasons{$reason} ++; die "$reason\n" unless $q1line; return "$q1line\t$reason\n";}
