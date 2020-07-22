#! /usr/bin/perl
use strict; use warnings;
use File::Spec;
#ToDo: refannot; enable search for group I&II introns using low minsize (150 bp?); annotate via uninterrupted; annotate broken target pfam

# OPTIONS
my $invocation = "Called \"$0 @ARGV\" on " . localtime . "\n";
my $dir = File::Spec->rel2abs($0); $dir =~ s/\/([^\/]+)$//; my $scriptname = $1;
my ($verbose, $nickname, $force) = ('', 'NoNick', '');
my ($outfolder, $db, $prefix, $inDna, $complete, %cts);
my $cpu = 1;
my $tax = 'B';
my $qlenIsland = 15000;
my $gencode = 11;
my $circularity = 'C';
my $search = 'island';
Options(); # see bottom of script; help and other messages, Getopt
print $invocation if $verbose;
$inDna = File::Spec->rel2abs($inDna);
if ($outfolder) {$prefix = $inDna; $prefix =~ s/.*\///} else {$outfolder = $inDna; $outfolder =~ s/([^\/]+)$//; $prefix = $1}
$prefix =~ s/\.[^\.]+$//;
my %qlen = (island => $qlenIsland, IS => 3000);

# PROGRAM
my %aalookup = qw/Ala A Arg R Asn N Asp D Cys C Glu E Gln Q Gly G His H Ile I Ile2 J Leu L Lys K Met M Phe F Pro P Ser S Thr T Trp W Tyr Y Val V Pyl O SeC U tmRNA Z iMet B fMet B Sup X Undet ?/;
my (%finals, %mobs, %genes, $ptnCt, $hypothCt, %sizes, %hskpEnrich, %fornEnrich, %is, %stats, %dnaRA, $taxonomy, %is607);
for (qw/IS1535 IS1536 IS1537 IS1538 IS1539 IS1602 IS1921 IS607 ISAfe10 ISAfe11 ISArma1 ISBce17 ISBlo12 ISC1904 ISC1913 ISC1926 ISCaje1 ISCaje3 ISCARN1 ISCARN56 ISCbe1 ISCbo10 ISCbo6 ISCbo9 ISCbt3 ISCbt4 ISCfe1 ISChh1 ISDka1 ISLhe60 ISLhe9 ISMae17 ISMae20 ISNma20 ISPfu4 ISSis5 ISSis7 ISSoc2 ISSto11 ISSto12 ISSto13 ISTko1 ISTko2 ISTsi1 ISTvo1 ISvAR158_1 ISvMimi_1 ISvMimi_2 ISvNY2A_1 ISvNY2A_2 ISvPBCV_1/)
{$is607{$_} ++}
chdir $outfolder;
StripHeaders() unless -e "$prefix.fa";
my ($extra) = (''); $extra = "-nickname $nickname" if $nickname;
ReadTax();
#open JOBS, ">>/home/brilau/jobs/IScore.jobs"; # For making core jobs list
RunCommand("perl $dir/tater.pl -tax $tax -gencode $gencode -circle $circularity $extra $verbose $force $prefix", "$prefix.gff"); # Prepare gff with tRNAs, integrases
if ($complete) {$extra = "-complete "} else {$extra = ' ';}
RunCommand("perl $dir/dnaStats.pl -circle $circularity $extra$force $verbose $prefix", "$prefix.stats"); # take stats and prepare DNA subfolders
ReadStats();
ReadGff();
print "$cts{mob} mobs\n" if $verbose;
Load_housekeep_enrich();
if ($search eq "IS") {IS()} else {Island()} # calls RunComp(tigercore.pl, ConvertEnds, MobGenes, Annotate, Crossover)
PrintGff($search); # calls MeanSd, Scores(DeltaInt, HousekeepIndex, HypothIndex, Bias, Foreignness), ResolveOverlaps(intersectBed)
RunCommand("perl $dir/bin/merge.pl $prefix.$search.gff $search $prefix $nickname", "$prefix.$search.nonoverlap.gff");
#close JOBS; # For making core jobs list

# SUBROUTINES
sub StripHeaders {
 open OUT, ">$prefix.fa";
 for (@{ReadFile("$inDna")}) {s/^>[a-z]+\|[0-9]+\|[a-z]+\|([^\|]+)\|\S*/>$1/; print OUT $_ . "\n"}
 close OUT;
}

sub RunCommand {
 my ($command, $checkfile) = @_;
 if (not $force and ($checkfile and -e $checkfile)) {print "Skipping command (product $checkfile exists): $command\n" if $verbose; return}
 print "Running command: $command\n" if $verbose;
 my $out = system($command);
 if ($out) {die "Command '$command' failed with error message $out\n"}
 else {print "Command '$command' succeeded\n"}
}

sub ReadFile {
 my @ret;
 print "Reading file $_[0]\n" if $verbose;
 open IN, $_[0] or die "Can't open $_[0]\n";
 while (<IN>) {last if /^##FASTA/; next if /^#/; push @ret, $_}
 close IN; chomp @ret; return \@ret;
}

sub ReadTax {
 $taxonomy = "division=;phylum=;order=;class=;family=;genus=;species=;org=;taxid=gencode=$gencode;";
 return unless -f "$prefix.tax";
 for (`cat $prefix.tax`) {
  chomp; my ($taxid, $org, $rank, $code, $nick) = split "\t";
  $gencode = $code;
  my ($div, $phy, $ord, $cla, $fam, $gen, $spp) = split ';', $rank;
  $taxonomy = "division=$div;phylum=$phy;order=$ord;class=$cla;family=$fam;genus=$gen;species=$spp;org=$org;taxid=$taxid;gencode=$gencode;";
  $tax = 'B'; $tax = 'A' if $div eq 'Archaea'; $tax = 'M' if $gencode == 4; $tax = 'G' if $gencode == 25;
  $nickname = $nick if $nick;
 }
}

sub ReadStats {
 for (@{ReadFile("$prefix.stats")}) {
  next unless s/^(\S+)\t//;
  my $dna = $1;
  for my $cat (qw/len circ replicon cds pfam hypoth hskp forn/) {$stats{$dna}{$cat} = $1 if s/^(\S+)\t//}
  @{$dnaRA{$dna}} = ($dna, split "\t"); # relative abundances of C and key dinucleotides
 }
}

sub ReadGff {
 for (@{ReadFile("$prefix.gff")}) {
  my @f = split "\t";
  my ($dna, $type, $L, $R, $orient, $info) = @f[0,2,3,4,6,8];
  $orient = '+' if $orient eq '.';
  my $id = ''; $id = $1 if $f[8] =~ /ID=(.*?);/;
  $id = $1 if /annot=([^;]+)/;
  push @{$genes{$dna}}, {type => $type, L => $L, R => $R, id => $id, info => $info, orient => $orient, line => $_};
  $genes{$dna}[-1]{questionable} ++ if /questionable=[^;]/; #/\ttRNA\t.*(und|pseudogene|not called by both|low cove score|possibly problematic gene)/i;
  if ($type eq 'tmRNA') {$genes{$dna}[-1]{product} = $type} elsif ($type =~ /^tRNA/ and /aa=([^;]+).*ac=([^;]+)/) {$genes{$dna}[-1]{product} = $type . "-$1($2)"} #; die "$genes{$dna}[-1]{product}\n"}
  next unless $type eq "CDS";
  my $call = 'hypothetical'; if (/pfam1=([^;]+)/) {$call = $1}
  $genes{$dna}[-1]{pfam} = $call;
  next unless /annot=([^;]+)/;
  $cts{mob} ++;
  my $mob = $1;
  $genes{$dna}[-1]{annot} = $mob;
  my $mid = int(($L+$R)/2);
  push @{$mobs{$dna}}, {dna => $dna, L => $L, R => $R, orient => $orient, id => $mob, mid => $mid};
 }
}

sub Load_housekeep_enrich {
 for (@{ReadFile("$dir/db/housekeep_enrich.txt")}) {$hskpEnrich{$1} = $2 if /(\S+)\t(\S+)/}
 for (@{ReadFile("$dir/db/foreign_enrich.txt")})   {$fornEnrich{$1} = $2 if /(\S+)\t(\S+)/}
}

sub Island { # run TIGER to find islands
 print "Island mode initiated\n" if $verbose;
 my ($islemin, $islemax, $islemaxover, $islemobtype) = (2000, 200000, 250, "island");
 for my $dna (keys %mobs) {
  for my $mobile (@{$mobs{$dna}}) {
   next unless $$mobile{id} =~ /Y-Int|S-Int|S-Core/;
   print "Island test for $$mobile{id}\n" if $verbose;
   RunComp($islemobtype, $mobile, $islemin, $islemax, $islemaxover, $dna);
  }
 }
 my ($ISmin, $ISmax, $ISmaxover, $ISmobtype) = (500, 15000, 30, "IS");
 for my $dna (keys %mobs) {
  for my $IS (@{$mobs{$dna}}) {
   next unless $$IS{IS};
   print "IS test for $$IS{IS}\n" if $verbose;
   RunComp($ISmobtype, $IS, $ISmin, $ISmax, $ISmaxover, $dna) if $$IS{IS};
  }
 }
 PrintGff("ISartifact");
 RunCommand("perl $dir/bin/merge.pl $prefix.ISartifact.gff ISartifact $prefix $nickname", "$prefix.ISartifact.nonoverlap.gff");
 open IN, "$prefix.ISartifact.nonoverlap.gff";
 while (<IN>) { # Island / IS relationships
  my @flds = split /\t/;
  $flds[8] =~ /;int=([^:;]+)/; 
  my $IS = $1; # TransposaseID
  for (keys %finals) { # Both islands and ISartifacts are in final
   my $final = $finals{$_};
   next unless $flds[0] eq $$final{dna};
   if (($$final{L} > $flds[3] and $$final{L} < $flds[4]) or ($$final{R} > $flds[3] and $$final{R} < $flds[4])) { # One end of final is within IS
    if ($$final{ISoverlap} =~ $IS) {$$final{transposon} = "$IS"}  # Idea is that twice-overlapped by same IS type has same IS at both ends, like a transposon
    $$final{ISoverlap} .= "$flds[3]-$flds[4].score:$flds[5].name:$IS,";
   }
   if (abs($flds[3]-$$final{L}) < 20 and abs($flds[4]-$$final{R}) < 20) {
    $$final{ISidentical} .= "$flds[3]-$flds[4].score:$flds[5].name:$IS,";
   }
   next unless $$final{type} eq "island"; # Reject island when q1 is entirely within IS
   if ($$final{mobQ1} =~ /$IS/ and $$final{q1convert} =~ /(\d+)-(\d+)/) {
    my ($qL, $qR)= ($1, $2); if ($qL > $qR) {$qL = $2; $qR = $1}
    if ($flds[3] - 20 < $qL and $flds[4] + 20 > $qR) {$$final{IS} .= "q1:$flds[3]-$flds[4].score:$flds[5].name:$IS,"}
   } 
   if ($$final{mobQ2} =~ /$IS/ and $$final{q2} =~ /\S+\>(\d+)\-(\d+)/) {
    my ($qL, $qR)= ($1, $2); if ($qL > $qR) {$qL = $2; $qR = $1}
    if ($flds[3] - 20 < $qL and $flds[4] + 20 > $qR) {$$final{IS} .= "q2:$flds[3]-$flds[4].score:$flds[5].name:$IS,"}
   }
  }
 }
 close IN;
 RunCommand("blastn -db $dir/db/is -query $prefix.fa -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen' > $prefix.IS.blast", "$prefix.IS.blast");
 open IN, "$prefix.IS.blast";
 while (<IN>) { # Same treatment as above for ISFinder matches
  my @flds = split /\t/;
  next if $is607{$flds[1]};
  for (keys %finals) {
   my $final = $finals{$_};
   next unless $flds[0] eq $$final{dna};
   next unless $flds[3] > 100;
   my ($sstart, $send) = ($flds[8], $flds[9]);
   if ($flds[8] > $flds[9]) {$sstart = $flds[9]; $send = $flds[8]}
   if (($$final{L} > $flds[6] and $$final{L} < $flds[7]) or ($$final{R} > $flds[6] and $$final{R} < $flds[7])) { # IS overlap either side?
    if ($$final{ISoverlap} =~ $flds[1]) {
     $$final{transposon} = "$flds[1]";
     $$final{ISoverlap} =~ /score:(\d+)-(\d+)/;
     unless ($1 > $send -20 or $sstart > $2 - 20) {$$final{transposon} .= ".match"}
    }
    $$final{ISoverlap} .= "$flds[6]-$flds[7].score:$sstart-$send.name:$flds[1],";
   }
   if (abs($flds[6]-$$final{L}) < 20 and abs($flds[7]-$$final{R}) < 20) {
    $$final{ISidentical} .= "$flds[6]-$flds[7].score:$flds[11].name:$flds[1],";
   }
   $$final{q1convert} =~ /(\d+)-(\d+)/;
   my ($q1L, $q1R)= ($1, $2); if ($q1L > $q1R) {$q1L = $2; $q1R = $1}
   if ($flds[6] - 20 < $q1L and $flds[7] + 20 > $q1R) {$$final{IS} .= "q1:$flds[6]-$flds[7].score:$flds[11].name:$flds[1],"}
   $$final{q2} =~ /\S+\>(\d+)\-(\d+)/;
   my ($q2L, $q2R)= ($1, $2); if ($q2L > $q2R) {$q2L = $2; $q2R = $1}
   if ($flds[6] - 20 < $q2L and $flds[7] + 20 > $q2R) {$$final{IS} .= "q2:$flds[6]-$flds[7].score:$flds[11].name:$flds[1],"} 
  } 
 }
}

sub IS { # run TIGER to find IS
 my ($ISmin, $ISmax, $maxover, $mobtype) = (500, 15000, 30, "IS");
 for my $dna (keys %mobs) {
  for my $mobile (@{$mobs{$dna}}) {
   RunComp($mobtype, $mobile, $ISmin, $ISmax, $maxover, $dna) if $$mobile{id} =~ /Tnp/;
  }
 }
}

sub RunComp { # run tigercore.pl script and parse uninterrupted.fa
 my ($mobtype, $mobile_ref, $min, $max, $over, $dna) = @_;
 #next unless $dna eq 'NC_011416.1';
 my %mobile = %$mobile_ref;
 #print "$mobtype $mobile{id}: $dna.$mobile{mid}\n" if $verbose;
 if ($mobile{orient} eq '+') {$mobile{id}  .= ":$mobile{L}-$mobile{R}"} else {$mobile{id} .= ":$mobile{R}-$mobile{L}"}
 my $cmd = "perl $dir/bin/tigercore.pl $prefix.fa $dna $mobile{mid} $db $min $max $over $qlen{$mobtype} 3000 $mobtype";
 #print JOBS "cd $prefix; $cmd >& $dna.$mobile{mid}.IS/core.log;\n"; next; # For making core jobs list
 my $subdir = "$dna.$mobile{mid}.$mobtype";
 unless (-d $subdir and -f "$subdir/uninterrupteds.txt") {print "Running: $cmd\n" if $verbose; system $cmd}
 open IN, "$subdir/uninterrupteds.txt" or die "No output file $subdir/uninterrupteds.txt\n";
 my (%hits, $q1, @align, $q2, @side); 
 while ($q1 = <IN>) {
  if ($q1 =~ /^# (.) side of coordinate: (\S+)/) {@side = ($1, $2); $side[1] = 0 if $side[1] eq 'no'; next}
  for (0..3) {$align[$_] = <IN>; chomp $align[$_]}
  $q2 = <IN>;
  die "parse $subdir q1 $q1" unless $q1 =~ s/^>(\S+);\s//;
  my $mob = "$mobile{id};$1";
  #$mob =~ /^\S+;(\d+):(\d+)-(\d+)$/;  next unless abs($2-$3) == 4253;
  die "parse q2 $q2" unless $q2 =~ s/^crossover=(\d+); bitsum=([0-9\.]+); hit=//;
  my ($crossover, $bitsum) = ($1, $2);
  if ($finals{$mob} and $finals{$mob}{bitsum} >= $bitsum) {
   $finals{$mob}{ct}++;
   $finals{$mob}{bitlist} .= ",$bitsum";
   next;
  }
  die "parse isle $mob\n" unless $mob =~ /^\S+;(\d+):(\d+)-(\d+)$/;
  my ($len, $L, $R, $orient, $refdir) = ($1, $2, $3, '+', 1);
  ($L, $R, $orient, $refdir) = ($R, $L, '-', -1) if $L > $R;
  my $ct = 1;
  $ct += $finals{$mob}{ct} if $finals{$mob}{ct};
  my (      $pctid1, $qS1, $qE1, $sS1, $sE1) = (split "\t", $q1)[2,6,7,8,9];
  my ($gnm, $pctid2, $qS2, $qE2, $sS2, $sE2) = (split "\t", $q2)[0,2,6,7,8,9];
  my ($qGnmL, $qGnmR) = ConvertEnds($qS1, $qE1, $mobile{mid}, $side[0], $mobtype, $stats{$dna}{len});
  #warn "$qGnmL, $qGnmR\n";
  my ($mobQ1, $mobQ2) = ('','');
  $mobQ1 = MobGenes($dna, $qGnmL, $qGnmR, $mobtype);
  $mobQ2 = MobGenes($dna, $sS2, $sE2, $mobtype);
  my $ints = IntsWithin($mobtype, $L, $R, $dna, -1);
  die "$subdir\n" unless defined $sE1;
  my ($OLL, $OLR, $ORL, $ORR, $OUL, $OUR, $seqL, $seqU, $seqR, $contextsum, $reorient, $context, $prefCoords) =
   Annotate($qGnmL, $qGnmR, $sS2, $sE2, ($sE1-(1+$refdir*$crossover)), $sE1, @align[1..3], $crossover, $mobile{mid}, $orient, $dna, $ints); 
  %{$finals{$mob}} = (len => $len, L => $L, R => $R, orient => $reorient, ct => $ct, gnm => $gnm, prefCoords => $prefCoords,
   bitsum => $bitsum, bitlist => $bitsum, crossover => $crossover, int => $mobile{id} , mid => $mobile{mid}, dna => $dna, side => $side[0].$side[1],
   q1 => "$pctid1:$qS1-$qE1($qGnmL-$qGnmR)>$sS1-$sE1", q2 => "$pctid2:$qS2-$qE2>$sS2-$sE2", OL => "$OLL-$OLR", OR => "$ORL-$ORR", OU => "$OUL-$OUR",
   isleLseq => $seqL, unintSeq => $seqU, isleRseq => $seqR, mobQ1 => $mobQ1, mobQ2 => $mobQ2, type => $mobtype, IS => '', , ISoverlap => '', ISidentical => '', transposon => '', ints => $ints,
   context => $context, contextsum => $contextsum, q1convert => "$qGnmL-$qGnmR", q1identity => $pctid1, q2identity => $pctid2, origOrient => $orient);
 }
 close IN;
}

sub PrintGff {
 my $mobtype = shift;
 open OUT, ">$prefix.$mobtype.gff" or die "Can't write $prefix.$mobtype.gff\n";
 open REJ, ">$prefix.$mobtype.rejects.gff" or die "Can't write $prefix.$mobtype.rejects.gff\n";
 if ($mobtype eq "ISartifact") {$mobtype = "IS"}
 for (sort {$finals{$a}{dna} cmp $finals{$b}{dna} || $finals{$a}{L} <=> $finals{$b}{L} || $finals{$a}{R} <=> $finals{$b}{R} } keys %finals) {
  my $final = $finals{$_};
  if ($$final{type} eq $mobtype) {
   my @n = split(/,/, $$final{bitlist});
   my ($mean, $sd) = MeanSd(@n);
   my ($delta_int, $forn, $hskp, $hypoth, $delta_GC, $dinuc, $overall) = Scores($$final{ints}, $$final{L}, $$final{R}, $$final{dna}, -1, $$final{len});
   my $out = join("\t", $$final{dna}, 'TIGER', $$final{type}, $$final{L}, $$final{R}, $$final{ct}, $$final{orient}, '.', '');
   $$final{brief} = sprintf('%.0f', $$final{len}/1000) . '.' . $$final{contextsum};
   $$final{brief} =~ s/\|.{2,3}\|/\|/;
   $$final{contextsum} =~ s/\|//g;
   for (qw/brief len contextsum prefCoords bitsum gnm q1 q2 crossover int mid side OL OR OU mobQ1 mobQ2 IS ISoverlap transposon ISidentical context origOrient q1identity q2identity isleLseq unintSeq isleRseq/)
   {die "$_ $$final{L}, $$final{R}\n" unless defined $$final{$_}; $out .= "$_=$$final{$_};"}
   $out .= "mean=$mean;SD=$sd;deltaint=$delta_int;foreign=$forn;housekeep=$hskp;hypoth=$hypoth;delta_GC=$delta_GC;dinuc=$dinuc;FPscore=$overall;project=$prefix;$taxonomy"
    . "replicon=$stats{$$final{dna}}{replicon};qlen=$qlenIsland;refannot=;ints=" . join(',', keys %{$$final{ints}});
   if ($$final{ISidentical} =~ /.+name:IS/) {print(REJ "$out;reject=identical;\n"); next} # if $$final{ISidentical} =~ /.+name:IS/; # Removes any islands with nearly identical (20bp leniency) endpoints as an IS (only from ISFinder)
   elsif ($$final{IS} =~ /.+/) {print(REJ "$out;reject=IS;\n"); next} #, next if $$final{IS} =~ /.+/; # Removes any IS artifacts, in which a q1/q2 blast hit was to a known IS (either from TIGER or ISFinder)
   elsif ($$final{crossover} > 2000) {print(REJ "$out;reject=crossover;\n"); next} # if $$final{crossover} > 2000; # Removes islands with crossovers near 2251, our new too-long "max" in TIGERCore
   elsif ($$final{transposon} =~ /match/) {print(REJ "$out;reject=transposon;\n"); next} # if $$final{transposon} =~ /match/; # Removes islands flagged as transposons 
   print OUT "$out;reject=;\n";
  }
 }
 close OUT;
 close REJ;
}

sub MeanSd {
 my ($mean, $stdev, $n, @vals) = (0, 0, scalar(@_), @_);
 return (0, 0) if $n == 0;
 for (@vals) {$mean += $_}               $mean /= $n;
 for (@vals) {$stdev += ($mean - $_)**2} $stdev = ($stdev/$n)**0.5;
 return ($mean, $stdev);
}

sub ConvertEnds { # Genomic coordinates (L<R) for hit portion of q1
 my ($L, $R, $mid, $side, $mobtype, $dnalen) = @_;
 #warn "$L, $R, $mid, $side, $mobtype, $dnalen\n";
 my ($genomeL, $genomeR, $qlenActual) = ($mid, $mid, $qlen{$mobtype});
 if ($side eq 'L') {
  if ($mid < $qlenActual) {$qlenActual = $mid}
  $genomeL += $L - $qlenActual;
  $genomeR += $R - $qlenActual;
 } else {
  if ($mid+$qlenActual > $dnalen) {$qlenActual = $dnalen -$mid}
  $genomeL += $qlenActual - $R;
  $genomeR += $qlenActual - $L;
 }
 return ($genomeL, $genomeR);
}

sub Annotate {
 my ($q1L, $q1R, $q2L, $q2R, $OUL, $OUR, $seqL, $seqU, $seqR, $crossover, $mob, $orient, $dna, $ints) = @_;
 #warn "$q1L, $q1R, $q2L, $q2R, $OUL, $OUR, $seqL, $seqU, $seqR, $crossover, $mob, $orient, $dna, $ints\n";
 my ($OLR, $ORL) = ($q1R, $q2L); # Proximal coordinates
    ($OLR, $ORL) = ($q2L, $q1L) if $orient eq '-';  # q2L is always proximal
 my $OLL = $OLR - $crossover + 1; # Distal
 my $ORR = $ORL + $crossover - 1; # OLL,OLR,ORL,ORR are in genomic order for now, may flip later based on target gene orientation
 my $flipseq; if ($orient eq '-') {$flipseq = 1} # Sequences, and OUL/R are still ordered/oriented by int-att side 
 my ($highR, $Lig, $Rig, %overlaps) = (-1, -1, -1); # Lig and Rig are gff lines for left and right nonoverlapping (in case intergene search), %overlaps are gff lines overlapping OL and OR
 foreach my $i (0 .. $#{$genes{$dna}}) {
  next if $genes{$dna}[$i]{questionable}; # Skip questionable (island-split fragment) tRNAs
  next if $genes{$dna}[$i]{annot} and $$ints{$genes{$dna}[$i]{annot}}; # Skip int within island
  if ($genes{$dna}[$i]{R} < $OLL and $genes{$dna}[$i]{R} > $highR) {($highR, $Lig) = ($genes{$dna}[$i]{R}, $i)} # Find Lig: gene with closest R end to left of OL and not overlapping OL
  $overlaps{L}[0] = $i unless $genes{$dna}[$i]{L} >= $OLR or $genes{$dna}[$i]{R} <= $OLL or $overlaps{L}[0]; # Genes overlapping OL
  $overlaps{R}[0] = $i unless $genes{$dna}[$i]{L} >= $ORR or $genes{$dna}[$i]{R} <= $ORL; # Genes overlapping OR
  if ($genes{$dna}[$i]{L} > $ORR) {$Rig = $i; last} # Find Rig: gene with closest L end to right of OR and not overlapping OR
 }
 if ($Lig > -1) {unshift @{$overlaps{L}}, $Lig} # Process Lig together with OL overlaps
 if ($Rig > -1) {push    @{$overlaps{R}}, $Rig}
 my ($nonintergene, @calls, @context);
 for my $side (qw/L R/) {
  next unless $overlaps{$side};
  for my $i (@{$overlaps{$side}}) { # $i corresponds to element of full %genes{$dna} gff array
   my ($pfam, $gene, $product, $annot, $dupl, $len, $intergene, $orient) = ('noPfam', 'noGene', 'noProduct', '', '3prime', $OLL - $genes{$dna}[$i]{L}, '', $genes{$dna}[$i]{orient});
   $pfam = $1 if $genes{$dna}[$i]{info} =~ /pfam1=([^;]+)/;
   $nonintergene ++ unless $i == $Lig or $i == $Rig or $genes{$dna}[$i]{info} =~ /annot=/; 
   if ($genes{$dna}[$i]{info} =~ /annot=([^;]+)/) {$annot = $1; $gene = $1}
   elsif ($genes{$dna}[$i]{info} =~ /gene=([^;]+)/) {$gene = $1}
   if ($genes{$dna}[$i]{info} =~ /product=([^;]+)/) {$product = $1}
   if (not $product and $genes{$dna}[$i]{info} =~ /rpt_family=([^;]+)/) {$product = $1}
   if ($product =~ /([^;]+)/) {$product =~ s/\s+/_/g; $product =~ s/[^a-zA-Z0-9_\-\.\(\)]//g}
   $dupl = '5prime' if ($side eq 'L' and $genes{$dna}[$i]{orient} eq '-') or ($side eq 'R' and $genes{$dna}[$i]{orient} eq '+');
   $len = $genes{$dna}[$i]{R} - $ORR if $side eq 'R';
   $intergene = 'intergene' if $i == $Lig or $i == $Rig;
   push @calls, {name => ContextName($pfam, $gene, $product, $annot), side => $side.$intergene, dupl => $dupl, len => $len, orient => $orient};
   $gene    = '' if $gene eq 'noGene';
   $pfam    = '' if $pfam eq 'noPfam';
   $product = '' if $gene or $pfam; # Avoid long product name if gene or pfam available
   push @context, [$pfam, $gene, $product, $side.$intergene, $dupl, $len];
  }
 }
 if ($nonintergene) {for (my $i = $#context; $i >=0; $i--) {splice(@context, $i, 1) if $context[$i] =~ /intergene/}}
 my (@parts, $relation, $flip);
 for (@calls) {
  next if $nonintergene and $$_{side} =~ /intergene/; 
  next if not $nonintergene and $$_{side} !~ /intergene/;
  push @parts, $_;
 }
 my @out = ($OLL, $OLR, $ORL, $ORR, $OUL, $OUR, $seqL, $seqU, $seqR, $parts[0]{name}, '+');
 # Handle simple cases of single gene overlap and same gene overlapping both ends (which is reduced to former)
 if ($parts[1] and $parts[0]{name} eq $parts[1]{name} and $parts[0]{orient} eq $parts[1]{orient}) { # Reduce same-name pair to single
  if ($parts[0]{len} < $parts[1]{len}) {splice @parts, 0, 1} else {splice @parts, 1, 1} # Leave longest gene portion as "main" gene
 }
 if (scalar(@parts) == 1 and $nonintergene) {$flip ++ if $parts[0]{orient} eq '-'}
 else { # Handle two-gene relationships
  for (1..2) {push @parts, {name => '', orient => 0} if scalar(@parts) < 2} # Raise number of entries in @parts to 2
  for (@parts) {$$_{orient} =~ tr/+-/></}
  $relation = $parts[0]{orient}.$parts[1]{orient};
  # Flip unidirectional if backward; Alphabetize (con|di)vergent genes; Four remaining relation cases 00|>0|0>|>> are already OK
  if    ($relation =~ /<0|0<|<</) {$flip ++; $relation = reverse $relation; $relation =~ tr/</>/}
  elsif ($relation =~ /<>|></ and ($parts[0]{name} cmp $parts[1]{name}) == 1) {$flip ++}
  if (not $nonintergene) {substr $relation, 1, 0, '/'} # Mark intergene with a slash
 }
 if ($flip) {
  @parts = reverse @parts;
  @context = reverse @context; for (@context) {$$_[3] =~ tr/LR/RL/}
  @out[0..3,-1] = (@out[3,2,1,0], '-'); # Flip OL and OR, match sign of island to sign of target gene
  $flipseq --; # 
 }
 if ($flipseq) {for (6..8) {$out[$_] = Revcomp($out[$_])} @out[4..8] = @out[5,4,8,7,6]}
 my @prefCoords; # target-distal side L or R within OL and OR 
 if (scalar @parts == 2) {
  $out[-2] = "$parts[0]{name}\|$relation\|$parts[1]{name}";
  unless ($relation =~ /<>|></) {
   my $testcontext = 0; $testcontext = 1 if $relation =~ /^0/;
   if    ($parts[$testcontext]{dupl} eq '5prime') {@prefCoords = ($out[1],$out[3])}
   elsif ($parts[$testcontext]{dupl} eq '3prime') {@prefCoords = ($out[0],$out[2])}
  }
 }
 elsif ($parts[0]{dupl} eq '5prime') {@prefCoords = ($out[1],$out[3])}
 elsif ($parts[0]{dupl} eq '3prime') {@prefCoords = ($out[0],$out[2])}
 @prefCoords = sort {$a <=> $b} @prefCoords;
 my @context2 = map {join '/', @{$_}} @context;
 return @out, join(',', @context2), join(',', @prefCoords);
}

sub ContextName { 
 my ($pfam, $gene, $product, $annot) = @_;
 if ($annot) {$annot =~ s/\.[0-9]+$//; return $annot}
 #die "$product $1\n" if $product =~ /^tRNA.*-(.+)\(/;
 #die "$product\n" if $product =~ /^tRNA.*\(/;
 #die "$product\n" if $product =~ /^tRNA/ and length($product) <12;
 return "$aalookup{$1}" if $product =~ /^tRNA.*-(.+)\(...\)$/;
 return 'Z' if $product eq 'tmRNA';
 unless ($gene eq 'noGene') {$gene =~ s/_\d+$//; return $gene}
 unless ($pfam eq 'noPfam') {return $pfam}
 return 'HYP' if $product =~ /^hypothetical_protein$/;
 unless ($product eq 'noProduct') {return $product}
 return 'unknown';
}

sub Revcomp {my $ret = reverse $_[0]; $ret =~ tr/ACGTacgt/TGCAtgca/; return $ret}

sub MobGenes {
 my ($dna, $S, $E, $mobtype) = @_;
 my $L = $S; my $R = $E; my $genes = '';
 if ($S > $E) {$L = $E; $R = $S}
 print "mobzyme test for dna $dna interval $L-$R\n" if $verbose;
 for my $mobile (@{$mobs{$dna}}) {
  next if ($$mobile{L} > $R or $$mobile{R} < $L);
  $genes .= "$$mobile{id},";
  if ($mobtype eq "island") {
    next unless $$mobile{id} =~ /^Tnp_/;
    $$mobile{IS} = $$mobile{id}; #load $$mobile{IS} with mobile id only if there is a transposase hit in Q1 or Q2 -> runs TIGER.pl on these
  }  
 }
 print "  Success\n" if $genes and $verbose;
 return $genes;
}

sub IntsWithin {
 my ($mobtype, $L, $R, $dna, $origin) = @_;
 my %ints;
 for my $int (@{$mobs{$dna}}) { # int must be over half within island
  next if $mobtype eq 'IS'     and $$int{id} !~ /^Tnp_/;
  next if $mobtype eq 'island' and $$int{id} !~ /Y-Int|S-Int|S-Core|Xer|Integron/;
  my $within = 0;
  my $midpt = $$int{mid};
  if ($origin == -1) {$within = 1 if $midpt > $L && $midpt < $R}
  else               {$within = 1 if $midpt > $L || $midpt < $R}
  #print "IntsWithin:$L $R $dna $origin $$int{id} $$int{mid} within=$within\n"; # 728586 738185 NC_000913.2 1 PhageIntegrase.1 295582 1
  next unless $within;
  %{$ints{$$int{id}}} = (L => $$int{L}, R => $$int{R});
 }
 return \%ints;
}

sub Scores {
 my ($ints, $L, $R, $dna, $origin, $len) = @_;
 #print "scoring $L, $R, $origin\n";
 my $delta_int = DeltaInt($ints, $L, $R, $dna, $origin, $len);
 my ($hypoth, $forn, $hskp) = Foreign($ints, $L, $R, $dna, $origin);
 my ($delta_GC, $dinuc) = Bias($L, $R, $dna, $origin);
 my $overall = exp(-35.2 + 6.4*log10($len) + 1.2*log10($delta_int) + 0.006*$forn + -0.39*$hskp + -5.3*$hypoth + 7.0*abs($delta_GC) + -25.7*$dinuc);
 return ($delta_int, $forn, $hskp, $hypoth, $delta_GC, $dinuc, $overall);
}

sub log10 {return log(shift)/log(10);}

sub DeltaInt { # Shortest distance between any internal integrase gene end and an island end
 my ($ints, $L, $R, $dna, $origin, $ret) = @_;
 $R += $stats{$dna}{len} if $origin == 1; # virtual right end for origin-spanning island
 for (keys %{$ints}) {
  my ($intL, $intR) = ($$ints{$_}{L}, $$ints{$_}{R});
  if ($origin == 1 and $intL < $R) { for ($intL, $intR) {$_ += $stats{$dna}{len}} }
  $ret = $intL-$L if $intL-$L < $ret;
  $ret = $R-$intR if $R-$intR < $ret;
 }
 $ret = 1 if $ret < 1;
 return $ret
}

sub Foreign {
 my ($ints, $L, $R, $dna, $origin) = @_;
 my $x=''; for (keys %{$ints}) {$x .= $_}  
 my ($cds, $hypoth, $pfam, $forn, $hskp) = (0,0,0,0,0);
 for my $prot (@{$genes{$dna}}) {
  next unless $$prot{line} =~ /\tCDS\t/;
  #die "$$prot{id}, $L, $R, $origin, $$prot{L}, $$prot{R}\n" unless $$prot{call};# $$ints{$$prot{id}}\n";
  #print "$$prot{id}, $L, $R, $origin, $$prot{L}, $$prot{R}, x='$x'\n";
  unless ($$ints{$$prot{id}}) {
   if ($origin == 1) {next if $$prot{L} <= $L and $$prot{R} >= $R}
   else {next unless $$prot{L} >= $L and $$prot{R} <= $R; last if $$prot{L} > $R}
  }
  $cds ++;
  if ($$prot{pfam} eq 'hypothetical') {$hypoth ++}
  else {
   $pfam ++;
   $forn += $fornEnrich{$$prot{pfam}} if $fornEnrich{$$prot{pfam}};
   $hskp += $hskpEnrich{$$prot{pfam}} if $hskpEnrich{$$prot{pfam}};
  }
 }
 $hypoth /= $cds if $cds;   $hypoth -= $stats{$dna}{hypoth}; # Difference between hypoth Densities of island and its replicon
 $forn   /= $pfam if $pfam; $forn   -= $stats{$dna}{forn};
 $hskp   /= $pfam if $pfam; $hskp   -= $stats{$dna}{hskp};
 return $hypoth, $forn, $hskp;
}

sub Bias {
 my ($L, $R, $dna, $origin) = @_;
 if ($origin == -1) {system "perl $dir/bin/collectSeq.pl -i $prefix.fa -e $dna -L $L -R $R > test.fa"}
 else {
  system "perl $dir/bin/collectSeq.pl -i $prefix.fa -e $dna -L $L -R $stats{$dna}{len} > test.fa";
  system "perl $dir/bin/collectSeq.pl -i $prefix.fa -e $dna -L 1 -R $R >> test.fa";
 }
 #die "perl $dir/bin/collectSeq.pl -i $prefix.fa -e $dna -L $L -R $R > test.fa\n$L, $R, $dna, $origin, $prefix\n" unless -s "test.fa";
 my $out = `perl $dir/bin/relAbun.pl test.fa`; chomp $out; $out =~ s/\n.*//s; my @testRA = split "\t", $out;
 my ($delta_GC, $di) = (($testRA[1]-$dnaRA{$dna}[1])/2, 0);
 for (2,3,4,6,7,9) {$di += abs($testRA[$_]-$dnaRA{$dna}[$_])}
 $di *= 2; # Double the 6 asymmetrical dinucs above to account for their complements, but don't double the 4 symmetrical dinucs below
 for (5,8,10,11) {$di += abs($testRA[$_]-$dnaRA{$dna}[$_])}
 return ($delta_GC, $di/16);
}

sub BinarySearchAoH { # Array of hashes $aref numerically presorted by $$aref[$_]{$aItem}
 my ( $target , $aref , $aItem ) = @_ ;
 my ( $low , $hi , $try , $result ) = ( 0 , $#$aref , '' , '' );
 while ( 1 ) {
  $try = int ( ( $low + $hi ) / 2 ) ; # Midpoint between hi and low
  $result = $aref->[$try]{$aItem} <=> $target ; # -1, 0, or 1
  if ($result < 0) { # Try is too low
   $low = $hi , next if $try == $low && $hi != $low; # required to reach ceiling
   return $try + 0.5 if $try == $low;
   $low = $try;
  } elsif ($result > 0) { # Try is too high (floor reachable by int round-down in $mid calculation)
   return $try - 0.5 if $try == $hi ;
   $hi = $try ;
  } else { return $try } # Try is correct
 }
}

sub Options {
my $version = '1.0 (Nov 2016)';
#   '    |    '    |    '    |    '    |    '    |    '    |    '    |    '    |
my $help = <<END;
$scriptname version $version
Usage: perl $scriptname [options] -db <RefDatabase> -fasta <GenomicDNA>
  -fasta:    Genomic fasta DNA sequence file.
  -db:       Blast database of reference genomes, absolute path.
  -search:     Search type. Specify island or IS. Default: island.
  -tax:      Taxonomic info for query genome. Enter name of a file containing 
              NCBI taxonomy string, or use B for Bacteria, A for Archaea, M for 
              Mycoplasmatales/Entomoplasmatales, G for Gracilibacteria/candidate
              division SR1. Automatically sets -gencode. Default: B.
  -gencode:  Genetic code table to use (see NCBI). Default: 11.
  -nickname: Brief name for genome (as might be used to start a locus_tag).
  -circle:   Specify C if all genomic DNA sequences are circular, L if all DNAs 
              are linear, or a filename for a tab-delimited file of query and 
              circularity (eg. acc.vers[tab]circular/linear). Default: C.
  -complete: Consider genome complete and categorize replicons. Default:
              consider genome incomplete and call all entries contigs.
  -qlen:     Query length for islands. Default: $qlenIsland (3000 is always used
              in the second pass test for IS's that rule out island artifacts).
  -force:    Overwrite current output files. Default: leave existing files.
  -outDir:   Output directory. Default: same directory as GENOME_FASTA_FILE.
  -cpu:      Number of cpus to use. Default: 1.
  Additional options: -help, -version, -verbose, -authors, -license

Example: perl $scriptname -verbose -db [path]/refseq_genomic -fasta Ecoli.fa

END
my $authors = "AUTHORS: Kelly Williams (kpwilli\@sandia.gov), Britney Lau, Julian Wagner\n";
my $license = "COPYRIGHT AND LICENSE: Copyright 2018, National Technology & Engineering Solutions of Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS, the U.S. Government retains certain rights in this software.
Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS 'AS IS' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.\n";

#Copyright (2013) Sandia Corporation. Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains certain rights in this software.\n";

die $help if @ARGV == 0;
use Getopt::Long;
my $options_okay = GetOptions(
 'help' => sub {print $help; exit},
 'version' => sub {print "$scriptname version $version\n"; exit},
 'authors' => sub {print $authors; exit},
 'license' => sub {print $license; exit},
 'verbose' => sub {$verbose = "-verbose"},
 'force' => sub {$force = "-force"},
 'outDir=s' => \$outfolder,
 'gencode=i' => \$gencode,
 'fasta=s' => \$inDna,
 'qlen=i' => \$qlenIsland,
 'db=s' => \$db,
 'complete' => \$complete,
 'tax=s' => \$tax,
 'nickname=s' => \$nickname,
 'cpu=i' => \$cpu,
 'circle=s' => \$circularity,
 'search=s' => \$search,
);
die $help if !$options_okay;
die "-db and -fasta are required\n\n$help" unless $inDna and $db;
}
