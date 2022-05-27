#! /usr/bin/perl
use strict; use warnings;
use File::Spec;
#ToDo: refannot; enable search for group I&II introns using low minsize (150 bp?); annotate via uninterrupted; annotate broken target pfam

# OPTIONS
my $invocation = "Called \"$0 @ARGV\" on " . localtime . "\n";
my $dir = File::Spec->rel2abs($0); $dir =~ s/\/([^\/]+)$//; my $scriptname = $1;
my ($verbose, $nickname, $force) = ('', 'NoNick', '');
my ($outfolder, $db, $prefix, $inDna, $complete, %cts);
my ($cpu, $tax, $qlenIsland, $gencode, $circularity, $search, $cross) = (1, 'B', 15000, 11, 'C', 'island', 'intact');  # Defaults
my ($j, %rankgenetype); for (qw/tmRNA tRNA gene annot pfam product rpt HYP unknown empty/) {$j ++; $rankgenetype{$_} = $j}
Options(); # see bottom of script; help and other messages, Getopt
if ($cross eq 'cross') {($circularity, $complete) = ('L', '')}
elsif ($cross eq 'circleOrigin') {($circularity, $complete) = ('C', 1)}
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
print STDERR "Running perl $dir/tater.pl -tax $tax -gencode $gencode -circle $circularity $extra $verbose $force $prefix\n";
RunCommand("perl $dir/tater.pl -tax $tax -gencode $gencode -circle $circularity $extra $verbose $force $prefix", "$prefix.gff"); # Prepare gff with tRNAs, integrases
if ($complete) {$extra = "-complete "} else {$extra = ' ';}
print STDERR "Running perl $dir/dnaStats.pl -circle $circularity $extra$force $verbose $prefix";
RunCommand("perl $dir/dnaStats.pl -circle $circularity $extra$force $verbose $prefix", "$prefix.stats"); # take stats and prepare DNA subfolders
ReadStats();
ReadGff();
print "$cts{mob} mobs\n" if $verbose;
Load_housekeep_enrich();
if ($search eq "IS") {IS()} else {Island()} # calls RunComp(tigercore.pl, ConvertEnds, MobGenes, Annotate, Crossover)
PrintGff($search); # calls MeanSd, Scores(DeltaInt, HousekeepIndex, HypothIndex, Bias, Foreignness), ResolveOverlaps(intersectBed)
print STDERR "perl $dir/merge.pl $prefix.$search.gff $search $prefix $nickname";
RunCommand("perl $dir/merge.pl $prefix.$search.gff $search $prefix $nickname", "$prefix.$search.nonoverlap.gff");

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
 for (@{ReadFile("$dir/../db/housekeep_enrich.txt")}) {$hskpEnrich{$1} = $2 if /(\S+)\t(\S+)/}
 for (@{ReadFile("$dir/../db/foreign_enrich.txt")})   {$fornEnrich{$1} = $2 if /(\S+)\t(\S+)/}
}

sub Island { # run TIGER to find islands
 print "Island mode initiated\n" if $verbose;
 my ($islemin, $islemax, $islemaxover, $islemobtype, $intsTried) = (2000, 200000, 250, "island");
 for my $dna (keys %mobs) {
  for my $mobile (@{$mobs{$dna}}) {
   next unless $$mobile{id} =~ /Y-Int|S-Int|S-Core/;
   $intsTried ++;
   print "Island test for $$mobile{id}\n" if $verbose;
   RunComp($islemobtype, $mobile, $islemin, $islemax, $islemaxover, $dna);
  }
 }
 unless ($intsTried) {print "No integrases found in genome\n" if $verbose; exit}
 my $ct = scalar(keys %finals);
 print "$ct calls after island mode\n" if $verbose;
 my ($ISmin, $ISmax, $ISmaxover, $ISmobtype) = (500, 15000, 30, "IS");
 for my $dna (keys %mobs) {
  for my $IS (@{$mobs{$dna}}) {
   next unless $$IS{IS};
   print "IS test for $$IS{IS}\n" if $verbose;
   RunComp($ISmobtype, $IS, $ISmin, $ISmax, $ISmaxover, $dna) if $$IS{IS};
  }
 }
 $ct = scalar(keys %finals); print "$ct calls after ISartifact mode\n";
 PrintGff("ISartifact");
 print STDERR "perl $dir/merge.pl $prefix.ISartifact.gff ISartifact $prefix $nickname";
 RunCommand("perl $dir/merge.pl $prefix.ISartifact.gff ISartifact $prefix $nickname", "$prefix.ISartifact.nonoverlap.gff");
 open IN, "$prefix.ISartifact.nonoverlap.gff";
 while (<IN>) { # Island / IS relationships
  my @f = split /\t/;
  $f[8] =~ /coord=([^;]+).*;int=([^:;]+)/; 
  my ($coord, $IS, @ends) = ($1, $2); # TransposaseID
  $coord =~ /^([^\/]+)\/(\d+)-(\d+)/;
  push @ends, [$1, sort {$a <=> $b} ($2, $3)];
  if ($coord =~ /\+([^\/]+)\/(\d+)-(\d+)/) {push @ends, [$1, sort {$a <=> $b} ($2, $3)]}
  for (keys %finals) { # Both islands and ISartifacts are in final
   my ($final, $e) = ($finals{$_}, $finals{$_}{ends});
   for my $j (0,1) {
    last unless $ends[$j];
    my ($dna, $L, $R) = @{$ends[$j]};
    if ($$final{mobQ1} =~ /$IS/ and $$final{q1} =~ /\(([^:]+):(\d+)-(\d+)\)/) {
     my ($qdna, $qL, $qR)= ($1, sort {$a <=> $b} ($2, $3));
     if ($dna eq $qdna and $L - 20 < $qL and $R + 20 > $qR) {$$final{IS} .= "q1:$L-$R.score:$f[5].name:$IS,"}
    } 
    if ($$final{mobQ2} =~ /$IS/ and $$final{q2} =~ /\S+\>([^:]+):(\d+)\-(\d+)/) {
     my ($qdna, $qL, $qR)= ($1, sort {$a <=> $b} ($2, $3));
     if ($dna eq $qdna and $L - 20 < $qL and $R + 20 > $qR) {$$final{IS} .= "q2:$L-$R.score:$f[5].name:$IS,"}
    }
    for my $i (0,1) {
     next unless $dna eq $$e[$i]{dna};
     if (($$e[$i]{L} > $L and $$e[$i]{L} < $R) or ($$e[$i]{R} > $L and $$e[$i]{R} < $R)) { # One end of final is within IS
      if ($$final{ISoverlap} =~ $IS) {$$final{transposon} = "$IS"}  # Idea is that twice-overlapped by same IS type has same IS at both ends, like a transposon
      $$final{ISoverlap} .= "$dna:$L-$R.score:$f[5].name:$IS,";
     }
     if (abs($L-$$e[$i]{L}) < 20 and abs($R-$$e[$i]{R}) < 20) {
      $$final{ISidentical} .= "$dna:$L-$R.score:$f[5].name:$IS,";
     }
     #next unless $$final{type} eq "island"; # Reject island when q1 is entirely within IS
     last if $$final{compose} eq 'simple';
    }
   }
  }
 }
 close IN;
 print STDERR "blastn -db $dir/../db/is -query $prefix.fa -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen";
 RunCommand("blastn -db $dir/../db/is -query $prefix.fa -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen' > $prefix.IS.blast", "$prefix.IS.blast");
 open IN, "$prefix.IS.blast";
 while (<IN>) { # Same treatment as above for ISFinder matches
  my @f = split /\t/;
  next if $is607{$f[1]};
  for (keys %finals) {
   my ($final, $e) = ($finals{$_}, $finals{$_}{ends});
   $$final{q1} =~ /\(([^:]+):(\d+)-(\d+)\)/;
   my ($qdna, $q1L, $q1R) = ($1, sort {$a <=> $b} ($2, $3));
   #print "$$final{context} $qdna $f[0] $f[6] $q1L $f[7] $q1R\n";
   if ($qdna eq $f[0] and $f[6] - 20 < $q1L and $f[7] + 20 > $q1R) {$$final{IS} .= "q1:$f[6]-$f[7].score:$f[11].name:$f[1],"}
   $$final{q2} =~ /\S+\>([^:]+):(\d+)\-(\d+)/;
   my ($qdna2, $q2L, $q2R) = ($1, sort {$a <=> $b} ($2, $3));
   #print "$qdna2 $f[0] $f[6] $q2L $f[7] $q2R\n";
   if ($qdna2 eq $f[0] and $f[6] - 20 < $q2L and $f[7] + 20 > $q2R) {$$final{IS} .= "q2:$f[6]-$f[7].score:$f[11].name:$f[1],"} 
   for my $i (0,1) {
    #print "$$e[$i]{dna}\n";
    next unless $f[0] eq $$e[$i]{dna};
    next unless $f[3] > 100;
    my ($sstart, $send) = ($f[8], $f[9]);
    if ($f[8] > $f[9]) {$sstart = $f[9]; $send = $f[8]}
    if (($$e[$i]{L} > $f[6] and $$e[$i]{L} < $f[7]) or ($$e[$i]{R} > $f[6] and $$e[$i]{R} < $f[7])) { # IS overlap either side?
     if ($$final{ISoverlap} =~ $f[1]) {
      $$final{transposon} = "$f[1]";
      $$final{ISoverlap} =~ /score:(\d+)-(\d+)/;
      unless ($1 > $send -20 or $sstart > $2 - 20) {$$final{transposon} .= ".match"}
     }
     $$final{ISoverlap} .= "$f[6]-$f[7].score:$sstart-$send.name:$f[1],";
    }
    if (abs($f[6]-$$e[$i]{L}) < 20 and abs($f[7]-$$e[$i]{R}) < 20) {
     $$final{ISidentical} .= "$f[6]-$f[7].score:$f[11].name:$f[1],";
    }
    last if $$final{compose} eq 'simple';
   } 
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
 my %mobile = %$mobile_ref;
 #print "$mobtype $mobile{id}: $dna.$mobile{mid}\n" if $verbose;
 if ($mobile{orient} eq '+') {$mobile{id}  .= ":$mobile{L}-$mobile{R}"} else {$mobile{id} .= ":$mobile{R}-$mobile{L}"}
 my $cmd = "perl $dir/tigercore.pl $prefix.fa $dna $mobile{mid} $db $min $max $over $qlen{$mobtype} 3000 $mobtype $cross";
 my $subdir = "$dna.$mobile{mid}.$mobtype";
 unless (-d $subdir and -f "$subdir/uninterrupteds.txt") {print "Running: $cmd\n" if $verbose; system $cmd}
 open IN, "$subdir/uninterrupteds.txt" or die "No output file $subdir/uninterrupteds.txt\n";
 my (%hits, $q1, @align, $q2, @side); 
 while ($q1 = <IN>) {
  if ($q1 =~ /^# (.) side of coordinate: (\S+)/) {@side = ($1, $2, 1); $side[1] = 0 if $side[1] eq 'no'; $side[2] = -1 if $1 eq 'R'; next}
  # >4115:CYEF01000004.1/1602-1+CYEF01000029.1/1-2515; VDLP01000003.1; bitsum=31179(26179+5000); length=4115; crossover=19; cross
  # proximal:CYEF01000004.1,1543-1661 reference:VDLP01000003.1,51826-51944 distal:CYEF01000029.1,2456-2574
  # Three sequence lines
  $q2 = <IN>;
  for (0..2) {$align[$_] = <IN>; chomp $align[$_]}
  die "parse $subdir q1 $q1" unless $q1 =~ s/^>(\S+);\s.*bitsum=([0-9\.]+).*crossover=(\d+); compose=(\S+)//;
  my ($mob, $bitsum, $crossover, $compos) = ("$mobile{id};$1", $2, $3, $4);
  #print "$mob\n" if $verbose;
  if ($finals{$mob} and $finals{$mob}{bitsum} >= $bitsum) {  # Already saw this island, with better bitsum
   $finals{$mob}{ct}++;
   $finals{$mob}{bitlist} .= ",$bitsum";
   next;
  }
  my $ct = 1; $ct += $finals{$mob}{ct} if $finals{$mob}{ct};
  my $bitlist = $bitsum; $bitlist .= ",$finals{$mob}{bitlist}" if $finals{$mob}{bitlist};
  die "$q2\n" unless $q2 =~ /proximal:[^,]+,([^,]+).*reference:([^,]+).*distal:([^,]+),([^,]+)/;
  my ($hitsum1, $gnm, $hitdna, $hitsum2) = ($1, $2, $3, $4);
  my ($pctid1, $oqS1, $oqE1, $sS1, $sE1) = (split '/', $hitsum1);
  my ($pctid2, $qS2,  $qE2,  $sS2, $sE2) = (split '/', $hitsum2);
  my ($qS1, $qE1) = ConvertEnds($oqS1, $oqE1, $mobile{mid}, $side[0], $mobtype, $stats{$dna}{len});
  # q1=100.000:1-2535(4499512-4502046)>5658-8192;q2=100.000:26-3000>4498181-4495207;
  my $q1out = "$pctid1:$oqS1-$oqE1($dna:$qS1-$qE1)>$sS1-$sE1";
  my $q2out = "$pctid2:$qS2-$qE2>$hitdna:$sS2-$sE2";
  die "$mob parse\n" unless $mob =~ /^\S+;(\d+):([^\/]+)\/(\d+)-(\d+)/;
  my ($len, @ends) = ($1);
  %{$ends[0]} = (dna => $2, S => $3, E => $4, dir => 1, OS => $qE1, seq => $align[0], compos => $compos);
  %{$ends[1]} = (dna => $2, E => $3, S => $4, dir => 1, OS => $sS2, seq => $align[2]);
  @{$ends[1]}{qw/dna E S/} = ($1, $2, $3) if $mob =~ /\+([^\/]+)\/(\d+)-(\d+)$/;
  for my $i (0..$#ends) {
   if ($ends[$i]{S} > $ends[$i]{E}) {$ends[$i]{dir} *= -1}
   $ends[$i]{OE} = $ends[$i]{OS} - $ends[$i]{dir} * ($crossover - 1);
   @{$ends[$i]}{qw/L R/}   = sort {$a <=> $b} @{$ends[$i]}{qw/S E/};
   @{$ends[$i]}{qw/OL OR/} = sort {$a <=> $b} @{$ends[$i]}{qw/OS OE/};
  }
  my $mobQ1 = MobGenes($ends[0]{dna}, $qS1, $qE1, $mobtype);
  my $mobQ2 = MobGenes($ends[1]{dna}, $sS2, $sE2, $mobtype);
  my $ints = IntsWithin($mobtype, \@ends);
  my %ref = (seq => $align[1], dna => $gnm, E => $sE1, dir => 1);
  $ref{dir} = -1 if $sS1 > $sE1;
  $ref{S} = $sE1-(1+$ref{dir}*$crossover);
  die "$subdir\n" unless defined $sE1;
  my ($reend, $flanksum, $overlap, $flip) = Annotate($crossover, $mobile{mid}, \@ends, $ints); 
  @ends = @{$reend};
  if ($flip) {$ref{dir} *= -1; @ref{qw/seq S E/} = (Revcomp($ref{seq}), @ref{qw/E S/})}
  my $coord = "$ends[0]{dna}/$ends[0]{S}-$ends[0]{E}";
  $coord  .= "+$ends[1]{dna}/$ends[1]{E}-$ends[1]{S}" unless $compos eq 'simple';
  #print "$mob => $coord, '$compos'\n" if $verbose;
  my ($dnanow, $L, $R, $dir, $OLL, $OLR, $seqL) = @{$ends[0]}{qw/dna L R dir OL OR seq/};
  my $orient = '+'; $orient = '-' if $dir < 0;
  my ($ORL, $ORR, $seqR) = @{$ends[1]}{qw/OL OR seq/};
  my ($OUL, $OUR, $seqU) = @ref{qw/S E seq/};
  %{$finals{$mob}} = (len => $len, L => $L, R => $R, orient => $orient, ct => $ct, gnm => $gnm, compose => $compos,
   coord => $coord, bitsum => $bitsum, bitlist => $bitlist, crossover => $crossover, int => $mobile{id}, flip => $flip,
   mid => $mobile{mid}, dna => $dnanow, side => $side[0].$side[1], OL => "$OLL-$OLR", OR => "$ORL-$ORR", OU => "$OUL-$OUR",
   isleLseq => $seqL, unintSeq => $seqU, isleRseq => $seqR, mobQ1 => $mobQ1, mobQ2 => $mobQ2, type => $mobtype,
   IS => '', , ISoverlap => '', ISidentical => '', transposon => '', ints => $ints, q1 => $q1out, q2 => $q2out, ends => \@ends,
   flanks => $flanksum, context => $overlap, q1convert => "$qS1-$qE1", q1identity => $pctid1, q2identity => $pctid2);
 }
 close IN;
}

sub PrintGff {
 my $mobtype = shift;
 print "Writing $prefix.$mobtype.gff\n" if $verbose;
 open OUT, ">$prefix.$mobtype.gff" or die "Can't write $prefix.$mobtype.gff\n";
 open REJ, ">$prefix.$mobtype.rejects.gff" or die "Can't write $prefix.$mobtype.rejects.gff\n";
 if ($mobtype eq "ISartifact") {$mobtype = "IS"}
 for (sort {$finals{$a}{dna} cmp $finals{$b}{dna} || $finals{$a}{L} <=> $finals{$b}{L} || $finals{$a}{R} <=> $finals{$b}{R} } keys %finals) {
  my $final = $finals{$_};
  if ($$final{type} eq $mobtype) {
   my @n = split(/,/, $$final{bitlist});
   my ($mean, $sd) = MeanSd(@n);
   my ($delta_int, $forn, $hskp, $hypoth, $delta_GC, $dinuc, $overall) = Scores($$final{ints}, $$final{coord}, -1, $$final{len});
   my $out = join("\t", $$final{dna}, 'TIGER', $$final{type}, $$final{L}, $$final{R}, $$final{ct}, $$final{orient}, '.', '');
   $$final{brief} = sprintf('%.0f', $$final{len}/1000) . '.' . $$final{context};
   $$final{brief} =~ s/\|.{2,3}\|/\|/;
   for (qw/brief coord compose len context flanks flip bitsum gnm crossover int mid side OL OR OU mobQ1 mobQ2 IS ISoverlap transposon ISidentical q1 q2 q1identity q2identity isleLseq unintSeq isleRseq/)
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

sub ConvertEnds { # Genomic coordinates (Start>End) for hit portion of q1
 my ($S, $E, $mid, $side, $mobtype, $dnalen) = @_;
 #warn "$L, $R, $mid, $side, $mobtype, $dnalen\n";
 my ($genomeS, $genomeE, $qlenActual) = ($mid, $mid, $qlen{$mobtype});
 if ($side eq 'L') {
  if ($mid < $qlenActual) {$qlenActual = $mid}
  $genomeS += $S - $qlenActual;
  $genomeE += $E - $qlenActual;
 } else {
  if ($mid+$qlenActual > $dnalen) {$qlenActual = $dnalen - $mid}
  $genomeS += $qlenActual - $S;
  $genomeE += $qlenActual - $E;
 }
 return ($genomeS, $genomeE);
}

sub Annotate {
 my ($crossover, $mob, $preends, $ints) = @_;
 my @e = @{$preends};  # $ends[0] is int-proximal att, seq aiming into island 
 my ($flip, @flanks, @olaps) = (0);
 for my $j (0,1) {  # Gene dir is adjusted to end dir then distal is flipped
  (@{$flanks[$j]}{qw/name dir rank/}, @{$olaps[$j]}{qw/name dir rank/}) = GeneEdge(@{$e[$j]}{qw/dna OL OR dir/}, $j, $ints);
  #print "$j @{$flanks[$j]}{qw/name dir rank/}, @{$olaps[$j]}{qw/name dir rank/}; @{$e[$j]}{qw/dna OL OR dir/}\n"; #if $e[0]{L} == 1211524;
 }
 $olaps[1]{dir} *= -1;
 #print "hi1 '$olaps[1]{name}' $flip $flanks[0]{name}$flanks[1]{name} $olaps[0]{dir} $olaps[1]{dir} $flanks[0]{dir} $flanks[1]{dir} $e[0]{dir} $e[1]{dir}\n"; # if $e$
 $flanks[1]{dir} *= -1;
 $e[1]{dir} *= -1; # Set distal gene dir to context of proximal
 @olaps = sort {$$a{rank} <=> $$b{rank} || $$a{name} cmp $$b{name}} @olaps;
 my $name = $olaps[0]{name};
 if ($name) {$flip ++ if $olaps[0]{dir} < 0} #!= $e[0]{dir}}
 else {
  my $pref = (sort {$flanks[$a]{rank} <=> $flanks[$b]{rank} || $flanks[$a]{name} cmp $flanks[$a]{name}} 0,1)[0];
  $flip ++ if $flanks[$pref]{dir} < 1; #* $e[$pref]{dir}< 1;
  #print "flanks $pref $flip $flanks[$pref]{dir}\n";
 }
 #print "hi '$name' $flip $flanks[0]{name}$flanks[1]{name} $olaps[0]{dir} $olaps[1]{dir} $flanks[0]{dir} $flanks[1]{dir} $e[0]{dir} $e[1]{dir}\n"; # if $e[0]{R} == 1211524;
 if ($flip) {
  @e = reverse @e;
  for (@e) {$$_{dir} *= -1; $$_{seq} = Revcomp($$_{seq})}
  @flanks = reverse @flanks;
  for (@flanks) {$$_{dir} *= -1}
 }
 for ($name, $flanks[0]{name}, $flanks[1]{name}) {s/\|//g}
 for (0,1) {my $d = $flanks[$_]{dir}; if ($d > 0) {$d = '>'} elsif ($d < 0) {$d = '<'}; $flanks[$_]{arrow} = $d}
 $name = $flanks[0]{name} . '|' . $flanks[1]{name} unless $name;
 my $flanksum = join('', $flanks[0]{name}, $flanks[0]{arrow}, $flanks[1]{arrow}, $flanks[1]{name});
 return (\@e, $flanksum, $name, $flip);
}

sub GeneEdge {
 my ($dna, $L, $R, $dir, $j, $ints) = @_;
 die "dna=$dna, L=$L, R=$R, dir=$dir, j=$j\n" unless $L;
 my ($highR, $Lig, $Rig, @overlaps, @out) = (-1, -1, -1); # Lig and Rig are gff lines for left and right nonoverlapping (in case intergene search), %overlaps are gff lines overlapping OL and OR
 for my $i (0 .. $#{$genes{$dna}}) {
  next if $genes{$dna}[$i]{questionable}; # Skip questionable (island-split fragment) tRNAs
  next if $genes{$dna}[$i]{annot} and $$ints{$genes{$dna}[$i]{annot}}; # Skip int within island
  if ($genes{$dna}[$i]{R} < $L and $genes{$dna}[$i]{R} > $highR) {($highR, $Lig) = ($genes{$dna}[$i]{R}, $i)} # Find Lig: gene with closest R end to left of OL and not overlapping OL
  push @overlaps, [GeneName($dna, $dir, $i)] unless $genes{$dna}[$i]{L} >= $R or $genes{$dna}[$i]{R} <= $L;
  if ($genes{$dna}[$i]{L} > $R) {$Rig = $i; last} # Find Rig: gene with closest L end to right of OR and not overlapping OR
 }
 if ($dir == 1) {@out = ($Lig)} else {@out = ($Rig)}
 #$dir *= -1 if $j;  # Keeps overlap gene order same for both ends, as in original tiger
 push @overlaps, [GeneName($dna, $dir, -1)] unless @overlaps;  # Empty
 my $pick = (sort {$$a[2] <=> $$b[2] || $$a[0] cmp $$b[0]} @overlaps)[0];
 return (GeneName($dna, $dir, $out[0]), @{$pick});  # flank, overlap
}

sub GeneName {
 my ($dna, $dir, $i) = @_;
 return (0, 0, $rankgenetype{empty}) if $i < 0;
 my ($name, $type) = ('');
 $dir *= -1 if $genes{$dna}[$i]{orient} eq '-';
 if ($genes{$dna}[$i]{info} =~ /annot=([^;]+)/) {$name = $1; $name =~ s/\.[0-9]+$//; $type = 'annot'}
 elsif ($genes{$dna}[$i]{info} =~ /product=tRNA.*-(.+)\(...\)/) {$name = $aalookup{$1}; $type = 'tRNA'}
 elsif ($genes{$dna}[$i]{info} =~ /product=tmRNA;/) {$name = 'Z'; $type = 'tRNA'}
 elsif ($genes{$dna}[$i]{info} =~ /gene=([^;]+)/) {$name = $1; $name =~ s/_\d+$//; $type = 'gene'}
 elsif ($genes{$dna}[$i]{info} =~ /pfam1=([^;]+)/) {$name = $1; $type = 'pfam'}
 elsif ($genes{$dna}[$i]{info} =~ /product=([^;]+)/) {$name = $1; $type = 'product'}
 elsif ($genes{$dna}[$i]{info} =~ /rpt_family=([^;]+)/) {$name = $1; $type = 'rpt'}
 else {$name = 'unknown'; $type = 'unknown'}
 if ($name =~ /^hypothetical protein$/) {$name = 'HYP'; $type = 'HYP'}
 $name =~ s/\s+/_/g;
 $name =~ s/[^a-zA-Z0-9_\-\.\(\)]//g;
 return ($name, $dir, $rankgenetype{$type})
}

sub Revcomp {my $ret = reverse $_[0]; $ret =~ tr/ACGTacgt/TGCAtgca/; return $ret}

sub MobGenes {
 my ($dna, $L, $R, $mobtype) = @_;
 my $genes = '';
 ($L, $R) = ($R, $L) if $L > $R;
 #print "mobzyme test for dna $dna interval $L-$R\n" if $verbose;
 for my $mobile (@{$mobs{$dna}}) {
  next if ($$mobile{L} > $R or $$mobile{R} < $L);
  $genes .= "$$mobile{id},";
  if ($mobtype eq "island") {
   next unless $$mobile{id} =~ /^Tnp_/;
   $$mobile{IS} = $$mobile{id}; #load $$mobile{IS} with mobile id only if there is a transposase hit in Q1 or Q2 -> runs TIGER.pl on these
  }  
 }
 #print "  Success\n" if $genes and $verbose;
 return $genes;
}

sub IntsWithin {
 my ($mobtype, $e) = @_;
 my %ints;
 for my $i (0,1) {
  my ($dna, $L, $R) = @{$$e[$i]}{qw/dna L R/};
  for my $int (@{$mobs{$dna}}) { # int must be over half within island
   next if $mobtype eq 'IS'     and $$int{id} !~ /^Tnp_/;
   next if $mobtype eq 'island' and $$int{id} !~ /Y-Int|S-Int|S-Core/;  # Changed in CROSS |Xer|Integron/;
   my ($within, $midpt) = (0, $$int{mid});
   $within = 1 if $midpt > $$e[$i]{L} && $midpt < $$e[$i]{R};
   %{$ints{$$int{id}}} = (dna => $dna, L => $$int{L}, R => $$int{R}) if $within;
  }
  last if $$e[$i]{compos} and $$e[$i]{compos} eq 'simple';
 }
 return \%ints;
}

sub Scores {
 my ($ints, $coord, $origin, $len, @coords) = @_;
 die unless $coord =~ s/([^\/]+)\/(\d+)-(\d+)//;
 @{$coords[0]}{qw/dna S E Ecorr dir/} = ($1, $2, $3, $3, 1);
 if ($2 == $3) {$coords[0]{dir} = -1 if $2>1} else {$coords[0]{dir} = ($3-$2)/abs($3-$2)}
 my $otherlen = abs($3-$2)+1;
 if ($coord =~ s/\+([^\/]+)\/(\d+)-(\d+)//) {
  @{$coords[1]}{qw/dna E S Ecorr dir/} = ($1, $2, $3, $2, 1);
  if ($2 == $3) {$coords[0]{dir} = -1 if $2>1} else {$coords[0]{dir} = ($3-$2)/abs($3-$2)}
  $coords[1]{Ecorr} -= $coords[1]{dir} * $otherlen;
  $coords[0]{Ecorr} += $coords[0]{dir} * abs($3-$2)+1;
 }
 for my $i (0, 1) {next unless $coords[$i]; @{$coords[$i]}{qw/L R/} = sort {$a <=> $b} @{$coords[$i]}{qw/S E/}}
 my $delta_int = DeltaInt($ints, \@coords, $origin, $len);
 my ($hypoth, $forn, $hskp) = Foreign($ints, \@coords, $origin);
 my ($delta_GC, $dinuc) = Bias(\@coords, $origin);
 my $overall = exp(-35.2 + 6.4*log10($len) + 1.2*log10($delta_int) + 0.006*$forn + -0.39*$hskp + -5.3*$hypoth + 7.0*abs($delta_GC) + -25.7*$dinuc);
 return ($delta_int, $forn, $hskp, $hypoth, $delta_GC, $dinuc, $overall);
}

sub log10 {return log(shift)/log(10);}

sub DeltaInt { # Shortest distance between any internal integrase gene end and an island end
 my ($ints, $c, $ret) = @_;
 for my $i (0,1) {
  next unless $$c[$i];
  my ($dna, $L, $R) = ($$c[$i]{dna}, sort {$a <=> $b} @{$$c[$i]}{qw/S Ecorr/});
  for (keys %{$ints}) {
   next unless $dna eq $$ints{$_}{dna};
   my ($intL, $intR) = ($$ints{$_}{L}, $$ints{$_}{R});
   $ret = $intL-$L if $intL-$L < $ret;
   $ret = $R-$intR if $R-$intR < $ret;
  }
 }
 $ret = 1 if $ret < 1;
 return $ret
}

sub Foreign {
 my ($ints, $c) = @_;
 my ($cds, $hypoth, $pfam, $forn, $hskp) = (0,0,0,0,0);
 for my $i (0,1) {
  next unless $$c[$i];
  my ($dna, $L, $R) = ($$c[$i]{dna}, sort {$a <=> $b} @{$$c[$i]}{qw/S E/});
  for my $prot (@{$genes{$dna}}) {
   next unless $$prot{line} =~ /\tCDS\t/;
   unless ($$ints{$$prot{id}}) {
    next unless $$prot{L} >= $L and $$prot{R} <= $R; last if $$prot{L} > $R;
   }
   $cds ++;
   if ($$prot{pfam} eq 'hypothetical') {$hypoth ++}
   else {
    $pfam ++;
    $forn += $fornEnrich{$$prot{pfam}} if $fornEnrich{$$prot{pfam}};
    $hskp += $hskpEnrich{$$prot{pfam}} if $hskpEnrich{$$prot{pfam}};
   }
  }
 }
 $hypoth /= $cds if $cds;   $hypoth -= $stats{all}{hypoth}; # Difference between hypoth Densities of island and its replicon
 $forn   /= $pfam if $pfam; $forn   -= $stats{all}{forn};
 $hskp   /= $pfam if $pfam; $hskp   -= $stats{all}{hskp};
 return $hypoth, $forn, $hskp;
}

sub Bias {
 my ($c) = @_;
 print STDERR "perl $dir/collectSeq.pl -i $prefix.fa -e $$c[0]{dna} -L $$c[0]{L} -R $$c[0]{R} >  test.fa";
 system  "perl $dir/collectSeq.pl -i $prefix.fa -e $$c[0]{dna} -L $$c[0]{L} -R $$c[0]{R} >  test.fa";
 if ($$c[1]) {
  print STDERR "perl $dir/collectSeq.pl -i $prefix.fa -e $$c[1]{dna} -L $$c[1]{L} -R $$c[1]{R} >> test.fa";
  system "perl $dir/collectSeq.pl -i $prefix.fa -e $$c[1]{dna} -L $$c[1]{L} -R $$c[1]{R} >> test.fa";
 }
 print STDERR "perl $dir/relAbun.pl test.fa";
 my $out = `perl $dir/relAbun.pl test.fa`; chomp $out; $out =~ s/\n.*//s; my @testRA = split "\t", $out;
 my ($delta_GC, $di) = (($testRA[1]-$dnaRA{all}[1])/2, 0);
 for (2,3,4,6,7,9) {$di += abs($testRA[$_]-$dnaRA{all}[$_])}
 $di *= 2; # Double the 6 asymmetrical dinucs above to account for their complements, but don't double the 4 symmetrical dinucs below
 for (5,8,10,11) {$di += abs($testRA[$_]-$dnaRA{all}[$_])}
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
my $version = '2.0 (May 2022)';
#   '    |    '    |    '    |    '    |    '    |    '    |    '    |    '    |
my $help = <<END;
$scriptname version $version
Usage: perl $scriptname [options] -db <RefDatabase> -fasta <GenomicDNA>
  -fasta:    Genomic fasta DNA sequence file.
  -db:       Blast database of reference genomes, absolute path.
  -search:   Search type. Specify island or IS. Default: $search.
  -tax:      Taxonomic info for query genome. Enter name of a file containing 
              NCBI taxonomy string, or use B for Bacteria, A for Archaea, M for 
              Mycoplasmatales/Entomoplasmatales, G for Gracilibacteria/candidate
              division SR1. Automatically sets -gencode. Default: B.
  -gencode:  Genetic code table to use (see NCBI). Default: 11.
  -nickname: Brief name for genome (as might be used to start a locus_tag).
  -circle:   Specify C if all genomic DNA sequences are circular, L if all DNAs 
              are linear, or a filename for a tab-delimited file of query and 
              circularity (eg. acc.vers[tab]circular/linear). Default: $circularity.
  -cross:    Three options: intact, cross, or circleOrigin. Default: $cross. 
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
 'cross=s' => \$cross,
 'tax=s' => \$tax,
 'nickname=s' => \$nickname,
 'cpu=i' => \$cpu,
 'circle=s' => \$circularity,
 'search=s' => \$search,
);
die $help if !$options_okay;
die "-db and -fasta are required\n\n$help" unless $inDna and $db;
die "ERROR: illegal cross option\n" . $help unless $cross =~ /^intact|cross|circleOrigin$/;
}
