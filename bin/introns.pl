#! /usr/bin/perl
use strict; use warnings;
# autocorrect aragorn miscall of anticodon loop with lower-casing of spliced.struct

die "Usage: $0 infasta outdir [AMG taxa]\n" if @ARGV < 2;
use Cwd 'abs_path';
my ($infile, $outdir, $binpath) = ($ARGV[0], $ARGV[1], $0); 
for ($infile, $binpath) {$_ = abs_path($_); die "Can't find infasta\n" unless $_}
$binpath =~ s/\/([^\/]+)$//; my $scriptname = $1;
my $lib = $binpath; $lib =~ s/[^\/]*$/db/;
my ($tax, $gencode, %anticodon2aa, %one_letter, %final, %n, $id, %serials, $dna, $collect, @gpi, $tbed, %ct, %xcptns) = ('B', '');

LoadCode();
LoadExceptions();
$tax = $ARGV[2] if $ARGV[2];
if ($tax eq 'M') {($tax, $gencode, $anticodon2aa{CTA}) = ('B', "-g $lib/gencode/gcode.othmito", 'Trp')}  # gencode=4
if ($tax eq 'G') {($tax, $gencode, $anticodon2aa{CTA}) = ('B', "-g $lib/gencode/gcode.gracil" , 'Gly')}  # gencode=25
mkdir $outdir;   chdir $outdir;
#`rm -r introns`; 
#`rm introns.gff`; # DELETE
open FINAL, ">introns.gff" or die "Can't open outfile $outdir/introns.gff: $!\n";
mkdir 'introns'; chdir 'introns';
#`rm spliced*`;  # ADD
#for (`ls`) {print $_}
#for (`pwd`) {print $_}
#warn "No gpi.tbl\n" unless -f 'gpi.tbl'; exit;

Introns();
warn "Examining gpI Introns in $infile\n";
unless (-s 'gpi.bed') {warn "No gpI Introns detected\n"; exit}
my $flag = Trna();
if ($flag) {warn "$flag\n"; exit}
for my $id (sort keys %final) {
 $serials{final} ++;
 print FINAL join("\t", $n{$id}{dna}, qw/tfind tRNA_gpIintron/, @{$n{$id}}{qw/L R score ori/}, '.', "ID=intron.$serials{final};");
 for (qw/aa ac a_site j_site disc_fill disc_end iso iso_score introns note gpi_range gpi_L gpi_R gpi_score splice_site splice_score ivs_L ivs_R intron_length seq cca struct/)
 {if (defined $n{$id}{$_}) {print FINAL "$_=$n{$id}{$_};"} else {print FINAL "$_=;"}}
 print FINAL "splices=", join(',', @{$n{$id}{splices}}), ";\n";
}
close FINAL;

# SUBROUTINES
sub Introns {
 my $cmd = "cmsearch --tblout gpi.tbl $lib/cm/Intron_gpI.cm $infile &> /dev/null";
 unless (-s 'gpi.tbl') {warn "$cmd\n"; system $cmd}
 warn "Can't find gpi.tbl\n" unless -f 'gpi.tbl';
 my @gpi;
 for (`cat gpi.tbl`) {
  next if /^#/;  my @f = split /\s+/;
  ##target_name         accession query_name           accession mdl mdl_from   mdl_to seq_from   seq_to strand trunc pass   gc  bias  score   E-value inc description_of_target
  #Intron_gpI           RF00028   AM747720.1           -          cm       49      160  1222644  1222744      +    no    1 0.62   2.6   28.9   4.5e-06 !   -
  next if $f[14] < 20;  # CM score cutoff for Intron_gpI; cm "trusted cutoff" is 40.01
  $f[10] =~ s/\'/prime/g; $f[10] =~ s/5prime\&3prime/both/; my $range = "$f[10]$f[5]..$f[6]";
  ($f[7], $f[8]) = ($f[8], $f[7]) if $f[7] > $f[8];
  $serials{gpi} ++;
  push @gpi, [$f[0], $f[7]-1, $f[8], "gpi$serials{gpi}", $f[14], $f[9], $range];
 }
 return unless @gpi;
 open OUT, ">gpi.bed"; for (sort {$$a[0] cmp $$b[0] || $$a[1] <=> $$b[1]} @gpi) {print OUT join("\t", @{$_}), "\n"} close OUT;
}

sub Trna {
 warn "Finding split tRNA gene calls\n";
 system "$binpath/aragorn1.1 -w -t -l -i -seq -o introns.aragorn $infile" unless -f "introns.aragorn";  # aragorn1.1
 for (`cat introns.aragorn`) {
  if (/^>(\S+)/) {$dna = $1; $id = ''}
  elsif (/^T /) {$id = ''}
  elsif (/^([a-z\.]+)\s*$/) {next unless $id; $n{$id}{seq} .= $1; $n{$id}{aragorn} .= $_}
  elsif (/^TI/) {
   die "ERROR: Cannot parse aragorn1.1 $infile $_" unless /^TI\s+.*(.)\[([\-\d]+),(\d+)\]\s+(\d+)\s+\S+\s+(\d+)\s+(\d+)/;
   $serials{trna} ++; $id = "trna$serials{trna}";
   if ($6 < 50) {warn "$id: Intron too short\n"; $id = ''; next}
   %{$n{$id}} = (dna => $dna, dir => 1, ori => '+', L => $2, R => $3, a_site => $4, intron_coord => $5, intron_len => $6, aragorn => $_);
   if ($1 eq 'c') {@{$n{$id}}{qw/dir ori/} = (-1, '-'); $n{$id}{L} += 3} else {$n{$id}{R} -= 3}
   $tbed .= join("\t", $n{$id}{dna}, $n{$id}{L}-1, $n{$id}{R}, $id, '', $n{$id}{ori}) . "\n";
  } else {die "ERROR: Cannot parse aragorn1.1 $infile $_"}
 }
 return "No split trna calls" unless scalar(keys %n);
 warn scalar(keys %n), " original split trna calls\n";
 for my $id (keys %n) {  # Truncations
  $n{$id}{seq} = uc $n{$id}{seq};
  $n{$id}{seq} =~ s/(...)$//; $n{$id}{cca} = uc $1; $n{$id}{cca} =~ s/\./N/g; $n{$id}{cca} =~ s/N*$//;
  if ($n{$id}{seq} =~ s/^(\.+)//) {my $l = length $1; if ($n{$id}{dir} > 0) {$n{$id}{L} += $l; $n{$id}{trunc_start} = $l} else {$n{$id}{R} -= $l; $n{$id}{trunc_end} = $l}}
  if ($n{$id}{seq} =~ s/(\.+)$//) {my $l = length $1; if ($n{$id}{dir} > 0) {$n{$id}{R} -= $l; $n{$id}{trunc_end} = $l} else {$n{$id}{L} += $l; $n{$id}{trunc_start} = $l}}
  $n{$id}{seq} =~ s/\./N/g;
  $n{$id}{seq} .= 'N' x $n{$id}{trunc_end} if $n{$id}{trunc_end};  # Ns improve tSE calls on truncations
  warn "ERROR: $id: internal non-acgt\n" if $n{$id}{seq} =~ /\./;
 }
 warn "Finding group I intron-bearing split tRNA gene calls\n";
 for (`echo '$tbed' | bedtools intersect -s -wo -f 1 -a gpi.bed -b stdin`) {  # -f insists that intron be fully contained within split trna
  # AM747720.1	1222643	1222744	gpi1	28.9	+	no49..160	AM747720.1	1222570	1222903	trna1		+	101
  my @f = split "\t";
  push @{$n{$f[10]}{gpis}}, [$f[1]+1, @f[2,4,6]];
  $ct{intersect} ++;
 }
 return "No gpI Introns overlapped split trna calls" unless $ct{intersect};
 open OUT, ">spliced.fa";  # Prepare 8 anticodon loop splicings
 my (%seqs, %lengths);
 for my $id (sort keys %n) {  # Make all 8 splicings of given intron length from anticodon loop
  unless ($n{$id}{gpis}) {warn "$id: No Intron_gpI hit\n"; delete $n{$id}; next}
  my ($seq, $ac, $in, $len) = @{$n{$id}}{qw/seq a_site intron_coord intron_len/};
  if ($xcptns{$seq} and $xcptns{$seq}{intron_len}) {$n{$id}{intron_len} = $len = $xcptns{$seq}{intron_len}; warn "Exception detected\n"}
  #warn "$seq, $ac, $in, $len\n";
  my $pos = $ac -3;
  $pos -= $len if $ac > $in;
  for (0..7) {
   my ($exon1, $intron, $exon2) = (substr($seq, 0, $pos), lc (substr($seq, $pos, $len)), substr($seq, $pos+$len));
   my $sum = "$_" . substr($exon1, -1) . substr($intron, -1); $sum =~ tr/tT/uU/;
   print OUT ">$id.$sum\n$exon1${exon2}CCA\n";
   $lengths{$id} = length($exon1 . $exon2);
   $seqs{$id}{$sum} = "$exon1$intron$exon2";
   $pos ++;
  }
 }
 close OUT;

 warn "Choosing splice sites for the ", scalar(keys %n), " split trna calls with a gpI intron hit\n";
 my $spliced = Scan();
 return "No spliced trnas called for these gpI introns" unless keys %{$spliced};
 warn scalar(keys %{$spliced}), " trna calls for the 8 splice sites tested for each previous\n";
 my (%scores, %tops, %top_splice, %splice_scores, @all);
 for my $id (sort keys %{$spliced}) {  # Find top scoring and top splice site scores
  for my $sum (sort keys %{$$spliced{$id}}) {
   warn "ERROR: Cannot parse splice $sum $id\n" unless $sum =~ /^(\d)([A-Za-z])([A-Za-z])/;
   my ($pos, $pre, $omega) = ($1, $2, $3);
   my $splice = "$one_letter{$$spliced{$id}{$sum}{aa}}$sum" . 10*$$spliced{$id}{$sum}{iso_score} . "$one_letter{$$spliced{$id}{$sum}{iso}}";
   push @{$n{$id}{splices}}, $splice;
   my $score = 0; $score += 3 if $omega eq 'g'; $score += 2 if $pre eq 'U'; $score += 1 if $pre eq 'C';  # Max = 5 for Ug
   $$spliced{$id}{$sum}{splice_score} = $score;
  }
 }
 my $ct = 0;
 for my $id (sort keys %n) {  # Choose top splice site, then best trna score
  unless (keys %{$$spliced{$id}}) {warn "$id: No tSE hits\n"; next}
  my $pick = (sort {$$spliced{$id}{$b}{splice_score} <=> $$spliced{$id}{$a}{splice_score} ||
                    $$spliced{$id}{$b}{iso_score}    <=> $$spliced{$id}{$a}{iso_score}    ||  
                    $$spliced{$id}{$b}{IPD}          <=> $$spliced{$id}{$a}{IPD}            } keys %{$$spliced{$id}})[0];
  $n{$id}{splice_site} = $pick; #$n{$id}{splice_site} =~ s/.*\.(\d[A-Z][a-z]).*/$1/;
  $n{$id}{seq} = $seqs{$id}{$pick};
  #print "$id $pick $n{$id}{seq}\n";
  if ($$spliced{$id}{$pick}{L} != 1) {  # Aragorn may overextend L or R end; correct
   my $diff = $$spliced{$id}{$pick}{L} - 1;
   if ($n{$id}{dir} == 1) {$n{$id}{L} += $diff} else {$n{$id}{R} -= $diff}
   $n{$id}{seq} =~ s/^.{$diff}// if $diff;
   for ($n{$id}{a_site}, $n{$id}{gpI_L}, $n{$id}{gpI_R}) {$_ -= $diff}
  }
  if ($$spliced{$id}{$pick}{R} < $lengths{$id}) {  # Aragorn may overextend L or R end; correct, but may have added Ns
   my $diff = $lengths{$id} - $$spliced{$id}{$pick}{R};
   #if ($$spliced{$id}{$pick}{cca}) {
   for (1..$diff) {
    $n{$id}{seq} =~ s/(.)$//; next if $1 eq 'N'; $n{$id}{cca} = $1 . $n{$id}{cca} ; 
    if ($n{$id}{dir} == 1) {$n{$id}{R} --} else {$n{$id}{L} ++}
   }
   $n{$id}{cca} =~ s/^(...).*/$1/;
  }
  $n{$id}{seq} =~ s/N*$//;
  my $seq = $n{$id}{seq};
  $seq =~ s/[A-Z]+$//; $n{$id}{ivs_R} = length $seq;
  $seq =~ s/[a-z]+$//; $n{$id}{ivs_L} = length $seq;
  $n{$id}{intron_length} = $n{$id}{ivs_R} - $n{$id}{ivs_L};
  warn "$id: intron too short $n{$id}{intron_length}\n" if $n{$id}{intron_length} < 150;
  if ($$spliced{$id}{$pick}{introns}) {warn "ERROR: $id: tRNAscan-SE called post-splicing intron\n"; next}
  $ct ++;
  @{$n{$id}}{                   qw/struct ac aa a_site score iso iso_score disc_fill disc_end introns note splice_score trunc_start trunc_end/}
    = (@{$$spliced{$id}{$pick}}{qw/str    ac aa a_site score iso iso_score disc_fill disc_end introns note splice_score trunc_start trunc_end/});
  $n{$id}{a_site} += $n{$id}{ivs_R} - $n{$id}{ivs_L} if $n{$id}{a_site} > $n{$id}{ivs_L};
  $n{$id}{aa} = $anticodon2aa{$n{$id}{ac}} if $n{$id}{ac} eq 'CTA' and $anticodon2aa{CTA};  # Two alternate genetic codes
  Fix($id, \%{$$spliced{$id}}) if $n{$id}{aa} eq 'Undet';
  #$n{$id}{j_site_post} = $lengths{$id} - 7; $n{$id}{j_site_post} -- if $n{$id}{aa} eq 'SeC';
  $n{$id}{j_site} = $n{$id}{R} - $n{$id}{L} - 6; $n{$id}{j_site} -- if $n{$id}{aa} eq 'SeC';
  my ($L, $R, $gpi_score, $gpi_range) = @{(sort {$$b[2] <=> $$a[2]} @{$n{$id}{gpis}})[0]};  # Choose best if multiple gpI calls
  if ($n{$id}{ori} eq '+') {for ($L, $R) {$_ = $_-$n{$id}{L}+1}} else {for ($L, $R) {$_ = $n{$id}{R}-$_+1} ($L, $R) = ($R, $L)}
  @{$n{$id}}{qw/gpi_L gpi_R gpi_score gpi_range/} = ($L-1, $R, $gpi_score, $gpi_range);
  $final{$id} ++;
 }
 warn "$ct gpI intron-containing tRNAs survive\n";
 return '';
}

sub Fix {  # Aragorn undercalled intron length in 26 uniq tRNAs leading to lc segment in ac loop in spliced.struc
 #warn "Possible miscalled intron, try adjusting based on tRNAscan-SE ac loop excess\n";
 #my ($id, $spliced) = @_;
 #for my $sum (sort keys %{$spliced}) { }
 warn "Possible error in intron length, aa Undet\n";
}

sub Scan {
 my (%t, $splice);
 my $call = "tRNAscan-SE -Q -q -f spliced.struct -o spliced.out $gencode -$tax --detail --brief spliced.fa &> /dev/null";  # Also like -s $base.iso, but not strictly necessary
 unless (-f "spliced.out") {warn "$call\n"; system $call}
 for (`cat spliced.struct`) {
  next unless /^\S/;
  chomp;
  if (/^(\S+).trna\d* \(\d*-\d*/) {$id = $1; $id =~ s/\.([^\.]+)$//; $splice = $1}
  elsif (/^Type: \S+\s+Anticodon: \S+ at (\d+)\S+ \S+\s+Score: \S+/) {$t{$id}{$splice}{a_site} = $1 + 1}
  elsif (/^Possible intron: (\d+)-(\d+)/) {push @{$t{$id}{$splice}{intLs}}, $1; push @{$t{$id}{$splice}{intRs}}, $2}
  elsif (/^Possible/) {}
  elsif (s/^Pre: (\S+).*/$1/) {while (s/\[([^\]]*)\]/lc($1)/e) {} $t{$id}{$splice}{seq} = $_}  # Lower-case (multiple) introns
  elsif (/^Seq: (\S+)/) {$t{$id}{$splice}{mat} = $1}
  elsif (/^Str: (\S+)/) {$t{$id}{$splice}{str} = $1}
  else {warn "ERROR: Cannot parse tRNAscan-SE line $_\n"}
 }
 for (`cat spliced.out`) {
  #fake            13      2000494 2000565 Undet   NNN     0       0       20.3    Arg     20.2    pseudo,trunc_start:32
  chomp; s/\s+\d+//; my @f = split "\t"; $id = $f[0]; $id =~ s/\.([^\.]+)$//; $splice = $1; $f[10] = '' unless $f[10];
  die "Cannot parse dna $_\n" unless /^(\S+)\t/;
  @{$t{$id}{$splice}}{qw/dna dir L R aa ac score iso iso_score note IPD/} = ($1, 1, @f[1..4], @f[7..10], 0);
  if ($t{$id}{$splice}{L} == 0) {$t{$id}{$splice}{L} = 1; warn "ERROR: Hard-trunc unresolved $_\n" unless $t{$id}{$splice}{str} =~ s/^\.>/>/}  # Correctable error noted on truncated Ype23.trna11 et al
  if ($f[1] > $f[2]) {warn "ERROR: tRNA call in wrong orientation\n"; next}
  $t{$id}{$splice}{IPD} = $1 if $f[10] =~ /IPD:([0-9\.\-]+)/;
  $t{$id}{$splice}{aa} = $anticodon2aa{$t{$id}{$splice}{ac}} if $t{$id}{$splice}{ac} eq 'CTA' and $anticodon2aa{CTA};  # Two alternate genetic codes
  if ($t{$id}{$splice}{seq}) {
   my $test = $t{$id}{$splice}{seq}; $test =~ s/[a-z]+//g; warn "ERROR: spliced pre ne mature $id\n" if $test ne $t{$id}{$splice}{mat};
  } else {$t{$id}{$splice}{seq} = $t{$id}{$splice}{mat}}  # Seq is pre-sequence with introns, mat is mature spliced
  $t{$id}{$splice}{str} =~ /.*</;  # Find last position of acceptor stem
  $t{$id}{$splice}{disc_end} = $+[0] +1;
  my ($expAcc, $len, $toy, @loops) = (7, length($t{$id}{$splice}{mat}), $t{$id}{$splice}{str});
  $expAcc ++ if $t{$id}{$splice}{aa} eq 'SeC';
  while ($toy =~ s/>\.+<([\.<X]+)$/X$1/) {  # Replace loops, marking positions of last stem ntds
   unshift @loops, $-[1];
   while ($toy =~ s/>\.*X(\.*<)/X/) {$loops[-1] += length $1}
  }
  my $accLen = $toy =~ s/>/>/g; $accLen = 0 unless $accLen;  #  Count acceptor stem bp after removing loops
  $t{$id}{$splice}{j_site} = $loops[-1];
  $t{$id}{$splice}{disc_fill} = $t{$id}{$splice}{disc_end} + $expAcc - $accLen;  # Discriminator based on filling from end of (short?) acc
  $t{$id}{$splice}{disc_end} = $t{$id}{$splice}{j_site} + $expAcc + 1 unless $accLen;  # Discriminator based on last acc bp (otherwise filling from end of tstem)
  if ($len > $t{$id}{$splice}{disc_fill} and $t{$id}{$splice}{seq} =~ s/CCA$// and $t{$id}{$splice}{str} =~ s/\.{3}$//) {$t{$id}{$splice}{cca} = 'disc_fill'}
  elsif ($len == $t{$id}{$splice}{disc_end} +3 and $t{$id}{$splice}{seq} =~ s/CCA$// and $t{$id}{$splice}{str} =~ s/\.{3}$//) {$t{$id}{$splice}{cca} = 'disc_end'}
  #warn "disc_fill=$t{$id}{$splice}{disc_fill}, disc_end=$t{$id}{$splice}{disc_end}, len=$len, cca=$t{$id}{$splice}{cca}\n";
  if ($t{$id}{$splice}{cca}) {if ($t{$id}{$splice}{dir} > 0) {$t{$id}{$splice}{R} -= 3} else {$t{$id}{$splice}{L} += 3}}
  my ($inttot, @introns) = (0);
  for my $i (0 .. $#{$t{$id}{$splice}{intLs}}) {
   push @introns, "$t{$id}{$splice}{intLs}[$i]..$t{$id}{$splice}{intRs}[$i]";
   my $intlen = $t{$id}{$splice}{intRs}[$i] - $t{$id}{$splice}{intLs}[$i] +1;
   die "$id $splice a_site\n" unless $t{$id}{$splice}{a_site};
   $t{$id}{$splice}{a_site} += $intlen if $t{$id}{$splice}{a_site} >= $t{$id}{$splice}{intLs}[$i];
   $inttot += $intlen;
  }
  $t{$id}{$splice}{introns} = join(',', @introns);
  for (qw/j_site disc_fill disc_end/) {$t{$id}{$splice}{$_} += $inttot}  # These had been based on mature seq
 }
 return \%t;
}

sub LoadCode {
 %anticodon2aa =      qw/TGA Ser GGA Ser CGA Ser AGA Ser GAA Phe AAA Phe TAA Leu CAA Leu GTA Tyr ATA Tyr GCA Cys ACA Cys CCA Trp 
 TAG Leu GAG Leu CAG Leu AAG Leu TGG Pro GGG Pro CGG Pro AGG Pro GTG His ATG His TTG Gln CTG Gln TCG Arg GCG Arg CCG Arg ACG Arg 
 TAT Ile GAT Ile AAT Ile CAT Met TGT Thr GGT Thr CGT Thr AGT Thr GTT Asn ATT Asn TTT Lys CTT Lys GCT Ser ACT Ser TCT Arg CCT Arg 
 TAC Val GAC Val CAC Val AAC Val TGC Ala GGC Ala CGC Ala AGC Ala GTC Asp ATC Asp TTC Glu CTC Glu TCC Gly GCC Gly CCC Gly ACC Gly/;
 %one_letter = (qw/Ala A Arg R Asn N Asp D Cys C Gln Q Glu E Gly G His H Ile I Ile2 I Leu L Lys K Met M Phe F Pro P SeC U Ser S Thr T Trp W Tyr Y Val V fMet B iMet B Sup * Undet X None X/);
}

sub LoadExceptions {
 for (`cat $lib/introns.exceptions`) {next unless /(\S+)\s(\S+)/; my ($seq, $info) = ($1, $2); for (split ';', $info) {/([^=]+)=(.*)/; $xcptns{uc($seq)}{$1} = $2}}
}
