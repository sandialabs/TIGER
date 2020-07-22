#! /usr/bin/perl
use strict; use warnings;
# autocorrect aragorn miscall of anticodon loop with lower-casing of spliced.struct

die "Usage: $0 fasta [AMG taxa]\n" if @ARGV < 1;
use File::Spec;
my ($file, $binpath) = ($ARGV[0], $0); 
for ($file, $binpath) {$_ = File::Spec->rel2abs($_); die "Can't find infasta\n" unless $_}
my $dir = '.'; $dir = $1 if $file =~ /(.*)\/[^\/]$/;
$binpath =~ s/\/([^\/]+)$//; my $scriptname = $1;
my $lib = $binpath; $lib =~ s/[^\/]*$/db/;
my ($tax, $gencode, %anticodon2aa, %one_letter, %final, %n, $id, %serials, $dna, $collect, @gpi, $tbed, $tbedmid, %ct, %xcptns, %spliced) = ('B', '');
my @isos = qw/Ala Arg Asn Asp Cys Gln Glu Gly His Ile Ile2 Leu Lys Met Phe Pro SeC Ser Thr Trp Tyr Val fMet/;

LoadCode();
$tax = $ARGV[1] if $ARGV[1];
if ($tax eq 'M') {($tax, $gencode, $anticodon2aa{CTA}) = ('B', "-g $lib/gencode/gcode.othmito", 'Trp')}  # gencode=4
if ($tax eq 'G') {($tax, $gencode, $anticodon2aa{CTA}) = ('B', "-g $lib/gencode/gcode.gracil" , 'Gly')}  # gencode=25
$tax = 'B' if $tax eq 'A';  # Archaeal mode ruins bacterial tRNAs that have been miscalled as archaea like Par15
chdir $dir; mkdir 'trna'; chdir 'trna';
open FINAL, ">introns.gff" or die "Can't open outfile $dir/introns.gff: $!\n";
mkdir 'introns'; chdir 'introns';
system "rm -r spliced* pre* ghost* &> /dev/null" if -f "spliced.fa"; # XXXXXXXXXX
warn "Examining gpI Introns in $file\n";
Introns();
unless (-s 'gpi.bed') {warn "No gpI Introns detected\n"; exit}
my $flag = Trna();
if ($flag) {warn "$flag\n"; exit}
for my $id (sort keys %final) {
 #my $name = "$n{$id}{aa}.$n{$id}{splice_site}"; $serials{$name} ++; if $serials{$name}
 $serials{final} ++;
 $n{$id}{note} = join(',', @{$n{$id}{notes}}) if $n{$id}{notes};
 $n{$id}{seq} =~ s/$n{$id}{ttc}/$n{$id}{torig}/;
 $n{$id}{start} = 0 unless $n{$id}{start};
 #warn "$n{$id}{ttc}/$n{$id}{torig} $n{$id}{ttcold} $n{$id}{start} $id\n";
 if ($n{$id}{dir}  == 1) {$n{$id}{L} += $n{$id}{start}; $n{$id}{R} = $n{$id}{L} + length($n{$id}{seq}) -1}
 else           {$n{$id}{R} -= $n{$id}{start}; $n{$id}{L} = $n{$id}{R} - length($n{$id}{seq}) +1}
 $n{$id}{L} = 1 if $n{$id}{L} < 1;
 print FINAL join("\t", $n{$id}{dna}, qw/tfind tRNA_gpIintron/, @{$n{$id}}{qw/L R score ori/}, '.', "ID=$n{$id}{aa}.$n{$id}{splice_site};");
 for (qw/aa ac a_site j_site disc_fill disc_end iso iso_score introns note gpi_range gpi_L gpi_R gpi_score splice_site splice_score ivs_L ivs_R intron_length seq cca struct/)
 {if (defined $n{$id}{$_}) {print FINAL "$_=$n{$id}{$_};"} else {print FINAL "$_=;"}}
 print FINAL "ttc=$n{$id}{ttcold};splices=", join(',', @{$n{$id}{splices}}), ";\n";
}
close FINAL;

# SUBROUTINES
sub Trna {
 Aragorn1();
 return "No split trna calls" unless scalar(keys %n);
 warn "Finding group I intron-bearing split tRNA gene calls\n";
 for (`echo '$tbed' | bedtools intersect -s -wo -f 1 -a gpi.bed -b stdin`) {  # -f insists that intron be fully contained within split trna
  # AM747720.1	1222643	1222744	gpi1	28.9	+	no49..160	AM747720.1	1222570	1222903	trna1		+	101
  my @f = split "\t";
  push @{$n{$f[10]}{gpis}}, [$f[1]+1, @f[2,4,6]];
  $ct{intersect} ++;
 }
 for (`echo '$tbedmid' | bedtools intersect -s -wo -a gpi.bed -b stdin`) {  # allow if intron covers midpoint of split trna, even if intron call not fully within split trna
  my @f = split "\t";
  next if $n{$f[10]}{gpis};
  push @{$n{$f[10]}{gpis}}, [$f[1]+1, @f[2,4,6]];
  $ct{intersect} ++;
 }
 return "No gpI Introns within a split trna call" unless $ct{intersect};
 Ghost();
 Tloop();
 Presplice();
 open OUT, ">spliced.fa";  # Prepare 8 anticodon loop splicings
 my (%seqs, %lengths);
 for my $id (sort keys %n) {  # Make all 8 splicings of given intron length from anticodon loop
  #my ($seq, $ac, $in, $len, $ttc, $torig) = @{$n{$id}}{qw/seq a_site intron_coord intron_len ttc torig/};
  my ($seq, $ac, $len, $ttc, $torig) = @{$n{$id}}{qw/seq a_site intron_len ttc torig/};
  $seq =~ s/$torig/$ttc/;
  my $pos = $ac -4;
  #$pos -= $len if $ac > $in;
  #print "$seq, $ac, $in, $len, $ttc, $torig $pos\n";
  for (0..7) {
   my ($exon1, $intron, $exon2) = (substr($seq, 0, $pos), lc (substr($seq, $pos, $len)), substr($seq, $pos+$len));
   my $sum = "$_" . substr($exon1, -1) . substr($intron, -1); $sum =~ tr/tT/uU/;
   #print OUT ">$id.$sum\n$exon1${exon2}CCA\n";
   print OUT ">$id.$sum\n$exon1${exon2}\n";
   $lengths{$id} = length($exon1 . $exon2);
   $seqs{$id}{$sum} = "$exon1$intron$exon2";
   $pos ++;
  }
 }
 close OUT;

 warn "Choosing splice sites for the ", scalar(keys %n), " split trna calls with a gpI intron hit\n";
 my $scanout = Scan(); %spliced = %{$scanout};
 Aragorn2() unless scalar(keys %spliced) == scalar(keys %n);  # Backup if tRNAscan-SE misses tRNA
 for my $id (sort keys %spliced) {  # Find top scoring and top splice site scores
  for my $sum (sort keys %{$spliced{$id}}) {
   warn "ERROR: Cannot parse splice $sum $id\n" unless $sum =~ /^(\d)([A-Za-z])([A-Za-z])/;
   my ($pos, $pre, $omega) = ($1, $2, $3);
   my $splice = "$one_letter{$spliced{$id}{$sum}{aa}}$sum" . sprintf("%.0f", 10*$spliced{$id}{$sum}{iso_score}) . "$one_letter{$spliced{$id}{$sum}{iso}}";
   push @{$n{$id}{splices}}, $splice;
   my $score = 0; $score += 3 if $omega eq 'g'; $score += 2 if $pre eq 'U'; $score += 1 if $pre eq 'C';  # Max = 5 for Ug
   $spliced{$id}{$sum}{splice_score} = $score;
   #warn "$id: $sum sco=$score nom=$spliced{$id}{$sum}{aa} isco=$spliced{$id}{$sum}{iso_score} iso=$spliced{$id}{$sum}{iso} ipd=$spliced{$id}{$sum}{IPD}\n";
  }
 }
 my $ct = 0;
 for my $id (sort keys %n) {  # Choose top splice site, then best trna score
  unless (keys %{$spliced{$id}}) {warn "$id: No tSE hits\n"; next}
  my $pick = (sort {$spliced{$id}{$b}{splice_score} <=> $spliced{$id}{$a}{splice_score} ||
                    $spliced{$id}{$b}{iso_score}    <=> $spliced{$id}{$a}{iso_score}    # Best iso corrected by IPD = nominal iso
  } keys %{$spliced{$id}})[0];
  $n{$id}{splice_site} = $pick; #$n{$id}{splice_site} =~ s/.*\.(\d[A-Z][a-z]).*/$1/;
  $n{$id}{seq} = $seqs{$id}{$pick};
  #print "$id $pick $n{$id}{seq}\n";
  if ($spliced{$id}{$pick}{L} != 1) {  # Aragorn may overextend L or R end; correct
   my $diff = $spliced{$id}{$pick}{L} - 1;
   if ($n{$id}{dir} == 1) {$n{$id}{L} += $diff} else {$n{$id}{R} -= $diff}
   $n{$id}{seq} =~ s/^.{$diff}// if $diff;
   for ($n{$id}{a_site}, $n{$id}{gpI_L}, $n{$id}{gpI_R}) {$_ -= $diff}
  }
  if ($spliced{$id}{$pick}{R} < $lengths{$id}) {  # Aragorn may overextend L or R end; correct, but may have added Ns
   my $diff = $lengths{$id} - $spliced{$id}{$pick}{R};
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
  if ($spliced{$id}{$pick}{introns}) {warn "ERROR: $id: tRNAscan-SE called post-splicing intron\n"; next}
  $ct ++;
  push @{$n{$id}{notes}}, $spliced{$id}{$pick}{note} if $spliced{$id}{$pick}{note};
  @{$n{$id}}{                  qw/struct ac aa a_site score iso iso_score disc_fill disc_end introns splice_score trunc_start trunc_end/}
    = (@{$spliced{$id}{$pick}}{qw/str    ac aa a_site score iso iso_score disc_fill disc_end introns splice_score trunc_start trunc_end/});
  $n{$id}{a_site} += $n{$id}{ivs_R} - $n{$id}{ivs_L} if $n{$id}{a_site} > $n{$id}{ivs_L};
  $n{$id}{aa} = $anticodon2aa{$n{$id}{ac}} if $n{$id}{ac} eq 'CTA' and $anticodon2aa{CTA};  # Two alternate genetic codes
  $n{$id}{j_site} = $n{$id}{R} - $n{$id}{L} - 6; $n{$id}{j_site} -- if $n{$id}{aa} eq 'SeC';
  my ($L, $R, $gpi_score, $gpi_range) = @{(sort {$$b[2] <=> $$a[2]} @{$n{$id}{gpis}})[0]};  # Choose best if multiple gpI calls
  if ($n{$id}{ori} eq '+') {for ($L, $R) {$_ = $_-$n{$id}{L}+1}} else {for ($L, $R) {$_ = $n{$id}{R}-$_+1} ($L, $R) = ($R, $L)}
  @{$n{$id}}{qw/gpi_L gpi_R gpi_score gpi_range/} = ($L-1, $R, $gpi_score, $gpi_range);
  $final{$id} ++;
 }
 warn "$ct gpI intron-containing tRNAs survive\n";
 return '';
}

sub Aragorn1 {
 warn "Finding split tRNA gene calls\n";
 RunCommand("$binpath/aragorn1.1 -w -t -l -i -seq -o introns.aragorn $file", "introns.aragorn");  # aragorn1.1
 for (`cat introns.aragorn`) {
  if (/^>(\S+)/) {$dna = $1; $id = ''}
  elsif (/^T /) {$id = ''}
  elsif (/^([a-z\.]+)\s*$/) {next unless $id; $n{$id}{seq} .= $1; $n{$id}{aragorn} .= $_}
  elsif (/^TI/) {
   die "ERROR: Cannot parse aragorn1.1 $dir $_" unless /^TI\s+.*(.)\[([\-\d]+),(\d+)\]\s+(\d+)\s+\S+\s+(\d+)\s+(\d+)/;
   $serials{trna} ++; $id = "trna$serials{trna}";
   if ($6 < 100) {warn "$id: Intron too short\n"; $id = ''; next}
   %{$n{$id}} = (dna => $dna, dir => 1, ori => '+', L => $2, R => $3, stt => $2, stp => $3, a_site => $4+1, intron_coord => $5, intron_len => $6, start => 0);
   $n{$id}{L} = 1 if $n{$id}{L} < 1;
   my $mid = int(($2 + $3)/2);
   if ($1 eq 'c') {@{$n{$id}}{qw/dir ori stt stp/} = (-1, '-', $3, $2); $n{$id}{L} += 3} else {$n{$id}{R} -= 3}
   $tbed    .= join("\t", $n{$id}{dna}, $n{$id}{L}-1, $n{$id}{R}, $id, '', $n{$id}{ori}) . "\n";
   $tbedmid .= join("\t", $n{$id}{dna}, $mid-1,       $mid,       $id, '', $n{$id}{ori}) . "\n";
  } else {die "ERROR: Cannot parse aragorn1.1 $dir $_"}
 }
 return "No split trna calls" unless scalar(keys %n);
 warn scalar(keys %n), " original split trna calls\n";
  # Todo: Store original Tloop, swap in perfect one so trnascan-SE works, restore at end
 for my $id (keys %n) {  # Truncations
  $n{$id}{seq} = uc $n{$id}{seq};
  $n{$id}{seq} =~ s/(...)$//; $n{$id}{cca} = uc $1; $n{$id}{cca} =~ s/\./N/g; $n{$id}{cca} =~ s/N*$//;
  if ($n{$id}{seq} =~ s/^(\.+)//) {my $l = length $1; if ($n{$id}{dir} > 0) {$n{$id}{L} += $l; $n{$id}{trunc_start} = $l} else {$n{$id}{R} -= $l; $n{$id}{trunc_end} = $l}}
  if ($n{$id}{seq} =~ s/(\.+)$//) {my $l = length $1; if ($n{$id}{dir} > 0) {$n{$id}{R} -= $l; $n{$id}{trunc_end} = $l} else {$n{$id}{L} += $l; $n{$id}{trunc_start} = $l}}
  $n{$id}{seq} =~ s/\./N/g;
  $n{$id}{seq} .= 'N' x $n{$id}{trunc_end} if $n{$id}{trunc_end};  # Ns improve tSE calls on truncations
  warn "ERROR: $id: internal non-acgt\n" if $n{$id}{seq} =~ /\./;
 }
}

sub Presplice {
 warn "Re-evaluating acceptor stem\n";
 open OUT, ">prespliced.fa";  # Prepare candidates deleting all but 12 bp from each intron end, due to some incorrect Aragorn acceptor stem calls
 for my $id (sort keys %n) {
  #$n{$id}{seq} =~ s/\./N/g;
  my ($seq, $ac, $in, $len, $ttc, $torig) = @{$n{$id}}{qw/seq a_site intron_coord intron_len ttc torig/};
  #$n{$id}{start} = 0;
  my $pos = $ac;
  $pos -= $len if $ac > $in;
  my $presplice = substr($seq, 0, $pos) . lc(substr($seq, $pos, 12));
  $pos += $len;
  $presplice .= lc(substr($seq, $pos-12, 12)) . substr($seq, $pos);
  $presplice =~ s/$torig/$ttc/;
  print OUT ">$id.cca=$n{$id}{cca}.cut=0\n${presplice}CCA\n";
  #die "$seq $n{$id}{cca} $presplice\n";
  $presplice =~ s/(.)$//;  # Aragorn can add an 8th acc stem bp; try also shifting 3' end back one and see who scores best
  my $cca2 = $1 . $n{$id}{cca}; $cca2 =~ s/.$//;
  print OUT ">$id.cca=$cca2.cut=1\n${presplice}CCA\n";
 }
 close OUT;
 my ($start, $cca, $cut, $stop, $score, %p);
 my $call = "tRNAscan-SE -Q -q --thread 1 -f prespliced.struct $gencode -$tax --detail --brief prespliced.fa &> /dev/null" unless -f "prespliced.struct";
 RunCommand($call, "prespliced.struct");
 for (`cat prespliced.struct`) {
  next unless /^\S/;
  chomp;
  ($id, $cca, $cut, $start, $stop) = ($1, $2, $3, $4-1, $5) if /^(\S+)\.cca=(\S*)\.cut=(\d+)\.trna\d* \((\d*)-(\d*)/;
  $score = $1 if /Score: ([\d\.]+)/;
  next unless /^Str: (\S+)/;
  my $str = $1;
  next unless $str =~ /^(([^<]+<[^>]+[^<]+>)\.+<)/;  # Find anticodon stem closure
  next if $p{$id}{presplice_score} and $score <= $p{$id}{presplice_score};  # Choose better-scoring of two cca positions tested
  my ($acL, $acR) = (length($2), length($1));
  $acL ++, $acR -- if /<\.[\.>]{4}\.{20}/;  # Handles non-standard ac-closing base pair as in Aba1905
  my ($new_len, $new_a) = ($n{$id}{intron_len} + $acR - $acL - 32, $acL + 4);
  #warn "Presplice $id: old intron length $n{$id}{intron_len} now $new_len; a_site $n{$id}{a_site} now $new_a; acL=$acL, acR=$acR start=$start cut=$cut score=$score\n";
  @{$p{$id}}{qw/intron_len a_site presplice_score start cca cut/} = ($new_len, $new_a, $score, $start, $cca, $cut);   
 }
 for my $id (sort keys %p) {
  @{$n{$id}}{qw/intron_len a_site presplice_score start cca/} = @{$p{$id}}{qw/intron_len a_site presplice_score start cca/};
  #warn "$n{$id}{seq} $p{$id}{cut}\n";
  $n{$id}{seq} =~ s/.$// if $p{$id}{cut};
  $n{$id}{seq} =~ s/^.{$p{$id}{start}}// if $p{$id}{start};
  #warn "$n{$id}{seq} $p{$id}{cut}\n";
  $n{$id}{intron_coord} -= $p{$id}{start} if $p{$id}{start};
  if ($id =~ /(.*)\.(\dprime)/) {
   my ($orig, $ghost) = ($1, $2); $ghost = 'both' if $ghost eq '1prime';
   push @{$n{$orig}{notes}}, "ghost $ghost presplice score $p{$id}{presplice_score} (current $n{$orig}{presplice_score})";
   next if $p{$id}{presplice_score} <= $n{$orig}{presplice_score} or $n{$orig}{ghost};
   push @{$n{$orig}{notes}}, "correct to $ghost";
   $n{$orig}{ghost} = $ghost;
   $n{$id}{start} += $start;
   print "$id $orig $ghost $start\n";   
   @{$n{$orig}}{qw/start a_site intron_coord intron_len seq cca ghost presplice_score/} = (@{$n{$id}}{qw/start a_site intron_coord intron_len seq cca/}, $ghost, $score); 
  }
 }
 for (keys %n) {delete $n{$_} if /prime|both/}
}

sub Tloop {
 warn "Examining T-stem\n";
 open OUT, ">preT.fa";  # Prepare candidates deleting all but 12 bp from each intron end, due to some incorrect Aragorn acceptor stem calls
 for my $id (sort keys %n) {
  my ($seq, $ac, $in, $len) = @{$n{$id}}{qw/seq a_site intron_coord intron_len/};
  my $pos = $ac;
  $pos -= $len if $ac > $in;
  my $presplice = substr($seq, 0, $pos) . lc(substr($seq, $pos, 12));
  $pos += $len;
  $presplice .= lc(substr($seq, $pos-12, 12)) . substr($seq, $pos);
  print OUT ">$id\n${presplice}$n{$id}{cca}\n";
 }
 close OUT;
 my $call = "$binpath/aragorn1.2.40 -w -br -t -l -i -e -o preT.aragorn preT.fa";
 RunCommand($call, "preT.aragorn");
 my ($seq, $ttc, $str, $torig, $ttcold);
 for (`cat preT.aragorn`) {
  chomp;
  next if /found/;
  if (/^>(\S+)/) {$id = $1}
  elsif ($id and /^[a-z]/) {$seq = uc($_);}
  next unless $id and /\({2}/;
  $seq =~ s/\./N/g;
  /(.*[^\(])(\(+t+\)+)/;
  $str = $2;
  $torig = $ttc = substr $seq, length($1), length($2);
  $str =~ /([^t]+)/;
  $ttcold = substr $ttc, length($1), 3;
  substr $ttc, length($1), 3, 'TTC';
  #warn "$torig $ttc $ttcold\n";
  @{$n{$id}}{qw/ttcold torig ttc/} = ($ttcold, $torig, $ttc);
 }
}

sub Aragorn2 {
 my $call = "$binpath/aragorn1.2.40 -w -br -t -l -i -e -o spliced.aragorn spliced.fa";
 #>trna2.4Ug
 #1 gene found
 #1   tRNA-Asn                        [1,75]      113.6   33      (gtt)
 #tcccctgtagctcaatggtagagcggctcgctgttaacgagaaggttcgtggtccgagtccacgcgggggagcca
 #(((((((ss((((ddddddd))))s.((((ccAAAcc)))).vvvvv(((((ttttttt))))))))))))
 
 RunCommand($call, "spliced.aragorn");
 my (%t, $splice);
 for (`cat spliced.aragorn`) {
  chomp;
  next if /found/;
  if (/^>(\S+)\.(\S+)/) {($id, $splice) = ($1, $2); $id = '' if $spliced{$id} and $spliced{$id}{$splice};}  # Only tRNAs missed by tRNAscan-SE
  elsif ($id and /^\d+\s+tRNA-(\S+)\s+\[(\d+),(\d+)\]\s+(\S+)\s+(\d+)\s+\(([a-z]+)/) {
   @{$spliced{$id}{$splice}}{qw/dir aa L R score a_site ac iso iso_score IPD note/} = (1, $1, $2, $3, $4, 1+$5, uc($6), 'None', $4, 0, 'Aragorn2');
  }
  elsif ($id and /^[a-z]/) {@{$spliced{$id}{$splice}}{qw/seq mat/} = (uc($_), uc($_));}
  elsif ($id and /\({2}/) {
   s/[^\(\)]/./g; s/\(/</g; s/\)/>/g;
   $spliced{$id}{$splice}{seq} =~ s/$n{$id}{ttc}/$n{$id}{torig}/;
   $spliced{$id}{$splice}{str} = $_ . '.' x (length($spliced{$id}{$splice}{seq}) - length($_));
   my $retee = Acc($splice, \%spliced); %spliced = %{$retee};
  }
 }
}

sub Scan {
 my (%t, $splice, $ttc);
 my $call = "tRNAscan-SE -Q -q --thread 1 -s spliced.iso -f spliced.struct -o spliced.out $gencode -$tax --detail --brief spliced.fa &> /dev/null";
 RunCommand($call, "spliced.iso");
 for (`cat spliced.struct`) {
  next unless /^\S/;
  chomp;
  if (/^(\S+)\.trna\d* \(\d*-\d*/) {$id = $1; $id =~ s/\.([^\.]+)$//; $splice = $1}
  elsif (/^Type: \S+\s+Anticodon: \S+ at (\d+)\S+ \S+\s+Score: \S+/) {$t{$id}{$splice}{a_site} = $1 + 1}
  elsif (/^Possible intron: (\d+)-(\d+)/) {push @{$t{$id}{$splice}{intLs}}, $1; push @{$t{$id}{$splice}{intRs}}, $2}
  elsif (/^Possible/) {}
  elsif (s/^Pre: (\S+).*/$1/) {while (s/\[([^\]]*)\]/lc($1)/e) {} $t{$id}{$splice}{seq} = $_}  # Lower-case (multiple) introns
  elsif (/^Seq: (\S+)/) {$t{$id}{$splice}{mat} = $1}
  elsif (/^Str: (\S+)/) {$t{$id}{$splice}{str} = $1}
  else {warn "ERROR: Cannot parse tRNAscan-SE line $_\n"}
 }
 for (`cat spliced.iso`) {
  chomp; my @f = split "\t";
  my ($id, $nom) = (shift(@f), shift(@f)); $id =~ s/\.trna\d*$//; $id =~ s/\.([^\.]+)$//; my $splice = $1;
  for (@isos) {$t{$id}{$splice}{isos}{$_} = shift @f; $t{$id}{$splice}{isos}{$_} = 0 if $t{$id}{$splice}{isos}{$_} < 0}
 }
 for (`cat spliced.out`) {
  #fake            13      2000494 2000565 Undet   NNN     0       0       20.3    Arg     20.2    pseudo,trunc_start:32
  chomp; s/\s+\d+//; my @f = split "\t"; $id = $f[0]; $id =~ s/\.([^\.]+)$//; $splice = $1; $f[10] = '' unless $f[10];
  die "Cannot parse dna $_\n" unless /^(\S+)\t/;
  @{$t{$id}{$splice}}{qw/dna dir L R aa ac score iso iso_score note IPD/} = ($1, 1, @f[1..4], @f[7..10], 0);
  if ($t{$id}{$splice}{L} == 0) {$t{$id}{$splice}{L} = 1; warn "ERROR: Hard-trunc unresolved $_\n" unless $t{$id}{$splice}{str} =~ s/^\.>/>/}  # Correctable error noted on truncated Ype23.trna11 et al
  if ($f[1] > $f[2]) {warn "ERROR: tRNA call in wrong orientation\n"; next}
  $t{$id}{$splice}{iso_score} += $1 if $f[10] =~ /IPD:([0-9\.\-]+)/;
  my $aa = $t{$id}{$splice}{aa};
  $t{$id}{$splice}{aa} = $anticodon2aa{$t{$id}{$splice}{ac}} if $t{$id}{$splice}{ac} eq 'CTA' and $anticodon2aa{CTA};  # Two alternate genetic codes
  @{$t{$id}{$splice}}{qw/aa iso_score/} = ('fMet', $t{$id}{$splice}{isos}{fMet}) if $t{$id}{$splice}{isos}{fMet} > $t{$id}{$splice}{isos}{Met} and $aa eq 'Met';  # best CAT if other iso wins
  @{$t{$id}{$splice}}{qw/aa iso_score/} = ('Ile2', $t{$id}{$splice}{isos}{Ile2}) if $t{$id}{$splice}{isos}{Ile2} > $t{$id}{$splice}{isos}{Met} and $t{$id}{$splice}{isos}{Ile2} > $t{$id}{$splice}{isos}{fMet} and $aa eq 'Met';
  $t{$id}{$splice}{note} .= ",original call Met" if $aa eq 'Met' and $t{$id}{$splice}{aa} ne 'Met';
  if ($t{$id}{$splice}{seq}) {
   my $test = $t{$id}{$splice}{seq}; $test =~ s/[a-z]+//g; warn "ERROR: spliced pre ne mature $id\n" if $test ne $t{$id}{$splice}{mat};
  } else {$t{$id}{$splice}{seq} = $t{$id}{$splice}{mat}}  # Seq is pre-sequence with introns, mat is mature spliced
  #$t{$id}{$splice}{seq} =~ s/$n{$id}{ttc}/$n{$id}{torig}/;
  my $retee = Acc($splice, \%t); %t = %{$retee}; #next;
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

sub Acc {
 my ($splice, %t) = ($_[0], %{$_[1]});
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
 return \%t;
}

sub Introns {
 #my $cmd = "cmsearch --cpu 0 --tblout gpi.tbl $lib/cm/Intron_gpI.cm $file &> /dev/null";
 my $cmd = "cmsearch --cpu 0 --tblout gpi.tbl $lib/cm/gpi_bact_trna.cm $file &> /dev/null";
 RunCommand($cmd, 'gpi.tbl');
 warn "Can't find completed gpi.tbl\n" unless -s 'gpi.tbl';
 my @gpi;
 for (`cat gpi.tbl`) {
  next if /^#/;  my @f = split /\s+/;
  ##target_name         accession query_name           accession mdl mdl_from   mdl_to seq_from   seq_to strand trunc pass   gc  bias  score   E-value inc description_of_target
  #Intron_gpI           RF00028   AM747720.1           -          cm       49      160  1222644  1222744      +    no    1 0.62   2.6   28.9   4.5e-06 !   -
  next if $f[14] < 18;  # CM score cutoff for Intron_gpI; cm "trusted cutoff" is 40.01
  $f[10] =~ s/\'/prime/g; $f[10] =~ s/5prime\&3prime/both/; my $range = "$f[10]$f[5]..$f[6]";
  ($f[7], $f[8]) = ($f[8], $f[7]) if $f[7] > $f[8];
  $serials{gpi} ++;
  push @gpi, [$f[0], $f[7]-1, $f[8], "gpi$serials{gpi}", $f[14], $f[9], $range];
 }
 return unless @gpi;
 open OUT, ">gpi.bed"; for (sort {$$a[0] cmp $$b[0] || $$a[1] <=> $$b[1]} @gpi) {print OUT join("\t", @{$_}), "\n"} close OUT;
}

sub LoadCode {
 %anticodon2aa =      qw/TGA Ser GGA Ser CGA Ser AGA Ser GAA Phe AAA Phe TAA Leu CAA Leu GTA Tyr ATA Tyr GCA Cys ACA Cys CCA Trp 
 TAG Leu GAG Leu CAG Leu AAG Leu TGG Pro GGG Pro CGG Pro AGG Pro GTG His ATG His TTG Gln CTG Gln TCG Arg GCG Arg CCG Arg ACG Arg 
 TAT Ile GAT Ile AAT Ile CAT Met TGT Thr GGT Thr CGT Thr AGT Thr GTT Asn ATT Asn TTT Lys CTT Lys GCT Ser ACT Ser TCT Arg CCT Arg 
 TAC Val GAC Val CAC Val AAC Val TGC Ala GGC Ala CGC Ala AGC Ala GTC Asp ATC Asp TTC Glu CTC Glu TCC Gly GCC Gly CCC Gly ACC Gly/;
 %one_letter = (qw/Ala A Arg R Asn N Asp D Cys C Gln Q Glu E Gly G His H Ile I Ile2 I Leu L Lys K Met M Phe F Pro P SeC U Ser S Thr T Trp W Tyr Y Val V fMet B iMet B Sup * Undet X None X/);
}

sub Ghost {
 warn "Checking for ghost exons\n";
 open OUT, ">ghost.fa";  # Prepare candidates deleting all but 12 bp from each intron end, due to some incorrect Aragorn acceptor stem calls
 my ($ghost, %seqs, %starts);
 for my $id (sort keys %n) {
  unless ($n{$id}{gpis}) {warn "$id: ghost No Intron_gpI hit\n"; delete $n{$id}; next}
  my ($seq, $in, $len) = @{$n{$id}}{qw/seq intron_coord intron_len/};
  $starts{"$id.5prime"} = $starts{"$id.1prime"} = $in-1;
  $starts{"$id.3prime"} = 0;
  $seqs{"$id.5prime"} = substr($seq, $in-1) . $n{$id}{cca};
  $seqs{"$id.3prime"} = substr($seq, 0, $in+$len-1);
  $seqs{"$id.1prime"} = substr($seq, $in-1, $len);
  print OUT ">$id.5prime\n", $seqs{"$id.5prime"}, "\n",
            ">$id.3prime\n", $seqs{"$id.3prime"}, "\n",
            ">$id.1prime\n", $seqs{"$id.1prime"}, "\n"
 }
 close OUT;
 my $call = "$binpath/aragorn1.1 -w -t -l -i -seq -o ghost.aragorn ghost.fa" unless -f "ghost.aragorn";  # aragorn1.1
 RunCommand($call, "ghost.aragorn");
 for (`cat ghost.aragorn`) {
  $id = $1 if /^>(\S+)/;
  next unless /^TI/;
  die "ERROR: Cannot parse aragorn1.1 ghost.aragorn $_" unless /^TI\s+.*(.)\[([\-\d]+),(\d+)\]\s+(\d+)\s+\S+\s+(\d+)\s+(\d+)/;
  next if $1 eq 'c' or $6 < 50;  # Wrong orientation or intron too short
  @{$n{$id}}{qw/start a_site intron_coord intron_len seq/} = ($starts{$id}+$2-1, $4+1, $5, $6, uc substr($seqs{$id}, $2-1, $3-$2+1));
  $n{$id}{seq} =~ s/(...)$//; $n{$id}{cca} = $1; $n{$id}{cca} =~ s/\./N/g; $n{$id}{cca} =~ s/N*$//;
  if ($n{$id}{seq} =~ s/^(\.+)//) {my $l = length $1; if ($n{$id}{dir} > 0) {$n{$id}{L} += $l; $n{$id}{trunc_start} = $l} else {$n{$id}{R} -= $l; $n{$id}{trunc_end} = $l}}
  if ($n{$id}{seq} =~ s/(\.+)$//) {my $l = length $1; if ($n{$id}{dir} > 0) {$n{$id}{R} -= $l; $n{$id}{trunc_end} = $l} else {$n{$id}{L} += $l; $n{$id}{trunc_start} = $l}}
  $n{$id}{seq} =~ s/\./N/g;
  $n{$id}{seq} .= 'N' x $n{$id}{trunc_end} if $n{$id}{trunc_end};  # Ns improve tSE calls on truncations
 }
}

sub RunCommand {
 my ($command, $checkfile) = @_;
 if ($checkfile and -e $checkfile) {print "Skipping command: $command\n"; return}
 print "Running command: $command\n";
 my $out = system($command);
 if ($out) {print "Command '$command' failed with error message $out\n"; exit}
 else {print "Command '$command' succeeded\n"}
}

