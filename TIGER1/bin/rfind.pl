#! /usr/bin/perl
use strict; use warnings;
use IPC::Run3 qw/run3/;
use File::Spec;
# todo: exclude curated pseudos, place in alignments, opp orient hard mask, insist on cod/acc hit

die "Usage: $0 db infasta outdirectory\n" unless @ARGV == 3;
my ($db, $infile, $outdir, $binpath) = ($ARGV[0], $ARGV[1], $ARGV[2], $0); for ($infile, $outdir, $binpath) {$_ = File::Spec->rel2abs($_)}
$binpath =~ s/\/([^\/]+)$//; my $scriptname = $1;
my $lib = $binpath; $lib =~ s/[^\/]*$/db/; my $dbfiles = "$lib/$db/$db";
die "No blast nucleotide database $db available in $dbfiles\n" unless -f "$dbfiles.nin";
my (%meta, %inseqs, %insizes, $id, $seq, @regions, %serials, %t, %seen, %gencode);
my ($round, $flank) = (1, 250);

Initial();
#exit;
chdir $outdir;
open FINAL, ">rfind.gff" or die "Cannot open outfile $db.rfind\n";
mkdir 'rfind'; chdir 'rfind';
#my @blastCmd = (qw/blastn -soft_masking false -lcase_masking -db/, $dbfiles, qw/-outfmt 6/);  # Run3 wants an array of its command words!
my @blastCmd = (qw/blastn -dust no -soft_masking false -lcase_masking -db/, $dbfiles, qw/-outfmt 6/);  # Run3 wants an array of its command words!
my @mergeCmd = (qw/bedtools merge -s -c/, '4,6', qw/-o collapse -i stdin/);
#warn scalar(keys %inseqs), " inseqs, first of size ", length($inseqs{(keys %inseqs)[0]}), "; db=$dbfiles\n";
for $round (1..10) {
 unless (-f "$db.$round.merge") {
  my ($bed, $gnm, @hits, @beds, @merges) = ('', '');
  for (keys %inseqs) {$gnm .= ">$_\n$inseqs{$_}\n"}
  #open OUT, ">genome.$round.fa"; print OUT $gnm; close OUT;
  my $cmd = join(' ', @blastCmd); warn "$cmd\n"; 
  eval { run3 \@blastCmd, \$gnm, \@hits };
  if    ( $@        ) { warn "Error: $@\n";                     }
  elsif ( $? & 0x7F ) { warn "Killed by signal ".( $? & 0x7F ) . "\n"; }
  elsif ( $? >> 8   ) { warn "Exited with error ".( $? >> 8 ) . "\n";  }
  else                { warn "Completed successfully\n";          }
  #run3(\@blastCmd, \$gnm, \@hits);
  warn scalar(@hits), " hits in round $round\n";
  #open OUT, ">hits.$round.blast"; print OUT join('', @hits); close OUT;
  for (@hits) {
   chomp; my @f = split "\t"; my $ori = '+'; ($f[8], $f[9], $ori) = ($f[9], $f[8], '-') if $f[8] > $f[9]; for (@f[6,8]) {$_ --}
   my ($len, $form) = @{$meta{$f[1]}}{qw/len form/};
   warn "no len for $infile $_\n" unless $len;
   my $score = sprintf "%.2f", $f[2]*($f[9]-$f[8])/$len;  # Use 2 signif digits
   #warn "$f[2],$f[9],$f[8],$len  $_\n" if $score == 121.33;
   push @beds, [$f[0], $f[6], $f[7], join(':', @f[1,6,7,8,9,2,11], $score, $len, $form), '', $ori];
  }
  for (sort {$$a[0] cmp $$b[0] || $$a[1] <=> $$b[1]} @beds) {$bed .= join("\t", @{$_}) . "\n";}
  open OUT, ">$round.bed"; print OUT $bed; close OUT;
  eval { run3 \@mergeCmd, \$bed, \@merges };
  if    ( $@        ) { warn "Error: $@\n";                     }
  elsif ( $? & 0x7F ) { warn "Killed by signal ".( $? & 0x7F ) . "\n"; }
  elsif ( $? >> 8   ) { warn "Exited with error ".( $? >> 8 ) . "\n";  }
  else                { warn "Completed successfully\n";          }
  #run3(\@mergeCmd, \$bed, \@merges) if @hits;
  open OUT, ">$db.$round.merge"; print OUT join('', @merges); close OUT;
 }
 my $ct = 0;
 for (`cat $db.$round.merge`) {  # Choose best hit for each merged region, mask genome at that hit
  chomp; my @f = split "\t";
  #warn "f0=$f[0]\n\nf1=$f[1]\n\nf2=$f[2]\n\nf3=$f[3]\n\nf4=$f[4]\n";
  my (@hits, @modelparts, $score, $model);
  for (split ',', $f[3]) {push @hits, [split(':', $_)]; print join(',', @{$hits[-1]}), " hit without all data\n" unless $hits[-1][6]}  # Array of arrays [0]model [1]gnmL [2]gnmR [3]modelL [4]modelR [5]pctMatch [6]bits [7]score [8]lenModel [9]form
  for my $hit (sort {$$b[6] <=> $$a[6]} @hits) {if ($$hit[9] eq 'int' and $$hit[4] > 0.9 * $$hit[8]) {$model = $$hit[0]; last}}  # Int-bearing tmRNA bitscore can be lower for int model than std
  for my $hit (sort {$$b[6] <=> $$a[6] || $$b[7] <=> $$a[7]} @hits) {  # Tile hits to same model within region
   $model = $$hit[0] unless $model;
   next unless $$hit[0] eq $model;
   for my $prev (@modelparts) {
    if    ($$hit[3] < $$prev[3]) {$$hit[4] = $$prev[3] if $$hit[4]>$$prev[3]}
    elsif ($$hit[4] > $$prev[4]) {$$hit[3] = $$prev[4] if $$hit[3]<$$prev[4]}
    else {$$hit[4] = 0}
   }
   next unless $$hit[4];
   push @modelparts, [@{$hit}, "$$hit[3]-$$hit[4]:$$hit[5]:$$hit[6]"];
   $score += ($$hit[4]-$$hit[3])*$$hit[5];  # pctMatch * trimmedLength
  }
  my $L  = (sort {$a <=> $b} map $_->[1], @modelparts)[0];
  my $R  = (sort {$a <=> $b} map $_->[2], @modelparts)[-1];
  my $mL = (sort {$a <=> $b} map $_->[3], @modelparts)[0];
  my $mR = (sort {$a <=> $b} map $_->[4], @modelparts)[-1];
  die "no length for model $model $L $R\n" unless $meta{$model}{len};
  $score = sprintf "%.2f", $score/$meta{$model}{len};
  my $hitsum = join(',', map($_->[10], @modelparts));
  my ($dna, $len) = ($f[0], $R-$L);
  #warn "model=$model, score=$score, L=$L, R=$R, mL=$mL, mR=$mR, len=$len hitsum=$hitsum\n";
  $f[4] =~ /^(.)/; my $ori = $1;
  push @regions, {dna=>$dna, L=>$L, R=>$R, mL=>$mL, mR=> $mR, ori=> $ori, model=>$model, score=>$score, hitsum=>$hitsum} unless $seen{"$dna.$L.$R"};
  $seen{"$dna.$L.$R"} ++;
  $ct ++;
  $seq = $inseqs{$dna}; $inseqs{$dna} = (substr($seq, 0, $L) . lc(substr($seq, $L, $len)) . substr($seq, $L+$len));  # Mask hit for next round of blast
 }
 warn "$ct new regions after Round $round\n";
 if (scalar @regions >= 10) {@regions = sort {$$b{score} <=> $$a{score}} @regions; splice @regions, 10; warn "Truncating to 10\n"; last}
 #die scalar(@regions), " regions\n";
 last unless $ct;
}
@regions = sort {$$a{dna} cmp $$b{dna} || $$a{L} <=> $$b{L}} @regions;
my $ct = 0;
for my $f (@regions) {
 my ($dna, $L, $R, $model, $mL, $mR, $ori, $gnmOverL, $gnmOverR) = ($$f{dna}, $$f{L}, $$f{R}, $$f{model}, $$f{mL}, $$f{mR}, $$f{ori}, 0, 0);
 if ($meta{$model}{ivsL} and $mL >= $meta{$model}{ivsL}-5 and $mR <= $meta{$model}{ivsR}+5) {next}  # Skip hits within IVS (such as RPEs) 
 my $core = uc(substr $inseqs{$dna}, $L, $R-$L);
 my ($posC, $lenC) = ($L-$mL, $R + $meta{$model}{len}-$mR - $L+$mL);  # Len calculated same in either orientn
 $posC = $L-$meta{$model}{len}+$mR if $ori eq '-';
 my ($left, $right, $posL, $lenL, $posR, $lenR) = ('', '', $posC-$flank, $flank, $posC+$lenC, $flank);
 if ($posC < 0) {$gnmOverL = -1 * $posC; $lenC += $posC; $posC = 0}  # Truncated at left end of contig
 else {if ($posL < 0) {$lenL += $posL; $posL = 0} $left = lc(substr($inseqs{$dna}, $posL, $lenL))}
 if ($posC+$lenC > $insizes{$dna}) {$gnmOverR = $posC+$lenC-$insizes{$dna}; $lenC -= $gnmOverR}
 else {if ($posR+$lenR > $insizes{$dna}) {$lenR -= $posR+$lenR-$insizes{$dna}} $right = lc(substr($inseqs{$dna}, $posR, $lenR))}
 $seq = uc(substr($inseqs{$dna}, $posC, $lenC));
 my ($gnmTruncUp, $gnmTruncDn) = ($gnmOverL, $gnmOverR);
 #my @core = ($mL, $mR);
 if ($ori eq '-') {for ($left, $seq, $right, $core) {$_ = Revcomp($_)} ($left, $right) = ($right, $left); ($gnmTruncUp, $gnmTruncDn) = ($gnmOverR, $gnmOverL)}
 my $tag = '';
 if ($meta{$model}{tagR} and $meta{$model}{tagR} <= $mR+3 and $meta{$model}{tagL} >= $mL) {
  my ($L, $len) = ($meta{$model}{tagL}-$gnmTruncUp, $meta{$model}{tagR}-$meta{$model}{tagL}+3);
  if (length $seq >= $L+$len) {my $orf = substr $seq, $L, $len; $tag = Translate($orf)}
 }
 #print "posC=$posC posR=$posR lenC=$lenC lenL=$lenL modellen=$meta{$model}{len} $seq\tmodeltagR=$meta{$model}{tagR} mR=$mR and modeltagL=$meta{$model}{tagL} mL=$mL gnmTruncUp=$gnmTruncUp orf=$orf tag=$tag\n";
 $serials{final} ++;
 print FINAL join("\t", $dna, qw/rfind tmRNA/, $L+1, $R, $$f{score}, $ori, '.', "ID=rfind.$serials{final};") .
  "model=$model;modelL=$mL;modelR=$mR;modelLen=$meta{$model}{len};gnmTruncUp=$gnmTruncUp;gnmTruncDn=$gnmTruncDn;hitsum=$$f{hitsum};tag=$tag;";
 for (qw/tag ivs int/) {
  next unless $meta{$model}{$_ .'R'};
  my ($L, $R) = ($meta{$model}{$_ .'L'}-$gnmTruncUp, $meta{$model}{$_ .'R'}-$gnmTruncUp);
  next if $L < 0 or $R > length($seq);
  print FINAL "${_}L=$L;${_}R=$R;";
  my $fill = lc substr($seq, $L, $R-$L); substr($seq, $L, $R-$L, $fill);
 }
 $seq =~ /$core/i; 
 $seq = substr($seq, 0, $-[0]) . ',' . substr($seq, $-[0], $+[0]-$-[0]) . ',' . substr($seq, $+[0]); 
 for (qw/cca form/) {if ($meta{$model}{$_}) {print FINAL "$_=$meta{$model}{$_};"} else {print FINAL "$_=;"}}
 print FINAL "seq=$left,$seq,$right;\n";
 $ct ++;
}
warn "$ct regions reported for $outdir\n";

sub Initial {
 LoadCode();
 for (`cat $infile`) {if (s/>(\S+)\s*//) {chomp; $id = $1} else {chomp; $inseqs{$id} .= $_}}
 for (keys %inseqs) {$insizes{$_} = length $inseqs{$_}; $inseqs{$_} = uc $inseqs{$_}}
 for (`cat $dbfiles.fa`) {
  if (s/^>(\S+)//) {$id = $1; while (s/^ ([^=]+)=([^;]+);/ /) {$meta{$id}{$1} = $2} die "$id\n" unless $meta{$id}{form}; next}  #die $_ unless /form=([^;]+);/; $meta{$id}{form} = $1; next}
  chomp;
  @{$meta{$id}}{qw/seq len/} = ($_, length($_));
  if ($meta{$id}{form} eq 'std' and $meta{$id}{tag} =~ /[a-z]/) {}
  elsif ($meta{$id}{form} eq 'std') {die "$id $_\n" unless /^[A-Z]+([a-z]+)[A-Z]+[a-z]{3}$/; @{$meta{$id}}{qw/tagL tagR/} = ($-[1], $+[1])}
  if ($meta{$id}{form} eq 'prm' and $meta{$id}{tag} =~ /[a-z]/) {die "$id\n" unless /^[A-Z]+([a-z]+)[A-Z]+$/; @{$meta{$id}}{qw/ivsL ivsR/} = ($-[1], $+[1])}
  elsif ($meta{$id}{form} eq 'prm') {die "$id\n" unless /^[A-Z]+([a-z]+)[A-Z]+([a-z]+)[A-Z]+$/; @{$meta{$id}}{qw/ivsL ivsR tagL tagR/} = ($-[1], $+[1], $-[2], $+[2])}
  if ($meta{$id}{form} eq 'int') {die "$id\n" unless /^[A-Z]+([a-z]+)[A-Z]+([a-z]+)[A-Z]+[a-z]{3}$/; @{$meta{$id}}{qw/tagL tagR intL intR/} = ($-[1], $+[1], $-[2], $+[2])}
  if ($meta{$id}{form} eq 'short') {die "$id\n" unless /^[A-Z]+$/}
  if ($meta{$id}{tagR}) {my $orf = substr $_, $meta{$id}{tagL}, $meta{$id}{tagR}-$meta{$id}{tagL}; my $tag = Translate(uc $orf); die "$meta{$id}{tag} ne $tag for $id\n" unless $meta{$id}{tag} eq $tag}
 }
}

sub Revcomp {my $ret = reverse $_[0]; $ret =~ tr/ACGTacgt/TGCAtgca/; return $ret}

sub LoadCode {
 %gencode = qw(GCT A GCC A GCA A GCG A TGT C TGC C GAT D GAC D GAA E GAG E TTT F TTC F GGT G GGC G GGA G GGG G CAT H CAC H ATT I ATC I
   ATA I AAA K AAG K CTT L CTC L CTA L CTG L TTA L TTG L ATG M AAT N AAC N CCT P CCC P CCA P CCG P CAA Q CAG Q CGT R CGC R CGA R CGG R
   AGA R AGG R TCT S TCC S TCA S TCG S AGT S AGC S ACT T ACC T ACA T ACG T GTT V GTC V GTA V GTG V TGG W TAT Y TAC Y TAA * TAG * TGA *);
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

