#! /usr/bin/perl
use strict; use warnings;

die "Usage: $0 infasta outdir [AMG taxa]\n" if @ARGV < 1;
use File::Spec;
my ($infile, $outdir, $binpath) = ($ARGV[0], $ARGV[1], $0); for ($infile, $outdir, $binpath) {$_ = File::Spec->rel2abs($_)}
$binpath =~ s/\/([^\/]+)$//; my $scriptname = $1;
my $lib = $binpath; $lib =~ s/[^\/]*$/lib/;
my ($tax, $gencode, %anticodon2aa, %one_letter, %final, %t, $id, %serials, $dna, $collect) = ('B', '');

LoadCode();
$tax = $ARGV[2] if $ARGV[2];
if ($tax eq 'M') {($tax, $gencode, $anticodon2aa{CTA}) = ('B', "-g $lib/gencode/gcode.othmito", 'Trp')}  # ncbi gencode=4
if ($tax eq 'G') {($tax, $gencode, $anticodon2aa{CTA}) = ('B', "-g $lib/gencode/gcode.gracil" , 'Gly')}  # ncbi gencode=25
chdir $outdir;
open FINAL, ">trna.gff" or die "Can't open outfile $outdir/trna.gff: $!\n";
mkdir 'trna';
chdir 'trna';

Trna();
for my $dna (sort keys %final) {
 for my $L (sort {$a <=> $b} keys %{$final{$dna}}) {
  for my $R (sort {$b <=> $a} keys %{$final{$dna}{$L}}) {
   my $id = $final{$dna}{$L}{$R};
   $serials{final} ++;
   print FINAL join("\t", $t{$id}{dna}, qw/tRNAscan-SE2.0 tRNA/, @{$t{$id}}{qw/L R score ori/}, '.', "ID=trna.$serials{final};");
   for (qw/aa ac a_site j_site disc_fill disc_end iso iso_score introns trunc_start trunc_end note struct seq cca/)
   {if (defined $t{$id}{$_}) {print FINAL "$_=$t{$id}{$_};"} else {print FINAL "$_=;"}}
   print FINAL "\n";
  }
 }
}
close FINAL;

# SUBROUTINES
sub Trna {
 my $call = "tRNAscan-SE -Q -q --thread 1 -f trna.struct -o trna.out $gencode -$tax --detail --brief $infile &> /dev/null";  # Also want -s base.iso, but not strictly necessary
 RunCommand($call, 'trna.out');
 for (`cat trna.struct`) {
  next unless /^\S/;
  chomp;
  if (/^(\S+.trna\d*) \(\d*-\d*/) {$id = $1}
  elsif (/^Type: \S+\s+Anticodon: \S+ at (\d+)\S+ \S+\s+Score: \S+/) {$t{$id}{a_site_post} = $t{$id}{a_site} = $1 + 1}
  elsif (s/^Pre: (\S+).*/$1/) {while (s/\[([^\]]*)\]/lc($1)/e) {} $t{$id}{seq} = $_}  # Lower-case (multiple) introns
  elsif (/^Possible intron: (\d+)-(\d+)/) {push @{$t{$id}{intLs}}, $1; push @{$t{$id}{intRs}}, $2}
  elsif (/^Possible/) {}
  elsif (/^Seq: (\S+)/) {$t{$id}{mat} = uc $1}
  elsif (/^Str: (\S+)/) {$t{$id}{struct} = $1}
  else {warn "ERROR: Cannot parse tRNAscan-SE line $_\n"}
 }
 warn scalar(keys(%t)), " calls\n";
 for (`cat trna.out`) {
  #fake            13      2000494 2000565 Undet   NNN     0       0       20.3    Arg     20.2    pseudo,trunc_start:32
  chomp;
  die "Cannot parse dna $_\n" unless s/^(\S+)\s+/${1}.trna/; my @f = split /\s+/; $id = $f[0]; $f[10] = '' unless $f[10];
  @{$t{$id}}{qw/dna dir ori L R aa ac score iso iso_score note/} = ($1, 1, '+', @f[1..4], @f[7..10]);
  if ($t{$id}{L} == 0) {$t{$id}{L} = 1; warn "ERROR: Hard-trunc unresolved $_\n" unless $t{$id}{struct} =~ s/^\.>/>/}  # Correctable error noted on truncated Ype23.trna11 et al
  @{$t{$id}}{qw/dir ori L R/} = (-1, '-', @f[2,1]) if $f[1] > $f[2];
  $t{$id}{aa} = $anticodon2aa{$t{$id}{ac}} if $t{$id}{ac} eq 'CTA' and $anticodon2aa{CTA};  # Two alternate genetic codes
  if ($t{$id}{seq}) {
   my $test = $t{$id}{seq}; $test =~ s/[a-z]+//g; warn "ERROR: spliced pre ne mature $id\n" if $test ne $t{$id}{mat};
  } else {$t{$id}{seq} = $t{$id}{mat}}  # Seq is pre-sequence with introns, mat is mature spliced
  if ($t{$id}{note} =~ s/,*trunc_start:(\d+)//) {$t{$id}{trunc_start} = $1}  
  if ($t{$id}{note} =~ s/,*trunc_end:(\d+)//  ) {$t{$id}{trunc_end} = $1}  

  # TRIM CCA
  $t{$id}{struct} =~ /.*</;  # Find last position of acceptor stem
  $t{$id}{disc_end} = $+[0] +1;
  my ($expAcc, $len, $toy, @loops) = (7, length($t{$id}{mat}), $t{$id}{struct}); $expAcc ++ if $t{$id}{aa} eq 'SeC';
  warn "0 $toy $t{$id}{disc_end}\n";
  while ($toy =~ s/>\.+<([\.<X]+)$/X$1/) {  # Replace loops, marking positions of last stem ntds
   unshift @loops, $-[1];
   warn "1 $toy @loops\n";
   while ($toy =~ s/>\.*X(\.*<)/X/) {$loops[0] += length $1; warn "2 $toy @loops\n";}
  }
  my $accLen = $toy =~ s/>/>/g; $accLen = 0 unless $accLen;  #  Count acceptor stem bp after removing loops
  $t{$id}{j_site} = $loops[-1];
  $t{$id}{disc_fill} = $t{$id}{disc_end} + $expAcc - $accLen;  # Discriminator based on filling from end of (short?) acc
  $t{$id}{disc_end} = $t{$id}{j_site} + $expAcc + 1 unless $accLen;  # Discriminator based on last acc bp (otherwise filling from end of tstem)
  if ($len > $t{$id}{disc_fill} and $t{$id}{seq} =~ s/CCA$// and $t{$id}{struct} =~ s/\.{3}$//) {$t{$id}{cca} = 'disc_fill'}
  elsif ($len == $t{$id}{disc_end} +3 and $t{$id}{seq} =~ s/CCA$// and $t{$id}{struct} =~ s/\.{3}$//) {$t{$id}{cca} = 'disc_end'}
  if ($t{$id}{cca}) {if ($t{$id}{dir} > 0) {$t{$id}{R} -= 3} else {$t{$id}{L} += 3}}

  my ($inttot, @introns) = (0);
  for my $i (0 .. $#{$t{$id}{intLs}}) {
   push @introns, "$t{$id}{intLs}[$i]..$t{$id}{intRs}[$i]";
   my $intlen = $t{$id}{intRs}[$i] - $t{$id}{intLs}[$i] +1;
   $t{$id}{a_site} += $intlen if $t{$id}{a_site} >= $t{$id}{intLs}[$i];
   $inttot += $intlen;
  }
  $t{$id}{introns} = join(',', @introns);
  for (qw/j_site disc_fill disc_end/) {$t{$id}{$_} += $inttot}  # These had been based on mature seq
  $final{$t{$id}{dna}}{$t{$id}{L}}{$t{$id}{R}} = $id;
 }
}

sub LoadCode {
 %anticodon2aa =      qw/TGA Ser GGA Ser CGA Ser AGA Ser GAA Phe AAA Phe TAA Leu CAA Leu GTA Tyr ATA Tyr GCA Cys ACA Cys CCA Trp 
 TAG Leu GAG Leu CAG Leu AAG Leu TGG Pro GGG Pro CGG Pro AGG Pro GTG His ATG His TTG Gln CTG Gln TCG Arg GCG Arg CCG Arg ACG Arg 
 TAT Ile GAT Ile AAT Ile CAT Met TGT Thr GGT Thr CGT Thr AGT Thr GTT Asn ATT Asn TTT Lys CTT Lys GCT Ser ACT Ser TCT Arg CCT Arg 
 TAC Val GAC Val CAC Val AAC Val TGC Ala GGC Ala CGC Ala AGC Ala GTC Asp ATC Asp TTC Glu CTC Glu TCC Gly GCC Gly CCC Gly ACC Gly/;
 %one_letter = (qw/Ala A Arg R Asn N Asp D Cys C Gln Q Glu E Gly G His H Ile I Ile2 I Leu L Lys K Met M Phe F Pro P SeC U Ser S
  Thr T Trp W Tyr Y Val V fMet B iMet B Sup * Undet X/);
}

sub RunCommand {
 my ($command, $checkfile) = @_;
 if ($checkfile and -e $checkfile) {print "Skipping command: $command\n"; return}
 print "Running command: $command\n";
 my $out = system($command);
 if ($out) {print "Command '$command' failed with error message $out\n"; exit}
 else {print "Command '$command' succeeded\n"}
}

