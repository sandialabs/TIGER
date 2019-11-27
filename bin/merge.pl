#! /usr/bin/perl
use strict; use warnings;

die "Called \"$0 @ARGV\" on " . localtime . "\nUsage: $0 file type outfile_prefix nickname\ne.g., perl $0 57779.island.gff island 57779 Eco837\n" unless @ARGV == 4;
my ($file, $type, $prefix, $nick) = @ARGV;
my (%comp, %isles, %serials);
my $tol = 10; # bp allowed to fall outside of group overlap regions; tol=0 gave 12843 compmerge.gff lines, tol=10 gave 12541
my ($compCutoff) = (3); 

open IN, "$file";
while (<IN>) {
 # NC_017765.1	comparator	island	9930995	9936042	3	+	.	brief=5HYP|ykoV;len=5047;contextsum=HYP>/>ykoV;prefCoords=9930995,9936041;bitsum=5923;gnm=NZ_CP013219.1:8562516-8565515;q1=94.105:11502-14268(9928230-9930996)>8559973-8562766;q2=93.739:250-1398>9936041-9937167;crossover=2;int=PhageIntegrase.9:9931135-9932322;mid=9931728;side=L612;OL=9930995-9930996;OR=9936041-9936042;OU=8562763-8562766;mobQ1=;mobQ2=;IS=;ISoverlap=;transposon=;ISidentical=;context=//hypothetical_protein/Lintergene/3prime/752,/ykoV_2//Rintergene/5prime/697;origOrient=+;q1identity=94.105;q2identity=93.739;isleLseq=ATGCCCCGAACGATCTGGTCCGGCGCGATCTCCTTCGGCCTGGTCACGGTgcTTATCTAGATCTGCATCACGGGTGTGAGGCACGTGAGGTCAGGTTTCATT;unintSeq=ATGCCCCGAACGATCTGGTCCGGCGCGATCTCCTTCGACCTGGTCACGGTgcCGATCAATGTGGTCGGCGCGACTGAGGACCACAGCATCCACTTCCACCAG;isleRseq=AGGAGAGCGCTCCTGACCTGCACATTCGCGGACGACATCTAGATAAGCACgcCGATCAATGTGATCGGCACCACCGAGGACCACAGCATCCACTTCCACCAG;mean=3255.66666666667;SD=1886.08948768492;deltaint=140;foreign=3.45631636123481;housekeep=13.1361382343199;hypoth=0.240129377648896;delta_GC=-0.07055;dinuc=0.06655;FPscore=6.71216010418801e-08;project=89409;division=Bacteria;phylum=Actinobacteria;order=Actinobacteria;class=Streptomycetales;family=Streptomycetaceae;genus=Streptomyces;species=Streptomyces hygroscopicus;org=Streptomyces hygroscopicus subsp. jinggangensis 5008;taxid=1133850;gencode=11;replicon=Chr;qlen=15000;ints=PhageIntegrase.9,PhageIntegrase.10;
 #next if /replicon=Plm/;
 chomp;
 next unless /\S/;
 my @f = split "\t";
 die "scores for $_\n" unless /(?:len|size)=([^;]+).*;side=(.).*delta_*int=(.*);foreign=(.*);housekeep=(.*);hypoth=(.*);delta_GC=(.*);dinuc=([^;]+);/;
 # 'side' as captured is single letter L/R only
 my ($scoComp, $scoIslr) = Score($1,$3,$4,$5,$6,$7,$8);
 my %l = (dna => $f[0], source => $f[1], L => $f[3], R => $f[4], supp => $f[5], orient => $f[6], line => $_, side => $2, f8 => $f[8], compScore => $scoComp, islrScore => $scoIslr, positive => 'true');
 for (split ';', $f[8]) {$l{$1} = $2 if /^(.+)=([^;]+)/}
 unless ($f[8] =~ /OL=(\d+)-(\d+);OR=(\d+)-(\d+);/) {print "no OL or OR for $_\n"; next}
 ($l{OLL}, $l{OLR}, $l{ORL}, $l{ORR}) = sort {$a <=> $b} ($1, $2, $3, $4);
 $l{positive} = 'false' if $scoComp >= $compCutoff;
 $l{line} .= "positive=$l{positive};"; 
 push @{$comp{$f[0]}}, \%l;
}
close IN;

open MRG, ">$prefix.$type.merge.gff";
open NOV, ">$prefix.$type.nonoverlap.gff";
for my $dna (keys %comp) {
 my %groups;
 @{$comp{$dna}} = sort {$$b{supp} <=> $$a{supp}} @{$comp{$dna}};
 my ($index) = (0); 
 for my $i (@{$comp{$dna}}) { # Make groups based on OL, OR
  $$i{index} = $index;
  my $group  = $index;
  if (defined $$i{group}) {
   $group = $$i{group};
   for (qw/OLL ORL/) {$groups{$group}{$_} = $$i{$_} if $$i{$_} < $groups{$group}{$_}} # Always expanding overlap
   for (qw/OLR ORR/) {$groups{$group}{$_} = $$i{$_} if $$i{$_} > $groups{$group}{$_}}
  } else {
   $$i{group} = $group;
   for (qw/OLL OLR ORL ORR/) {$groups{$group}{$_} = $$i{$_}}
  }
  $groups{$group}{members}{$index} ++;
  my $index2 = $index;
  $index ++;
  for my $j (@{$comp{$dna}}[$index..$#{$comp{$dna}}]) {
   $index2 ++;
   die "$#{$comp{$dna}}, $dna, $index, $index2, $$j{L}, $$j{R}, $groups{$group}{OLL}\n" unless $$j{L} and $groups{$group}{OLL};
   next if $$j{L} < $groups{$group}{OLL}-$tol or $$j{L} > $groups{$group}{OLR}+$tol;
   next if $$j{R} < $groups{$group}{ORL}-$tol or $$j{R} > $groups{$group}{ORR}+$tol;
   if ($$j{group} and $$j{group} ne $group and $groups{$$j{group}}) { # Two groups need merging
    my ($good, $bad) = sort {$a <=> $b} ($group, $$j{group});
    for (keys %{$groups{$bad}{members}}) {$groups{$good}{members}{$_} ++; $comp{$dna}[$_]{group} = $good}
    for (qw/OLL ORL/) {$groups{$good}{$_} = $groups{$bad}{$_} if $groups{$bad}{$_} < $groups{$good}{$_}} # Always expanding overlap
    for (qw/OLR ORR/) {$groups{$good}{$_} = $groups{$bad}{$_} if $groups{$bad}{$_} > $groups{$good}{$_}}
    delete $groups{$bad};
    $group = $good;
   }
   for (qw/OLL ORL/) {$groups{$group}{$_} = $$j{$_} if $$j{$_} < $groups{$group}{$_}} # Always expanding overlap
   for (qw/OLR ORR/) {$groups{$group}{$_} = $$j{$_} if $$j{$_} > $groups{$group}{$_}}
   $groups{$group}{members}{$index2} ++;
   $$j{group} = $group; 
  }
 }
 $index = 0;
 my @merge;
 for my $i (@{$comp{$dna}}) { # Transfer labels, add supports, for mergers
  $index ++;
  for (qw/islander brief contextsum IS/) {my $val = ''; $val = $$i{$_} if $$i{$_}; if ($_ eq 'IS') {$val =~ s/q[12],//g; $val =~ s/,score:[0-9]+//g}; $groups{$$i{group}}{$_}{$val} ++}
  next if $$i{reject};
  my (%support);
  for (keys %$i) {$support{"$$i{int}$$i{side}"}{$_} = $$i{$_}} # Combining supports for virtually same island
  for my $j (@{$comp{$dna}}[$index..$#{$comp{$dna}}]) {
   next if $$j{reject};
   next if $$j{group} ne $$i{group};
   for (qw/islander brief contextsum IS/) {my $val = ''; $val = $$j{$_} if $$j{$_}; if ($_ eq 'IS') {$val =~ s/q[12],//g; $val =~ s/,score:[0-9]+//g}; $groups{$$i{group}}{$_}{$val} ++}
   $$j{reject} ++;
   my $cat = "$$j{int}$$j{side}"; # Prevents overcounts from short island:qlen (found in both orientations) or from multiple mobility zymes (different ints)
   unless ($support{$cat}) {
    for (keys %$j) {$support{$cat}{$_} = $$j{$_}}
    $support{$cat}{supp} = 0; # Will add back in next line
   }
   $support{$cat}{supp} += $$j{supp};
  }
  push @merge, {};
  my $top = (sort {$support{$b}{supp} <=> $support{$a}{supp}} keys %support)[0];
  $support{$top}{line} =~ s/^(\S+\t\S+\t\S+\t\S+\t\S+\t)(\S+)/$1$support{$top}{supp}/;
  for (keys %{$support{$top}}) {next if /^reject$/; $merge[-1]{$_} = $support{$top}{$_}}
  for (qw/islander brief contextsum IS/) { # Did these vary within the group?
   my @vals = sort {$groups{$$i{group}}{$_}{$b} <=> $groups{$$i{group}}{$_}{$a}} keys %{$groups{$$i{group}}{$_}};
   for (@vals) {s/,/\//g}
   $merge[-1]{$_ . 's'} = join ',', @vals;
   $merge[-1]{$_ . 's'} =~ s/^,//; $merge[-1]{$_ . 's'} =~ s/,$//;
  }
 }
 @merge = sort {$$b{supp} <=> $$a{supp}} @merge;
 $index = 0;
 for my $i (@merge) {  # Better-supported first, no previous rejects remain
  $index ++;
  ($$i{overlap}, $$i{overlap2}) = (0,0) unless $$i{overlap};
  $$i{line} =~ s/FPscore=[^;]+/islrScore=$$i{islrScore};compScore=$$i{compScore}/; #$$i{FPscore} = $$i{score};
  my $diffs; for (qw/brief contextsum IS/) {$diffs .= $_ . "s=$$i{$_ . 's'};"}
  $$i{line} =~ s/IS=[^;]+/IS=/ if $$i{ISs} =~ /^,/ or $$i{ISs} =~ /,$/ or $$i{ISs} =~ /,,/; # When IS calls/noncalls among a merged group are mixed the noncall seems to usually be correct
  die "brief parse $_\n" unless $$i{brief} =~ s/^(\d+)\.//;
  my $len = $1;
  $$i{brief} =~ s/^\|//; $$i{brief} =~ s/\|$//;
  $$i{brief} =~ s/\|HYP//; $$i{brief} =~ s/HYP\|//;
  $$i{brief} =~ s/^(.)\|..+/$1/; $$i{brief} =~ s/^..+\|(.)$/$1/;
  $$i{name} = "$nick.$len.$$i{brief}";
  if ($serials{$$i{name}}) {$serials{$$i{name}} ++; $$i{name} .= ".$serials{$$i{name}}"} else {$serials{$$i{name}} ++}
  $$i{line} =~ s/brief=/ID=$$i{name};brief=/;
  $$i{line} =~ s/;OL=/;OLL=$groups{$$i{group}}{OLL};OLR=$groups{$$i{group}}{OLR};ORL=$groups{$$i{group}}{ORL};ORR=$groups{$$i{group}}{ORR};OL=/;
  if ($$i{overlap}) {print MRG $$i{line} . "status=overlapped;\n"; next} else {$$i{line} .= "status=nonoverlap;\n"; print MRG $$i{line}}
  for my $j (@merge[$index..$#merge]) {
   next if $$j{R} < $groups{$$i{group}}{OLR} or $$j{L} > $groups{$$i{group}}{ORL}; # Skip if no overlap
   my @overlapends = (sort {$a <=> $b} $$i{L}, $$i{R}, $$j{L}, $$j{R})[1,2];
   my $overlap = $overlapends[1] - $overlapends[0];
   ($$j{overlap}, $$j{overlap2}) = ($overlap, $$i{len});
  }
  print NOV $$i{line};
 }
}
close MRG;
close NOV;

sub Score {
 my ($len, $deltaint, $foreign, $housekeep, $hypoth, $deltaGC, $dinuc) = @_;
 my $scoComp = 2.71828182846 ** (-26.67594652 + 4.05867565 * log10($len) + 2.50735718 * log10($deltaint) + 0.08897193 * $foreign -
   0.57223169 * $housekeep - 8.13469540 * $hypoth + 0.96718617 * log10(abs($deltaGC)+0.01) + 14.69015423 * $dinuc);
 my $scoIslr = 2.71828182846 ** (-48.9141549 + 8.8517226 * log10($len) + 1.1004900 * log10($deltaint) + -0.6165349 * $foreign +
   -0.1554525 * $housekeep + -5.2098695 * $hypoth + -1.0174611 * log10(abs($deltaGC)+0.01) + -32.4246985 * $dinuc);
 return ($scoComp, $scoIslr);
}
sub log10 {return log(shift)/log(10)}

