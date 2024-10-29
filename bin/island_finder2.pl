#! /usr/bin/perl
use strict; use warnings;

# OPTIONS
use File::Spec;
use List::Util qw(shuffle);
my $path = File::Spec->rel2abs($0); $path =~ s/\/([^\/]+)$//; my $scriptname = $1;
my $dbpath = $path; $dbpath =~ s/[^\/]+$/db/;
my $verbose;
my $maxSize = 200000; # Maximum island size
my $minSize =   2000;
my $nickname = '';
my $criterion = 'score';
my $cross = 'cross';
Options(); # See bottom of script; help and other messages, Getopt
my $prefix = $ARGV[0];

# PROGRAM # Read infiles, Filter candidates, Treat tandems/overlaps, Write
my %aalookup = qw/Ala A Arg R Asn N Asp D Cys C Glu E Gln Q Gly G His H Ile I Ile2 J Leu L Lys K Met M Phe F Pro P Ser S Thr T Trp W Tyr Y Val V Pyl O SeC U tmRNA Z iMet B fMet B Sup X Undet ?/;
my (%domains, %proteins, %dnaSeq, %dnaRA, %integrases, %trnas, %reasons,
 $interval, @intervals, @isles, %tandem_gps, %intPfam, %hskpEnrich, %fornEnrich, %stats);
Read_stats();
Read_dom();
Read_integrase_pfam_ref();
Read_gff();
Read_fa();
Load_housekeep_enrich();
Check_intervals();
#TandemTile();
#ResolveOverlaps();
Print_islands();

# SUBROUTINES
sub Read_stats {
 for (@{ReadFile("$prefix.stats")}) {
  my $dna = $1 if /^(\S+)/;
  for my $cat (qw/dna len circ replicon cds pfam hypoth hskp forn/) {$stats{$dna}{$cat} = $1 if s/^(\S+)\t//}
  @{$dnaRA{$dna}} = ($stats{$dna}{dna}, split "\t"); # relative abundances of C and key dinucleotides
 }
}

sub ReadFile {
 my @ret;
 print "Reading file $_[0]\n" if $verbose;
 open IN, $_[0] or die "Can't open $_[0]\n";
 while (<IN>) {next if /^#/; push @ret, $_}
 close IN; chomp @ret; return \@ret;
}

sub Read_dom { # Read domains gff file, store the domain accession, coordinates, size
 for (@{ReadFile("$prefix.domains.gff")}) {
  my @f = split "\t";
  /ID=([^;]+);pfam=([^;]+)/;
  push @{$domains{$f[0]}}, {dna => $f[0], L => $f[3], R => $f[4], id => $1, fam => $2};
 }
}

sub Read_integrase_pfam_ref { # Allow int-associated pfam domains to overlap fragment
 for (@{ReadFile("$dbpath/integrase_domain_pfams.txt")}) {$intPfam{$1} = $1 if /^\S+\t(\S+)/}
}

sub Read_gff { # Read protein gff file, store the domain accession, coordinates, size
 for (@{ReadFile("$prefix.gff")}) {
  my @f = split "\t";
  my ($type, $L, $R) = @f[2,3,4];
  my $id; $id = $1 if $f[8] =~ /ID=(.*?);/;
  if ($type eq "CDS") {
   my $call = 'hypothetical';
   if (/pfam1=([^;]+)/) {$call = $1}
   if (/annot=(Y-Int[^;]+)/) {
    push @{$integrases{$f[0]}}, {dna => $f[0], L => $f[3], R => $f[4], id => $id, annot => $1};
   }
   push @{$proteins{$f[0]}}, {dna => $f[0], L => $L, R => $R, call => $call, id => $id};
  }
  #next unless /ID=([^;]+).*product=([^;]+).*j_site=([^;]+).*sequence=([^;]+).*cca=([^;]+)/;
  next unless $type =~ /^(tRNA|tmRNA|tRNA_gpIintron)$/;
  %{$trnas{$id}} = (dna => $f[0], L => $f[3], R => $f[4], orient => $f[6], dir => $f[6].'1', a_site => '', oneletter => '?') ;
  for (split ';', $f[8]) {/^([^=]+)=(.*)/; $trnas{$id}{$1} = $2}
  if ($type eq 'tmRNA') {$trnas{$id}{product} = $type . '_' . $trnas{$id}{form}; $trnas{$id}{aa} = $type}
  else {$trnas{$id}{product} = "${type}_$trnas{$id}{aa}($trnas{$id}{ac})"}
  $trnas{$id}{seq} .= $trnas{$id}{cca} unless $trnas{$id}{form} and $trnas{$id}{form} eq 'prm';
  #my ($aa) = ($2);
  #if ($aa =~ /tRNA-(.)(..)/) {$aa = uc($1) . lc($2)}
  #$trnas{$id}{aa} = 'tmRNA;
  #if (/undetermined|pseudogene|not called by both|low cove score|possibly problematic gene/i or $aa eq "Und") {$trnas{$id}{questionable} ='yes'} else {$trnas{$id}{questionable} =''}
  #$trnas{$id}{a_site} = $1 if /anticodon_center=(\d+)/;
  $trnas{$id}{len} = length $trnas{$id}{seq};
  $trnas{$id}{oneletter} = $aalookup{$trnas{$id}{aa}} if $aalookup{$trnas{$id}{aa}};
 }
}

sub Read_fa { # dna accession and size # use blastdbcmd for size
 my $accession;
 for (@{ReadFile("$prefix.fa")}) {if (/^>(\S+)/) {$accession = $1} else {$dnaSeq{$accession} .= $_} }
}

sub Load_housekeep_enrich {
 for (@{ReadFile("$dbpath/housekeep_enrich.txt")}) {$hskpEnrich{$1} = $2 if /(\S+)\t(\S+)/}
 for (@{ReadFile("$dbpath/foreign_enrich.txt")})   {$fornEnrich{$1} = $2 if /(\S+)\t(\S+)/}
}

sub Check_intervals {
 my $ct = 0;
 print "Checking intervals\n" if $verbose;
 my %hits;
 open (HITS, "$prefix.trnablast") or die "Can't open $prefix.trnabast\n";
 while (<HITS>) {
  chomp; $ct ++; my @f = split "\t";
  my ($trna, $dna, $qStart, $qEnd, $hitStart, $hitEnd) = @f[0,1,6..9];
  unless ($trnas{$trna}) {die "tRNA identifier \"$trna\" in blast results not found in ttm.fa file\n"}
  unless ($cross eq 'cross') { next unless $trnas{$trna}{dna} eq $dna}
  $hits{$trna}{$ct}{id} = "$trna.$ct";
  my $hit = $hits{$trna}{$ct};
  @{$hit}{qw\trna dna percent_id mismatches qStart qEnd hitStart hitEnd bit_score\} = @f[0,1,2,4,6..9,11];
  $$hit{hitDir} = abs($hitEnd - $hitStart)/($hitEnd - $hitStart);
  $$hit{hitOrient} = '+'; $$hit{hitOrient} = '-' if $$hit{hitDir} < 0;
  $$hit{dir} = $trnas{$trna}{dir};
  $$hit{qStart_genomic} = ($trnas{$trna}{L}, $trnas{$trna}{R})[0.5-0.5*$$hit{dir}] + $$hit{dir} * ($$hit{qStart} -1);
  $$hit{dupe} = 0;
  if ($qStart == 1) {$$hit{dupe} = -1}
  elsif ($qEnd >= $trnas{$trna}{len} - 7) {$$hit{dupe} = '+1'}
  if ($trnas{$trna}{a_site} and $qStart < $trnas{$trna}{a_site}) {$$hit{int_type} = 'a_site'}
  elsif ($qStart < $trnas{$trna}{j_site}) {$$hit{int_type} = 'j_site'}
  $$hit{int_site} = $trnas{$trna}{$$hit{int_type}} if $$hit{int_type};
 }
 for my $trna (keys %hits) {
  print "$trna\n";
  MAIN: for (keys %{$hits{$trna}}) {
   print "$_\n";
   my $hit = $hits{$trna}{$_};
   $interval = "$$hit{dna} trna:$trna $trnas{$trna}{orient}$$hit{qStart_genomic} frag:$$hit{hitOrient}$$hit{hitStart}";
   unless ($$hit{dupe}) {Reason('0. not from tRNA terminus'); next}
   $$hit{questionable} = '';
   for my $id (keys %trnas) {
    if ($$hit{hitStart} == $$hit{qStart_genomic}) {Reason('1. hit same or already-known tRNA'); next MAIN}
    if ($$hit{hitStart} >= $trnas{$id}{L} && $$hit{hitStart} <= $trnas{$id}{R} or $$hit{hitEnd} >= $trnas{$id}{L} && $$hit{hitStart} <= $trnas{$id}{R}) {
     if ($trnas{$id}{questionable}) {$$hit{questionable} .= "$id,"; next}
     Reason('1. hit same or already-known tRNA'); next MAIN;
    }
   }
   if ($$hit{dna} eq $trnas{$trna}{dna} and $$hit{hitDir} != $$hit{dir}) {Reason('2. opposite tRNA and hit orientations on same contig'); next}
   unless ($$hit{int_type}) {
    Reason('3. covers neither A nor J site'); next;
   }
   $interval .= ' dupe:' . (4-$$hit{dupe}) . 'prime';
   my ($coord, $compos);
   my $start_to_site = $$hit{int_site} - $$hit{qStart};
   my ($proximal, $distal) = ($$hit{qStart_genomic}, $$hit{hitStart});
   $proximal += $$hit{dir} * $start_to_site; $distal += $$hit{hitDir} * $start_to_site;
   my ($isleStart, $isleEnd) = ($proximal, $distal);
   ($isleStart, $isleEnd) = ($distal, $proximal) if $$hit{dupe} == -1;
   my $endHi = abs($isleEnd - $isleStart)/($isleEnd - $isleStart);
   my $origin;
   my @island = sort {$a <=> $b} ($isleStart, $isleEnd);
   @island = reverse @island if $origin == 1;
   if ($$hit{dna} eq $trnas{$trna}{dna}) {
    for (-1, 1) {
     $origin = $_ if $_ != $$hit{dir} * $endHi;
     last unless $stats{$$hit{dna}}{circ} eq 'Cir' or $cross ne 'simple';
     last if $origin;
    }
    unless ($origin) {Reason('4. tRNA duplication on wrong side of linear DNA'); next}
    if ($origin == 1) {
     if ($cross ne 'simple' and $stats{$$hit{dna}}{circ} ne 'Cir') {
      $coord = "$$hit{dna}/$isleEnd-$stats{$$hit{dna}}{len}+$$hit{dna}/1-$isleStart";
      $compos = "circleJxn";
     } elsif ($stats{$$hit{dna}}{circ} eq 'Cir') {
      $coord = "$$hit{dna}/$isleEnd-$isleStart";
      $compos = "simple";
     }
    } else {
     $coord = "$$hit{dna}/$isleStart-$isleEnd";
     $compos = "simple";
    }
   } else {
    $origin = -1;
    $compos = "cross";
    my $trnaEnd = (1, $stats{$trnas{$trna}{dna}}{len})[0.5-0.5*$$hit{dir}*$$hit{dupe}];
    my $hitEnd = ($stats{$$hit{dna}}{len}, 1)[0.5-0.5*$$hit{hitDir}*$$hit{dupe}];
    my ($trnaStart, $hitStart) = ($isleStart, $isleEnd);
    ($trnaStart, $hitStart) = ($isleEnd, $isleStart) if $$hit{dupe} == -1;
    my $trnaCoord = "$trnas{$trna}{dna}/$trnaStart-$trnaEnd";
    my $hitCoord = "$$hit{dna}/$hitStart-$hitEnd";
    my $coord = ($$hit{dupe} == 1) ? "$trnaCoord+$hitCoord" : "$hitCoord+$trnaCoord";
    @island = sort {$a <=> $b} ($trnaStart, $trnaEnd);
   }
   my $isleSize;
   unless ($compos eq "cross") {
    $isleSize = abs($isleStart - $isleEnd);
    $isleSize = $stats{$$hit{dna}}{len} - $isleSize if $origin == 1;
   } else {
    die unless $coord =~ /[^\/]+\/(\d+)-(\d+)\+[^\/]+\/(\d+)-(\d+)/;
    $isleSize = abs($1 - $2) + abs($3 - $4);
   }
   if ($isleSize <= $minSize) {Reason('5. too short'); next}
   die unless $coord =~ /^([^\/]+)\/(\d+)-(\d+)/;
   my @ts;
   push @ts, {dna => $1, L => $2, R => $3};
   if ($compos ne "simple") {
    die unless $coord =~ /\+([^\/]+)\/(\d+)-(\d+)/;
    push @ts, {dna => $1, L => $2, R => $3};
   }
   my $ints = IntsWithin(\@ts);
   unless (keys %$ints) {Reason('6. contains no integrase gene'); next}
   my $overlap = 0;
   DOMAIN: for my $dom (@{$domains{$$hit{dna}}}) {
    next if $$ints{$$dom{id}} or $intPfam{$$dom{fam}};
    for ($$hit{hitStart}, $$hit{hitEnd}) {
     next unless $_ >= $$dom{L} && $_ <= $$dom{R};
     $overlap ++;
     last DOMAIN;
    }
   }
   if ($overlap) {Reason('7. non-integrase protein gene core overlap'); next}
   my ($Lseq, $Rseq) = getSeq($coord, $origin, $$hit{dupe});
   my $id = scalar(@isles);
   my $group = (4-$$hit{dupe})."prime.$trna";
   push @isles, {
    id => $id,
    dna => $$hit{dna},
    trna => $trna,
    group => $group,
    dir => $$hit{dir},
    strand => substr ($$hit{dir}, 0, 1),
    L => $island[0],
    R => $island[1],
    coord => $coord,
    compos => $compos,
    proximal => $proximal,
    distal => $distal,
    size => $isleSize,
    isleLseq => $Lseq,
    isleRseq => $Rseq,
    int_site => $$hit{int_site},
    int_site_type => $$hit{int_type},
    trna_dupe => $$hit{dupe},
    tRNA_L => $trnas{$trna}{L},
    tRNA_R => $trnas{$trna}{R},
    tRNA_len => $trnas{$trna}{R}-$trnas{$trna}{L}+1,
    tRNA_aaa => $trnas{$trna}{aa},
    tRNA_aa => $trnas{$trna}{oneletter},
    A_site => $trnas{$trna}{a_site},
    J_site => $trnas{$trna}{j_site},
    mismatches => $$hit{mismatches},
    origin => $origin,
    qStart => $$hit{qStart},
    qEnd => $$hit{qEnd},
    hitStart => $$hit{hitStart},
    hitEnd => $$hit{hitEnd},
    percent_id => $$hit{percent_id},
    bit_score => $$hit{bit_score},
    questionable => $trnas{$trna}{questionable},
    overall => 0,
   };
   #push @{$tandem_gps{$group}{orig}}, $id if $$hit{hitDir} == $$hit{dir};
   if ($isleSize >= $maxSize) {Reason('8. too long (either way around circle)'); next}
   Reason('9. Accepted after initial interval check');
   @{$isles[-1]}{qw/delta_int foreign housekeep hypoth delta_GC dinuc overall/} = Scores ($coord, $origin);
  }
 }
 if ($verbose) {
  print "$ct candidate islands (blast hits) initially\n";
  my $sum = 0;
  for (sort keys %reasons) {$sum += $reasons{$_}; print "$reasons{$_} got $_: runningTotal=$sum\n"}
  print scalar(@isles), " candidate islands (including too-longs for now)\n";
 }
 %reasons = (); 
}

sub Reason { # Handle reasons for rejection from Check_intervals
 $reasons{$_[0]} ++;
 print "$interval $_[0]\n" if $verbose;
}

sub IntsWithin {
 my ($e) = @_;
 my %ints;
 for my $i (0,1) {
  next unless $$e[$i];
  #warn "$$e[$i]{dna}";
  my ($dna, $L, $R) = @{$$e[$i]}{qw/dna L R/};
  for my $int (@{$integrases{$dna}}) { # int must be over half within island
   my ($within, $midpt) = (0, $$int{mid});
   $within = 1 if $midpt > $$e[$i]{L} && $midpt < $$e[$i]{R};
   %{$ints{"$$int{annot}:$$int{L}-$$int{R}"}} = (dna => $dna, L => $$int{L}, R => $$int{R}) if $within;
  }
 }
 #warn keys %ints;
 return \%ints;
}

sub TandemTile { # tile spans within tandem group, until too-long int-less break
 for my $gp (keys %tandem_gps) {
  my ($break, $last, $origin, @segments); my $segment = 0;
  for my $i (sort {$isles[$a]{size} <=> $isles[$b]{size}} @{$tandem_gps{$gp}{orig}}) { # isle ids, working from proximal toward distal
   $segment ++;
   if ($break) {$isles[$i]{reject} = '11. Broken chain'; $tandem_gps{$gp}{sum} .= ",$segment"; next}
   if (not $origin) {$last = $isles[$i]{proximal}; $origin = $isles[$i]{origin}} # Initialize
   elsif ($isles[$i]{origin} == 1) { # Only one segment can cross origin
    if ($origin == -1) {$origin = 1} # First encounter of origin-spanner
    else {$origin = -2; $isles[$i]{origin} = -1} # Have passed origin-spanning, can't return to 1; new segment will inherit non-span from $isle
   }
   my $prox = $last; $last = $isles[$i]{distal};
   my $size = abs($last - $prox); $size = $stats{$isles[$i]{dna}}{len} - $size if $origin == 1;
   if ($size >= $maxSize) {$isles[$i]{reject} = '10. Too long'; $break ++; $tandem_gps{$gp}{sum} .= ",B$segment";}
   my ($intTest, $finalTest) = ('', 0);   # Final segment
   if ($break) {$intTest = 'i' if @segments and scalar keys %{$segments[-1]{ints}}; $finalTest = 1}
   else {
    my @ts; push @ts, $isles[$i];
    my $ints = IntsWithin(@ts);
    $intTest = 'i' if scalar keys %{$ints}; # New segment has an int
    @{$isles[$i]}{qw/size segment ints proximal/} = ($size, $segment.$intTest, $ints, $prox);
    push @segments, {}; for (keys %{$isles[$i]}) {$segments[-1]{$_} = $isles[$i]{$_}}
    $segments[-1]{segment} .= 'o' if $origin == 1;
    push @{$tandem_gps{$gp}{frag}}, "$isles[$i]{hitStart}-$isles[$i]{hitEnd}";
    for (sort keys %{$ints}) {push @{$tandem_gps{$gp}{ints}}, ${$isles[$i]{ints}}{$_}{annot}}
    $finalTest = 1 if $segment == scalar @{$tandem_gps{$gp}{orig}};
   }
   if ($intTest or $finalTest) { # Skip if internal segment without int; intTest only: split only, finalTest only: proximal only, Both: split if possible else whole unit
    my (@best, @best1); 
    SEG: for my $s (0 .. $#segments) {
     my @pair;
     last if $s == $#segments and not $finalTest; # Prepare full concat only if final segment
     my $score = 0; my $min = ($s + 1, $#segments)[$s + 1 > $#segments];
     my @pairfinal = ({proximal => $segments[0]{proximal}, distal=> $segments[$s]{distal}, origin => $segments[0]{origin}},
                   {proximal => $segments[$min]{proximal} , distal=> $last,                 origin => $segments[-1]{origin}});
     my $oripair = -1;
     my $oridna = '';
     for my $s2 (0 .. $#segments) {
      if ($s2 <= $s) {
       push @{$pair[0]}, {}; for (keys %{$segments[$s2]}) {$pair[0][-1]{$_} = $segments[$s2]{$_};}
       if ($segments[$s2]{origin} == 1) { $oripair = 0; $oridna = $segments[$s2]{dna}; $pairfinal[0]{origin} = 1;}
      } else {
       push @{$pair[1]}, {}; for (keys %{$segments[$s2]}) {$pair[1][-1]{$_} = $segments[$s2]{$_};}
       if ($segments[$s2]{origin} == 1) { $oripair = 1; $oridna = $segments[$s2]{dna}; $pairfinal[1]{origin} = 1;}
      }
     }
     unless ($oripair == -1) {
      for my $k (0 .. $#{$pair[$oripair]}) {
       $pair[$oripair][$k]{origin} = 1 if $pair[$oripair][$k]{dna} eq $oridna;
      }
     }
     if ($segments[0]{dna} ne $segments[$s]{dna}) {
      $pairfinal[0]{cross} = 'cross';
      $pairfinal[0]{coord} = "$segments[0]{dna}/$segments[0]{proximal}-$segments[0]{distal}+$segments[$s]{dna}/$segments[$s]{proximal}-$segments[$s]{distal}";
     } else {
      $pairfinal[0]{cross} = 'simple';
      $pairfinal[0]{coord} = "$segments[0]{dna}/$pairfinal[0]{proximal}-$pairfinal[0]{distal}";
     }
     if ($segments[$min]{dna} ne $segments[-1]{dna}) {
      $pairfinal[1]{cross} = 'cross';
      $pairfinal[1]{coord} = "$segments[$min]{dna}/$segments[$min]{proximal}-$segments[$min]{distal}+$segments[-1]{dna}/$segments[-1]{proximal}-$segments[$s]{distal}";
     } else {
      $pairfinal[1]{cross} = 'simple';
      $pairfinal[1]{coord} = "$segments[$min]{dna}/$pairfinal[1]{proximal}-$pairfinal[1]{distal}";
     }
     for my $j (0, 1) {
      my ($L, $R) = ($pairfinal[$j]{proximal}, $pairfinal[$j]{distal}); ($L, $R) = sort {$a <=> $b} $L, $R;
      my $size = 0;
      for my $l (0 .. $#{$pair[$j]}) {
       if ($pair[$j][$l]{origin} == 1) {
        $size += $stats{$pair[$j][$l]{dna}}{len} - abs($pair[$j][$l]{proximal} - $pair[$j][$l]{distal}); ($L, $R) = ($R, $L);
       } else {
        $size += abs($pair[$j][$l]{proximal} - $pair[$j][$l]{distal});
       }
      }
      my $ints2 = IntsWithin(\@{$pair[$j]});
      next SEG if $size >= $maxSize or $size <= $minSize or not $ints2; # allow distal to be short, may elongate?
      #warn ref($pairfinal[$j]);
      @{$pairfinal[$j]}{qw/L R size ints delta_int foreign housekeep hypoth delta_GC dinuc overall/} =
       ($L, $R, $size, $ints2, Scores($ints2, $size, \@{$pair[$j]}));
      $score += $pairfinal[$j]{overall};  # !!!! Scores not additive !!!!
      @best1 = ($score, $s, \@pairfinal) if $j == 0 and $finalTest and (not @best1 or $best1[0] > $score); # Final may exclude terminal segmen
      last if $j==0 and ($s == $#segments or not $intTest); # Last distal is empty; only testing proximals if finalTest onl
	     }
     last if $s == $#segments and $finalTest and $intTest and @best; # Consider full concat when final has int only if no split available
     @best = ($score, $s, \@pairfinal) if not @best or $best[0] > $score;
    }
    @best = @best1 unless @best; # Best1 is among proximal-only's when finalTest, this may occur when intTest also
    for my $j (0,1) {
     next unless @best;
     my $id = scalar @isles;
     my $bestidx = $best[1];
     warn "$best[1]";
     push @{$tandem_gps{$gp}{tiled}}, $id;
     $segments[$bestidx]{id} = $id;
     for (keys %{$segments[$bestidx]}) {$isles[$id]{$_} = $segments[$bestidx]{$_}; warn "$_: $segments[$bestidx]{$_}\n"}
     for (keys %{${$best[2]}[$j]}) {$isles[$id]{$_} = ${$best[2]}[$j]{$_}; warn "$_: ${$best[2]}[$j]{$_}\n"}
     my @parts; for my $j (0 .. $best[1]) {push @parts, $segments[0]{segment}; shift @segments}
     $isles[$id]{segments} = '[' . join(',', @parts) . ']';
     $tandem_gps{$gp}{sum} .= $isles[$id]{segments};
     $isles[$id]{segment} = scalar @{$tandem_gps{$gp}{tiled}};
     $best[1] = $#segments; last if $best[1] == -1; # Final island (if called) includes last segment
    }
   }
  }
  for (@{$tandem_gps{$gp}{orig}}) {$isles[$_]{reject} = "Retiled" unless $isles[$_]{reject};}
 }
 #exit;
}

sub Scores {
 my ($ints, $len, $e) = @_;
 #($L, $R) = sort {$a <=> $b} $L, $R; ($L, $R) = ($R, $L) if $origin == 1;
 #print "scoring $L, $R, $origin\n";
 my $delta_int = DeltaInt($ints, $stats{all}{len}, \$e);
 my ($hypoth, $forn, $hskp) = Foreign($ints, \$e);
 my ($delta_GC, $dinuc) = Bias(\$e);
 #my $overall = exp(-31.8612 + 5.4279*log10($len) + 0.9279*log10($delta_int) + -0.4803*$forn + -0.0586*$hskp + -3.2921*$hypoth + -2.1053*log10(abs($delta_GC)+0.01) + 21.3018*$dinuc); #FP score as of 6/14/17; cutoff= 1.123398
 ##sixpercent##
 #noscore#my $overall = exp(-46.9060924 + 8.8548917*log10($len) + 1.2061817*log10($delta_int) + -0.5071771*$forn + -0.2149491*$hskp + -4.7611535*$hypoth + -0.6902204*log10(abs($delta_GC)+0.01) + -72.1607107*$dinuc); #FP score as of 6/14/17 iteration 1; cutoff= 0.6329089
 #iteration1# my $overall = exp(-138.6464535 + 17.4400839*log10($len) + 13.8959223*log10($delta_int) + -2.6253524*$forn + -0.3931175*$hskp + -17.3815851*$hypoth + 7.0123013*log10(abs($delta_GC)+0.01) + 133.9646467*$dinuc); #FP score as of 6/14/17 iteration 1; cutoff= 0.6329089
 #iteration2#my $overall = exp(-59.0850486 + 7.2883969*log10($len) + 5.6108163*log10($delta_int) + 0.4253299*$forn + -0.3019590*$hskp + -12.7197672*$hypoth + 1.6160732*log10(abs($delta_GC)+0.01) + 27.4692643*$dinuc); #FP score as of 6/14/17 iteration 2; cutoff= 0.7948245
 ##sevenpercent##
 #noscore#my $overall = exp(-41.6759067 + 7.7320839*log10($len) + 1.1147119*log10($delta_int) + -0.6388200*$forn + -0.2571841*$hskp + -5.7317222*$hypoth + -1.4778523*log10(abs($delta_GC)+0.01) + -46.2233896*$dinuc); #FP score as of 6/14/17 noscore; cutoff= 1.374656
 #my $overall = exp(-53.19720867 + 9.79511470*log10($len) + 1.11791950*log10($delta_int) + -0.73160503*$forn + -0.09266934*$hskp + -6.14852295*$hypoth +  -0.59046929*log10(abs($delta_GC)+0.01) + -36.04755375 *$dinuc); #FP score as of 6/14/17 noscore; cutoff= 1.230556
 ##final island list##
 #sevenandahalf_percent#
 #original#my $overall = exp(-43.7650175 + 8.0473167*log10($len) + 1.2673848*log10($delta_int) + -0.7526923*$forn + -0.2293208*$hskp + -3.7405608*$hypoth + -0.6821956*log10(abs($delta_GC)+0.01) + -49.1119053*$dinuc); #cutoff= 0.8151459; 269 FPs; lost double positives = 4
 my $overall = exp(-48.9464877 + 8.8456270*log10($len) + 1.1003517*log10($delta_int) + -0.6153832*$forn + -0.1563892*$hskp + -5.2095476*$hypoth + -1.0500543*log10(abs($delta_GC)+0.01) + -31.9852781*$dinuc); #cutoff= 1.361073; 233 FPs; lost double positives = 2
 #BEST9/6/18 my $overall = exp(30.6077118+ -8.2975315*log10($len) + -0.3310587*log10($delta_int) + -3.2087991*$forn + -1.7145020 *$hskp + -22.2829392*$hypoth + -16.0497881*log10(abs($delta_GC)+0.01) + -149.2433560*$dinuc); #cutoff=0.51
 return ($delta_int, $forn, $hskp, $hypoth, $delta_GC, $dinuc, $overall);
}

sub log10 {return log(shift)/log(10);}

sub DeltaInt { # Shortest distance between any internal integrase gene end and an island end
 my ($ints, $ret, $c) = @_;
 if (ref($$c) eq 'ARRAY') {
  for my $i (0 .. scalar @{$$c} - 1) {
   my ($dna, $L, $R, $origin) = @{$$$c[$i]}{qw/dna L R origin/};
   $R += $stats{$dna}{len} if $origin == 1; # virtual right end for origin-spanning island
   for (keys %$ints) {
    next unless $dna eq $$ints{$_}{dna};
    my ($intL, $intR) = ($$ints{$_}{L}, $$ints{$_}{R});
    if ($origin == 1 and $intL < $R) { for ($intL, $intR) {$_ += $stats{$dna}{len}} }
    $ret = $intL-$L if $intL-$L < $ret;
    $ret = $R-$intR if $R-$intR < $ret;
   }
  }
 } else {
  my ($dna, $L, $R, $origin) = @{$$c}{qw/dna L R origin/};
  $R += $stats{len} if $origin == 1; # virtual right end for origin-spanning island
  for (keys %$ints) {
   my ($intL, $intR) = ($$ints{$_}{L}, $$ints{$_}{R});
   if ($origin == 1 and $intL < $R) { for ($intL, $intR) {$_ += $stats{len}} }
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
 if (ref($$c) eq 'ARRAY') {
  for my $i (0 .. scalar @{$$c} - 1) {
   my ($dna, $L, $R, $origin) = @{$$$c[$i]}{qw/dna L R origin/};
   for my $prot (@{$proteins{$dna}}) {
    unless ($$ints{$$prot{id}}) {
     if ($origin == 1) {
      next if $$prot{L} <= $L and $$prot{R} >= $R;
     } else {
      next unless $$prot{L} >= $L and $$prot{R} <= $R;
      last if $$prot{L} > $R
     }
    }
    $cds ++;
    if ($$prot{call} eq 'hypothetical') {$hypoth ++}
    else {
     $pfam ++;
     $forn += $fornEnrich{$$prot{call}} if $fornEnrich{$$prot{call}};
     $hskp += $hskpEnrich{$$prot{call}} if $hskpEnrich{$$prot{call}};
    }
   }
  }
 } else {
  my ($dna, $L, $R, $origin) = @{$$c}{qw/dna L R origin/};
  for my $prot (@{$proteins{$dna}}) {
  	#die "$L $R $origin $x $$prot{L} $$prot{R}\n" unless $$prot{id};# $$ints{$$prot{id}}\n";
   unless ($$ints{$$prot{id}}) {
    if ($origin == 1) {next if $$prot{L} <= $L and $$prot{R} >= $R}
    else {next unless $$prot{L} >= $L and $$prot{R} <= $R; last if $$prot{L} > $R}
   }
   $cds ++;
   if ($$prot{call} eq 'hypothetical') {$hypoth ++}
   else {
    $pfam ++;
    $forn += $fornEnrich{$$prot{call}} if $fornEnrich{$$prot{call}};
    $hskp += $hskpEnrich{$$prot{call}} if $hskpEnrich{$$prot{call}};
   }
  }
 }
 $hypoth /= $cds;  $hypoth -= $stats{all}{hypoth};
 $forn   /= $pfam; $forn   -= $stats{all}{forn};
 $hskp   /= $pfam; $hskp   -= $stats{all}{hskp};
 return $hypoth, $forn, $hskp;
}

sub Bias {
 my ($c) = @_;
 system "echo \"\" > test.fa";
 if (ref($$c) eq 'ARRAY') {
  for my $i (0 .. scalar @{$$c} - 1) {
   my ($dna, $L, $R, $origin) = @{$$$c[$i]}{qw/dna L R origin/};
   if ($origin == -1) {system "perl $path/collectSeq.pl -i $prefix.fa -e $dna -L $L -R $R >> test.fa"}
   else {
    system "perl $path/collectSeq.pl -i $prefix.fa -e $dna -L $L -R $stats{$dna}{len} >> test.fa";
    system "perl $path/collectSeq.pl -i $prefix.fa -e $dna -L 1 -R $R >> test.fa";
   }
  }
  my $out = `perl $path/relAbun.pl test.fa`; chomp $out; $out =~ s/\n.*//s; my @testRA = split "\t", $out;
  my ($delta_GC, $di) = (($testRA[1]-$dnaRA{all}[1])/2, 0);
  for (2,3,4,6,7,9) {$di += abs($testRA[$_]-$dnaRA{all}[$_])}
  $di *= 2; # Double the 6 asymmetrical dinucs above to account for their complements, but don't double the 4 symmetrical dinucs below
  for (5,8,10,11) {$di += abs($testRA[$_]-$dnaRA{all}[$_])}
  return ($delta_GC, $di/16);
 } else {
  my ($dna, $L, $R, $origin) = @{$$c}{qw/dna L R origin/};
  if ($origin == -1) {system "perl $path/collectSeq.pl -i $prefix.fa -e $dna -L $L -R $R >> test.fa"}
  else {
   system "perl $path/collectSeq.pl -i $prefix.fa -e $dna -L $L -R $stats{$dna}{len} >> test.fa";
   system "perl $path/collectSeq.pl -i $prefix.fa -e $dna -L 1 -R $R >> test.fa";
  }
  my $out = `perl $path/relAbun.pl test.fa`; chomp $out; $out =~ s/\n.*//s; my @testRA = split "\t", $out;
  my ($delta_GC, $di) = (($testRA[1]-$dnaRA{all}[1])/2, 0);
  for (2,3,4,6,7,9) {$di += abs($testRA[$_]-$dnaRA{all}[$_])}
  $di *= 2; # Double the 6 asymmetrical dinucs above to account for their complements, but don't double the 4 symmetrical dinucs below
  for (5,8,10,11) {$di += abs($testRA[$_]-$dnaRA{all}[$_])}
  return ($delta_GC, $di/16);
 }
}

sub ResolveOverlaps {
 for my $isle (@isles) { warn "$$isle{overall}\n"; }
 if ($criterion eq 'score') {@isles = sort { $$a{overall} <=> $$b{overall} } @isles} # lower FP scores = better islands (reverse if using bitscore)
 elsif ($criterion eq 'deltaGC') {@isles = sort {$$b{delta_GC} <=> $$a{delta_GC} } @isles}
 elsif ($criterion eq 'random') {@isles = shuffle @isles}
 for my $i (0 .. $#isles - 1) { # Find winner for overlap from different tandem group
  next if $isles[$i]{reject};
  for my $j ($i+1 .. $#isles) {
   next if $isles[$j]{reject};
   next if $isles[$i]{group} eq $isles[$j]{group}; # Only different tRNA group can reject at this stage
   if     ($isles[$i]{origin} == -1 and $isles[$j]{origin} == -1) { # Neither crosses origin
    next if $isles[$i]{L} > $isles[$j]{R} or  $isles[$j]{L} > $isles[$i]{R}  # Skip if no overlap
   } elsif ($isles[$i]{origin} == -1 xor $isles[$j]{origin} == -1) { # Only one crosses origin
    next if $isles[$i]{L} > $isles[$j]{R} and $isles[$j]{L} > $isles[$i]{R} # Skip if no overlap
   } # Remaining pair overlaps (automatically if both cross origin); determine who's better: 1. questionable loses, else 2. by overall score
   if ($isles[$j]{questionable} and not $isles[$i]{questionable}) {
    push @{$isles[$j]{better}}, $isles[$i]{id};
    print "$isles[$i]{group}:$isles[$i]{proximal}-$isles[$i]{distal} beat questionable $isles[$j]{group}:$isles[$j]{proximal}-$isles[$j]{distal}\n";
   } elsif ($isles[$i]{questionable} and not $isles[$j]{questionable}) {
    push @{$isles[$i]{better}}, $isles[$j]{id};
    print "$isles[$j]{group}:$isles[$j]{proximal}-$isles[$j]{distal} beat questionable $isles[$i]{group}:$isles[$i]{proximal}-$isles[$i]{distal}\n";
   } elsif ($isles[$i]{overall} > $isles[$j]{overall}) {
    push @{$isles[$j]{better}}, $isles[$i]{id};
    print "$isles[$i]{group}:$isles[$i]{proximal}-$isles[$i]{distal} beat $isles[$j]{group}:$isles[$j]{proximal}-$isles[$j]{distal}\n";
   } else {
    push @{$isles[$i]{better}}, $isles[$j]{id};
    print "$isles[$j]{group}:$isles[$j]{proximal}-$isles[$j]{distal} beat $isles[$i]{group}:$isles[$i]{proximal}-$isles[$i]{distal}\n";
   }
  }
 }
 # Now we have hierarchies of overlap winners: isle A need not be rejected by B, if B was rejected by C; hope there are no circular relationships (A < B < C < A)?
 @isles = sort {$$a{id} <=> $$b{id}} @isles;
 for my $i (0 .. $#isles) {$isles[$i]{on} = 1 unless $isles[$i]{better} or $isles[$i]{reject} }
 for my $i (0 .. $#isles) {
  $isles[$i]{tiling_order} = '';
  next if defined $isles[$i]{on} or $isles[$i]{reject};
  IntegrateBettersSignals($i, @{$isles[$i]{better}});
 }
 for my $gp (keys %tandem_gps) { # Prepare all tilings; most numerous, then highest summed score wins
  my ($break, @tilings, @win);
  #print "Group $gp: ", join(',', @{$tandem_gps{$gp}}), "\n";
  for my $i (@{$tandem_gps{$gp}{orig}}) { # Seed tilings with tRNA-touching islands only
   next if $isles[$i]{reject} or not $isles[$i]{trna}; # tRNA-touching islands
   push @tilings, {isles => {$i => 1}, distal => $isles[$i]{distal}, ct => 1, score => $isles[$i]{overall} };
   @win = (1, $isles[$i]{overall}, $#tilings) if @win == 0 or $win[1] < $isles[$i]{overall};
   print "Tiling: gp $gp; isles=$i,; ct=1; distal=$isles[$i]{distal}; score=$isles[$i]{overall}; winner=$win[1]\n" if $verbose; 
  }
  for my $i (@{$tandem_gps{$gp}{orig}}) { # Append segment when its proximal end matches tiling's distal end
   next if $isles[$i]{reject} or $isles[$i]{trna}; # Segments
   for my $oldtiling (@tilings) {
    #print "old{distal}=$$oldtiling{distal} new(isle $i){proximal}=$isles[$i]{proximal}\n";
    next unless $$oldtiling{distal} == $isles[$i]{proximal};
    push @tilings, {distal => $isles[$i]{distal}, ct => $$oldtiling{ct}+1, score => $$oldtiling{score}+$isles[$i]{overall} };
    my ($serial, $out);
    for ((sort {$$oldtiling{isles}{$a} <=> $$oldtiling{isles}{$b}} keys(%{$$oldtiling{isles}})), $i) {
     $serial ++;
     $tilings[-1]{isles}{$_} = $serial;
     $out .= "$_,";
    }
    @win = ($tilings[-1]{ct}, $tilings[-1]{score}, $#tilings) if $win[0] < $tilings[-1]{ct} or $win[1] < $tilings[-1]{score};
    print "Tiling: gp $gp; isles=$out; ct=$tilings[-1]{ct}; distal=$tilings[-1]{distal}; score=$tilings[-1]{score}; winner=$win[1]\n" if $verbose; 
   }
  }
  for my $i (@{$tandem_gps{$gp}{orig}}) {
   next if $isles[$i]{reject};
   unless ($tilings[$win[2]]{isles}{$i}) {$isles[$i]{reject} = '13. Not part of best tiling'; next};
   $isles[$i]{tiling_order} = $tilings[$win[2]]{isles}{$i};
  }
 }
}

sub IntegrateBettersSignals { # integrate on/off signals from better-scoring overlapping islands, hierarchically
 my ($isle, $signal, $overlapper) = (shift, 1, '');
 for my $better (@_) {
  IntegrateBettersSignals($better, @{$isles[$better]{better}}) unless defined $isles[$better]{on};
  $signal *= 1-$isles[$better]{on}; # Switch toggle (ON better sends OFF signal and vice versa); OFF signal (multiplication by zero) dominates
  $overlapper = $better if $signal == 0;
 }
 $isles[$isle]{on} = $signal;
 return if $signal; # Done if island on, else reject island due to overlap by better-scoring island
 my $trna = 'different'; $trna = 'same' if $isles[$overlapper]{group} eq $isles[$isle]{group};
 $isles[$isle]{reject} = "12. Better-scoring ON $isles[$overlapper]{strand} $isles[$overlapper]{L}-$isles[$overlapper]{R} island " .
  "in $trna $isles[$overlapper]{group} at $isles[$overlapper]{tRNA_R} overlaps";
 print "Better-scoring ON island $overlapper overlaps $isles[$isle]{group}:$isles[$isle]{proximal}-$isles[$isle]{distal}\n" if $verbose;
}

sub Print_islands{ #This outputs the final results of the program (overlaps still unresolved)
 print "Outputting islands\n" if $verbose;
 open GFF, ">islander.island.gff" || die "cannot write to islander.island.gff\n";
 #open OVP, ">islander.islandOverlapped.gff" || die "cannot write to islander.islandOverlapped.gff\n"; # Potential islands/tandems rejected due to overlaps
 my $serial = 0;
 for (sort {$$a{L} <=> $$b{L} || $$a{R} <=> $$b{R}} @isles) {
  my ($overlapped, $out);
  if ($$_{reject}) {$interval = "$$_{group}:$$_{L}..$$_{R}"; Reason($$_{reject}); $overlapped ++; next unless $$_{reject} =~ /^1[23]\./}
  $$_{L} ++ if $$_{L} == $$_{distal};
  $$_{R} -- if $$_{R} == $$_{distal};
  $$_{name} = "$nickname." . sprintf '%.0f', $$_{size}/1000;
  if ($$_{trna}) {$$_{name} .= $$_{tRNA_aa}} else {$$_{name} .= $isles[$$_{distIsle}]{tRNA_aa}}
  $$_{pot_frags} = join(',', @{$tandem_gps{$$_{group}}{frag}});
  $$_{pot_ints}  = join(',', @{$tandem_gps{$$_{group}}{ints}});
  $$_{group} .= $tandem_gps{$$_{group}}{sum};
  $$_{side} = "";
  $serial++;
  for my $int (keys %{$$_{ints}}) {$$_{intCoords} .= "$$_{ints}{$int}{L}..$$_{ints}{$int}{R},"}
  for my $int (keys %{$$_{ints}}) {$$_{intList} .= "$$_{ints}{$int}{annot},"}
  $out = "$$_{dna}\tisland_finder\tgenomic_island\t$$_{L}\t$$_{R}\t$$_{overall}\t$$_{strand}\t.\tID=$$_{dna}.$serial;";
  for my $attribute (qw/trna int_site int_site_type trna_dupe tRNA_L tRNA_R tRNA_len tRNA_aaa tRNA_aa A_site J_site questionable qStart qEnd hitStart hitEnd percent_id bit_score/) {
   if ($$_{trna}) {$out .= "$attribute=$$_{$attribute};"}
   else {$out .= "$attribute=$isles[$$_{distIsle}]{$attribute};"}
  }
  for my $attribute (qw/size side group segment origin proximal distal isleLseq isleRseq delta_int foreign housekeep hypoth delta_GC dinuc intList intCoords pot_ints pot_frags name/) {
   $out .= "$attribute=$$_{$attribute};";
  }
  print GFF "$out\n";
  #if ($overlapped) {print OVP $out . "reject=$$_{reject};\n"} else {print GFF "$out\n"}
 }
 close GFF;
 #close OVP;
 if ($verbose) {
  print "Later stage rejections\n";
  my $sum = 0;
  for (sort keys %reasons) {$sum += $reasons{$_}; print "$reasons{$_} got $_: runningTotal=$sum\n"}
  print "$serial final islands\n";
 }
}

sub getSeq {
 my ($coord, $ori, $dupe) = @_;
 my ($acc1, $L1, $R1) = ($1, $2, $3) if $coord =~ /^([^\/]+)\/(\d+)-(\d+)/;
 my ($acc2, $L2, $R2) = ($1, $2, $3) if $coord =~ /([^\/]+)\/(\d+)-(\d+)$/;
 my ($L, $R) = ($L1, $R2);
 if ($dupe == 1) { if ($ori eq '+') { $R ++ } else { $L -- } }
 else { if ($ori eq '+') { $L -- } else { $R ++ } }
 my $Lseq = `perl $path/collectSeq.pl -i $prefix.fa -e $acc1 -L $L -R $L -f 100 -s`;
 my $Rseq = `perl $path/collectSeq.pl -i $prefix.fa -e $acc2 -L $R -R $R -f 100 -s`;
 for ($Lseq, $Rseq) {chomp; $_ = Revcomp($_) if $ori eq '-'; $_ = uc $_}
 if ($ori eq '-') {($Lseq, $Rseq) = ($Rseq, $Lseq)}
 my @Lf = split //, substr($Lseq, 0, 100, '');
 my @Rf = split //, substr($Rseq, 0, 100, '');
 my @Lb = split //, $Lseq;
 my @Rb = split //, $Rseq;
 for my $i (reverse 0..99) {last if $Lf[$i] ne $Rf[$i]; $Lf[$i] = lc $Lf[$i]; $Rf[$i] = lc $Rf[$i]}
 for my $i (0..100) {last unless $Lb[$i] and $Rb[$i]; last if $Lb[$i] ne $Rb[$i]; $Lb[$i] = lc $Lb[$i]; $Rb[$i] = lc $Rb[$i]}
 my ($isleL, $isleR) = (join('', @Lf, @Lb), join('', @Rf, @Rb));
 #warn "isleLseq: $isleL\nisleRseq: $isleR\n";
 return ($isleL, $isleR);
}

sub Revcomp {my $seq = reverse $_[0]; $seq =~ tr/acgtACGTbdhkmrvyBDHKMRVY/tgcaTGCAvhdmkybrVHDMKYBR/; return $seq}

sub Options {
 my $version = '0.2 (Nov 2016)';
#   '    |    '    |    '    |    '    |    '    |    '    |    '    |    '    |
my $help = <<END;
$scriptname version $version
Usage: perl $scriptname [options] SINGLE_FASTA_FILE
  -maxSize:  Maximum island size. Default: $maxSize.
  -minSize:  Minimum island size. Default: $minSize.
  -nickname: Brief name for genome (as might be used to start a locus_tag).
  -criterion: Basis for overlap resolution, 3 options: random, score (7-test false positive formula), deltaGC. Default = $criterion.
  -cross: Three options: intact, cross, or circleOrigin. Default: $cross.
  Additional options: -help, -version, -verbose, -authors, -license

Example: perl $scriptname -verbose SINGLE_FASTA_FILE

END
 my $authors = "AUTHORS: Kelly Williams (kpwilli\@sandia.gov), Corey Hudson, Britney Lau, Owen Solberg\n";
 my $license = "LICENSE AND COPYRIGHT: Copyright (2013) Sandia Corporation. Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains certain rights in this software.\n";

 die $help if @ARGV == 0;
 use Getopt::Long;
 my $options_okay = GetOptions(
  'help' => sub {print $help; exit},
  'version' => sub {print "$scriptname version $version\n"; exit},
  'authors' => sub {print $authors; exit},
  'license' => sub {print $license; exit},
  'verbose' => sub {$verbose = "--verbose"},
  'criterion=s' => \$criterion,
  'nickname=s' => \$nickname,
  'cross=s' => \$cross,
  'maxSize=i' => \$maxSize,
  'minSize=i' => \$minSize,
 ) ;
 die $help if !$options_okay;
 die "Only criteria allowed: random, score, deltaGC\n" unless $criterion =~ /^random|score|deltaGC$/;
}
