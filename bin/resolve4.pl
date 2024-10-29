#! /usr/bin/perl
use strict; use warnings;
use List::Util qw(shuffle);
use File::Spec;
use Getopt::Long;

# INITIALIZATIONS
my ($cutoff, $maxIsle, $criterion, $verbose, $logic, $prefix, $tigerfile, $islanderfile, $output) = (0.51, 200000, 'bestSupp', 1, 'lenient', 'genome');
warn "Called '$0 @ARGV' on " . localtime . "\n" if $verbose;
my $dir = File::Spec->rel2abs($0); $dir =~ s/\/([^\/]+)$//; my $scriptname = $1;
Options();
my $minIsle = ($logic) ? 5000 : 2000;
my ($nick, %isl, %ends, $endorder, %tandems, %seen, %tandemtypes, %ct, %rejects, %splits, %dnas, %draws, %segCts, %trnas, %uniqIsles, %contexts, $org, %gnm, %tandCur, %serials, %rrna);
my (@dnaRA, %hskpEnrich, %fornEnrich, %stats, $in, %prevScores, %oldprevs, @order, %intList);
my (%is607, %supconfl, %skipInt, %litProv);
my %aalookup = qw/Ala A Arg R Asn N Asp D Cys C Glu E Gln Q Gly G His H Ile I Ile2 J Leu L Lys K Met M Phe F Pro P Ser S Thr T Trp W Tyr Y Val V Pyl O SeC U SeC(p) U tmRNA Z iMet B fMet B Sup X Undet ?/;

# LOAD DATA
IntList();
LoadPrevScores();
LoadDNA();
LoadTax();
LoadHousekeepEnrich();

# RESOLUTION
if ($islanderfile) {
 LoadIslands($islanderfile, "Islander"); # Load Islander output into isl hash
} # Creating the Islands for Islander data

if ($tigerfile) {
 LoadIslands($tigerfile, "TIGER"); # Load TIGER output into isl hash
} # Creating the Islands for TIGER data
Ends(); # Loading in ends from all sources
if ($islanderfile) {
 FreshTandem();
}

if ($tigerfile and $islanderfile) {
 while (1) {my $change = TandemBuildTrna(); last if $change == 0}
} # Combining Tandems that are overlapping from TIGER and Islander

if ($tigerfile) {
 TandemBuildCompOnly(); # Create new tandems for TIGER islands that cannot fit in the already existing Islander tandems
}  # Treat islands that didn't fit into Islander Tandem

WriteTandems(); exit();

ReInt(); # Loading TIGER Islands into prevScores
if ($tigerfile and $islanderfile) {
 EndTrna(); # For all ends, add id of trna that the end is inside of
 my $ttot = 0; for (keys %tandems) {$ttot += keys %{$tandems{$_}}} warn "$ttot tandems\n" if $verbose;
 TandemSplit(); # During tandem ends into candidate spans
 TandemDeoverlap('multi'); # Reject lower support tandem (the whole tandem) if it is overlapped. ISSUE: We are questioning if this is the right thing to do
} #Tandem Resolution when for TIGER,ISLANDER Islands

# Resolution of Tandems for all data
TandemSplit(); # During tandem ends into candidate spans
TandemDeoverlap('multi'); # # Reject lower support tandem (the whole tandem) if it is overlapped. ISSUE: We are questioning if this is the right thing to do
SpanDeoverlap(); # Deoverlap spans within tandems, may be able to take over TandemDeoverlap function
Write(); # Write the output and prevScores
warn "Resolve Succeeded!\n";

# SUBROUTINES
sub LoadHousekeepEnrich {
 for (@{ReadFile("$dir/../db/housekeep_enrich.txt")}) {$hskpEnrich{$1} = $2 if /(\S+)\t(\S+)/}
 for (@{ReadFile("$dir/../db/foreign_enrich.txt")})   {$fornEnrich{$1} = $2 if /(\S+)\t(\S+)/}
} # Load housekeeping and foreign gene data

sub LoadTax {
 $org = "division=;phylum=;order=;class=;family=;genus=;species=;prg=;";
 return unless -f "genome.tax";
 $org = "";
 my @taxa;
 for (`cat genome.tax`) {
  chomp; my ($nickname, $spp, $rank, $code) = split "\t";
  @taxa = split ';', $rank;
  $nick = $nickname if $nickname and not $nick;
 }
 for (qw/division phylum order class family genus species org/) {$org .= "$_=" . shift(@taxa) . ";"}
} # Load tax data if available

sub LoadTrna {  # Change to marking oneLetter+q for questionable?
 open IN, "trna/ttm.gff" or die "Can't open ttm.gff\n";
 while (<IN>) {
  my @f = split "\t";
  my $id = 'Z';
  $id = $aalookup{$1} if /aa=([^;]+)/;
  $id = '?' if /questionable=[^;]/;
  push @{$trnas{$f[0]}}, {L => $f[3], R => $f[4], id => 1 + scalar(@{$trnas{$f[0]}}) . $id};
 }
 close IN;
} # Load trna data from ttm.gff
# ISSUE: New Pseudo type of trna gene classification is not properly treated

sub LoadDNA {
 for (@{ReadFile("$prefix.stats")}) {
  next unless s/^(\S+)\t//;
  my $dna = $1;
  for my $cat (qw/len circ replicon cds pfam hypoth hskp forn/) {$stats{$dna}{$cat} = $1 if s/^(\S+)\t//}
  @{$stats{$dna}{dnaRA}} = ($_[0], split "\t"); # relative abundances of C and key dinucleotides
  @{$stats{$dna}{annots}} = ();
 }
 for (@{ReadFile("$prefix.gff")}) {
  my @f = split "\t";
  my %l = (dna => $f[0], source => $f[1], type => $f[2], L => $f[3], R => $f[4], supp => $f[5], orient => $f[6], line => $_);
  for (split ';', $f[8]) {$l{$1} = $2 if /^([^=]+)=([^;]+)/}
  $l{id} = $l{ID}; $l{id} = $l{annot} if $l{annot};
  $l{pfam} = 'hypothetical'; $l{pfam} = $l{pfam1} if $l{pfam1};
  push @{$stats{$f[0]}{annots}}, \%l;
 }
 my $dna = '';
 for (@{ReadFile("genome.fa")}) {
  $dna = $1 , next if /^>(\S+)/;
  chomp;
  $stats{$dna}{seq} .= $_;
 }
} # Load in data from genome.gff, including ints and sequences

sub LoadPrevScores {  # Because score calculation is slow
 return;
 return unless (-e "prevScores.gff");
 for (@{ReadFile("prevScores.gff")}) {
  my @f = split "\t";
  #die unless /coord=([^;]+);/;
  my $coord = $1;
  my $p = \%{$prevScores{$coord}};
  for my $cat (qw/coord deltaside len deltaint foreign housekeep hypoth delta_GC dinuc overall logscore/)
  {die "$f[8]" unless $f[8] =~ /$cat=([^;]*);/; $$p{$cat} = $1}
  die unless ($f[8] =~ /ints=([^;]*);/); my $ints = $1;
  for (split(",", $ints)) {
   $$p{ints}{$_} = 1;
  }
  ($$p{overall}, $$p{logscore}) = Score($$p{len}, $$p{deltaint}, $$p{foreign}, $$p{housekeep}, $$p{hypoth}, $$p{delta_GC}, $$p{dinuc});
  $$p{line} = join("\t", @f);
 }
} # Load previously calculated scores
# ISSUE: probably need to delete old prevScores with this rewrite

sub ReadFile {
 my @ret;
 print "Reading file $_[0]\n" if $verbose;
 unless (open IN, $_[0]) {warn "Can't open $_[0]\n"; return \@ret;}
 while (<IN>) {next if /^#/; push @ret, $_}
 close IN; chomp @ret; return \@ret;
}

sub IntList {
 for (`cat genome.gff`) {
  my @f = split "\t";
  push @{$rrna{$f[0]}}, [$f[3], $f[4], $1] if /\trRNA\t.*product=(23S|16S) ribosomal RNA/;
  my $annot;
  if ($logic) {
   next unless $f[8] =~ /annot=(Y-Int[^;]*|S-Int[^;]*|Tnp[^;]*);/; $annot = $1;
  } else {
   next unless $f[8] =~ /annot=([^;]+);/; $annot = $1;
  }
  next if $annot =~ /Xer[^;]*|Integron-Int[^;]*|Hom_end_hint[^;]*/;
  @{$intList{$f[0]}{$annot}}{qw/id L R mid/} = ($annot, $f[3], $f[4], int(($f[3]+$f[4])/2));  # Gene midpoint
 }
} # Create a list of annots for all dna

sub EndCt {
 my $ct = 0;
 for my $dna (sort keys %ends) {$ct += scalar keys %{$ends{$dna}}}
 warn "$ct ends\n" if $verbose;
} # Count ends

sub EndTrna {
 for my $dna (sort keys %ends) {
  for my $end (keys %{$ends{$dna}}) {
   for (@{$trnas{$dna}}) {next if $$_{R} < $end; last if $$_{L} > $end; $ends{$dna}{$end}{trna}{$$_{id}} =1}
  }
 }
} # For all ends, add id of trna that the end is inside of

sub LoadIslands {
 my ($file, $mode) = @_;
 warn "$mode Input File: $file\n" if $verbose;
 for (@{ReadFile($_[0])}) {
  chomp; next if /replicon=Plm/;
  my $line = $_;
  my $coord; $coord = $1 if /coord=([^;]+);/;
  my @f = split "\t";
  $f[8] .= "compose=simple;" unless /compose=([^;]+);/; $coord = ($f[6] eq '+') ? "$f[0]/$f[3]-$f[4]" : "$f[0]/$f[4]-$f[3]" unless $coord; # Backwards compatibility
  my %l = (source => 'TIGER', supp => $f[5], line => $line, f8 => $f[8], rRNA => '', trnaSide => '');
  if ($mode eq "Islander") {
   @l{qw/source score supp trnaSide reject_type/} = ('Islander', $f[5], 0, 'L', 'final');
  }
  my $inverter = 1;
  for (split(/\+/, $coord)) {
   die unless /([^\/]+)\/(\d+)-(\d+)/;
   my ($dna, $start, $stop) = ($1, $2, $3); #warn "$dna, $start, $stop";
   my $orient = $inverter * ($stop-$start);
   $orient = ($orient > 0) ? '+' : '-';
   my ($L, $R) = sort {$a <=> $b} ($start, $stop);
   push @{$l{contigs}}, {dna => $dna, start => $start, stop => $stop, L => $L, R => $R, orient => $orient};
   $inverter = -1; # Flip orientation of second half so it points to terminal end
  } # Storing dna, orient, L, and R so they are specific to each contig
  my %flag;
  $l{len} = 0;
  for my $contig (@{$l{contigs}}) {
   for (@{$rrna{$$contig{dna}}}) {
    last if $$_[0] > $$contig{R}; next if $$_[1] < $$contig{L};
   $flag{$$_[2]} ++;
   }
   $l{len} += $$contig{R} - $$contig{L} + 1;
  }
  print "Rejecting $coord, rRNA operon\n", next if $flag{'23S'} and $flag{'16S'};
  for (split ';', $f[8]) {$l{$1} = $2 if /^([^=]+)=([^;]+)/}
  if ($mode eq "Islander") {
   next if $l{origin} == 1; #ISSUE: Old code says to manually deal with this in question marks, not sure if it every did
   $l{brief} = sprintf('%0.f', $l{len}/1000) . '.' . $l{tRNA_aa};
   $l{db} = 'NA';
   $l{gnmok} = $l{ID}; $l{gnmok} =~ s/\.[^\.]+$//;
   $l{target} = $l{tRNA_aa};
   $nick = $l{gnmok} unless $nick;
   $l{reject} = '';
  } else {
   $l{score} = $l{compScore};
   $serials{$l{brief}} ++ if $mode eq "TIGER";
   $nick = $l{ID} unless $nick; $nick =~ /^([^\.]+)/; $nick = $1;
   $l{ID} = "$nick.$l{brief}.$serials{$l{brief}}";
   $l{target} = $l{brief}; $l{target} =~ s/[^\.]+\.//;
   $gnm{bestSupp} = $f[5] unless $gnm{bestSupp} and $gnm{bestSupp} > $f[5];
  }
  next if $rejects{$l{ID}};
  if ($l{compose} eq "circleJxn" and $l{contigs}[0]{L} == $l{contigs}[-1]{R}) {
   print "Rejecting $coord, cicleJxn across whole contig\n", next;
  }
  if ($mode eq "TIGER") {
   unless ($f[8] =~ /OL=(\d+)-(\d+);OR=(\d+)-(\d+);/) {print "no OL or OR for $_\n"; next}
   ($l{OLL}, $l{OLR}) = ($1, $2);
   ($l{ORL}, $l{ORR}) = ($3, $4);
  } else {
   my ($coordside, $dir, $contig, $otherside) = ($l{contigs}[0]{L}, 1, 0, 'R'); # Start calculating tRNA query hit genomic coordinates
   if (abs($l{hitStart} - $l{contigs}[0]{start}) < abs($l{hitStart} - $l{contigs}[-1]{stop})) {$coordside = $l{contigs}[-1]{R}; $l{trnaSide} = 'R'; $otherside = 'L'; $contig = -1} # Hit closer to L (so tRNA closer to R)
   $dir = -1 if $l{contigs}[$contig]{orient} eq '-';
   ($l{'O'.$l{trnaSide}.'L'}, $l{'O'.$l{trnaSide}.'R'}) = sort {$a <=> $b} ($coordside + $dir*($l{qEnd} - $l{int_site}), $coordside + $dir*($l{qStart} -$l{int_site}));
   ($l{'O'.$otherside.'L'}, $l{'O'.$otherside.'R'}) = sort {$a <=> $b} ($l{hitStart}, $l{hitEnd});
  }
  $l{coord} = $coord unless $l{coord};
  $ct{islct}++;
  $isl{$ct{islct}} = \%l;
 }
} # Load Tiger and Islander data

sub Ends {
 for my $isle (sort {$isl{$a}{source} cmp $isl{$b}{source} || $isl{$b}{supp} <=> $isl{$a}{supp} || $isl{$b}{score} <=> $isl{$a}{score}} keys %isl) {
  my %i; @i{qw/id L R dnaL dnaR supp LL LR RL RR/} = ($isle, $isl{$isle}{contigs}[0]{start}, $isl{$isle}{contigs}[-1]{stop}, $isl{$isle}{contigs}[0]{dna}, $isl{$isle}{contigs}[-1]{dna}, $isl{$isle}{supp}, $isl{$isle}{OLL}, $isl{$isle}{OLR}, $isl{$isle}{ORL}, $isl{$isle}{ORR});
  for my $side (qw/L R/) {
   my $dna = $i{'dna'.$side};
   my ($expand, $done);
   $i{$side} = $i{$side.'L'} if $i{$side} < $i{$side.'L'};  # Adjust occasional cases (eg, Tsa1.50D) when end called outside overlap window
   $i{$side} = $i{$side.'R'} if $i{$side} > $i{$side.'R'};
   for (sort {$a <=> $b} keys %{$ends{$dna}}) {
    my $end = $ends{$dna}{$_};
    next if $i{$side.'L'} > $$end{R};
    last if $i{$side.'R'} < $$end{L};
    if ($i{$side.'L'} < $$end{L}) {$$end{L} = $i{$side.'L'}; $expand = $_}
    if ($i{$side.'R'} > $$end{R}) {$$end{R} = $i{$side.'R'}; $expand = $_}
    $$end{members}{$i{id}.$side} = $dna;
    $$end{supp} += $i{supp};
    $i{$side} = $$end{coord};
    $done ++;
   }
   unless ($done) { # new end
    $endorder ++;
    %{$ends{$dna}{$i{$side}}} = (L => $i{$side.'L'}, R => $i{$side.'R'}, order => $endorder, nominal => $i{$side}, coord => $i{$side}, founder => $i{id}.$side, members => {$i{id}.$side => $dna}, supp => $i{supp});
   }
   next unless $expand;
   my $to = $ends{$dna}{$expand};
   for (sort {$a <=> $b} keys %{$ends{$dna}}) {
    next unless $ends{$dna}{$_};
    my $from = $ends{$dna}{$_};
    next if $$to{L} > $$from{R};
    last if $$to{R} < $$from{L};
    next if $$from{coord} eq $expand;  # Skip self
    if ($$from{order} < $$to{order}) { # Merge into older end
     ($from, $to) = ($to, $from);
     $expand = $$to{coord};
     $i{$side} = $$to{coord};
    }
    for (keys %{$$from{members}}) {
     $$to{members}{$_} = $$from{members}{$_};
     die "$dna $_ parse" unless /(.*)([LR])$/;
     $isl{$1}{'end'.$2} = $$to{coord};
    }
    if ($$from{L} < $$to{L}) {$$to{L} = $$from{L}}
    if ($$from{R} > $$to{R}) {$$to{R} = $$from{R}}
    delete $ends{$dna}{$$from{coord}};
   }
  }
  if ($isl{$isle}{compose} eq "simple") {
   @{$isl{$isle}{contigs}[0]}{qw/endL endR/} = sort {$a <=> $b} ($i{L}, $i{R});
  } else {
   @{$isl{$isle}{contigs}[0]}{qw/endL endR/} = ($i{L}, $isl{$isle}{contigs}[0]{stop});
   @{$isl{$isle}{contigs}[1]}{qw/endL endR/} = ($isl{$isle}{contigs}[1]{start}, $i{R});
  }
 }
 EndCt();
} # Load all ends from islands, excluding terminal ones, with their overlapping regio. If any overlapping regions overlap, merge them together. Return the final L and R coordinate from the highest supported island of the overlap region. NOTE: islands must be processed by Islander || support || score priority.

sub FreshTandem {
 warn "Fresh Tandem\n" if $verbose;
 for my $isle (keys %isl) {
  my $i = \%{$isl{$isle}};
  next unless $$i{source} eq "Islander";
  #@{$i}{qw/len deltaint deltaside ints foreign housekeep hypoth delta_GC dinuc overall logscore/} = Span($$i{coord}); ISSUE: Span not written yet, disabling for now
  my $t = \%{$tandems{$$i{group}}};
  $$t{sources}{Islander} = 1; $$t{members}{$isle} = 1;
  for my $contig (@{$$i{contigs}}) {
   $$t{ends}{$$contig{dna}}{$$contig{endL}} = 1;
   $$t{ends}{$$contig{dna}}{$$contig{endR}} = 1;
   $$t{dna}{$$contig{dna}} = 1;
  }
  $$t{questionable} = ($$i{questionable}) ? 1 : 0;
  $$t{bestSupp} = 0; $$t{trnaSide} = $$i{trnaSide}; $$t{id} = 'I' . $$i{group}; $$t{target} = $$i{tRNA_aa};
  $$t{trnaStepoverTest} = sub {return 1 if shift() < $$i{contigs}[0]{endL}};
  $$t{trnaStepoverTest} = sub {return 1 if shift() > $$i{contigs}[-1]{endR}} if $$t{trnaSide} eq 'R';
  $$t{bestScore} = $$i{logscore} unless $$t{bestScore} and $$t{bestScore} < $$i{logscore};
  $$i{tandem} = $$i{group};
 }
} # Build tandems for islander data

sub TandemBuildTrna {
 my $change = 0;
 for my $tand (keys %tandems) {
  my $t = \%{$tandems{$tand}};
  for my $dna (keys %{$$t{dna}}) {
   for my $end (keys %{$$t{ends}{$dna}}) {
    my $j = $ends{$dna}{$end};
    for (keys %{$$j{members}}) {
     die "$dna $_ $$j{L} $$j{R}" unless/(.*)([LR])$/;  # $1=island $2=side
     my ($isle, $side, $otherside) = ($1, $2, $2); $otherside =~ tr/LR/RL/;
     my $contig = ($otherside eq "L") ? 0 : -1;
     $otherside = $isl{$isle}{contigs}[$contig]{'end'.$otherside};
     next if &{$$t{trnaStepoverTest}}($otherside);  # Prevents tRNA stepover
     next if $$t{members}{$isle};
     next if $isl{$isle}{group} and $isl{$dna}{$isle}{group} ne $tand;
     $change++;
     $isl{$isle}{'end'.$side} = $end;
     $$t{bestSupp} = $isl{$isle}{supp} if $$t{bestSupp} < $isl{$isle}{supp};
     $$t{totSupp} = $isl{$isle}{supp} if $$t{bestSupp} < $isl{$isle}{supp};
     $$t{sources}{$isl{$isle}{source}} = 1;
     $$t{members}{$isle} = 1;
     for my $contig (@{$isl{$isle}{contigs}}) {
      $$t{dna}{$$contig{dna}} = 1;
      $$t{ends}{$$contig{dna}}{$$contig{endL}} = 1;
      $$t{ends}{$$contig{dna}}{$$contig{endR}} = 1;
     }
     $isl{$isle}{tandem} = $tand;
    }
   }
  }
 }
 warn "TandemBuildTrna, $change changes\n" if $verbose;
 return $change;
} # Take Tandems already created by Islander data and insert any TIGER data that matches as well

sub TandemBuildCompOnly {
 warn "TandemBuildCompOnly\n" if $verbose;
 my ($mode) = ($_[0]);
 for my $dna (keys %ends) {
  my %islands;
  for (keys %{$ends{$dna}}) {
   my ($end, $j) = ($_, $ends{$dna}{$_});
   for (@{$trnas{$dna}}) {next if $$_{R} < $$j{L}; last if $$_{L} > $$j{R}; $$j{trna}{$$_{id}} =1}
   my ($tand, %jmembers, %jtands);
   my ($bestSupp, $totSupp, %flag) = (0,0);
   for (keys %{$$j{members}}) {
    die "ERROR parsing island name $_ for end $end: $dna $_ $$j{L} $$j{R}\n" unless/(.*)([LR])$/;  # $1=island $2=side
    my ($isle, $side) = ($1, $2);
    next if $isl{$isle}{rejected_by} or $isl{$isle}{tandem};  # Skip since previously evaluated as part of Islander tandem
    $jmembers{$isle} ++;
    $jtands{$islands{$isle}} ++ if $islands{$isle};  # Tandem call for other island end?
    $isl{$isle}{'end'.$side} = $$j{coord};  # Record ends for isls
    if ($bestSupp < $isl{$isle}{supp}) {$bestSupp = $isl{$isle}{supp};}
    $totSupp += $isl{$isle}{supp}
   }
   next unless keys %jmembers;
   if (keys %jtands == 0) {
    $$j{founder} =~ /(.*)([LR])$/; $tand = $1;
   } else {
    $tand = (keys %jtands)[0]; #Choose one at random
   }
   if (keys %jtands > 1) {
    warn "Merging " . join(", ", keys %jtands) . "into $tand\n";
    for my $old (keys %jtands) {
     $totSupp += $tandems{$old}{totSupp};
     $bestSupp = $tandems{$old}{bestSupp} if $bestSupp < $tandems{$old}{bestSupp};
     next if $old eq $tand;
     for (keys %{$tandems{$old}{members}}) {
      $tandems{$tand}{members}{$_} = 1; $islands{$_} = $tand;
     }
     for my $td (keys %{$tandems{$old}{ends}}) {
      for (keys %{$tandems{$old}{ends}{$_}}) {
       $tandems{$tand}{ends}{$td}{$_} = 1;
       $tandems{$tand}{dna}{$td} = 1;
      }
     }
     warn "Deleting $old\n";
     delete $tandems{$old};
    }
   } # Add info to chosen tandem about this end (j)
   $tandems{$tand}{ends}{$dna}{$$j{coord}} = 1;
   $tandems{$tand}{bestSupp} = $bestSupp if not defined $tandems{$tand}{bestSupp} or $bestSupp > $tandems{$tand}{bestSupp};
   $tandems{$tand}{totSupp} += $totSupp;
   $tandems{$tand}{dna}{$dna} = 1;
   for my $isle (keys %jmembers) {
    $tandems{$tand}{sources}{$isl{$isle}{source}} = 1;
    $tandems{$tand}{members}{$isle} = 1;
    for my $contig (@{$isl{$isle}{contigs}}) {
     $tandems{$tand}{dna}{$$contig{dna}} = 1;
     $tandems{$tand}{ends}{$$contig{dna}}{$$contig{endL}} = 1;
     $tandems{$tand}{ends}{$$contig{dna}}{$$contig{endR}} = 1;
    }
    $islands{$isle} = $tand;
   }
  }
  for (keys %tandems) {
   my $t = $tandems{$_};
   my @dnas = sort {$a cmp $b} keys %{$$t{ends}};
   my @ends1 = sort {$a <=> $b} keys %{$$t{ends}{$dnas[0]}};
   my @ends2 = sort {$a <=> $b} keys %{$$t{ends}{$dnas[-1]}};
   $$t{id} = "C$ends1[0]-$ends2[-1]($$t{bestSupp})" unless $$t{id};
   if ($$t{sources}{TIGER} and $$t{id} =~ s/^I/CI/) {};
   $$t{questionable} = 0;
  } # ID for all of the tandems
 }
} # Build Tandems that contain TIGER data only

sub WriteTandems {
 open TAND, ">tandem_struct.txt";
 for my $tand (sort {$tandems{$b}{bestSupp} <=> $tandems{$a}{bestSupp}} keys %tandems) {
  my $t = \%{$tandems{$tand}};
  my $islnum = keys %{$$t{members}};
  print TAND ">$$t{id}\t$tand\t$$t{bestSupp}\t$islnum\n"; # 
  for my $isle (sort {$isl{$a}{source} cmp $isl{$b}{source} || $isl{$b}{supp} <=> $isl{$a}{supp} || $isl{$b}{score} <=> $isl{$a}{score}} keys %{$$t{members}}) {
   my $i = \%{$isl{$isle}};
   my $co = (defined $$i{crossover}) ? $$i{crossover} : "Islander";
   my $ints = (defined $$i{ints}) ? $$i{ints} : "NA";
   print TAND "$isle $$i{ID}:$$i{coord}|supp=$$i{supp};score=$$i{score};crossover=$co;ints=$ints\t";
  }
  print TAND "\n";
  for my $dna (sort {$a cmp $b} keys %{$$t{dna}}) {
   print TAND "$dna\t";
   for my $end (sort {$a <=> $b} keys %{$$t{ends}{$dna}}) {
    my $endpt = (defined $ends{$dna}{$end}{founder}) ? $ends{$dna}{$end}{founder} : "terminus";
    my $mem = ""; $mem = join(",",keys %{$ends{$dna}{$end}{members}}) if (defined $ends{$dna}{$end}{founder});
    print TAND "$endpt:$end($mem)\t";
   }
   print TAND "\n";
  }
 }
} # Write out the Tandem structures



sub Options {
 my $help = "Usage: $0 -logic [strict/lenient] -prefix [prefix_name] -tiger TIGER_file -islander ISLANDER_file -o output_file\nLaunch from within genome directory\n";
 my $version = "3.0 (April 2024)";
 my $authors = ""; # ISSUE: Someone fill in authors here
 my $license = "COPYRIGHT AND LICENSE: Copyright 2018, National Technology & Engineering Solutions of Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS, the U.S. Government retains certain rights in this software.
Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS 'AS IS' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.\n";
 die $help if @ARGV == 0;
 my $options_okay = GetOptions(
  'help' => sub {print $help; exit},
  'version' => sub {print "$scriptname version $version\n"; exit},
  'authors' => sub {print $authors; exit},
  'license' => sub {print $license; exit},
  'verbose' => \$verbose,
  'prefix=s' => \$prefix,
  'logic=s' => \$logic,
  'tiger=s' => \$tigerfile,
  'islander=s' => \$islanderfile,
  'o=s' => \$output
 );
 die "Usage: $0 -logic [strict/lenient] -prefix [prefix_name] -tiger TIGER_file -islander ISLANDER_file\nLaunch from within genome directory\n" if !$options_okay;
 die "ERROR: illegal logic option\n" . $help unless $logic =~ /^strict|lenient|$/;
 die "ERROR: prefix required\n" . $help unless $prefix;
}
