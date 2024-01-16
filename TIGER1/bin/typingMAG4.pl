#! /usr/bin/perl
use strict; use warnings;
use 5.10.0;
use File::Spec;

die "Usage: perl $0 Island File\n" unless @ARGV == 1;
my $study = $ARGV[0]; 
my $dir = File::Spec->rel2abs($0); $dir =~ s/\/[^\/]+$//;

my (%isles, %cats, %phage, %ice, %transposases, %invRep, $topSupp, $prefix, $stripfix);

system "(mkdir Isles)"; 
system "perl $dir/gff2faMAG.pl $study > isles.fa";

for (`cat isles.fa`) {
 if (/>([^\s\|\+\']+)/) {#close OUT; 
  my $nick = $1;
  my $head = substr($nick, 0, 18);
 system ("mkdir Isles/$nick"); 
 open OUT, ">Isles/$nick/$nick.fa"}
 print OUT $_;
}


open OUT, ">job";
for (glob "Isles/*"){
 next unless /(\S+)/;
 $prefix = $1;
 $prefix =~ s/Isles\///; 
 $stripfix = $prefix; $stripfix =~ s/\.x$//; 
 my $head = substr($prefix, 0, 15); 
 print OUT "cd Isles/$prefix; perl $dir/iceMAG.pl $stripfix.fa; cd ../../;\n"; 
 print OUT "cd Isles/$prefix; perl $dir/phageMAG.pl $stripfix.fa; cd ../../;\n"; 
 if ($prefix =~ /\.x$/) {next} else {system "perl -pi -e 's/>\\S+/>$head/' Isles/$prefix/$prefix.fa; \n";}
}
close OUT; 

system "bash job"; 
system "cat Isles/*/phage.txt > phage.txt"; 
system "cat Isles/*/ICE.txt > ICE.txt"; 


for (`cat phage.txt`) {chomp; my @f = split "\t"; $phage{$f[0]} = $f[3]}
for (`cat ICE.txt`) {chomp; my @f = split "\t"; $ice{$f[0]} = $f[9]}

my (%uniqIsles, %gnmnames);
my (%intfams, @outIsles, $lastgnm);
open OUT, ">islesFinal.gff";
for (`cat $study`) {
 my @f = split "\t";
 my %i; while ($f[8] =~ s/^([^=;]+)=([^;]*);//g) {$i{$1} = $2}
 die unless $i{ID} =~ /^([^\.]+)/; my $nick = $1;		
 my ($name, %intcats, $intsum, $int, $longest) = ($i{ID});
 for (split ',', $i{int}) {
  die unless s/:(\d+)-(\d+)$//;
  my $len = (abs($1-$2)+1)/3;
  $longest = $len unless $longest and $longest > $len;
  if (/Y-Int/) {$intcats{Tyr}++; $int='Tyr'} else {s/\..*//; $intcats{$_}++; $int = $_}
 }
 $i{intLongest} = $longest;
 $name =~ s/\.\d+$//;
 $name =~ s/^(([^\.]+)\.[^\.]+\.)//;
 my $id = $i{ID};
 $id =~ s/\|//g;
 $id =~ s/\'//g;
 $id =~ s/fa//;
 $id =~ s/\?/trna/;
 $id =~ s/\s+//;
 $id =~ s/\,/x/g;
 $id =~ s/\(/x/g;
 $id =~ s/\)/x/g;
 $lastgnm = $2;
 $uniqIsles{$id} ++;
 if ($uniqIsles{$id} == 2) {for (@outIsles) {last if s/\tID=\Q$id\E;/\tID=$id.1;/}}
 if ($uniqIsles{$id} > 1) {$id .= ".$uniqIsles{$id}"}
 for (qw/target gnm OL OU OR/) {if ($isles{$i{ID}}{$_ . 'new'}) {$i{$_ . 'old'} = $i{$_}; $i{$_} = $isles{$i{ID}}{$_ . 'new'}}}
 $i{dirOnU} = $isles{$i{ID}}{dirOnU} if $isles{$i{ID}}{dirOnU};
 $i{source} = $f[1];
 $i{ice} = ''; $i{ice} = $ice{$id} if defined $ice{$id}; 
 $i{phage} = ''; $i{phage} = $phage{$id} if defined $phage{$id};
 $id = $i{ID};  
 my $stripParenth = $id;
 $stripParenth =~ s/\([^\)]*\)//g;
 push @outIsles, join("\t", @f[0..7]) . "\tID=$stripParenth;";
 $i{type} = ''; 
   if ($i{phage} == '1' and $i{ice} == '1') {$i{type} = 'PhageICE'}
   elsif ($i{phage} == '1' and $i{ice} == '2') {$i{type} = 'PhageICE'}
   elsif (($i{phage} == '2' or $i{phage} == '3' or $i{phage} == '4') and ($i{ice} == '1')){$i{type} = 'ICE1'}
   elsif (($i{phage} == '2' or $i{phage} == '3' or $i{phage} == '4') and ($i{ice} == '2')){$i{type} = 'ICE2'}
   elsif ($i{phage} == '1'){$i{type} = 'Phage1'}
   elsif ($i{phage} == '2' or $i{phage} == '3'){$i{type} = 'Phage2'}
   elsif ($i{phage} == '4'){$i{type} = 'PhageFil'}
   elsif ($i{ice} == '1'){$i{type} = 'ICE1'}
   elsif ($i{ice} == '2'){$i{type} = 'ICE2'}
   else {$i{type} = 'other'}
 for (qw/supp intLongest compose coord crossover type phage ice source OLL OLR ORL ORR splice_site splice_site_type trna_dupe tRNA_len qStart qEnd bit_score
  context contextsum isleLseq unintSeq isleRseq gnm dirOnU OL OU OR division phylum order class family
 genus species org ints deltaside len deltaint foreign housekeep hypoth delta_GC dinuc/) {
  if (defined $i{$_}) {$outIsles[-1] .= "$_=$i{$_};"} else {$outIsles[-1] .= "$_=;"}
 }
 $outIsles[-1] .= "\n";
}
print OUT join('', @outIsles); @outIsles = ();
close OUT;
