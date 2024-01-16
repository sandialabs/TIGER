#!/usr/bin/perl
use strict; use warnings;

# OPTIONS
my ($verbose, $gffin, $protgff, $domfiles, $domgff);
Options(); # See bottom of script; help and other messages, Getopt

# PROGRAM
my %prots;
open IN, $gffin or die "No gff file $gffin\n";
while (<IN>) {
 last if /^##FASTA/; next if /^#/; chomp; 
 next unless /\tCDS\t.*ID=([^;]+)/;
 my $id = $1;
 my @f = split "\t";
 %{$prots{$id}} = (L => $f[3], R => $f[4], orient => $f[6], dna => $f[0]);
}
close IN;

for my $domfile (split ',', $domfiles) {
 open IN, $domfile or die "No domain file $domfile\n";
 while (<IN>) {
  chomp;
  next if /^#/;
  my ($id, $fam, $hStart, $hEnd, $aStart, $aEnd, $score) = (split /\s+/)[0,3,15,16,17,18,13]; # Parse the HMMER domtbl line
  my $dna = $id; $dna =~ s/_[0-9]+$//;
  my ($L, $R) = ($prots{$id}{L} + 3*$aStart -3, $prots{$id}{L} + 3*$aEnd -1); # Genome coords for domain (differing if negative orientation, below)
  ($L, $R) = ($prots{$id}{R} - 3*$aEnd +1, $prots{$id}{R} - 3*$aStart +3) if $prots{$id}{orient} eq '-';
  push @{$prots{$id}{doms}}, {fam => $fam, score => $score, range => "$aStart..$aEnd", hrange => "$hStart..$hEnd", L => $L, R => $R};
 }
}

my @doms;
open OUT, ">$protgff"; # Per-protein
for my $id (sort keys %prots){
 print OUT join("\t", $prots{$id}{dna}, 'pfamA', 'CDS', $prots{$id}{L}, $prots{$id}{R}, '.', $prots{$id}{orient}, '.', "ID=$id;");
 if ($prots{$id}{doms}){
  my $i = 0;
  for my $d (sort {$$b{score} <=> $$a{score}} @{$prots{$id}{doms}}){ # pfam1 has highest score
   $i++;
   print OUT "pfam$i=$$d{fam};score$i=$$d{score};hmmRange$i=$$d{hrange};range$i=$$d{range};";
   push @doms, {id => $id."_dom$i", dna => $prots{$id}{dna}, L => $$d{L}, R => $$d{R}, fam => $$d{fam}, score => $$d{score}, orient => $prots{$id}{orient}};
  }
 }
 print OUT "\n";
}
close OUT;

open OUT, ">$domgff"; # Per-domain
for my $d (sort {$$a{dna} cmp $$b{dna} || $$a{L} <=> $$b{L}} @doms) {
 print OUT join("\t", $$d{dna}, 'pfamA', 'CDS', $$d{L}, $$d{R}, $$d{score}, $$d{orient}, '.', "ID=$$d{id};pfam=$$d{fam};"), "\n";
}
close OUT;

sub Options {
 my $scriptname = $0; $scriptname =~ s/.*\///;
 my $version = '1.0 (Nov 2016)';
#   '    |    '    |    '    |    '    |    '    |    '    |    '    |    '    |
 my $help = <<END;
$scriptname version $version
Usage: perl $scriptname -gff <file> -dom <file> -protgff <file> -domgff <file>
 -gffin:    Prokka .gff file.
 -dom:      Pfam domain table.
 -protgff:  Output protein gff file.
 -domgff:   Output domain gff file.
 Additional options: -help, -version, -verbose, -authors, -license

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
  'gffin=s' => \$gffin,
  'dom=s' => \$domfiles,
  'protgff=s' => \$protgff,
  'domgff=s' => \$domgff,
 );
}
