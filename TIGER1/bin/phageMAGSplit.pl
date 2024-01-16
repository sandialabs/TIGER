#! /usr/bin/perl
use strict; use warnings;
use 5.10.0;
use File::Spec;
use List::MoreUtils 'true'; 

my (%r, @phage_structure, @phage_funct, %l, @phage, @lfunct, @lstruct, $lStot, $lFtot, @rfunct, @rstruct, $rStot, $rFtot, $nick, @lsname, @lfname, @zot,  @rsname, @rfname,);
my ($lphagecall, $rphagecall) = (0,0);
die "Usage: perl $0 FastaFile\n" unless @ARGV >= 1;
my $dir = File::Spec->rel2abs($0); $dir =~ s/\/([^\/]+)$//;
my $file = $ARGV[0]; 
die unless $file =~ /(\S+)\.[^\.]+$/;
$nick = $1; 


for (`cat $dir/../db/PhageStructure.txt`){
 chomp; 
 next if /^#/;
 push @phage_structure,  $_; 
}
for (`cat $dir/../db/OtherPhage.txt`){
 chomp;
 next if /^#/;
 push @phage_funct, $_; 
}

for (`cat L.$nick.gff`){
 chomp; 
 next unless /\S+\t\S+\t\S+\t\S+\t\S+\t\S+\t\S+\t\S+/;
 my @f = split "\t";  
 %l = (dna => $f[0], type => $f[2], L => $f[3], R => $f[4], dir => $f[6]);
 while ($f[8] =~ s/^([^=;]+)=([^;]+);//g) {$l{$1} = $2};
 if ($l{pfam1} ~~ \@phage_structure) {push @lstruct, $l{ID}, push @lsname, $l{pfam1}}
 elsif ($l{pfam1} ~~ @phage_funct) {push @lfunct, $l{ID}, push @lfname, $l{pfam1}}
 $lStot = scalar @lsname; 
 $lFtot = scalar @lfname;
}

if (grep {$_ eq 'Zot'} @lfname) {$lphagecall = 4}
elsif ($lStot >=1 and $lFtot >= 1) {$lphagecall = 1}
elsif ($lStot >=1) {$lphagecall = 2}
elsif ($lFtot >= 1) {$lphagecall = 3}

for (`cat R.$nick.gff`){
 chomp;
 next unless /\S+\t\S+\t\S+\t\S+\t\S+\t\S+\t\S+\t\S+/;
 my @f = split "\t";
 %r = (dna => $f[0], type => $f[2], L => $f[3], R => $f[4], dir => $f[6]);
 while ($f[8] =~ s/^([^=;]+)=([^;]+);//g) {$r{$1} = $2};
 if ($r{pfam1} ~~ \@phage_structure) {push @rstruct, $r{ID}, push @rsname, $r{pfam1}}
 elsif ($r{pfam1} ~~ @phage_funct) {push @rfunct, $r{ID}, push @rfname, $r{pfam1}}
 $rStot = scalar @rsname;
 $rFtot = scalar @rfname;
}


if (grep {$_ eq 'Zot'} @rfname) {$rphagecall = 4} 
elsif ($rStot >=1 and $rFtot >= 1) {$rphagecall = 1}
elsif ($rStot >=1) {$rphagecall = 2}
elsif ($rFtot >= 1) {$rphagecall = 3}


open (OUT, ">> phage.txt");
#print OUT "#nick\tstructTot\tnonStructTot\tphagecall\tstructural\tnonStruct\n";
$nick =~ s/.*\///;
print OUT "L.${nick}\t$lphagecall\t", "\n";
print OUT "R.${nick}\t$rphagecall\t", "\n";
close (OUT);



