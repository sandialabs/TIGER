#! /usr/bin/perl
use strict; use warnings;
use 5.10.0;
use File::Spec;
use List::MoreUtils 'true'; 

my ($filTot, @filname, @fil, $pici, @pici, @pici_pfam, @phage_structure, @phage_funct, %l, @phage, @funct, @struct, $Stot, $Ftot, $nick, @sname, @fname, @zot, $Ptot, @pname);
my ($phagecall) = (0);
die "Usage: perl $0 FastaFile\n" unless @ARGV >= 1;
my $dir = File::Spec->rel2abs($0); $dir =~ s/\/([^\/]+)$//;
my $file = $ARGV[0]; 
die unless $file =~ /(\S+)\.[^\.]+$/;
$nick = $1; 

system("perl $dir/../tater.pl  $file") unless -f "$nick.gff";

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
for (`cat $dir/../db/PICI.txt`){
 chomp;
 next if /^#/;
 push @pici_pfam,  $_;
}

my @fil_pfam = qw(CTX_RstB Phage_Coat_A Phage_Coat_B Phage_Coat_Gp8 Zot); 


for (`cat $nick.gff`){
 chomp; 
 next unless /\S+\t\S+\t\S+\t\S+\t\S+\t\S+\t\S+\t\S+/;
 my @f = split "\t";  
 next if ($f[3] eq 'misc_RNA');
 %l = (dna => $f[0], type => $f[2], L => $f[3], R => $f[4], dir => $f[6]);
 while ($f[8] =~ s/^([^=;]+)=([^;]+);//g) {$l{$1} = $2};
 if ($l{pfam1} ~~ \@phage_structure) {push @struct, $l{ID}; push @sname, $l{pfam1}}
 elsif ($l{pfam1} ~~ @phage_funct) {push @funct, $l{ID}; push @fname, $l{pfam1}}
 if ($l{pfam1} ~~ @pici_pfam) {push @pici, $l{ID}; push @pname, $l{pfam1}}
 if ($l{pfam1} ~~ @fil_pfam){push @fil, $l{ID}; push @filname, $l{pfam1}}
}
$Stot = scalar @struct; 
$Ftot = scalar @funct;
$Ptot = scalar @pici;  
$filTot = scalar @fil;

my $len = 0; for (`cat $file`) {next if /^>/; chomp; $len += length $_;}
if ($filTot >= 1 and $len < 13000) {$phagecall = 4}
elsif ($len < 20500 and $Ptot >= 3){$phagecall = 5} 
elsif ($Stot >=1 and $Ftot >= 1) {$phagecall = 1}
elsif ($Stot >=1) {$phagecall = 2}
elsif ($Ftot >= 1) {$phagecall = 3}

open (OUT, "> FilPICIphage.txt");
print OUT "#nick\tstructTot\tnonStructTot\tPICItot\tphagecall\tstructural\tnonStruct\tPICI\n";
$nick =~ s/.*\///;
print OUT "$nick\t$Stot\t$Ftot\t$Ptot\t$phagecall\t", join(',', @sname), "\t", join(',', @fname), "\t", join(',', @pname), "\n";
close (OUT);



