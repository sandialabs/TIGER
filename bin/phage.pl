#! /usr/bin/perl
use strict; use warnings;
use 5.10.0;
use File::Spec;
use List::MoreUtils 'true'; 

my (@phage_structure, @phage_funct, %l, @phage, @funct, @struct, $Stot, $Ftot, $nick, $name, @sname, @fname, @zot);
my ($phagecall) = (0);
die "Usage: perl $0 FastaFile\n" unless @ARGV >= 1;
my $dir = File::Spec->rel2abs($0); $dir =~ s/\/([^\/]+)$//;
my $file = $ARGV[0];
$file = File::Spec->rel2abs($file);
die unless $file =~ /([^\/]+)\/([^\/]+)\.[^\.]+$/;
$name = $1;
$nick = $2; 

system("perl $dir/tater.pl $file") unless -f "$nick.gff";

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

for (`cat $nick.gff`){
 chomp; 
 next unless /\S+\t\S+\t\S+\t\S+\t\S+\t\S+\t\S+\t\S+/;
 my @f = split "\t";  
 next unless $f[2] eq 'CDS';
 %l = (dna => $f[0], type => $f[2], L => $f[3], R => $f[4], dir => $f[6]);
 while ($f[8] =~ s/^([^=;]+)=([^;]+);//g) {$l{$1} = $2};
 if ($l{pfam1} ~~ \@phage_structure) {push @struct, $l{ID}, push @sname, $l{pfam1}}
 elsif ($l{pfam1} ~~ @phage_funct) {push @funct, $l{ID}, push @fname, $l{pfam1}}
 $Stot = scalar @sname; 
 $Ftot = scalar @fname;
}


my $len = 0; for (`cat $file`) {next if /^>/; chomp; $len += length $_;}
if (grep {$_ eq 'Zot'} @fname and $len < 13000) {$phagecall = 4} 
elsif ($Stot >=1 and $Ftot >= 1) {$phagecall = 1}
elsif ($Stot >=1) {$phagecall = 2}
elsif ($Ftot >= 1) {$phagecall = 3}


open (OUT, "> phage.txt");
print OUT "#nick\tstructTot\tnonStructTot\tphagecall\tstructural\tnonStruct\n";
$name =~ s/.*\///;
print OUT "$name\t$Stot\t$Ftot\t$phagecall\t", join(',', @sname), "\t", join(',', @fname), "\n";
close (OUT);



