#! /usr/bin/perl
use strict; use warnings;
use Cwd 'abs_path';
my $dir = abs_path($0); $dir =~ s/\/[^\/]+$//;

die "Usage: perl $0 Island File\n" unless @ARGV == 1;
for (`cat $ARGV[0]`) {
 chomp;
 my @f = split "\t";
 /ID=([^;]+);/;
 my $name = $1;
 $name =~ s/\|//g;
($f[3], $f[4]) = ($1, $2) if /endL=([^;]+).*endR=([^;]+)/;
 my $cmd = "$dir/collectSeq.pl -h -i genome.fa -e $f[0] -L $f[3] -R $f[4] -n $name";
 warn "$cmd\n";
 my $fa = `$cmd`;
 print $fa;
}
