#! /usr/bin/perl
use strict; use warnings;
# For setup, see bottom of script
my $path = $0;
$path =~ s/\/*[^\/]+$//;
my @releases = `ls $path | grep '^gtdb'`;
die Usage() unless @ARGV == 2;
my ($keeps_file, $gtdb_release) = @ARGV;
my @keeps = `cat $keeps_file`;
$path .= "/$gtdb_release";
die "No working directory $path found, check available gtdb releases below\n\n", Usage() unless -d $path;
chdir $path;
for (qw/labels d__Archaea.tree d__Bacteria.tree/) {die "Missing $_ file in working directory $path\n" unless -f $_}
my (%labels, %divs, %rev);
for (`cat labels`) {chomp; my ($rep, $sp, $div) = split "\t"; %{$labels{$sp}} = (rep => $rep, div => $div);}
for my $keep (@keeps) {
 chomp $keep;
 unless ($labels{$keep}) {warn "$keep not available\n"; next}
 #die "$keep $labels{$keep}{rep}\n";
 $divs{$labels{$keep}{div}}{$keep} = $labels{$keep}{rep};
}
unless (keys %divs) {die "No keepers were available, output tree empty\n"}
my $div = (sort {scalar keys %{$divs{$b}} <=> scalar keys %{$divs{$a}} } keys %divs)[0];
open OUT, ">keeps";
for my $sp (keys %{$divs{$div}}) {my $rep = $labels{$sp}{rep}; print OUT "$rep\n"; $rev{$rep} = $sp}
close OUT;
my $treefile = "$div.tree";
#die "java -jar ../../bin/PareTree.jar -t O -keep keeps -f $treefile";
`java -jar ../PareTree.jar -t O -keep keeps -f $treefile 2> /dev/null`;
$treefile =~ s/\./_pared./;
my $tree = do {local(@ARGV, $/) = $treefile; <>};
for my $rep (keys %rev) {
 $tree =~ s/([\(,])$rep([:\),])/$1$rev{$rep}$2/;
}
print $tree;


sub Usage {
 for (@releases) {chomp} my $releaselist = join(', ', @releases);
#   '    |    '    |    '    |    '    |    '    |    '    |    '    |    '    |
my $help = <<END;
This script requires two arguments:
 1. keeps_file: list of species to retain in tree. Each line's format: Escherichia__flexneri
 2. gtdb_release: currently available: $releaselist
Call like this: \$ $0 keeps_file gtdb_release > outtree
Note: if your list mixes Archaea with Bacteria, output will have only one (largest) group
Note: generic intermediate files are collected, so don't try multiple jobs in parallel
END
 return $help;
}

__END__
Place this script in a folder 'trees' along with the PareTree.jar file and make GTDB release subfolders like this:
$ mkdir gtdb202; cd gtdb202
$ perl -pe 's/:[gfocpd]__[^:]+//g' ../../gtdb/bac120_r202.tree | perl -pe "s/\'//g" > d__Bacteria.tree
$ perl -pe 's/:[gfocpd]__[^:]+//g' ../../gtdb/ar122_r202.tree | perl -pe "s/\'//g" > d__Archaea.tree
$ grep -v '^Repre' ../../gtdb/sp_clusters_r202.tsv | perl -pe 's/s__(\S+) /$1__/' | perl -pe 's/;.*//' > labels

