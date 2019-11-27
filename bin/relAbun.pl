#!/usr/bin/perl
use strict; use warnings;
use Getopt::Long;

# MEASURES RELATIVE ABUNDANCES OF KEY MONO- and DI-NUCLEOTIDES (C=G AA=TT AC=GT AG=CT AT CA=TG CC=GG CG GA=TC GC TA)
# reference: Karlin, Mrazek, Campbell JBact179:3899
# INPUT: FastA
# AUTHOR: Kelly Williams, kpwilli@sandia.gov
# for bias between two DNAs, take average absolute difference of relative abundance values for each mono- or di-nucleotide

# MESSAGES
my $version = "$0 1.0 (Nov 2016)\n";
my $usage = "usage: perl $0 [options] -- fastaFile(s)\n";
my $help = "$version\n$usage\n".
 "-L Don't treat DNA as circle\n".
 "-e Entry within fasta file\n".
 "-x Relaxed search for entry name\n".
 "-coord Coordinates in entry (format: start-stop)\n".
 "-o outfile (default=stdout)\n".
 "3 example usages:\n".
 "an island: $0 -L -e NC -x -coord 1446290-1461401 NC_016795.fna\n".
 "a replicon: $0 NC_016795.fna\n".
 "a genome: $0 NC_016795.fna NC_016777.fna\n".
 "Additional options: -version, -usage, -help\n";

# ARGUMENTS
die ("$help\n") unless (@ARGV > 0);
my ($outfile, $entry, $relaxed, $coord, $linear) = ('', 0, 0, 0, 0);
my ($verbose, $stdCommand) = (0, 0);
my $options_okay = GetOptions (
 'o=s' => \$outfile,
 'L' => \$linear,
 'e=s' => \$entry,
 'x' => \$relaxed,
 'coord=s' => \$coord,
 'verbose' => \$verbose,
 'version' => sub { $stdCommand = $version },
 'usage' => sub { $stdCommand = $usage },
 'help' => sub { $stdCommand = $help },
);
die $stdCommand if $stdCommand;
die $help if !$options_okay;
die ("Please list input files after double-dash\n$help") unless @ARGV ;
die "-coord must be in format start-stop\n" if $coord and $coord !~ /^(\d+)-(\d+)$/;

# MAIN PROGRAM
if ($outfile) {open my $OUT, ">$outfile"; select $OUT}
my @monos = qw/A C G T/;
my @dis = qw/AA AC AG AT CA CC CG CT GA GC GG GT TA TC TG TT/;
my @keyDinucs = qw/AA AC AG AT CA CC CG GA GC TA/;
my %allCts; 
for my $infile (@ARGV) {
 die "No infile $infile\n" unless -s $infile;
 my $text = do { local( @ARGV, $/ ) = "$infile"; <> };
 my $dnaCt = 0;
 while ( $text =~ s/>(\S+)[^\n]*\n([^>]+)// ) {
  my ($record, $dna) = ($1, uc $2);
  if ($entry) {
   next if not $relaxed and $entry ne $record;
   next unless $record =~ /$entry/;
  }
  $dna =~ s/[0-9\s]//g;
  $dnaCt ++;
  if ($coord and $coord =~ /^(\d+)-(\d+)$/) {
   my $dnalen = length $dna;
   my ($L, $len) = ($1-1, 1+abs($1-$2));
   $L = $2-1 if $1 > $2;
   die "dna $entry shorter than $L + $len\n" if $dnalen < $L + $len;
   $dna = substr $dna, $L, $len;
   Write("$record:$coord", RAs(CountNuc($dna))); exit;
  }
  Write($record, RAs(CountNuc($dna)));
 }
 Write('all', RAs(\%allCts)); # if $dnaCt > 1;
}

# SUBROUTINES
sub Write {my %ra = %{$_[1]}; print "$_[0]\t$ra{C}"; for (@keyDinucs) {print "\t$ra{$_}"} print "\n"}

sub CountNuc {
 my $dna = $_[0];
 my (%cts); for (@monos, @dis) {$cts{$_} = 0}
 unless ($linear) {$dna = $dna . substr $dna, 0, 1} # permute circle to complete last dinuc
 for ($dna =~ /(?=(.{2}))/g) { $cts{$_} ++ } # DINUC COUNTS: the slowest step
 for my $di (@dis) {
  $cts{tot} += 2 * $cts{$di};
  $di =~ /^(.)/;
  $cts{$1} += $cts{$di} ; # MONONUC COUNTS
 }
 for (qw/A C/) {$cts{$_} = $cts{Revcom($_)} += $cts{$_}}
 for (@keyDinucs) {$cts{$_} += $cts{Revcom($_)}}
 for (@monos, @keyDinucs, 'tot') {$allCts{$_} += $cts{$_}}
 return \%cts;
}

sub RAs { # FREQUENCIES, RELATIVE ABUNDANCES: factor out lower-order biases
 my %cts = %{$_[0]};
 my (%abuns); 
 for (@monos, @keyDinucs) {$abuns{$_} = $cts{$_}/$cts{tot}} # Frequencies at this stage
 for (@keyDinucs) {
  /^(.)(.)$/;
  $abuns{$_} /= $abuns{$1} * $abuns{$2} if $abuns{$1} and $abuns{$2};
 }
 for (@monos) {$abuns{$_} /= 0.25}
 for (keys %abuns) {$abuns{$_} = sprintf ( "%.4f", $abuns{$_} )}
 return \%abuns;
}

sub Revcom {my $ret = reverse $_[0]; $ret =~ tr/ACGT/TGCA/; return $ret}
