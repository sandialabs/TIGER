#!/usr/bin/perl
use strict; use warnings;

# OPTIONS
use Cwd 'abs_path';
my $dir = abs_path($0); $dir =~ s/\/([^\/]+)$//; my $scriptname = $1;
my ($verbose, $nickname, $complete, $force, $virus, $taxonomy) = ('', '', '', '', '', '');
my ($outfolder, $prefix);
my $invocation = "Called \"$0 @ARGV\" on " . localtime . "\n";
my $cpu = 1;
my $tax = 'B';
my $gencode = 11;
my $criterion = 'score';
Options(); # see bottom of script; help and other messages, Getopt
print $invocation if $verbose;
my $inDna = abs_path($ARGV[0]);
die "No path for GENOME_FASTA_FILE $ARGV[0]\n" unless $inDna;
if ($outfolder) {$prefix = $inDna; $prefix =~ s/.*\///}
else {$outfolder = $inDna; $outfolder =~ s/([^\/]+)$//; $prefix = $1}
$prefix =~ s/\.f(ast|n)*a(\.gz)*$//;
#die "$outfolder $prefix\n";

# PROGRAM
chdir $outfolder;
my (%stats); 
ReadTax();
#die "$nickname $taxonomy\n";
$nickname = "-nickname " . $nickname if $nickname;
my ($extra) = (' '); $extra = $nickname if $nickname; 
RunCommand("perl $dir/tater.pl -tax $tax -gencode $gencode $extra $verbose $force $prefix", "$prefix.gff"); # Prepare gff with tRNAs, integrases
if ($complete) {$extra = "-complete "} else {$extra = ' '}
if ($virus) {$extra .= "-virus $virus "}
RunCommand("perl $dir/dnaStats.pl $extra$force $verbose $prefix", "$prefix.stats"); # take stats and prepare DNA subfolders
for (@{ReadFile("$prefix.stats")}) {$stats{$1} = $_ if /(\S+)/}
for (keys %stats) {PerDna($_) unless $_ eq 'all'}
PrintFinal();

# SUBROUTINES
sub ReadTax {
 $taxonomy = "division=;phylum=;order=;class=;family=;genus=;species=;org=;taxid=gencode=$gencode;";
 return unless -f "$prefix.tax";
 for (`cat $prefix.tax`) {
  chomp; my ($taxid, $org, $rank, $code, $nick) = split "\t";
  ($nickname, $gencode) = ($nick, $code);
  my ($div, $phy, $ord, $cla, $fam, $gen, $spp) = split ';', $rank;
  $taxonomy = "division=$div;phylum=$phy;order=$ord;class=$cla;family=$fam;genus=$gen;species=$spp;org=$org;taxid=$taxid;gencode=$gencode;";
  $tax = 'B'; $tax = 'A' if $div eq 'Archaea'; $tax = 'M' if $gencode == 4; $tax = 'G' if $gencode == 25;
 }
}

sub RunCommand {
 my ($command, $checkfile) = @_;
 if (not $force and ($checkfile and -e $checkfile)) {print "Skipping command: $command\n" if $verbose; return}
 print "Running command: $command\n" if $verbose;
 my $out = system($command);
 if ($out) {print "Command '$command' failed with error message $out\n"; exit}
}

sub ReadFile {
 my @ret;
 print "Reading file $_[0]\n" if $verbose;
 if ($_[0] =~ /\.gz$/ and -f $_[0]) {open IN, "zcat $_[0] |"} elsif ($_[0] =~ /\.gz$/) {die "No file $_[0]\n"} else {open IN, $_[0] or die "Can't open $_[0]\n"};
 while (<IN>) {last if /^##FASTA/; next if /^#/; push @ret, $_}
 close IN; chomp @ret; return \@ret;
}

sub PerDna { # Find islands
 my ($entry) = @_;
 chdir($entry); print "Changing to $entry directory\n" if $verbose;
 RunCommand("makeblastdb -in in.fa -out in -dbtype nucl > /dev/null", "in.nsq");
 RunCommand("blastn -query in.trna.fa -db in -outfmt 6 -dust no -task  blastn > in.trnablast", "in.trnablast");
 RunCommand("perl $dir/bin/island_finder.pl in $prefix $nickname $verbose -criterion $criterion > in.islander.log", "in.islander.log");
 chdir('..');
}

sub PrintFinal {
 my %serials;
 open OUT, ">islander.gff"; 
 for (glob "*/in.island.gff") {
  open IN, "$_";
  while (<IN>) {
   chomp;
   next unless /;name=([^;]+)/;
   my $name = $1;
   $serials{$name} ++;
   $name .= $serials{$name} if $serials{$name} > 1;
   s/ID=[^;]+/ID=$name/;
   s/;name=([^;]+)//;
   print OUT $_ . $taxonomy . "\n";
  }
  close IN;
 }
 close OUT;
}

sub log10 {return log(shift)/log(10);}

sub Options {
my $version = '2.1 (Nov 2016)';
#   '    |    '    |    '    |    '    |    '    |    '    |    '    |    '    |
my $help = <<END;
$scriptname version $version
Usage: perl $scriptname [options] GENOME_FASTA_FILE
  * NOTE: gzipped GENOME_FASTA_FILE accepted, autodetecting .gz suffix
  * NOTE: -tax, -gencode, -nickname unnecessary if proper .tax file available, 
     with same prefix in same directory as GENOME_FASTA_FILE
  -outDir:   Output directory. Default: same directory as GENOME_FASTA_FILE.
  -tax:      Taxonomic info for query genome. Enter file in outDir containing 
              NCBI taxonomy string, or use B for Bacteria, A for Archaea, M for 
              Mycoplasmatales/Entomoplasmatales, G for Gracilibacteria/candidate
              division SR1. Automatically sets -gencode. Default: $tax.
  -gencode:  Genetic code table (see NCBI). Default: $gencode.
  -nickname: Brief name for genome (as might be used to start a locus_tag).
  -criterion: Basis for overlap resolution, 3 options: random, score (7-test
              false positive formula), deltaGC. Default = $criterion.
  -virus:    Comma-separated list of entries assigned as viruses.
  -complete: Consider genome complete and categorize replicons. Default:
              consider genome incomplete and call all entries contigs.
  -force:    Overwrite current output files. Default: leave existing files.
  -cpu:      Number of cpus to use. Default: $cpu.
  Additional options: -help, -version, -verbose, -authors, -license

Example: perl $scriptname -verbose GENOME_FASTA_FILE

END
my $authors = "AUTHORS: Kelly Williams (kpwilli\@sandia.gov), Corey Hudson, Britney Lau, Owen Solberg\n";
my $license = "COPYRIGHT AND LICENSE: Copyright (2013) Sandia Corporation. Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains certain rights in this software.
Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS 'AS IS' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.\n";

die $help if @ARGV == 0;
use Getopt::Long;
my $options_okay = GetOptions(
 'help' => sub {print $help; exit},
 'version' => sub {print "$scriptname version $version\n"; exit},
 'authors' => sub {print $authors; exit},
 'license' => sub {print $license; exit},
 'verbose' => sub {$verbose = "-verbose"},
 'outDir=s' => \$outfolder,
 'force' => sub {$force = "-force"},
 'gencode=i' => \$gencode,
 'tax=s' => \$tax,
 'nickname=s' => \$nickname,
 'criterion=s' => \$criterion,
 'virus' => \$virus,
 'complete' => \$complete,
 'cpu=i' => \$cpu,
);
die $help if !$options_okay;
}
