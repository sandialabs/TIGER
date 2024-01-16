#!/usr/bin/perl
use strict; use warnings;

# OPTIONS
use File::Spec;
my $dir = File::Spec->rel2abs($0); $dir =~ s/\/([^\/]+)$//; my $scriptname = $1;
my ($verbose, $faa, $gff, $intDoms, $integron, $pfam, $xer, $force);
my $extraDoms = '';
my $cpu = 0;
my $intRatio = 0.85;
Options(); # See bottom of script; help and other messages, Getopt

# PROGRAM
my (%ints, %intAssoc, @order);
my $base = $faa; $base =~ s/\.faa$//g;
my $doms = join (',', "$base.domtbl", (split ',', $extraDoms));
RunCommand("$dir/hmmsearch --domtblout $base.domtbl --cpu 0 --cut_tc $pfam $base.faa > /dev/null", "$base.domtbl");
RunCommand("perl $dir/dom_to_gff.pl -gffin $gff -dom $doms -protgff $base.pfam.gff -domgff $base.domains.gff", "$base.pfam.gff");
my $gffct = LoadTRs("$base.pfam.gff");
print "$gffct gffs\n",  scalar(keys %ints), " tyrosine recombinase family members\n";
my $recombinases = Read_faa();
my $xerct = Remove_xers($recombinases);
print scalar(keys %ints), " integrases after removal of $xerct xers\n";
my $integrons = IntegronsIn();
for (@{Read_file($intDoms)}) {next unless /^\S+\s+(\S+)/; $intAssoc{$1} = 1} # other domains associated with Ints
print scalar(keys %intAssoc), " reference Int-associated domains\n";
my ($phageCt, $integronCt, $topCt) = Call();
print "$phageCt phage integrases after removal of $integronCt integrons and $topCt with better non-int Pfam calls\n";

# SUBROUTINES
sub RunCommand {
 my ($command, $checkfile) = @_;
 if (not $force and ($checkfile and -s $checkfile)) {print "Skipping command: $command\n" if $verbose; return}
 print "Running command: $command\n" if $verbose;
 my $out = system($command);
 die "Command '$command' failed with error message $out\n" if $out;
}

sub LoadTRs {
 my %TRs;
 my $gffct = 0;
 for (@{Read_file($_[0])}) {
  $gffct ++;
  next unless /ID=([^;]+);pfam1=([^;]+);score1=([^;]+)/;
  my ($id, $top, $topscore) = ($1, $2, $3);
  next unless /=Phage_integrase;score\d+=([^;]+)/;
  %{$ints{$id}} = (phage => $1, top => $top, topscore => $topscore, line => $_ . "\n", integron => 0);
  push @order, $id;
 }
 return $gffct;
}

sub Read_file {
 open IN, $_[0] or die "Can't open $_[0]\n"; my @ret = <IN>; close IN;
 chomp @ret; return \@ret;
}

sub Read_faa {
 my (%prots, $prot);
 for (@{Read_file($faa)}) {
  if (/^>(\S+)/){
   $prot = $1;
   $prot = '' unless $ints{$prot};
  } elsif ($prot) {$prots{$prot} .= $_}
 }
 return \%prots;
}

sub Remove_xers {
 my %fasta = %{$_[0]};
 my $xerct = 0;
 open OUT, ">$base.xer";
 for my $prot (keys %fasta) {
  my $run = "echo '>$prot\n$fasta{$prot}' | $dir/pfscan -f - $xer";
  next unless `$run`; # To Do: collect score
  print OUT "$ints{$prot}{line}"; $xerct++; delete $ints{$prot};
 }
 close OUT;
 return $xerct;
}

sub IntegronsIn {
 RunCommand("$dir/hmmsearch --acc --tblout $base.integron.tbl --cpu 0 --cut_tc $integron $faa > /dev/null", "$base.integron.tbl");
 for (@{Read_file("$base.integron.tbl")}) {
  next if /^#/;
  my @f = split /\s+/;
  my $id = $f[0];
  $id =~ s/.*id\|([0-9]*).*/$1/; # Required?
  next unless $ints{$id};
  my $score = sprintf('%G', $f[4]);
  #next if $score > sprintf('%G', '1.2e-24'); # For famint9
  next if $score > sprintf('%G', '1e-3'); # For intI_Cterm
  #$ints{$id}{integron} = $score;
  $ints{$id}{integron} = $f[8];
 }
}

sub Call {
 open PHAGE, ">$base.phage";
 open INTEGRON, ">$base.integrons";
 my ($phageCt, $integronCt, $topCt) = (0, 0, 0);
 for my $id (@order) {
  next unless $ints{$id};
  # print "$id $ints{$id}{topscore} $ints{$id}{phage} $ints{$id}{integron}\n";
  if (not $intAssoc{$ints{$id}{top}} and $ints{$id}{topscore} > $ints{$id}{phage} and $ints{$id}{topscore} > $ints{$id}{integron}) {
   if ($ints{$id}{top}) {
    $topCt ++; next;
   }
  }
  unless ($ints{$id}{integron}) {print PHAGE $ints{$id}{line}; $phageCt ++; next}
  #if ($ints{$id}{phage}/$ints{$id}{integron} < $intRatio) {print INTEGRON $ints{$id}{line}; $integronCt ++} # show integron score # For famint9
  if ($ints{$id}{integron}) {print INTEGRON $ints{$id}{line}; $integronCt ++} # show integron score # For intI_Cterm
  else {print PHAGE $ints{$id}{line}; $phageCt ++}
 }
 return ($phageCt, $integronCt, $topCt);
}

sub Options {
 my $version = '0.3 (Nov 2016)';
#   '    |    '    |    '    |    '    |    '    |    '    |    '    |    '    |
 my $help = <<END;
$scriptname version $version
Usage: perl $scriptname [options] -faa <file> -gff <file> -pfamDb <file> \
      -intDoms <file> -integron <file> -xer <file>
 Options:
 -faa:       Protein sequence .faa file.
 -gff:       Prokka annotation .gff file.
 -pfam:      Pfam hmm database file.
 -extraDoms: Additional .domtbl files, comma-separated
 -intDoms:   File listing Pfams associated with phage integrases. 
 -integron:  Hmm file for integron integrases.
 -xer:       Profile for Xer proteins.
 -force:     Overwrite current output files. Default: leave existing files.
 -cpu:       Number of cpus to use. Default: $cpu.
 Additional options: -help, -version, -verbose, -authors, -license

Example: perl $scriptname -faa genome.faa -gff genome.gff -integron DB/famint9.hmm -pfamDb DB/Pfam-A.hmm -xer DB/xers.prf -intDoms DB/integrase_domain_pfams.txt -verbose

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
  'faa=s' => \$faa,
  'gff=s' => \$gff,
  'cpu=i' => \$cpu,
  'pfam=s' => \$pfam,
  'extraDoms=s' => \$extraDoms,
  'intDoms=s' => \$intDoms,
  'integron=s' => \$integron,
  'xer=s' => \$xer,
  'force' => \$force,
 ) ;
 die $help if !$options_okay;
}
