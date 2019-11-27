#! /usr/bin/perl
use strict; use warnings;
use IPC::Run3 qw/run3/;
use Cwd 'abs_path';

die "Usage: $0 infasta outdirectory\n" unless @ARGV == 2;
my ($db, $infile, $outdir, $dbfiles) = ('tmrna', $ARGV[0], $ARGV[1], $0); for ($infile, $outdir, $dbfiles) {$_ = abs_path($_)}
$dbfiles =~ s/[^\/]*\/([^\/]+)$/db\/cm\/$db.cm/;
# my $dbfiles = "$lib/cm/$db.cm";

chdir $outdir;
open FINAL, ">rfam.gff" or die "Cannot open outfile rfam.gff\n";
mkdir 'rfam'; chdir 'rfam';
system "cmscan -o /dev/null --cpu 0 --tblout tmrna.tbl --oskip --fmt 2 $dbfiles $infile" unless -f 'tmrna.tbl';

my (@out, $serial);
my %cm = qw/RF00023 tmRNA RF01849 alpha_tmRNA RF01850 beta_tmRNA RF01851 cyano_tmRNA/;
my %clen = qw/RF00023 354 RF01849 338 RF01850 331 RF01851 288/;
for (`cat tmrna.tbl`) {
 next if /^#/;
 my @f = split /\s+/;
 ##idx target_name          accession query_name           accession clan_name mdl mdl_from   mdl_to seq_from   seq_to strand trunc pass   gc  bias  score   E-value inc olp anyidx afrct1 afrct2 winidx wfrct1 wfrct2 description_of_target
 #1    tmRNA                RF00023   GG666849.1           -         -          cm        1      354   192875   192950      +    no    1 0.63   0.0   31.5   2.8e-07  !   ^       -      -      -      -      -      - transfer-messenger RNA
 #next if $f[16] < 20;  # CM score cutoff
 $f[12] =~ s/\'/prime/g; $f[12] =~ s/5prime\&3prime/both/;
 my ($tail, $form) = ("model=$cm{$f[2]};modelL=$f[7];modelR=$f[8];modelLen=$clen{$f[2]};trunc=$f[12];form=", 'std;'); $form = 'prm;' if $cm{$f[2]} =~ /_/;
 ($f[9], $f[10]) = ($f[10], $f[9]) if $f[9] > $f[10];
 push @out, [$f[3], 'cmscan', 'tmRNA', @f[9,10,16,11], '.', $tail.$form];
}
@out = sort {$$a[0] cmp $$b[0] || $$a[3] <=> $$b[3] || $$b[4] <=> $$a[4]} @out;
for (@out) {$serial ++; $$_[8] = "ID=rfam.$serial;" . $$_[8]; print FINAL join("\t", @{$_}), "\n"}

