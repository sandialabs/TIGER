#! /usr/bin/perl -w
use strict ;
use Getopt::Long ;

# ARGUMENTS
my ( $query , $flank , $supp , $case ) = ( '' , 'r' , 0 , 'u' );
my $version = "$0 1.1 (Jan 2013)" ;
my $usage = "usage: perl $0 [options] -q query seqfiles\n" ;
my $help = "$version\n$usage\nextends perfectly-matched query in both orientations from set of read sequences\n" .
   "test sequence files should have single line for each sequence as in fastq\n" .
   "note that query must be short and fit entirely within read\n" .
   "implements tests for palindromic query and for query duplication in a read\n" .
   "-q     query sequence\n" .
   "-f     flank to extend: l=left, r=right (default=$flank)\n" .
   "-s     suppress extensions below this length (default=$supp)\n" .
   "-c     case of read sequences: u=upper-case, l=lower-case (default=$case)\n" .
   "Additional options: -version, -usage, -verbose, -help" ;

my ( $verbose , $stdCommand ) ;
die "\n$help\n\n" unless ( @ARGV > 0 ) ;
my $options_okay = GetOptions (
 'q=s' => \$query ,
 'f=s' => \$flank ,
 's=i' => \$supp ,
 'c=s' => \$case ,
 'verbose' => \$verbose ,
 'version' => sub { $stdCommand = $version } ,
 'usage' => sub { $stdCommand = $usage } ,
 'help' => sub { $stdCommand = $help } ,
);
die "\n$stdCommand\n\n" if ( $stdCommand ) ;
die "\n$help\n\n" if !$options_okay ;
die "\n$help\n\n" if !@ARGV || !$query || $query =~ /[^ACGTRYMKVBHD]/i || $supp =~ /[^0-9]/ || $case !~ /^[ul]$/ || $f
lank !~ /^[lr]$/ ;

my @files = ( @ARGV ) ;
for ( @files ) { die "No file $_\n" unless -f $_ }
$query = uc $query ; $query = lc $query if $case eq 'l' ; 
my $rquery = Revcomp ( $query ) ; 
if ( $query eq $rquery ) { warn "WARNING: palindromic query\n" }

my ( @hits1 , @hits2 , %ct ) ;
my $longest = 0 ;
for ( @files ) { # GET FORWARD QUERY
 if ($_ =~ /\.gz$/) {push @hits1 , `zcat $_ | grep $query`}
 else {push @hits1 , `grep $query $_` }
}
for ( @hits2 ) { 
 chomp ;
 my $dupetest = $_ ;
 if ( $dupetest =~ s/$query//g > 1 ) { warn "query $query duplicated in read $_\n" }
 if ( $flank eq 'r' ) { s/$rquery.*// } 
 else { 
  s/.*$rquery// ; 
  if ( length $_ > $longest ) { $longest = length $_ }
 }
 $_ = Revcomp ( $_ ) if $flank eq 'r' ;
}
#push @hits1 , ( @hits2 ) ;

$ct{hit} = 0;
for ( sort (@hits1, @hits2) ) { 
 my $len = length $_ ;
 if ( $len < $supp ) { $ct{supp} ++ ; next }
 $ct{hit} ++;
 if ( $flank eq 'l' ) { $_ = ' ' x ( $longest - $len ) . Revcomp ( $_ ) }
 print "$_\n" ; 
}
$ct{supp} = 0 unless $ct{supp} ;
my $sum = $ct{supp} + $ct{hit};
print "$sum hits ($ct{hit} after omitting any < $supp nt) to $flank of $query in @ARGV\n";

sub Revcomp { # Compute the reverse complement of DNA sequence, even if degenerate
 my $dna = $_[0];
 my $revcomp = reverse $dna;
 $revcomp =~ tr/ACGTRYMKVBHDacgtrymkvbhd/TGCAYRKMBVDHtgcayrkmbvdh/; #N,S,W unchanged
 return $revcomp;
}