#! /usr/bin/perl -w 
use strict ;
use Getopt::Long;

my $version = "$0 1.1 (Sep 2009)" ;
my $usage = "usage: perl $0 -i infile [options]";
my $help = "$version\n$usage\n".
	"-i <file> (Multi-) FastA DNA sequence file: path/filename\n".
	"-e        Entry (record) in file (default=first record)\n".
	"-x        Relaxed search for 'entry' within header line\n".
	"-L        Left coordinate\n".
	"-R        Right coordinate\n".
	"-d        Direction (F or R; default=F)\n".
	"-c        Coordinates in 'start-stop' format (overrides -L, -R and -d)\n".
	"-f        Flank (length of flanks added in lowercase)\n" .
	"-fL       Left Flank (override above)\n" .
	"-fR       Right Flank (override above)\n" .
	"-s        suppress header line (sequence only)\n" .
	"-h        one-word header\n" .
	"-ULU      Convert flankSEQflank to FLANKseqFLANK\n" .
	"-n        Name\n" .
	"Additional options: -version, -usage, -help";

# ARGUMENTS
my ( $infile, $record , $flank , $flankL , $flankR , $L , $R , $coord, $orient , $ulu , $name , $relax, $nohead, $singlehead);
die ("\n$help\n\n") unless (@ARGV > 0);
my $verbose = 0 ;
my $stdCommand = 0 ;
my $options_okay = GetOptions (
	'i=s' => \$infile,
	'e=s' => \$record,
	'x' => \$relax,
	's' => \$nohead,
	'h' => \$singlehead,
	'L=i' => \$L,
	'R=i' => \$R,
	'd=s' => \$orient,
	'c=s' => \$coord,
	'n=s' => \$name,
	'f=i' => \$flank ,
	'fL=i' => \$flankL ,
	'fR=i' => \$flankR ,
	'ULU' => \$ulu ,
	'verbose' => \$verbose,
	'version' => sub { $stdCommand = $version },
	'usage' => sub { $stdCommand = $usage },
	'help' => sub { $stdCommand = $help },
);
die "\n$stdCommand\n\n" if ($stdCommand);
die "\nUnknown arguments @ARGV\n\n$help\n\n" if (@ARGV);
die ("\n$help\n\n") if !$options_okay;
if ($coord) {
 unless ($coord =~ /^(\d+)-(\d+)$/) {die "Coordinate $coord must be in 'x-y' format where x and y are integers\n"}
 if ($1 <= $2) {($L, $R, $orient) = ($1, $2, 'F')} else {($L, $R, $orient) = ($2, $1, 'R')}
}
if ( $L and $R and $L > $R ) { die "Left coordinate cannot be greater than right\n" }
if ( $flank and not $flankL ) { $flankL = $flank }
if ( $flank and not $flankR ) { $flankR = $flank }

my ( $collect , $header , $seq ) ;
#die "$infile $relax $record\n";
open IN , "$infile" or die "No $infile\n" ;
while ( <IN> ) {
	chomp ;
	if ( /^>(\S+)/ ) {
		if ( $collect ) { last }
		if ( $relax and $record and not /$record/ ) { next }
		if ( not $relax and $record and $1 ne $record ) { next }
		s/^>(\S+)// ;
		$header = $_ ;
		#if ( $name ) { $header = $name }
		$collect = $1 ;
		next ;
	}
	unless ( $collect ) { next }
	$seq .= $_ ;
}
unless ( $seq ) { die "No sequence collected\n" }
$L = 1 unless $L;
$R = length $seq unless $R;
#if ( $orient and $orient =~ /^R/ ) {( $flankL , $flankR ) = ( $flankR , $flankL )}
$seq =~ s/[\s0-9]// ;
unless ( $collect ) { die "No record $record found in $infile\n" }
unless ( $seq ) { die "No sequence collected\n" }
$seq = uc $seq ;
print length ( $seq ) , " characters in sequence\n" if $verbose ;
my $rightFlank = substr $seq , $R , length ( $seq ) , '' ;
print length ( $seq ) , " characters in sequence after right trim\n" if $verbose ;
my $leftFlank = substr $seq , 0 , $L - 1 , '' ;
print length ( $seq ) , " characters in sequence after trim\n" if $verbose ;
if ( $flankR ) {
        if ($flankR > length ($rightFlank)) { $flankR = length ($rightFlank) }
	$seq .= lc substr ( $rightFlank , 0 , $flankR ) ;
}
if ( $flankL ) {
        if ($flankL > length ($leftFlank)) { $flankL = length ($leftFlank) }
	$seq = lc ( substr ( $leftFlank , -1 * $flankL ) ) . $seq ;
}
print length ( $seq ) , " characters in final sequence\n" if $verbose ;
if ( $ulu ) { $seq =~ tr/acgtACGTbdhkmrvyBDHKMRVYnN/ACGTacgtBDHKMRVYbdhkmrvyNn/ } 
my $head = ">$collect:$L\-$R$header";
if ( $orient and $orient =~ /^R/ ) {
	$seq = Revcom ( $seq ) ;
	$head = ">$collect:$R\-$L$header<\n" ;
}
$head =~ s/ .*// if $singlehead;
$head = ">$name" if $name;
unless ($nohead) {print "$head\n"}
print "$seq\n" ;

sub Revcom {
	my $seq = reverse $_[0] ;
	$seq =~ tr/acgtACGTbdhkmrvyBDHKMRVY/tgcaTGCAvhdmkybrVHDMKYBR/ ;
	return $seq ;
}

