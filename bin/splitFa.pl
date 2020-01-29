for (`cat isles.fa`) {
 if (/>([^\s\|\+]+)/) {#close OUT; 
  my $nick = $1;
 system ("mkdir Isles/$nick"); 
 open OUT, ">Isles/$nick/$nick.fa"}
 print OUT $_;
}
