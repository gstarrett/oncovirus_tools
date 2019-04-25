#!/usr/bin/perl -w
use strict;
use Bio::DB::Fasta;

my @chrs = (1..22,"X","Y");
my $ref = shift;
my $window = 1;
my $db = Bio::DB::Fasta->new($ref);

for (my $i=0;$i<1000;$i++) {
  my $chr = $chrs[int(rand($#chrs))];
  my $chrLen = $db->length($chr)-$window;
  my $chrStart = int(rand($chrLen));
  my $chrEnd = $chrStart+$window;
  print join("\t",$chr,$chrStart,$chrEnd), "\n";
}
