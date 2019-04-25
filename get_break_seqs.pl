#!/usr/perl -w
use strict;
use BIO::DB::Fasta;

my $in = shift;
my $reffile = shift;
my $db = Bio::DB::Fasta->new($reffile);

open(IN, "< $in");

while (my $line = <IN>) {
  chomp($line);
  my @f = split("\t", $line);
  my ($seq1, $seq2, $seq3) = ('') x 3;
  if ($f[1] ne "") {
    my $edge = 0;
    if ($f[7] < -0.1) {
      $edge = $f[2];
    } elsif ($f[7] > 0.1) {
      $edge = $f[1]
    } else {
      next;
    }
    my $start1 = $edge - 20;
    my $end1 = $edge + 20;
    $seq1 = $db->seq($f[0], $start1 => $end1);
  }
  if ($f[4] ne "") {
    my $start2 = $f[4] - 20;
    my $end2 = $f[4] + 20;
    $seq2 = $db->seq($f[3], $start2 => $end2);
  }
  if ($f[5] ne "") {
    my $start3 = $f[5] - 20;
    my $end3 = $f[5] + 20;
    $seq3 = $db->seq($f[3], $start3 => $end3);
  }
  print join("\t", $f[9], $seq1, $seq2, $seq3), "\n";
}
