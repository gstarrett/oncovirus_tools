#!/usr/bin/perl -w
use strict;
use Statistics::Basic qw(:all);

my $in = shift;
my $virus = shift;
chomp($virus);
my %intHash;

open(SAM, "samtools view $in |");
open(OUT, "> $in.$virus.txt");
open(FQ, "> $in.$virus.fastq");
open(SUM, "> $in.$virus.summary.txt");
while (<SAM>) {
  next unless $_ =~ /$virus/;
  chomp($_);
  my @f = split("\t", $_);
  if ($f[2] eq $virus) {
    if ($f[6] ne "=") {
      push(@{$intHash{$f[6]}}, $f[7]);
    }
  } elsif ($f[6] eq $virus) {
    print OUT join("\t", $in, @f[2,3]), "\n";
  }
  print FQ "@" . $f[0] . "-" . $f[2] . ":" . $f[3] . "\n" . $f[9] . "\n+\n" . $f[10] . "\n";
}

for my $key (keys %intHash) {
  my $modal = mode($intHash{$key});
  print SUM join("\t", $in, $key, median($intHash{$key}), stddev($intHash{$key}), scalar(@{$intHash{$key}}), $modal->is_multimodal), "\n";
}

close(SAM);
close(OUT);
close(FQ);
close(SUM);
