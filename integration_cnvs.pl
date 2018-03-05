#!/usr/bin/perl -w
use strict;

my $index = shift;
open(IDX, "< $index");
while (my $line = <IDX>) {
  my @f = split("\t", $line);
  # get integration sites
  my $bam = $f[0];
  my $chr = $f[1];
  my $pos = $f[2];
  my $start = $pos - 1000000;
  my $end = $pos + 1000000;
  my $region = "$chr:$start-$end";
  my $out = "$bam.$region.bed";
  system("samtools view -b $bam $region | bedtools genomecov -bg -ibam stdin -g /Users/starrettgj/Documents/REFs/grch37_virus/DFCI_GRCh37_virus.fa > $out");
}
