#!/usr/bin/perl -w
use strict;
use Statistics::Basic qw(:all);

my $in = shift;
my $virus = shift;
chomp($virus);
my %intHash;
my %virHash;

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
      push(@{$virHash{$f[6]}}, $f[3]);
    }
  } elsif ($f[6] eq $virus) {
    print OUT join("\t", $in, @f[2,3]), "\n";
  }
  print FQ "@" . $f[0] . "-" . $f[2] . ":" . $f[3] . "\n" . $f[9] . "\n+\n" . $f[10] . "\n";
}

print SUM join("\t","#filename","chromosome","1st_quartile_pos","2nd_quartile_pos","3rd_quartile_pos","interquartile_range","1st_quartile_read_count","3rd_quartile_read_count","total_viral_pairs_on_chromosome","virus_1st_quartile_pos","virus_2nd_quartile_pos","virus_3rd_quartile_pos","virus_interquartile_range"),"\n";

for my $key (sort keys %intHash) {
  my @qrts =  &quartiles($intHash{$key});
  my $q1Cnt = &count($intHash{$key}, $qrts[0]);
  my $q3Cnt = &count($intHash{$key}, $qrts[2]);
  my @virQrts =  &quartiles($virHash{$key});
  print SUM join("\t", $in, $key, @qrts, $q1Cnt, $q3Cnt, scalar(@{$intHash{$key}}),@virQrts), "\n";
  #print SUM join("\t", $in, $key, median($intHash{$key}), stddev($intHash{$key}), @qrts, $q1Cnt, $q3Cnt, scalar(@{$intHash{$key}}), $modal->is_multimodal), "\n";
}

close(SAM);
close(OUT);
close(FQ);
close(SUM);

sub quartiles {
  my $vector = shift;
  my $n = scalar @$vector;
  my @sort = sort { $a <=> $b } @$vector;
  my $q1 = int(.25 * ($n + 1)) - 1;
  my $q2 = int(.5 * ($n + 1)) - 1;
  my $q3 = int(.75 * ($n + 1)) - 1;
  return ($sort[$q1], $sort[$q2], $sort[$q3], $sort[$q3] - $sort[$q1])
}

sub count {
  my $vector = shift;
  my $site = shift;
  my $start = $site - 500;
  my $end = $site + 500;
  my $count = 0;
  for my $entry (@$vector) {
    if ($start <= $entry && $end >= $entry) {
      $count++;
    }
  }
  return $count;
}
