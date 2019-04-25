#!/usr/bin/perl -w
use strict;
use POSIX qw(ceil);

my $in = shift;
my $window = 50;
my %edges;
open(IN, "< $in");
while (my $line = <IN>) {
  if ($line =~ /^>(.+);/) {
    my @nodes = split(":", $1);
    my $seqName = "";
    if ($nodes[0] =~ /(.+)'$/) {
      $seqName = $1;
    } else {
      $seqName = $nodes[0];
    }
    my @data = split("_", $nodes[0]);
    unless (exists $edges{$data[1] . "_0"}) {
      my $length = ceil($data[3]/$window);
      for (my $i=0;$i<$length;$i++) {
        if ($i==0) {
          my $from = $data[1] . "_0";
          my $to = $data[1] . "_0" . $i;
          if ($i+1>=$length) {
            $to = $data[1] . "_1";
          }
          my $weight = 1;
          if ($data[5] =~ /(.+)'/) {
            $weight = $1;
          } else {
            $weight = $data[5];
          }
          print join("\t", $from, $to, $weight, $data[1], $data[3], "$seqName", 1, ($i+1)*$window), "\n";
          $edges{$from} = 1;
          $edges{$to} = 1;
        } elsif ($i+1>=$length) {
          my $from = $data[1] . "_0". ($i-1);
          my $to = $data[1] . "_1";
          my $weight = 1;
          if ($data[5] =~ /(.+)'/) {
            $weight = $1;
          } else {
            $weight = $data[5];
          }
          print join("\t", $from, $to, $weight, $data[1], $data[3], "$seqName", $i*$window+1, $data[3]), "\n";
          $edges{$to} = 1;
        } else {
          my $from = $data[1] . "_0" . ($i-1);
          my $to = $data[1] . "_0" . $i;
          my $weight = 1;
          if ($data[5] =~ /(.+)'/) {
            $weight = $1;
          } else {
            $weight = $data[5];
          }
          print join("\t", $from, $to, $weight, $data[1], $data[3], "$seqName", $i*$window+1, (($i+1)*$window)), "\n";
          $edges{$to} = 1;
        }
      }
    }
    if (scalar @nodes > 1) {
      my $dir1 = 1;
      if ($nodes[0] =~ /.+'$/) {
        $dir1 = 0;
      }
      my @start = split("_", $nodes[0]);
      my @ends = split(",", $nodes[1]);
      for my $end (@ends) {
        my $dir2 = 0;
        if ($end =~ /.+'$/) {
          $dir2 = 1;
        }
        my @f = split("_", $end);
        my $from = $start[1] . "_" . $dir1;
        my $to = $f[1] . "_" . $dir2;
        print join("\t", $from, $to, 0.5, "NA", 1, "", "", ""), "\n";
      }
    }
  }
}
#print "\n";
#for my $key (sort keys %edges) {
#  print $key . "\n";
#}
