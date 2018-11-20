#!/usr/bin/perl -w
use strict;

my $in = shift;
my $virus = shift;
chomp($virus);
my %intHash;
# my %virHash;
my $threshold = 0.05; # minimum fraction of one-end, virus-mapped pairs needed for integration call
my $readLen = 100; # arbitrary value for viral read length
my $totalPairs = 0;
open(SAM, "samtools view $in |"); # read in bam/sam
open(OUT, "> $in.$virus.bedpe");
open(FQ, "> $in.$virus.fastq");

while (<SAM>) {
  next unless $_ =~ /$virus/;
  chomp($_);
  my @f = split("\t", $_);
  if ($f[6] eq $virus) {
    my $alnLen = 0;
    my $strand = "+"; # default is positive strand for read mapping
    if ($f[1] & hex("0x10")) { # find SAM flag for reverse complemented read
      $strand = "-"; # set as reverse strand
    }
    while ($f[5] =~ /(\d+)([MIDNSHPX=])/g) { # functional mapping length relative to reference is calculated by cigar string
      my $alpha = $2;
      my $num = $1;
      if ($alpha =~ /[MDN=X]/) {
          $alnLen += $num;
      }
    }
    push( @{$intHash{$f[2]}}, [$f[3],$f[3]+$alnLen, $f[7]]);
    print OUT join("\t", @f[2,3], $f[3]+$alnLen, ".", ".", $strand, @f[6,7], $f[7]+$readLen, ".", ".", "."), "\n"; # Not actually bedpe format, but structured to work with bedtools
    $totalPairs++; # counting only host read with viral mate pair
  }
  print FQ "@" . $f[0] . "-" . $f[2] . ":" . $f[3] . "\n" . $f[9] . "\n+\n" . $f[10] . "\n"; # write out all reads in pairs containing viral sequence
}
if ($totalPairs == 0) {
  print "No discordant pairs containing $virus found!\n\n";
  exit 0;
}
open(SUM, "> $in.$virus.int.txt");
print SUM join("\t","#hostChrom","hostStart","hostEnd","virChrom","virStart","virEnd","numPairs","hostSkew", "virSkew"),"\n";

# for my $key (sort keys %intHash) {
#   my @qrts =  &quartiles($intHash{$key});
#   my $q1Cnt = &count($intHash{$key}, $qrts[0]);
#   my $q3Cnt = &count($intHash{$key}, $qrts[2]);
#   my @virQrts =  &quartiles($virHash{$key});
#   print SUM join("\t", $in, $key, @qrts, $q1Cnt, $q3Cnt, scalar(@{$intHash{$key}}),@virQrts), "\n";
#   #print SUM join("\t", $in, $key, median($intHash{$key}), stddev($intHash{$key}), @qrts, $q1Cnt, $q3Cnt, scalar(@{$intHash{$key}}), $modal->is_multimodal), "\n";
# }

close(SAM);
close(OUT);
# close(FQ);

open(MERGE, "bedtools merge -s -i $in.$virus.bedpe -c 7,8,9 -o count,min,max |");
while (<MERGE>) {
  chomp($_);
  my @f = split("\t",$_);
  my $strand = splice(@f,3,1);
  if($f[3]/$totalPairs>$threshold){
    # print $_, "\n";
    my $hostMid = ($f[1] + $f[2])/2;
    my $hostLeft = 0;
    my $hostRight = 0;
    my $virMid = ($f[4] + $f[5])/2; # midpoint for skewness calculation
    # Using start site biases skewness towards 5', but simplifies code and doesn't have a major effect on observed results
    my $virLeft = 0;
    my $virRight = 0;
    for my $entry (@{$intHash{$f[0]}}) {
      # print join("\t", @$entry), "\n";
      if ($$entry[0] >= $f[1] && $$entry[0] < $hostMid) { # count host read start sites upstream of midpoint
        $hostLeft++;
      } elsif ($$entry[0] > $hostMid && $$entry[0] <= $f[2]) { # count host read start sites downstream of midpoint
        $hostRight++;
      }
      if ($$entry[2] >= $f[1] && $$entry[2] < $virMid) { # count viral read start sites upstream of midpoint
        $virLeft++;
      } elsif ($$entry[2] >= $virMid && $$entry[2] <= $f[2]) { # count viral read start sites upstream of midpoint
        $virRight++;
      }
    }
    my $virCov = ($virLeft+$virRight); #total reads evaluated for integration event
    my $hostCov = ($hostLeft+$hostRight);
    my $hostSkew = 0; # set default skew value of 0 in case of div/0 errors
    my $virSkew = 0;
    $hostSkew = ($hostLeft-$hostRight)/$hostCov if $hostCov > 0;
    $virSkew = ($virLeft-$virRight)/$virCov if $virCov > 0;
    print SUM join("\t", @f[0..2], $virus, @f[4,5], $f[3], $hostSkew, $virSkew), "\n"; # final output
  }
}
close(SUM);

my $spades_out = $in . "_spades";
my $spades_opt = "-t 4";
system("spades.py -s $in.$virus.fastq -o $spades_out $spades_opt");

# sub quartiles {
#   my $vector = shift;
#   my $n = scalar @$vector;
#   my @sort = sort { $a <=> $b } @$vector;
#   my $q1 = int(.25 * ($n + 1)) - 1;
#   my $q2 = int(.5 * ($n + 1)) - 1;
#   my $q3 = int(.75 * ($n + 1)) - 1;
#   return ($sort[$q1], $sort[$q2], $sort[$q3], $sort[$q3] - $sort[$q1])
# }
#
# sub count {
#   my $vector = shift;
#   my $site = shift;
#   my $start = $site - 500;
#   my $end = $site + 500;
#   my $count = 0;
#   for my $entry (@$vector) {
#     if ($start <= $entry && $end >= $entry) {
#       $count++;
#     }
#   }
#   return $count;
# }
