#!/usr/bin/perl -w
use strict;
use Bio::Seq;
use Bio::SeqIO;
use Bio::DB::Fasta;
use Algorithm::BloomFilter;
use Text::Levenshtein qw(distance);
use Bio::Perl;
use Data::Dumper;

my $makeblastdb = qx(which makeblastdb);
chomp($makeblastdb);

my $in = shift;
my $virus = shift;
my $ref = shift;
my $k = 32;
my $timestamp = time();
chomp($virus);
my %intHash;
# my %virHash;
my $threshold = 0.05; # minimum fraction of one-end, virus-mapped pairs needed for integration call
my $readLen = 100; # arbitrary value for viral read length
my $totalPairs = 0;

# read in fasta files
print STDERR "Reading in reference fasta\n";
#my @refFiles = split(",", $ref);
my $seqs = Bio::DB::Fasta->new($ref);
# my @ids = $seqs->get_all_primary_ids;
# print Dumper(@ids);
# goto BLAST;

open(SAM, "samtools view $in |"); # read in bam/sam
open(OUT, "> $in.$virus.$timestamp.bedpe");
open(FQ, "> $in.$virus.$timestamp.fastq");

# make virus kmer bloom filter
print STDERR "Finding unique kmers in $virus and making Bloom filter\n";
my %seen;
my $filter = Algorithm::BloomFilter->new(1000000, 129);
#print "$seqs{$virus}\n";
#exit;
my $sequence = $seqs->seq("$virus");
my $length = $seqs->length("$virus");
my $revsequence = revcom($sequence)->seq();
for my $i (0 .. $length - $k) {
  my $kmer = substr($sequence, $i, $k);
  my $revKmer = substr($revsequence, $i, $k);
  $filter->add($kmer) unless $seen{$kmer}++;
  $filter->add($revKmer) unless $seen{$revKmer}++;
  #print "$kmer\t$revKmer\n";
}

#print serialize($filter), "\n";

# process bam
print STDERR "Processing $in for sequences matching $virus kmers\n";
while (<SAM>) {
  chomp($_);
  my @f = split("\t", $_);
  if ($_ =~ /$virus/) {
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
  } else {
    for my $i (0 .. length($f[9]) - $k) {
      my $kmer = substr($f[9], $i, $k);
      if ($filter->test($kmer)) {
        $totalPairs++; # counting because it is likely a read spanning the junction
        print join("\t", @f[2,3], $f[3], ".", ".", "+", @f[6,7], $f[7]+$readLen, ".", ".", "."), "\n"; # Not actually bedpe format, but structured to work with bedtools
        print FQ "@" . $f[0] . "-" . $f[2] . ":" . $f[3] . "\n" . $f[9] . "\n+\n" . $f[10] . "\n"; # write out all reads in pairs containing viral sequence
        last;
      } else {
        next;
      }
    }
  }
}
if ($totalPairs == 0) {
  print "No discordant pairs containing $virus found!\n\n";
  exit 0;
}
open(SUM, "> $in.$virus.$timestamp.int.txt");
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

open(MERGE, "bedtools merge -s -i $in.$virus.$timestamp.bedpe -c 7,8,9 -o count,min,max |");
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

my $spades_out = $in . "_$timestamp" . "_spades";
my $spades_opt = "-t 4";
system("spades.py -s $in.$virus.$timestamp.fastq -o $spades_out $spades_opt");

# BLAST:
# $spades_out = shift;
# write out saved merge fasta and create blat reference
system("$makeblastdb -in $ref -dbtype nucl -out $spades_out/tmpblastdb");
# blat contigs against human virus references
system("blastn -query $spades_out/assembly_graph.fastg -out $spades_out/assembly_graph.fastg.blastn.txt -evalue \"1e-10\" -db $spades_out/tmpblastdb -outfmt \"6 qseqid sseqid qstart qend sstart send length mismatch gapope sstrand evalue bitscore qseq sseq\"");

#BLAST:
#$spades_out = shift;
my $fastg = Bio::DB::Fasta->new("$spades_out/assembly_graph.fastg");
open(BLAST, "< $spades_out/assembly_graph.fastg.blastn.txt");
my %blastHits;
while(my $line = <BLAST>) {
  chomp($line);
  my @f = split("\t", $line);
  #print $f[9], "\n";
  if($f[9]<1e-18){
    my @nodes = split(":", $f[0]);
    if ($line =~ /$virus/) {
      push(@{$blastHits{$f[0]}[0]}, \@f);
    } else {
      push(@{$blastHits{$f[0]}[1]}, \@f);
    }
  }
}
my $window = 10;
open(MHSEQ, "> $spades_out/mh.$window.seq.txt");
open(MHOUT, "> $spades_out/mh.$window.txt");

for my $key (keys %blastHits) {
    for my $virHit (@{$blastHits{$key}[0]}) {
      for my $hostHit (@{$blastHits{$key}[1]}) {
        if ($$virHit[2] <= $$hostHit[3]+30 && $$virHit[2] >= $$hostHit[2]-30) {
          my $overlap = $$hostHit[3] - $$virHit[2];
          my ($virStart, $virEnd);
          if ($$virHit[4] < $$virHit[5]) {
            $virStart = $$virHit[4] - $window;
            $virEnd = $$virHit[4] + $window + $overlap;
          } else {
            $virStart = $$virHit[4] + $window;
            $virEnd = $$virHit[4] - $window - $overlap;
          }

          my ($hostStart,$hostEnd);
          if ($$hostHit[4] < $$hostHit[5]) {
            $hostStart = $$hostHit[5] - $window - $overlap;
            $hostEnd = $$hostHit[5] + $window;
          } else {
            $hostStart = $$hostHit[5] + $window + $overlap;
            $hostEnd = $$hostHit[5] - $window;
          }
          my $contigStart = $$virHit[2] - $window;
          my $contigEnd = $$virHit[2] + $window + $overlap;
          my $contigSeq = $fastg->seq($key,$contigStart,$contigEnd);
          #print("$$hostHit[1]\t$$virHit[1]\n");
          my $hostSeq = $seqs->seq($$hostHit[1],$hostStart,$hostEnd);;
          my $virSeq = $seqs->seq($$virHit[1],$virStart,$virEnd);
          print MHSEQ join("\n",$key,$virSeq,$contigSeq,$hostSeq),"\n";
          my @mh = (0,0,0,0,0,0,0,0);
          my $prev = 1;
          for (my $i = 0; $i<length($hostSeq)-5; $i++) {
            if (substr($hostSeq,$i,2) eq substr($virSeq,$i,2)) {
              $prev++;
            } else {
              if ($prev > 1) {
                if ($prev > 8) {
                  $mh[7]++;
                } else {
                  $mh[$prev-2]++;
                }
                $prev = 1;
              }
            }
          }
          my $Ldist = distance($hostSeq, $virSeq);
          my $dist = (length($hostSeq)-$Ldist)/length($hostSeq);
          print MHOUT join("\t", $key, @$hostHit[1..5,9], @$virHit[1..5,9],@mh, $dist),"\n";
          last;
        } elsif ($$virHit[3] >= $$hostHit[2]-30 && $$virHit[3] <= $$hostHit[3]+30){
          my $overlap = $$virHit[3] - $$hostHit[2];
          my ($virStart, $virEnd);
          if ($$virHit[4] < $$virHit[5]) {
            $virStart = $$virHit[5] - $window - $overlap;
            $virEnd = $$virHit[5] + $window;
          } else {
            $virStart = $$virHit[5] + $window + $overlap;
            $virEnd = $$virHit[5] - $window;
          }

          my ($hostStart,$hostEnd);
          if ($$hostHit[4] < $$hostHit[5]) {
            $hostStart = $$hostHit[4] - $window;
            $hostEnd = $$hostHit[4] + $window + $overlap;
          } else {
            $hostStart = $$hostHit[4] + $window;
            $hostEnd = $$hostHit[4] - $window - $overlap;
          }
          my $contigStart = $$virHit[3] - $window - $overlap;
          my $contigEnd = $$virHit[3] + $window;
          my $contigSeq = $fastg->seq($key,$contigStart,$contigEnd);
          #print("$$hostHit[1]\t$$virHit[1]\n");
          my $hostSeq = $seqs->seq($$hostHit[1],$hostStart,$hostEnd);;
          my $virSeq = $seqs->seq($$virHit[1],$virStart,$virEnd);
          print MHSEQ join("\n",$key,$virSeq,$contigSeq,$hostSeq),"\n";
          my @mh = (0,0,0,0,0,0,0,0);
          my $prev = 1;
          for (my $i = 0; $i<length($hostSeq)-1; $i++) {
            if (substr($hostSeq,$i,2) eq substr($virSeq,$i,2)) {
              $prev++;
            } else {
              if ($prev > 1) {
                if ($prev > 8) {
                  $mh[7]++;
                } else {
                  $mh[$prev-2]++;
                }
                $prev = 1;
              }
            }
          }
          my $Ldist = distance($hostSeq, $virSeq);
          my $dist = (length($hostSeq)-$Ldist)/length($hostSeq);
          print MHOUT join("\t", $key, @$virHit[1..5,9], @$hostHit[1..5,9],@mh, $dist),"\n";
          last;
        }
      }
    }
}

# my %contigs;
# my $fastg = Bio::SeqIO->new(-format => 'fasta', -file => "$spades_out/assembly_graph.fastg");
# while( my $seq = $fastg->next_seq() ) {
#   my $virusJunction;
#   my $hostJunction;
#   my $contigName = $seq->display_id();
#   my @nodes = split(":", $contigName);
#   unless ($nodes[0] =~ /.+'$/) {
#     my @data = split("_", $nodes[0]);
#     unless (exists $contigs{$data[1]}) {
#       my $contigLen = $nodes[3];
#       for (my $i = 1; $i<$contigLen; $i++) {
#         if (exists $blastHits{$data[1]}) {
#           for my $hit (@$blastHits{$data[1]}) {
#               if ($i =< $$hit[3] && $i >= $$hit[2])
#
#           }
#         }
#       }
#     }
#   }
# }

# for my $key (sort keys %contigs) {
#   if (defined $contigs{$key}[0] && defined $contigs{$key}[1]) {
#     for my $hit (@${$contigs{$key}[0]}) {
#       print join("\t", @$hit), "\n";
#     }
#     print join("\t", @${$contigs{$key}[1]}), "\n";
#
#   } elsif (defined $contigs{$key}[0]) {
#     for my $hit (@${$contigs{$key}[0]}) {
#       print join("\t", @$hit), "\n";
#     }
#   } elsif (defined $contigs{$key}[1]) {
#     print join("\t", @${$contigs{$key}[1]}), "\n";
#   }
# }

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
