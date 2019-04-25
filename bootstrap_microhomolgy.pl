#!/usr/bin/perl -w
use strict;
use Bio::DB::Fasta;
use Text::Levenshtein qw(distance);

my @chrs = (1..22,"X","Y");
my $ref = shift;
my $virus = shift;
my $window = 20;
my $db = Bio::DB::Fasta->new($ref);

for (my $i=0;$i<1000;$i++) {
  my $virLen = $db->length($virus)-$window;
  my $virStart = int(rand($virLen));
  my $virEnd = $virStart+$window;
  my $virStrand = int(rand(1));
  my $virSeq = $db->seq($virus,$virStart,$virEnd,$virStrand);

  my $chrSeq = "N";
  until ($chrSeq !~ /N/) {
    my $chr = $chrs[int(rand($#chrs))];
    my $chrLen = $db->length($chr)-$window;
    my $chrStart = int(rand($chrLen));
    my $chrEnd = $chrStart+$window;
    my $chrStrand = int(rand(1));
    $chrSeq = $db->seq($chr,$chrStart,$chrEnd,$chrStrand);
  }

  my @mh = (0,0,0,0,0,0,0,0);
  my $prev = 1;
  for (my $i = 0; $i<length($chrSeq)-1; $i++) {
    if (substr($chrSeq,$i,2) eq substr($virSeq,$i,2)) {
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
  my $Ldist = distance($chrSeq, $virSeq);
  my $dist = (length($chrSeq)-$Ldist)/length($chrSeq);
  print join("\t", @mh, $dist),"\n";
}
