#!/usr/bin/perl -w
use strict;
use Data::Dumper;

# Input blast data format
# find /Users/starrettgj/Documents/MCC_assemblies/181024 -name "assembly_graph.fastg" | xargs -n1 -P6 -I{} blastn -query {} -out {}.blastn.txt -evalue "1e-10" -db /Volumes/data/REFS/Homo_sapiens/hg19_virus/VirusHg19_combined_Update -outfmt "6 qseqid sseqid qstart qend sstart send length mismatch gapope sstrand evalue bitscore qseq sseq"

my $fastgData = shift;
my $blast = shift;

my %blastAnn;
open(BLAST, "< $blast");
while (my $line = <BLAST>) {
  my @f = split("\t", $line);
  my $seqName = "";
  if ($f[0] =~ /^([^':]+)(:.+)?;/) {
    $seqName = $1;
    my $data = [$seqName, @f[2,3], "$f[1]:$f[4]-$f[5]"];
    push(@{$blastAnn{$seqName}}, $data);
  }
}
#print Dumper(%blastAnn);
close(BLAST);
open(FASTG, "< $fastgData");
while(my $line = <FASTG>){
  chomp($line);
  if ($line =~ /NA/) {
    print $line, "\t", "NA", "\n";
  } else {
    my @f = split("\t", $line);
    my $found = 0;
    if(exists $blastAnn{$f[5]}) {
      for my $entry (@{$blastAnn{$f[5]}}) {
        if (($$entry[1] <= $f[6] && $$entry[2] >= $f[7]) || ($$entry[1] > $f[6] && $$entry[1] < $f[7]) ||  ($$entry[2] > $f[6] && $$entry[2] < $f[7])) {
          print $line, "\t", $$entry[3], "\n";
          $found = 1;
          last;
        }
      }
      if ($found == 0) {
        print $line, "\t", "NA", "\n";
      }
    }
  }
}
