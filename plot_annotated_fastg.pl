#!/usr/bin/perl -w
use strict;

my $file = shift;

my $n = time;

my $Rscript = "
library(\"ggraph\")
library(\"tidygraph\")
graph.data <- read.table(\"$file\", sep=\"\t\")
routes_tidy <- tbl_graph(edges = graph.data, directed = F)
pdf(\"$file.pdf\", width=8, height=6)
ggraph(routes_tidy, layout = \"kk\") + geom_edge_link(aes(width = V3, color=as.factor(V9))) + theme_graph(base_family = 'Helvetica')
dev.off()
";
open(OUT, "> $file.R");
print OUT $Rscript;
close(OUT);
system("Rscript $file.R");
