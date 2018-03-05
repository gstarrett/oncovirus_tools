library("reshape2")
oncopanel.df <- read.table("", sep="\t", header=T, row.names=1)

 <- apply(oncopanel.df, 1, function(x) apply(oncopanel.df ,1 ,function(y) return(fisher.test(matrix(c(length(which(x > 0 & y > 0)), length(which(x > 0 & y == 0)), length(which(x == 0 & y > 0)), length(which(x == 0 & y == 0))), nrow=2, ncol=2, byrow=TRUE ))$p.value) ))
