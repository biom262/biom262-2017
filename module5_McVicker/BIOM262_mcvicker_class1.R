
gwas.tab <- read.table("class1/BIOM262/RA_GWAS.txt", header=T)

##
## Make a histogram of p-values
##
hist(gwas.tab$p.val, breaks=50)

## notice that the distribution looks faily uniform, but with an excess of low
## p-values

##
## We can better see this by using a quantile-quantile plot, which is
## used to compare two distributions. Under the null hypothesis,
## p-values should be uniformly distributed. We can use a qq-plot to
## compare the observed p-values to what would be expected under the
## null
##
expect.p <- 1:nrow(gwas.tab) / nrow(gwas.tab)
obs.p <- gwas.tab$p.val
qqplot(-log10(expect.p), -log10(obs.p),
       xlab="observed -log10(p-values)",
       ylab="observed -log10(p-values)")

## add a line with a slope of 1, to show what the p-values would look
## like if the followed the null distribution
abline(a=0, b=1, col="red")

## there are some very low p-values that make it difficult
## to see the relationship between the two distributions,
## so threshold p-values below 1e-20 to 1e-20

obs.p[obs.p < 1e-20] <- 1e-20
qqplot(-log10(expect.p), -log10(obs.p),
       xlab="observed -log10(p-values)",
       ylab="observed -log10(p-values)")
abline(a=0, b=1, col="red")


##
## Make a Manhattan plot!
##

## assign colors to each chromosome
chrom <- unique(gwas.tab$chrom)
chrom.color <- rep(c("slategray3", "slategray4"), (length(chrom)+1)/2)
gwas.tab["color"] <- rep("slategray3", nrow(gwas.tab))
for(i in seq(1,length(chrom))) {
  gwas.tab$color[gwas.tab$chrom == chrom[i]] <- chrom.color[i]
}

plot(gwas.tab$genome.pos.hg19, -log10(gwas.tab$p.val),
     col=gwas.tab$color, xlab="genome position", ylab="-log10(P-value)")

##
## Q: Where is the region with the extremely low p-values?
##
gwas.tab[gwas.tab$p.val < 1e-200,]

##
## A: chromosome 6, 32,300,000-32,600,000
##
## Make a zoomed-in plot of the p-values just from chromosome 6 32-33MB
##

f <- which(gwas.tab$chrom=="6" &
           gwas.tab$chrom.pos.hg19 > 32e6 &
           gwas.tab$chrom.pos.hg19 < 33e6)

plot(gwas.tab$chrom.pos.hg19[f], -log10(gwas.tab$p.val[f]),
     xlab="chr6 position", ylab="-log10(P-value)")

##
## Q: Does anyone know what region this is? Can look in UCSC genome browser
##
## A: MHC Class II region
##

##
## The MHC region has a super low p-value, which obscures other signals.
## Let's make the plot again, this time thresholding p-values to 1e-20

p.val <- gwas.tab$p.val
p.val[p.val < 1e-20] <- 1e-20

plot(gwas.tab$genome.pos.hg19, -log10(p.val),
     col=gwas.tab$color, xlab="genome position", ylab="-log10(P-value)",
     ylim=c(0, 20))


## also draw a line indicating the threshold
## for genome-wide significance (5e-8)
lines(x=c(0, max(gwas.tab$genome.pos.hg19)), y=rep(-log10(5e-8),2), col="red")


