#!/usr/bin/env Rscript
library(edgeR)
data_dir<-c("exon_dir", "intron_dir")
count.matrix<-list()
length.matrix<-list()
for (i in 1:length(data_dir)) {
print(paste0("I am saving count matrices to ",data_dir[i]))
fns <- system2('ls', args = paste0("featureCounts_dir/",data_dir[i],"/*.counts.txt"), stdout = T)
counts.df <- NULL
lengths.df <- NULL
for (fn in fns) {
df <- read.table(fn, check.names=F, sep='\t', header=T)
sname<-strsplit(fn[1], '\\.')[[1]][1]
colnames(df)[3] <- strsplit(sname, '/')[[1]][3]
if (is.null(counts.df)) {
lengths.df <- subset(df, select=c(GeneID, Length))
df <- subset(df, select=-c(Length))
counts.df <- df
} else {
df <- subset(df, select=-c(Length))
counts.df <- merge(counts.df, df, by="GeneID", all.x=T, sort=F)
    }
  }
featurename<-strsplit(data_dir[i],'_')[[1]][1]
write.csv(counts.df, file=paste0("featureCounts_dir/",data_dir[i],'/featureCounts_raw.',featurename,'.csv'), row.names=F)
write.csv(lengths.df, file=paste0("featureCounts_dir/",data_dir[i],'/featureCounts_gene_lengths.',featurename,'.csv'), row.names=F)

cat('\nCount matrix file written to "featureCounts_raw.csv"')

genes <- as.vector(counts.df[[1]])
gene.lengths <- lengths.df[[2]]

# Making a matrix of counts
count.mat <- as.matrix(counts.df[2:length(counts.df)])

## Computing CPM
d <- DGEList(counts=count.mat)
cpm.mat <- cpm(d)
rownames(cpm.mat) <- genes
write.csv(cpm.mat, file=paste0("featureCounts_dir/",data_dir[i],'/featureCounts_cpm.',featurename,'.csv'), row.names=T)
cat('\nCPM matrix file written to "featureCounts_cpm.csv"')

## Computing RPKM
d$genes$Length <- gene.lengths
rpkm.mat <- rpkm(d, gene.lengths)
rownames(rpkm.mat) <- genes
write.csv(rpkm.mat, file=paste0("featureCounts_dir/",data_dir[i],'/featureCounts_rpkm.',featurename,'.csv'), row.names=T)
cat('\nRPKM matrix written to "featureCounts_rpkm.csv"\n')
}
