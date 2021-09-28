##### Sample preparation ##########

library(MuSiC)
library(dplyr)
load('dn.marker2.Rda')
top50 <- markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC)
celltypes<-levels(dn.integrated)
markers.list<-list()
for (i in 1:length(levels(dn.integrated))){
  cell1<-top50[top50$cluster==celltypes[i],]
  cell1<-cell1$gene
  markers.list[[i]]<-cell1
}
names(markers.list)<-celltypes

bulk<-read.csv('Raw_counts_tub.csv')
bulk<-bulk[!duplicated(bulk$Gene.symbol),]
rownames(bulk)<-bulk$Gene.symbol
bulk<-bulk[,-1]
bulk<-bulk[rowSums(bulk)>4,]

firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}

names(markers.list)<-c('PCT','PST','Inj.PT')
save(markers.list, file = 'markers.list.Rda')

rownames(bulk)<-firstup(tolower(rownames(bulk)))
all_pt<-SetIdent(dn.integrated, value = 'group')
all_pt<-subset(all_pt, idents=c('Grp3b'))
all_pt<-SetIdent(all_pt, value = 'celltype')
#all_pt<-subset(all_pt, idents=c('PCT','PST','Inj.PT'))
pheno.matrix<-all_pt@meta.data
pheno.matrix$cell<-rownames(pheno.matrix)
gene_exprs.matrix<-all_pt@assays$RNA@counts
gene_exprs.matrix<-gene_exprs.matrix[rowSums(gene_exprs.matrix)>0,]
gene_exprs.matrix<-as.matrix(gene_exprs.matrix)
metadata <- data.frame(labelDescription= colnames(pheno.matrix), row.names=colnames(pheno.matrix))
SC.eset = ExpressionSet(assayData = gene_exprs.matrix, phenoData =  new("AnnotatedDataFrame", data = pheno.matrix, varMetadata = metadata) )
save(SC.eset, file = 'sc.eset.Rda')

pheno.matrix_bulk<-data.frame(Sample=colnames(bulk), Group=c(rep('DN',19),rep('LD',20)), row.names = colnames(bulk))
metadata <- data.frame(labelDescription= colnames(pheno.matrix_bulk), row.names=colnames(pheno.matrix_bulk))
bulk.eset = ExpressionSet(assayData = data.matrix(bulk), phenoData =  new("AnnotatedDataFrame", data = pheno.matrix_bulk, varMetadata = metadata) )
save(bulk.eset, file = 'bulk.dn.Rda')

save(markers.list, file = 'markers.list.Rda')

