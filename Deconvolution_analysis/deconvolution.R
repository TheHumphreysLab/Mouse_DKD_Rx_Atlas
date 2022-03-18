##### Sample preparation ##########
library(BisqueRNA)
library(dplyr)
library(plot1cell)
#### deconvolution mouse bulk ######
dn<-SetIdent(janssen_dn, value = 'group')
dn<-subset(dn, idents=c("Grp3b"))
dn<-SetIdent(dn, value = 'celltype')
sc.eset <- BisqueRNA::SeuratToExpressionSet(dn, delimiter="_", position=1, version="v3")
bulk<-read.csv('all_counts.csv')
rownames(bulk)<-bulk$X
bulk<-bulk[,-1]
metadata<-pbmc.integrated@meta.data
metadata$mouse<-substr(metadata$orig.ident, 2, 5)
metadata<-metadata[!duplicated(metadata$orig.ident),]
metadata<-metadata[order(metadata$mouse),]
metadata<-metadata[!grepl('Rosi', metadata$group2),]
metadata<-metadata[!grepl('A3029', metadata$orig.ident),]
colnames(bulk)<-metadata$orig.ident
bulk2<-as.matrix(bulk)
bulk.eset <- Biobase::ExpressionSet(assayData = bulk.matrix2)
save(bulk.eset, sc.eset, file = 'mousedata_for_bisque.Rda' )
load("~/janss.marker2.Rda")
markers %>%
  group_by(cluster) %>%
  top_n(n = 100, wt = avg_log2FC) -> topgenes
topgenes<-unique(topgenes$gene)
res <- BisqueRNA::ReferenceBasedDecomposition(bulk.eset, sc.eset, markers=topgenes, use.overlap=F)
save(res, file = 'bisque_mousedb_full.Rda' )

######### deconvolution Fan et al ########
bulk<-read.csv('/mnt/sdc/Human_DN_diabetes/fastq/SRR10691631/featureCounts_count.csv')
rownames(bulk)<-bulk$GeneID
bulk<-bulk[,-1]
colnames(bulk)<-c(paste0('AD',1:21),paste0('ED',1:6),paste0('Norm',1:9))
bulk<-convert_geneid(bulk, species = 'human')
rownames(bulk)<-firstup(rownames(bulk))
common_gene<-intersect(rownames(sc.eset), rownames(bulk))
bulk<-bulk[common_gene,]
pheno.matrix_bulk<-data.frame(Sample=colnames(bulk), Group=c(rep('AD',21),rep('ED',6), rep('Norm',9)), row.names = colnames(bulk))
metadata <- data.frame(labelDescription= colnames(pheno.matrix_bulk), row.names=colnames(pheno.matrix_bulk))
bulk.eset = ExpressionSet(assayData = data.matrix(bulk), phenoData =  new("AnnotatedDataFrame", data = pheno.matrix_bulk, varMetadata = metadata) )
save(bulk.eset, file = 'bulk.dn1.Rda')
res <- BisqueRNA::ReferenceBasedDecomposition(bulk.eset, sc.eset, markers=topgenes, use.overlap=F)
save(res, file = 'bisque_humandn1_full.Rda' )

######### deconvolution Levin et al ########
bulk<-read.csv('/mnt/sdc/Human_DN_diabetes/Raw_counts_tub.csv')
bulk<-bulk[!duplicated(bulk$Gene.symbol),]
rownames(bulk)<-bulk$Gene.symbol
bulk<-bulk[,-1]
bulk<-bulk[rowSums(bulk)>4,]
rownames(bulk)<-firstup(rownames(bulk))
common_gene<-intersect(rownames(sc.eset), rownames(bulk))
bulk<-bulk[common_gene,]
pheno.matrix_bulk<-data.frame(Sample=colnames(bulk), Group=c(rep('DN',19),rep('LD',20)), row.names = colnames(bulk))
metadata <- data.frame(labelDescription= colnames(pheno.matrix_bulk), row.names=colnames(pheno.matrix_bulk))
bulk.eset = ExpressionSet(assayData = data.matrix(bulk), phenoData =  new("AnnotatedDataFrame", data = pheno.matrix_bulk, varMetadata = metadata) )
save(bulk.eset, file = 'bulk.dn2.Rda')
res <- BisqueRNA::ReferenceBasedDecomposition(bulk.eset, sc.eset, markers=topgenes, use.overlap=F)
save(res, file = 'bisque_humandn2_full.Rda' )




