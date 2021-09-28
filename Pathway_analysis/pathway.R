library(fgsea)

run_gsea <- function(cell_deg){
gsc<-gmtPathways('Mouse_GOBP_AllPathways_with_GO_iea_November_01_2020_symbol.gmt')
cell_deg<-cell_deg[order(cell_deg$avg_log2FC, decreasing = F),]
ranks<-cell_deg$avg_log2FC
names(ranks)<-rownames(cell_deg)
cell_gsea <- fgsea(pathways = gsc,
                 stats = ranks,
                 minSize=5,
                 maxSize=500)
return(cell_gsea)
}

run_topGO <- function(geneset,bgset, topgo){
library(org.Mm.eg.db)
library(KEGGREST)
library(topGO)
library(grid)
library(pheatmap)
library(reshape2)
sigGenes<-geneset
anno <- AnnotationDbi::select(org.Mm.eg.db, 
                              keys=bgset, 
                              columns=c("ENSEMBL","SYMBOL", "GENENAME",'ENTREZID'),
                              keytype="SYMBOL")
anSig <- as.data.frame(subset(anno, SYMBOL %in% sigGenes))
backG <- setdiff(bgset,  anSig$SYMBOL)
onts = c( "MF", "BP", "CC" )
geneIDs = bgset
inUniverse = geneIDs %in% c(anSig$SYMBOL,  backG) 
inSelection =  geneIDs %in% anSig$SYMBOL 
alg <- factor( as.integer( inSelection[inUniverse] ) )
names(alg) <- geneIDs[inUniverse]
tgd <- new( "topGOdata", ontology='BP', allGenes = alg, nodeSize=5,
            annot=annFUN.org, mapping="org.Mm.eg.db", ID = "SYMBOL" )
resultTopGO.elim <- runTest(tgd, algorithm = "elim", statistic = "fisher" )
resultTopGO.classic <- runTest(tgd, algorithm = "classic", statistic = "fisher" )
tab <- GenTable( tgd, Fisher.elim = resultTopGO.elim, 
                 Fisher.classic = resultTopGO.classic,
                 orderBy = "Fisher.classic" , topNodes = topgo)
return(tab)
     
}
