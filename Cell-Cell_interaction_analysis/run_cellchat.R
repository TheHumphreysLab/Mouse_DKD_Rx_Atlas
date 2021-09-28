run_cellchat<-function(group, n_core, signal_type){
library(CellChat)
print(paste("Group:", group, ";", "n_core:", n_core, ";", "Signal_type:",signal_type))
load('dn.integrated.Rda')
dn.integrated@meta.data$lineage<-as.character(dn.integrated@meta.data$celltype)
dn.integrated@meta.data$lineage[dn.integrated@meta.data$lineage %in% c("PEC","PCT" ,   "PST"  ,  "Inj.PT", "DTL" ,   "TAL" ,   "MD","DCT")]<-"Nephron"
dn.integrated@meta.data$lineage[dn.integrated@meta.data$lineage %in% c("CNT","PC"   ,  "ICA"   , "ICB" )]<-"Ureteric"
dn.integrated@meta.data$lineage[dn.integrated@meta.data$lineage %in% c("EC" )]<-"Vascular"
dn.integrated@meta.data$lineage[dn.integrated@meta.data$lineage %in% c("Fib","JGA" )]<-"Interstitial"
dn.integrated@meta.data$lineage[dn.integrated@meta.data$lineage %in% c("Immune" )]<-"Immune"
dn.integrated@meta.data$lineage[dn.integrated@meta.data$lineage %in% c("Podo" )]<-"Glomerular"
dn.integrated<-SetIdent(dn.integrated, value = "group")
ctrl<-subset(dn.integrated, idents=group)
ctrl<-ctrl@assays$RNA@data
ctrl<-ctrl[rowSums(ctrl)>0,]
meta_data<-dn.integrated@meta.data
ctrl_meta<-meta_data[meta_data$group==group,]
ctrl_meta = data.frame(labels = ctrl_meta$lineage, row.names = colnames(ctrl))
cellchat <- createCellChat(object = ctrl, meta = ctrl_meta, group.by = "labels")
cellchat <- addMeta(cellchat, meta = ctrl_meta)
cellchat <- setIdent(cellchat, ident.use = "labels")
groupSize <- as.numeric(table(cellchat@idents))
CellChatDB <- CellChatDB.mouse 
showDatabaseCategory(CellChatDB)
glimpse(CellChatDB$interaction)
CellChatDB.use <- subsetDB(CellChatDB, search = signal_type) 
cellchat@DB <- CellChatDB.use
cellchat <- subsetData(cellchat)
plan("multiprocess", workers = n_core)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.mouse)
cellchat <- computeCommunProb(cellchat,type = "truncatedMean",trim = 0.1,population.size = TRUE)
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
save(cellchat, file = paste0(group,'_cellchat.Rda'))
}
