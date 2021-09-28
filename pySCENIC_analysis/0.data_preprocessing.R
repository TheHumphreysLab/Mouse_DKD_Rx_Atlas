load("janssen_seurat.Rda")
injPT<-subset(dn.integrated, idents="Inj.PT")
injPT_count<-injPT@assays$RNA@counts
injPT_count<-injPT_count[rowSums(injPT_count)>4,]
injPT<-subset(injPT, features=rownames(injPT_count))
SaveH5Seurat(injPT, filename = ".h5Seurat")
Convert("injPT.h5Seurat", dest = "h5ad")
