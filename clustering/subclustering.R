#### TAL subclustering analysis ######
tal <- subset(dn.integrated, idents="TAL"))
tal.list<-SplitObject(tal, split.by = 'orig.ident')
for (i in names(tal.list)) {
  tal.list[[i]] <- SCTransform(tal.list[[i]], verbose = F)
  print(paste(i, 'done!'))
}

tal.features <- SelectIntegrationFeatures(object.list = tal.list, nfeatures = 1000)

tal.list <- PrepSCTIntegration(object.list = tal.list, anchor.features = tal.features)
reference_dataset <- which(names(tal.list) == "A3020")
tal.anchors <- FindIntegrationAnchors(object.list = tal.list,
                                       anchor.features = tal.features, reference = reference_dataset, normalization.method = "SCT")
Sys.time()
tal.integrated <- IntegrateData(anchorset = tal.anchors)
Sys.time()
tal.integrated<-ScaleData(tal.integrated)
tal.integrated <- RunPCA(object = tal.integrated, npcs = 30)
tal.integrated <- RunUMAP(object = tal.integrated, dims = 1:30)
tal.integrated <- FindNeighbors(tal.integrated, dims = 1:30)
tal.integrated <- FindClusters(tal.integrated, resolution = 0.3)
save(tal.integrated, file = 'tal.integrated.Rda')
Sys.time()

#### EC subclustering analysis ######
ec <- subset(dn.integrated, idents="EC"))
ec.list<-SplitObject(ec, split.by = 'orig.ident')
for (i in names(ec.list)) {
  ec.list[[i]] <- SCTransform(ec.list[[i]], verbose = F)
  print(paste(i, 'done!'))
}

ec.features <- SelectIntegrationFeatures(object.list = ec.list, nfeatures = 1000)

ec.list <- PrepSCTIntegration(object.list = ec.list, anchor.features = ec.features)
reference_dataset <- which(names(ec.list) == "A3020")
ec.anchors <- FindIntegrationAnchors(object.list = ec.list,
                                      anchor.features = ec.features, reference = reference_dataset, normalization.method = "SCT")
Sys.time()
ec.integrated <- IntegrateData(anchorset = ec.anchors)
Sys.time()
ec.integrated<-ScaleData(ec.integrated)
ec.integrated <- RunPCA(object = ec.integrated, npcs = 30)
ec.integrated <- RunUMAP(object = ec.integrated, dims = 1:10)
ec.integrated <- FindNeighbors(ec.integrated, dims = 1:10)
ec.integrated <- FindClusters(ec.integrated, resolution = 0.3)
save(ec.integrated, file = 'ec.integrated.Rda')
Sys.time()

#### Fib subclustering analysis ######
fib <- subset(dn.integrated, idents="Fib"))
fib.list<-SplitObjfibt(fib, split.by = 'orig.ident')
for (i in names(fib.list)) {
  fib.list[[i]] <- SCTransform(fib.list[[i]], verbose = F)
  print(paste(i, 'done!'))
}

fib.features <- SelfibtIntegrationFeatures(objfibt.list = fib.list, nfeatures = 1000)

fib.list <- PrepSCTIntegration(objfibt.list = fib.list, anchor.features = fib.features)
reference_dataset <- which(names(fib.list) == "A3020")
fib.anchors <- FindIntegrationAnchors(objfibt.list = fib.list,
                                     anchor.features = fib.features, reference = reference_dataset, normalization.method = "SCT")
Sys.time()
fib.integrated <- IntegrateData(anchorset = fib.anchors)
Sys.time()
fib.integrated<-ScaleData(fib.integrated)
fib.integrated <- RunPCA(objfibt = fib.integrated, npcs = 30)
fib.integrated <- RunUMAP(objfibt = fib.integrated, dims = 1:10)
fib.integrated <- FindNeighbors(fib.integrated, dims = 1:10)
fib.integrated <- FindClusters(fib.integrated, resolution = 0.3)
save(fib.integrated, file = 'fib.integrated.Rda')
Sys.time()

#### Immune subclustering analysis ######
immune <- subset(dn.integrated, idents="Immune"))
### Immune cells subclustering was directly performed on the intial integration data due to the fact that the cells are sparse in some samples
### which prevents a new integration analysis https://github.com/satijalab/seurat/issues/1209
immune.integrated <- RunPCA(objfibt = immune.integrated, npcs = 30)
immune.integrated <- RunUMAP(objfibt = immune.integrated, dims = 1:10)
immune.integrated <- FindNeighbors(immune.integrated, dims = 1:10)
immune.integrated <- FindClusters(immune.integrated, resolution = 0.3)
save(immune.integrated, file = 'immune.integrated.Rda')
Sys.time()

###### Glomerulus subclustering analysis ########
load("glom.Rda")
glom.list<-SplitObject(glom, split.by = 'orig.ident')
for (i in names(glom.list)) {
  glom.list[[i]] <- SCTransform(glom.list[[i]], verbose = F)
  print(paste(i, 'done!'))
}
glom.features <- SelectIntegrationFeatures(object.list = glom.list, nfeatures = 1000)

glom.list <- PrepSCTIntegration(object.list = glom.list, anchor.features = glom.features)

reference_dataset <- which(names(glom.list) == "A3020")

glom.anchors <- FindIntegrationAnchors(object.list = glom.list,
                                       anchor.features = glom.features, reference = reference_dataset, normalization.method = "SCT")Sys.time()
rm(glom.list)
gc()

for (i in 1:length(names(glom.anchors@object.list))){
  glom.anchors@object.list[[i]][['integrated']]<-NULL
}
glom <- IntegrateData(anchorset = glom.anchors, k.weight = 40)
Sys.time()
rm(glom.anchors)
gc()
glom<-ScaleData(glom)
glom <- RunPCA(object = glom, npcs = 30)
glom <- RunUMAP(object = glom, dims = 1:20)
glom <- FindNeighbors(glom, dims = 1:20)
glom <- FindClusters(glom, resolution = 0.1)
save(glom, file = 'glom.Rda')
Sys.time()


