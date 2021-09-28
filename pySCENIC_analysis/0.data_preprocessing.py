import scanpy as sc
import numpy as np
import loompy as lp
adata = sc.read("injPT.h5ad")
row_attrs = {
  "Gene": np.array(adata.var_names),
}
col_attrs = {
  "CellID": np.array(adata.obs_names),
  "nGene": np.array(np.sum(adata.X.transpose()>0, axis=0)).flatten(),
  "nUMI": np.array(np.sum(adata.X.transpose(),axis=0)).flatten(),
}
lp.create("injPT_filtered.loom",adata.X.transpose(),row_attrs,
          col_attrs)
