docker run -it --rm \
    -v /data:/data \
    aertslab/pyscenic:0.10.0 pyscenic aucell \
        /data/injPT_filtered.loom \
        /data/inj_regulons.csv \
        -o /data/inj_auc_mtx.csv \
        --num_workers 12
 
