docker run -it --rm \
    -v /data:/data \
    aertslab/pyscenic:0.10.0 pyscenic ctx \
        /data/injPT_expr_mat.adjacencies.tsv \
        /data/mm9-tss-centered-5kb-7species.mc9nr.feather \
        /data/mm9-tss-centered-10kb-7species.mc9nr.feather \
        --annotations_fname /data/motifs-v9-nr.mgi-m0.001-o0.0.tbl \
        --expression_mtx_fname /data/injPT_filtered.loom \
        --mode "dask_multiprocessing" \
        --output /data/inj_regulons.csv \
        --num_workers 12
