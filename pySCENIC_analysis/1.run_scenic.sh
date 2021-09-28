docker run -it --rm \
    -v /data:/data \
    aertslab/pyscenic:0.10.0 pyscenic grn \
        --num_workers 12 \
        -o /data/injPT_expr_mat.adjacencies.tsv \
        /data/injPT_filtered.loom \
        /data/mm_mgi_tfs.txt
