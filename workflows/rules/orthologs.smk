# Ortholog rules — build cross-species ortholog tables


rule build_orthologs:
    output:
        backbone="results/orthologs/backbone.parquet",
    shell:
        "wombat orthologs build"
