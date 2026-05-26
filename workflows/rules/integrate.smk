# Integration rules — cross-species cell-state integration


rule integrate_stromal:
    input:
        backbone="results/orthologs/backbone.parquet",
        qc=expand("results/qc/{acc}.h5ad", acc=SCRNA_ACCESSIONS),
    output:
        h5ad="results/integrated/stromal_cross_species.h5ad",
    shell:
        "wombat integrate --mode stromal --method harmony"
