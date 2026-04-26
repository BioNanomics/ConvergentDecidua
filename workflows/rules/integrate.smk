# Integration rules — cross-species cell-state integration


rule integrate_stromal:
    input:
        backbone="results/orthologs/backbone.parquet",
        qc=expand("results/qc/{acc}.h5ad", acc=[
            d["accession"] for d in config
            if "scrna" in d.get("assay", "").lower() or "snrna" in d.get("assay", "").lower()
        ]),
    output:
        h5ad="results/integrated/stromal_harmony.h5ad",
    shell:
        "wombat integrate --mode stromal --method harmony"
