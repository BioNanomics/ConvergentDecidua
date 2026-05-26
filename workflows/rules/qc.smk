# QC rules — filter, normalize, and quality-check datasets


rule qc_dataset:
    input:
        h5ad="results/processed/{accession}.h5ad",
    output:
        h5ad="results/qc/{accession}.h5ad",
    params:
        species=lambda wc: SPECIES_OF[wc.accession],
    shell:
        "wombat qc --species {params.species}"
