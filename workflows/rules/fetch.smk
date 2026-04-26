# Fetch rules — download datasets and convert to standardized h5ad


rule fetch_dataset:
    """Download and convert a single dataset to h5ad."""
    output:
        "results/processed/{accession}.h5ad",
    shell:
        "wombat fetch --dataset {wildcards.accession}"
