# Report rules — generate automated reports and manifests


rule generate_reports:
    input:
        integrated="results/integrated/stromal_cross_species.h5ad",
        backbone="results/orthologs/backbone.parquet",
    output:
        methods="results/reports/methods.md",
        coverage="results/reports/coverage.md",
        qc_summary="results/reports/qc_summary.md",
        orthologs="results/reports/orthologs.md",
        integration_qc="results/reports/integration_qc.md",
        manifest="results/reports/manifest.md",
    shell:
        "wombat generate-reports"
