# Report rules — generate automated reports and manifests


rule generate_reports:
    input:
        scored="results/scored/stromal_scored.h5ad",
        backbone="results/orthologs/backbone.parquet",
    output:
        methods="results/reports/methods.md",
        coverage="results/reports/coverage.md",
        manifest="results/reports/manifest.md",
    shell:
        "wombat generate-reports"
