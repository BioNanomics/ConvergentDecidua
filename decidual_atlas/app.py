"""DecidualAtlas — interactive Streamlit visualization app.

Multi-page app for exploring the comparative decidualization atlas.
Launch via: ``wombat serve-atlas``
"""

from __future__ import annotations

from pathlib import Path

import streamlit as st

# Project paths
PROJECT_ROOT = Path(__file__).resolve().parent.parent
RESULTS_DIR = PROJECT_ROOT / "results"

st.set_page_config(
    page_title="DecidualAtlas",
    page_icon="🔬",
    layout="wide",
    initial_sidebar_state="expanded",
)


def main() -> None:
    """Main app entry point with sidebar navigation."""
    st.sidebar.title("🔬 DecidualAtlas")
    st.sidebar.markdown("*Comparative Decidualization Atlas*")
    st.sidebar.divider()

    page = st.sidebar.radio(
        "Navigate",
        [
            "🏠 Overview",
            "📋 Dataset Browser",
            "🧬 Species Comparison",
            "🔎 Cell-State Viewer",
            "🧪 Gene Explorer",
        ],
    )

    if page == "🏠 Overview":
        _page_overview()
    elif page == "📋 Dataset Browser":
        _page_dataset_browser()
    elif page == "🧬 Species Comparison":
        _page_species_comparison()
    elif page == "🔎 Cell-State Viewer":
        _page_cell_state_viewer()
    elif page == "🧪 Gene Explorer":
        _page_gene_explorer()


# ---------------------------------------------------------------------------
# Pages
# ---------------------------------------------------------------------------


def _page_overview() -> None:
    st.title("DecidualAtlas — Overview")
    st.markdown(
        """
        **ConvergentDecidua** builds a comparative atlas of decidualization
        across species. This viewer lets you explore:

        - **Dataset Browser**: Registry of all ingested datasets
        - **Species Comparison**: Side-by-side UMAPs and score distributions
        - **Cell-State Viewer**: Interactive embedding colored by type/species/score
        - **Gene Explorer**: Search expression by gene symbol
        """
    )

    # Show summary stats if data exists
    registry_path = RESULTS_DIR / "registry.parquet"
    if registry_path.exists():
        import pandas as pd

        df = pd.read_parquet(registry_path)
        col1, col2, col3 = st.columns(3)
        col1.metric("Datasets", len(df))
        col2.metric("Species", df["species"].nunique())
        col3.metric("Assays", df["assay"].nunique())
    else:
        st.info("No registry found. Run `wombat build-registry` first.")


def _page_dataset_browser() -> None:
    st.title("📋 Dataset Browser")

    registry_path = RESULTS_DIR / "registry.parquet"
    if not registry_path.exists():
        st.warning("Registry not found. Run `wombat build-registry`.")
        return

    import pandas as pd

    df = pd.read_parquet(registry_path)

    # Filters
    col1, col2 = st.columns(2)
    with col1:
        species_filter = st.multiselect(
            "Species", df["species"].unique(), default=list(df["species"].unique())
        )
    with col2:
        assay_filter = st.multiselect(
            "Assay", df["assay"].unique(), default=list(df["assay"].unique())
        )

    filtered = df[df["species"].isin(species_filter) & df["assay"].isin(assay_filter)]
    st.dataframe(filtered, use_container_width=True)


def _page_species_comparison() -> None:
    st.title("🧬 Species Comparison")

    integrated_path = RESULTS_DIR / "integrated" / "stromal_harmony.h5ad"
    if not integrated_path.exists():
        st.warning("Integrated data not found. Run `wombat integrate` first.")
        return

    adata = _load_h5ad_cached(str(integrated_path))

    if "X_umap" not in adata.obsm:
        st.warning("UMAP not computed in integrated data.")
        return

    import plotly.express as px

    umap_df = _get_umap_df(adata)

    # Side-by-side UMAPs
    col1, col2 = st.columns(2)
    for sp, col in zip(["human", "mouse"], [col1, col2], strict=False):
        subset = umap_df[umap_df["species"] == sp]
        with col:
            st.subheader(sp.capitalize())
            if len(subset) > 0:
                fig = px.scatter(
                    subset,
                    x="UMAP1",
                    y="UMAP2",
                    color="cell_type",
                    title=f"{sp.capitalize()} ({len(subset)} cells)",
                    opacity=0.6,
                    width=500,
                    height=450,
                )
                fig.update_traces(marker_size=3)
                st.plotly_chart(fig, use_container_width=True)
            else:
                st.info(f"No {sp} cells in integrated data")

    # Score comparison
    score_cols = [c for c in adata.obs.columns if c.endswith("_score")]
    if score_cols:
        st.subheader("Score Distributions by Species")
        selected_score = st.selectbox("Select score module", score_cols)
        fig = px.violin(
            umap_df,
            x="species",
            y=selected_score,
            color="species",
            box=True,
            title=selected_score,
        )
        st.plotly_chart(fig, use_container_width=True)


def _page_cell_state_viewer() -> None:
    st.title("🔎 Cell-State Viewer")

    scored_path = RESULTS_DIR / "scored" / "stromal_scored.h5ad"
    integrated_path = RESULTS_DIR / "integrated" / "stromal_harmony.h5ad"
    data_path = scored_path if scored_path.exists() else integrated_path

    if not data_path.exists():
        st.warning("No integrated/scored data found.")
        return

    import plotly.express as px

    adata = _load_h5ad_cached(str(data_path))

    if "X_umap" not in adata.obsm:
        st.warning("UMAP not computed.")
        return

    umap_df = _get_umap_df(adata)

    # Color-by selector
    color_options = ["cell_type", "species", "dataset"]
    score_cols = [c for c in adata.obs.columns if c.endswith("_score")]
    color_options.extend(score_cols)
    color_by = st.selectbox("Color by", color_options)

    fig = px.scatter(
        umap_df,
        x="UMAP1",
        y="UMAP2",
        color=color_by,
        title=f"UMAP — colored by {color_by}",
        opacity=0.6,
        width=900,
        height=600,
    )
    fig.update_traces(marker_size=3)
    st.plotly_chart(fig, use_container_width=True)


def _page_gene_explorer() -> None:
    st.title("🧪 Gene Explorer")

    scored_path = RESULTS_DIR / "scored" / "stromal_scored.h5ad"
    integrated_path = RESULTS_DIR / "integrated" / "stromal_harmony.h5ad"
    data_path = scored_path if scored_path.exists() else integrated_path

    if not data_path.exists():
        st.warning("No data found.")
        return

    adata = _load_h5ad_cached(str(data_path))

    gene = st.text_input("Gene symbol (human)", value="PRL")
    if gene and gene in adata.var_names:
        import numpy as np
        import plotly.express as px

        expr = (
            np.asarray(adata[:, gene].X.todense()).flatten()
            if hasattr(adata[:, gene].X, "todense")
            else np.asarray(adata[:, gene].X).flatten()
        )

        umap_df = _get_umap_df(adata)
        umap_df["expression"] = expr

        fig = px.scatter(
            umap_df,
            x="UMAP1",
            y="UMAP2",
            color="expression",
            color_continuous_scale="Viridis",
            title=f"{gene} expression",
            opacity=0.7,
            width=900,
            height=600,
        )
        fig.update_traces(marker_size=3)
        st.plotly_chart(fig, use_container_width=True)

        # Also show violin by cell type
        if "cell_type" in umap_df.columns:
            fig2 = px.violin(
                umap_df,
                x="cell_type",
                y="expression",
                color="species" if "species" in umap_df.columns else None,
                title=f"{gene} by cell type",
                box=True,
            )
            st.plotly_chart(fig2, use_container_width=True)

        # DuckDB query on backbone
        _gene_ortholog_lookup(gene)

    elif gene:
        st.warning(f"Gene '{gene}' not found in dataset. Available: {adata.n_vars} genes.")


def _gene_ortholog_lookup(gene: str) -> None:
    """Look up orthologs via DuckDB on backbone parquet."""
    backbone_path = RESULTS_DIR / "orthologs" / "backbone.parquet"
    if not backbone_path.exists():
        return

    try:
        import duckdb

        conn = duckdb.connect(":memory:")
        result = conn.execute(
            f"SELECT * FROM read_parquet('{backbone_path}') "
            f"WHERE source_symbol = ? OR target_symbol = ?",
            [gene, gene],
        ).fetchdf()
        conn.close()

        if len(result) > 0:
            st.subheader("Ortholog Mapping")
            st.dataframe(result, use_container_width=True)
    except ImportError:
        pass


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


@st.cache_data(show_spinner="Loading data...")
def _load_h5ad_cached(path: str):
    """Load h5ad with Streamlit caching."""
    import anndata as ad

    return ad.read_h5ad(path)


def _get_umap_df(adata):
    """Extract UMAP coordinates + obs into a DataFrame."""
    import pandas as pd

    df = pd.DataFrame(
        adata.obsm["X_umap"],
        columns=["UMAP1", "UMAP2"],
        index=adata.obs.index,
    )
    for col in adata.obs.columns:
        df[col] = adata.obs[col].values
    return df


if __name__ == "__main__":
    main()
