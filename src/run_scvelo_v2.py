import scvelo as scv
import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scipy.sparse as sp
import glob
import re

scv.settings.verbosity = 3
scv.settings.presenter_view = True
scv.set_figure_params('scvelo')

# Base directory for loom files
LOOM_DIR = '/gscratch/kawaldorflab/jcorn427/placenta/velocyto_out'


def extract_sample_id(filename):
    """Extract animal ID and tissue from filename"""
    animal_match = re.search(r'[JLRZAM]\d+', filename)

    if 'CV' in filename or '_cv_' in filename.lower():
        tissue = 'CV'
    elif 'decidua' in filename.lower() or 'dedidua' in filename.lower():
        tissue = 'Decidua'
    elif 'meminoc' in filename.lower():
        tissue = 'MemInoc'
    else:
        tissue = None

    if animal_match and tissue:
        return f"{animal_match.group()}_{tissue}"
    return None


def load_all_loom_files():
    """Load all loom files once and return concatenated AnnData object"""
    print("="*60)
    print("LOADING ALL LOOM FILES")
    print("="*60)

    loom_files = glob.glob(f'{LOOM_DIR}/*.loom')
    print(f"\nFound {len(loom_files)} loom files")

    adatas = []
    for loom_file in loom_files:
        sample_name = loom_file.split('/')[-1].replace('.loom', '')
        sample_id = extract_sample_id(loom_file)

        print(f"Loading {sample_name} -> {sample_id}...")
        adata_loom = sc.read_loom(loom_file)

        adata_loom.obs['loom_sample'] = sample_name
        adata_loom.obs['sample_id'] = sample_id
        adata_loom.var_names_make_unique()

        adatas.append(adata_loom)

    print("\nConcatenating all loom files...")
    if len(adatas) == 1:
        adata_all = adatas[0]
    else:
        adata_all = sc.concat(adatas, axis=0, join='outer')

    print(f"Combined loom data: {adata_all.shape[0]} cells, {adata_all.shape[1]} genes")

    # Extract barcode from format: "sample_name:BARCODEx"
    def extract_barcode(full_name):
        parts = full_name.split(':')
        if len(parts) == 2:
            barcode = parts[1].rstrip('x') + '-1'
            return barcode
        return full_name

    adata_all.obs['barcode_clean'] = [extract_barcode(x) for x in adata_all.obs_names]

    print(f"\nExample barcode conversion:")
    print(f"  Original: {adata_all.obs_names[0]}")
    print(f"  Cleaned:  {adata_all.obs['barcode_clean'].iloc[0]}")

    return adata_all


def run_velocity_analysis(adata_all, metadata_file, subset_name, cell_type_col='cell_type_v4',
                          use_dynamical=True, n_jobs=4):
    """
    Run scVelo analysis on a subset.

    Parameters:
    -----------
    adata_all : AnnData
        Pre-loaded concatenated loom data
    metadata_file : str
        Path to CSV with metadata and UMAP coordinates
    subset_name : str
        Name for saving outputs
    cell_type_col : str
        Column name for cell type annotations
    use_dynamical : bool
        If True, attempt dynamical mode first with stochastic fallback.
        If False, use stochastic mode only.
    n_jobs : int
        Number of parallel jobs for dynamical mode (recover_dynamics is slow)
    """

    print(f"\n{'='*60}")
    print(f"Processing {subset_name}")
    print(f"{'='*60}\n")

    # Load metadata
    metadata = pd.read_csv(metadata_file, index_col='barcode')
    print(f"Loaded metadata: {metadata.shape[0]} cells")

    # Extract barcode from metadata index
    metadata['barcode_clean'] = metadata.index.str.split('_').str[-1]
    metadata['sample_id'] = metadata['sample']

    print(f"Example metadata barcode: {metadata.index[0]}")
    print(f"  Extracted barcode: {metadata['barcode_clean'].iloc[0]}")
    print(f"  Sample ID: {metadata['sample_id'].iloc[0]}")

    print("\nMatching cells between loom and Seurat data...")

    # Vectorized matching
    metadata['match_key'] = metadata['sample_id'] + '_' + metadata['barcode_clean']
    adata_all.obs['match_key'] = adata_all.obs['sample_id'] + '_' + adata_all.obs['barcode_clean']

    common_keys = set(metadata['match_key']).intersection(set(adata_all.obs['match_key']))
    print(f"Found {len(common_keys)} matching cells")

    metadata_matched = metadata[metadata['match_key'].isin(common_keys)].copy()

    match_key_to_obs = dict(zip(adata_all.obs['match_key'], adata_all.obs_names))
    adata_obs_names = [match_key_to_obs[key] for key in metadata_matched['match_key']]

    adata = adata_all[adata_obs_names, :].copy()

    assert len(adata) == len(metadata_matched), \
        f"Cell count mismatch: adata has {len(adata)}, metadata has {len(metadata_matched)}"

    print(f"Successfully matched {len(metadata_matched)} out of {len(metadata)} cells "
          f"({100*len(metadata_matched)/len(metadata):.1f}%)")

    if len(metadata_matched) == 0:
        print("ERROR: No cells matched!")
        return None

    # Set obs_names to match Seurat barcodes
    adata.obs_names = metadata_matched.index.values

    # Add all metadata
    for col in metadata_matched.columns:
        if col not in ['barcode_clean', 'sample_id', 'match_key']:
            adata.obs[col] = metadata_matched[col].values

    # Add UMAP coordinates
    if 'umap_1' in metadata_matched.columns and 'umap_2' in metadata_matched.columns:
        adata.obsm['X_umap'] = metadata_matched[['umap_1', 'umap_2']].values
        print("Added UMAP coordinates from Seurat")

    # -------------------------------------------------------------------------
    # PRE-NORMALIZATION DIAGNOSTIC
    # Helps catch issues like pre-normalized counts or failed barcode matching.
    # Spliced/unspliced values should be raw integer counts before normalization.
    # -------------------------------------------------------------------------
    spliced_raw = adata.layers['spliced']
    unspliced_raw = adata.layers['unspliced']
    spliced_vals = spliced_raw.data if sp.issparse(spliced_raw) else spliced_raw.flatten()
    unspliced_vals = unspliced_raw.data if sp.issparse(unspliced_raw) else unspliced_raw.flatten()
    total_entries = adata.shape[0] * adata.shape[1]

    print(f"\n--- PRE-NORMALIZATION DIAGNOSTIC ---")
    print(f"Cells: {adata.shape[0]}, Genes: {adata.shape[1]}")
    print(f"Spliced   — nonzero mean: {spliced_vals.mean():.2f}, max: {spliced_vals.max():.0f}, "
          f"% nonzero: {100*len(spliced_vals)/total_entries:.2f}%, "
          f"integer counts: {np.all(spliced_vals == spliced_vals.astype(int))}")
    print(f"Unspliced — nonzero mean: {unspliced_vals.mean():.2f}, max: {unspliced_vals.max():.0f}, "
          f"% nonzero: {100*len(unspliced_vals)/total_entries:.2f}%, "
          f"integer counts: {np.all(unspliced_vals == unspliced_vals.astype(int))}")
    print(f"Sample spliced values (first 10 nonzero): {spliced_vals[:10]}")
    print(f"--- END DIAGNOSTIC ---\n")

    # -------------------------------------------------------------------------
    # PREPROCESS
    # -------------------------------------------------------------------------
    print("Preprocessing...")
    scv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=2000)
    scv.pp.moments(adata, n_pcs=30, n_neighbors=30)

    # -------------------------------------------------------------------------
    # VELOCITY COMPUTATION
    # Dynamical mode models full splicing kinetics and is preferred for myeloid
    # cells. Falls back to stochastic if the dynamical fit is degenerate.
    # -------------------------------------------------------------------------
    def _run_stochastic(adata):
        scv.tl.velocity(adata, mode='stochastic')
        scv.tl.velocity_graph(adata)
        return 'stochastic'

    velocity_mode_used = None

    if use_dynamical:
        try:
            print(f"\nAttempting dynamical mode (n_jobs={n_jobs})...")
            print("Note: recover_dynamics can take 10-30 min on large subsets.")
            scv.tl.recover_dynamics(adata, n_jobs=n_jobs)
            scv.tl.velocity(adata, mode='dynamical')
            scv.tl.velocity_graph(adata)

            # Check for degenerate fit — dynamical mode can succeed but produce
            # all-NaN velocities if splicing kinetics are unidentifiable
            frac_finite = np.isfinite(adata.layers['velocity']).mean()
            print(f"Fraction of finite velocity values: {frac_finite:.3f}")

            if frac_finite < 0.1:
                print("WARNING: >90% of velocity values are NaN — dynamical fit is degenerate.")
                print("Falling back to stochastic mode...")
                velocity_mode_used = _run_stochastic(adata)
            else:
                velocity_mode_used = 'dynamical'
                print("Dynamical mode successful.")
        except Exception as e:
            print(f"Dynamical mode failed: {e}")
            print("Falling back to stochastic mode...")
            velocity_mode_used = _run_stochastic(adata)
    else:
        print("\nUsing stochastic mode...")
        velocity_mode_used = _run_stochastic(adata)

    print(f"Velocity computed using: {velocity_mode_used} mode")

    # -------------------------------------------------------------------------
    # NaN FILTERING
    #
    # KEY INSIGHT: velocity is only computed for velocity genes — a subset of
    # all 2000 genes. All other columns in the velocity layer are NaN by design.
    # Using .all(axis=1) would filter out every single cell because no cell
    # has finite values across all gene columns.
    #
    # Fix: identify velocity genes (columns with any finite value), then only
    # remove cells that have NO finite velocity values whatsoever.
    # -------------------------------------------------------------------------
    print("\nFiltering cells with non-finite velocity values...")

    velocity_layer = adata.layers['velocity']
    vel_dense = velocity_layer.toarray() if sp.issparse(velocity_layer) else np.array(velocity_layer)

    # Velocity genes = columns that have at least one finite value
    velocity_gene_mask = np.isfinite(vel_dense).any(axis=0)
    n_velocity_genes = velocity_gene_mask.sum()
    print(f"Velocity genes identified: {n_velocity_genes} / {vel_dense.shape[1]}")

    if n_velocity_genes == 0:
        print(f"ERROR: No velocity genes found for {subset_name}. "
              f"Splicing signal may be too weak.")
        return None

    # Per-cell diagnostics across velocity genes only
    vel_subset = vel_dense[:, velocity_gene_mask]
    frac_finite_per_cell = np.isfinite(vel_subset).mean(axis=1)
    print(f"Per-cell finite velocity fraction — "
          f"mean: {frac_finite_per_cell.mean():.3f}, "
          f"median: {np.median(frac_finite_per_cell):.3f}, "
          f"min: {frac_finite_per_cell.min():.3f}")

    # Only remove cells with zero finite velocity values (truly empty cells)
    valid_mask = np.isfinite(vel_subset).any(axis=1)
    n_invalid = (~valid_mask).sum()

    if n_invalid > 0:
        print(f"Removing {n_invalid} cells with no finite velocity values "
              f"({100*n_invalid/len(adata):.1f}%)")
        adata = adata[valid_mask, :].copy()
        print(f"Remaining cells: {len(adata)}")
    else:
        print("All cells have valid velocity values")

    if len(adata) == 0:
        print(f"ERROR: No valid cells remain after filtering for {subset_name}.")
        print("Check the pre-normalization diagnostic above for clues.")
        return None

    # -------------------------------------------------------------------------
    # PLOTTING
    # -------------------------------------------------------------------------
    print("\nGenerating plots...")

    umap_min_x = adata.obsm['X_umap'][:, 0].min() - 1
    umap_max_x = adata.obsm['X_umap'][:, 0].max() + 1
    umap_min_y = adata.obsm['X_umap'][:, 1].min() - 1
    umap_max_y = adata.obsm['X_umap'][:, 1].max() + 1

    print(f"Plot limits: x=({umap_min_x:.1f}, {umap_max_x:.1f}), "
          f"y=({umap_min_y:.1f}, {umap_max_y:.1f})")

    print("Recomputing velocity embedding for visualization...")
    scv.tl.velocity_embedding(adata, basis='umap')

    # --- STREAM PLOTS (SVG) ---
    print("\n--- Creating stream plots (SVG format) ---")

    # 1. Stream by cell type
    try:
        print("Creating velocity stream plot by cell type...")
        fig, ax = plt.subplots(figsize=(10, 8))
        scv.pl.velocity_embedding_stream(
            adata,
            basis='umap',
            color=cell_type_col,
            legend_loc='right margin',
            title=f'{subset_name} ({velocity_mode_used}) - Velocity Stream by Cell Type',
            ax=ax,
            show=False
        )
        ax.set_xlim(umap_min_x, umap_max_x)
        ax.set_ylim(umap_min_y, umap_max_y)
        plt.tight_layout()
        plt.savefig(f'{subset_name}_stream_celltype.svg', dpi=300, bbox_inches='tight')
        plt.close()
        print(f"Saved: {subset_name}_stream_celltype.svg")
    except Exception as e:
        print(f"Error creating stream plot by cell type: {e}")
        plt.close()

    # 2. Stream by tissue
    if 'tissue' in adata.obs.columns:
        try:
            print("Creating velocity stream plot by tissue...")
            fig, ax = plt.subplots(figsize=(10, 8))
            scv.pl.velocity_embedding_stream(
                adata,
                basis='umap',
                color='tissue',
                legend_loc='right margin',
                title=f'{subset_name} ({velocity_mode_used}) - Velocity Stream by Tissue',
                ax=ax,
                show=False
            )
            ax.set_xlim(umap_min_x, umap_max_x)
            ax.set_ylim(umap_min_y, umap_max_y)
            plt.tight_layout()
            plt.savefig(f'{subset_name}_stream_tissue.svg', dpi=300, bbox_inches='tight')
            plt.close()
            print(f"Saved: {subset_name}_stream_tissue.svg")
        except Exception as e:
            print(f"Error creating stream plot by tissue: {e}")
            plt.close()

    # --- ARROW PLOTS (PDF) ---
    print("\n--- Creating arrow plots (PDF format) ---")

    # 3. Arrows by cell type
    try:
        print("Creating velocity arrows plot by cell type...")
        fig, ax = plt.subplots(figsize=(10, 8))
        scv.pl.velocity_embedding(
            adata,
            basis='umap',
            color=cell_type_col,
            arrow_length=3,
            arrow_size=2,
            legend_loc='right margin',
            title=f'{subset_name} ({velocity_mode_used}) - Velocity Arrows by Cell Type',
            ax=ax,
            show=False
        )
        ax.set_xlim(umap_min_x, umap_max_x)
        ax.set_ylim(umap_min_y, umap_max_y)
        plt.tight_layout()
        plt.savefig(f'{subset_name}_arrows_celltype.pdf', dpi=300, bbox_inches='tight')
        plt.close()
        print(f"Saved: {subset_name}_arrows_celltype.pdf")
    except Exception as e:
        print(f"Error creating arrows plot by cell type: {e}")
        plt.close()

    # 4. Arrows by tissue
    if 'tissue' in adata.obs.columns:
        try:
            print("Creating velocity arrows plot by tissue...")
            fig, ax = plt.subplots(figsize=(10, 8))
            scv.pl.velocity_embedding(
                adata,
                basis='umap',
                color='tissue',
                arrow_length=3,
                arrow_size=2,
                legend_loc='right margin',
                title=f'{subset_name} ({velocity_mode_used}) - Velocity Arrows by Tissue',
                ax=ax,
                show=False
            )
            ax.set_xlim(umap_min_x, umap_max_x)
            ax.set_ylim(umap_min_y, umap_max_y)
            plt.tight_layout()
            plt.savefig(f'{subset_name}_arrows_tissue.pdf', dpi=300, bbox_inches='tight')
            plt.close()
            print(f"Saved: {subset_name}_arrows_tissue.pdf")
        except Exception as e:
            print(f"Error creating arrows plot by tissue: {e}")
            plt.close()

    # --- SCATTER PLOTS (PDF) ---
    print("\n--- Creating scatter plots (PDF format) ---")

    # 5. Velocity confidence
    scv.tl.velocity_confidence(adata)
    try:
        print("Creating velocity confidence plot...")
        fig, ax = plt.subplots(figsize=(10, 8))
        umap_coords = adata.obsm['X_umap']
        confidence = adata.obs['velocity_confidence'].values
        scatter = ax.scatter(umap_coords[:, 0], umap_coords[:, 1],
                             c=confidence, cmap='coolwarm', s=5, alpha=0.8,
                             vmin=np.percentile(confidence, 5),
                             vmax=np.percentile(confidence, 95))
        ax.set_xlim(umap_min_x, umap_max_x)
        ax.set_ylim(umap_min_y, umap_max_y)
        ax.set_xlabel('UMAP_1')
        ax.set_ylabel('UMAP_2')
        ax.set_title(f'{subset_name} ({velocity_mode_used}) - Velocity Confidence')
        cbar = plt.colorbar(scatter, ax=ax)
        cbar.set_label('velocity_confidence')
        plt.tight_layout()
        plt.savefig(f'{subset_name}_confidence.pdf', dpi=300, bbox_inches='tight')
        plt.close()
        print(f"Saved: {subset_name}_confidence.pdf")
    except Exception as e:
        print(f"Error creating confidence plot: {e}")
        plt.close()

    # 6. Velocity pseudotime
    # Call terminal_states explicitly to avoid hanging inside velocity_pseudotime
    print("Computing terminal states explicitly...")
    try:
        scv.tl.terminal_states(adata, self_transitions=True, n_jobs=n_jobs)
    except Exception as e:
        print(f"Warning: terminal_states failed ({e}), pseudotime may be less accurate")

    scv.tl.velocity_pseudotime(adata)
    try:
        print("Creating velocity pseudotime plot...")
        fig, ax = plt.subplots(figsize=(10, 8))
        umap_coords = adata.obsm['X_umap']
        pseudotime = adata.obs['velocity_pseudotime'].values
        scatter = ax.scatter(umap_coords[:, 0], umap_coords[:, 1],
                             c=pseudotime, cmap='gnuplot', s=5, alpha=0.8)
        ax.set_xlim(umap_min_x, umap_max_x)
        ax.set_ylim(umap_min_y, umap_max_y)
        ax.set_xlabel('UMAP_1')
        ax.set_ylabel('UMAP_2')
        ax.set_title(f'{subset_name} ({velocity_mode_used}) - Velocity Pseudotime')
        cbar = plt.colorbar(scatter, ax=ax)
        cbar.set_label('velocity_pseudotime')
        plt.tight_layout()
        plt.savefig(f'{subset_name}_pseudotime.pdf', dpi=300, bbox_inches='tight')
        plt.close()
        print(f"Saved: {subset_name}_pseudotime.pdf")
    except Exception as e:
        print(f"Error creating pseudotime plot: {e}")
        plt.close()

    # 7. Dynamical mode extras: latent time and top likelihood genes
    if velocity_mode_used == 'dynamical':
        try:
            print("\nComputing latent time (dynamical mode only)...")
            scv.tl.latent_time(adata)
            fig, ax = plt.subplots(figsize=(10, 8))
            umap_coords = adata.obsm['X_umap']
            latent_time = adata.obs['latent_time'].values
            scatter = ax.scatter(umap_coords[:, 0], umap_coords[:, 1],
                                 c=latent_time, cmap='gnuplot', s=5, alpha=0.8)
            ax.set_xlim(umap_min_x, umap_max_x)
            ax.set_ylim(umap_min_y, umap_max_y)
            ax.set_xlabel('UMAP_1')
            ax.set_ylabel('UMAP_2')
            ax.set_title(f'{subset_name} - Latent Time')
            cbar = plt.colorbar(scatter, ax=ax)
            cbar.set_label('latent_time')
            plt.tight_layout()
            plt.savefig(f'{subset_name}_latent_time.pdf', dpi=300, bbox_inches='tight')
            plt.close()
            print(f"Saved: {subset_name}_latent_time.pdf")
        except Exception as e:
            print(f"Error computing latent time: {e}")
            plt.close()

        try:
            print("Plotting top likelihood genes (dynamical mode only)...")
            top_genes = adata.var['fit_likelihood'].sort_values(ascending=False).head(10).index
            scv.pl.scatter(adata, basis=top_genes, ncols=5, frameon=False,
                           save=f'{subset_name}_top_likelihood_genes.pdf')
            print(f"Saved: {subset_name}_top_likelihood_genes.pdf")
        except Exception as e:
            print(f"Error plotting likelihood genes: {e}")

    # 8. Top velocity genes
    if len(adata.obs[cell_type_col].unique()) > 1:
        try:
            print("Ranking velocity genes...")
            scv.tl.rank_velocity_genes(adata, groupby=cell_type_col, min_corr=.3)
            df = pd.DataFrame(adata.uns['rank_velocity_genes']['names'])
            print(f"\nTop velocity genes by {cell_type_col}:")
            print(df.head(10))
            df.to_csv(f'{subset_name}_top_velocity_genes.csv')
            print(f"Saved: {subset_name}_top_velocity_genes.csv")
        except Exception as e:
            print(f"Error ranking velocity genes: {e}")

    # Save
    adata.write(f'{subset_name}_velocity.h5ad')
    print(f"\nSaved results to {subset_name}_velocity.h5ad")
    print(f"Velocity mode used: {velocity_mode_used}")

    return adata


# =============================================================================
# MAIN
# =============================================================================
if __name__ == "__main__":

    print("Starting scVelo trajectory analysis on SUBSETS...")
    print(f"Loom files directory: {LOOM_DIR}\n")

    # Load all loom files once
    adata_all = load_all_loom_files()

    print("\n" + "="*60)
    print("LOOM FILES LOADED - NOW PROCESSING SUBSETS")
    print("="*60)

    # # DSC subset
    # dsc_adata = run_velocity_analysis(
    #     adata_all=adata_all,
    #     metadata_file='dsc_metadata.csv',
    #     subset_name='DSC',
    #     cell_type_col='cell_type_v5',
    #     use_dynamical=True,
    #     n_jobs=4
    # )

    # MAC subset — dynamical mode preferred for myeloid cells
    #mac_adata = run_velocity_analysis(
    #   adata_all=adata_all,
    #    metadata_file='mac_metadata.csv',
    #    subset_name='MAC',
    #    cell_type_col='cell_type_v5',
    #    use_dynamical=True,
    #    n_jobs=4
    #)

    # # T cell subset
    # t_adata = run_velocity_analysis(
    #     adata_all=adata_all,
    #     metadata_file='t_cell_metadata.csv',
    #     subset_name='t_cell',
    #     cell_type_col='cell_type_v5',
    #     use_dynamical=True,
    #     n_jobs=4
    # )

    # TB subset
    ctb_adata = run_velocity_analysis(
        adata_all=adata_all,
        metadata_file='tb_traj_metadata.csv',
        subset_name='TB',
        cell_type_col='cell_type_v5',
        use_dynamical=True,
        n_jobs=4
    )

    # # NK subset
    # unk_adata = run_velocity_analysis(
    #     adata_all=adata_all,
    #     metadata_file='nk_metadata.csv',
    #     subset_name='uNK',
    #     cell_type_col='cell_type_v5',
    #     use_dynamical=True,
    #     n_jobs=4
    # )

    print("\n" + "="*60)
    print("ALL ANALYSES COMPLETE!")
    print("="*60)