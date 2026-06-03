import scvelo as scv
import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
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


def run_velocity_full_dataset(adata_all, metadata_file, cell_type_col='predicted.celltype'):
    """
    Run scVelo analysis on the full dataset
    
    Parameters:
    -----------
    adata_all : AnnData
        Pre-loaded concatenated loom data
    metadata_file : str
        Path to CSV with metadata and UMAP coordinates
    cell_type_col : str
        Column name for cell type annotations
    """
    
    print(f"\n{'='*60}")
    print(f"Processing FULL DATASET")
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
    
    # OPTIMIZED VECTORIZED MATCHING
    # Create matching keys for both datasets
    metadata['match_key'] = metadata['sample_id'] + '_' + metadata['barcode_clean']
    adata_all.obs['match_key'] = adata_all.obs['sample_id'] + '_' + adata_all.obs['barcode_clean']
    
    # Find intersection of match keys
    common_keys = set(metadata['match_key']).intersection(set(adata_all.obs['match_key']))
    print(f"Found {len(common_keys)} matching cells")
    
    # Filter both datasets to common cells
    metadata_matched = metadata[metadata['match_key'].isin(common_keys)].copy()
    
    # Create a dictionary mapping match_key to adata obs_names for faster lookup
    match_key_to_obs = dict(zip(adata_all.obs['match_key'], adata_all.obs_names))
    
    # Get the adata obs_names corresponding to the matched metadata
    adata_obs_names = [match_key_to_obs[key] for key in metadata_matched['match_key']]
    
    # Subset adata using the matched obs_names
    adata = adata_all[adata_obs_names, :].copy()
    
    # Verify we have the right number of cells
    assert len(adata) == len(metadata_matched), f"Cell count mismatch: adata has {len(adata)}, metadata has {len(metadata_matched)}"
    
    print(f"Successfully matched {len(metadata_matched)} out of {len(metadata)} cells ({100*len(metadata_matched)/len(metadata):.1f}%)")
    
    if len(metadata_matched) == 0:
        print("ERROR: No cells matched!")
        return None
    
    # Set the obs_names to match Seurat barcodes
    adata.obs_names = metadata_matched.index.values
    
    # Add all metadata
    for col in metadata_matched.columns:
        if col not in ['barcode_clean', 'sample_id', 'match_key']:
            adata.obs[col] = metadata_matched[col].values
    
    print("Columns after metadata transfer:", adata.obs.columns.tolist())
    print("cell_type_v1 sample:", adata.obs['cell_type_v1'].head())
    
    # Add UMAP coordinates
    if 'umap_1' in metadata_matched.columns and 'umap_2' in metadata_matched.columns:
        adata.obsm['X_umap'] = metadata_matched[['umap_1', 'umap_2']].values
        print("Added UMAP coordinates from Seurat")
    
    # Preprocess
    print("\nPreprocessing...")
    scv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=2000)
    scv.pp.moments(adata, n_pcs=30, n_neighbors=30)
    
    # Compute velocity
    print("Computing RNA velocity...")
    scv.tl.velocity(adata, mode='stochastic')
    scv.tl.velocity_graph(adata)
    
    # COMPREHENSIVE NaN FILTERING
    print("\nFiltering cells with non-finite values...")
    
    # Check velocity
    valid_velocity = np.isfinite(adata.layers['velocity']).all(axis=1)
    
    # Check velocity_umap if it exists
    if 'velocity_umap' in adata.obsm.keys():
        valid_velocity_umap = np.isfinite(adata.obsm['velocity_umap']).all(axis=1)
        valid_mask = valid_velocity & valid_velocity_umap
    else:
        valid_mask = valid_velocity
    
    n_invalid = (~valid_mask).sum()
    
    if n_invalid > 0:
        print(f"Removing {n_invalid} cells with non-finite values ({100*n_invalid/len(adata):.1f}%)")
        adata = adata[valid_mask, :].copy()
        print("Columns after NaN filtering:", adata.obs.columns.tolist())
        print(f"Remaining cells: {len(adata)}")
    else:
        print("All cells have valid values")
    
    # Generate plots
    print("\nGenerating plots...")
    
    # Get UMAP coordinate ranges for proper axis limits
    umap_min_x = adata.obsm['X_umap'][:, 0].min() - 1
    umap_max_x = adata.obsm['X_umap'][:, 0].max() + 1
    umap_min_y = adata.obsm['X_umap'][:, 1].min() - 1
    umap_max_y = adata.obsm['X_umap'][:, 1].max() + 1
    
    print(f"Setting plot limits: x=({umap_min_x:.1f}, {umap_max_x:.1f}), y=({umap_min_y:.1f}, {umap_max_y:.1f})")
    
    # Recompute velocity embedding after filtering
    print("Recomputing velocity embedding for visualization...")
    scv.tl.velocity_embedding(adata, basis='umap')
    
    # STREAM PLOTS (as SVG)
    print("\n--- Creating stream plots (SVG format) ---")
    
    # 1. Velocity stream plot colored by cell type
    try:
        print("Creating velocity stream plot by cell type...")
        fig, ax = plt.subplots(figsize=(10, 8))
        scv.pl.velocity_embedding_stream(
            adata, 
            basis='umap', 
            color=cell_type_col,
            legend_loc='right margin',
            title='Full Dataset - Velocity Stream by Cell Type',
            ax=ax,
            show=False
        )
        ax.set_xlim(umap_min_x, umap_max_x)
        ax.set_ylim(umap_min_y, umap_max_y)
        plt.tight_layout()
        plt.savefig('full_dataset_stream_celltype.svg', dpi=300, bbox_inches='tight')
        plt.close()
        print("Saved: full_dataset_stream_celltype.svg")
    except Exception as e:
        print(f"Error creating stream plot by cell type: {e}")
        plt.close()
    
    # 2. Velocity stream colored by tissue
    if 'tissue' in adata.obs.columns:
        try:
            print("Creating velocity stream plot by tissue...")
            fig, ax = plt.subplots(figsize=(10, 8))
            scv.pl.velocity_embedding_stream(
                adata, 
                basis='umap', 
                color='tissue',
                legend_loc='right margin',
                title='Full Dataset - Velocity Stream by Tissue',
                ax=ax,
                show=False
            )
            ax.set_xlim(umap_min_x, umap_max_x)
            ax.set_ylim(umap_min_y, umap_max_y)
            plt.tight_layout()
            plt.savefig('full_dataset_stream_tissue.svg', dpi=300, bbox_inches='tight')
            plt.close()
            print("Saved: full_dataset_stream_tissue.svg")
        except Exception as e:
            print(f"Error creating stream plot by tissue: {e}")
            plt.close()
    
    # 3. Velocity stream colored by sample
    try:
        print("Creating velocity stream plot by sample...")
        fig, ax = plt.subplots(figsize=(10, 8))
        scv.pl.velocity_embedding_stream(
            adata, 
            basis='umap', 
            color='sample',
            legend_loc='right margin',
            title='Full Dataset - Velocity Stream by Sample',
            legend_fontsize=6,
            ax=ax,
            show=False
        )
        ax.set_xlim(umap_min_x, umap_max_x)
        ax.set_ylim(umap_min_y, umap_max_y)
        plt.tight_layout()
        plt.savefig('full_dataset_stream_sample.svg', dpi=300, bbox_inches='tight')
        plt.close()
        print("Saved: full_dataset_stream_sample.svg")
    except Exception as e:
        print(f"Error creating stream plot by sample: {e}")
        plt.close()
    
    # ARROW PLOTS (as PDF)
    print("\n--- Creating arrow plots (PDF format) ---")
    
    # 4. Velocity arrows colored by cell type
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
            title='Full Dataset - Velocity Arrows by Cell Type',
            ax=ax,
            show=False
        )
        ax.set_xlim(umap_min_x, umap_max_x)
        ax.set_ylim(umap_min_y, umap_max_y)
        plt.tight_layout()
        plt.savefig('full_dataset_arrows_celltype.pdf', dpi=300, bbox_inches='tight')
        plt.close()
        print("Saved: full_dataset_arrows_celltype.pdf")
    except Exception as e:
        print(f"Error creating arrows plot by cell type: {e}")
        plt.close()
    
    # 5. Velocity arrows colored by tissue
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
                title='Full Dataset - Velocity Arrows by Tissue',
                ax=ax,
                show=False
            )
            ax.set_xlim(umap_min_x, umap_max_x)
            ax.set_ylim(umap_min_y, umap_max_y)
            plt.tight_layout()
            plt.savefig('full_dataset_arrows_tissue.pdf', dpi=300, bbox_inches='tight')
            plt.close()
            print("Saved: full_dataset_arrows_tissue.pdf")
        except Exception as e:
            print(f"Error creating arrows plot by tissue: {e}")
            plt.close()
    
    # 6. Velocity arrows colored by sample
    try:
        print("Creating velocity arrows plot by sample...")
        fig, ax = plt.subplots(figsize=(10, 8))
        scv.pl.velocity_embedding(
            adata, 
            basis='umap', 
            color='sample',
            arrow_length=3,
            arrow_size=2,
            legend_loc='right margin',
            title='Full Dataset - Velocity Arrows by Sample',
            legend_fontsize=6,
            ax=ax,
            show=False
        )
        ax.set_xlim(umap_min_x, umap_max_x)
        ax.set_ylim(umap_min_y, umap_max_y)
        plt.tight_layout()
        plt.savefig('full_dataset_arrows_sample.pdf', dpi=300, bbox_inches='tight')
        plt.close()
        print("Saved: full_dataset_arrows_sample.pdf")
    except Exception as e:
        print(f"Error creating arrows plot by sample: {e}")
        plt.close()
    
    # SCATTER PLOTS (as PDF)
    print("\n--- Creating scatter plots (PDF format) ---")
    
    # 7. Velocity confidence - using manual matplotlib scatter
    scv.tl.velocity_confidence(adata)
    try:
        print("Creating velocity confidence plot...")
        fig, ax = plt.subplots(figsize=(10, 8))
        
        # Get UMAP coordinates and confidence values
        umap_coords = adata.obsm['X_umap']
        confidence = adata.obs['velocity_confidence'].values
        
        # Create scatter plot manually
        scatter = ax.scatter(umap_coords[:, 0], umap_coords[:, 1], 
                           c=confidence, cmap='coolwarm', s=5, alpha=0.8,
                           vmin=np.percentile(confidence, 5), 
                           vmax=np.percentile(confidence, 95))
        
        ax.set_xlim(umap_min_x, umap_max_x)
        ax.set_ylim(umap_min_y, umap_max_y)
        ax.set_xlabel('UMAP_1')
        ax.set_ylabel('UMAP_2')
        ax.set_title('Full Dataset - Velocity Confidence')
        
        # Add colorbar
        cbar = plt.colorbar(scatter, ax=ax)
        cbar.set_label('velocity_confidence')
        
        plt.tight_layout()
        plt.savefig('full_dataset_confidence.pdf', dpi=300, bbox_inches='tight')
        plt.close()
        print("Saved: full_dataset_confidence.pdf")
    except Exception as e:
        print(f"Error creating confidence plot: {e}")
        plt.close()
    
    # 8. Velocity pseudotime - using manual matplotlib scatter
    scv.tl.velocity_pseudotime(adata)
    try:
        print("Creating velocity pseudotime plot...")
        fig, ax = plt.subplots(figsize=(10, 8))
        
        # Get UMAP coordinates and pseudotime values
        umap_coords = adata.obsm['X_umap']
        pseudotime = adata.obs['velocity_pseudotime'].values
        
        # Create scatter plot manually
        scatter = ax.scatter(umap_coords[:, 0], umap_coords[:, 1], 
                           c=pseudotime, cmap='gnuplot', s=5, alpha=0.8)
        
        ax.set_xlim(umap_min_x, umap_max_x)
        ax.set_ylim(umap_min_y, umap_max_y)
        ax.set_xlabel('UMAP_1')
        ax.set_ylabel('UMAP_2')
        ax.set_title('Full Dataset - Velocity Pseudotime')
        
        # Add colorbar
        cbar = plt.colorbar(scatter, ax=ax)
        cbar.set_label('velocity_pseudotime')
        
        plt.tight_layout()
        plt.savefig('full_dataset_pseudotime.pdf', dpi=300, bbox_inches='tight')
        plt.close()
        print("Saved: full_dataset_pseudotime.pdf")
    except Exception as e:
        print(f"Error creating pseudotime plot: {e}")
        plt.close()
    
    # 9. Top genes driving velocity
    if len(adata.obs[cell_type_col].unique()) > 1:
        try:
            print("\n--- Ranking velocity genes ---")
            scv.tl.rank_velocity_genes(adata, groupby=cell_type_col, min_corr=.3)
            
            df = pd.DataFrame(adata.uns['rank_velocity_genes']['names'])
            print(f"\nTop velocity genes by {cell_type_col}:")
            print(df.head(10))
            
            # Save top genes
            df.to_csv('full_dataset_top_velocity_genes.csv')
            print("Saved: full_dataset_top_velocity_genes.csv")
        except Exception as e:
            print(f"Error ranking velocity genes: {e}")
    
    # Save results
    adata.write('full_dataset_velocity.h5ad')
    print(f"\nSaved results to full_dataset_velocity.h5ad")
    
    return adata


# Main execution
if __name__ == "__main__":
    
    print("Starting scVelo trajectory analysis on FULL DATASET...")
    print(f"Loom files directory: {LOOM_DIR}\n")
    
    # LOAD ALL LOOM FILES ONCE
    adata_all = load_all_loom_files()
    
    print("\n" + "="*60)
    print("LOOM FILES LOADED - NOW PROCESSING FULL DATASET")
    print("="*60)
    
    # Run velocity on full dataset
    full_adata = run_velocity_full_dataset(
        adata_all=adata_all,
        metadata_file='full_placenta_metadata.csv',
        cell_type_col='cell_type_v1'
    )
    
    print("\n" + "="*60)
    print("ANALYSIS COMPLETE!")
    print("="*60)