# normalize_ihs.py
import pandas as pd
import numpy as np
import os
import zipfile
import gzip
from glob import glob
import logging
import traceback

# --- Constants ---
NUM_DAF_BINS = 50
NUM_RECOMB_BINS = 10
MIN_SNPS_PER_BIN = 2
COMPRESS_LEVEL_GZIP = 6
DEFAULT_FILL_VALUE = 0.0

# --- Column Names (Reference Only) ---
HAPBIN_EXPECTED_CONTENT_ORDER = [ # Reflects content, not exact header string
    'Index', 'ID', 'Freq', 'iHH_0', 'iHH_1', 'iHS', 'Std_iHS' # Using underscore for internal ref
]
SELSCAN_IHS_COLUMN_NAMES = [
    'ID', 'POS', 'Freq', 'iHH_1', 'iHH_0', 'iHS'
]

# --- Helper Functions ---
# (load_map_file, calculate_global_norm, JV500bin_norm - remain the same)
def load_map_file(map_filepath):
    """Loads .map file to get cM and POS values."""
    logging.info(f"    Loading map file: {os.path.basename(map_filepath)}")
    try:
        df = pd.read_csv(
            map_filepath, sep=r'\s+', header=None, usecols=[1, 2, 3],
            names=['ID', 'cM', 'POS'],
            dtype={'ID': 'string', 'cM': 'float32', 'POS': 'uint32'},
            engine='c', on_bad_lines='warn'
        )
        df['cM'] = pd.to_numeric(df['cM'], errors='coerce').fillna(0.0).clip(lower=0.0)
        df['POS'] = pd.to_numeric(df['POS'], errors='coerce')
        df.dropna(subset=['ID', 'POS'], inplace=True)
        if not df.empty: df['POS'] = df['POS'].astype('uint32')
        logging.info(f"      Loaded {len(df):,} SNPs from map file.")
        return df[['ID', 'cM', 'POS']]
    except FileNotFoundError:
        logging.error(f"    Map file not found: {map_filepath}")
        return pd.DataFrame(columns=['ID', 'cM', 'POS']).astype({'ID': 'string', 'cM': 'float32', 'POS': 'uint32'})
    except Exception as e:
        logging.error(f"    Error loading map file {map_filepath}: {e}", exc_info=True)
        return pd.DataFrame(columns=['ID', 'cM', 'POS']).astype({'ID': 'string', 'cM': 'float32', 'POS': 'uint32'})

def calculate_global_norm(full_df):
    """Genome-wide DAF bin normalization (Robust version)"""
    logging.info("  Calculating global DAF normalization...")
    full_df['iHS_numeric'] = pd.to_numeric(full_df['iHS'], errors='coerce')
    valid_ihs = full_df['iHS_numeric'].dropna()
    if valid_ihs.empty:
         logging.warning("  ⚠️ No valid iHS values found for global normalization.")
         full_df['GW_DAF_bin_iHS'] = DEFAULT_FILL_VALUE
         return full_df.drop(columns=['iHS_numeric'], errors='ignore')
    overall_mean_ihs = valid_ihs.mean()
    overall_std_ihs = valid_ihs.std(ddof=1)
    logging.info(f"    Overall iHS: mean={overall_mean_ihs:.4f}, std={overall_std_ihs:.4f}")
    if pd.isna(overall_std_ihs) or overall_std_ihs <= 0:
        logging.warning("  ⚠️ Overall iHS standard deviation is zero, NaN, or negative. Cannot perform DAF bin normalization.")
        full_df['GW_DAF_bin_iHS'] = DEFAULT_FILL_VALUE
        return full_df.drop(columns=['iHS_numeric'], errors='ignore')
    bins = np.linspace(0, 1, NUM_DAF_BINS + 1)
    full_df['DAF_bin_global'] = pd.cut(full_df['Freq'], bins, right=False, include_lowest=True, labels=False)
    global_stats = full_df.dropna(subset=['iHS_numeric']).groupby('DAF_bin_global', observed=False)['iHS_numeric'].agg(
        global_mean=lambda x: np.mean(x) if not x.empty else np.nan,
        global_std=lambda x: np.std(x, ddof=1) if len(x) > 1 else np.nan,
        count='size'
    ).reset_index()
    global_stats['adj_global_mean'] = np.where(
        (global_stats['count'] >= MIN_SNPS_PER_BIN) & (~global_stats['global_mean'].isna()),
        global_stats['global_mean'], overall_mean_ihs)
    global_stats['adj_global_std'] = np.where(
        (global_stats['count'] >= MIN_SNPS_PER_BIN) & (~global_stats['global_std'].isna()) & (global_stats['global_std'] > 0),
        global_stats['global_std'], overall_std_ihs)
    global_stats['adj_global_std'] = global_stats['adj_global_std'].fillna(overall_std_ihs).clip(lower=np.finfo(float).eps)
    merged = pd.merge(full_df, global_stats[['DAF_bin_global', 'adj_global_mean', 'adj_global_std']], on='DAF_bin_global', how='left')
    merged['adj_global_mean'] = merged['adj_global_mean'].fillna(overall_mean_ihs)
    merged['adj_global_std'] = merged['adj_global_std'].fillna(overall_std_ihs).clip(lower=np.finfo(float).eps)
    merged['GW_DAF_bin_iHS'] = ((merged['iHS_numeric'] - merged['adj_global_mean']) / merged['adj_global_std']).fillna(DEFAULT_FILL_VALUE)
    logging.info("  Finished global DAF normalization.")
    return merged.drop(columns=['DAF_bin_global', 'adj_global_mean', 'adj_global_std', 'iHS_numeric'], errors='ignore')

def JV500bin_norm(full_df):
    """DAF + Recombination Rate bin normalization (Robust version)"""
    logging.info("  Calculating DAF + Recombination Rate normalization...")
    full_df['iHS_numeric'] = pd.to_numeric(full_df['iHS'], errors='coerce')
    full_df['cM'] = pd.to_numeric(full_df['cM'], errors='coerce').fillna(0.0).clip(lower=0.0)
    full_df['Freq'] = pd.to_numeric(full_df['Freq'], errors='coerce')
    valid_ihs = full_df['iHS_numeric'].dropna()
    if valid_ihs.empty:
        logging.warning("  ⚠️ No valid iHS values found for JV500bin normalization.")
        full_df['JV_rmap_adj_norm_iHS'] = DEFAULT_FILL_VALUE
        return full_df.drop(columns=['iHS_numeric'], errors='ignore')
    global_mean = valid_ihs.mean()
    global_std = valid_ihs.std(ddof=1)
    if pd.isna(global_std) or global_std <= 0:
        logging.warning("  ⚠️ Overall iHS standard deviation is zero, NaN, or negative. Cannot perform JV500bin normalization.")
        full_df['JV_rmap_adj_norm_iHS'] = DEFAULT_FILL_VALUE
        return full_df.drop(columns=['iHS_numeric'], errors='ignore')
    daf_bins = np.linspace(0, 1, NUM_DAF_BINS + 1)
    full_df['DAF_bin'] = pd.cut(full_df['Freq'], bins=daf_bins, right=False, include_lowest=True, labels=False)
    full_df['DAF_bin'] = 'D' + full_df['DAF_bin'].astype(str)
    min_cM, max_cM = full_df['cM'].min(), full_df['cM'].max()
    if pd.isna(min_cM) or pd.isna(max_cM) or min_cM >= max_cM:
        logging.warning(f"  ⚠️ Recombination rate range invalid (min={min_cM}, max={max_cM}). Using single recombination bin.")
        full_df['recomb_bin'] = 'R_single'
    else:
        recomb_bins = np.linspace(min_cM, max_cM, NUM_RECOMB_BINS + 1)
        recomb_bins = np.unique(recomb_bins)
        if len(recomb_bins) < 2:
            logging.warning(f"  ⚠️ Could not create unique recombination bins. Using single recombination bin.")
            full_df['recomb_bin'] = 'R_single'
        else:
            full_df['recomb_bin'] = pd.cut(full_df['cM'], bins=recomb_bins, right=False, include_lowest=True, labels=False, duplicates='drop')
            full_df['recomb_bin'] = 'R' + full_df['recomb_bin'].astype(str)
    full_df['DAF_bin'] = full_df['DAF_bin'].fillna('D_NaN')
    full_df['recomb_bin'] = full_df['recomb_bin'].fillna('R_NaN')
    full_df['combined_bin'] = full_df['DAF_bin'] + '_' + full_df['recomb_bin']
    bin_stats = full_df.dropna(subset=['iHS_numeric']).groupby('combined_bin', observed=False)['iHS_numeric'].agg(
        mean_adj=lambda x: np.mean(x) if not x.empty else np.nan,
        std_adj=lambda x: np.std(x, ddof=1) if len(x) > 1 else np.nan,
        count='size'
    ).reset_index()
    bin_stats['mean_adj'] = np.where((bin_stats['count'] >= MIN_SNPS_PER_BIN) & (~bin_stats['mean_adj'].isna()), bin_stats['mean_adj'], global_mean)
    bin_stats['std_adj'] = np.where((bin_stats['count'] >= MIN_SNPS_PER_BIN) & (~bin_stats['std_adj'].isna()) & (bin_stats['std_adj'] > 0), bin_stats['std_adj'], global_std)
    bin_stats['mean_adj'] = bin_stats['mean_adj'].fillna(global_mean)
    bin_stats['std_adj'] = bin_stats['std_adj'].fillna(global_std).clip(lower=np.finfo(float).eps)
    merged = pd.merge(full_df, bin_stats[['combined_bin', 'mean_adj', 'std_adj']], on='combined_bin', how='left')
    merged['mean_adj'] = merged['mean_adj'].fillna(global_mean)
    merged['std_adj'] = merged['std_adj'].fillna(global_std).clip(lower=np.finfo(float).eps)
    merged['JV_rmap_adj_norm_iHS'] = ((merged['iHS_numeric'] - merged['mean_adj']) / merged['std_adj']).fillna(DEFAULT_FILL_VALUE)
    logging.info("  Finished DAF + Recombination Rate normalization.")
    return merged.drop(columns=['DAF_bin', 'recomb_bin', 'combined_bin', 'mean_adj', 'std_adj', 'iHS_numeric'], errors='ignore')


# === Main Normalization Function ===
def run_normalization(pop, aaref, rmap, software, alt_na, maf, chromosomes_to_process, base_ihs_dir, base_map_dir, output_dir, output_prefix):
    """
    Main function for Step 4: Normalizing iHS results.
    Processes only the chromosomes specified in chromosomes_to_process.
    """
    logging.info(f"Running Step 4: Normalize iHS Results for Pop={pop}, Software={software}, MAF=0.{maf}, AARef={aaref}, RMap={rmap}, Chromosomes={chromosomes_to_process}")
    os.makedirs(output_dir, exist_ok=True)

    all_chr_dfs = []
    files_found = False
    overall_success = True

    # --- Loop through specified chromosomes ---
    for chrom_num in chromosomes_to_process:
        chrom = str(chrom_num)
        logging.info(f"  Processing chromosome {chrom}...")

        # --- Determine expected raw IHS input filename ---
        raw_ihs_filepath = None
        selscan_expected_cols = len(SELSCAN_IHS_COLUMN_NAMES)
        hapbin_expected_cols = len(HAPBIN_EXPECTED_CONTENT_ORDER) + 1 # +1 because 'Std iHS' gets split
        header_ihs = None

        if software == 'selscan':
            raw_ihs_filename = f"{pop}.{software}.{alt_na}.{aaref}.{rmap}.maf{maf}.chr{chrom}.ihs.out"
            raw_ihs_filepath = os.path.join(base_ihs_dir, raw_ihs_filename)
            header_ihs = None # No header
        elif software == 'hapbin':
            maf_suffix_norm = f"maf{int(maf):d}"
            raw_ihs_filename = f"{pop}.{chrom}.{maf_suffix_norm}.ihs"
            raw_ihs_filepath = os.path.join(base_ihs_dir, raw_ihs_filename)
            header_ihs = 0 # Has header
        else:
            logging.error(f"Unsupported software '{software}' for normalization.")
            raise ValueError(f"Unsupported software '{software}'")

        # --- Load Raw IHS Data ---
        if not os.path.exists(raw_ihs_filepath):
            logging.warning(f"    Raw iHS file not found for Chr {chrom}: {raw_ihs_filepath}. Skipping chromosome.")
            continue

        try:
            logging.info(f"    Reading raw iHS file: {os.path.basename(raw_ihs_filepath)}")
            # Read without assigning names initially if header exists
            ihs_df = pd.read_csv(
                raw_ihs_filepath,
                sep=r'\s+',
                header=header_ihs,
                dtype={'ID': 'string'}, # Base dtype, others inferred or set later
                engine='c', comment='#', on_bad_lines='warn'
            )
            files_found = True
            logging.info(f"      Read {len(ihs_df)} records. Columns: {list(ihs_df.columns)}")
            if ihs_df.empty:
                logging.warning(f"      Raw iHS data frame is empty for Chr {chrom}. Skipping.")
                continue

            # --- Data Cleaning/Prep & Column Renaming ---
            if software == 'hapbin':
                # Check number of columns read vs expected based on header splitting
                if len(ihs_df.columns) < hapbin_expected_cols: # Expecting 8 columns due to split header
                     logging.error(f"      Hapbin output {raw_ihs_filepath} has fewer columns ({len(ihs_df.columns)}) than expected ({hapbin_expected_cols}) after read. Skipping.")
                     continue

                # Assign standard names based on column *index*, handling the split 'Std iHS'
                # Index 0: Index, 1: ID, 2: Freq, 3: iHH_0, 4: iHH_1, 5: iHS, 6: 'Std', 7: 'iHS.1' (value)
                # We need: ID, Freq, iHH_0, iHH_1, iHS, iHS_chr_std
                try:
                    rename_dict = {
                        ihs_df.columns[1]: 'ID',
                        ihs_df.columns[2]: 'Freq',
                        ihs_df.columns[3]: 'iHH_0',
                        ihs_df.columns[4]: 'iHH_1',
                        ihs_df.columns[5]: 'iHS',        # Unstandardized iHS
                        ihs_df.columns[7]: 'iHS_chr_std' # Chromosome-standardized iHS
                    }
                    # Select only the necessary columns by index and rename
                    needed_indices = [1, 2, 3, 4, 5, 7]
                    ihs_df = ihs_df.iloc[:, needed_indices].rename(columns=rename_dict)
                    logging.info(f"      Renamed Hapbin columns: {list(ihs_df.columns)}")
                except IndexError:
                     logging.error(f"      Error accessing expected columns by index in Hapbin output {raw_ihs_filepath}. Skipping.")
                     continue

            elif software == 'selscan':
                 # Assign standard names as there was no header
                 if len(ihs_df.columns) == len(SELSCAN_IHS_COLUMN_NAMES):
                     ihs_df.columns = SELSCAN_IHS_COLUMN_NAMES
                 else:
                      logging.error(f"      Selscan output {raw_ihs_filepath} has unexpected number of columns ({len(ihs_df.columns)} vs {len(SELSCAN_IHS_COLUMN_NAMES)} expected). Skipping.")
                      continue

            # Ensure essential columns are numeric after name assignment/check
            essential_cols_numeric = ['Freq', 'iHS', 'iHH_0', 'iHH_1']
            if software == 'hapbin': essential_cols_numeric.append('iHS_chr_std')
            if software == 'selscan': essential_cols_numeric.append('POS') # POS is essential for selscan input

            for col in essential_cols_numeric:
                 if col in ihs_df.columns:
                     ihs_df[col] = pd.to_numeric(ihs_df[col], errors='coerce')
                 else:
                      logging.error(f"      Essential column '{col}' missing before numeric conversion in Chr {chrom}. Skipping.")
                      overall_success = False; continue # Skip chromosome if essential col missing

            # Drop rows if ID or essential numeric cols became NaN
            essential_check_cols = ['ID', 'Freq', 'iHS']
            if software == 'selscan': essential_check_cols.append('POS')
            ihs_df.dropna(subset=essential_check_cols, inplace=True)
            if ihs_df.empty:
                 logging.warning(f"      No valid data after initial cleaning/numeric conversion for Chr {chrom}. Skipping.")
                 continue

            # --- Load Corresponding Map File ---
            map_filename = f"{pop}.{chrom}.{aaref}.{rmap}.map"
            map_filepath = os.path.join(base_map_dir, map_filename)
            map_df = load_map_file(map_filepath) # Gets ID, cM, POS

            # --- Merge IHS and Map Data ---
            if map_df.empty:
                logging.warning(f"    Map data unavailable for Chr {chrom}. Merging failed. Skipping chromosome.")
                # We absolutely need POS and cM for normalization/output.
                continue
            else:
                logging.info(f"    Merging IHS ({len(ihs_df):,}) and Map ({len(map_df):,}) data...")
                # Check for ID columns before merge
                if 'ID' not in ihs_df.columns or 'ID' not in map_df.columns:
                    logging.error("    'ID' column missing in IHS or Map data. Cannot merge.")
                    continue

                # Use left merge to keep all valid IHS results
                merged_df = pd.merge(ihs_df, map_df, on='ID', how='left')
                logging.info(f"      Merged {len(merged_df):,} SNPs.")

                # Check for essential columns POST-merge (POS and cM must come from map)
                if 'POS' not in merged_df.columns or 'cM' not in merged_df.columns:
                     logging.error(f"    POS or cM column missing after merge for Chr {chrom}. Cannot proceed.")
                     continue

                # Handle missing map info (shouldn't happen with how='left' on IHS if map was loaded, but check)
                merged_df['cM'] = merged_df['cM'].fillna(0.0).clip(lower=0.0)
                merged_df['POS'] = pd.to_numeric(merged_df['POS'], errors='coerce')
                pos_na_before = merged_df['POS'].isna().sum()
                if pos_na_before > 0:
                     logging.warning(f"    {pos_na_before} SNPs missing POS after merge for Chr {chrom}. Will be dropped.")
                     merged_df.dropna(subset=['POS'], inplace=True)
                if not merged_df.empty:
                     merged_df['POS'] = merged_df['POS'].astype('uint32')
                else:
                     logging.warning(f"     DataFrame empty after dropping missing POS for Chr {chrom}. Skipping.")
                     continue

            # Add chromosome column
            merged_df['chr'] = chrom
            all_chr_dfs.append(merged_df)

        except Exception as e:
            logging.error(f"    Error processing data for Chr {chrom}: {e}", exc_info=True)
            overall_success = False
            continue

    # --- Combine data and Normalize ---
    if not files_found:
        logging.error(f"No raw iHS input files were found or processed successfully for Pop={pop}. Normalization cannot proceed.")
        raise FileNotFoundError("No raw iHS input files found.")
    if not all_chr_dfs:
        logging.error(f"No valid data collected across chromosomes for Pop={pop}. Normalization cannot proceed.")
        raise ValueError("No data available for normalization after processing chromosomes.")

    logging.info("Concatenating data from all processed chromosomes...")
    full_df = pd.concat(all_chr_dfs, ignore_index=True)
    logging.info(f"Total SNPs across all processed chromosomes: {len(full_df):,}")
    if full_df.empty:
        logging.error("Concatenated DataFrame is empty. Normalization cannot proceed.")
        raise ValueError("Concatenated DataFrame is empty.")

    # --- Data Cleaning before Normalization ---
    initial_rows = len(full_df)
    logging.info("Cleaning combined data before normalization...")
    critical_cols = ['POS', 'Freq', 'iHS', 'cM']
    if not all(col in full_df.columns for col in critical_cols):
        missing_crit = [col for col in critical_cols if col not in full_df.columns]
        logging.error(f"Critical columns missing in combined DataFrame: {missing_crit}. Aborting.")
        raise KeyError(f"Missing critical columns in combined data: {missing_crit}")
    for col in critical_cols: full_df[col] = pd.to_numeric(full_df[col], errors='coerce')
    full_df['cM'] = full_df['cM'].fillna(0.0).clip(lower=0.0)
    full_df.dropna(subset=critical_cols, inplace=True)
    rows_after_dropna = len(full_df)
    if rows_after_dropna < initial_rows: logging.info(f"  Dropped {initial_rows - rows_after_dropna:,} rows with NaN in critical columns.")
    if full_df.empty:
        logging.error("No valid data remaining after cleaning combined data. Normalization cannot proceed.")
        raise ValueError("No valid data remaining after cleaning combined data.")

    # --- Perform Normalizations ---
    full_df = calculate_global_norm(full_df)
    full_df = JV500bin_norm(full_df)

    # --- Prepare Final Output Columns ---
    logging.info("Preparing final output columns...")
    # Define base columns expected after merge/prep
    base_cols = ['chr', 'POS', 'ID', 'Freq', 'iHH_0', 'iHH_1', 'iHS'] # Use standard internal names
    norm_cols = ['GW_DAF_bin_iHS', 'JV_rmap_adj_norm_iHS']
    optional_cols = []
    if software == 'hapbin' and 'iHS_chr_std' in full_df.columns: # Check internal name
        optional_cols.append('iHS_chr_std')
    final_expected_cols = base_cols + optional_cols + norm_cols
    output_df = full_df.copy()
    missing_final_cols = [col for col in final_expected_cols if col not in output_df.columns]
    if missing_final_cols:
        logging.warning(f"⚠️ Final output columns missing: {missing_final_cols}. Filling with {DEFAULT_FILL_VALUE}.")
        for col in missing_final_cols: output_df[col] = DEFAULT_FILL_VALUE
    # Rename columns for final output
    rename_map = {
        'chr': 'CHR', 'ID': 'rsID', 'Freq': 'DAF', 'iHH_0': 'iHH0', 'iHH_1': 'iHH1',
        'iHS': 'unstd_iHS', 'iHS_chr_std': 'chr_std_iHS', # Rename internal name to output name
        'GW_DAF_bin_iHS': 'gw_std_iHS', 'JV_rmap_adj_norm_iHS': 'std_iHS'}
    output_df = output_df[[col for col in final_expected_cols if col in output_df.columns]].rename(columns=rename_map)
    # Define standard output order
    standard_order = ['CHR', 'POS', 'rsID', 'DAF', 'iHH0', 'iHH1', 'unstd_iHS']
    if 'chr_std_iHS' in output_df.columns: standard_order.append('chr_std_iHS')
    standard_order.extend(['gw_std_iHS', 'std_iHS'])
    output_df = output_df[[col for col in standard_order if col in output_df.columns]]

    # Sort final output
    logging.info("Sorting final output by CHR and POS...")
    output_df['CHR'] = pd.to_numeric(output_df['CHR'], errors='coerce')
    output_df.dropna(subset=['CHR','POS'], inplace=True)
    output_df['CHR'] = output_df['CHR'].astype(int)
    output_df['POS'] = pd.to_numeric(output_df['POS'], errors='coerce').astype('Int64')
    output_df.dropna(subset=['POS'], inplace=True)
    output_df['POS'] = output_df['POS'].astype(int)
    output_df = output_df.sort_values(by=['CHR', 'POS'])
    output_df['CHR'] = output_df['CHR'].astype(str)
    output_df['POS'] = output_df['POS'].astype(str)

    # --- Save Final Compressed Output ---
    final_output_path = os.path.join(output_dir, f"{output_prefix}.norm.ihs.gz")
    logging.info(f"Saving final normalized results to: {final_output_path}")
    try:
        output_df.to_csv(
            final_output_path, sep='\t', index=False, float_format='%.6f', encoding='utf-8', na_rep='NA',
            compression={'method': 'gzip', 'compresslevel': COMPRESS_LEVEL_GZIP}
        )
        logging.info(f"✅ Successfully saved {len(output_df):,} normalized SNPs to {final_output_path}")
    except Exception as e:
        logging.error(f"❌ Error writing final output file {final_output_path}: {e}", exc_info=True)
        raise

    if not overall_success:
        logging.warning("Normalization completed, but errors occurred during processing of some chromosomes.")
    else:
        logging.info(f"Step 4 completed successfully for {pop}.")
