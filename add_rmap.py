# add_rmap.py
import os
import sys
import glob
import re
import time
import bisect
import shutil
import pandas as pd
import logging
import numpy as np

# --- Constants ---
# Column indices (0-based) for input AA-TPED
TPED_CHR_COL = 0
TPED_RSID_COL = 1
TPED_CM_COL = 2
TPED_POS_COL = 3
TPED_GENO_START_COL = 4

# Column indices for different map formats (used by loaders)
# PHV map format (chr pos rate cM)
MAP_CHR_COL_PHV = 0
MAP_POS_COL_PHV = 1
MAP_RATE_COL_PHV = 2
MAP_CM_COL_PHV = 3
# PYRHO/JV/TIAN map format (pos chr cM or similar)
MAP_POS_COL_OTHER = 0
MAP_CHR_COL_OTHER = 1
MAP_CM_COL_OTHER = 2


# === Helper Functions for Loading Recombination Maps ===

# --- PHV Map Loader (for interpolation) ---
def load_rmap_phv(map_file_path):
    """
    Loads genetic map data (position, cM) from the PHV format file.
    Expects format: chr pos rate cM (with header)
    Returns: list of tuples [(position (int), cM (float))] sorted by position, or None.
    """
    logging.info(f"  Loading PHV recombination map from: {os.path.basename(map_file_path)}...")
    map_data = []
    line_count = 0
    loaded_count = 0
    skipped_malformed = 0
    try:
        with open(map_file_path, 'r', encoding='utf-8') as infile:
            header_line = infile.readline()
            if not header_line:
                logging.error(f"    PHV map file is empty: {os.path.basename(map_file_path)}")
                return None
            for line in infile:
                line_count += 1
                try:
                    stripped_line = line.strip()
                    if not stripped_line: continue
                    fields = re.split(r'\s+', stripped_line)
                    if len(fields) <= max(MAP_POS_COL_PHV, MAP_CM_COL_PHV):
                        skipped_malformed += 1
                        continue
                    pos_str = fields[MAP_POS_COL_PHV]
                    cm_str = fields[MAP_CM_COL_PHV]
                    if not pos_str.isdigit():
                        skipped_malformed += 1
                        continue
                    try: cm_val = float(cm_str)
                    except ValueError:
                        skipped_malformed += 1
                        continue
                    map_data.append((int(pos_str), cm_val))
                    loaded_count += 1
                except (ValueError, IndexError) as e_inner:
                    logging.warning(f"    Warning: Error processing PHV map line {line_count + 1}: {e_inner}. Line: '{stripped_line}'")
                    skipped_malformed += 1
                    continue
        if skipped_malformed > 0:
             logging.warning(f"    Skipped {skipped_malformed:,} malformed lines during PHV map loading.")
        if not map_data:
             logging.error(f"    Error: No valid map data loaded from {os.path.basename(map_file_path)}.")
             return None
        map_data.sort(key=lambda x: x[0]) # Sort
        # Deduplicate positions
        seen_pos = set()
        dedup_map_data = []
        duplicates_found = False
        for pos, cm in map_data:
            if pos not in seen_pos:
                dedup_map_data.append((pos, cm))
                seen_pos.add(pos)
            else: duplicates_found = True
        if duplicates_found:
             original_count = len(map_data)
             map_data = dedup_map_data
             logging.warning(f"    Removed {original_count - len(map_data):,} duplicate positions from PHV map data.")
        if not map_data:
             logging.error(f"    Error: No map data remaining after deduplication for {os.path.basename(map_file_path)}.")
             return None
        if len(map_data) < 2:
            logging.warning(f"    Warning: Fewer than 2 unique map points ({len(map_data)}) loaded from {os.path.basename(map_file_path)}. Interpolation might not be effective.")
        logging.info(f"    Successfully loaded and prepared {len(map_data):,} unique PHV map points.")
        return map_data # List of (int_pos, float_cM) tuples
    except FileNotFoundError:
        logging.error(f"    PHV map file not found: {map_file_path}")
        return None
    except Exception as e:
        logging.error(f"    Error reading PHV map file {map_file_path}: {e}", exc_info=True)
        return None

# --- PYRHO/JV/TIAN Map Loader (for last point lookup) ---
def load_rmap_other(rmap_file, rmap_type):
    """
    Reads the recombination map file (PYRHO/JV/TIAN format: pos chr cM or similar).
    Returns positions (list of ints) and cm_values (list of strings). Handles header.
    """
    logging.info(f"  Loading {rmap_type.upper()} recombination map from: {os.path.basename(rmap_file)}...")
    positions = []
    cm_values = [] # Store cM as string
    try:
        with open(rmap_file, 'r', encoding='utf-8') as f:
            first_line = f.readline().strip()
            is_header = any(c.isalpha() for c in first_line)
            if not is_header and first_line:
                 logging.info(f"    Info: Assuming no header found in {os.path.basename(rmap_file)}. Processing first line as data.")
                 f.seek(0)
            elif is_header:
                 logging.info(f"    Info: Detected header in {os.path.basename(rmap_file)}: '{first_line}'. Skipping.")
            line_count = 1 if is_header else 0
            skipped_malformed = 0
            loaded_count = 0
            for line in f:
                line_count += 1
                stripped_line = line.strip()
                if not stripped_line: continue
                try:
                    parts = re.split(r'\s+', stripped_line)
                    # Assume format pos chr cM for these types
                    if len(parts) < 3:
                         skipped_malformed += 1
                         continue
                    pos_str = parts[MAP_POS_COL_OTHER]
                    cm_str = parts[MAP_CM_COL_OTHER]
                    if not pos_str.isdigit():
                        skipped_malformed += 1
                        continue
                    try: float(cm_str) # Check cM is numeric
                    except ValueError:
                         skipped_malformed += 1
                         continue
                    positions.append(int(pos_str))
                    cm_values.append(cm_str) # Keep as string
                    loaded_count += 1
                except (ValueError, IndexError) as e:
                    logging.warning(f"    Warning: Skipping malformed line {line_count} in {rmap_type.upper()} map file {os.path.basename(rmap_file)}: {line.strip()} - Error: {e}")
                    skipped_malformed += 1
                    continue
            if skipped_malformed > 0:
                logging.warning(f"    Warning: Skipped {skipped_malformed:,} malformed lines in {os.path.basename(rmap_file)}.")
            if loaded_count == 0:
                 logging.error(f"    Error: No valid data loaded from recombination map {os.path.basename(rmap_file)}.")
                 return None, None
            logging.info(f"    Successfully loaded {loaded_count:,} {rmap_type.upper()} map points.")
            # Return positions (int list) and cm_values (string list)
            return positions, cm_values
    except FileNotFoundError:
        logging.warning(f"    Recombination map file not found: {rmap_file}. Will use cM=0.")
        return None, None
    except Exception as e:
        logging.error(f"    Error reading recombination map file {rmap_file}: {e}", exc_info=True)
        return None, None


# === Helper Functions for Processing AA-TPED data ===

# --- Interpolation function (used only for rmap='phv') ---
# (Function interpolate_cm remains the same as previous version)
def interpolate_cm(target_pos, map_data):
    """
    Interpolates genetic map position (cM) for a target base pair position.
    map_data: List of (position, cM) tuples, SORTED BY POSITION.
    Returns interpolated cM (float) or None if interpolation is not possible.
    """
    if not map_data or len(map_data) < 2:
        return None # Cannot interpolate with fewer than 2 points

    positions = [p for p, c in map_data]
    cms = [c for p, c in map_data]

    if target_pos < positions[0] or target_pos > positions[-1]:
        return None
    if target_pos == positions[-1]:
        return cms[-1]

    idx = bisect.bisect_left(positions, target_pos)

    if idx < len(positions) and positions[idx] == target_pos:
        return cms[idx]
    if idx == 0:
         logging.warning(f"Interpolation logic error for pos {target_pos}. idx=0 unexpectedly.")
         return None

    pos_lower = positions[idx-1]
    pos_upper = positions[idx]
    cm_lower = cms[idx-1]
    cm_upper = cms[idx]

    if pos_upper == pos_lower:
         logging.warning(f"    Warning: Found duplicate positions ({pos_lower}) in map data during interpolation. Using lower cM value.")
         return cm_lower

    fraction = (target_pos - pos_lower) / (pos_upper - pos_lower)
    interpolated_cM = cm_lower + fraction * (cm_upper - cm_lower)
    return interpolated_cM


# --- Core Processing Function (handles both interpolation and lookup) ---
# (Function process_aa_tped_with_rmap remains largely the same, ensures robust column handling)
def process_aa_tped_with_rmap(input_aa_tped_path, rmap_type, map_data_or_positions, map_cm_values, expected_chrom_str):
    """
    Reads AA-TPED, updates cM based on rmap_type, filters, sorts, and returns DataFrame.
    Returns processed DataFrame or None on error. Also returns counts.
    """
    logging.info(f"  Processing AA-TPED data from: {os.path.basename(input_aa_tped_path)} using rmap type '{rmap_type}'")
    variants_read_count = 0
    variants_after_chrom_filter_count = 0
    variants_after_cm_update_count = 0

    try:
        # Read AA-TPED using pandas
        map_df = pd.read_csv(input_aa_tped_path, delim_whitespace=True, header=None, dtype=object)
        variants_read_count = len(map_df)

        num_cols = map_df.shape[1]
        if num_cols < TPED_GENO_START_COL + 1:
             logging.error(f"Error: AA-TPED file {os.path.basename(input_aa_tped_path)} has < {TPED_GENO_START_COL + 1} columns.")
             return None, variants_read_count, 0, 0
        col_names = ['chr', 'rsid', 'cM_old', 'pos'] + [f'geno_{i}' for i in range(num_cols - TPED_GENO_START_COL)]
        map_df.columns = col_names

        # Convert position to integer
        map_df['pos'] = pd.to_numeric(map_df['pos'], errors='coerce')
        initial_rows_pos = len(map_df)
        map_df.dropna(subset=['pos'], inplace=True)
        map_df['pos'] = map_df['pos'].astype('Int64') # Nullable integer
        if len(map_df) < initial_rows_pos:
            logging.warning(f"    Dropped {initial_rows_pos - len(map_df)} rows with invalid 'pos' values.")
            variants_read_count = len(map_df) # Update count
        if map_df.empty: return None, variants_read_count, 0, 0 # Stop if no valid positions

        # --- Chromosome Filtering ---
        original_variant_count_before_chrom = len(map_df)
        map_df['chr_norm'] = map_df['chr'].astype(str).str.upper().str.replace('CHR', '')
        expected_chrom_norm = expected_chrom_str.upper().replace('CHR','')
        map_df = map_df[map_df['chr_norm'] == expected_chrom_norm]
        variants_after_chrom_filter_count = len(map_df)
        if variants_after_chrom_filter_count < original_variant_count_before_chrom:
            logging.info(f"    Filtered out {original_variant_count_before_chrom - variants_after_chrom_filter_count:,} variants due to chromosome mismatch (expected chr {expected_chrom_str}).")
        if variants_after_chrom_filter_count == 0:
            logging.warning(f"    Warning: No variants remaining after chromosome filtering for chr {expected_chrom_str}.")
            return None, variants_read_count, 0, 0
        map_df = map_df.drop(columns=['chr_norm'])

        # --- cM Update Logic ---
        map_available = (rmap_type == 'phv' and map_data_or_positions is not None and len(map_data_or_positions)>0) or \
                        (rmap_type != 'phv' and map_data_or_positions is not None and map_cm_values is not None and len(map_data_or_positions)>0)

        if not map_available:
            logging.warning(f"    Recombination map data unavailable/empty for rmap type '{rmap_type}'. Assigning cM = 0.0.")
            map_df['cM'] = "0.00000000"
            variants_after_cm_update_count = len(map_df)
        elif rmap_type == 'phv':
            # --- Interpolation Logic (PHV) ---
            logging.info(f"    Interpolating cM values using PHV map...")
            genetic_map_data = map_data_or_positions
            map_df['cM_new_float'] = map_df['pos'].apply(lambda p: interpolate_cm(p, genetic_map_data))
            original_count_before_cm_filter = len(map_df)
            map_df = map_df.dropna(subset=['cM_new_float'])
            variants_after_cm_update_count = len(map_df)
            skipped_cm_interp = original_count_before_cm_filter - variants_after_cm_update_count
            if skipped_cm_interp > 0: logging.info(f"    Filtered out {skipped_cm_interp:,} variants where cM interpolation failed.")
            if variants_after_cm_update_count == 0:
                 logging.warning(f"    Warning: No variants remaining after cM interpolation filtering.")
                 return None, variants_read_count, variants_after_chrom_filter_count, 0
            map_df['cM'] = map_df['cM_new_float'].apply(lambda cm: f"{cm:.8f}")
            map_df = map_df.drop(columns=['cM_new_float'])
        else:
            # --- Last Point Lookup Logic (PYRHO/JV/TIAN) ---
            logging.info(f"    Assigning cM values using {rmap_type.upper()} map (method: last point <= pos)...")
            rmap_positions = map_data_or_positions
            # Ensure map is sorted
            rmap_combined = list(zip(rmap_positions, map_cm_values))
            if not all(rmap_combined[i][0] <= rmap_combined[i+1][0] for i in range(len(rmap_combined)-1)):
                logging.warning(f"    Map positions for {rmap_type.upper()} were not sorted. Sorting now.")
                try:
                    rmap_combined.sort(key=lambda item: item[0])
                    rmap_positions = [item[0] for item in rmap_combined]
                    map_cm_values = [item[1] for item in rmap_combined]
                except Exception as sort_e:
                     logging.error(f"    Error sorting {rmap_type.upper()} map positions: {sort_e}. Assigning cM = 0.")
                     map_df['cM'] = "0"
                     variants_after_cm_update_count = len(map_df)
                     rmap_positions = None # Prevent further use
            if rmap_positions: # Proceed if map data is valid
                new_cm_list = []
                for pos_val in map_df['pos']: # Use the integer position
                    idx = bisect.bisect_right(rmap_positions, pos_val)
                    if idx > 0: new_cm_list.append(map_cm_values[idx-1])
                    else: new_cm_list.append("0")
                map_df['cM'] = new_cm_list
                variants_after_cm_update_count = len(map_df)
                logging.info(f"    Assigned cM for {variants_after_cm_update_count:,} variants.")

        # --- Genotype Filtering ('0' or '1' only) ---
        logging.info(f"  Filtering genotypes (keeping only 0/1)...")
        genotype_cols = [col for col in map_df.columns if col.startswith('geno_')]
        if not genotype_cols:
            logging.error("Error: No genotype columns found after processing.")
            return None, variants_read_count, variants_after_chrom_filter_count, variants_after_cm_update_count

        def check_genotypes(row):
            return row.isin(['0', '1']).all()

        valid_genotype_mask = map_df[genotype_cols].apply(check_genotypes, axis=1)
        original_count_before_geno_filter = len(map_df)
        map_df = map_df[valid_genotype_mask]
        variants_after_genotype_filter = len(map_df)
        skipped_genotype_count = original_count_before_geno_filter - variants_after_genotype_filter
        if skipped_genotype_count > 0: logging.info(f"    Filtered out {skipped_genotype_count:,} variants with invalid genotypes.")
        logging.info(f"    Kept {variants_after_genotype_filter:,} variants after genotype filtering.")
        if variants_after_genotype_filter == 0:
             logging.warning(f"    Warning: No variants remaining after genotype filtering.")
             return None, variants_read_count, variants_after_chrom_filter_count, variants_after_cm_update_count

        # --- Sort by Position ---
        logging.info(f"    Sorting {len(map_df):,} variants by position...")
        map_df.sort_values(by='pos', inplace=True)
        map_df.reset_index(drop=True, inplace=True)

        # --- Final Column Selection ---
        final_genotype_cols = [col for col in map_df.columns if col.startswith('geno_')]
        final_columns = ['chr', 'rsid', 'cM', 'pos'] + final_genotype_cols
        if 'cM' not in map_df.columns: map_df['cM'] = "0.00000000" # Safety
        map_df = map_df[final_columns]

        # Convert 'pos' back to string for output if needed, or keep as int? Keep as int for now.
        # map_df['pos'] = map_df['pos'].astype(str)

        return map_df, variants_read_count, variants_after_chrom_filter_count, variants_after_cm_update_count

    except FileNotFoundError:
        logging.error(f"Error: Input AA-TPED file not found: {input_aa_tped_path}")
        return None, 0, 0, 0
    except pd.errors.EmptyDataError:
        logging.warning(f"Warning: Input AA-TPED file is empty: {input_aa_tped_path}")
        return None, 0, 0, 0
    except Exception as e:
        logging.error(f"Error processing AA-TPED file {os.path.basename(input_aa_tped_path)} with pandas: {e}", exc_info=True)
        return None, 0, 0, 0


# === Main Function for Adding Recombination Map ===
def run_add_rmap(pop, aaref, rmap, chromosomes_to_process, base_aa_tped_dir, base_rmap_dir_pyrho, base_rmap_dir_jv, base_rmap_dir_phv, base_rmap_dir_tian, output_dir):
    """
    Main function for Step 2: Adding recombination map info to AA-TPED files.
    Processes only the chromosomes specified in chromosomes_to_process.
    Uses 'pyrho' for the rmap previously called 'arg'.
    """
    logging.info(f"Running Step 2: Add Recombination Map for Pop={pop}, AARef={aaref}, RMap={rmap}, Chromosomes={chromosomes_to_process}")
    os.makedirs(output_dir, exist_ok=True)

    # Select rmap directory, loading function, and filename pattern
    if rmap == 'phv':
        rmap_dir = base_rmap_dir_phv
        load_func = load_rmap_phv
        rmap_filename_func = lambda chrom: f"sexavg_chr{chrom}.txt"
        needs_cm_values_list = False # load_func returns list of tuples
    elif rmap == 'pyrho': # Changed from 'arg'
        rmap_dir = base_rmap_dir_pyrho # Use the renamed variable
        load_func = lambda f: load_rmap_other(f, rmap_type='pyrho') # Pass rmap type to loader
        rmap_filename_func = lambda chrom: f"PJL_map_chr_{chrom}.txt"
        needs_cm_values_list = True # load_func returns positions, cm_values
    elif rmap == 'jv':
        rmap_dir = base_rmap_dir_jv
        load_func = lambda f: load_rmap_other(f, rmap_type='jv')
        rmap_filename_func = lambda chrom: f"PJL.{chrom}.map.txt"
        needs_cm_values_list = True
    elif rmap == 'tian':
        rmap_dir = base_rmap_dir_tian
        load_func = lambda f: load_rmap_other(f, rmap_type='tian')
        rmap_filename_func = lambda chrom: f"PJL.{chrom}.map.txt"
        needs_cm_values_list = True
    else:
        raise ValueError(f"Invalid rmap specified: {rmap}")

    overall_success = True
    processed_files_count = 0
    # Add counters specific to this step if needed

    # --- Loop through specified chromosomes ---
    for chrom_num in chromosomes_to_process:
        chrom_str = str(chrom_num).upper()
        # Find input AA-TPED file (output of Step 1)
        aa_tped_file_pattern = os.path.join(base_aa_tped_dir, f"{pop}.{chrom_str}.aa.tped")
        aa_tped_files_found = glob.glob(aa_tped_file_pattern)

        if not aa_tped_files_found:
            logging.warning(f"Input AA-TPED file not found for chromosome {chrom_str} (pattern: {aa_tped_file_pattern}). Skipping.")
            continue
        aa_tped_file = aa_tped_files_found[0] # Assume only one match
        base_aa_tped = os.path.basename(aa_tped_file)

        logging.info(f"Processing Chr {chrom_str} from AA-TPED file: {base_aa_tped}")
        start_time_chr = time.time()

        # --- Load Recombination Map ---
        rmap_filename = rmap_filename_func(chrom_str)
        rmap_filepath = os.path.join(rmap_dir, rmap_filename)

        map_positions = None
        map_cm_values = None
        map_data_tuples = None

        if needs_cm_values_list:
            map_positions, map_cm_values = load_func(rmap_filepath)
        else:
            map_data_tuples = load_func(rmap_filepath)

        # --- Process the AA-TPED file with the loaded map ---
        processed_df, read_count, chrom_filt_count, cm_upd_count = process_aa_tped_with_rmap(
            input_aa_tped_path=aa_tped_file,
            rmap_type=rmap,
            map_data_or_positions=map_data_tuples if rmap == 'phv' else map_positions,
            map_cm_values=map_cm_values,
            expected_chrom_str=chrom_str
        )

        if processed_df is None:
            logging.error(f"  Skipping chromosome {chrom_str} due to processing errors in process_aa_tped_with_rmap.")
            overall_success = False # Mark failure if critical processing fails
            continue

        if processed_df.empty:
             logging.warning(f"  No variants remained after processing for {base_aa_tped}. No output files generated for Chr {chrom_str}.")
             continue

        processed_files_count += 1 # Increment count of successfully processed files

        # --- Define Output File Paths ---
        output_base = f"{pop}.{chrom_str}.{aaref}.{rmap}" # Base name for outputs
        final_tped_filepath = os.path.join(output_dir, f"{output_base}.tped")
        map_out_filepath = os.path.join(output_dir, f"{output_base}.map")
        hap_out_filepath = os.path.join(output_dir, f"{output_base}.hap")

        # --- Save Outputs ---
        try:
            logging.info(f"  Saving final outputs for chromosome {chrom_str}...")
            genotype_cols = [col for col in processed_df.columns if col.startswith('geno_')]
            map_cols = ['chr', 'rsid', 'cM', 'pos']

            # Save Final TPED (.tped) - space separated
            # Convert pos back to string for standard TPED if needed, or keep integer?
            # Sticking to pandas default string conversion for now.
            processed_df.to_csv(final_tped_filepath, sep=' ', header=False, index=False, lineterminator='\n', encoding='utf-8')
            logging.info(f"    Successfully wrote {len(processed_df):,} variants to {os.path.basename(final_tped_filepath)}.")

            # Save MAP file (.map) - space separated
            processed_df[map_cols].to_csv(map_out_filepath, sep=' ', header=False, index=False, lineterminator='\n', encoding='utf-8')
            logging.info(f"    Successfully wrote {len(processed_df):,} variants to {os.path.basename(map_out_filepath)}.")

            # Save HAP file (.hap) - space separated
            processed_df[genotype_cols].to_csv(hap_out_filepath, sep=' ', header=False, index=False, lineterminator='\n', encoding='utf-8')
            logging.info(f"    Successfully wrote {len(processed_df):,} variants to {os.path.basename(hap_out_filepath)}.")

            end_time_chr = time.time()
            logging.info(f"  Finished chromosome {chrom_str} processing in {end_time_chr - start_time_chr:.2f} seconds.")

        except Exception as e_save:
            logging.error(f"  Error saving output files for {base_aa_tped}: {e_save}", exc_info=True)
            overall_success = False
            for f_path in [final_tped_filepath, map_out_filepath, hap_out_filepath]:
                if os.path.exists(f_path):
                    try: os.remove(f_path)
                    except OSError: pass

    # End of loop through specified chromosomes
    logging.info("-" * 20)
    logging.info(f"Population {pop} Summary (Step 2 - {len(chromosomes_to_process)} chr requested, {processed_files_count} processed):")
    # Add more summary stats if needed

    if processed_files_count == 0 and chromosomes_to_process:
         logging.error(f"Step 2 failed: No input AA-TPED files were found or processed for the requested chromosome(s).")
         raise FileNotFoundError("No input AA-TPED files processed in Step 2.")
    if not overall_success:
        raise RuntimeError(f"Step 2 failed for population {pop}, aaref {aaref}, rmap {rmap}. Check logs.")
    elif processed_files_count < len(chromosomes_to_process):
         logging.warning(f"Step 2 completed for population {pop}, but some requested chromosomes may have been skipped due to missing input or errors.")
    else:
        logging.info(f"Step 2 completed successfully for population {pop}.")
