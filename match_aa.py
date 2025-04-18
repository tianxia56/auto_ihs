# match_aa.py
import os
import sys
import glob
import re
import time
import gzip
import logging

# --- Constants ---
TPED_MISSING_ALLELE = '0' # Standard missing allele in TPED format

# === Helper Functions for Loading Ancestral Alleles ===
# (Loaders remain the same as in the previous version - load_ancestral_alleles_arg and load_ancestral_alleles_jv_phv)
# It's assumed load_ancestral_alleles_arg is used for aaref='arg' (even if rmap is 'pyrho')

# --- ARG Ancestral Allele Loader ---
def load_ancestral_alleles_arg(aa_ref_file_gz):
    """
    Loads ancestral alleles from the gzipped ErinG TSV format file (ARG).
    Assumes columns: CHR, POS, ErinG_ARG_AA, ... (tab-separated)
    Uses POS (column index 1) and ErinG_ARG_AA (column index 2).
    Returns: dict { position (int): allele (str) } or None on failure.
    Stores alleles as uppercase.
    """
    logging.info(f"  Loading ARG ancestral alleles from: {os.path.basename(aa_ref_file_gz)}...")
    aa_map = {}
    line_count = 0
    loaded_count = 0
    skipped_malformed = 0
    try:
        with gzip.open(aa_ref_file_gz, 'rt', encoding='utf-8') as infile: # Specify encoding
            header_line = infile.readline()
            if not header_line:
                logging.error(f"    Ancestral allele file is empty: {os.path.basename(aa_ref_file_gz)}")
                return None

            header = header_line.strip().split('\t')
            expected_header_start = ['CHR', 'POS', 'ErinG_ARG_AA']
            if len(header) < len(expected_header_start) or header[:len(expected_header_start)] != expected_header_start:
                 logging.error(f"    Unexpected header start in {os.path.basename(aa_ref_file_gz)}. Expected: {expected_header_start}, Got: {header}")
                 return None

            for line in infile:
                line_count += 1
                try:
                    stripped_line = line.strip()
                    if not stripped_line: continue
                    fields = stripped_line.split('\t')
                    if len(fields) < 3:
                        skipped_malformed += 1
                        continue

                    pos_str = fields[1]
                    allele = fields[2]
                    allele_upper = allele.upper()

                    if not pos_str.isdigit() or len(allele_upper) != 1 or allele_upper not in 'ACGT':
                         skipped_malformed += 1
                         continue
                    aa_map[int(pos_str)] = allele_upper
                    loaded_count +=1
                except (ValueError, IndexError) as e:
                    logging.warning(f"    Skipping invalid line format {line_count+1} in {os.path.basename(aa_ref_file_gz)}: {stripped_line} - Error: {e}")
                    skipped_malformed += 1
                    continue

        if skipped_malformed > 0:
             logging.warning(f"    Skipped {skipped_malformed:,} malformed or invalid lines during ARG AA loading.")
        if loaded_count == 0:
             logging.error(f"    No valid ancestral alleles loaded from {os.path.basename(aa_ref_file_gz)}")
             return None
        logging.info(f"    Successfully loaded {loaded_count:,} ARG ancestral alleles.")
        return aa_map
    except FileNotFoundError:
        logging.error(f"    Ancestral allele file not found: {aa_ref_file_gz}")
        return None
    except Exception as e:
        logging.error(f"    Error reading ancestral allele file {aa_ref_file_gz}: {e}", exc_info=True)
        return None

# --- JV/PHV Ancestral Allele Loader (shared logic) ---
def load_ancestral_alleles_jv_phv(readable_aa_file_gz):
    """
    Loads ancestral alleles from the gzipped readable format file (JV/PHV).
    Assumes format: chr\tposition\tancestral_allele (tab-separated, with header)
    Returns: dict { position (int): allele (str) } or None on failure.
    Stores alleles as uppercase.
    """
    logging.info(f"  Loading JV/PHV ancestral alleles from: {os.path.basename(readable_aa_file_gz)}...")
    aa_map = {}
    line_count = 0
    loaded_count = 0
    skipped_malformed = 0
    try:
        with gzip.open(readable_aa_file_gz, 'rt', encoding='utf-8') as infile: # Specify encoding
            header_line = infile.readline()
            if not header_line:
                logging.error(f"    Ancestral allele file is empty: {os.path.basename(readable_aa_file_gz)}")
                return None

            header = header_line.strip().split('\t')
            expected_header = ['chr', 'position', 'ancestral_allele']
            if header != expected_header:
                 logging.error(f"    Unexpected header in {os.path.basename(readable_aa_file_gz)}. Expected: {expected_header}, Got: {header}")
                 return None

            for line in infile:
                line_count += 1
                try:
                    stripped_line = line.strip()
                    if not stripped_line: continue
                    fields = stripped_line.split('\t')
                    if len(fields) != 3:
                        skipped_malformed += 1
                        continue

                    chrom, pos_str, allele = fields
                    allele_upper = allele.upper()

                    # Validation
                    if not pos_str.isdigit() or len(allele_upper) != 1 or allele_upper not in 'ACGT':
                         skipped_malformed += 1
                         continue
                    aa_map[int(pos_str)] = allele_upper
                    loaded_count +=1
                except (ValueError, IndexError) as e:
                    logging.warning(f"    Skipping invalid line format {line_count+1} in {os.path.basename(readable_aa_file_gz)}: {stripped_line} - Error: {e}")
                    skipped_malformed += 1
                    continue

        if skipped_malformed > 0:
             logging.warning(f"    Skipped {skipped_malformed:,} malformed or invalid lines during JV/PHV AA loading.")
        if loaded_count == 0:
             logging.error(f"    No valid ancestral alleles loaded from {os.path.basename(readable_aa_file_gz)}")
             return None
        logging.info(f"    Successfully loaded {loaded_count:,} JV/PHV ancestral alleles.")
        return aa_map
    except FileNotFoundError:
        logging.error(f"    Ancestral allele file not found: {readable_aa_file_gz}")
        return None
    except Exception as e:
        logging.error(f"    Error reading ancestral allele file {readable_aa_file_gz}: {e}", exc_info=True)
        return None


# === Main AA Matching Function ===
def run_match_aa(pop, aaref, chromosomes_to_process, base_tped_dir, base_aa_ref_dir_arg, base_aa_ref_dir_jv, base_aa_ref_dir_phv, output_dir):
    """
    Main function for Step 1: Matching ancestral alleles.
    Processes only the chromosomes specified in chromosomes_to_process.
    """
    logging.info(f"Running Step 1: Match Ancestral Alleles for Pop={pop}, AARef={aaref}, Chromosomes={chromosomes_to_process}")
    os.makedirs(output_dir, exist_ok=True)

    # Select AA ref directory and loading function based on aaref
    if aaref == 'arg':
        aa_ref_dir = base_aa_ref_dir_arg
        load_func = load_ancestral_alleles_arg
        # ARG uses chr prefix generally, even for X/Y? Verify if needed. Assuming yes for now.
        aa_filename_func = lambda chrom: f"chr{chrom}_ErinG_AA_hg19.tsv.gz"
    elif aaref == 'jv':
        aa_ref_dir = base_aa_ref_dir_jv
        load_func = load_ancestral_alleles_jv_phv
        aa_filename_func = lambda chrom: f"human_ancestor_{chrom}_readable.txt.gz"
    elif aaref == 'phv':
        aa_ref_dir = base_aa_ref_dir_phv
        load_func = load_ancestral_alleles_jv_phv
        aa_filename_func = lambda chrom: f"human_ancestor_{chrom}_readable.txt.gz"
    else:
        raise ValueError(f"Invalid aaref specified: {aaref}")

    overall_success = True
    total_variants_processed_pop = 0
    total_variants_written_pop = 0
    total_ancestral_missing_pop = 0

    # --- Loop through specified chromosomes ---
    processed_files_count = 0
    for chrom_num in chromosomes_to_process:
        chrom_str = str(chrom_num).upper() # Ensure consistency (e.g., "22", "X")
        # Find the corresponding input TPED file
        # Handle potential 'chr' prefix in input filename but not in chrom_str
        tped_file_pattern1 = os.path.join(base_tped_dir, f"{pop}.{chrom_str}.tped")
        tped_file_pattern2 = os.path.join(base_tped_dir, f"{pop}.chr{chrom_str}.tped")
        tped_files_found = glob.glob(tped_file_pattern1) + glob.glob(tped_file_pattern2)

        if not tped_files_found:
            logging.warning(f"No input TPED file found for chromosome {chrom_str} (patterns: {tped_file_pattern1}, {tped_file_pattern2}). Skipping.")
            continue # Skip this chromosome if TPED is missing

        if len(tped_files_found) > 1:
             logging.warning(f"Found multiple potential input TPED files for chromosome {chrom_str}: {tped_files_found}. Using the first one: {tped_files_found[0]}.")
        tped_file = tped_files_found[0]
        base_tped = os.path.basename(tped_file)

        logging.info(f"Processing Chr {chrom_str} from TPED file: {base_tped}")
        start_time_chr = time.time()

        # Construct path to the corresponding ancestral allele file
        aa_ref_filename = aa_filename_func(chrom_str)
        aa_ref_filepath = os.path.join(aa_ref_dir, aa_ref_filename)

        # Load ancestral alleles for this chromosome
        ancestral_allele_map = load_func(aa_ref_filepath)
        if ancestral_allele_map is None:
            logging.error(f"  Skipping chromosome {chrom_str} due to issues loading ancestral data from {aa_ref_filename}.")
            overall_success = False # Critical failure if AA ref is missing for a requested chromosome
            continue

        # Define output AA-TPED file path
        aa_tped_filename = f"{pop}.{chrom_str}.aa.tped" # Consistent output naming
        aa_tped_filepath = os.path.join(output_dir, aa_tped_filename)

        variants_in_file = 0
        variants_written_file = 0
        ancestral_missing_file = 0
        file_had_errors = False

        # Process single TPED file
        try:
            # Use UTF-8 encoding for reading and writing
            with open(tped_file, 'r', encoding='utf-8') as infile, \
                 open(aa_tped_filepath, 'w', encoding='utf-8') as outfile:
                processed_files_count += 1 # Count successfully opened files
                # Loop through variants in the TPED file
                for line_num, line in enumerate(infile, 1):
                    variants_in_file += 1
                    try:
                        stripped_line = line.strip()
                        if not stripped_line: continue
                        fields = stripped_line.split()
                        if len(fields) < 5:
                            logging.warning(f"  Skipping malformed TPED line {line_num} in {base_tped} (fields < 5): {stripped_line}")
                            continue

                        tped_chrom_raw = fields[0]
                        rsid = fields[1]
                        cm = fields[2]
                        pos_str = fields[3]
                        genotypes = fields[4:]

                        # Consistency check (normalize TPED chrom from file line)
                        tped_chrom_norm = tped_chrom_raw.upper().replace("CHR", "")
                        if tped_chrom_norm != chrom_str:
                             logging.warning(f"  Chromosome mismatch line {line_num} in {base_tped}. Expected '{chrom_str}', TPED has '{tped_chrom_raw}'. Skipping line.")
                             continue

                        pos = int(pos_str)

                        # Lookup Ancestral Allele
                        ancestral_al = ancestral_allele_map.get(pos)

                        if ancestral_al is None:
                            ancestral_missing_file += 1
                            continue

                        # Convert genotypes
                        aa_genotypes = []
                        valid_tped_alleles = True
                        # ancestral_al is already uppercase

                        for allele in genotypes:
                            allele_upper = allele.upper()
                            if allele_upper == TPED_MISSING_ALLELE:
                                aa_genotypes.append('0')
                            elif allele_upper == ancestral_al:
                                aa_genotypes.append('0')
                            elif allele_upper in 'ACGT':
                                aa_genotypes.append('1')
                            else:
                                logging.warning(f"  Unexpected TPED allele '{allele}' at {chrom_str}:{pos} line {line_num} in {base_tped}. Skipping variant.")
                                valid_tped_alleles = False
                                break

                        if not valid_tped_alleles:
                            continue

                        # Write the new AA-TPED line using original TPED chr format, TAB separated
                        outfile.write(f"{tped_chrom_raw}\t{rsid}\t{cm}\t{pos_str}\t{' '.join(aa_genotypes)}\n")
                        variants_written_file += 1

                    except ValueError as ve:
                        logging.warning(f"  Skipping invalid data on TPED line {line_num} in {base_tped}: {ve}. Line: {line.strip()}")
                        continue
                    except Exception as e_inner:
                        logging.error(f"  ERROR processing TPED line {line_num} in {base_tped}: {e_inner}. Line: {line.strip()}", exc_info=True)
                        file_had_errors = True
                        continue

            # Finished processing lines in one TPED file
            end_time_chr = time.time()
            total_variants_processed_pop += variants_in_file
            total_variants_written_pop += variants_written_file
            total_ancestral_missing_pop += ancestral_missing_file

            logging.info(f"  Finished Chr {chrom_str} ({base_tped}) ({end_time_chr - start_time_chr:.2f}s):")
            logging.info(f"    Variants processed              : {variants_in_file:,}")
            logging.info(f"    Variants skipped (pos not found): {ancestral_missing_file:,}")
            logging.info(f"    Variants written to AA-TPED     : {variants_written_file:,}")
            if variants_in_file > 0 and variants_written_file == 0 and ancestral_missing_file == variants_in_file:
                 logging.info(f"    Note: No variants written; all positions might have been missing from the '{aaref}' reference for Chr {chrom_str}.")
            if file_had_errors:
                 logging.warning(f"    Errors occurred while processing lines in {base_tped}.")
                 overall_success = False

        except Exception as e_outer:
            logging.error(f"Error processing file {base_tped}: {e_outer}", exc_info=True)
            overall_success = False
            # Clean up potentially incomplete output file
            if os.path.exists(aa_tped_filepath):
                try:
                    os.remove(aa_tped_filepath)
                    logging.info(f"  Removed potentially incomplete output file: {aa_tped_filepath}")
                except OSError as e_rem:
                    logging.warning(f"  Could not remove incomplete output file {aa_tped_filepath}: {e_rem}")

    # End of loop through specified chromosomes
    logging.info("-" * 20)
    logging.info(f"Population {pop} Summary (Step 1 - {len(chromosomes_to_process)} chr requested, {processed_files_count} processed):")
    logging.info(f"  Total variants processed         : {total_variants_processed_pop:,}")
    logging.info(f"  Total skipped (AA missing)       : {total_ancestral_missing_pop:,}")
    logging.info(f"  Total written to AA-TPED         : {total_variants_written_pop:,}")

    if processed_files_count == 0 and chromosomes_to_process:
        logging.error(f"Step 1 failed: No input TPED files were found or processed for the requested chromosome(s).")
        raise FileNotFoundError("No input TPED files processed.")
    if not overall_success:
        raise RuntimeError(f"Step 1 failed for population {pop} with aaref {aaref}. Check logs.")
    elif processed_files_count < len(chromosomes_to_process):
        logging.warning(f"Step 1 completed for population {pop}, but some requested chromosomes may have been skipped due to missing input or errors.")
    else:
        logging.info(f"Step 1 completed successfully for population {pop}.")
