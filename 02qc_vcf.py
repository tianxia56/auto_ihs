import os
import glob
import re
from collections import defaultdict
import sys # For sys.stderr

def load_excluded_rsids(filepath):
    """Loads a set of rsIDs to exclude from a file (one ID per line)."""
    excluded_ids = set()
    if not filepath or not os.path.exists(filepath):
        print(f"Warning: Exclude rsID file not provided or not found: {filepath}", file=sys.stderr)
        return excluded_ids # Return empty set

    print(f"Loading outdated rsIDs to exclude from: {filepath}")
    try:
        with open(filepath, 'r') as f:
            for line in f:
                rsid = line.strip()
                if rsid.startswith("rs"): # Ensure we are adding valid looking rsIDs
                    excluded_ids.add(rsid)
        print(f"  Loaded {len(excluded_ids):,} rsIDs to exclude.")
    except Exception as e:
        print(f"Error loading excluded rsID file {filepath}: {e}", file=sys.stderr)
        return set()
    return excluded_ids

def process_vcf_files(input_dir, output_dir, report_file, exclude_rsid_file=None):
    """
    Processes VCF files according to the specified rules:
    0. Remove rows with rsIDs listed in exclude_rsid_file.
    1. Standardize remaining rsIDs (rename if not starting with "rs").
    2. Remove rows with duplicate rsIDs (after potential rename).
    3. Remove rows with duplicate positions (CHROM + POS).
    4. Remove rows where any genotype allele code is not strictly 0 or 1 (removes '.' and others).
    5. Save cleaned files to the output directory.
    6. Generate a summary report including counts for all steps and in table format.
    """
    print(f"Input VCF directory: {input_dir}")
    print(f"Output QC VCF directory: {output_dir}")
    print(f"Summary report file: {report_file}")
    if exclude_rsid_file:
        print(f"Exclude rsID file: {exclude_rsid_file}")

    excluded_rsids = load_excluded_rsids(exclude_rsid_file)

    os.makedirs(output_dir, exist_ok=True)
    print(f"Ensured output directory exists: {output_dir}")

    vcf_pattern = os.path.join(input_dir, "PJL.chr*.recode.vcf")
    input_files = sorted(glob.glob(vcf_pattern))

    if not input_files:
        print(f"Error: No VCF files found matching pattern '{vcf_pattern}'")
        return

    print(f"Found {len(input_files)} VCF files to process.")

    overall_stats = defaultdict(int)
    chromosome_stats = {}
    # --- Step 4 Definition: Define the ONLY valid genotype codes ---
    # --- This set EXCLUDES '.' ---
    valid_genotype_codes = {'0', '1'}
    print(f"Strict genotype filtering enabled: Only codes {valid_genotype_codes} allowed (variants with '.' or other codes will be removed).")

    for vcf_file_path in input_files:
        filename = os.path.basename(vcf_file_path)
        print(f"\nProcessing file: {filename}...")

        match = re.search(r'chr(\d+|X|Y|M)', filename, re.IGNORECASE)
        if not match:
            print(f"  Warning: Could not extract chromosome number from filename '{filename}'. Skipping.")
            continue
        chrom_str = match.group(1)

        output_filename = f"PJL.{chrom_str}.qc.vcf"
        output_file_path = os.path.join(output_dir, output_filename)

        initial_variants = 0
        outdated_rsids_removed = 0
        ids_renamed = 0
        duplicate_ids_removed = 0
        duplicate_positions_removed = 0
        invalid_genotypes_removed = 0 # Counter now includes variants with '.'
        final_variants_written = 0
        seen_ids = set()
        seen_positions = set()

        try:
            with open(vcf_file_path, 'r') as infile, open(output_file_path, 'w') as outfile:
                header_written = False
                gt_index = -1

                for line_num, line in enumerate(infile, 1):
                    if line.startswith('#'):
                        outfile.write(line)
                        header_written = True
                        continue

                    if not header_written:
                         print(f"  Error: Data line encountered before VCF header in {filename}. Aborting file.", file=sys.stderr)
                         raise ValueError("VCF format error: Missing or misplaced header.")

                    initial_variants += 1
                    fields = line.strip().split('\t')

                    if len(fields) < 9:
                        print(f"  Warning: Skipping line {line_num} in {filename} - insufficient fields (<9): {line.strip()}", file=sys.stderr)
                        continue

                    chrom = fields[0]
                    pos = fields[1]
                    rsid_original = fields[2]
                    format_field = fields[8]
                    genotype_fields = fields[9:]

                    # --- 0. Check for outdated/excluded rsID ---
                    if rsid_original != '.' and rsid_original in excluded_rsids:
                        outdated_rsids_removed += 1
                        continue

                    # Determine GT index on first data line if needed
                    if gt_index == -1 and genotype_fields:
                        try:
                            format_parts = format_field.split(':')
                            gt_index = format_parts.index('GT')
                        except ValueError:
                            print(f"  Error: Cannot find 'GT' in FORMAT field ('{format_field}') on line {line_num}. Cannot process genotypes for this file. Skipping file.", file=sys.stderr)
                            raise ValueError("GT field missing in FORMAT")

                    # --- 1. Check and potentially rename rsID ---
                    if not rsid_original.startswith("rs") and rsid_original != '.':
                        current_id = f"rspos{pos}"
                        if current_id != rsid_original:
                             ids_renamed += 1
                    else:
                        current_id = rsid_original

                    # --- 2. Check for duplicate rsID ---
                    if current_id != '.' and current_id in seen_ids:
                        duplicate_ids_removed += 1
                        continue
                    if current_id != '.':
                        seen_ids.add(current_id)

                    # --- 3. Check for duplicate position ---
                    position_key = (chrom, pos)
                    if position_key in seen_positions:
                        duplicate_positions_removed += 1
                        if current_id != '.' and current_id in seen_ids:
                            try: seen_ids.remove(current_id)
                            except KeyError: pass
                        continue
                    seen_positions.add(position_key)

                    # --- 4. Check for invalid/missing genotype codes (strict: only 0 or 1 allowed) ---
                    # --- This section ensures '.' leads to removal ---
                    invalid_or_missing_genotype_found = False
                    if not genotype_fields:
                         invalid_or_missing_genotype_found = True # Skip if no samples
                    elif gt_index == -1:
                         invalid_or_missing_genotype_found = True # Skip if GT index unknown
                    else:
                        for sample_field in genotype_fields:
                            try:
                                sample_parts = sample_field.split(':')
                                if gt_index >= len(sample_parts):
                                    gt_call = "./." # Treat truncated field as having missing GT
                                else:
                                    gt_call = sample_parts[gt_index]

                                # --- Quick check: If GT call contains '.', it's invalid for this rule ---
                                if '.' in gt_call:
                                    invalid_or_missing_genotype_found = True
                                    break # Found missing data, no need to check further codes/samples

                                allele_codes = re.split(r'[/|]', gt_call)
                                for code in allele_codes:
                                    # --- The core check: is the code '0' or '1'? ---
                                    if code not in valid_genotype_codes: # valid_genotype_codes is {'0', '1'}
                                        invalid_or_missing_genotype_found = True
                                        break # Found invalid code (not '0' or '1'), stop checking this sample
                            except Exception as e_gt:
                                print(f"  Warning: Error parsing GT field '{sample_field}' at {chrom}:{pos} line {line_num}. Treating as invalid. Error: {e_gt}", file=sys.stderr)
                                invalid_or_missing_genotype_found = True
                                break # Stop checking this variant if parsing fails
                            if invalid_or_missing_genotype_found:
                                break # Stop checking other samples if invalid found in this one

                    # --- Step 4 Outcome: If any code wasn't '0' or '1' (including '.'), skip the row ---
                    if invalid_or_missing_genotype_found:
                        invalid_genotypes_removed += 1 # Increment the counter
                        # Clean up potentially added ID and Position
                        if current_id != '.' and current_id in seen_ids:
                            try: seen_ids.remove(current_id)
                            except KeyError: pass
                        if position_key in seen_positions:
                            try: seen_positions.remove(position_key)
                            except KeyError: pass
                        continue # Skip writing this row

                    # --- If passed ALL checks, write the modified line ---
                    fields[2] = current_id
                    outfile.write("\t".join(fields) + "\n")
                    final_variants_written += 1

        except FileNotFoundError:
             print(f"  Error: Input VCF file not found: {vcf_file_path}", file=sys.stderr)
             continue
        except ValueError as ve:
             print(f"  Skipping file {filename} due to critical error: {ve}", file=sys.stderr)
             if os.path.exists(output_file_path):
                 try: os.remove(output_file_path)
                 except OSError: pass
             continue
        except Exception as e:
            print(f"  Unexpected error processing file {filename}: {e}", file=sys.stderr)
            if os.path.exists(output_file_path):
                 try: os.remove(output_file_path)
                 except OSError: pass
            continue

        # --- Store stats ---
        stats = {
            'initial': initial_variants,
            'outdated': outdated_rsids_removed,
            'ids_renamed': ids_renamed,
            'dup_ids': duplicate_ids_removed,
            'dup_pos': duplicate_positions_removed,
            'inv_gt': invalid_genotypes_removed, # This count now includes variants with '.'
            'final': final_variants_written
        }
        chromosome_stats[chrom_str] = stats

        print(f"  Finished {filename}:")
        print(f"    Initial variants            : {initial_variants:,}")
        print(f"    Outdated rsIDs removed      : {outdated_rsids_removed:,}")
        print(f"    IDs renamed                 : {ids_renamed:,}")
        print(f"    Duplicate IDs removed       : {duplicate_ids_removed:,}")
        print(f"    Duplicate Positions removed : {duplicate_positions_removed:,}")
        print(f"    Invalid/Missing GT removed  : {invalid_genotypes_removed:,}") # Label reflects '.' removal
        print(f"    Final variants written      : {final_variants_written:,}")

        # --- Update overall stats ---
        overall_stats['initial'] += initial_variants
        overall_stats['outdated'] += outdated_rsids_removed
        overall_stats['ids_renamed'] += ids_renamed
        overall_stats['dup_ids'] += duplicate_ids_removed
        overall_stats['dup_pos'] += duplicate_positions_removed
        overall_stats['inv_gt'] += invalid_genotypes_removed
        overall_stats['final'] += final_variants_written

    # --- Write summary report ---
    print(f"\nWriting summary report to: {report_file}...")
    successfully_processed_count = len(chromosome_stats)
    total_files_found = len(input_files)
    try:
        with open(report_file, 'w') as report:
            report.write("VCF QC Summary Report\n")
            report.write("======================\n\n")
            report.write(f"Input Directory: {input_dir}\n")
            report.write(f"Output Directory: {output_dir}\n")
            if exclude_rsid_file:
                 report.write(f"Exclude rsID File: {exclude_rsid_file} ({len(excluded_rsids):,} IDs loaded)\n")
            report.write(f"Processed {successfully_processed_count} files successfully (out of {total_files_found} found).\n")
            # --- Explicitly state the filter rule ---
            report.write(f"Strict Genotype Filter: Kept only variants where ALL genotype codes were '0' or '1'.\n\n")


            report.write("--- QC Statistics Summary ---\n")

            # Adjusted column header slightly
            header = "{:<6} {:>15} {:>11} {:>12} {:>10} {:>10} {:>10} {:>15}".format(
                "chr", "initial_var", "OutdatedID", "ID_renamed", "Dup_ID", "Dup_pos", "Inv/MisGT", "final_var"
            )
            separator = "-" * len(header)
            row_fmt = "{:<6} {:>15,} {:>11,} {:>12,} {:>10,} {:>10,} {:>10,} {:>15,}"

            report.write(header + "\n")
            report.write(separator + "\n")

            def sort_key(chrom_key):
                if chrom_key.isdigit(): return (0, int(chrom_key))
                else: return (1, chrom_key)

            sorted_chroms = sorted(chromosome_stats.keys(), key=sort_key)

            for chrom in sorted_chroms:
                stats = chromosome_stats[chrom]
                report.write(row_fmt.format(
                    chrom,
                    stats['initial'],
                    stats['outdated'],
                    stats['ids_renamed'],
                    stats['dup_ids'],
                    stats['dup_pos'],
                    stats['inv_gt'], # This column reflects removal of non '0'/'1' codes (incl '.')
                    stats['final']
                ) + "\n")

            report.write(separator + "\n")
            report.write(row_fmt.format(
                "Sum",
                overall_stats['initial'],
                overall_stats['outdated'],
                overall_stats['ids_renamed'],
                overall_stats['dup_ids'],
                overall_stats['dup_pos'],
                overall_stats['inv_gt'], # This sum reflects removal of non '0'/'1' codes (incl '.')
                overall_stats['final']
            ) + "\n")

        print("Summary report written successfully.")

    except Exception as e:
        print(f"Error writing summary report {report_file}: {e}", file=sys.stderr)

    print("\nProcessing complete.")
    if successfully_processed_count < total_files_found:
         print(f"\nWarning: {total_files_found - successfully_processed_count} file(s) encountered errors and were skipped or incomplete. Check logs above.", file=sys.stderr)


# --- Configuration ---
INPUT_VCF_DIR = "/home/tx56/palmer_scratch/100kga/pjl_ihs/vcf"
OUTPUT_QC_DIR = "/home/tx56/palmer_scratch/100kga/pjl_ihs/qc_vcf"
REPORT_FILE = "/home/tx56/palmer_scratch/100kga/pjl_ihs/qc_vcf.txt"
EXCLUDE_RSID_FILE = "/home/tx56/1KGP/exclude_rsid.txt"

# --- Run the processing ---
if __name__ == "__main__":
    process_vcf_files(INPUT_VCF_DIR, OUTPUT_QC_DIR, REPORT_FILE, exclude_rsid_file=EXCLUDE_RSID_FILE)
