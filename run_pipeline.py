# run_pipeline.py
import os
import sys
import argparse
import logging
import time
import re

# Import functions from other modules
import match_aa
import add_rmap
import run_ihs_calc
import normalize_ihs

# --- Configuration ---
# Base directories
BASE_TPED_DIR = "/home/tx56/palmer_scratch/100kga/pjl_ihs/tped"
BASE_AA_REF_DIR_ARG = "/vast/palmer/pi/reilly/jfa38/datasets_for_annotation/ErinG_ARG_hg19_AncestralAlleles/tidied_AA_tsvs"
BASE_AA_REF_DIR_JV = "/home/tx56/palmer_scratch/100kga/jv_inputs/aa_ref"
BASE_AA_REF_DIR_PHV = "/home/tx56/palmer_scratch/100kga/phv_inputs/aa_ref"

# Renamed BASE_RMAP_DIR_ARG to BASE_RMAP_DIR_PYRHO
BASE_RMAP_DIR_PYRHO = "/home/tx56/selscan_maps/PJL"
BASE_RMAP_DIR_JV = "/home/tx56/palmer_scratch/100kga/jv_inputs/rmap"
BASE_RMAP_DIR_PHV = "/home/tx56/Downloads/pophum_genetic_map_b37"
BASE_RMAP_DIR_TIAN = "/home/tx56/palmer_scratch/100kga/ihs_auto/tian_rmap"

# Intermediate/Output Base Directories
STEP1_OUTPUT_BASE = "/home/tx56/palmer_scratch/100kga/ihs_auto/aa_tped"
STEP2_OUTPUT_BASE = "/home/tx56/palmer_scratch/100kga/ihs_auto/qced_input_maphap"
STEP3_OUTPUT_BASE = "/home/tx56/palmer_scratch/100kga/ihs_auto/ihs"
STEP4_OUTPUT_BASE = "/home/tx56/palmer_scratch/100kga/ihs_auto"

LOG_DIR = "/home/tx56/palmer_scratch/100kga/ihs_auto/logs"

# --- Logging Setup ---
os.makedirs(LOG_DIR, exist_ok=True)
log_filename = os.path.join(LOG_DIR, f"pipeline_run_{time.strftime('%Y%m%d_%H%M%S')}.log")
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler(log_filename),
        logging.StreamHandler(sys.stdout)
    ]
)

# --- Argument Parsing ---
def parse_args():
    parser = argparse.ArgumentParser(description="Run the full iHS analysis pipeline.")
    parser.add_argument('--pop', required=True, help="Population identifier (e.g., PJL)")
    parser.add_argument('--aaref', required=True, choices=['arg', 'jv', 'phv'], help="Ancestral allele reference set")
    # Updated rmap choices
    parser.add_argument('--rmap', required=True, choices=['pyrho', 'jv', 'phv', 'tian'], help="Recombination map source ('pyrho' uses the ARG map)")
    parser.add_argument('--software', required=True, choices=['selscan', 'hapbin'], help="iHS calculation software")
    parser.add_argument('--alt_na', required=True, choices=['alt', 'na'], help="'alt' for selscan --alt flag, 'na' otherwise")
    parser.add_argument('--maf', required=True, choices=['01', '05'], help="Minor allele frequency cutoff (e.g., 05 for 0.05)")
    # Added optional --chr argument
    parser.add_argument('--chr', type=int, choices=range(1, 23), default=None, help="Optional: Specify a single chromosome (1-22) to process.")

    args = parser.parse_args()

    # --- Parameter Validation ---
    if args.alt_na == 'alt' and args.software != 'selscan':
        parser.error("--alt_na 'alt' is only compatible with --software 'selscan'")
    if args.software == 'hapbin' and args.alt_na != 'na':
         logging.warning("--alt_na is ignored when --software is 'hapbin'. Setting to 'na'.")
         args.alt_na = 'na'

    return args

# --- Main Pipeline Function ---
def main():
    args = parse_args()
    pipeline_start_time = time.time()

    logging.info("=" * 50)
    logging.info(f"Starting iHS Pipeline Run")
    logging.info(f"Parameters:")
    logging.info(f"  Population: {args.pop}")
    logging.info(f"  AA Reference: {args.aaref}")
    logging.info(f"  Rec. Map: {args.rmap}")
    logging.info(f"  Software: {args.software}")
    logging.info(f"  Alt/NA: {args.alt_na}")
    logging.info(f"  MAF: 0.{args.maf}")
    if args.chr:
        logging.info(f"  Chromosome: {args.chr}")
        chromosomes_to_process = [args.chr]
    else:
        logging.info(f"  Chromosome: ALL (1-22)")
        chromosomes_to_process = list(range(1, 23))
    logging.info("=" * 50)

    # --- Define Dynamic Paths ---
    # These directory paths remain the same, filenames inside will depend on chromosome list
    step1_output_dir = os.path.join(STEP1_OUTPUT_BASE, f"{args.pop}_{args.aaref}")
    step2_output_dir = os.path.join(STEP2_OUTPUT_BASE, f"{args.pop}_{args.aaref}_{args.rmap}")
    norm_map_dir = step2_output_dir # Directory where map files reside for norm step
    step3_output_dir = os.path.join(STEP3_OUTPUT_BASE, f"{args.pop}_{args.software}_{args.alt_na}_{args.aaref}_{args.rmap}_maf{args.maf}")
    final_output_prefix = f"{args.pop}.{args.software}.{args.alt_na}.{args.aaref}.{args.rmap}.{args.maf}"
    # Final Output Zip file path - adjust if needed based on actual norm output
    # If processing single chr, final output name might need adjustment, but let's keep it genome-wide for now.
    # Normalization combines results *before* saving the final file.
    final_gz_path = os.path.join(STEP4_OUTPUT_BASE, f"{final_output_prefix}.norm.ihs.gz")

    # --- Execute Pipeline Steps ---
    try:
        # Step 1: Match AA Ref
        logging.info("-" * 30)
        logging.info(f"Step 1: Matching Ancestral Alleles (aaref={args.aaref})")
        step1_start = time.time()
        match_aa.run_match_aa(
            pop=args.pop,
            aaref=args.aaref,
            chromosomes_to_process=chromosomes_to_process, # Pass the list
            base_tped_dir=BASE_TPED_DIR,
            base_aa_ref_dir_arg=BASE_AA_REF_DIR_ARG,
            base_aa_ref_dir_jv=BASE_AA_REF_DIR_JV,
            base_aa_ref_dir_phv=BASE_AA_REF_DIR_PHV,
            output_dir=step1_output_dir
        )
        logging.info(f"Step 1 completed in {time.time() - step1_start:.2f} seconds.")

        # Step 2: Add Recombination Map
        logging.info("-" * 30)
        logging.info(f"Step 2: Adding Recombination Map (rmap={args.rmap})")
        step2_start = time.time()
        add_rmap.run_add_rmap(
            pop=args.pop,
            aaref=args.aaref,
            rmap=args.rmap,
            chromosomes_to_process=chromosomes_to_process, # Pass the list
            base_aa_tped_dir=step1_output_dir,
            # Pass renamed base map dir
            base_rmap_dir_pyrho=BASE_RMAP_DIR_PYRHO, # Renamed
            base_rmap_dir_jv=BASE_RMAP_DIR_JV,
            base_rmap_dir_phv=BASE_RMAP_DIR_PHV,
            base_rmap_dir_tian=BASE_RMAP_DIR_TIAN,
            output_dir=step2_output_dir
        )
        logging.info(f"Step 2 completed in {time.time() - step2_start:.2f} seconds.")

        # Step 3: Run iHS Calculation
        logging.info("-" * 30)
        logging.info(f"Step 3: Running iHS Calculation (software={args.software}, alt_na={args.alt_na}, maf=0.{args.maf})")
        step3_start = time.time()
        run_ihs_calc.run_ihs_calculation(
            pop=args.pop,
            aaref=args.aaref,
            rmap=args.rmap,
            software=args.software,
            alt_na=args.alt_na,
            maf=args.maf,
            chromosomes_to_process=chromosomes_to_process, # Pass the list
            base_maphap_dir=step2_output_dir,
            output_dir=step3_output_dir
        )
        logging.info(f"Step 3 completed in {time.time() - step3_start:.2f} seconds.")


        # Step 4: Normalize iHS Results
        logging.info("-" * 30)
        logging.info(f"Step 4: Normalizing iHS Results (software={args.software}, maf=0.{args.maf})")
        step4_start = time.time()
        normalize_ihs.run_normalization(
             pop=args.pop,
             aaref=args.aaref,
             rmap=args.rmap,
             software=args.software,
             alt_na = args.alt_na,
             maf=args.maf,
             chromosomes_to_process=chromosomes_to_process, # Pass the list
             base_ihs_dir=step3_output_dir,
             base_map_dir=norm_map_dir,
             output_dir=STEP4_OUTPUT_BASE,
             output_prefix=final_output_prefix
        )
        logging.info(f"Step 4 completed in {time.time() - step4_start:.2f} seconds.")


        logging.info("=" * 50)
        logging.info(f"Pipeline completed successfully in {time.time() - pipeline_start_time:.2f} seconds.")
        logging.info(f"Final normalized output expected at: {final_gz_path}")
        logging.info("=" * 50)

    except Exception as e:
        logging.error(f"Pipeline failed: {e}", exc_info=True)
        sys.exit(1)

if __name__ == "__main__":
    main()
