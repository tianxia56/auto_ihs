# run_ihs_calc.py
import os
import sys
import subprocess
import logging
import time
import shutil # To check for executables

# --- Constants ---
HAPBIN_EXECUTABLE = "/home/tx56/hapbin/build/ihsbin"
# Define SELSCAN executable if not in PATH, otherwise assume 'selscan' works
SELSCAN_EXECUTABLE = "selscan" # Or provide full path: "/path/to/selscan"

# --- Function to get available CPUs (from SLURM or default) ---
def get_available_cpus(default=4):
    """Gets the number of CPUs allocated by SLURM, or returns a default."""
    try:
        cpus = int(os.environ.get('SLURM_CPUS_PER_TASK', default))
        logging.debug(f"Detected {cpus} CPUs available (SLURM_CPUS_PER_TASK or default).")
        return max(1, cpus) # Ensure at least 1 CPU
    except (ValueError, TypeError):
        logging.warning(f"Could not parse SLURM_CPUS_PER_TASK. Using default: {default}")
        return max(1, default)

# --- Main iHS Calculation Function ---
def run_ihs_calculation(pop, aaref, rmap, software, alt_na, maf, chromosomes_to_process, base_maphap_dir, output_dir):
    """
    Executes iHS calculations directly using subprocess (Step 3).
    Processes only the chromosomes specified in chromosomes_to_process.
    """
    logging.info(f"Running Step 3: Execute iHS Calculation for Pop={pop}, AARef={aaref}, RMap={rmap}, Software={software}, Alt/NA={alt_na}, MAF=0.{maf}, Chromosomes={chromosomes_to_process}")
    os.makedirs(output_dir, exist_ok=True)

    # Check if executables exist
    if software == 'selscan' and not shutil.which(SELSCAN_EXECUTABLE):
         logging.error(f"{software.capitalize()} executable '{SELSCAN_EXECUTABLE}' not found or not executable. Please check installation and PATH.")
         raise FileNotFoundError(f"{SELSCAN_EXECUTABLE} not found")
    elif software == 'hapbin' and not shutil.which(HAPBIN_EXECUTABLE):
         logging.error(f"{software.capitalize()} executable '{HAPBIN_EXECUTABLE}' not found or not executable. Please check path.")
         raise FileNotFoundError(f"{HAPBIN_EXECUTABLE} not found")

    overall_success = True
    calculations_run = 0
    available_threads = get_available_cpus() # Get CPU count for selscan

    # Base name components for output files
    base_output_prefix = f"{pop}.{software}.{alt_na}.{aaref}.{rmap}.maf{maf}"

    # --- Loop through specified chromosomes ---
    for chrom_num in chromosomes_to_process:
        chrom = str(chrom_num)
        logging.info(f"  Executing iHS calculation for chromosome {chrom}...")
        start_time_chr = time.time()

        # Define input file basenames (output from Step 2)
        input_base = f"{pop}.{chrom}.{aaref}.{rmap}"
        tped_filepath = os.path.join(base_maphap_dir, f"{input_base}.tped")
        map_filepath = os.path.join(base_maphap_dir, f"{input_base}.map")
        hap_filepath = os.path.join(base_maphap_dir, f"{input_base}.hap")

        # Construct the command list
        cmd = []
        output_filepath_raw = None # Track the expected output file

        if software == 'selscan':
            # Check input file
            if not os.path.exists(tped_filepath):
                 logging.error(f"  Input TPED file not found for Chr {chrom}: {tped_filepath}. Skipping calculation.")
                 overall_success = False
                 continue
            # Define output base for selscan
            selscan_out_base = os.path.join(output_dir, f"{base_output_prefix}.chr{chrom}")
            output_filepath_raw = f"{selscan_out_base}.ihs.out" # Expected output
            cmd = [
                SELSCAN_EXECUTABLE, "--ihs", "--tped", tped_filepath,
                "--out", selscan_out_base,
                "--threads", str(available_threads), # Use detected threads
                "--maf", f"0.{maf}"
            ]
            if alt_na == 'alt':
                cmd.append("--alt")

        elif software == 'hapbin':
            # Check input files
            if not os.path.exists(map_filepath):
                 logging.error(f"  Input MAP file not found for Chr {chrom}: {map_filepath}. Skipping calculation.")
                 overall_success = False
                 continue
            if not os.path.exists(hap_filepath):
                 logging.error(f"  Input HAP file not found for Chr {chrom}: {hap_filepath}. Skipping calculation.")
                 overall_success = False
                 continue
            # Define output path for hapbin (matches norm expectation)
            maf_suffix_norm = f"maf{int(maf):d}"
            output_filepath_raw = os.path.join(output_dir, f"{pop}.{chrom}.{maf_suffix_norm}.ihs")
            cmd = [
                 HAPBIN_EXECUTABLE,
                 "--hap", hap_filepath,
                 "--map", map_filepath,
                 "--out", output_filepath_raw,
                 "--minmaf", f"0.{maf}"
            ]
        else:
            # This case should not be reached due to argument parsing, but good practice
            logging.error(f"  Unsupported software '{software}' specified.")
            overall_success = False
            continue

        # --- Execute the command directly ---
        logging.info(f"    Executing command: {' '.join(cmd)}")
        try:
            # Run the command, wait for completion, capture output
            process = subprocess.run(
                cmd,
                capture_output=True,
                text=True,            # Decode stdout/stderr as text
                check=True,           # Raise CalledProcessError on non-zero exit
                encoding='utf-8'      # Specify encoding
            )
            # Log stdout/stderr even on success for debugging (can be verbose)
            if process.stdout:
                 logging.info(f"      {software} stdout:\n{process.stdout.strip()}")
            if process.stderr:
                 logging.info(f"      {software} stderr:\n{process.stderr.strip()}") # Use info level for stderr as tools might print info there

            calculations_run += 1
            logging.info(f"    Finished Chr {chrom} calculation in {time.time() - start_time_chr:.2f} seconds.")

            # Optional: Check if the expected output file was actually created
            if output_filepath_raw and not os.path.exists(output_filepath_raw):
                 logging.warning(f"    {software} completed for Chr {chrom}, but expected output file '{output_filepath_raw}' was not found!")
                 overall_success = False # Treat missing output as failure

        except FileNotFoundError as e:
             # This catches if the executable (selscan/hapbin) itself isn't found
             logging.error(f"  Executable not found for command: {' '.join(cmd)}")
             logging.error(f"  Error details: {e}")
             overall_success = False
             break # Stop if executable is missing
        except subprocess.CalledProcessError as e:
            # Log detailed error if the command failed (non-zero exit status)
            logging.error(f"  {software} command failed for Chr {chrom} with exit code {e.returncode}.")
            logging.error(f"  Command: {' '.join(e.cmd)}")
            logging.error(f"  Stdout:\n{e.stdout.strip()}")
            logging.error(f"  Stderr:\n{e.stderr.strip()}")
            overall_success = False
            # Decide whether to continue with other chromosomes or stop
            # continue # Let's try to continue for now
        except Exception as e:
            # Catch any other unexpected errors during execution
            logging.error(f"  An unexpected error occurred during {software} execution for Chr {chrom}: {e}", exc_info=True)
            overall_success = False
            # continue # Let's try to continue

    # End of chromosome loop
    logging.info("-" * 20)
    logging.info(f"Step 3 execution summary:")
    logging.info(f"  Requested chromosomes: {len(chromosomes_to_process)}")
    logging.info(f"  Calculations executed: {calculations_run}")

    if calculations_run == 0 and chromosomes_to_process:
         logging.error("Step 3 failed: No iHS calculations were successfully executed (check for input file/executable errors).")
         raise RuntimeError("No iHS calculations executed in Step 3.")
    if not overall_success:
        raise RuntimeError(f"Step 3 encountered errors during iHS calculation. Check logs in {output_dir} and the main pipeline log.")
    elif calculations_run < len(chromosomes_to_process):
         logging.warning(f"Step 3 completed, but only {calculations_run}/{len(chromosomes_to_process)} calculations were executed successfully. Check logs.")
    else:
        logging.info("Step 3 iHS calculations completed successfully.")
        logging.info(f"  Raw iHS results can be found in: {output_dir}")
