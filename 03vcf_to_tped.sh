#!/bin/bash

# --- Configuration ---
VCF_INPUT_DIR="/home/tx56/palmer_scratch/100kga/pjl_ihs/qc_vcf"  # Directory with PJL.*.qc.vcf files
TPED_OUTPUT_DIR="/home/tx56/palmer_scratch/100kga/pjl_ihs/tped" # Output directory for .tped files
LOG_DIR="/home/tx56/palmer_scratch/100kga/pjl_ihs/logs" # Optional: Directory for logs
# --- End Configuration ---

# Create output and log directories if they don't exist
mkdir -p "$TPED_OUTPUT_DIR"
mkdir -p "$LOG_DIR"

# Get a timestamp for the log file (optional)
TIMESTAMP=$(date +%Y%m%d_%H%M%S)
LOG_FILE="${LOG_DIR}/bash_vcf_to_tped_${TIMESTAMP}.log"

echo "Starting VCF to TPED conversion using Bash/AWK..."
echo "Input VCF Directory: $VCF_INPUT_DIR"
echo "Output TPED Directory: $TPED_OUTPUT_DIR"
# echo "Detailed log (includes warnings): $LOG_FILE" # Uncomment if using log file redirection
echo "--------------------------------------"

# Check if input directory exists
if [ ! -d "$VCF_INPUT_DIR" ]; then
    echo "Error: Input directory not found: $VCF_INPUT_DIR" >&2
    exit 1
fi

# Find VCF files
shopt -s nullglob # Prevent loop from running if no files match
vcf_files=("$VCF_INPUT_DIR"/PJL.*.qc.vcf)
shopt -u nullglob # Turn off nullglob

if [ ${#vcf_files[@]} -eq 0 ]; then
    echo "Error: No VCF files found matching pattern '$VCF_INPUT_DIR/PJL.*.qc.vcf'" >&2
    exit 1
fi

echo "Found ${#vcf_files[@]} VCF files to convert."

conversion_errors=0

# Loop through each VCF file
for vcf_file in "${vcf_files[@]}"; do
    base_vcf=$(basename "$vcf_file")

    # Extract chromosome number using parameter expansion or sed
    # Using sed for robustness in case of complex filenames later
    chrom=$(echo "$base_vcf" | sed -E 's/^PJL\.([0-9]+)\.qc\.vcf$/\1/')

    if [ -z "$chrom" ]; then
        echo "Warning: Could not extract chromosome number from '$base_vcf'. Skipping." >&2
        ((conversion_errors++))
        continue
    fi

    # Define output TPED file path
    tped_file="${TPED_OUTPUT_DIR}/PJL.${chrom}.tped"

    echo "Processing: $base_vcf  ==>  $(basename "$tped_file")"

    # Use AWK to perform the conversion
    # Explanation of awk script:
    # BEGIN { FS="\t"; OFS=" "; TPED_MISSING="0" } # Set separators, define TPED missing allele
    # /^#/ { next } # Skip header lines
    # {
    #   chrom = $1; pos = $2; id = $3; ref = $4; alt_full = $5; format_field = $9
    #   split(alt_full, alt_arr, ","); alt = alt_arr[1] # Handle multi-allelic, take first ALT
    #
    #   # Find GT index in FORMAT field
    #   gt_idx = -1; split(format_field, fmt_arr, ":");
    #   for(i=1; i<=length(fmt_arr); i++) { if(fmt_arr[i]=="GT") {gt_idx=i; break} }
    #   if(gt_idx == -1) { print "Warning: No GT field in line: " $0 > "/dev/stderr"; next }
    #
    #   tped_line = chrom OFS id OFS "0" OFS pos # Start TPED line
    #
    #   # Process genotypes for each sample (fields 10 to NF)
    #   for(sample_idx=10; sample_idx<=NF; sample_idx++) {
    #     split($sample_idx, sample_data, ":");
    #     gt_call = (gt_idx <= length(sample_data)) ? sample_data[gt_idx] : "./." # Get GT, default to missing if error
    #
    #     n_alleles = split(gt_call, vcf_codes, /[/|]/); # Split GT code (e.g., 0/1 -> vcf_codes[1]=0, vcf_codes[2]=1)
    #
    #     # Output two alleles for TPED format
    #     for (allele_pos=1; allele_pos<=2; allele_pos++) {
    #         allele_out = TPED_MISSING # Default to missing
    #         if (allele_pos <= n_alleles) { # If VCF provided info for this allele position
    #             code = vcf_codes[allele_pos]
    #             if (code == "0") { allele_out = (ref==".") ? TPED_MISSING : ref } # Map 0 to REF (handle missing REF)
    #             else if (code == "1") { allele_out = (alt==".") ? TPED_MISSING : alt } # Map 1 to ALT (handle missing ALT)
    #             else if (code != ".") { print "Warning: Unexpected VCF code '" code "' in GT '" gt_call "' line: " $0 > "/dev/stderr" }
    #             # if code == ".", allele_out remains TPED_MISSING
    #         }
    #         tped_line = tped_line OFS allele_out
    #     }
    #   }
    #   print tped_line # Print the completed TPED line
    # }

    awk '
    BEGIN { FS="\t"; OFS=" "; TPED_MISSING="0" }
    /^#/ { next } # Skip header lines that start with #
    {
      # Basic VCF field extraction
      chrom = $1; pos = $2; id = $3; ref = $4; alt_full = $5; format_field = $9

      # Handle multi-allelic sites: take the first ALT allele
      split(alt_full, alt_arr, ",");
      alt = alt_arr[1]

      # Find the index of the GT field within the FORMAT string
      gt_idx = -1;
      split(format_field, fmt_arr, ":");
      for(i=1; i<=length(fmt_arr); i++) {
        if(fmt_arr[i]=="GT") {
          gt_idx=i;
          break
        }
      }

      # Skip line if GT field is not found in FORMAT
      if(gt_idx == -1) {
        # print "Warning: No GT field found in FORMAT for line: " $0 > "/dev/stderr" # Optional: Log warning
        next
      }

      # Start building the TPED output line: Chr, SNP_ID, Gen_Dist(0), Position
      tped_line = chrom OFS id OFS "0" OFS pos

      # Iterate through sample genotype fields (from field 10 to the end)
      for(sample_idx=10; sample_idx<=NF; sample_idx++) {
        # Split the sample field (e.g., "0/1:30:...") by ":"
        split($sample_idx, sample_data, ":");

        # Extract the GT call using the index found earlier; default to missing if not present
        gt_call = (gt_idx <= length(sample_data) && sample_data[gt_idx] != "") ? sample_data[gt_idx] : "./."

        # Split the GT call (e.g., "0/1" or "0|1") into individual allele codes
        n_alleles = split(gt_call, vcf_codes, /[/|]/);

        # TPED requires diploid representation (two alleles per sample)
        for (allele_pos=1; allele_pos<=2; allele_pos++) {
            allele_out = TPED_MISSING # Default allele is the TPED missing code '0'

            # Check if the VCF genotype call provides information for this allele position
            if (allele_pos <= n_alleles) {
                code = vcf_codes[allele_pos]
                current_ref = (ref==".") ? TPED_MISSING : ref # Use TPED missing if VCF REF is '.'
                current_alt = (alt==".") ? TPED_MISSING : alt # Use TPED missing if VCF ALT is '.'

                if (code == "0") { allele_out = current_ref }         # VCF '0' maps to REF allele
                else if (code == "1") { allele_out = current_alt }    # VCF '1' maps to first ALT allele
                else if (code != ".") {
                    # Optionally warn about unexpected codes (e.g., '2' for second ALT)
                    # print "Warning: Unexpected VCF code '" code "' in GT '" gt_call "' line: " $0 > "/dev/stderr"
                    allele_out = TPED_MISSING # Treat unexpected codes as missing
                }
                # If code is '.', allele_out remains TPED_MISSING
            }
            # Append the determined allele (REF, ALT, or TPED_MISSING) to the line
            tped_line = tped_line OFS allele_out
        }
      }
      # Print the complete TPED line for this variant
      print tped_line

    }' "$vcf_file" > "$tped_file"

    # Check awk exit status (optional but good practice)
    if [ $? -ne 0 ]; then
        echo "Error: AWK processing failed for $base_vcf" >&2
        ((conversion_errors++))
        # Optional: Remove potentially incomplete output file
        # rm -f "$tped_file"
    fi

done

# --- Final Status Report ---
echo "--------------------------------------"
if [ $conversion_errors -eq 0 ]; then
    echo "VCF to TPED conversion completed successfully."
    echo "Output files are in: $TPED_OUTPUT_DIR"
    exit 0
else
    echo "VCF to TPED conversion finished with $conversion_errors errors/warnings." >&2
    echo "Please check the output directory and potential warnings above." >&2
    exit 1
fi
