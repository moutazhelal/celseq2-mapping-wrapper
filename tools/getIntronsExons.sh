#!/bin/bash

# Derived from: Mouse gastruloids transcriptomics analysis
# Upstream: https://github.com/anna-alemany/mouseGastruloids_scRNAseq_tomoseq
# Original authors (c) upstream contributors; see upstream history
# Modifications (c) 2025 Moutaz Helal
# License: GPL-3.0
# SPDX-License-Identifier: GPL-3.0-only

if [ $# -ne 7 ]
then
    echo "Please, give 7 input parameters:"
    echo "1) sorted bam file (map with star)"
    echo "2) file with introns"
    echo "3) file with exons"
    echo "4) output root files"
    echo "5) path to bedtools"
    echo "6) path to samtools"
    echo "7) path to scripts directory (containing countExonsIntrons.py)"
    exit 1
fi

inbam=$1
intron=$2
exon=$3
out=$4
p2b_raw=$5
p2s_samtools=$6
p2s_scripts=$7

# ====================================================================
# ROBUST BEDTOOLS PATH HANDLING
# ====================================================================
# Handle different ways users might specify bedtools path:
# 1. Just "bedtools" (command name)
# 2. "/usr/bin/bedtools" (full path to executable)
# 3. "/usr/bin/" (directory path - incorrect but handle gracefully)
# 4. "/conda/envs/env/bin/bedtools" (conda environment)

# Function to determine the correct bedtools command
get_bedtools_cmd() {
    local bedtools_input="$1"

    # Remove trailing slash if present
    bedtools_input="${bedtools_input%/}"

    # Test different scenarios
    if command -v "$bedtools_input" >/dev/null 2>&1; then
        # Direct executable (bedtools or /path/to/bedtools)
        echo "$bedtools_input"
    elif [[ -x "$bedtools_input" ]]; then
        # Full path to executable file
        echo "$bedtools_input"
    elif [[ -d "$bedtools_input" ]] && [[ -x "$bedtools_input/bedtools" ]]; then
        # Directory path containing bedtools executable
        echo "$bedtools_input/bedtools"
    elif [[ "$bedtools_input" == *"/" ]]; then
        # Looks like a directory path, try appending bedtools
        local dir_path="${bedtools_input%/}"
        if [[ -x "$dir_path/bedtools" ]]; then
            echo "$dir_path/bedtools"
        else
            echo "ERROR: bedtools not found at $dir_path/bedtools" >&2
            return 1
        fi
    else
        # Try as-is and hope for the best
        echo "$bedtools_input"
    fi
}

# Get the proper bedtools command
p2b=$(get_bedtools_cmd "$p2b_raw")
if [[ $? -ne 0 ]]; then
    echo "FATAL: Cannot determine correct bedtools command from: $p2b_raw"
    exit 1
fi

# Verify bedtools works
if ! "$p2b" --version >/dev/null 2>&1; then
    echo "FATAL: bedtools command not working: $p2b"
    echo "Input was: $p2b_raw"
    exit 1
fi

# Similarly handle samtools (though this is less commonly problematic)
if ! "$p2s_samtools" --version >/dev/null 2>&1; then
    echo "FATAL: samtools command not working: $p2s_samtools"
    exit 1
fi

echo "Processing: $(basename $inbam)"
echo "  Introns BED: $intron"
echo "  Exons BED: $exon"
echo "  Output prefix: $out"
echo "  Bedtools: $p2b (from input: $p2b_raw)"
echo "  Samtools: $p2s_samtools"
echo "  Scripts dir: $p2s_scripts"

# Verify input files exist
[[ -f "$inbam" ]] || { echo "ERROR: BAM file not found: $inbam"; exit 1; }
[[ -f "$intron" ]] || { echo "ERROR: Introns BED not found: $intron"; exit 1; }
[[ -f "$exon" ]] || { echo "ERROR: Exons BED not found: $exon"; exit 1; }

# ====================================================================
# MAIN PROCESSING (using corrected bedtools subcommand syntax)
# ====================================================================

# Convert BAM to BED using bedtools bamtobed subcommand
echo "Converting BAM to BED..."
"$p2b" bamtobed -i "$inbam" -split > "${out}_bam2bed.bed"
if [[ $? -ne 0 ]] || [[ ! -s "${out}_bam2bed.bed" ]]; then
    echo "ERROR: Failed to convert BAM to BED"
    exit 1
fi

# Use bedtools intersect subcommand (not separate executable)
echo "Finding intron intersections..."
"$p2b" intersect -a "${out}_bam2bed.bed" -b "$intron" -wb | awk '{if (($6==$10) && ($5==255)) print $1"\t"$2"\t"$3"\t"$4"\t"$6"\t"$NF}' | uniq > "${out}_intron1.bed" &

echo "Finding exon intersections..."
"$p2b" intersect -a "${out}_bam2bed.bed" -b "$exon" -wb | awk '{if (($6==$10) && ($5==255)) print $1"\t"$2"\t"$3"\t"$4"\t"$6"\t"$NF}' | uniq > "${out}_exon1.bed"
wait

# Check that intersect operations produced output
if [[ ! -f "${out}_intron1.bed" ]] || [[ ! -f "${out}_exon1.bed" ]]; then
    echo "ERROR: Intersect operations failed to produce output files"
    exit 1
fi

echo "Getting chromosome names..."
chroms=$("$p2s_samtools" view -H "$inbam" | grep '^@SQ' | awk -F 'SN:|\t' '{print $3}')

if [[ -z "$chroms" ]]; then
    echo "ERROR: No chromosomes found in BAM header"
    exit 1
fi

# Clean up existing output files
[[ -f "${out}_intron.bed" ]] && rm "${out}_intron.bed"
[[ -f "${out}_exon.bed" ]] && rm "${out}_exon.bed"

echo "Processing chromosomes: $chroms"
for i in $(echo $chroms | awk '{for (i=1; i<=NF; i++) print $i}')
do
    echo "  Processing chromosome: $i"
    awk -v i=$i '{if ($1==i) print $0}' "${out}_intron1.bed" > "${out}_TMP_intron1_${i}.bed" &
    awk -v i=$i '{if ($1==i) print $0}' "${out}_exon1.bed" > "${out}_TMP_exon1_${i}.bed"
    wait

    "$p2b" intersect -a "${out}_TMP_intron1_${i}.bed" -b "${out}_TMP_exon1_${i}.bed" -v >> "${out}_intron.bed" &
    "$p2b" intersect -a "${out}_TMP_exon1_${i}.bed" -b "${out}_TMP_intron1_${i}.bed" -v >> "${out}_exon.bed"
    wait

    rm "${out}_TMP_exon1_${i}.bed" "${out}_TMP_intron1_${i}.bed"
done

echo "Cleaning up intermediate files..."
rm "${out}_bam2bed.bed" "${out}_intron1.bed" "${out}_exon1.bed"

echo "Running count script..."
if [[ -f "${p2s_scripts}/countExonsIntrons.py" ]]; then
    python3 "${p2s_scripts}/countExonsIntrons.py" "${out}_intron.bed" "${out}_exon.bed" "$out"
    if [[ $? -eq 0 ]]; then
        echo "Count script completed successfully"
    else
        echo "WARNING: Count script failed, creating basic summaries..."
    fi
else
    echo "WARNING: countExonsIntrons.py not found at ${p2s_scripts}/countExonsIntrons.py"
    echo "Creating basic count summaries instead..."
fi

# Create basic summaries regardless (as backup or primary method)
if [[ -f "${out}_intron.bed" ]] && [[ -f "${out}_exon.bed" ]]; then
    {
        echo "intron_reads: $(wc -l < "${out}_intron.bed")"
        echo "exon_reads: $(wc -l < "${out}_exon.bed")"
        echo "total_reads: $(($(wc -l < "${out}_intron.bed") + $(wc -l < "${out}_exon.bed")))"
    } > "${out}_counts_summary.txt"

    # Create gene-level counts if gene info is available (column 6)
    if awk 'NF>=6 {print $6; exit}' "${out}_intron.bed" | grep -q .; then
        awk '{gene[$6]++} END{for(g in gene) print g"\t"gene[g]}' "${out}_intron.bed" | sort -k2,2nr > "${out}_intron_gene_counts.txt"
        echo "Created intron gene counts: ${out}_intron_gene_counts.txt"
    fi

    if awk 'NF>=6 {print $6; exit}' "${out}_exon.bed" | grep -q .; then
        awk '{gene[$6]++} END{for(g in gene) print g"\t"gene[g]}' "${out}_exon.bed" | sort -k2,2nr > "${out}_exon_gene_counts.txt"
        echo "Created exon gene counts: ${out}_exon_gene_counts.txt"
    fi

    echo "Basic count summaries created: ${out}_counts_summary.txt"
fi

echo "getIntronsExons.sh completed successfully"
echo "Output files:"
[[ -f "${out}_intron.bed" ]] && echo "  - ${out}_intron.bed ($(wc -l < "${out}_intron.bed") reads)"
[[ -f "${out}_exon.bed" ]] && echo "  - ${out}_exon.bed ($(wc -l < "${out}_exon.bed") reads)"
[[ -f "${out}_counts_summary.txt" ]] && echo "  - ${out}_counts_summary.txt"
[[ -f "${out}_intron_gene_counts.txt" ]] && echo "  - ${out}_intron_gene_counts.txt"
[[ -f "${out}_exon_gene_counts.txt" ]] && echo "  - ${out}_exon_gene_counts.txt"
