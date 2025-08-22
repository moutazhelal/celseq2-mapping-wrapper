#!/usr/bin/env bash
set -Eeuo pipefail

# =====================================================================
# CEL-seq2 mapping Pipeline (portable wrapper)
# ---------------------------------------------------------------------
# This script wraps your existing pipeline to make it portable for users.
# - Uses a manifest CSV with columns: sample,raw_fq_dir,cbc_file
# - Generates STAR genome index and intron/exon BEDs if missing
# - Locates bundled helper scripts relative to this file by default
# - Exposes clear CLI flags; prints helpful usage on --help
# =====================================================================

VERSION="1.0.0"

# -------- Defaults (can be overridden by CLI flags) --------
THREADS="${THREADS:-10}"
RESULTS_DIR="${RESULTS_DIR:-results}"
MANIFEST_CSV="${MANIFEST_CSV:-data/manifest.csv}"   # header: sample,raw_fq_dir,cbc_file

# Genome resources
GENOME_INDEX="${GENOME_INDEX:-}"                     # optional; will be created next to FASTA if not provided
GENOME_FASTA="${GENOME_FASTA:-}"
GTF="${GTF:-}"
SJDB_OVERHANG="${SJDB_OVERHANG:-99}"

# Intron/Exon resources
EXINT_DIR="${EXINT_DIR:-data/exint}"
EXINT_BASENAME="${EXINT_BASENAME:-EXINt}"   # will create {basename}_{introns,exons}.bed

# Tool locations (use PATH by default)
TRIM_GALORE="${TRIM_GALORE:-trim_galore}"
CUTADAPT_BIN="${CUTADAPT_BIN:-cutadapt}"
STAR_BIN="${STAR_BIN:-STAR}"
BEDTOOLS_BIN="${BEDTOOLS_BIN:-bedtools}"
SAMTOOLS_BIN="${SAMTOOLS_BIN:-samtools}"

# Helper scripts (bundled in repo's tools/ by default)
SCRIPT_DIR="$(cd -- "$(dirname "${BASH_SOURCE[0]}")" >/dev/null 2>&1 && pwd)"
SCRIPTS_DIR_DEFAULT="${SCRIPT_DIR}/tools"
SCRIPTS_DIR="${SCRIPTS_DIR:-${SCRIPTS_DIR_DEFAULT}}"
CONCATENATOR_PY="${CONCATENATOR_PY:-${SCRIPTS_DIR}/concatenator.py}"
EXINT_MAKER_PY="${EXINT_MAKER_PY:-${SCRIPTS_DIR}/get_IntronsExons_fromGTF_mod.py}"
GET_INTRONSEXONS="${GET_INTRONSEXONS:-${SCRIPTS_DIR}/getIntronsExons.sh}"

# concatenator.py parameters
BCR_READ="${BCR_READ:-R1}"
BIO_READ="${BIO_READ:-R2}"
LEN_CBC="${LEN_CBC:-6}"
LEN_UMI="${LEN_UMI:-6}"
UMI_FIRST="${UMI_FIRST:-true}"    # if "true", pass --umifirst
CBC_HD="${CBC_HD:-0}"

# -------- Pretty logging & error handling --------
banner () {
  echo "============================================================"
  echo "[$(date '+%F %T')] $*"
  echo "============================================================"
}
step () { echo "[$(date '+%F %T')]   → $*"; }
die () { echo "FATAL: $*" >&2; exit 1; }

trap 'ec=$?; echo "Error on line ${BASH_LINENO[0]} (cmd: ${BASH_COMMAND})" >&2; exit $ec' ERR

usage () {
  cat <<EOF
Mouse CEL-seq Pipeline v${VERSION}

USAGE:
  $(basename "$0") --manifest data/manifest.csv --genome-fasta path.fa --gtf path.gtf [options]

Required:
  --genome-fasta PATH         Genome FASTA (e.g., Mus_musculus.GRCm39.dna.primary_assembly.fa)
  --gtf PATH                  Annotation GTF
  --manifest PATH             CSV with columns: sample,raw_fq_dir,cbc_file

Optional:
  --results-dir DIR           Output directory (default: ${RESULTS_DIR})
  --genome-index DIR          Existing/target STAR index directory (default: created next to FASTA)
  --sjdb-overhang N           STAR sjdbOverhang (default: ${SJDB_OVERHANG})
  --threads N                 Threads for STAR and other tools (default: ${THREADS})

  # Helper resources
  --exint-dir DIR             Directory for EXINt BEDs (default: ${EXINT_DIR})
  --exint-basename NAME       Basename for EXINt BEDs (default: ${EXINT_BASENAME})

  # Tools (set if not on PATH)
  --trim-galore PATH          trim_galore executable (default: ${TRIM_GALORE})
  --cutadapt PATH             cutadapt executable (default: ${CUTADAPT_BIN})
  --star PATH                 STAR executable (default: ${STAR_BIN})
  --bedtools PATH             bedtools executable (default: ${BEDTOOLS_BIN})
  --samtools PATH             samtools executable (default: ${SAMTOOLS_BIN})
  --multiqc PATH             MultiQC executable (default: ${MULTIQC_BIN})
  --multiqc-config PATH      Optional MultiQC config YAML
  --no-multiqc               Disable MultiQC report generation

  # Scripts (override if stored elsewhere)
  --scripts-dir DIR           Folder containing helper scripts (default: ${SCRIPTS_DIR_DEFAULT})
  --concatenator PATH         concatenator.py (default: tools/concatenator.py)
  --exint-maker PATH          get_IntronsExons_fromGTF_mod.py (default: tools/get_IntronsExons_fromGTF_mod.py)
  --get-intronsexons PATH     getIntronsExons.sh (default: tools/getIntronsExons.sh)

  # concatenator.py parameters
  --bcread STR                Barcode read name (default: ${BCR_READ})
  --bioread STR               Bio read name (default: ${BIO_READ})
  --len-cbc N                 Barcode length (default: ${LEN_CBC})
  --len-umi N                 UMI length (default: ${LEN_UMI})
  --umi-first [true|false]    Whether UMI precedes CBC (default: ${UMI_FIRST})
  --cbc-hd N                  Hamming distance for CBC (default: ${CBC_HD})

  -h, --help                  Show this help and exit

Examples:
  $(basename "$0") --manifest data/manifest.csv \
    --genome-fasta /ref/Mus_musculus.GRCm39.dna.primary_assembly.fa \
    --gtf /ref/Mus_musculus.GRCm39.112.gtf \
    --threads 16

EOF
}

# -------- Arg parsing --------

# ---- Pre-parse tool override flags (do this before the main parser) ----
STAR_BIN="${STAR_BIN:-}"
BEDTOOLS_BIN="${BEDTOOLS_BIN:-}"
SAMTOOLS_BIN="${SAMTOOLS_BIN:-}"
TRIM_GALORE="${TRIM_GALORE:-}"
CUTADAPT_BIN="${CUTADAPT_BIN:-}"

if [[ $# -gt 0 ]]; then
  __newargs=()
  while [[ $# -gt 0 ]]; do
    case "$1" in
      --star) STAR_BIN="$2"; shift 2;;
      --bedtools) BEDTOOLS_BIN="$2"; shift 2;;
      --samtools) SAMTOOLS_BIN="$2"; shift 2;;
      --trim-galore) TRIM_GALORE="$2"; shift 2;;
      --cutadapt) CUTADAPT_BIN="$2"; shift 2;;
      --multiqc) MULTIQC_BIN="$2"; shift 2;;
      --multiqc-config) MULTIQC_CONFIG="$2"; shift 2;;
      --no-multiqc) RUN_MULTIQC="false"; shift 1;;
      *) __newargs+=("$1"); shift;;
    esac
  done

# ---- Final quality summary ----
run_multiqc
  set -- "${__newargs[@]}"
fi

# Defaults to environment (PATH) unless explicitly overridden above
: "${STAR_BIN:=STAR}"
: "${BEDTOOLS_BIN:=bedtools}"
: "${SAMTOOLS_BIN:=samtools}"
: "${TRIM_GALORE:=trim_galore}"
: "${CUTADAPT_BIN:=cutadapt}"

# ---- Tool resolution helpers ----
die () { echo "FATAL: $*" >&2; exit 1; }
__resolve_bin () {
  local var_name="$1" default_name="$2"
  local current; eval "current=\"\${${var_name}:-}\""
  local candidate="${current:-$default_name}"
  if ! command -v "$candidate" >/dev/null 2>&1; then
    die "Required tool '${var_name}' not found (tried: '${candidate}'). Install it in your environment or pass the CLI flag."
  fi
  local abs; abs="$(command -v "$candidate")"
  eval "${var_name}=\"${abs}\""
}
resolve_all_tools () {
  __resolve_bin STAR_BIN STAR
  __resolve_bin BEDTOOLS_BIN bedtools
  __resolve_bin SAMTOOLS_BIN samtools
  __resolve_bin TRIM_GALORE trim_galore
  __resolve_bin CUTADAPT_BIN cutadapt
  if [[ "${RUN_MULTIQC}" == "true" ]]; then
    __resolve_bin MULTIQC_BIN multiqc
  fi
}
print_tool_versions () {
  echo "Tools detected:"
  echo "  STAR:        ${STAR_BIN}  ($("$STAR_BIN" --version 2>/dev/null | head -n1 || true))"
  echo "  bedtools:    ${BEDTOOLS_BIN}  ($("$BEDTOOLS_BIN" --version 2>/dev/null | head -n1 || true))"
  echo "  samtools:    ${SAMTOOLS_BIN}  ($("$SAMTOOLS_BIN" --version 2>/dev/null | head -n1 || true))"
  echo "  trim_galore: ${TRIM_GALORE}  ($("$TRIM_GALORE" --version 2>/dev/null | head -n1 || true))"
  if ver="$("$CUTADAPT_BIN" --version 2>/dev/null)"; then
    echo "  cutadapt:    ${CUTADAPT_BIN}  (v${ver})"
  else
    echo "  cutadapt:    ${CUTADAPT_BIN}"
  fi
  if [[ "${RUN_MULTIQC}" == "true" ]]; then
    if ver="$("$MULTIQC_BIN" --version 2>/dev/null | head -n1)"; then
      echo "  multiqc:     ${MULTIQC_BIN}  (${ver})"
    else
      echo "  multiqc:     ${MULTIQC_BIN}"
    fi
  fi
}
# -----------------------------------------------------------------------
ARGS=()
while [[ $# -gt 0 ]]; do
  case "$1" in
    --manifest) MANIFEST_CSV="$2"; shift 2;;
    --results-dir) RESULTS_DIR="$2"; shift 2;;
    --genome-index) GENOME_INDEX="$2"; shift 2;;
    --genome-fasta) GENOME_FASTA="$2"; shift 2;;
    --gtf) GTF="$2"; shift 2;;
    --sjdb-overhang) SJDB_OVERHANG="$2"; shift 2;;
    --threads) THREADS="$2"; shift 2;;
    --exint-dir) EXINT_DIR="$2"; shift 2;;
    --exint-basename) EXINT_BASENAME="$2"; shift 2;;
    --trim-galore) TRIM_GALORE="$2"; shift 2;;
    --cutadapt) CUTADAPT_BIN="$2"; shift 2;;
    --star) STAR_BIN="$2"; shift 2;;
    --bedtools) BEDTOOLS_BIN="$2"; shift 2;;
    --samtools) SAMTOOLS_BIN="$2"; shift 2;;
    --scripts-dir) SCRIPTS_DIR="$2"; shift 2;;
    --concatenator) CONCATENATOR_PY="$2"; shift 2;;
    --exint-maker) EXINT_MAKER_PY="$2"; shift 2;;
    --get-intronsexons) GET_INTRONSEXONS="$2"; shift 2;;
    --bcread) BCR_READ="$2"; shift 2;;
    --bioread) BIO_READ="$2"; shift 2;;
    --len-cbc) LEN_CBC="$2"; shift 2;;
    --len-umi) LEN_UMI="$2"; shift 2;;
    --umi-first) UMI_FIRST="$2"; shift 2;;
    --cbc-hd) CBC_HD="$2"; shift 2;;
    -h|--help) usage; exit 0;;
    *) ARGS+=("$1"); shift;;
  esac
done
set -- "${ARGS[@]:-}"

# Required checks
[[ -n "${GENOME_FASTA}" ]] || { echo "ERROR: --genome-fasta is required"; usage; exit 2; }
[[ -n "${GTF}" ]] || { echo "ERROR: --gtf is required"; usage; exit 2; }
[[ -s "${MANIFEST_CSV}" ]] || { echo "ERROR: manifest not found: ${MANIFEST_CSV}"; usage; exit 2; }

# Create output dirs
mkdir -p "${RESULTS_DIR}" "${EXINT_DIR}"

# -------- Helpers --------
check_writable_dir () {
  local d="$1"
  mkdir -p "$d" || die "Cannot create directory: $d"
  [[ -w "$d" ]] || die "No write permission in: $d"
}
valid_star_index () {
  local idx="$1"
  [[ -d "$idx" && -s "$idx/Genome" && -s "$idx/SA" && -s "$idx/SAindex" ]]
}
ensure_star_index () {
  if [[ -n "${GENOME_INDEX}" ]] && valid_star_index "${GENOME_INDEX}"; then
    step "Found existing STAR index: ${GENOME_INDEX}"
    return 0
  fi
  if [[ -z "${GENOME_INDEX}" ]]; then
    GENOME_INDEX="$(dirname "$GENOME_FASTA")/STAR_Index"
  fi
  banner "STAR index not found → generating at: ${GENOME_INDEX}"
  check_writable_dir "${GENOME_INDEX}"
  [[ -s "${GENOME_FASTA}" ]] || die "GENOME_FASTA not found: ${GENOME_FASTA}"
  [[ -s "${GTF}" ]] || die "GTF not found: ${GTF}"
  "${STAR_BIN}" \
    --runThreadN "${THREADS}" \
    --runMode genomeGenerate \
    --genomeDir "${GENOME_INDEX}" \
    --genomeFastaFiles "${GENOME_FASTA}" \
    --sjdbGTFfile "${GTF}" \
    --sjdbOverhang "${SJDB_OVERHANG}"
  valid_star_index "${GENOME_INDEX}" || die "STAR index generation failed for ${GENOME_INDEX}"
  step "STAR index ready: ${GENOME_INDEX}"
}

ensure_exint () {
  local introns="${EXINT_DIR}/${EXINT_BASENAME}_introns.bed"
  local exons="${EXINT_DIR}/${EXINT_BASENAME}_exons.bed"
  if [[ -s "$introns" && -s "$exons" ]]; then
    step "Found EXINt BEDs: ${introns} , ${exons}"
    return 0
  fi
  banner "EXINt BEDs not found → generating in: ${EXINT_DIR}"
  check_writable_dir "${EXINT_DIR}"
  [[ -s "${GTF}" ]] || die "GTF not found: ${GTF}"
  [[ -s "${EXINT_MAKER_PY}" ]] || die "Generator script not found: ${EXINT_MAKER_PY}"
  ( cd "${EXINT_DIR}" && python3 "${EXINT_MAKER_PY}" "${GTF}" "${EXINT_BASENAME}" )
  [[ -s "$introns" && -s "$exons" ]] || die "Failed to generate EXINt BEDs in ${EXINT_DIR}"
  step "EXINt BEDs ready."
}

# Resolve tool paths now (prefers overrides; falls back to PATH)
resolve_all_tools
print_tool_versions

check_tools () {
  command -v python3 >/dev/null || die "python3 not found in PATH"
  command -v "${CUTADAPT_BIN}" >/dev/null || die "cutadapt not found: ${CUTADAPT_BIN}"
  command -v "${TRIM_GALORE}" >/dev/null || die "trim_galore not found: ${TRIM_GALORE}"
  command -v "${STAR_BIN}" >/dev/null || die "STAR not found: ${STAR_BIN}"
  command -v "${BEDTOOLS_BIN}" >/dev/null || die "bedtools not found: ${BEDTOOLS_BIN}"
  command -v "${SAMTOOLS_BIN}" >/dev/null || die "samtools not found: ${SAMTOOLS_BIN}"
  if [[ "${RUN_MULTIQC}" == "true" ]]; then command -v "${MULTIQC_BIN}" >/dev/null || die "multiqc not found: ${MULTIQC_BIN}"; fi
  [[ -s "${CONCATENATOR_PY}" ]] || die "concatenator.py not found: ${CONCATENATOR_PY}"
  [[ -s "${GET_INTRONSEXONS}" ]] || die "getIntronsExons.sh not found: ${GET_INTRONSEXONS}"
}

run_sample () {
  local SAMPLE="$1" RAW_FQ_DIR="$2" CBC_FILE="$3"
  local SAMPLE_DIR="${RESULTS_DIR}/${SAMPLE}"
  local LOG_DIR="${SAMPLE_DIR}/logs"

  banner "START sample: ${SAMPLE}"
  check_writable_dir "${SAMPLE_DIR}"
  check_writable_dir "${LOG_DIR}"

  local IN_PREFIX="${RAW_FQ_DIR}/${SAMPLE}"
  if [[ ! -f "${IN_PREFIX}_R1.fq.gz" || ! -f "${IN_PREFIX}_R2.fq.gz" ]]; then
    echo "ERROR: Missing ${IN_PREFIX}_R1.fq.gz or ${IN_PREFIX}_R2.fq.gz" | tee -a "${LOG_DIR}/00_inputs.log"
    return 1
  fi
  local OUT_PREFIX="${SAMPLE_DIR}/${SAMPLE}"

  # 0) Ensure prerequisites (index + EXINt)
  ensure_star_index
  ensure_exint

  # 1) Extract barcodes/UMIs → *_cbc.fastq
  step "(${SAMPLE}) Extracting cell barcodes and UMIs"
  python3 "${CONCATENATOR_PY}" \
    --fqf "${IN_PREFIX}" \
    --bcread "${BCR_READ}" \
    --bioread "${BIO_READ}" \
    --lencbc "${LEN_CBC}" \
    --lenumi "${LEN_UMI}" \
    $([[ "${UMI_FIRST}" == "true" ]] && echo "--umifirst") \
    --cbcfile "${CBC_FILE}" \
    --cbchd "${CBC_HD}" \
    --out "${OUT_PREFIX}" \
    2>&1 | tee "${LOG_DIR}/01_concatenator.log"

  # 2) Compress for trimming (concatenator outputs *_cbc.fastq)
  step "(${SAMPLE}) Compressing barcode FASTQ"
  gzip -f "${OUT_PREFIX}_cbc.fastq"
  local CONCAT_OUT_GZ="${OUT_PREFIX}_cbc.fastq.gz"

  # 3) Trim adapters/low-quality bases
  step "(${SAMPLE}) Trimming adapters and low-quality bases with Trim Galore"
  "${TRIM_GALORE}" \
    --path_to_cutadapt "${CUTADAPT_BIN}" \
    --output_dir "${SAMPLE_DIR}" \
    "${CONCAT_OUT_GZ}" \
    2>&1 | tee "${LOG_DIR}/02_trim_galore.log"

  local TRIMMED="${OUT_PREFIX}_cbc_trimmed.fq.gz"
  if [[ ! -f "${TRIMMED}" ]]; then
    echo "ERROR: Expected Trim Galore output not found: ${TRIMMED}" | tee -a "${LOG_DIR}/02_trim_galore.log"
    return 1
  fi

  # 4) Map with STAR (STAR manages its own temp dir)
  step "(${SAMPLE}) Mapping reads to reference genome with STAR"
  "${STAR_BIN}" \
    --runThreadN "${THREADS}" \
    --genomeDir "${GENOME_INDEX}" \
    --readFilesIn "${TRIMMED}" \
    --readFilesCommand zcat \
    --outFileNamePrefix "${OUT_PREFIX}_" \
    --outSAMtype BAM SortedByCoordinate \
    --outSAMattributes All \
    --outSAMstrandField intronMotif \
    --outFilterMultimapNmax 1 \
    --quantMode GeneCounts \
    2>&1 | tee "${LOG_DIR}/03_star.log"

  # 5) Introns/Exons counting (your helper script)
  local BAM="${OUT_PREFIX}_Aligned.sortedByCoord.out.bam"
  local INTRONS_BED="${EXINT_DIR}/${EXINT_BASENAME}_introns.bed"
  local EXONS_BED="${EXINT_DIR}/${EXINT_BASENAME}_exons.bed"
  if [[ ! -s "$BAM" ]]; then
    echo "ERROR: Expected BAM not found: ${BAM}" | tee -a "${LOG_DIR}/03_star.log"
    return 1
  fi
  step "(${SAMPLE}) Computing intron/exon coverage"
  bash "${GET_INTRONSEXONS}" \
    "${BAM}" \
    "${INTRONS_BED}" \
    "${EXONS_BED}" \
    "${OUT_PREFIX}" \
    "${BEDTOOLS_BIN}" \
    "${SAMTOOLS_BIN}" \
    "${SCRIPTS_DIR}" \
    > "${LOG_DIR}/04_introns_exons.log" 2>&1

  banner "FINISHED sample: ${SAMPLE} → ${SAMPLE_DIR}"
  echo
}

run_multiqc () {
  if [[ "${RUN_MULTIQC}" != "true" ]]; then
    step "Skipping MultiQC (disabled)."
    return 0
  fi
  banner "Running MultiQC summary"
  local outdir="${RESULTS_DIR}/multiqc"
  mkdir -p "${outdir}"
  # Search the entire results dir for logs/metrics from all steps
  if [[ -n "${MULTIQC_CONFIG}" ]]; then
    "${MULTIQC_BIN}" "${RESULTS_DIR}" -o "${outdir}" -c "${MULTIQC_CONFIG}" 2>&1 | tee "${outdir}/multiqc.log"
  else
    "${MULTIQC_BIN}" "${RESULTS_DIR}" -o "${outdir}" 2>&1 | tee "${outdir}/multiqc.log"
  fi
  step "MultiQC report: ${outdir}/multiqc_report.html"
}

# -------- Tool sanity checks --------
check_tools

# -------- Read manifest & run --------
python3 - "${MANIFEST_CSV}" <<'PY' | while IFS=$'\t' read -r SAMPLE RAW_FQ_DIR CBC_FILE; do
import csv, sys
path = sys.argv[1]
with open(path, newline='') as f:
    r = csv.DictReader(f)
    required = {"sample","raw_fq_dir","cbc_file"}
    missing = required - set(r.fieldnames or [])
    if missing:
        raise SystemExit(f"Missing CSV columns: {', '.join(sorted(missing))}")
    for row in r:
        print("\t".join([row["sample"], row["raw_fq_dir"], row["cbc_file"]]))
PY
  run_sample "$SAMPLE" "$RAW_FQ_DIR" "$CBC_FILE"
done
