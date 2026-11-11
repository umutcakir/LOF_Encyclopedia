#!/bin/bash
#SBATCH -p soeding
#SBATCH -N 1
#SBATCH -n 18
#SBATCH -t 2-00:00:00
#SBATCH --mem=350G
#SBATCH -J sra_download
#SBATCH -o slurm_logs/download_%j.out
#SBATCH -e slurm_logs/download_%j.err
#SBATCH --mail-type=END,FAIL

# =============================================================================
# download_fastq.slurm
#
# Downloads the 50 RNA-Seq runs from Supplementary Dataset S17 as gzipped
# FASTQ, using parallel workers.
#
# Submit with:
#     mkdir -p slurm_logs
#     sbatch download_fastq.slurm
#
# Monitor with:
#     squeue -u $USER
#     tail -f slurm_logs/download_<jobid>.out
#     ./download_fastq.slurm --check      # (runs locally, no job needed)
#
# Resumable: cancel and resubmit at any time; completed runs are skipped.
# =============================================================================

set -uo pipefail

# =============================================================================
# CONFIG
# =============================================================================
WORKDIR="${SLURM_SUBMIT_DIR:-$PWD}"
OUTDIR="${WORKDIR}/fastq"
LOGDIR="${WORKDIR}/logs_per_run"
SRA_CACHE="${WORKDIR}/sra_cache"

# Parallelism: downloads are network-limited, not CPU-limited.
# NCBI throttles heavily above ~12 concurrent connections, so more workers
# here makes things slower, not faster.
N_PARALLEL=4
THREADS_PER_JOB=4          # 4 x 4 = 16 cores

MAX_RETRIES=3

# --- Scratch location -----------------------------------------------------
# Everything stays under the working directory. Node-local /tmp and /dev/shm
# are deliberately NOT used: fasterq-dump needs several times the final output
# size as temp space, and with 12 workers running at once that easily exceeds
# a small /tmp or a RAM-disk cap, killing the job partway through.
#
# Set USE_SHM=true only if you have confirmed /dev/shm is large enough
# (check with: df -h /dev/shm). It is faster, but riskier.
USE_SHM=false

if [[ "${USE_SHM}" == "true" && -d /dev/shm && -w /dev/shm ]]; then
    SCRATCH="/dev/shm/sra_${SLURM_JOB_ID:-$$}"
else
    SCRATCH="${WORKDIR}/scratch_${SLURM_JOB_ID:-$$}"
fi

# Redirect anything that reads TMPDIR into our scratch too, so no stray
# tool writes into a node-local /tmp we do not control.
export TMPDIR="${SCRATCH}"

# =============================================================================
# The 50 runs (accession|group)
# =============================================================================
SRR_LIST=(
    "SRR5740870|N10"  "SRR5740864|N10"  "SRR5740861|N10"
    "SRR949248|N10"   "SRR949249|N10"   "SRR949250|N10"
    "SRR949251|N10"   "SRR949252|N10"   "SRR949253|N10"
    "SRR5740862|N14"  "SRR5740866|N14"  "SRR5740869|N14"
    "SRR18299116|N14" "SRR18299115|N14" "SRR18299114|N14"
    "SRR11637999|B"   "SRR11637998|B"   "SRR11637997|B"
    "SRR350519|F"
    "SRR7589436|L"    "SRR6979196|L"    "SRR6979187|L"
    "SRR6979188|L"    "SRR6979192|L"    "SRR6979193|L"
    "SRR949254|R"     "SRR949255|R"     "SRR949256|R"
    "SRR949257|R"     "SRR949258|R"     "SRR949259|R"
    "SRR3997852|R"    "SRR3997853|R"    "SRR3997864|R"
    "SRR12971740|Se"  "SRR12971741|Se"  "SRR12971749|Se"
    "SRR12971750|Se"  "SRR12971753|Se"  "SRR12971754|Se"
    "SRR18944067|St"  "SRR18944066|St"
    "SRR10058814|W"   "SRR10058818|W"   "SRR10058822|W"
    "SRR2060975|W"
    "SRR2060979|RoD"
    "SRR15439239|PC"  "SRR15439238|PC"  "SRR15439237|PC"
)

# =============================================================================
# --check mode (run directly, not via sbatch)
# =============================================================================
if [[ "${1:-}" == "--check" ]]; then
    done_n=0; missing_n=0
    printf "%-14s %-6s %-10s %s\n" "ACCESSION" "GROUP" "STATUS" "SIZE"
    for entry in "${SRR_LIST[@]}"; do
        srr="${entry%%|*}"; grp="${entry##*|}"
        if [[ -f "${OUTDIR}/${srr}.done" ]]; then
            sz=$(du -shc "${OUTDIR}/${srr}"*.fastq.gz 2>/dev/null | tail -1 | cut -f1)
            printf "%-14s %-6s %-10s %s\n" "${srr}" "${grp}" "DONE" "${sz}"
            ((done_n++))
        else
            printf "%-14s %-6s %-10s\n" "${srr}" "${grp}" "MISSING"
            ((missing_n++))
        fi
    done
    echo ""
    echo "Complete: ${done_n}/${#SRR_LIST[@]}   Missing: ${missing_n}"
    [[ -d "${OUTDIR}" ]] && echo "Total on disk: $(du -sh "${OUTDIR}" 2>/dev/null | cut -f1)"
    exit 0
fi

# =============================================================================
# Environment setup
# =============================================================================
echo "=========================================================="
echo "Job ID       : ${SLURM_JOB_ID:-interactive}"
echo "Node         : $(hostname)"
echo "Started      : $(date)"
echo "Work dir     : ${WORKDIR}"
echo "Output dir   : ${OUTDIR}"
echo "Scratch      : ${SCRATCH}"
echo "Parallelism  : ${N_PARALLEL} workers x ${THREADS_PER_JOB} threads"
echo "=========================================================="

# --- Activate conda environment -------------------------------------------
# `conda activate` is a shell function that does not exist in a
# non-interactive batch shell, so the conda profile must be sourced first.
# Without this you get: "CommandNotFoundError: Your shell has not been
# properly configured to use 'conda activate'".
CONDA_ENV="sra_env"

if command -v conda >/dev/null 2>&1; then
    CONDA_BASE="$(conda info --base)"
elif [[ -f "${HOME}/miniconda3/etc/profile.d/conda.sh" ]]; then
    CONDA_BASE="${HOME}/miniconda3"
elif [[ -f "${HOME}/anaconda3/etc/profile.d/conda.sh" ]]; then
    CONDA_BASE="${HOME}/anaconda3"
else
    echo "ERROR: could not locate a conda installation."
    echo "       Set CONDA_BASE manually near the top of this script."
    exit 1
fi

# shellcheck disable=SC1091
source "${CONDA_BASE}/etc/profile.d/conda.sh"

if ! conda activate "${CONDA_ENV}"; then
    echo "ERROR: failed to activate conda environment '${CONDA_ENV}'."
    echo "       Available environments:"
    conda env list
    exit 1
fi

echo "conda env    : ${CONDA_DEFAULT_ENV}"
echo "conda prefix : ${CONDA_PREFIX}"

if ! command -v prefetch >/dev/null 2>&1 || ! command -v fasterq-dump >/dev/null 2>&1; then
    echo ""
    echo "ERROR: prefetch / fasterq-dump not found even after activating '${CONDA_ENV}'."
    echo "       Install them with:"
    echo "           conda install -n ${CONDA_ENV} -c bioconda sra-tools"
    exit 1
fi

echo "prefetch     : $(command -v prefetch)"
echo "fasterq-dump : $(command -v fasterq-dump)"
echo "sra-tools    : $(fasterq-dump --version 2>&1 | tr -d '\n')"

# --- Point the SRA cache at our workdir, not $HOME -------------------------
# Default cache is ~/.ncbi and will blow through a home quota fast.
export NCBI_SETTINGS="${SRA_CACHE}/user-settings.mkfg"
mkdir -p "${SRA_CACHE}"
if [[ ! -f "${NCBI_SETTINGS}" ]]; then
    cat > "${NCBI_SETTINGS}" <<EOF
/config/default = "false"
/repository/user/main/public/root = "${SRA_CACHE}"
/repository/user/cache-disabled = "false"
EOF
fi

# --- Compressor -----------------------------------------------------------
if command -v pigz >/dev/null 2>&1; then
    ZIPCMD="pigz -p ${THREADS_PER_JOB}"
else
    ZIPCMD="gzip"
fi
echo "compressor   : ${ZIPCMD}"

# --- Directories ----------------------------------------------------------
mkdir -p "${OUTDIR}" "${LOGDIR}" "${SCRATCH}"

# Always clean RAM scratch on exit, including on scancel or timeout
cleanup() {
    echo "Cleaning scratch ${SCRATCH} ..."
    rm -rf "${SCRATCH}"
}
trap cleanup EXIT INT TERM

# --- Disk space checks ----------------------------------------------------
# Note: the human-readable line goes to stderr, so that capturing this
# function's stdout with $(...) does not swallow the display output.
report_space() {
    local label="$1" path="$2"
    local avail fs
    avail=$(df -BG "${path}" 2>/dev/null | awk 'NR==2 {gsub("G","",$4); print $4}')
    fs=$(df -h "${path}" 2>/dev/null | awk 'NR==2 {print $1" ("$2" total)"}')
    printf "  %-10s %-45s %s GB free\n" "${label}" "${fs}" "${avail}" >&2
    echo "${avail}"
}

echo "Disk space:"
AVAIL_OUT=$(report_space "output:"  "${OUTDIR}")
AVAIL_SCR=$(report_space "scratch:" "${SCRATCH}")
echo ""

# Output needs ~300-600 GB for 50 compressed runs.
if [[ "${AVAIL_OUT}" =~ ^[0-9]+$ ]] && [[ ${AVAIL_OUT} -lt 400 ]]; then
    echo "WARNING: under 400 GB free at output location."
fi

# Scratch holds uncompressed FASTQ for N_PARALLEL runs simultaneously.
# Budget roughly 40 GB per concurrent worker to be safe.
NEED_SCR=$(( N_PARALLEL * 40 ))
if [[ "${AVAIL_SCR}" =~ ^[0-9]+$ ]] && [[ ${AVAIL_SCR} -lt ${NEED_SCR} ]]; then
    echo "WARNING: scratch has ${AVAIL_SCR} GB but ${N_PARALLEL} parallel workers"
    echo "         need roughly ${NEED_SCR} GB. Reducing N_PARALLEL to fit."
    N_PARALLEL=$(( AVAIL_SCR / 40 ))
    [[ ${N_PARALLEL} -lt 1 ]] && N_PARALLEL=1
    THREADS_PER_JOB=$(( 96 / N_PARALLEL ))
    [[ ${THREADS_PER_JOB} -gt 16 ]] && THREADS_PER_JOB=16
    echo "         -> now ${N_PARALLEL} workers x ${THREADS_PER_JOB} threads"
fi
echo ""

# =============================================================================
# Worker function — downloads one accession
# Exported so xargs subshells can call it.
# =============================================================================
download_one() {
    local entry="$1"
    local srr="${entry%%|*}"
    local grp="${entry##*|}"
    local runlog="${LOGDIR}/${srr}.log"
    local tmp="${SCRATCH}/${srr}"
    local attempt=1

    if [[ -f "${OUTDIR}/${srr}.done" ]]; then
        echo "[SKIP ] ${srr} (${grp}) already complete"
        return 0
    fi

    while [[ ${attempt} -le ${MAX_RETRIES} ]]; do
        echo "[START] ${srr} (${grp}) attempt ${attempt}/${MAX_RETRIES}"
        local t0=$SECONDS

        rm -rf "${tmp}"; mkdir -p "${tmp}"
        rm -f "${OUTDIR}/${srr}"*.fastq.gz

        # --- prefetch the .sra container ---
        if ! prefetch "${srr}" \
                --output-directory "${SRA_CACHE}" \
                --max-size u \
                --progress \
                >> "${runlog}" 2>&1; then
            echo "[WARN ] ${srr} prefetch failed (attempt ${attempt})"
            ((attempt++)); sleep 60; continue
        fi

        # --- validate the .sra before converting ---
        # prefetch already does an md5 check, but vdb-validate additionally
        # verifies the internal column checksums and blob integrity. This is
        # the only point where a meaningful checksum comparison against NCBI
        # is possible -- once converted to FASTQ, no upstream reference exists.
        if ! vdb-validate "${SRA_CACHE}/${srr}" >> "${runlog}" 2>&1; then
            echo "[WARN ] ${srr} FAILED vdb-validate - corrupt download (attempt ${attempt})"
            rm -rf "${SRA_CACHE:?}/${srr}"
            ((attempt++)); sleep 30; continue
        fi
        echo "[VALID] ${srr} passed vdb-validate"

        # Record the .sra md5 now, while the file still exists, so the
        # download can be audited later without re-downloading.
        SRA_FILE=$(find "${SRA_CACHE}/${srr}" -name '*.sra' | head -1)
        SRA_MD5=""
        SRA_BYTES=""
        if [[ -n "${SRA_FILE}" ]]; then
            SRA_MD5=$(md5sum "${SRA_FILE}" | awk '{print $1}')
            SRA_BYTES=$(stat -c%s "${SRA_FILE}")
        fi

        # --- convert to FASTQ ---
        # --split-3 handles both PAIRED and SINGLE layouts correctly,
        # which matters here since the 50 runs are a mix of both.
        if ! fasterq-dump "${SRA_CACHE}/${srr}" \
                --outdir "${tmp}" \
                --temp "${tmp}" \
                --threads "${THREADS_PER_JOB}" \
                --split-3 \
                --skip-technical \
                >> "${runlog}" 2>&1; then
            echo "[WARN ] ${srr} fasterq-dump failed (attempt ${attempt})"
            ((attempt++)); sleep 60; continue
        fi

        local nfq
        nfq=$(find "${tmp}" -maxdepth 1 -name '*.fastq' | wc -l)
        if [[ ${nfq} -eq 0 ]]; then
            echo "[WARN ] ${srr} produced no FASTQ files"
            ((attempt++)); continue
        fi

        # --- compress into final location ---
        local ok=1
        for fq in "${tmp}"/*.fastq; do
            if ! ${ZIPCMD} -c "${fq}" > "${OUTDIR}/$(basename "${fq}").gz" 2>> "${runlog}"; then
                echo "[WARN ] ${srr} compression failed"
                ok=0; break
            fi
        done
        [[ ${ok} -eq 0 ]] && { ((attempt++)); continue; }

        # --- record provenance, checksums, and read counts ------------------
        local firstfq nreads
        firstfq=$(find "${tmp}" -maxdepth 1 -name '*.fastq' | head -1)
        nreads=$(( $(wc -l < "${firstfq}") / 4 ))
        {
            echo "accession=${srr}"
            echo "group=${grp}"
            echo "n_fastq_files=${nfq}"
            echo "reads_first_file=${nreads}"
            echo "vdb_validate=PASS"
            echo "sra_md5=${SRA_MD5}"
            echo "sra_bytes=${SRA_BYTES}"
            echo "completed=$(date -Iseconds)"
            echo "slurm_job=${SLURM_JOB_ID:-NA}"
            echo "sratools=$(fasterq-dump --version 2>&1 | tr -d '\n')"
            # md5 of each delivered .fastq.gz, for detecting later corruption
            # or incomplete transfer. These are OUR checksums, not NCBI's --
            # they verify the files have not changed since download, and
            # cannot be compared against any upstream reference.
            for gz in "${OUTDIR}/${srr}"*.fastq.gz; do
                [[ -e "${gz}" ]] || continue
                echo "output_md5=$(md5sum "${gz}" | awk '{print $1}')  $(basename "${gz}")"
            done
        } > "${OUTDIR}/${srr}.done"

        rm -rf "${tmp}" "${SRA_CACHE}/${srr}"

        local sz elapsed
        sz=$(du -shc "${OUTDIR}/${srr}"*.fastq.gz 2>/dev/null | tail -1 | cut -f1)
        elapsed=$((SECONDS - t0))
        echo "[DONE ] ${srr} (${grp}) ${nfq} file(s), ${nreads} reads, ${sz}, ${elapsed}s"
        return 0
    done

    echo "[FAIL ] ${srr} (${grp}) after ${MAX_RETRIES} attempts - see ${runlog}"
    rm -rf "${tmp}"
    return 1
}

export -f download_one
export OUTDIR LOGDIR SCRATCH SRA_CACHE ZIPCMD THREADS_PER_JOB MAX_RETRIES
export NCBI_SETTINGS SLURM_JOB_ID

# =============================================================================
# Single-accession mode
# If accessions are given as arguments, download only those and exit.
# Useful for testing the pipeline or retrying individual failures:
#     ./download_fastq.slurm SRR350519
#     ./download_fastq.slurm SRR350519 SRR7589436
# =============================================================================
if [[ $# -gt 0 ]]; then
    echo "SINGLE-ACCESSION MODE: $# accession(s) requested"
    echo "(the other $(( ${#SRR_LIST[@]} - $# )) runs in the set will NOT be downloaded)"
    echo ""

    single_fail=0
    for want in "$@"; do
        # Look up the group label for this accession
        entry=""
        for e in "${SRR_LIST[@]}"; do
            if [[ "${e%%|*}" == "${want}" ]]; then
                entry="${e}"
                break
            fi
        done
        if [[ -z "${entry}" ]]; then
            echo "NOTE: ${want} is not in Supplementary Dataset S17; labelling it 'manual'"
            entry="${want}|manual"
        fi
        download_one "${entry}" || single_fail=1
    done

    echo ""
    echo "=========================================================="
    echo "Single-accession mode finished: $(date)"
    echo "Output: ${OUTDIR}"
    du -sh "${OUTDIR}" 2>/dev/null
    echo "=========================================================="
    exit ${single_fail}
fi

# =============================================================================
# Run all accessions in parallel
# =============================================================================
echo "Launching ${N_PARALLEL} parallel workers over ${#SRR_LIST[@]} accessions..."
echo ""

printf '%s\n' "${SRR_LIST[@]}" | \
    xargs -P "${N_PARALLEL}" -I{} bash -c 'download_one "$@"' _ {}

# =============================================================================
# Summary
# =============================================================================
echo ""
echo "=========================================================="
echo "Download phase finished: $(date)"

done_n=0; fail_n=0; failed=()
for entry in "${SRR_LIST[@]}"; do
    srr="${entry%%|*}"
    if [[ -f "${OUTDIR}/${srr}.done" ]]; then
        ((done_n++))
    else
        ((fail_n++)); failed+=("${srr}")
    fi
done

echo "Succeeded : ${done_n}/${#SRR_LIST[@]}"
echo "Failed    : ${fail_n}"
if [[ ${fail_n} -gt 0 ]]; then
    echo "Failed accessions: ${failed[*]}"
    echo ""
    echo "Per-run logs are in ${LOGDIR}/"
    echo "Resubmit this same script to retry only the failures."
fi
echo "Total size: $(du -sh "${OUTDIR}" 2>/dev/null | cut -f1)"
echo "Output    : ${OUTDIR}"
echo "=========================================================="

# Non-zero exit if anything failed, so SLURM marks the job appropriately
[[ ${fail_n} -eq 0 ]] || exit 1
exit 0
