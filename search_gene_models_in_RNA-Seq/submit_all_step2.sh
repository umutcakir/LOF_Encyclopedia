#!/usr/bin/env bash
# =============================================================================
# 02_submit_all.sh  --  submit one search job per sample
#     ./02_submit_all.sh                 submit all 
#     ./02_submit_all.sh SRR350519       submit one
#     ./02_submit_all.sh --check         status only, submits nothing
#     ./02_submit_all.sh --collect       submit all, then chain 03_collect.sh
# =============================================================================
set -uo pipefail

OUTDIR="search/per_sample"
FASTQDIR="fastq"
mkdir -p slurm_logs "${OUTDIR}"

if [[ ! -f search/ref/transcripts.nin && ! -f search/ref/transcripts.nal ]]; then
    echo "ERROR: BLAST database missing. Run ./00_prepare.sh first."
    exit 1
fi

# Sample list comes from the downloads that actually completed
mapfile -t ALL < <(find "${FASTQDIR}" -maxdepth 1 -name '*.done' -printf '%f\n' \
                   | sed 's/\.done$//' | sort)

CHAIN=false
TARGETS=()
for a in "$@"; do
    case "${a}" in
        --check)   MODE=check ;;
        --collect) CHAIN=true ;;
        *)         TARGETS+=("${a}") ;;
    esac
done
[[ ${#TARGETS[@]} -eq 0 ]] && TARGETS=("${ALL[@]}")

# --- status only -----------------------------------------------------------
if [[ "${MODE:-}" == "check" ]]; then
    ndone=0; npend=0
    printf "%-14s %-8s %s\n" "SAMPLE" "STATUS" "HSPs"
    for acc in "${ALL[@]}"; do
        if [[ -f "${OUTDIR}/${acc}.done" ]]; then
            h=$(grep '^kept_hsps=' "${OUTDIR}/${acc}.done" | cut -d= -f2)
            printf "%-14s %-8s %s\n" "${acc}" "DONE" "${h}"
            ((ndone++))
        else
            printf "%-14s %-8s\n" "${acc}" "PENDING"
            ((npend++))
        fi
    done
    echo ""
    echo "Done: ${ndone}/${#ALL[@]}   Pending: ${npend}"
    exit 0
fi

# --- submit ----------------------------------------------------------------
JOBIDS=()
nsub=0; nskip=0
for acc in "${TARGETS[@]}"; do
    if [[ -f "${OUTDIR}/${acc}.done" ]]; then
        echo "SKIP   ${acc} (already complete)"
        ((nskip++)); continue
    fi
    jid=$(sbatch --parsable -J "srch_${acc}" 01_search_sample.slurm "${acc}")
    echo "SUBMIT ${acc}  ->  job ${jid}"
    JOBIDS+=("${jid}")
    ((nsub++))
done

echo ""
echo "Submitted ${nsub}, skipped ${nskip}"

# --- optionally chain the collection step ---------------------------------
if [[ "${CHAIN}" == "true" && ${#JOBIDS[@]} -gt 0 ]]; then
    DEP=$(IFS=:; echo "${JOBIDS[*]}")
    cid=$(sbatch --parsable --dependency=afterok:"${DEP}" \
                 -J collect -p soeding -N 1 -n 8 -t 14-00:00:00 --mem=64G \
                 -o slurm_logs/collect_%j.out -e slurm_logs/collect_%j.err \
                 03_collect.sh)
    echo "Collection job ${cid} will run once all searches succeed."
else
    echo "When all jobs finish, run:  ./03_collect.sh"
fi


