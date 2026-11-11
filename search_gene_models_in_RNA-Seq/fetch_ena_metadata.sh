#!/usr/bin/env bash
# =============================================================================
# fetch_ena_metadata.sh
#
# Builds ena_metadata.tsv by querying ENA one accession at a time.
# The comma-separated multi-accession form of the filereport endpoint is
# unreliable and often returns only a header, so this loops instead.
#
# Falls back to NCBI's runinfo endpoint if ENA has no record for a run.
#
# Run this on the LOGIN NODE (compute nodes are usually firewalled), then:
#     sbatch verify_downloads.slurm --quick
#
# Usage:  ./fetch_ena_metadata.sh
# =============================================================================

set -uo pipefail

OUTDIR="${OUTDIR:-fastq}"
OUT="ena_metadata.tsv"
ENA="https://www.ebi.ac.uk/ena/portal/api/filereport"
NCBI="https://trace.ncbi.nlm.nih.gov/Traces/sra-db-be/runinfo"

# Collect accessions from the .done markers
mapfile -t ACCS < <(find "${OUTDIR}" -maxdepth 1 -name '*.done' -printf '%f\n' \
                    | sed 's/\.done$//' | sort)

if [[ ${#ACCS[@]} -eq 0 ]]; then
    echo "No .done files found in ${OUTDIR}/"
    exit 1
fi

echo "Fetching metadata for ${#ACCS[@]} accessions..."
printf 'run_accession\tread_count\tbase_count\tlibrary_layout\n' > "${OUT}"

n_ena=0; n_ncbi=0; n_miss=0; missing=()

for acc in "${ACCS[@]}"; do
    row=""

    # --- Try ENA first -----------------------------------------------------
    resp=$(curl -sS --max-time 30 --retry 2 \
        "${ENA}?accession=${acc}&result=read_run&fields=run_accession,read_count,base_count,library_layout&format=tsv" \
        2>/dev/null)
    # Take the first data line (skip header), require a non-empty read_count
    row=$(echo "${resp}" | awk -F'\t' 'NR>1 && $2 != "" {print; exit}')

    if [[ -n "${row}" ]]; then
        echo "${row}" >> "${OUT}"
        printf "  %-14s ENA   %s reads\n" "${acc}" "$(echo "${row}" | cut -f2)"
        ((n_ena++))
        continue
    fi

    # --- Fall back to NCBI runinfo -----------------------------------------
    # runinfo is CSV; spots = read count, bases = base count,
    # LibraryLayout = SINGLE/PAIRED
    resp=$(curl -sS --max-time 30 --retry 2 "${NCBI}?acc=${acc}" 2>/dev/null)
    row=$(echo "${resp}" | awk -F',' -v a="${acc}" '
        NR==1 {
            for (i=1; i<=NF; i++) {
                if ($i=="spots")         c_spots=i
                if ($i=="bases")         c_bases=i
                if ($i=="LibraryLayout") c_layout=i
            }
            next
        }
        $1==a && c_spots {
            printf "%s\t%s\t%s\t%s\n", a, $c_spots, $c_bases, toupper($c_layout)
            exit
        }')

    if [[ -n "${row}" ]]; then
        echo "${row}" >> "${OUT}"
        printf "  %-14s NCBI  %s reads\n" "${acc}" "$(echo "${row}" | cut -f2)"
        ((n_ncbi++))
    else
        printf "  %-14s NOT FOUND\n" "${acc}"
        ((n_miss++)); missing+=("${acc}")
    fi
done

echo ""
echo "=========================================================="
echo "Wrote ${OUT}"
echo "  from ENA  : ${n_ena}"
echo "  from NCBI : ${n_ncbi}"
echo "  missing   : ${n_miss}"
[[ ${n_miss} -gt 0 ]] && echo "  ${missing[*]}"
echo ""
echo "Next:  sbatch verify_downloads.slurm --quick"
echo "=========================================================="
