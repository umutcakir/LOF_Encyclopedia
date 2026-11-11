#!/bin/bash
#SBATCH -p soeding
#SBATCH -N 1
#SBATCH -n 24
#SBATCH -t 14-00:00:00
#SBATCH --mem=150G
#SBATCH -o slurm_logs/search_%j.out
#SBATCH -e slurm_logs/search_%j.err

# =============================================================================
# 01_search_sample.slurm  --  search one sample for matches to the 40 queries
#
#     sbatch 01_search_sample.slurm SRR350519
#
# Two stages: reduce 2.7 billion reads to a candidate set, then BLAST the
# candidates for exact HSP statistics. BLASTing every read directly would
# take weeks; the prefilter is what makes this tractable.
#
#   EXACT MODE (bbduk, k-mer)     -- used when identity is high enough that
#                                    a usable exact k-mer is guaranteed.
#                                    Provably lossless: no true match can
#                                    be missed.
#
#   SENSITIVE MODE (bbmap, align) -- used when the thresholds are too
#                                    permissive for any exact-seed guarantee.
#                                    Highly sensitive in practice, but NOT
#                                    provably exhaustive. See the note below.
# =============================================================================

set -uo pipefail

ACC="${1:?Usage: sbatch 01_search_sample.slurm <SRR_ACCESSION>}"

# =============================================================================
# CONFIG
# =============================================================================
WORKDIR="${SLURM_SUBMIT_DIR:-$PWD}"
FASTQDIR="${WORKDIR}/fastq"
REFDIR="${WORKDIR}/search/ref"
OUTDIR="${WORKDIR}/search/per_sample"
TMPBASE="${WORKDIR}/search/tmp"

# --- Match thresholds ------------------------------------------------------
# Set these to the MOST permissive values you might want. 03_collect.sh can
# filter up from here in minutes; filtering down means redoing the search.
#
# Read lengths in this dataset range 36-150 bp:
#     SRR350519             36 bp   <- binding constraint on MIN_ALNLEN
#     SRR12971xxx (seeds)   50 bp
#     SRR949xxx             50-51 bp per mate
#     others                100-150 bp per mate
MIN_PIDENT=95        # minimum percent identity of an HSP
MIN_ALNLEN=30        # minimum alignment length in nt
#
# 30 nt at 95% permits 1 mismatch, giving a guaranteed 14-mer -> exact mode,
# provably lossless. Chance hits are negligible at this length: a 14-mer
# recurs often, but a 30 nt match at 95% essentially does not occur by chance
# in a 270 Gb database, so what comes back is real transcript evidence.
#
# Read lengths here span 36-150 bp, so 30 nt means 83% coverage of a
# SRR350519 read but only 20% of a bud read. Use perc_read_coverage in
# all_matches.tsv if you need stringency comparable across samples.

# --- Prefilter selection (COMPUTED -- do not hand-edit) --------------------
# A match of length L at identity P has at most  m = floor(L*(100-P)/100)
# mismatches, which split it into m+1 exact segments. The longest is at
# least (L-m)/(m+1). An exact k-mer prefilter is lossless only if k fits
# in that guaranteed segment:
#
#     k_max = floor( (MIN_ALNLEN - m) / (m + 1) )
#
#   30 nt @ 100%  -> m=0, k=30   exact mode
#   30 nt @  97%  -> m=0, k=30   exact mode
#   30 nt @  95%  -> m=1, k=14   exact mode
#   30 nt @  90%  -> m=3, k=6    too small -> sensitive mode
#   30 nt @  85%  -> m=4, k=5    too small -> sensitive mode
#
# A 5-mer occurs by chance roughly every kilobase, so an exact filter at
# that size would pass essentially every read and save nothing.
MAXMM=$(( MIN_ALNLEN * (100 - MIN_PIDENT) / 100 ))
KMER=$(( (MIN_ALNLEN - MAXMM) / (MAXMM + 1) ))
[[ ${KMER} -gt 31 ]] && KMER=31

if [[ ${KMER} -ge 12 ]]; then
    MODE="exact"
else
    MODE="sensitive"
fi

# For sensitive mode: align with a margin below the reporting threshold so
# borderline hits are not lost at the prefilter stage. BLAST applies the
# real cutoff afterwards.
BBMAP_MINID=$(awk -v p=${MIN_PIDENT} 'BEGIN{v=(p-10)/100; if(v<0.70)v=0.70; printf "%.2f", v}')
BBMAP_K=9            # seed length; 8 is BBMap's minimum, 9 balances speed

# --- Low-complexity handling ----------------------------------------------
# NONE. No entropy filtering, no dust masking, no reference masking.
# Every read that shares a k-mer with the queries is kept and reported,
# including reads that match only on poly-A/T or simple-repeat sequence.
# This makes the search exhaustive over the full query sequence, UTRs
# included, at the cost of a large low-complexity artefact fraction
# (~17% of reads matched in the SRR350519 test vs 0.23% expected by chance).

THREADS=${SLURM_NTASKS:-24}

# =============================================================================
# Environment
# =============================================================================
echo "=========================================================="
echo "Sample    : ${ACC}"
echo "Job       : ${SLURM_JOB_ID:-interactive}  on $(hostname)"
echo "Started   : $(date)"
echo "Identity  : >= ${MIN_PIDENT}%"
echo "Min match : >= ${MIN_ALNLEN} nt  (allows ${MAXMM} mismatch(es))"
echo "Prefilter : ${MODE} mode"
if [[ "${MODE}" == "exact" ]]; then
echo "            bbduk k=${KMER}, lossless"
else
echo "            bbmap k=${BBMAP_K} minid=${BBMAP_MINID}, heuristic"
echo ""
echo "  NOTE: at ${MIN_PIDENT}% identity over ${MIN_ALNLEN} nt, no seed-based"
echo "  method can guarantee exhaustive recall at this scale. Results are"
echo "  highly sensitive but not provably complete. For guaranteed-complete"
echo "  results, raise identity to >=95%, which restores exact mode."
fi
echo "=========================================================="

if command -v conda >/dev/null 2>&1; then
    CONDA_BASE="$(conda info --base)"
else
    CONDA_BASE="${HOME}/miniconda3"
fi
# shellcheck disable=SC1091
source "${CONDA_BASE}/etc/profile.d/conda.sh"
conda activate sra_env || { echo "ERROR: cannot activate sra_env"; exit 1; }

NEED=(blastn seqtk)
if [[ "${MODE}" == "exact" ]]; then NEED+=(bbduk.sh); else NEED+=(bbmap.sh); fi
for tool in "${NEED[@]}"; do
    if ! command -v "${tool}" >/dev/null 2>&1; then
        echo "ERROR: ${tool} not found. Install with:"
        echo "    conda install -n sra_env -c bioconda bbmap blast seqtk"
        exit 1
    fi
done

mkdir -p "${OUTDIR}" "${TMPBASE}"
TMP="${TMPBASE}/${ACC}_${SLURM_JOB_ID:-$$}"
mkdir -p "${TMP}"
trap 'rm -rf "${TMP}"' EXIT INT TERM

# =============================================================================
# Gather this sample's FASTQ files
# NOTE the glob: ${ACC}*.fastq.gz picks up _1, _2 AND the unpaired file that
# --split-3 writes. SRR7589436 has all three, with >half its data unpaired --
# a _1/_2 glob would silently discard it.
# =============================================================================
mapfile -t FQ < <(find "${FASTQDIR}" -maxdepth 1 -name "${ACC}*.fastq.gz" | sort)
if [[ ${#FQ[@]} -eq 0 ]]; then
    echo "ERROR: no FASTQ files found for ${ACC} in ${FASTQDIR}"
    exit 1
fi
echo "Input files (${#FQ[@]}):"
printf '  %s\n' "${FQ[@]}"
echo ""

# =============================================================================
# STAGE 1 -- reduce to candidate reads
# =============================================================================
echo "--- Stage 1: prefilter (${MODE}) ---"
HITFQ="${TMP}/candidates.fq"
: > "${HITFQ}"

TOTAL_IN=0
for fq in "${FQ[@]}"; do
    base=$(basename "${fq}" .fastq.gz)
    echo "  ${base} ..."
    HITS="${TMP}/${base}.hits.fq"

    if [[ "${MODE}" == "exact" ]]; then
        # outm = reads sharing an exact k-mer with the reference.
        # rcomp=t checks both strands; hdist=0 means no mismatch in the seed.
        #
        # -da disables Java assertions. BBTools 39.80's vectorized FASTQ
        # parser throws a spurious "missing plus" AssertionError on files
        # where the third line carries the full header (which is what
        # fasterq-dump writes) rather than a bare "+". Both forms are valid
        # FASTQ. If this causes trouble, pin an older build instead:
        #     conda install -n sra_env -c bioconda 'bbmap=39.06'
        bbduk.sh \
            in="${fq}" \
            outm="${HITS}" \
            ref="${REFDIR}/transcripts.seed.fasta" \
            k=${KMER} hdist=0 rcomp=t mm=f \
            -da \
            threads=${THREADS} -Xmx120g overwrite=t \
            2> "${TMP}/${base}.filter.log"
    else
        # Align against the 40 transcripts and keep anything that maps.
        #   local=t     allow soft-clipping, so a 30 nt match inside a
        #               150 bp read is found rather than requiring the
        #               whole read to align
        #   ambig=all   a read matching several splice variants is kept
        #               for all of them; Igor asked for ALL matches
        #   nodisk=t    build the index in memory, not on BeeGFS
        bbmap.sh \
            in="${fq}" \
            outm="${HITS}" \
            ref="${REFDIR}/transcripts.seed.fasta" \
            nodisk=t local=t ambig=all \
            minid=${BBMAP_MINID} k=${BBMAP_K} \
            -da \
            threads=${THREADS} -Xmx120g overwrite=t \
            2> "${TMP}/${base}.filter.log"
    fi

    if [[ ! -s "${HITS}" ]]; then
        echo "    no candidates"
        rm -f "${HITS}"
        continue
    fi

    # Tag headers with provenance: @ACC|MATE|ORIGINAL_ID
    # MATE is a short label, not the filename -- for a single-end sample the
    # filename equals the accession, and SRA read IDs already start with the
    # accession, so using the filename here printed it three times over.
    case "${base}" in
        *_1) MATE=R1  ;;
        *_2) MATE=R2  ;;
        *)   if [[ ${#FQ[@]} -gt 1 ]]; then MATE=unp; else MATE=se; fi ;;
    esac
    awk -v acc="${ACC}" -v mate="${MATE}" '
        NR%4==1 { sub(/^@/,""); print "@" acc "|" mate "|" $0; next }
        { print }
    ' "${HITS}" >> "${HITFQ}"

    n=$(( $(wc -l < "${HITS}") / 4 ))
    echo "    ${n} candidate reads"
    TOTAL_IN=$(( TOTAL_IN + n ))
    rm -f "${HITS}"
done

echo "  total candidates: ${TOTAL_IN}"
echo ""

if [[ ${TOTAL_IN} -eq 0 ]]; then
    echo "No candidates in ${ACC}. Writing empty result."
    : > "${OUTDIR}/${ACC}.tsv"
    : > "${OUTDIR}/${ACC}.reads.fasta"
    {
        echo "sample=${ACC}"; echo "candidates=0"; echo "kept_hsps=0"
        echo "mode=${MODE}"; echo "date=$(date -Iseconds)"
    } > "${OUTDIR}/${ACC}.done"
    exit 0
fi

# =============================================================================
# STAGE 2 -- blastn the candidates for exact HSP statistics
# =============================================================================
echo "--- Stage 2: blastn ---"
CANDFA="${TMP}/candidates.fasta"
seqtk seq -A "${HITFQ}" > "${CANDFA}"

# query = READ, subject = TRANSCRIPT
# -dust no / -soft_masking false : never mask; we want matches everywhere
# -word_size 7 in sensitive mode : BLAST's own seed must shrink too, or it
#                                  misses the divergent matches bbmap found
# -max_target_seqs 40            : the 40 queries are splice variants of six
#                                  genes and overlap heavily, so one read
#                                  legitimately matches many
if [[ "${MODE}" == "exact" ]]; then
    WORD=11
else
    WORD=7
fi

blastn \
    -query "${CANDFA}" \
    -db "${REFDIR}/transcripts" \
    -out "${TMP}/blast.tsv" \
    -outfmt '6 qseqid sseqid pident length qlen slen qstart qend sstart send evalue bitscore sstrand' \
    -perc_identity ${MIN_PIDENT} \
    -dust no -soft_masking false \
    -word_size ${WORD} \
    -evalue 1000 \
    -max_target_seqs 40 \
    -num_threads ${THREADS} \
    2>> "${TMP}/blast.log"

RAW=$(wc -l < "${TMP}/blast.tsv")
echo "  raw HSPs: ${RAW}"

awk -F'\t' -v m=${MIN_ALNLEN} '$4 >= m' "${TMP}/blast.tsv" > "${OUTDIR}/${ACC}.tsv"
KEPT=$(wc -l < "${OUTDIR}/${ACC}.tsv")
echo "  HSPs >= ${MIN_ALNLEN} nt: ${KEPT}"

# 03_collect.sh needs full read sequences to build the per-transcript FASTAs
cp "${CANDFA}" "${OUTDIR}/${ACC}.reads.fasta"

if [[ "${MODE}" == "exact" ]]; then LOSSLESS=yes; else LOSSLESS=no; fi
{
    echo "sample=${ACC}"
    echo "input_files=${#FQ[@]}"
    echo "candidates=${TOTAL_IN}"
    echo "raw_hsps=${RAW}"
    echo "kept_hsps=${KEPT}"
    echo "mode=${MODE}"
    echo "kmer=${KMER}"
    echo "bbmap_minid=${BBMAP_MINID}"
    echo "blast_word_size=${WORD}"
    echo "max_mismatches_at_min_len=${MAXMM}"
    echo "min_pident=${MIN_PIDENT}"
    echo "min_alnlen=${MIN_ALNLEN}"
    echo "lossless=${LOSSLESS}"
    echo "date=$(date -Iseconds)"
} > "${OUTDIR}/${ACC}.done"

echo ""
echo "=========================================================="
echo "${ACC} finished: $(date)"
echo "  ${OUTDIR}/${ACC}.tsv           (${KEPT} HSPs)"
echo "  ${OUTDIR}/${ACC}.reads.fasta   (${TOTAL_IN} reads)"
echo "=========================================================="
