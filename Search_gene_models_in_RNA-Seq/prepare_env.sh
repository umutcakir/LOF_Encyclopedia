# 1. Strips CRLF from the query FASTA (it arrives with Windows line endings,
#    which would otherwise become part of every transcript ID).
# 2. Builds a small BLAST nucleotide database from the 40 transcripts.
# 3. Reports query stats so the k-mer threshold can be chosen sensibly.

# Usage:  ./00_prepare.sh
# =============================================================================
set -euo pipefail

QUERY_RAW="40_documents_from_Models_for_RNA_Seq.fasta"
REFDIR="search/ref"
QUERY="${REFDIR}/transcripts.fasta"

mkdir -p "${REFDIR}"

# --- 1. Normalise line endings and wrap ------------------------------------
if [[ ! -f "${QUERY_RAW}" ]]; then
    echo "ERROR: ${QUERY_RAW} not found in $(pwd)"; exit 1
fi
tr -d '\r' < "${QUERY_RAW}" > "${QUERY}"
echo "Cleaned query FASTA -> ${QUERY}"

# --- 2. Sanity checks ------------------------------------------------------
NSEQ=$(grep -c '^>' "${QUERY}")
echo "Sequences: ${NSEQ}"
if [[ ${NSEQ} -ne 40 ]]; then
    echo "WARNING: expected 40 sequences, found ${NSEQ}"
fi

# Reject anything non-nucleotide (would silently break the k-mer prefilter)
BAD=$(grep -v '^>' "${QUERY}" | tr -d 'ACGTNacgtn\n' | wc -c)
if [[ ${BAD} -gt 0 ]]; then
    echo "ERROR: query contains ${BAD} non-ACGTN characters."
    echo "       This pipeline is nucleotide-only."
    exit 1
fi
echo "All sequences are nucleotide (ACGT only) -- OK"

# --- 3. Length distribution ------------------------------------------------
echo ""
echo "Query lengths:"
awk '/^>/{if(id)print n"\t"id; id=substr($0,2); n=0; next}
     {n+=length($0)}
     END{print n"\t"id}' "${QUERY}" | sort -n | \
    awk '{printf "  %6d  %s\n",$1,$2}' | (head -3; echo "  ..."; tail -3)

TOTAL=$(grep -v '^>' "${QUERY}" | tr -d '\n' | wc -c)
echo "  Total query bases: ${TOTAL}"

# --- 4. Seed reference for the prefilter ----------------------------------
# The prefilter seeds on this file. It is an exact copy of the query FASTA:
# NO low-complexity masking is applied, so poly-A/T tracts and simple repeats
# in the UTRs remain fully searchable and reads covering them are reported.
#
# Consequence to be aware of: the queries are ~65% AT and contain
# homopolymer runs, so a substantial fraction of reads will match on
# low-complexity sequence alone (~17% in the SRR350519 test). Those are real
# alignments, but they carry little information about WHICH transcript a
# read came from. Filter them downstream if needed -- see the note in
# 03_collect.sh.
cp "${QUERY}" "${REFDIR}/transcripts.seed.fasta"
echo ""
echo "Seed reference: ${REFDIR}/transcripts.seed.fasta (unmasked)"

# --- 5. Build BLAST database from the UNMASKED transcripts ----------------
# The database is the 40 transcripts; reads are streamed in as queries.
# This keeps the database ~100 kB instead of hundreds of GB.
if command -v makeblastdb >/dev/null 2>&1; then
    makeblastdb -in "${QUERY}" -dbtype nucl \
        -out "${REFDIR}/transcripts" \
        -title "Mt_transcript_models_40" \
        -parse_seqids > "${REFDIR}/makeblastdb.log" 2>&1
    echo ""
    echo "BLAST database built: ${REFDIR}/transcripts"
else
    echo ""
    echo "ERROR: makeblastdb not found. Install with:"
    echo "    conda install -n sra_env -c bioconda blast"
    exit 1
fi

# --- 6. Sample -> tissue group map ----------------------------------------
cat > "${REFDIR}/sample_groups.tsv" <<'EOF'
SRR5740870	10-day_nodules
SRR5740864	10-day_nodules
SRR5740861	10-day_nodules
SRR949248	10-day_nodules
SRR949249	10-day_nodules
SRR949250	10-day_nodules
SRR949251	10-day_nodules
SRR949252	10-day_nodules
SRR949253	10-day_nodules
SRR5740862	14-day_nodules
SRR5740866	14-day_nodules
SRR5740869	14-day_nodules
SRR18299116	14-day_nodules
SRR18299115	14-day_nodules
SRR18299114	14-day_nodules
SRR11637999	Buds
SRR11637998	Buds
SRR11637997	Buds
SRR350519	Flowers
SRR7589436	Leaves
SRR6979196	Leaves
SRR6979187	Leaves
SRR6979188	Leaves
SRR6979192	Leaves
SRR6979193	Leaves
SRR949254	Roots
SRR949255	Roots
SRR949256	Roots
SRR949257	Roots
SRR949258	Roots
SRR949259	Roots
SRR3997852	Roots
SRR3997853	Roots
SRR3997864	Roots
SRR12971740	Seeds
SRR12971741	Seeds
SRR12971749	Seeds
SRR12971750	Seeds
SRR12971753	Seeds
SRR12971754	Seeds
SRR18944067	Stems
SRR18944066	Stems
SRR10058814	Whole_plant
SRR10058818	Whole_plant
SRR10058822	Whole_plant
SRR2060975	Whole_plant
SRR2060979	Root_drought
SRR15439239	Petiole_control
SRR15439238	Petiole_control
SRR15439237	Petiole_control
EOF
echo "Sample-to-tissue map: ${REFDIR}/sample_groups.tsv"

echo ""
echo "=========================================================="
echo "Ready. Next:  ./02_submit_all.sh"
echo "=========================================================="
