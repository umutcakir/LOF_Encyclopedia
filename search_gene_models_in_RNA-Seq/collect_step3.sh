#!/usr/bin/env bash
#SBATCH -p soeding
#SBATCH -N 1
#SBATCH -n 8
#SBATCH -t 14-00:00:00
#SBATCH --mem=64G
#SBATCH -J collect
#SBATCH -o slurm_logs/collect_%j.out
#SBATCH -e slurm_logs/collect_%j.err

# =============================================================================
# 03_collect.sh  --  merge per-sample results
#
# Produces:
#   search/results/all_matches.tsv            one row per HSP, all samples
#   search/results/per_transcript/*.fasta     one FASTA per transcript,
#                                             ready to load into Geneious
#   search/results/summary_by_transcript.tsv  read counts per transcript
#   search/results/threshold_scan.tsv         hits vs. minimum match length
#
# Usage:  ./03_collect.sh            (or sbatch 03_collect.sh)
#         MIN_LEN=50 ./03_collect.sh  to change the length cutoff
# =============================================================================
set -uo pipefail

WORKDIR="${SLURM_SUBMIT_DIR:-$PWD}"
cd "${WORKDIR}"

PERSAMPLE="search/per_sample"
REFDIR="search/ref"
RESULTS="search/results"
FADIR="${RESULTS}/per_transcript"

# Reporting thresholds. These filter BLAST output that already exists, so
# changing them and re-running takes minutes -- no need to redo the search,
# PROVIDED the new values are not more permissive than what 01_search used
# (MIN_ALNLEN=32, MIN_PIDENT=95). Going below those gives incomplete results,
# because the k-mer prefilter was only lossless down to that point.
#
# Note that read length varies 36-150 bp across these samples, so a fixed
# length cutoff means different fractional coverage per sample. At 32 nt:
#   36 bp reads  -> 89% of the read
#   50 bp reads  -> 64%
#   150 bp reads -> 21%
# Use perc_read_coverage in all_matches.tsv if you need consistent stringency.
MIN_LEN="${MIN_LEN:-30}"
MIN_PID="${MIN_PID:-95}"

# Minimum percent of the READ covered by the alignment. This is the right
# knob for cutting result size, because it is length-normalised: a fixed
# nt cutoff like 150 would silently exclude the 36, 50 and 101 bp samples
# entirely and bias everything toward the deep 150 bp libraries.
#   0  = keep everything (default, matches what Igor asked for)
#   80 = reads aligning over most of their length -- far smaller files
MIN_COV="${MIN_COV:-0}"

mkdir -p "${RESULTS}" "${FADIR}"

echo "=========================================================="
echo "Collecting results"
echo "  min alignment length : ${MIN_LEN} nt"
echo "  min identity         : ${MIN_PID}%"
echo "  min read coverage    : ${MIN_COV}%"
echo "=========================================================="

mapfile -t SAMPLES < <(find "${PERSAMPLE}" -name '*.done' -printf '%f\n' \
                       | sed 's/\.done$//' | sort)
if [[ ${#SAMPLES[@]} -eq 0 ]]; then
    echo "ERROR: no completed samples in ${PERSAMPLE}"; exit 1
fi
echo "Samples with results: ${#SAMPLES[@]}"

python3 - "${MIN_LEN}" "${MIN_PID}" "${MIN_COV}" <<'PYEOF'
import os, sys, csv
from collections import defaultdict

MIN_LEN = int(sys.argv[1])
MIN_PID = float(sys.argv[2])
MIN_COV = float(sys.argv[3]) if len(sys.argv) > 3 else 0.0

PERSAMPLE = "search/per_sample"
REFDIR    = "search/ref"
RESULTS   = "search/results"
FADIR     = os.path.join(RESULTS, "per_transcript")

# --- sample -> tissue -------------------------------------------------------
groups = {}
with open(os.path.join(REFDIR, "sample_groups.tsv")) as fh:
    for line in fh:
        if line.strip():
            a, g = line.rstrip("\n").split("\t")
            groups[a] = g

# --- transcript list (preserve file order for stable output) ---------------
transcripts = []
with open(os.path.join(REFDIR, "transcripts.fasta")) as fh:
    for line in fh:
        if line.startswith(">"):
            transcripts.append(line[1:].strip().split()[0])

samples = sorted(
    f[:-5] for f in os.listdir(PERSAMPLE) if f.endswith(".done")
)

# --- load candidate read sequences -----------------------------------------
# Headers are  ACC|SOURCEFILE|ORIGINAL_ID
print("Loading candidate reads ...")
readseq = {}
for acc in samples:
    fa = os.path.join(PERSAMPLE, f"{acc}.reads.fasta")
    if not os.path.exists(fa):
        continue
    rid, buf = None, []
    with open(fa) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if rid:
                    readseq[rid] = "".join(buf)
                rid = line[1:].split()[0]
                buf = []
            else:
                buf.append(line)
        if rid:
            readseq[rid] = "".join(buf)
print(f"  {len(readseq):,} candidate reads loaded")

# --- sequence complexity helper --------------------------------------------
# Nothing is filtered on this: it is reported as a column so low-complexity
# matches can be separated downstream without redoing the search. Shannon
# entropy over 4-mers, normalised 0-1. Poly-A and simple repeats score near
# 0; ordinary coding sequence scores 0.8-1.0.
import math
def complexity(seq, k=4):
    if len(seq) < k:
        return 0.0
    counts = {}
    total = 0
    for i in range(len(seq) - k + 1):
        km = seq[i:i+k]
        counts[km] = counts.get(km, 0) + 1
        total += 1
    if total == 0:
        return 0.0
    h = -sum((c/total) * math.log2(c/total) for c in counts.values())
    hmax = math.log2(min(total, 4**k))
    return h / hmax if hmax > 0 else 0.0

# --- parse HSPs -------------------------------------------------------------
print("Parsing BLAST results ...")
rows = []
per_transcript = defaultdict(dict)      # transcript -> {read_id: seq}
length_hist    = defaultdict(int)       # alignment length -> count
per_tr_sample  = defaultdict(lambda: defaultdict(set))  # tr -> sample -> reads

n_raw = 0
for acc in samples:
    tsv = os.path.join(PERSAMPLE, f"{acc}.tsv")
    if not os.path.exists(tsv):
        continue
    with open(tsv) as fh:
        for row in csv.reader(fh, delimiter="\t"):
            if len(row) < 12:
                continue
            n_raw += 1
            (qseqid, sseqid, pident, length, qlen, slen,
             qstart, qend, sstart, send, evalue, bitscore) = row[:12]
            strand = row[12] if len(row) > 12 else ""

            aln  = int(length)
            pid  = float(pident)
            length_hist[(aln, pid)] += 1

            # Read span, not alignment length: BLAST counts gap positions in
            # `length`, which can exceed the read and give >100% coverage.
            span = abs(int(qend) - int(qstart)) + 1
            cov  = span / int(qlen) * 100.0

            if aln < MIN_LEN or pid < MIN_PID or cov < MIN_COV:
                continue

            # read id:  ACC|SOURCE|ORIGINAL
            parts   = qseqid.split("|", 2)
            sample  = parts[0] if parts else acc
            origid  = parts[2] if len(parts) > 2 else qseqid
            tissue  = groups.get(sample, "unknown")

            seq = readseq.get(qseqid, "")

            rows.append([
                sseqid,                       # transcript ID  (first column)
                origid,                       # read ID
                sample,                       # sample code
                tissue,                       # sample name / tissue
                aln,                          # length of match
                f"{pid:.2f}",                 # % nucleotide identity
                f"{cov:.2f}",                 # % coverage of the read
                qlen,                         # read length
                sstart, send, strand,         # position on the transcript
                qstart, qend,                 # position on the read
                evalue, bitscore,
                f"{complexity(seq):.3f}" if seq else "",  # 0=poly-A, 1=complex
                seq,                          # full read sequence
            ])

            # Geneious-friendly header: tissue|sample|read, e.g.
            #   Flowers|SRR350519|SRR350519.5934
            #   Roots|SRR949256|SRR949256.77/R1
            mate_sfx = f"/{parts[1]}" if len(parts) > 1 and parts[1] in ("R1","R2","unp") else ""
            per_transcript[sseqid][qseqid] = (
                f"{tissue}|{sample}|{origid}{mate_sfx}", seq)
            per_tr_sample[sseqid][sample].add(qseqid)

print(f"  {n_raw:,} HSPs read, {len(rows):,} passed filters")

# --- 1. all_matches.tsv -----------------------------------------------------
HEADER = ["transcript_id","read_id","sample_code","sample_name",
          "match_length_nt","perc_identity","perc_read_coverage","read_length",
          "transcript_start","transcript_end","strand",
          "read_start","read_end","evalue","bitscore",
          "read_complexity","read_sequence"]

out = os.path.join(RESULTS, "all_matches.tsv")
rows.sort(key=lambda r: (r[0], r[2], -int(r[4])))
with open(out, "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t")
    w.writerow(HEADER)
    w.writerows(rows)
print(f"  -> {out}")

# --- 2. per-transcript FASTA ------------------------------------------------
# One file per transcript, every matching read, deduplicated.
# This is what gets loaded into Geneious.
print("Writing per-transcript FASTA files ...")
GENEIOUS_WARN = 50000     # reads above which manual analysis stops being viable
oversized = []
for tr in transcripts:
    reads = per_transcript.get(tr, {})
    safe  = tr.replace("/", "_").replace(" ", "_")
    path  = os.path.join(FADIR, f"{safe}.fasta")
    with open(path, "w") as fh:
        # Values are (display_header, sequence) tuples, built during parsing.
        # The header is already TISSUE|ACCESSION|READ_ID[/R1|/R2|/unp], so
        # nothing needs rebuilding here -- unpack and write.
        for rid, (label, seq) in reads.items():
            fh.write(f">{label}\n{seq}\n")
    n = len(reads)
    mark = "  <-- very large" if n > GENEIOUS_WARN else ""
    if n > GENEIOUS_WARN:
        oversized.append((tr, n))
    print(f"  {safe}.fasta  {n} reads{mark}")

if oversized:
    print("")
    print(f"  WARNING: {len(oversized)} transcript file(s) exceed "
          f"{GENEIOUS_WARN:,} reads.")
    print(f"  At MIN_LEN={MIN_LEN} nt, most of these are chance matches rather")
    print( "  than real transcript evidence -- a 15-mer recurs about every")
    print( "  270 Mb, and this database is ~270 Gb. Files this size are not")
    print( "  practical to inspect by hand in Geneious.")
    print( "  Re-run with a longer minimum, e.g.:")
    print( "      MIN_LEN=30 ./03_collect.sh")
    print( "  and compare against threshold_scan.tsv before deciding.")

# --- 3. summary by transcript ----------------------------------------------
out = os.path.join(RESULTS, "summary_by_transcript.tsv")
with open(out, "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t")
    w.writerow(["transcript_id","total_reads","n_samples_with_hits",
                "tissues_with_hits"])
    for tr in transcripts:
        smap    = per_tr_sample.get(tr, {})
        tissues = sorted({groups.get(s,"unknown") for s in smap})
        w.writerow([tr, len(per_transcript.get(tr, {})), len(smap),
                    ",".join(tissues) if tissues else "-"])
print(f"  -> {out}")

# --- 4. threshold scan ------------------------------------------------------
# How many HSPs survive at each combination of minimum length and identity?
# This is the evidence for choosing cutoffs rather than guessing at them.
#
# IMPORTANT: BLAST was run with -perc_identity set at search time, so no HSP
# below that value exists in the data. Rows below it are marked INCOMPLETE
# rather than omitted, because a bare number there would look like a real
# count of zero additional hits when it actually means "never searched for".
search_pid = None
for acc in samples:
    d = os.path.join(PERSAMPLE, f"{acc}.done")
    if os.path.exists(d):
        for line in open(d):
            if line.startswith("min_pident="):
                search_pid = float(line.strip().split("=")[1])
                break
    if search_pid is not None:
        break

out = os.path.join(RESULTS, "threshold_scan.tsv")
with open(out, "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t")
    w.writerow(["min_alignment_length","min_identity",
                "n_HSPs_at_or_above","completeness"])
    for pid_thr in (95.0, 97.0, 99.0, 100.0):
        for len_thr in (15, 20, 25, 30, 40, 50, 60, 75, 100, 125, 150):
            n = sum(c for (L, P), c in length_hist.items()
                    if L >= len_thr and P >= pid_thr)
            if search_pid is not None and pid_thr < search_pid:
                flag = f"INCOMPLETE_searched_at_{search_pid:.0f}pct"
            else:
                flag = "complete"
            w.writerow([len_thr, f"{pid_thr:.0f}", n, flag])
print(f"  -> {out}")
if search_pid is not None:
    print(f"     (search ran at {search_pid:.0f}% identity; "
          f"lower rows are marked INCOMPLETE)")
PYEOF

echo ""
echo "=========================================================="
echo "Deliverables:"
echo "  ${RESULTS}/all_matches.tsv"
echo "  ${FADIR}/                    (40 FASTA files for Geneious)"
echo "  ${RESULTS}/summary_by_transcript.tsv"
echo "  ${RESULTS}/threshold_scan.tsv"
