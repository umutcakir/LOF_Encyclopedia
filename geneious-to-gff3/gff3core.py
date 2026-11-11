"""
gff3core.py - the logic behind the Geneious -> GFF3 tool (no web code in here).

    read_geneious()  finds annotated nucleotide sequences inside a .geneious file
    build_model()    turns the 'exon' and 'CDS' annotations of one sequence into a
                     transcript model (exons, CDS parts with phases, UTRs) and checks it.
                     A sequence with exons but no CDS annotation becomes a non-coding
                     transcript: exons only, no CDS, no UTR features and no protein.
    Placement / overlap_offsets() / placement_from_offset()
                     put a sequence onto the chromosome (typed coordinates, or found
                     automatically by matching a splice variant to another transcript)
    build_gene()     writes one gene with any number of mRNAs as GFF3 + FASTA + report
    merge_gff3()     cleans, validates, sorts and merges several GFF3 files

Everything this module produces uses Unix (LF) line endings.
"""
from __future__ import annotations

import gzip
import io
import re
import zipfile
import xml.etree.ElementTree as ET
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

MAX_UNZIPPED_BYTES = 500 * 1024 * 1024
MIN_MATCH_OVERLAP = 200   # identical overlap (nt) required to place a splice variant automatically
SHORT_INTRON = 30         # introns shorter than this are flagged

ID_OK = re.compile(r"^[A-Za-z0-9._:|\-]+$")
SEQID_OK = re.compile(r"^[A-Za-z0-9.:^*$@!+_?|\-]+$")


# --------------------------------------------------------------------------- messages
class Issues(list):
    """A list of (level, text) tuples; level is 'error', 'warning' or 'info'."""

    def error(self, text: str) -> None:
        self.append(("error", text))

    def warning(self, text: str) -> None:
        self.append(("warning", text))

    def info(self, text: str) -> None:
        self.append(("info", text))

    def of(self, level: str) -> List[str]:
        return [t for lvl, t in self if lvl == level]

    @property
    def has_errors(self) -> bool:
        return any(lvl == "error" for lvl, _ in self)


# --------------------------------------------------------------------------- sequences
_COMP = str.maketrans("ACGTRYMKSWBDHVNacgtrymkswbdhvn", "TGCAYRKMSWVHDBNtgcayrkmswvhdbn")
_B = "TCAG"
_AA = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
CODONS = {a + b + c: _AA[16 * i + 4 * j + k]
          for i, a in enumerate(_B) for j, b in enumerate(_B) for k, c in enumerate(_B)}
STOP_CODONS = {"TAA", "TAG", "TGA"}
SPLICE_KINDS = (("GT", "AG", "GT-AG"), ("GC", "AG", "GC-AG"), ("AT", "AC", "AT-AC"))


def revcomp(seq: str) -> str:
    return seq.translate(_COMP)[::-1]


def translate(nt: str) -> str:
    nt = nt.upper()
    return "".join(CODONS.get(nt[i:i + 3], "X") for i in range(0, len(nt) - 2, 3))


def wrap(seq: str, width: int = 60) -> str:
    return "\n".join(seq[i:i + width] for i in range(0, len(seq), width))


# --------------------------------------------------------------------------- Geneious files
@dataclass
class GInterval:
    start: int
    end: int
    direction: str          # leftToRight, rightToLeft or none


@dataclass
class GAnnotation:
    type: str
    name: str
    intervals: List[GInterval]


@dataclass
class GRecord:
    """One nucleotide sequence document found in a .geneious file."""
    name: str
    seq: str
    annotations: List[GAnnotation]
    file_name: str

    def of_type(self, wanted: str) -> List[GAnnotation]:
        return [a for a in self.annotations if a.type.strip().lower() == wanted]

    def ignored_types(self) -> Dict[str, int]:
        counts: Dict[str, int] = {}
        for a in self.annotations:
            t = a.type.strip() or "(no type)"
            if t.lower() not in ("exon", "cds"):
                counts[t] = counts.get(t, 0) + 1
        return counts

    @property
    def is_annotated(self) -> bool:
        return bool(self.of_type("exon") or self.of_type("cds"))


def _xml_payloads(data: bytes):
    buf = io.BytesIO(data)
    if zipfile.is_zipfile(buf):
        with zipfile.ZipFile(buf) as z:
            total = 0
            for info in z.infolist():
                if info.is_dir():
                    continue
                total += info.file_size
                if total > MAX_UNZIPPED_BYTES:
                    raise ValueError("the archive is too large")
                yield z.read(info)
    elif data[:2] == b"\x1f\x8b":
        yield gzip.decompress(data)
    else:
        yield data


def read_geneious(data: bytes, file_name: str) -> Tuple[List[GRecord], Issues]:
    """Return every nucleotide sequence (with its annotations) stored in a .geneious file."""
    issues = Issues()
    records: List[GRecord] = []
    seen = set()
    xml_found = False
    try:
        payloads = list(_xml_payloads(data))
    except (zipfile.BadZipFile, OSError, ValueError, EOFError) as exc:
        issues.error(f"{file_name}: the file cannot be opened ({exc}).")
        return records, issues
    for payload in payloads:
        try:
            root = ET.fromstring(payload)
        except ET.ParseError:
            continue
        xml_found = True
        for el in root.iter():
            if el.tag == "gappedSequence":          # copies of sequences inside alignments
                continue
            cs = el.find("charSequence")
            if cs is None or not (cs.text or "").strip():
                continue
            seq = re.sub(r"\s+", "", cs.text).upper()
            if "-" in seq:
                continue
            name = (el.findtext("name") or "").strip() or f"sequence {len(records) + 1}"
            anns = []
            for a in el.findall("sequenceAnnotations/annotation"):
                ivs = []
                for iv in a.findall("intervals/interval"):
                    try:
                        s, e = int(iv.findtext("minimumIndex")), int(iv.findtext("maximumIndex"))
                    except (TypeError, ValueError):
                        continue
                    ivs.append(GInterval(min(s, e), max(s, e), (iv.findtext("direction") or "none").strip()))
                anns.append(GAnnotation((a.findtext("type") or "").strip(),
                                        (a.findtext("description") or "").strip(), ivs))
            key = (name, seq, tuple((a.type, tuple((i.start, i.end) for i in a.intervals)) for a in anns))
            if key not in seen:
                seen.add(key)
                records.append(GRecord(name, seq, anns, file_name))
    if not xml_found:
        issues.error(f"{file_name}: this does not look like a Geneious file (no Geneious XML inside).")
    elif not records:
        issues.error(f"{file_name}: no nucleotide sequence was found in this file.")
    else:
        bare = [r.name for r in records if not r.is_annotated]
        if bare:
            issues.info(f"{file_name}: sequence(s) without exon/CDS annotations were skipped: {', '.join(bare)}.")
        if not any(r.is_annotated for r in records):
            issues.error(f"{file_name}: none of its sequences carries annotations of type 'exon' or 'CDS'.")
    return records, issues


# --------------------------------------------------------------------------- transcript model
@dataclass
class Intron:
    number: int
    start: int
    end: int
    donor: str
    acceptor: str
    kind: str               # GT-AG, GC-AG, AT-AC or non-canonical
    hint: str = ""

    @property
    def length(self) -> int:
        return self.end - self.start + 1


@dataclass
class TranscriptModel:
    """Positions are 1-based on the sequence read 5'->3' along the transcript."""
    name: str
    file_name: str
    seq: str
    flipped: bool                                   # Geneious sequence was reverse-complemented
    exons: List[Tuple[int, int]]
    cds_start: Optional[int]
    cds_end: Optional[int]
    cds_parts: List[Tuple[int, int, int, int]] = field(default_factory=list)   # start, end, phase, exon no.
    utr5: List[Tuple[int, int, int]] = field(default_factory=list)             # start, end, exon no.
    utr3: List[Tuple[int, int, int]] = field(default_factory=list)
    introns: List[Intron] = field(default_factory=list)
    issues: Issues = field(default_factory=Issues)

    @property
    def length(self) -> int:
        return len(self.seq)

    @property
    def coding(self) -> bool:
        """False when the sequence carried no CDS annotation (a non-coding transcript)."""
        return self.cds_start is not None

    @property
    def mrna(self) -> str:
        return "".join(self.seq[s - 1:e] for s, e in self.exons)

    @property
    def cds(self) -> str:
        return "".join(self.seq[s - 1:e] for s, e, _, _ in self.cds_parts)

    @property
    def protein(self) -> str:
        p = translate(self.cds)
        return p[:-1] if p.endswith("*") else p


def _intron_kind(seq: str, a: int, b: int) -> Tuple[str, str, str]:
    donor, acceptor = seq[a - 1:a + 1], seq[b - 2:b]
    for d, ac, kind in SPLICE_KINDS:
        if donor == d and acceptor == ac:
            return donor, acceptor, kind
    return donor, acceptor, "non-canonical"


def _canonical_shift(seq: str, left: Tuple[int, int], right: Tuple[int, int], max_shift: int = 30):
    """If the intron between two exons can slide to a canonical position without
    changing the spliced sequence, return (shift, kind); else None."""
    e1, s2 = left[1], right[0]
    options = []
    k = 1
    while k <= max_shift and s2 + k <= right[1] and seq[e1 + k - 1] == seq[s2 + k - 2]:
        options.append(k)
        k += 1
    k = 1
    while k <= max_shift and e1 - k >= left[0] and seq[e1 - k] == seq[s2 - k - 1]:
        options.append(-k)
        k += 1
    rank = {kind: i for i, (_, _, kind) in enumerate(SPLICE_KINDS)}
    best = None
    for sh in options:
        kind = _intron_kind(seq, e1 + sh + 1, s2 + sh - 1)[2]
        if kind in rank:
            cand = (rank[kind], abs(sh), sh, kind)
            best = cand if best is None or cand < best else best
    return (best[2], best[3]) if best else None


def _tx_pos(exons: List[Tuple[int, int]], p: int) -> int:
    """Position of sequence base p within the spliced transcript (1-based)."""
    n = 0
    for s, e in exons:
        if p > e:
            n += e - s + 1
        elif p >= s:
            return n + p - s + 1
    return n


def longest_orf(spliced: str) -> Tuple[int, int, int]:
    """Longest ATG-to-stop open reading frame in a spliced transcript, read on this strand.
    Returns (codons, start, end) as 1-based positions in the spliced transcript; (0, 0, 0) if none."""
    best = (0, 0, 0)
    for frame in range(3):
        i = frame
        while i + 2 < len(spliced):
            if spliced[i:i + 3] == "ATG":
                j = i
                while j + 2 < len(spliced):
                    if spliced[j:j + 3] in STOP_CODONS:
                        codons = (j + 3 - i) // 3
                        if codons > best[0]:
                            best = (codons, i + 1, j + 3)
                        i = j
                        break
                    j += 3
                else:
                    break
            i += 3
    return best


def _reads_correctly_reversed(seq, exons, cs, ce) -> bool:
    L = len(seq)
    rs = revcomp(seq)
    rex = sorted((L - e + 1, L - s + 1) for s, e in exons)
    rcs, rce = L - ce + 1, L - cs + 1
    cds = "".join(rs[max(s, rcs) - 1:min(e, rce)] for s, e in rex if max(s, rcs) <= min(e, rce))
    if not cds or len(cds) % 3:
        return False
    p = translate(cds)
    return cds.startswith("ATG") and p.endswith("*") and "*" not in p[:-1]


def build_model(rec: GRecord) -> TranscriptModel:
    """Build and check a transcript model from the 'exon' and 'CDS' annotations of one sequence."""
    iss = Issues()
    L = len(rec.seq)
    exon_anns, cds_anns = rec.of_type("exon"), rec.of_type("cds")
    ignored = rec.ignored_types()
    if ignored:
        iss.info("Ignored annotations (only the types 'exon' and 'CDS' are used): "
                 + ", ".join(f"{n} x {t}" for t, n in sorted(ignored.items())) + ".")
    if not exon_anns:
        iss.error("No annotation of type 'exon' was found on this sequence.")
    if not cds_anns:
        iss.info("No annotation of type 'CDS' was found, so this is treated as a non-coding transcript "
                 "(a lncRNA, for example). The GFF3 will hold the exons only: no CDS, no UTR features "
                 "and no protein. If the transcript does code for a protein, annotate its CDS in Geneious "
                 "and upload the file again.")

    directions = {iv.direction for a in exon_anns + cds_anns for iv in a.intervals} - {"none", ""}
    flipped = directions == {"rightToLeft"}
    if len(directions) > 1:
        iss.error("Some exon/CDS annotations point left-to-right and others right-to-left. "
                  "All of them must point in the direction of transcription.")
    seq = revcomp(rec.seq) if flipped else rec.seq
    if flipped:
        iss.info("The exon/CDS annotations point right-to-left, so the gene runs on the reverse strand of this "
                 "Geneious sequence. The tool reverse-complemented it internally; positions in this report are "
                 "counted from the 5' end of the transcript.")

    def norm(iv: GInterval) -> Tuple[int, int]:
        return (L - iv.end + 1, L - iv.start + 1) if flipped else (iv.start, iv.end)

    raw = []
    for a in exon_anns:
        if len(a.intervals) > 1:
            iss.info(f"The exon annotation '{a.name or 'unnamed'}' has {len(a.intervals)} parts; "
                     "each part is treated as one exon.")
        raw += [norm(iv) for iv in a.intervals]
    exons = sorted(set(raw))
    if len(exons) < len(raw):
        iss.warning("Some exons are annotated twice at identical positions; the duplicates were ignored.")
    for s, e in exons:
        if s < 1 or e > L:
            iss.error(f"Exon {s}-{e} lies outside the sequence (1-{L}).")
    for i in range(len(exons) - 1):
        (s1, e1), (s2, e2) = exons[i], exons[i + 1]
        if s2 <= e1:
            iss.error(f"Exons {s1}-{e1} and {s2}-{e2} overlap.")
        elif s2 == e1 + 1:
            iss.error(f"Exons {s1}-{e1} and {s2}-{e2} touch, with no intron between them. "
                      "If they are one exon, annotate it once.")
        elif s2 - e1 - 1 < SHORT_INTRON:
            iss.warning(f"The intron between exons {i + 1} and {i + 2} is only {s2 - e1 - 1} nt long.")

    cds_ivs = sorted(norm(iv) for a in cds_anns for iv in a.intervals)
    cs = cds_ivs[0][0] if cds_ivs else None
    ce = max(e for _, e in cds_ivs) if cds_ivs else None
    if len(cds_anns) > 1:
        iss.warning(f"{len(cds_anns)} CDS annotations were found; the CDS is taken to run from the first "
                    f"start ({cs}) to the last end ({ce}). Keeping a single CDS annotation avoids ambiguity.")

    model = TranscriptModel(rec.name, rec.file_name, seq, flipped, exons, cs, ce, issues=iss)
    if iss.has_errors or not exons:
        return model

    def exon_no(p: int) -> Optional[int]:
        return next((n for n, (s, e) in enumerate(exons, 1) if s <= p <= e), None)

    if model.coding:
        if exon_no(cs) is None:
            iss.error(f"The CDS starts at {cs}, which is not inside any exon.")
        if exon_no(ce) is None:
            iss.error(f"The CDS ends at {ce}, which is not inside any exon.")
        if iss.has_errors:
            return model

        done = 0
        for n, (s, e) in enumerate(exons, 1):
            lo, hi = max(s, cs), min(e, ce)
            if lo <= hi:
                model.cds_parts.append((lo, hi, (3 - done % 3) % 3, n))
                done += hi - lo + 1
            if s < cs:
                model.utr5.append((s, min(e, cs - 1), n))
            if e > ce:
                model.utr3.append((max(s, ce + 1), e, n))

        for a in cds_anns:
            if len(a.intervals) > 1 and sorted(norm(iv) for iv in a.intervals) != [(s, e) for s, e, _, _ in model.cds_parts]:
                iss.warning("The parts of the joined CDS annotation do not coincide with the exon borders. The GFF3 "
                            "uses the exons plus the CDS start and end; check the borders in Geneious.")
                break

    for n in range(1, len(exons)):
        left, right = exons[n - 1], exons[n]
        a, b = left[1] + 1, right[0] - 1
        donor, acceptor, kind = _intron_kind(seq, a, b)
        intron = Intron(n, a, b, donor, acceptor, kind)
        if kind == "non-canonical":
            msg = f"Intron {n} ({a}-{b}) has non-canonical borders {donor}...{acceptor}."
            alt = _canonical_shift(seq, left, right)
            if alt:
                sh, akind = alt
                intron.hint = (f"Moving it by {sh:+d} nt gives {akind} with an identical spliced sequence: "
                               f"exon {n} end = {left[1] + sh}, exon {n + 1} start = {right[0] + sh}.")
                msg += " " + intron.hint
            iss.warning(msg)
        elif kind != "GT-AG":
            iss.info(f"Intron {n} ({a}-{b}) is a minor-class {kind} intron: legitimate, but worth a second look.")
        model.introns.append(intron)

    if not model.coding:
        codons, o_start, o_end = longest_orf(model.mrna)
        if codons:
            iss.info(f"For reference only: the longest ATG-to-stop reading frame in this spliced transcript is "
                     f"{codons} codons ({codons - 1} amino acids), at positions {o_start}-{o_end} of the "
                     f"{len(model.mrna):,} nt transcript. It is not written to the GFF3.")
        if exons[0][0] > 1:
            iss.info(f"The sequence has {exons[0][0] - 1} nt before exon 1; they are not part of the transcript.")
        if exons[-1][1] < L:
            iss.info(f"The sequence has {L - exons[-1][1]} nt after the last exon; they are not part of the transcript.")
        return model

    cds = model.cds
    if len(cds) % 3:
        iss.error(f"The spliced CDS is {len(cds)} nt long, which is not a multiple of 3. "
                  "An exon border or the CDS start/end is off.")
    else:
        prot = translate(cds)
        trouble = False
        if cds[:3] != "ATG":
            iss.warning(f"The CDS starts with {cds[:3]}, not ATG.")
            trouble = True
        if not prot.endswith("*"):
            nxt = model.mrna[_tx_pos(exons, ce):_tx_pos(exons, ce) + 3]
            msg = f"The CDS does not end with a stop codon (its last codon is {cds[-3:]})."
            if nxt in STOP_CODONS:
                msg += (f" The next codon is {nxt}: the CDS annotation probably ends one codon early. "
                        "In GFF3 the CDS includes the stop codon, so extend the CDS annotation by 3 nt.")
            iss.warning(msg)
            trouble = True
        internal = [i + 1 for i, aa in enumerate(prot[:-1]) if aa == "*"]
        if internal:
            iss.error("In-frame stop codon(s) inside the CDS at codon "
                      + ", ".join(map(str, internal[:10])) + (" ..." if len(internal) > 10 else "")
                      + ". Check the exon borders.")
            trouble = True
        if trouble and _reads_correctly_reversed(seq, exons, cs, ce):
            iss.warning("This CDS reads correctly (ATG ... stop, no internal stops) only in the opposite "
                        "direction: the annotation arrows in Geneious probably point the wrong way.")

    if exons[0][0] > 1:
        iss.info(f"The sequence has {exons[0][0] - 1} nt before exon 1; they are not part of the transcript.")
    if exons[-1][1] < L:
        iss.info(f"The sequence has {L - exons[-1][1]} nt after the last exon; they are not part of the transcript.")
    return model


# --------------------------------------------------------------------------- chromosome placement
@dataclass
class Placement:
    start: int      # smaller chromosome coordinate of the whole Geneious sequence
    end: int        # larger chromosome coordinate
    how: str


def to_chrom(s: int, e: int, pl: Placement, strand: str) -> Tuple[int, int]:
    """Transcript-oriented sequence positions -> chromosome coordinates (start <= end)."""
    if strand == "+":
        return pl.start + s - 1, pl.start + e - 1
    return pl.end - e + 1, pl.end - s + 1


def overlap_offsets(ref: str, qry: str, k: int = 40, min_overlap: int = MIN_MATCH_OVERLAP) -> List[int]:
    """Offsets d such that qry position p corresponds to ref position p + d over an
    identical overlap of at least min_overlap nt (or the whole shorter sequence)."""
    k = min(k, len(ref), len(qry))
    cands = set()
    for needle, hay, sign in ((qry[:k], ref, 1), (ref[:k], qry, -1)):
        i = hay.find(needle)
        while i != -1:
            cands.add(sign * i)
            i = hay.find(needle, i + 1)
    need = min(min_overlap, len(ref), len(qry))
    good = []
    for d in sorted(cands):
        lo, hi = max(0, -d), min(len(qry), len(ref) - d)
        if hi - lo >= need and qry[lo:hi] == ref[lo + d:hi + d]:
            good.append(d)
    return good


def placement_from_offset(d: int, ref_pl: Placement, qry_len: int, strand: str, how: str) -> Placement:
    if strand == "+":
        start = ref_pl.start + d
        return Placement(start, start + qry_len - 1, how)
    end = ref_pl.end - d
    return Placement(end - qry_len + 1, end, how)


def offset_between(ref_pl: Placement, qry_pl: Placement, strand: str) -> int:
    return qry_pl.start - ref_pl.start if strand == "+" else ref_pl.end - qry_pl.end


def sequences_agree(ref: str, qry: str, d: int) -> Optional[bool]:
    """True/False = the sequences are identical/different where they overlap; None = no overlap."""
    lo, hi = max(0, -d), min(len(qry), len(ref) - d)
    if hi <= lo:
        return None
    return qry[lo:hi] == ref[lo + d:hi + d]


# --------------------------------------------------------------------------- GFF3 writing
# Feature types (Sequence Ontology terms) offered for the transcript line.
CODING_TX_TYPE = "mRNA"
NONCODING_TX_TYPES = (
    ("lnc_RNA", "long non-coding RNA (lncRNA): the usual choice for a spliced, non-coding transcript"),
    ("antisense_RNA", "transcribed from the opposite strand of another gene and overlapping it"),
    ("ncRNA", "non-coding RNA, when the class is not known or none of the others fits"),
    ("pre_miRNA", "primary transcript of a microRNA"),
    ("snoRNA", "small nucleolar RNA"),
    ("tRNA", "transfer RNA"),
    ("rRNA", "ribosomal RNA"),
    ("transcript", "a transcript of no stated class"),
)
GENE_TYPES = (
    ("gene", "accepted everywhere, for coding and non-coding genes alike"),
    ("ncRNA_gene", "says in the gene line itself that the gene is non-coding; some browsers prefer it"),
)
SO_TYPE_OK = re.compile(r"^[A-Za-z][A-Za-z0-9_]{1,40}$")
GENE_LINE_TYPES = frozenset([t for t, _ in GENE_TYPES] + ["pseudogene"])
TRANSCRIPT_TYPES = frozenset([CODING_TX_TYPE] + [t for t, _ in NONCODING_TX_TYPES]
                             + ["transcript", "miRNA", "snRNA", "lncRNA", "mRNA_region"])


@dataclass
class GeneSpec:
    seqid: str
    strand: str
    gene_id: str
    gene_name: str = ""
    alias: str = ""
    note: str = ""
    source: str = "manual_curation"
    gene_type: str = "gene"


@dataclass
class TranscriptEntry:
    model: TranscriptModel
    placement: Placement
    mrna_id: str
    mrna_name: str = ""
    tx_type: str = ""           # empty: mRNA for a coding model, lnc_RNA for a non-coding one

    @property
    def type(self) -> str:
        return self.tx_type or (CODING_TX_TYPE if self.model.coding else NONCODING_TX_TYPES[0][0])


@dataclass
class GeneResult:
    gff3: str
    rows: List[List[str]]
    proteins_fa: str
    cds_fa: str
    mrna_fa: str
    report: str
    transcripts: List[dict]
    issues: Issues
    span: Tuple[int, int]


def escape_value(text: str, keep_commas: bool = False) -> str:
    """URL-escape the characters GFF3 reserves inside attribute values."""
    out = []
    for ch in text:
        if ch in ";=&%" or (ch == "," and not keep_commas) or ord(ch) < 32:
            out.append("%{:02X}".format(ord(ch)))
        else:
            out.append(ch)
    return "".join(out)


def _attrs(pairs: List[Tuple[str, str]]) -> str:
    return ";".join(f"{k}={v}" for k, v in pairs if v)


def build_gene(spec: GeneSpec, entries: List[TranscriptEntry]) -> GeneResult:
    """Write one gene (one or more mRNAs) as GFF3 with chromosome coordinates."""
    issues = Issues()
    st = spec.strand
    blocks, spans, tables = [], [], []
    prot_fa, cds_fa, mrna_fa = [], [], []
    for ent in entries:
        m, pl, mid = ent.model, ent.placement, ent.mrna_id

        def conv(s, e, pl=pl):
            return to_chrom(s, e, pl, st)

        ex_chrom = [conv(s, e) for s, e in m.exons]
        t_start, t_end = min(a for a, _ in ex_chrom), max(b for _, b in ex_chrom)
        spans.append((t_start, t_end))
        base = [spec.seqid, spec.source]
        mrna_row = base + [ent.type, t_start, t_end, ".", st, ".",
                           _attrs([("ID", mid), ("Parent", spec.gene_id),
                                   ("Name", escape_value(ent.mrna_name) if ent.mrna_name else "")])]
        children = []   # (chromosome start of the exon, 0 = exon line / 1 = part, own start, row)
        for n, (a, b) in enumerate(ex_chrom, 1):
            children.append((a, 0, a, base + ["exon", a, b, ".", st, ".", f"ID={mid}.exon{n};Parent={mid}"]))
        for s, e, n in m.utr5:
            a, b = conv(s, e)
            children.append((ex_chrom[n - 1][0], 1, a,
                             base + ["five_prime_UTR", a, b, ".", st, ".", f"ID={mid}.utr5p{n};Parent={mid}"]))
        for s, e, ph, n in m.cds_parts:
            a, b = conv(s, e)
            children.append((ex_chrom[n - 1][0], 1, a,
                             base + ["CDS", a, b, ".", st, str(ph), f"ID={mid}.cds;Parent={mid}"]))
        for s, e, n in m.utr3:
            a, b = conv(s, e)
            children.append((ex_chrom[n - 1][0], 1, a,
                             base + ["three_prime_UTR", a, b, ".", st, ".", f"ID={mid}.utr3p{n};Parent={mid}"]))
        children.sort(key=lambda c: (c[0], c[1], c[2]))
        blocks.append([mrna_row] + [c[3] for c in children])

        prot, cds_seq, mrna_seq = m.protein, m.cds, m.mrna
        if m.coding:
            prot_fa.append(f">{mid} gene={spec.gene_id} length={len(prot)}aa\n{wrap(prot)}")
            cds_fa.append(f">{mid} gene={spec.gene_id} length={len(cds_seq)}nt\n{wrap(cds_seq)}")
        mrna_fa.append(f">{mid} gene={spec.gene_id} type={ent.type} length={len(mrna_seq)}nt\n{wrap(mrna_seq)}")

        cds_by_exon = {n: (s, e, ph) for s, e, ph, n in m.cds_parts}
        rows = []
        for n, ((s, e), (a, b)) in enumerate(zip(m.exons, ex_chrom), 1):
            intr = m.introns[n - 1] if n <= len(m.introns) else None
            cp = cds_by_exon.get(n)
            rows.append({
                "n": n, "seq": [s, e], "chrom": [a, b], "len": e - s + 1,
                "cds": [cp[0], cp[1]] if cp else None, "phase": cp[2] if cp else None,
                "intron": None if intr is None else {
                    "len": intr.length, "borders": f"{intr.donor}...{intr.acceptor}", "kind": intr.kind,
                    "hint": intr.hint},
            })
        tables.append({
            "mrna_id": mid, "mrna_name": ent.mrna_name, "name": m.name, "file": m.file_name,
            "coding": m.coding, "tx_type": ent.type,
            "length": m.length, "flipped": m.flipped,
            "placement": {"start": pl.start, "end": pl.end, "how": pl.how},
            "span": [t_start, t_end], "exons": rows,
            "cds_seq_range": [m.cds_start, m.cds_end], "cds_len": len(cds_seq),
            "start_codon": cds_seq[:3], "stop_codon": cds_seq[-3:],
            "protein": prot, "protein_len": len(prot), "mrna_len": len(mrna_seq),
            "issues": list(m.issues),
        })

    g_start, g_end = min(a for a, _ in spans), max(b for _, b in spans)
    alias = ",".join(escape_value(x.strip()) for x in re.split(r"[,;]", spec.alias) if x.strip())
    gene_type = spec.gene_type or "gene"
    if gene_type != "gene" and any(e.model.coding for e in entries):
        issues.info(f"The gene line says 'gene' rather than '{gene_type}', because at least one of its "
                    "transcripts codes for a protein.")
        gene_type = "gene"
    gene_row = [spec.seqid, spec.source, gene_type, g_start, g_end, ".", st, ".",
                _attrs([("ID", spec.gene_id), ("Name", escape_value(spec.gene_name or spec.gene_id)),
                        ("Alias", alias), ("Note", escape_value(spec.note))])]
    all_rows = [gene_row] + [r for b in blocks for r in b]
    text = "##gff-version 3\n" + "\n".join("\t".join(str(c) for c in r) for r in all_rows) + "\n"

    structures = {}
    for t in tables:
        key = tuple(tuple(x["chrom"]) for x in t["exons"]), tuple(x["cds"] and tuple(x["cds"]) for x in t["exons"])
        if key in structures:
            issues.warning(f"Transcripts {structures[key]} and {t['mrna_id']} have identical exons and CDS.")
        structures[key] = t["mrna_id"]
    for i, (a1, b1) in enumerate(spans):
        if len(spans) > 1 and all(b1 < a2 or b2 < a1 for j, (a2, b2) in enumerate(spans) if j != i):
            issues.warning(f"Transcript {entries[i].mrna_id} does not overlap any other transcript of this gene.")
    check = parse_gff3("generated", text.encode("utf-8"))
    issues.extend(check.issues)
    issues.extend(validate_features(check.features))

    report = _gene_report(spec, tables, issues, (g_start, g_end), len(all_rows), gene_type)
    fasta = ["\n".join(block) + "\n" if block else "" for block in (prot_fa, cds_fa, mrna_fa)]
    return GeneResult(text, [[str(c) for c in r] for r in all_rows], fasta[0], fasta[1], fasta[2],
                      report, tables, issues, (g_start, g_end))


def _gene_report(spec, tables, issues, span, n_lines, gene_type="gene") -> str:
    out = ["Geneious -> GFF3 report", "=" * 23, "",
           f"{gene_type} {spec.gene_id} ({spec.gene_name or spec.gene_id}) on {spec.seqid}, "
           f"strand {spec.strand}, {span[0]}-{span[1]}; {len(tables)} transcript(s); "
           f"{n_lines} GFF3 feature lines.", ""]
    for t in tables:
        pl = t["placement"]
        out += [f"{t['tx_type']} {t['mrna_id']}  <-  Geneious sequence '{t['name']}' "
                f"({t['file']}, {t['length']} bp)",
                f"  Sequence placed at {spec.seqid}:{pl['start']}-{pl['end']} ({pl['how']})",
                (f"  transcript {t['mrna_len']} nt, CDS {t['cds_len']} nt "
                 f"({t['start_codon']}...{t['stop_codon']}), protein {t['protein_len']} aa"
                 if t["coding"] else
                 f"  transcript {t['mrna_len']} nt, non-coding: no CDS, no UTR features, no protein"),
                "  exon  on sequence      on chromosome           length  CDS part         phase  intron after"]
        for x in t["exons"]:
            cds = f"{x['cds'][0]}-{x['cds'][1]}" if x["cds"] else "-"
            ph = "" if x["phase"] is None else str(x["phase"])
            intr = "" if not x["intron"] else f"{x['intron']['len']} nt {x['intron']['borders']} {x['intron']['kind']}"
            out.append(f"  {x['n']:>4}  {x['seq'][0]:>6}-{x['seq'][1]:<7}  {x['chrom'][0]:>10}-{x['chrom'][1]:<11}"
                       f"  {x['len']:>6}  {cds:<15}  {ph:^5}  {intr}")
        for lvl, txt in t["issues"]:
            out.append(f"  {lvl.upper()}: {txt}")
        out += [f"  Protein: {t['protein']}" if t["coding"] else "  (non-coding: no protein)", ""]
    for lvl, txt in issues:
        out.append(f"{lvl.upper()}: {txt}")
    return "\n".join(out).rstrip() + "\n"


# --------------------------------------------------------------------------- GFF3 reading, checking, merging
@dataclass
class GffFeature:
    file: str
    line_no: int
    cols: List[str]
    attrs: List[Tuple[str, str]]
    src: int = 0            # index of the input file (names can repeat)

    seqid = property(lambda self: self.cols[0])
    type = property(lambda self: self.cols[2])
    start = property(lambda self: int(self.cols[3]))
    end = property(lambda self: int(self.cols[4]))
    strand = property(lambda self: self.cols[6])
    phase = property(lambda self: self.cols[7])

    def get(self, key: str) -> Optional[str]:
        return next((v for k, v in self.attrs if k == key), None)

    @property
    def fid(self) -> Optional[str]:
        return self.get("ID")

    @property
    def parents(self) -> List[str]:
        v = self.get("Parent")
        return [p for p in v.split(",") if p] if v else []

    @property
    def where(self) -> str:
        return f"{self.file}, line {self.line_no}"


@dataclass
class ParsedGff:
    name: str
    features: List[GffFeature]
    directives: List[str]
    regions: Dict[str, Tuple[int, int]]
    fixes: List[str]
    issues: Issues


def parse_gff3(name: str, data: bytes, src: int = 0) -> ParsedGff:
    """Read GFF3 bytes, repairing harmless formatting problems (and listing each repair)."""
    fixes, issues = [], Issues()
    if data.startswith(b"\xef\xbb\xbf"):
        data = data[3:]
        fixes.append("removed the byte-order mark (BOM) at the start of the file")
    if data.startswith((b"\xff\xfe", b"\xfe\xff")):
        text = data.decode("utf-16")
        fixes.append("converted UTF-16 text (Excel 'Unicode Text') to UTF-8")
    else:
        try:
            text = data.decode("utf-8")
        except UnicodeDecodeError:
            text = data.decode("cp1252", errors="replace")
            fixes.append("the file was not UTF-8; it was read as Windows-1252 text")
    crlf = text.count("\r\n")
    lone_cr = text.count("\r") - crlf
    if crlf:
        fixes.append(f"converted Windows line endings (CRLF) to Unix (LF) on {crlf} line(s)")
    if lone_cr:
        fixes.append(f"converted old Mac line endings (CR) to Unix (LF) on {lone_cr} line(s)")
    lines = text.replace("\r\n", "\n").replace("\r", "\n").split("\n")

    feats, directives, regions = [], [], {}
    blank = trailing = quoted = comments = 0
    for i, raw in enumerate(lines, 1):
        if not raw.strip():
            blank += 1 if raw else 0      # lines holding only tabs/spaces (e.g. empty Excel rows)
            continue
        if raw.startswith("##"):
            d = raw.strip()
            if d.startswith("##gff-version") or d == "###":
                continue
            if d.startswith("##FASTA"):
                issues.warning(f"{name}, line {i}: an embedded ##FASTA section was found; "
                               "the sequences after it were not copied.")
                break
            if d.startswith("##sequence-region"):
                parts = d.split()
                try:
                    regions[parts[1]] = (int(parts[2]), int(parts[3]))
                except (IndexError, ValueError):
                    issues.warning(f"{name}, line {i}: unreadable ##sequence-region line was dropped.")
                continue
            directives.append(d)
            continue
        if raw.startswith("#"):
            comments += 1
            continue
        line = raw.rstrip(" \t")
        if line != raw:
            trailing += 1
        cols = line.split("\t")
        for j, c in enumerate(cols):
            if len(c) >= 2 and c[0] == '"' and c[-1] == '"':
                cols[j] = c[1:-1].replace('""', '"')
                quoted += 1
        if len(cols) != 9:
            hint = (" The columns seem to be separated by spaces instead of tabs."
                    if len(cols) < 9 and len(line.split()) >= 9 else "")
            issues.error(f"{name}, line {i}: {len(cols)} column(s) instead of 9.{hint}")
            continue
        bad = False
        for j, label in ((3, "start"), (4, "end")):
            if not cols[j].isdigit() or int(cols[j]) < 1:
                issues.error(f"{name}, line {i}: the {label} '{cols[j]}' is not a positive whole number.")
                bad = True
        if not bad and int(cols[3]) > int(cols[4]):
            issues.error(f"{name}, line {i}: start {cols[3]} is larger than end {cols[4]} "
                         "(in GFF3 start is always the smaller number, on both strands).")
            bad = True
        if cols[6] not in ("+", "-", ".", "?"):
            issues.error(f"{name}, line {i}: strand must be + or -, not '{cols[6]}'.")
            bad = True
        if cols[2] == "CDS" and cols[7] not in ("0", "1", "2"):
            issues.error(f"{name}, line {i}: a CDS line needs phase 0, 1 or 2, not '{cols[7]}'.")
            bad = True
        elif cols[2] != "CDS" and cols[7] not in (".", "0", "1", "2"):
            issues.error(f"{name}, line {i}: phase must be '.', not '{cols[7]}'.")
            bad = True
        attrs = []
        for part in cols[8].split(";"):
            part = part.strip()
            if not part or cols[8] == ".":
                continue
            if "=" not in part:
                issues.error(f"{name}, line {i}: attribute '{part}' has no '=' (expected tag=value).")
                bad = True
                continue
            k, v = part.split("=", 1)
            attrs.append((k.strip(), v.strip()))
        if not bad:
            feats.append(GffFeature(name, i, cols, attrs, src))
    if blank:
        fixes.append(f"removed {blank} line(s) containing only tabs/spaces (e.g. empty Excel rows)")
    if trailing:
        fixes.append(f"removed stray spaces/tabs at the end of {trailing} line(s)")
    if quoted:
        fixes.append(f"removed quote marks around {quoted} field(s) (an Excel artefact)")
    if comments:
        fixes.append(f"left out {comments} comment line(s) starting with #")
    if not feats and not issues.has_errors:
        issues.error(f"{name}: no feature lines were found.")
    return ParsedGff(name, feats, directives, regions, fixes, issues)


def validate_features(features: List[GffFeature]) -> Issues:
    """Cross-line checks: unique IDs, Parents, spans, CDS phases, exons of each mRNA."""
    iss = Issues()
    by_id: Dict[str, List[GffFeature]] = {}
    for f in features:
        if f.fid:
            by_id.setdefault(f.fid, []).append(f)
    for fid, fs in by_id.items():
        if len(fs) < 2:
            continue
        where = "; ".join(f.where for f in fs[:4])
        if len({(f.src, f.type, f.seqid, f.get("Parent")) for f in fs}) > 1:
            iss.error(f"ID '{fid}' is used by different features ({where}). IDs must be unique in the whole file.")
            continue
        parts = sorted((f.start, f.end) for f in fs)     # a legitimate multi-line feature (e.g. one CDS)
        if any(s2 <= e1 for (_, e1), (s2, _) in zip(parts, parts[1:])):
            iss.error(f"ID '{fid}' appears on overlapping or repeated lines ({where}); "
                      "was a line or a file included twice?")
    kids: Dict[str, List[GffFeature]] = {}
    for f in features:
        for p in f.parents:
            if p not in by_id:
                iss.error(f"{f.where}: Parent '{p}' is not defined by any ID.")
                continue
            kids.setdefault(p, []).append(f)
            par = by_id[p]
            ps, pe = min(x.start for x in par), max(x.end for x in par)
            if par[0].seqid != f.seqid:
                iss.error(f"{f.where}: on '{f.seqid}' but its parent '{p}' is on '{par[0].seqid}'.")
            elif f.start < ps or f.end > pe:
                iss.warning(f"{f.where}: {f.type} {f.start}-{f.end} extends beyond its parent '{p}' ({ps}-{pe}); "
                            "a parent's start and end should cover all of its children.")
            if par[0].strand in "+-" and f.strand in "+-" and par[0].strand != f.strand:
                iss.warning(f"{f.where}: strand {f.strand} differs from its parent '{p}' ({par[0].strand}).")
    for pid, ch in kids.items():
        cds = [k for k in ch if k.type == "CDS"]
        exons = sorted((k.start, k.end) for k in ch if k.type == "exon")
        if cds and cds[0].phase in ("0", "1", "2"):
            order = sorted(cds, key=lambda k: k.start, reverse=(cds[0].strand == "-"))
            p0, done = int(order[0].phase), 0
            for k in order:
                expected = (3 - ((done - p0) % 3)) % 3
                if k.phase != str(expected):
                    iss.error(f"{k.where}: CDS phase is {k.phase} but should be {expected}, given the lengths of "
                              f"the preceding CDS parts of '{pid}'.")
                    break
                done += k.end - k.start + 1
        for (s1, e1), (s2, e2) in zip(exons, exons[1:]):
            if s2 <= e1:
                iss.error(f"Exons {s1}-{e1} and {s2}-{e2} of '{pid}' overlap.")
        if exons:
            for k in cds:
                if not any(s <= k.start and k.end <= e for s, e in exons):
                    iss.warning(f"{k.where}: CDS {k.start}-{k.end} of '{pid}' is not inside one of its exons.")
    return iss


def _natural(s: str):
    return [int(t) if t.isdigit() else t.lower() for t in re.split(r"(\d+)", s)]


@dataclass
class MergeResult:
    text: Optional[str]
    issues: Issues
    files: List[dict]
    n_genes: int
    n_features: int


def merge_gff3(inputs: List[Tuple[str, bytes]], sort: bool = True) -> MergeResult:
    """Merge GFF3 files into one: one header, LF line endings, genes sorted by position."""
    parsed = [parse_gff3(n, d, i) for i, (n, d) in enumerate(inputs)]
    issues = Issues()
    feats: List[GffFeature] = []
    regions: Dict[str, Tuple[int, int]] = {}
    directives: List[str] = []
    files = []
    for p in parsed:
        issues.extend(p.issues)
        feats += p.features
        for sid, (a, b) in p.regions.items():
            if sid in regions and regions[sid] != (a, b):
                issues.warning(f"##sequence-region for '{sid}' differs between files; the combined range was kept.")
                a, b = min(a, regions[sid][0]), max(b, regions[sid][1])
            regions[sid] = (a, b)
        for d in p.directives:
            if d not in directives:
                directives.append(d)
        files.append({"name": p.name, "features": len(p.features),
                      "genes": sum(1 for f in p.features if f.type in GENE_LINE_TYPES),
                      "mrnas": sum(1 for f in p.features if f.type in TRANSCRIPT_TYPES),
                      "fixes": p.fixes})
    issues.extend(validate_features(feats))

    by_id: Dict[str, GffFeature] = {}
    for f in feats:
        if f.fid and f.fid not in by_id:
            by_id[f.fid] = f

    def root(f: GffFeature) -> GffFeature:
        seen = set()
        while f.parents and f.parents[0] in by_id and f.parents[0] not in seen:
            seen.add(f.parents[0])
            f = by_id[f.parents[0]]
        return f

    blocks: Dict[int, List[GffFeature]] = {}
    roots: Dict[int, GffFeature] = {}
    for f in feats:
        r = root(f)
        blocks.setdefault(id(r), []).append(f)
        roots[id(r)] = r
    order = list(blocks)
    if sort:
        order.sort(key=lambda k: (_natural(roots[k].seqid), roots[k].start, roots[k].end, roots[k].fid or ""))
    head = ["##gff-version 3"] + [f"##sequence-region {s} {a} {b}" for s, (a, b) in
                                  sorted(regions.items(), key=lambda x: _natural(x[0]))] + directives
    body = ["\t".join(f.cols) for k in order for f in blocks[k]]
    n_genes = sum(1 for k in order if roots[k].type in GENE_LINE_TYPES)
    text = None if issues.has_errors else "\n".join(head + body) + "\n"
    return MergeResult(text, issues, files, n_genes, len(feats))
