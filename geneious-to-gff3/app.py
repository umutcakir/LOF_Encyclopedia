"""
app.py - local web interface of the Geneious -> GFF3 tool.

Start it with   python app.py   (or the start_windows.bat / start_mac.command launchers)
and your browser opens the tool. Everything runs on your own computer.
Flask and its helper packages are included in the 'vendor' folder: nothing needs installing.
"""
import json
import os
import re
import secrets
import shutil
import socket
import sys
import threading
import webbrowser
from datetime import datetime
from pathlib import Path

VERSION = "1.1"            # 1.0 coding genes only; 1.1 adds non-coding transcripts

APP_DIR = Path(__file__).resolve().parent
# Flask and its helper packages travel in one zip file: Python imports them straight out of it.
# A single file cannot be half-extracted, which matters on Windows, where a deep folder can push
# individual files past the 260-character path limit. An unpacked 'vendor' folder also still works.
for _bundle in (APP_DIR / "vendor.zip", APP_DIR / "vendor"):
    if _bundle.exists():
        sys.path.insert(0, str(_bundle))

try:
    from flask import Flask, abort, flash, redirect, render_template, request, send_file, url_for
    from markupsafe import Markup, escape
except ImportError as exc:
    sys.exit(f"\n  The tool could not load its web component ({exc}).\n\n"
             f"  The file 'vendor.zip' must sit next to app.py, in this folder:\n"
             f"    {APP_DIR}\n"
             f"  {'It is there, so it may have been damaged: unzip the download again.' if (APP_DIR / 'vendor.zip').exists() else 'It is missing: unzip the download again, keeping every file together.'}\n"
             f"  Python 3.9 or newer is needed; this is Python {sys.version.split()[0]}.\n")

import gff3core as core


def default_data_dir():
    """Jobs folder. Windows limits file paths to 260 characters, and this tool may sit deep inside
    other folders, so on Windows jobs are kept in a short place in the user's profile."""
    if os.name == "nt":
        base = os.environ.get("LOCALAPPDATA") or str(Path.home() / "AppData" / "Local")
        return Path(base) / "Geneious_to_GFF3" / "jobs"
    return APP_DIR / "jobs"


DATA_DIR = Path(os.environ.get("GFF3TOOL_DATA") or default_data_dir())
DATA_DIR.mkdir(parents=True, exist_ok=True)
JOB_ID = re.compile(r"^[a-f0-9]{12}$")
GENE_FIELDS = ("seqid", "strand", "gene_id", "gene_name", "alias", "note", "source", "gene_type")
TX_FIELDS = ("mrna_id", "mrna_name", "start", "end", "tx_type")

app = Flask(__name__)
app.config["MAX_CONTENT_LENGTH"] = 200 * 1024 * 1024
app.secret_key = secrets.token_hex(16)


@app.template_filter("thousands")
def thousands(n):
    try:
        return f"{int(n):,}"
    except (TypeError, ValueError):
        return n


def gene_glyph(tracks, lo, hi, reverse=False, ruler=True, width=960):
    """Inline SVG drawn like a genome-browser track: thin line = intron, low box = UTR,
    tall box = coding part, chevrons = direction of transcription.
    'width' is the drawing's own width; pick it close to the space it is shown in."""
    lab = (128 if width >= 700 else 70) if any(t[0] for t in tracks) else 0
    row, top = 34, 4
    height = top + row * len(tracks) + (20 if ruler else 2)
    span = max(1, hi - lo + 1)

    def sx(p):
        return round(lab + (p - lo) / span * (width - lab - 2), 2)

    out = [f'<svg class="glyph" viewBox="0 0 {width} {height}" role="img" '
           f'aria-label="Exon-intron structure">']
    for i, (label, exons, cds) in enumerate(tracks):
        y = top + row * i + row / 2
        if label:
            out.append(f'<text class="gl-label" x="0" y="{y + 4}">{escape(label)}</text>')
        if not exons:
            continue
        out.append(f'<line class="gl-line" x1="{sx(exons[0][0])}" x2="{sx(exons[-1][1] + 1)}" y1="{y}" y2="{y}"/>')
        for (_, e1), (s2, _) in zip(exons, exons[1:]):
            a, b = sx(e1 + 1), sx(s2)
            if b - a > 16:
                m, d = (a + b) / 2, (-3 if reverse else 3)
                out.append(f'<path class="gl-chev" d="M{m - d},{y - 4} L{m + d},{y} L{m - d},{y + 4}"/>')
        for n, (s, e) in enumerate(exons, 1):
            num = len(exons) - n + 1 if reverse else n
            out.append(f'<rect class="gl-utr" x="{sx(s)}" y="{y - 5}" width="{max(1.5, sx(e + 1) - sx(s))}" '
                       f'height="10"><title>Exon {num}: {s:,}-{e:,} ({e - s + 1:,} nt)</title></rect>')
        for s, e in cds:
            out.append(f'<rect class="gl-cds" x="{sx(s)}" y="{y - 9}" width="{max(1.5, sx(e + 1) - sx(s))}" '
                       f'height="18"><title>Coding part: {s:,}-{e:,}</title></rect>')
    if ruler:
        yb = height - 5
        out.append(f'<text class="gl-ruler" x="{lab}" y="{yb}">{lo:,}</text>')
        out.append(f'<text class="gl-ruler" x="{width - 2}" y="{yb}" text-anchor="end">{hi:,}</text>')
    out.append("</svg>")
    return Markup("".join(out))


def tracks_from_rows(rows):
    """One drawing track per transcript, whatever its feature type (mRNA, lnc_RNA, ...)."""
    parsed = [(r, dict(kv.split("=", 1) for kv in r[8].split(";") if "=" in kv)) for r in rows]
    genes = {a["ID"] for r, a in parsed if "ID" in a and "Parent" not in a}
    tx, order = {}, []
    for r, a in parsed:
        if a.get("Parent") in genes and "ID" in a:
            tx[a["ID"]] = ([], [])
            order.append(a["ID"])
        elif r[2] in ("exon", "CDS") and a.get("Parent") in tx:
            tx[a["Parent"]][0 if r[2] == "exon" else 1].append((int(r[3]), int(r[4])))
    return [(m, sorted(tx[m][0]), sorted(tx[m][1])) for m in order]


@app.context_processor
def template_helpers():
    return {"data_dir": DATA_DIR, "gene_glyph": gene_glyph, "version": VERSION}


# --------------------------------------------------------------------------- job storage
def now():
    return datetime.now().strftime("%Y-%m-%d %H:%M")


def job_path(job_id):
    if not JOB_ID.match(job_id or "") or not (DATA_DIR / job_id / "job.json").is_file():
        abort(404)
    return DATA_DIR / job_id


def load_job(job_id):
    with open(job_path(job_id) / "job.json", encoding="utf-8") as fh:
        return json.load(fh)


def save_job(job):
    job["updated"] = now()
    target = DATA_DIR / job["id"] / "job.json"
    tmp = target.with_suffix(".tmp")
    with open(tmp, "w", encoding="utf-8", newline="\n") as fh:
        json.dump(job, fh, indent=1)
    tmp.replace(target)


def new_job(kind):
    job_id = secrets.token_hex(6)
    d = DATA_DIR / job_id
    (d / "input").mkdir(parents=True)
    (d / "output").mkdir()
    job = {"id": job_id, "kind": kind, "created": now(), "updated": now(), "status": "draft",
           "files": [], "gene": {"source": "manual_curation"}, "tx": {}, "outputs": {}}
    save_job(job)
    return job


def list_jobs():
    jobs = []
    for d in DATA_DIR.iterdir():
        if d.is_dir() and JOB_ID.match(d.name) and (d / "job.json").is_file():
            try:
                with open(d / "job.json", encoding="utf-8") as fh:
                    jobs.append(json.load(fh))
            except (OSError, ValueError):
                continue
    return sorted(jobs, key=lambda j: j.get("updated", ""), reverse=True)


def write_output(job_id, name, text):
    # binary write: keeps Unix (LF) line endings even on Windows
    with open(DATA_DIR / job_id / "output" / name, "wb") as fh:
        fh.write(text.encode("utf-8"))


def store_uploads(job, uploads):
    added = []
    for fs in uploads:
        if not fs or not fs.filename:
            continue
        suffix = re.sub(r"[^A-Za-z0-9.]", "", Path(fs.filename).suffix)[:12] or ".geneious"
        stored = f"{len(job['files']) + 1:02d}{suffix}"      # short on purpose (path limits)
        fs.save(DATA_DIR / job["id"] / "input" / stored)
        job["files"].append({"stored": stored, "original": fs.filename})
        added.append(fs.filename)
    return added


def safe_base(text, limit=40):
    return (re.sub(r"[^A-Za-z0-9._-]+", "_", text).strip("._") or "gene")[:limit]


# --------------------------------------------------------------------------- building a gene
def candidates_of(job):
    """Every annotated sequence in the job's files, with its checked transcript model."""
    cands, messages = [], core.Issues()
    for fi, f in enumerate(job["files"]):
        data = (DATA_DIR / job["id"] / "input" / f["stored"]).read_bytes()
        records, iss = core.read_geneious(data, f["original"])
        messages.extend(iss)
        for ri, rec in enumerate(r for r in records if r.is_annotated):
            model = core.build_model(rec)
            cands.append({"key": f"{fi}-{ri}", "file": f["original"], "model": model,
                          "n_exons": len(model.exons), "ignored": rec.ignored_types(),
                          "coding": model.coding,
                          "glyph_cds": [(s, e) for s, e, _, _ in model.cds_parts],
                          "cds_by_exon": {n: (s, e, ph) for s, e, ph, n in model.cds_parts}})
    return cands, messages


def read_form(job, cands):
    for k in GENE_FIELDS:
        job["gene"][k] = request.form.get(k, "").strip()
    for c in cands:
        t = job["tx"].setdefault(c["key"], {})
        t["include"] = request.form.get(f"tx-{c['key']}-include") == "on"
        for k in TX_FIELDS:
            t[k] = request.form.get(f"tx-{c['key']}-{k}", "").strip()


def parse_coord(text):
    cleaned = re.sub(r"[\s,._']", "", text or "")
    return int(cleaned) if cleaned.isdigit() else None


def validate_and_build(job, cands):
    """Returns (GeneResult or None, field_errors, general_errors)."""
    errors, general = {}, []
    g = job["gene"]
    g["source"] = g.get("source") or "manual_curation"
    if not g.get("seqid"):
        errors["seqid"] = "Enter the chromosome name, for example MtrunA17Chr4."
    elif not core.SEQID_OK.match(g["seqid"]):
        errors["seqid"] = "Use only letters, digits and . : _ - (no spaces)."
    if g.get("strand") not in ("+", "-"):
        errors["strand"] = "Choose the strand: + or -."
    if not g.get("gene_id"):
        errors["gene_id"] = "Enter a gene ID."
    elif not core.ID_OK.match(g["gene_id"]):
        errors["gene_id"] = "Use only letters, digits and . _ - : (no spaces, and none of ; = , &)."
    if "\t" in g["source"] or ";" in g["source"]:
        errors["source"] = "The source must not contain tabs or semicolons."
    if g.get("gene_type") not in [t for t, _ in core.GENE_TYPES]:
        g["gene_type"] = "gene"

    included = [c for c in cands if job["tx"].get(c["key"], {}).get("include", True)]
    if not cands:
        general.append("No sequence with exon/CDS annotations was found in the uploaded file(s).")
    elif not included:
        general.append("Tick at least one transcript to include.")
    for c in included:
        if c["model"].issues.has_errors:
            general.append(f"Transcript '{c['model'].name}' ({c['file']}) has problems that must be fixed in "
                           "Geneious first (listed in red under it). Fix them, upload the corrected file, "
                           "or untick this transcript.")

    used = set()
    for n, c in enumerate(included, 1):
        t = job["tx"].setdefault(c["key"], {})
        mid = t.get("mrna_id") or (f"{g['gene_id']}.{n}" if g.get("gene_id") else "")
        field = f"tx-{c['key']}-mrna_id"
        if mid and not core.ID_OK.match(mid):
            errors[field] = "Use only letters, digits and . _ - : (no spaces)."
        elif mid and mid == g.get("gene_id"):
            errors[field] = "The mRNA ID must differ from the gene ID."
        elif mid in used:
            errors[field] = "This mRNA ID is used twice; each transcript needs its own."
        used.add(mid)
        c["mrna_id"] = mid

        tfield = f"tx-{c['key']}-tx_type"
        if c["model"].coding:
            c["tx_type"] = core.CODING_TX_TYPE      # it has a CDS, so it is an mRNA
        else:
            chosen = (t.get("tx_type") or "").strip() or core.NONCODING_TX_TYPES[0][0]
            if not core.SO_TYPE_OK.match(chosen):
                errors[tfield] = ("Use a Sequence Ontology term: letters, digits and underscores only, "
                                  "for example lnc_RNA.")
            elif chosen == core.CODING_TX_TYPE:
                errors[tfield] = ("This sequence has no CDS annotation, so it cannot be written as an mRNA. "
                                  "Choose a non-coding type, or annotate the CDS in Geneious and upload "
                                  "the file again.")
            c["tx_type"] = chosen

    placed = {}
    for c in included:
        key, t, L = c["key"], job["tx"][c["key"]], c["model"].length
        fs, fe = f"tx-{key}-start", f"tx-{key}-end"
        s, e = parse_coord(t.get("start")), parse_coord(t.get("end"))
        if t.get("start") and s is None:
            errors[fs] = "Enter a whole number (digits only)."
        if t.get("end") and e is None:
            errors[fe] = "Enter a whole number (digits only)."
        if s is None and e is not None:
            s = e - L + 1
        if s is not None and s < 1:
            errors[fs if t.get("start") else fe] = "The region would start before position 1 of the chromosome."
        elif s is not None and e is not None and e - s + 1 != L:
            errors[fe] = (f"END - START + 1 = {e - s + 1:,}, but this sequence is {L:,} bp long. "
                          f"If START is right, END should be {s + L - 1:,}.")
        elif s is not None:
            how = "coordinates entered" if t.get("start") else "START calculated from the END entered"
            placed[key] = core.Placement(s, s + L - 1, how)

    if included and not placed and not any(k.endswith(("-start", "-end")) for k in errors):
        errors[f"tx-{included[0]['key']}-start"] = "Enter the chromosome START of this sequence."
    strand = g.get("strand")
    ref = next((c for c in included if c["key"] in placed), None)
    if ref is not None and strand in ("+", "-"):
        ref_pl = placed[ref["key"]]
        for c in included:
            if c is ref or c["model"].issues.has_errors:
                continue
            fs = f"tx-{c['key']}-start"
            offs = core.overlap_offsets(ref["model"].seq, c["model"].seq)
            if c["key"] in placed:
                d = core.offset_between(ref_pl, placed[c["key"]], strand)
                if core.sequences_agree(ref["model"].seq, c["model"].seq, d) is False and fs not in errors:
                    msg = (f"With these coordinates this sequence does not match transcript "
                           f"{ref['mrna_id'] or ref['model'].name} where the two overlap.")
                    if len(offs) == 1:
                        alt = core.placement_from_offset(offs[0], ref_pl, c["model"].length, strand, "")
                        msg += f" Matching the sequences suggests START = {alt.start:,}."
                    errors[fs] = msg
            elif len(offs) == 1:
                placed[c["key"]] = core.placement_from_offset(
                    offs[0], ref_pl, c["model"].length, strand,
                    f"placed automatically: its sequence matches {ref['mrna_id']} with an offset of {offs[0]:+,} nt")
            else:
                errors[fs] = ("Could not place this variant automatically "
                              + ("(its sequence matches the first transcript at several positions). "
                                 if offs else "(no identical overlap of at least 200 nt with the first transcript). ")
                              + "Enter its START.")
    if errors or general:
        return None, errors, general
    spec = core.GeneSpec(g["seqid"], strand, g["gene_id"], g.get("gene_name", ""), g.get("alias", ""),
                         g.get("note", ""), g["source"], g.get("gene_type") or "gene")
    entries = [core.TranscriptEntry(c["model"], placed[c["key"]], c["mrna_id"],
                                    job["tx"][c["key"]].get("mrna_name", ""), c["tx_type"])
               for c in included]
    return core.build_gene(spec, entries), {}, []


def save_build(job, res):
    out = DATA_DIR / job["id"] / "output"
    for old in out.iterdir():
        old.unlink()
    base = safe_base(job["gene"]["gene_id"])
    coding = any(t["coding"] for t in res.transcripts)
    files = {"gff3": (f"{base}.gff3", res.gff3),
             "mrna": (f"{base}_mRNA.fa" if coding else f"{base}_transcripts.fa", res.mrna_fa),
             "report": (f"{base}_report.txt", res.report)}
    if res.proteins_fa:       # left out when no transcript of this gene codes for a protein
        files["proteins"] = (f"{base}_proteins.fa", res.proteins_fa)
        files["cds"] = (f"{base}_CDS.fa", res.cds_fa)
    for name, text in files.values():
        write_output(job["id"], name, text)
    job["outputs"] = {k: name for k, (name, _) in files.items()}
    job["status"] = "built"
    job["result"] = {"transcripts": res.transcripts, "issues": list(res.issues), "rows": res.rows,
                     "span": list(res.span), "coding": coding}


# --------------------------------------------------------------------------- routes
@app.route("/")
def index():
    return render_template("index.html", jobs=list_jobs(), data_dir=DATA_DIR)


@app.route("/build", methods=["GET", "POST"])
def build_start():
    if request.method == "POST":
        uploads = [f for f in request.files.getlist("files") if f and f.filename]
        if not uploads:
            flash("Choose at least one .geneious file.", "error")
            return redirect(url_for("build_start"))
        job = new_job("build")
        store_uploads(job, uploads)
        save_job(job)
        return redirect(url_for("build_details", job_id=job["id"]))
    return render_template("build_upload.html")


@app.route("/build/<job_id>", methods=["GET", "POST"])
def build_details(job_id):
    job = load_job(job_id)
    if job["kind"] != "build":
        abort(404)
    cands, messages = candidates_of(job)
    errors, general = {}, []
    if request.method == "POST":
        read_form(job, cands)
        added = store_uploads(job, request.files.getlist("more_files"))
        if request.form.get("action") == "add":
            save_job(job)
            flash(f"Added: {', '.join(added)}" if added else "Choose a file to add first.",
                  "info" if added else "error")
            return redirect(url_for("build_details", job_id=job_id) + "#transcripts")
        if added:
            cands, messages = candidates_of(job)
        res, errors, general = validate_and_build(job, cands)
        if res is not None:
            save_build(job, res)
            save_job(job)
            return redirect(url_for("build_result", job_id=job_id))
        job["status"] = "draft"
        save_job(job)
    for n, c in enumerate(cands, 1):
        c["form"] = job["tx"].get(c["key"], {"include": True})
        c["index"] = n
    return render_template("build_details.html", job=job, cands=cands, messages=messages,
                           errors=errors, general=general, gene_types=core.GENE_TYPES,
                           noncoding_types=core.NONCODING_TX_TYPES)


@app.route("/build/<job_id>/result")
def build_result(job_id):
    job = load_job(job_id)
    if job.get("status") != "built":
        return redirect(url_for("build_details", job_id=job_id))
    res = job["result"]
    n_warn = sum(1 for lvl, _ in res["issues"] if lvl == "warning") + sum(
        1 for t in res["transcripts"] for lvl, _ in t["issues"] if lvl == "warning")
    return render_template("build_result.html", job=job, res=res, tracks=tracks_from_rows(res["rows"]),
                           n_warn=n_warn)


@app.route("/job/<job_id>/download/<kind>")
def download(job_id, kind):
    job = load_job(job_id)
    name = job.get("outputs", {}).get(kind)
    if not name:
        abort(404)
    return send_file(job_path(job_id) / "output" / name, as_attachment=True, download_name=name,
                     mimetype="text/plain")


@app.route("/job/<job_id>/delete", methods=["POST"])
def delete_job(job_id):
    shutil.rmtree(job_path(job_id), ignore_errors=True)
    flash("Job deleted.", "info")
    return redirect(url_for("index"))


@app.route("/merge", methods=["GET", "POST"])
def merge():
    built = [j for j in list_jobs() if j.get("kind") == "build" and j.get("status") == "built"]
    if request.method == "POST":
        inputs = []
        for jid in request.form.getlist("jobs"):
            j = load_job(jid)
            name = j.get("outputs", {}).get("gff3")
            if name:
                inputs.append((name, (job_path(jid) / "output" / name).read_bytes()))
        for fs in request.files.getlist("files"):
            if fs and fs.filename:
                inputs.append((fs.filename, fs.read()))
        if not inputs:
            flash("Tick at least one gene or upload at least one GFF3 file.", "error")
            return redirect(url_for("merge"))
        out_name = safe_base(request.form.get("out_name") or "new_gene_models", 60)
        out_name = out_name if out_name.lower().endswith(".gff3") else out_name + ".gff3"
        res = core.merge_gff3(inputs, sort=request.form.get("sort") == "on")
        mj = new_job("merge")
        mj.update({"inputs": [n for n, _ in inputs], "status": "merged" if res.text else "failed",
                   "report": {"files": res.files, "issues": list(res.issues), "n_genes": res.n_genes,
                              "n_features": res.n_features}})
        if res.text:
            write_output(mj["id"], out_name, res.text)
            mj["outputs"] = {"gff3": out_name}
            mj["preview"] = res.text.splitlines()[:80]
        save_job(mj)
        return redirect(url_for("merge_result", job_id=mj["id"]))
    return render_template("merge.html", built=built)


@app.route("/merge/<job_id>")
def merge_result(job_id):
    job = load_job(job_id)
    if job["kind"] != "merge":
        abort(404)
    return render_template("merge_result.html", job=job, rep=job["report"])


@app.errorhandler(413)
def too_large(_):
    flash("That upload is too large (limit 200 MB).", "error")
    return redirect(url_for("index"))


# --------------------------------------------------------------------------- start-up
def free_port(start):
    for port in range(start, start + 50):
        with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
            if s.connect_ex(("127.0.0.1", port)) != 0:
                return port
    return start


def main():
    import logging
    import flask.cli
    logging.getLogger("werkzeug").setLevel(logging.ERROR)      # no developer warnings or request log
    flask.cli.show_server_banner = lambda *args, **kwargs: None
    host = os.environ.get("GFF3TOOL_HOST", "127.0.0.1")
    port = free_port(int(os.environ.get("GFF3TOOL_PORT", "8765")))
    url = f"http://127.0.0.1:{port}/"
    if os.environ.get("GFF3TOOL_NO_BROWSER") != "1":
        threading.Timer(1.5, lambda: webbrowser.open(url)).start()
    print(f"\n  Geneious to GFF3 {VERSION} is running at  {url}\n"
          f"  Your browser should open by itself. If it does not, type that address into it.\n"
          f"  Keep this window open while you work, and close it when you are done.\n", flush=True)
    app.run(host=host, port=port, debug=False, load_dotenv=False)


if __name__ == "__main__":
    main()
