"""Quick check that the tool works on this computer (no browser needed):  python selftest.py

Builds the MtSYT1 example from example/MtSYT1_auto.geneious and compares the result with
example/MtSYT1_expected.gff3, a model checked by hand, with GenomeTools, and accepted by the
Medicago truncatula v5 curators. Then builds the same sequence without its CDS annotation, to
check the non-coding (lncRNA) path."""
import io
import zipfile
import xml.etree.ElementTree as ET
from pathlib import Path

import gff3core as core

HERE = Path(__file__).resolve().parent / "example"
NOTE = ("Joins v5 loci MtrunA17_Chr4g0037111 and MtrunA17_Chr4g0037121 into a unified model identical to "
        "v4 locus Medtr4g073400. Evidence based on RNA-Seq and MS proteomic data.")
RAW = (HERE / "MtSYT1_auto.geneious").read_bytes()

# ---- coding: the model the curators accepted
records, _ = core.read_geneious(RAW, "MtSYT1_auto.geneious")
model = core.build_model([r for r in records if r.is_annotated][0])
spec = core.GeneSpec("MtrunA17Chr4", "+", "MtSYT1", "MtSYT1", "Medtr4g073400", NOTE)
placement = core.Placement(35768230, 35768230 + model.length - 1, "entered")
result = core.build_gene(spec, [core.TranscriptEntry(model, placement, "MtSYT1.1")])
expected = (HERE / "MtSYT1_expected.gff3").read_bytes().decode("utf-8")

# ---- non-coding: the same sequence with its CDS annotation taken out
with zipfile.ZipFile(io.BytesIO(RAW)) as z:
    entry = z.namelist()[0]
    root = ET.fromstring(z.read(entry))
for holder in root.iter():
    block = holder.find("sequenceAnnotations")
    if block is not None:
        for ann in list(block):
            if (ann.findtext("type") or "").strip().lower() == "cds":
                block.remove(ann)
buf = io.BytesIO()
with zipfile.ZipFile(buf, "w", zipfile.ZIP_DEFLATED) as z:
    z.writestr(entry, ET.tostring(root, encoding="utf-8", xml_declaration=True))
nc_records, _ = core.read_geneious(buf.getvalue(), "no_cds.geneious")
nc_model = core.build_model([r for r in nc_records if r.is_annotated][0])
nc = core.build_gene(core.GeneSpec("MtrunA17Chr4", "+", "LNC", "LNC", "", "", gene_type="ncRNA_gene"),
                     [core.TranscriptEntry(nc_model, placement, "LNC.1")])
nc_types = [line.split("\t")[2] for line in nc.gff3.splitlines() if not line.startswith("#")]

checks = {
    "12 exons and 1 CDS read by annotation type": len(model.exons) == 12 and model.cds_start == 383,
    "protein of 537 amino acids": len(model.protein) == 537,
    "GFF3 identical to the checked example": result.gff3 == expected,
    "Unix (LF) line endings only": "\r" not in result.gff3 and "\r" not in nc.gff3,
    "a sequence without a CDS is accepted as non-coding":
        not nc_model.coding and not nc_model.issues.has_errors,
    "the non-coding GFF3 holds ncRNA_gene, lnc_RNA and 12 exons":
        nc_types == ["ncRNA_gene", "lnc_RNA"] + ["exon"] * 12,
    "no protein or CDS is written for it": nc.proteins_fa == "" and nc.cds_fa == "",
}
for label, passed in checks.items():
    print(("  ok    " if passed else "  FAIL  ") + label)
print("\nSelf-test passed." if all(checks.values()) else "\nSelf-test FAILED - please report this.")
