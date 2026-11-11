# Geneious to GFF3, version 1.1

A small program that runs in your web browser, on your own computer. It turns Geneious
sequences annotated with **exon** and **CDS** into GFF3 gene models with chromosome
coordinates. It can also merge several GFF3 files into one file for the genome curators.

- Protein-coding and **non-coding** transcripts alike: a sequence with exons but no CDS annotation
  becomes a lncRNA or another non-coding transcript.
- No exon coordinates to type: they are read from the Geneious annotations.
- Only two numbers to type per gene: the chromosome START (and END, if you want a check) of the sequence.
- Every file it writes has Unix (LF) line endings, so no conversion in Notepad++, VS Code or BBEdit.
- Nothing to install except Python, and nothing is sent over the internet. Flask, the web component,
  travels inside `vendor.zip` and is read straight out of it.

## Starting the tool

You need **Python 3.9 or newer** (free from <https://www.python.org/downloads/>). On Windows, tick
**"Add python.exe to PATH"** during installation. Nothing else has to be installed, and no internet is needed.

1. Unzip the download. On Windows, right-click the zip and choose **Extract All**: it makes a folder
   with the files directly inside. Keep all of them together, `vendor.zip` included, and don't run the
   tool from inside the zip.
2. Double-click the launcher:
   - **Windows:** `start_windows.bat`
   - **Mac:** `start_mac.command`. The first time, macOS may refuse to open it. Right-click it, choose
     **Open**, then **Open** again.
3. Your browser opens the tool. If it doesn't, type the address shown in the black window
   (normally `http://127.0.0.1:8765/`) into the browser.
4. Keep the black window open while you work. Close it to stop the tool.

The folder can sit deep inside other folders, even on Windows with its 260-character path limit: the
bundled packages are kept inside `vendor.zip`, one file, so Windows cannot leave some of them behind.

To check the tool without the browser, run `python selftest.py` in this folder. It rebuilds the
MtSYT1 example and compares it with the checked GFF3 in `example/`.

From the command line, `python app.py` starts the tool too.

## What the Geneious file must contain

- **One nucleotide sequence per transcript**: the genomic sequence (exons and introns) from the first
  to the last base of that transcript.
- **Annotations of type `exon`**, one per exon, including the UTR parts. Names (Ex1, Ex2, ...) and their
  order do not matter; only the annotation *type* counts.
- **For a coding transcript, one annotation of type `CDS`**, from the A of the start codon to the last base
  of the stop codon. It may span the introns; a joined (multi-part) CDS also works.
- **For a non-coding transcript, no CDS annotation at all**: annotate the exons only. See below.
- Anything else in the file (other annotation types, other documents) is ignored and listed.

The gene should read left to right. A minus-strand gene kept in chromosome orientation also works: make its
annotation arrows point right to left, and the tool reverse-complements the sequence itself.

Several splice variants can sit in one Geneious file or in separate files.

## Building one gene (one job = one gene)

1. **Upload** the Geneious file(s) of the gene.
2. **Describe the gene.** The page first shows what was read from the file: a drawing of each transcript,
   the exon/intron table, and any problems. Then fill in:
   - chromosome name (for example `MtrunA17Chr4`) and the strand of the gene;
   - gene ID, and optionally name, alias and note;
   - for a transcript with no CDS, what kind of non-coding transcript it is (see below);
   - for the first transcript, the **chromosome START** of the Geneious sequence. This is always the
     **smaller** coordinate, also for minus-strand genes. END is optional; if given, the tool checks that
     END − START + 1 equals the sequence length.
   - For further splice variants, START can stay empty. The tool finds where the variant's sequence
     matches the first transcript and places it there, needing an identical overlap of at least 200 nt.
     If you do type a START for a variant, it is checked the same way.
   - More variants can be added on the same page with **Add files**. What you typed is kept.
3. **Build GFF3.** The result page shows the model on the chromosome, one table per transcript and the
   GFF3 lines. Download the GFF3, plus FASTA files of the proteins, CDSs and mRNAs, and a text report.
   You can go back, change anything and build again.

## What is checked

Problems that **stop** the build, because the GFF3 would be wrong:
overlapping or touching exons, a CDS start or end outside the exons, a spliced CDS whose length is not a
multiple of 3, an in-frame stop codon inside the CDS, a missing exon annotation, annotation arrows
pointing both ways, a START/END that doesn't fit the sequence length, and a variant whose typed
coordinates contradict its sequence. A missing CDS annotation is **not** an error: it means the transcript
is non-coding.

Points to **check**, which don't stop the build:
- a CDS not starting with ATG;
- a CDS not ending with a stop codon. If the next codon is a stop, the tool says the CDS annotation
  probably ends 3 nt early.
- a non-canonical intron. If sliding it gives GT-AG with an identical spliced sequence, the tool gives the
  corrected exon borders; this catches border shifts made by alignment programs.
- a very short intron;
- a CDS that would only read correctly with the arrows reversed.

**Notes**: minor-class GC-AG or AT-AC introns (legitimate, like intron 9 of MtSYT1), ignored annotation
types, and sequence before exon 1 or after the last exon.

Exon, intron and splice-site checks apply to non-coding transcripts exactly as they do to coding ones; only
the CDS checks are skipped.

## Non-coding transcripts (lncRNA and others)

Leave the CDS unannotated and the tool treats the sequence as a non-coding transcript. The GFF3 then holds
the exons only: no CDS, no UTR features, no protein, and no protein or CDS FASTA file.

On the details page you choose:

- **What this transcript is** — column 3 of the transcript line, in place of `mRNA`. These are Sequence
  Ontology terms that genome browsers understand: `lnc_RNA` (the usual choice for a spliced lncRNA),
  `antisense_RNA`, `ncRNA`, `pre_miRNA`, `snoRNA`, `tRNA`, `rRNA` or plain `transcript`.
- **Feature type of the gene line** — `gene`, which is accepted everywhere, or `ncRNA_gene`, which states in
  the gene line itself that the gene is non-coding. If the curator has no preference, keep `gene`.

As a convenience the tool reports the longest ATG-to-stop reading frame it can find in the spliced
transcript. This is information only, to help you judge whether a transcript really is non-coding; nothing
about it is written to the GFF3.

A gene may mix coding and non-coding transcripts. In that case each transcript keeps its own type, and the
gene line stays `gene`.

## Merging or checking GFF3 files

Open **Merge or check GFF3 files**. Tick the genes built with the tool and/or upload other GFF3 files,
such as ones made with the Excel builder. The merged file has one `##gff-version 3` header, LF line
endings, and genes sorted by chromosome (Chr2 before Chr10) and position.

Repaired automatically, with every repair listed in the report: Windows (CRLF) or old Mac (CR) line
endings, byte-order marks, empty Excel rows, quote marks added by Excel, stray tabs at line ends, and
UTF-16 text.

Refused, so that no invalid file is written: duplicate IDs (including the same gene included twice),
Parents that don't exist, wrong CDS phases, start > end, wrong column counts (spaces instead of tabs),
overlapping exons within an mRNA.

To clean and check a **single** file, merge just that one file.

## Keeping version 1.0 as well

Version 1.1 (this one) is made to sit beside the earlier version without disturbing it:

- Unzip it into a **folder of its own**; do not overwrite the old one.
- It uses **port 8775**, the old one 8765, so both can run at the same time.
- On Windows it keeps its jobs in `Geneious_to_GFF3_v1_1\jobs`, apart from the old version's jobs.
  On a Mac the jobs always live next to the program, so separate folders are separate anyway.
- The version is shown next to the name at the top of every page, so you can tell the two windows apart.

The genes you built with the old version stay in the old version. Gene models already downloaded are
unaffected, and both versions write the same GFF3 for a coding gene: the MtSYT1 example is byte-identical.
Once you are happy with 1.1 you can simply delete the old folder (and, on Windows, the old jobs folder).

## Where your jobs are kept

Each gene you build and each merge is a job, listed on the Jobs page, where you can also delete it. The
exact folder is shown at the bottom of every page:

- **Windows:** `C:\Users\<you>\AppData\Local\Geneious_to_GFF3_v1_1\jobs`. That short location keeps
  every path under Windows' 260-character limit, however deep the tool folder is.
- **Mac:** the `jobs` folder next to the program.

## If something goes wrong

- **"Python 3.9 or newer was not found"**: install Python from python.org (Windows: tick "Add python.exe
  to PATH"), then double-click the launcher again.
- **"The tool could not load its web component"**: `vendor.zip` is not next to `app.py`. Unzip the
  download again and keep every file together. On Windows, avoid extracting into a folder that already
  carries a long path, and avoid ending up with the folder twice (`...\Geneious_to_GFF3_v1_1\Geneious_to_GFF3_v1_1\`):
  files whose full path passes 260 characters are dropped without a warning.
- **The browser shows nothing**: look in the black window for the address. If port 8765 is taken, the
  tool picks the next free one.
- **A different port or jobs folder**: set the environment variables `GFF3TOOL_PORT` or `GFF3TOOL_DATA`.
  `GFF3TOOL_NO_BROWSER=1` stops the browser from opening.

## Files

| File | Purpose |
|---|---|
| `app.py` | the web interface (Flask) |
| `gff3core.py` | reading Geneious files, building and checking models, writing and merging GFF3; usable without the web interface |
| `templates/` | the pages |
| `start_windows.bat`, `start_mac.command` | double-click launchers |
| `selftest.py`, `example/` | check with the MtSYT1 example, coding and non-coding |
| `vendor.zip` | Flask and its helper packages, imported straight from the zip (unmodified, BSD/MIT licensed) |
| `LICENSES_bundled.txt` | the licenses of everything in `vendor.zip` |
