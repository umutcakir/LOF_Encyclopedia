## Abstract:

*Medicago truncatula* is an important model plant for studies on symbiotic nitrogen fixation and many fundamental biological processes. Despite its value for the global community of plant biologists, *M. truncatula* still has no consolidated resource listing all genes targeted by loss-of-function (LOF) approaches. We have filled this gap by collecting and analyzing information on 769 genes studied using LOF methods since the model was introduced more than 30 years ago (655 studies). Because many publications, including recent ones, reported outdated, incomplete, incorrect, or no information about the identity of loci affected in LOF studies, a considerable effort was made to ensure the accurate mapping of affected loci to the latest genome assembly (v. 5.1.9; MtrunA17r5.0-ANR). The dataset lists not only detailed phenotypes of monogenic mutants, but also phenotypes observed in polygenic mutant backgrounds, even when such phenotypes are not different from the wild type. This comprehensive compendium also points out numerous inconsistencies, ambiguities, and errors in naming of genetic loci targeted by LOF approaches. Lastly, we restructured 12 gene models based on RNA-Seq data, with some models also supported by mass spectrometry proteomics. This updatable resource is publicly available as a separate track in the *M. truncatula* genome browser.

## Extended summary:

*M. truncatula* has been the main legume model for symbiotic nitrogen fixation (SNF) and arbuscular mycorrhizal (AM) symbiosis for more than 30 years, yet it had no consolidated record of genes targeted by loss-of-function (LOF) approaches. Without one, researchers duplicate effort, give conflicting names to genes, and publish outdated or missing locus identifiers. The authors built a manually curated LOF compendium mapped to the current genome assembly (v. 5.1.9; MtrunA17r5.0-ANR) and released it as a dedicated track in the *M. truncatula* genome browser.

### Approach

Current locus identifiers were resolved from primer sequences, RNAi constructs, Tnt1 flanking sequence tags, BLAST searches and correspondence with dozens of authors; language models helped gather literature but failed at the core interpretation. Phenotypes are classed as monogenic altered (M), polygenic altered (P) (seen in lines with several disrupted genes, without assigning causation) or unaltered (U), and every potentially affected locus is listed.

### Main findings

#### 1. The genome browser under-reported LOF evidence

Before the new track, only 266 of the 769 genes (35%) carried correct LOF references. A further 228 (30%) had no functional annotation at all and 256 (33%) had references that did not point to the LOF work. More than half of the genes characterized since 2023 were unflagged, but missed genes date back to 2001, with no bias by method, phenotype or protein class.

#### 2. Four types of naming errors

- Several names for one gene: the locus first studied as MtAOC1 (2005) was renamed MtAOC2 in 2023; similar cases are MtLOX9/MtLOX24/MtLOX4 and MtMAX2/MtDW1.
- One name for several genes: MtCHS was used for three different loci; MtMYB1, MtbZIP60 and MtSymSCL1 each for two. Provisional curator suffixes (e.g. MtCHSa–c) separate them.
- Converging acronyms: MtCRA1/2 denote both COMPACT ROOT ARCHITECTURE and CHITINASE-RELATED AGGLUTININ genes.
- Misassigned names in browser v. 5.1.9: at least 11 published names are attached to the wrong loci (e.g. MtMYC2).

#### 3. Gene-model corrections

Six v. 4 genes are split in v. 5.1.9 into adjacent loci only 1-5 nt apart: MtGSHS1, MtLICK1, MtROP7, MtSYP132, MtSYT1 and the uncharacterized MtbZIP. Reads from 50 public RNA-Seq runs support unified transcripts, mostly formed by alternative splicing, and junction-spanning MS peptides confirm translation of unified MtLICK1 (10-day nodules) and MtSYT1 (flowers, seeds). The original qRT-PCR primers for MtROP7 and MtSYT1 also span the junction. MtbZIP may be tri-cistronic through alternative translation start sites. Because independent transcription of the short models is not excluded, the unified models were added to the "Expert annotation" track alongside the originals. A separate problem is illustrated by MtMYB14: it has two identical loci in A17 but a single copy in R108, so the copy hit by Tnt1 cannot be assigned.

#### 4. Methods, processes and protein classes

Tnt1 insertional mutagenesis dominates the field. RNAi remains popular despite its off-target ambiguity: of 191 genes first studied in 2023-2026, 96 used Tnt1, 62 RNAi and 61 CRISPR/Cas9. The most studied processes are SNF (444 genes), root development (193) and AM symbiosis (153); 127 genes were studied for seed traits. Enzymes (263) and transcription factors (TFs; 238) account for more than half of the genes; substrates were identified for 65 of the 81 transporters.

#### 5. Phenotype proportions as a "cheat sheet"

Overall, 63% of genes show an M-phenotype for at least one process. Only 54% of genes tested for SNF do so (59% for non-SNF), which may reflect redundancy, although RNAi, common in hairy-root SNF work, is strongly associated with P-phenotypes. Of 262 genes tested in both contexts, 41% were M in both and only 9% were altered for SNF alone, suggesting that purely SNF-specific genes are uncommon. Protein classes differ significantly: kinases have more M- and fewer U-phenotypes than oxidoreductases; MYB-HB-like TFs reach 77% M versus 50% for AP2-EREBP; transporters (67% M) outperform TFs (55%); channels have the most U-phenotypes (36%). Some extreme ratios (SBP, C2C2-Dof) stem from single studies.

### Take-home message

The LOF track turns 30 years of scattered *M. truncatula* mutant literature into a curated, citable, browser-integrated reference covering 769 loci (current, it may be slightly higher in the recent update). It corrects identifiers and names, validates six unified gene models, and provides family-level probabilities of obtaining a scorable phenotype.
