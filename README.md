# TE-Discovery-Pipeline
Pipeline for discovering novel/underannotated TEs in the human genome

**NB: Currently Incomplete**

## Step 1 — Genome preparation

Run from inside the chromosome subdirectory (e.g. `cd /home/lukeok/RepeatFinder/T2T/chr22`).

1. **Fix FASTA header** — rename the NCBI accession to a simple `>chrN` header.
2. **Index** — `samtools faidx` + generate `.genome` file (chromosome name + size).
3. **Extract known repeats** — parse `chm13v2.0_RepeatMasker_4.1.2p1.2022Apr14.out` for the target chromosome and write a BED file of known repeat intervals.
4. **Soft-mask** — `bedtools maskfasta -soft` masks known repeats before `BuildDatabase`.
5. **Build database** — `BuildDatabase chrNdb <chrN_softmasked.fa>`

> **Soft-masking does not suppress re-discovery of known repeats.** RepeatModeler's de novo engines read lowercase as ordinary nucleotides. Every soft-masked copy remains in the signal, so known repeats — especially old, diverged relics below RepeatMasker's threshold — are freely re-modelled. Steps 3b and 9a–9b are required regardless.

---

## Step 2 — De novo repeat modelling

```bash
RepeatModeler -database chrNdb -threads 8 -LTRStruct
```

Runs four iterative rounds (RepeatScout + RECON) then a separate LTR structural pipeline (LTR_Retriever). Output: `RM_*/consensi.fa.classified` and `families.stk`.

Classification uses RepeatClassifier (BLAST against `RepeatPeps.lib` + `RepeatMasker.lib`).

---

## Step 3 — Extract Unknown families

Families classified as `#Unknown` are extracted from `consensi.fa.classified` → `unknown_consensi.fa`.

---

## Step 3b — TRF / satellite pre-filter ⚠️ run before Step 4

Satellites and simple tandem repeats escape BLAST filtering because BLAST's default DUST filter silently masks low-complexity sequence before alignment, producing zero hits even against known counterparts. RepeatClassifier has the same blind spot. Without this step, satellite and VNTR consensi pass all downstream filters and appear as novel.

```bash
RepeatMasker -species human -pa 8 -dir rm_trf_check unknown_consensi.fa
```

Remove families where TRF / Simple_repeat / Satellite accounts for ≥50% of the consensus length. Record excluded families in `is_satellite_or_simple.txt`; remaining families go to `unknown_nonsatellite_consensi.fa`.

> Use 30% as the threshold if high satellite contamination is expected (e.g. acrocentric chromosomes).

---

## Step 4 — Filter against RepeatMasker library (BLASTN)

```bash
makeblastdb -in <RepeatMasker.lib> -dbtype nucl -out rmlib_db

blastn -query unknown_nonsatellite_consensi.fa \
  -db rmlib_db \
  -outfmt "6 qseqid sseqid pident length evalue bitscore" \
  -perc_identity 60 -evalue 1e-5 -word_size 7 \
  -dust no \
  -num_threads 8 \
  > unknown_vs_RMlib.blast.out
```

Families with hits → `has_hit.txt` (removed). Remaining → `novel_consensi.fa`.

> **Use `-dust no`.** With DUST on, low-complexity sequences produce zero hits — not because they are novel, but because the aligner refuses to align them. Step 3b removes clear-cut cases; `-dust no` here catches any residual.

> This BLAST checks against the curated RepeatMasker flat-file library only. Diverged HERVs, solo LTRs, and degenerate fragments are better covered by the full Dfam search in Step 9a.

---

## Step 5 — Genomic copy confirmation (BLASTN)

Confirm each candidate family is genuinely repetitive (multiple copies in the target chromosome).

```bash
blastn -query novel_consensi.fa -db chrNdb -outfmt 6 -num_threads 8 \
  > novel_vs_chrN.blast.out
```

Families with ≥1 hit → `has_genomic_hit.txt`. Families with no hits are excluded.

---

## Step 6 — Length filtering and NCBI nt BLAST (optional)

Extract long families (above a length threshold) to `novel_long_consensi.fa` for more sensitive characterisation.

```bash
blastn -query novel_long_consensi.fa -db nt -remote -outfmt 6 \
  > ncbi_nt_blast.out
```

Families with nt hits → `has_nt_hit.txt`.

---

## Step 7 — TE protein BLASTX

```bash
blastx -query novel_consensi.fa -db RepeatPeps.lib \
  -outfmt "6 qseqid sseqid pident length evalue" -num_threads 8 \
  > novel_vs_repeatpeps.blastx.out
```

Detects coding domains and assists structural classification.

---

## Step 8 — RepeatMasker re-annotation with novel library

Quantify how much of the chromosome the candidate families annotate.

```bash
RepeatMasker -lib novel_consensi.fa -pa 8 -nolow chrN.fa
```

Outputs in `rm_novel_check/`: `chrN.fa.out`, `chrN.fa.tbl`.

> **Interpret masking figures with caution.** LTR/ERV elements are ~8% of the entire human genome; a handful of candidate families masking double-digit percentages of a single chromosome is a strong signal of artefact (diverged known TEs or p-arm HOR/SD content), not genuine discovery. All candidates must pass Steps 9a–9d before any masking figure is reported.

---

## Steps 9a–9d — Validation (required before any novelty claim)

A family is a real candidate only if it survives all four tests.

### 9a — Full Dfam (human clade) BLAST

```bash
famdb.py -i "$RMLIB" families \
  --format fasta_name --include-class-in-name \
  --ancestors --descendants 'Homo sapiens' > human-dfam.fa
makeblastdb -in human-dfam.fa -dbtype nucl -out human_dfam_db

blastn -query novel_consensi.fa -db human_dfam_db \
  -task dc-megablast -evalue 1e-5 -max_target_seqs 5 \
  -outfmt '6 qseqid sseqid pident length evalue bitscore stitle' \
  > novel_vs_dfam.tsv
```

### 9b — Profile-HMM search (catches diverged LTRs BLAST misses)

```bash
famdb.py -i "$RMLIB" families \
  --format hmm --ancestors --descendants 'Homo sapiens' > human-dfam.hmm

nhmmer --cpu 8 --tblout novel_vs_dfam_hmm.tbl -o /dev/null \
  human-dfam.hmm novel_consensi.fa
```

### 9c — Overlap with existing annotation

```bash
RepeatMasker -pa 8 -lib novel_consensi.fa -e rmblast chrN.fa
mkdir std && RepeatMasker -pa 8 -species human -e rmblast std/chrN.fa

awk 'NR>3{print $5"\t"$6-1"\t"$7"\t"$10}' chrN.fa.out     > novel.bed
awk 'NR>3{print $5"\t"$6-1"\t"$7"\t"$10}' std/chrN.fa.out > standard.bed

bedtools intersect -a novel.bed -b standard.bed -wo > novel_overlap.tsv
```

### 9d — Genomic distribution and structural sanity

Bin `novel.bed` by position. Families clustering in the short arm / rDNA / satellite block → flag as satellite/SD-derived pending structural proof. For LTR-called families, verify a real pair of terminal direct repeats + flanking TSDs (+ PBS/PPT for internal regions).

> **Local tandem array check (failure mode 4, first confirmed chr19).** A family can pass Steps 9a–9c (no Dfam hit, low annotation overlap) and still not be a TE — it can be a locally amplified tandem sequence that RepeatModeler mistook for a dispersed repeat family. Compute the inter-copy spacing distribution from `novel.bed` per family: sort copies by position, take consecutive gaps. **If >90% of a family's RM copies fall in a single ~5 Mb bin AND the median inter-copy gap is <5 kb, the family is a local tandem array, not a TE** — discard regardless of other test results. Also check whether RM copy count is much greater (>3×) than the independent BLAST genomic-hit count (Step 5) — that signals RepeatMasker fragmenting single copies into multiple intervals, not independent insertions; inspect manually before counting copies.

---

## Step 10 — Segmental duplication overlap check

```bash
wget -q "https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/annotation/chm13v2.0_SD.bed"

awk 'NR>3{print $5"\t"$6-1"\t"$7"\t"$10}' chrN.fa.out > novel.bed
bedtools intersect -a novel.bed -b chm13v2.0_SD.bed -wo > novel_sd_overlap.tsv
```

Families with the majority of copies inside SD blocks are SD-derived, not genuine TE dispersal.

---

## Step 11 — Translation screen (required before any novelty claim) ⚠️ added after chr19

**Failure mode 5 (first confirmed chr19): host gene exon collapse.** Conserved exons shared across a large paralogous gene family (ZNF clusters, olfactory receptor clusters, immunoglobulin loci) are picked up as repeat families by RepeatModeler. RepeatClassifier returns `#Unknown` because it has no TE model for host-gene domains, and the consensus passes Steps 9a–9c cleanly because the domain simply isn't in any TE database. This is invisible to every nucleotide-level test above — it requires translating the consensus and checking against a protein database.

Run on every family still standing after Step 10 (i.e. the CANDIDATE_NOVEL set):

```bash
# 6-frame translation, longest ORFs first
esl-translate novel_consensi.fa > novel_translated.fa

# Check against curated proteins (remote NCBI; minutes for 5-10 candidates)
blastx -query novel_consensi.fa -db swissprot -remote \
  -outfmt '6 qseqid sseqid pident length evalue bitscore stitle' \
  > novel_vs_swissprot.tsv

# Already-available local check vs TE protein domains (same as Step 7, rerun on survivors)
blastx -query novel_consensi.fa -db RepeatPeps.lib \
  -outfmt "6 qseqid sseqid pident length evalue bitscore" -num_threads 8 \
  > novel_vs_repeatpeps_final.blastx.out
```

**Diagnostic shortcuts:**
- A full-length, uninterrupted ORF spanning (or nearly spanning) the entire consensus → almost certainly a gene exon, not a TE. Confirm with the SwissProt blastx.
- `DVMLENY` (or the broader `QRNLYRDVMLENY` core) anywhere in a 6-frame translation → KRAB-A box; the family is a ZNF gene fragment. Discard without waiting for the BLAST job.
- All Step 5/9d genomic hits concentrated in a 1–3 Mb window that overlaps a known tandemly-duplicated gene cluster (ZNF, OR, Ig V/D/J) → check gene annotation before proceeding regardless of ORF length.
- Chromosomes at elevated risk for this failure mode: chr19 (ZNF + OR clusters), chr17 (ZNF cluster), chr11 (largest OR cluster in the genome), chr1 (OR cluster), chr6 (MHC + OR), and any chromosome with a large tandemly duplicated gene family.

A family that has no SwissProt hit, no RepeatPeps hit, and no full-length spanning ORF clears Step 11.

---

## Step 12 — Record results in `chromosome_results.md` (final step)

**The pipeline is not complete until the outcome is recorded.** This is the last step of every per-chromosome run, executed whether or not any CANDIDATE_NOVEL families survived. After Step 11 (or directly after Step 10 if 0 CANDIDATE_NOVEL families exist, so Step 11 is skipped), append a `## ChrN` section to `chromosome_results.md` containing:

- Filtering funnel: RM families → Unknown → Step 3b removals → Step 4 removals → candidate input → Steps 9a–9d+10 validation outcome → final calls
- RepeatModeler round breakdown table
- Validation summary table (call counts: CONFIRMED_KNOWN, SD_DERIVED, structural discard, AMBIGUOUS, CANDIDATE_NOVEL, etc.)
- Per-family detail table for every non-CONFIRMED_KNOWN family (BLAST/HMM/overlap/SD/tandem evidence)
- Masking figure from Step 8 (`rm_novel_check/chrN.fa.tbl`)
- Any structural follow-up notes for CANDIDATE_NOVEL survivors
- Note any deviation from the documented pipeline (e.g. a step run informationally rather than as a hard gate)

Then update the **Cross-chromosome summary** table and the running chromosome list in `chromosome_results.md` with the new chromosome's numbers.

No chromosome run should be reported as "done," and no claim language should be used, until this step is complete.

---

## Decision criteria

| Outcome | Criteria |
|---|---|
| **Satellite/simple repeat** (discard) | ≥50% of consensus masked by TRF / Simple_repeat / Satellite in Step 3b |
| **Confirmed known** (discard) | BLAST hit in full Dfam (Step 9a, E ≤ 1e-5) **or** significant nhmmer hit (Step 9b) **or** ≥50% genomic footprint overlaps existing annotation (Step 9c) **or** RepeatPeps protein-domain hit at E ≤ 1e-5 (Step 11 — protein-level homology is more sensitive than nucleotide BLAST for old, diverged TE remnants, so a RepeatPeps hit here is corroborating TE evidence, not host-gene evidence) |
| **SD-derived** (discard) | ≥50% of genomic copies fall inside segmental-duplication blocks (Step 10) |
| **Local tandem array** (discard) | >90% of RM copies in a single ~5 Mb bin **and** median inter-copy gap <5 kb (Step 9d) |
| **Host gene exon collapse** (discard) | Full-length ORF spanning consensus with a SwissProt hit (E ≤ 1e-5), **or** a recognisable domain signature (e.g. KRAB-A `DVMLENY`) in 6-frame translation (Step 11) — distinct from a RepeatPeps hit, which argues the opposite (real TE) |
| **Candidate novel** | Passes Step 3b; no Dfam BLAST or HMM hit; no RepBase hit; <10% overlap with existing annotation; <50% copies in SD; not a local tandem array; no SwissProt hit and no full-length spanning ORF; no RepeatPeps hit (Step 11); consistent TE structure (real LTR pair + TSDs) where applicable; multiple genomic copies; consensus length >100 bp |

