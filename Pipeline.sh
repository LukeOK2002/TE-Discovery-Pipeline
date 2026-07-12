#!/usr/bin/env bash
# T2T per-chromosome de novo repeat discovery pipeline (see pipeline.md)
# Usage: bash run_pipeline.sh <chrom>   (launch inside tmux — this is long-running)
set -euo pipefail

CHROM="${1:?Usage: run_pipeline.sh <chrom> (e.g. chr9)}"
WORKDIR="/home/lukeok/RepeatFinder/T2T/${CHROM}"
cd "$WORKDIR"

source /home/lukeok/miniforge3/etc/profile.d/conda.sh
conda activate repeatmodeler_env

# ── paths ──────────────────────────────────────────────────────────────────────
FAMDB=/home/lukeok/miniforge3/envs/repeatmodeler_env/share/RepeatMasker/famdb.py
RMLIB_DIR=/home/lukeok/miniforge3/envs/repeatmodeler_env/share/RepeatMasker/Libraries/famdb
RMLIB=/home/lukeok/miniforge3/envs/repeatmodeler_env/share/RepeatMasker/Libraries/RepeatMasker.lib
REPEATPEPS=/home/lukeok/miniforge3/envs/repeatmodeler_env/share/RepeatMasker/Libraries/RepeatPeps.lib
T2T_OUT=/home/lukeok/RepeatFinder/T2T/chm13v2.0_RepeatMasker_4.1.2p1.2022Apr14.out
SD_BED=/home/lukeok/RepeatFinder/T2T/chr22/validation/chm13v2.0_SD.bed

CPUS=8
BIN_SIZE=5000000
PAR_END=0          # p-arm/q-arm split coord for Step 9d; 0 if non-acrocentric

# ── centromere-exclusion fallback ───────────────────────────────────────────
# Set CEN_EXCLUDE=1 on chromosomes where soft-masked satellite alone still lets
# RECON's all-vs-all explode (e.g. under-annotated peri/centromeric HOR arrays).
# Hard-masks [CEN_START, CEN_END) to N in addition to Step 1d2; coordinates are
# otherwise preserved so downstream BED/coordinate steps are unaffected.
CEN_EXCLUDE=0
CEN_START=0
CEN_END=0

FNA="${WORKDIR}/${CHROM}.fna"
FA="${WORKDIR}/${CHROM}.fa"
MASKED_FA="${WORKDIR}/${CHROM}_masked.fa"
RMSK_BED="${WORKDIR}/${CHROM}_t2t_rmsk.bed"
GENOME="${WORKDIR}/${CHROM}.genome"
DB="${WORKDIR}/${CHROM}db"
LOG="${WORKDIR}/${CHROM}_RM.log"

# per-chromosome intermediate filenames
UNKNOWN_FA="${CHROM}_unknown_consensi.fa"
SAT_LIST="${CHROM}_is_satellite_or_simple.txt"
NONSAT_FA="${CHROM}_unknown_nonsatellite_consensi.fa"
RMLIB_BLAST="${CHROM}_unknown_vs_rmlib.blast"
HAS_HIT_TXT="${CHROM}_has_hit.txt"
NOVEL_FA="${CHROM}_novel_consensi.fa"
GENOMIC_BLAST="${CHROM}_novel_vs_${CHROM}.blast"
HAS_GENOMIC_HIT_TXT="${CHROM}_has_genomic_hit.txt"
REPEATPEPS_PRELIM="${CHROM}_novel_vs_repeatpeps.blastx.out"
CANDIDATE_FA="${CHROM}_candidate_novel_consensi.fa"
RM_NOVEL_OUT="rm_novel_check/${CHROM}.fa.out"

log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*"; }

# ──────────────────────────────────────────────────────────────────────────────
# Step 1 — Genome preparation
# ──────────────────────────────────────────────────────────────────────────────
log "Step 1a: Fixing FASTA header → >${CHROM}"
sed "1s/^>.*/>${CHROM}/" "$FNA" > "$FA"

log "Step 1b: Indexing FASTA and creating .genome"
samtools faidx "$FA"
cut -f1,2 "${FA}.fai" > "$GENOME"

log "Step 1c: Extracting ${CHROM} repeats from T2T .out → BED"
awk -v chrom="$CHROM" '
  NF >= 15 && $5 == chrom {
    start = $6 - 1; end = $7
    if (start < 0) start = 0
    print chrom "\t" start "\t" end
  }
' "$T2T_OUT" | sort -k1,1 -k2,2n > "$RMSK_BED"
log "  $(wc -l < "$RMSK_BED") repeat intervals extracted"

log "Step 1d: Soft-masking ${CHROM} with T2T repeats"
bedtools maskfasta -soft -fi "$FA" -bed "$RMSK_BED" -fo "$MASKED_FA"

# ── Step 1d2 — hard-mask tandem/low-complexity classes (→ N) ──────────────────
# RepeatModeler reads soft-masked (lowercase) sequence as ordinary input, so visible
# satellite/simple-repeat arrays can drive an O(n^2) HSP blowup in RECON. These
# classes aren't discovery targets anyway (Step 3b drops them), so hard-mask to N;
# dispersed TEs stay soft-masked and visible.
HARDMASK_BED="${WORKDIR}/${CHROM}_hardmask.bed"
log "Step 1d2: Hard-masking Satellite*/Simple_repeat/Low_complexity → N"
awk -v chrom="$CHROM" '
  NF >= 15 && $5 == chrom && ($11 ~ /Satellite/ || $11 == "Simple_repeat" || $11 == "Low_complexity") {
    start = $6 - 1; end = $7
    if (start < 0) start = 0
    print chrom "\t" start "\t" end
  }
' "$T2T_OUT" | sort -k1,1 -k2,2n > "$HARDMASK_BED"

if (( CEN_EXCLUDE == 1 )); then
  log "  CEN_EXCLUDE=1 → also hard-masking centromere ${CHROM}:${CEN_START}-${CEN_END}"
  printf '%s\t%s\t%s\n' "$CHROM" "$CEN_START" "$CEN_END" >> "$HARDMASK_BED"
  sort -k1,1 -k2,2n "$HARDMASK_BED" | bedtools merge -i - > "${HARDMASK_BED}.merged"
  mv "${HARDMASK_BED}.merged" "$HARDMASK_BED"
fi

log "  $(wc -l < "$HARDMASK_BED") intervals, $(awk '{s+=$3-$2} END{print s}' "$HARDMASK_BED") bp → N"
# hard mask (default = replace with N); preserves existing soft-mask case elsewhere
bedtools maskfasta -fi "$MASKED_FA" -bed "$HARDMASK_BED" -fo "${MASKED_FA}.hm"
mv "${MASKED_FA}.hm" "$MASKED_FA"
log "  hard-masked N content now: $(awk '!/^>/{for(i=1;i<=length($0);i++){c=substr($0,i,1); if(c==\"N\"||c==\"n\")n++; t++}} END{printf \"%.2f%%\", 100*n/t}' "$MASKED_FA")"

# ── disk guard ──────────────────────────────────────────────────────────────
# RECON's element DB can blow up to 100s of GB on satellite-dense chromosomes.
# Abort early rather than filling the disk.
AVAIL_GB=$(df -BG --output=avail "$WORKDIR" | tail -1 | tr -dc '0-9')
MIN_FREE_GB=200
if (( AVAIL_GB < MIN_FREE_GB )); then
  log "ERROR: <${MIN_FREE_GB} GB free — RECON risk. Free space or narrow input, then re-run."
  exit 1
fi
log "Free disk: ${AVAIL_GB} GB. During RECON, monitor: du -sh ${WORKDIR}/RM_*/round-*/ and df -h ${WORKDIR}."

log "Step 1e: Building RepeatModeler database"
BuildDatabase -name "$DB" "$MASKED_FA"
log "Step 1: done."

# ──────────────────────────────────────────────────────────────────────────────
# Step 2 — RepeatModeler (long-running; tmux keeps this alive)
# ──────────────────────────────────────────────────────────────────────────────
log "Step 2: Starting RepeatModeler (threads=$CPUS, -LTRStruct)..."
RepeatModeler -database "$DB" -threads "$CPUS" -LTRStruct 2>&1 | tee "$LOG"
log "Step 2: RepeatModeler done."

FAMILIES="${WORKDIR}/${CHROM}db-families.fa"
[[ -f "$FAMILIES" ]] || { log "ERROR: $FAMILIES not found. Check RM_* output."; exit 1; }

# ──────────────────────────────────────────────────────────────────────────────
# Step 3 — Extract Unknown families
# ──────────────────────────────────────────────────────────────────────────────
log "Step 3: Extracting #Unknown families from $FAMILIES..."
awk '/^>/{keep=/#Unknown/} keep{print}' "$FAMILIES" > "$UNKNOWN_FA"
NUNKNOWN=$(grep -c "^>" "$UNKNOWN_FA")
log "Step 3: $NUNKNOWN Unknown families extracted."

# ──────────────────────────────────────────────────────────────────────────────
# Step 3b — TRF / satellite pre-filter
# ──────────────────────────────────────────────────────────────────────────────
log "Step 3b: Running RepeatMasker on Unknown consensi (TRF/satellite detection)..."
mkdir -p rm_trf_check
RepeatMasker -species human -pa "$CPUS" -dir rm_trf_check "$UNKNOWN_FA"

log "Step 3b: Identifying satellite/simple-repeat families (>=50% masked)..."
python3 - <<PYEOF
import os

out_file = "rm_trf_check/${UNKNOWN_FA}.out"
fasta_file = "${UNKNOWN_FA}"

seq_len = {}
name = None
buf = []
with open(fasta_file) as f:
    for line in f:
        line = line.rstrip()
        if line.startswith(">"):
            if name:
                seq_len[name] = sum(len(s) for s in buf)
            name = line[1:].split()[0]
            buf = []
        else:
            buf.append(line)
    if name:
        seq_len[name] = sum(len(s) for s in buf)

masked_bp = {}
if os.path.exists(out_file):
    with open(out_file) as f:
        for i, line in enumerate(f):
            if i < 3: continue
            parts = line.split()
            if len(parts) < 11: continue
            qname = parts[4]
            try:
                qstart = int(parts[5]); qend = int(parts[6])
            except ValueError:
                continue
            repeat_class = parts[10]
            if any(c in repeat_class for c in ("Satellite", "Simple_repeat", "Low_complexity")):
                masked_bp[qname] = masked_bp.get(qname, 0) + (qend - qstart + 1)

satellite_families = [fam for fam, length in seq_len.items()
                      if length > 0 and masked_bp.get(fam, 0) / length >= 0.50]

with open("${SAT_LIST}", "w") as f:
    for fam in satellite_families:
        f.write(fam + "\n")
print(f"  {len(satellite_families)} families flagged as satellite/simple repeat (>=50% masked)")

sat_set = set(satellite_families)
kept = 0
with open("${UNKNOWN_FA}") as fin, \
     open("${NONSAT_FA}", "w") as fout:
    keep = False
    for line in fin:
        if line.startswith(">"):
            name = line[1:].split()[0]
            keep = name not in sat_set
        if keep:
            fout.write(line)
            if line.startswith(">"): kept += 1
print(f"  {kept} non-satellite families written to ${NONSAT_FA}")
PYEOF
log "Step 3b: done."

# ──────────────────────────────────────────────────────────────────────────────
# Step 4 — BLAST against RepeatMasker library (-dust no)
# ──────────────────────────────────────────────────────────────────────────────
log "Step 4: BLASTn against RepeatMasker library (dust no, perc_identity 60, evalue 1e-5)..."
blastn -query "$NONSAT_FA" \
  -db "$RMLIB" \
  -outfmt "6 qseqid sseqid pident length evalue bitscore" \
  -perc_identity 60 -evalue 1e-5 -word_size 7 \
  -dust no -num_threads "$CPUS" \
  > "$RMLIB_BLAST"

log "Step 4: $(wc -l < "$RMLIB_BLAST") hits. Partitioning..."
cut -f1 "$RMLIB_BLAST" | sort -u > "$HAS_HIT_TXT"

python3 - <<PYEOF
hit_set = set(l.strip() for l in open("${HAS_HIT_TXT}"))
kept = 0
with open("${NONSAT_FA}") as fin, \
     open("${NOVEL_FA}", "w") as fout:
    keep = False
    for line in fin:
        if line.startswith(">"):
            name = line[1:].split()[0]
            keep = name not in hit_set
        if keep:
            fout.write(line)
            if line.startswith(">"): kept += 1
print(f"  {len(hit_set)} families removed (RM library hit)")
print(f"  {kept} families in ${NOVEL_FA}")
PYEOF
log "Step 4: done."

# ──────────────────────────────────────────────────────────────────────────────
# Step 5 — Genomic copy confirmation
# ──────────────────────────────────────────────────────────────────────────────
log "Step 5: BLASTn novel families vs ${CHROM}db (genomic copy confirmation)..."
blastn -query "$NOVEL_FA" \
  -db "$DB" -outfmt 6 -num_threads "$CPUS" \
  > "$GENOMIC_BLAST"

cut -f1 "$GENOMIC_BLAST" | sort -u > "$HAS_GENOMIC_HIT_TXT"

python3 - <<PYEOF
genomic = set(l.strip() for l in open("${HAS_GENOMIC_HIT_TXT}"))
total = sum(1 for l in open("${NOVEL_FA}") if l.startswith(">"))
print(f"  {len(genomic)}/{total} families have genomic copies in ${CHROM}")
print(f"  {total - len(genomic)} families excluded (no genomic copies)")
PYEOF
log "Step 5: done."

# ──────────────────────────────────────────────────────────────────────────────
# Step 7 — TE protein BLASTX (preliminary; also run on survivors in Step 11)
# ──────────────────────────────────────────────────────────────────────────────
log "Step 7: blastx vs RepeatPeps.lib..."
blastx -query "$NOVEL_FA" -db "$REPEATPEPS" \
  -outfmt "6 qseqid sseqid pident length evalue bitscore" -num_threads "$CPUS" \
  > "$REPEATPEPS_PRELIM"
log "Step 7: $(wc -l < "$REPEATPEPS_PRELIM") hits."

# ──────────────────────────────────────────────────────────────────────────────
# Step 8 — RepeatMasker re-annotation with novel library
# ──────────────────────────────────────────────────────────────────────────────
log "Step 8: RepeatMasker with novel library..."
mkdir -p rm_novel_check
RepeatMasker -lib "$NOVEL_FA" -pa "$CPUS" -nolow -dir rm_novel_check "$FA"
log "Step 8: done."

# ──────────────────────────────────────────────────────────────────────────────
# Steps 9a–9d, 10 — Validation
# ──────────────────────────────────────────────────────────────────────────────
mkdir -p validation

log "9a: Exporting human-clade Dfam FASTA..."
python "$FAMDB" -i "$RMLIB_DIR" families \
  --format fasta_name --include-class-in-name \
  --ancestors --descendants 'Homo sapiens' \
  > validation/human-dfam.fa
log "9a: $(grep -c '>' validation/human-dfam.fa) sequences exported."

log "9a: Building BLAST database..."
makeblastdb -in validation/human-dfam.fa -dbtype nucl \
  -out validation/human_dfam_db -logfile validation/makeblastdb.log

log "9a: BLAST dc-megablast (evalue 1e-5)..."
blastn -query "$NOVEL_FA" -db validation/human_dfam_db \
  -task dc-megablast -evalue 1e-5 -max_target_seqs 5 -num_threads "$CPUS" \
  -outfmt '6 qseqid sseqid pident length evalue bitscore stitle' \
  > validation/novel_vs_dfam.tsv
log "9a: $(wc -l < validation/novel_vs_dfam.tsv) hits (dc-megablast)."

log "9a: Sensitive BLAST fallback (blastn, word_size 7, evalue 1e-3)..."
blastn -query "$NOVEL_FA" -db validation/human_dfam_db \
  -task blastn -word_size 7 -evalue 1e-3 -max_target_seqs 5 -num_threads "$CPUS" \
  -outfmt '6 qseqid sseqid pident length evalue bitscore stitle' \
  > validation/novel_vs_dfam_sensitive.tsv
log "9a: $(wc -l < validation/novel_vs_dfam_sensitive.tsv) hits (sensitive)."

log "9b: Exporting human-clade Dfam HMMs..."
python "$FAMDB" -i "$RMLIB_DIR" families \
  --format hmm --ancestors --descendants 'Homo sapiens' \
  > validation/human-dfam.hmm
log "9b: HMM export done."

log "9b: nhmmer (cpu=$CPUS)..."
nhmmer --cpu "$CPUS" --tblout validation/novel_vs_dfam_hmm.tbl -o /dev/null \
  validation/human-dfam.hmm "$NOVEL_FA"
log "9b: $(grep -v '^#' validation/novel_vs_dfam_hmm.tbl | wc -l) nhmmer hits."

log "9c: Converting novel RM output to BED..."
awk 'NR>3 && NF>=10 {print $5"\t"$6-1"\t"$7"\t"$10}' \
  "$RM_NOVEL_OUT" > validation/novel.bed
log "9c: $(wc -l < validation/novel.bed) novel intervals."

log "9c: Extracting ${CHROM} from T2T standard annotation..."
awk -v chrom="$CHROM" 'NR>3 && NF>=10 && $5==chrom {print $5"\t"$6-1"\t"$7"\t"$10}' \
  "$T2T_OUT" > validation/standard.bed
log "9c: $(wc -l < validation/standard.bed) standard ${CHROM} intervals."

log "9c: bedtools intersect..."
bedtools intersect -a validation/novel.bed -b validation/standard.bed -wo \
  > validation/novel_overlap.tsv
log "9c: $(wc -l < validation/novel_overlap.tsv) overlapping intervals."

log "9d: Genomic distribution + tandem-array spacing check..."
python3 - <<PYEOF
import csv
from collections import defaultdict
from statistics import median

PAR_END = $PAR_END
BIN_SIZE = $BIN_SIZE

intervals = defaultdict(list)
with open("validation/novel.bed") as f:
    for row in csv.reader(f, delimiter="\t"):
        if len(row) < 4: continue
        fam = row[3]
        start, end = int(row[1]), int(row[2])
        intervals[fam].append((start, end))

rows = []
for fam, ivs in intervals.items():
    ivs.sort()
    n = len(ivs)
    parm_bp = sum(e - s for s, e in ivs if e <= PAR_END)
    qarm_bp = sum(e - s for s, e in ivs if e > PAR_END)

    bin_counts = defaultdict(int)
    for s, e in ivs:
        bin_counts[s // BIN_SIZE] += 1
    max_bin_frac = max(bin_counts.values()) / n if n > 0 else 0.0

    gaps = [max(0, ivs[i][0] - ivs[i-1][1]) for i in range(1, n)]
    median_gap = median(gaps) if gaps else -1

    tandem_flag = (n >= 3) and (max_bin_frac > 0.90) and (median_gap >= 0) and (median_gap < 5000)

    rows.append([fam, n, parm_bp, qarm_bp, round(max_bin_frac, 3),
                 int(median_gap) if median_gap >= 0 else "NA", int(tandem_flag)])

rows.sort(key=lambda r: -r[2])
with open("validation/novel_distribution.tsv", "w", newline="") as f:
    w = csv.writer(f, delimiter="\t")
    w.writerow(["family", "interval_count", "bp_parm", "bp_qarm",
                "max_bin_frac", "median_gap_bp", "tandem_array_flag"])
    w.writerows(rows)

n_tandem = sum(r[6] for r in rows)
print(f"  {len(rows)} families summarised; {n_tandem} flagged as local tandem arrays")
PYEOF
log "9d: Distribution + tandem-array summary written."

log "10: Segmental duplication overlap..."
bedtools intersect -a validation/novel.bed -b "$SD_BED" -wo \
  > validation/novel_sd_overlap.tsv
bedtools intersect -a validation/novel.bed -b "$SD_BED" -u \
  > validation/novel_sd_hit_intervals.bed
log "10: $(wc -l < validation/novel_sd_overlap.tsv) overlapping records; $(wc -l < validation/novel_sd_hit_intervals.bed) distinct SD-overlapping intervals."

# ──────────────────────────────────────────────────────────────────────────────
# Per-family decision table v1 (pre-translation-screen)
# ──────────────────────────────────────────────────────────────────────────────
log "Producing per-family decision table (v1)..."
python3 - <<PYEOF
import csv, os
from collections import defaultdict, Counter

outdir = "validation"

blast_best = {}
for fname in (f"{outdir}/novel_vs_dfam.tsv", f"{outdir}/novel_vs_dfam_sensitive.tsv"):
    if os.path.exists(fname):
        with open(fname) as f:
            for row in csv.reader(f, delimiter="\t"):
                if len(row) < 6: continue
                qid = row[0]
                if qid not in blast_best:
                    blast_best[qid] = (row[1], float(row[4]), float(row[5]))

hmm_best = {}
if os.path.exists(f"{outdir}/novel_vs_dfam_hmm.tbl"):
    with open(f"{outdir}/novel_vs_dfam_hmm.tbl") as f:
        for line in f:
            if line.startswith("#"): continue
            parts = line.split()
            if len(parts) < 14: continue
            tgt = parts[0]
            try:
                evalue = float(parts[12]); score = float(parts[13])  # nhmmer tblout: E-value=col13, score=col14 (1-indexed)
            except (ValueError, IndexError):
                continue
            if tgt not in hmm_best or evalue < hmm_best[tgt][1]:
                hmm_best[tgt] = (parts[2], evalue, score)

family_bp = defaultdict(int)
family_olap = defaultdict(int)
if os.path.exists(f"{outdir}/novel.bed"):
    with open(f"{outdir}/novel.bed") as f:
        for row in csv.reader(f, delimiter="\t"):
            if len(row) < 4: continue
            family_bp[row[3]] += int(row[2]) - int(row[1])
if os.path.exists(f"{outdir}/novel_overlap.tsv"):
    with open(f"{outdir}/novel_overlap.tsv") as f:
        for row in csv.reader(f, delimiter="\t"):
            if len(row) < 9: continue
            family_olap[row[3]] += int(row[-1])

dist = {}; tandem = {}
if os.path.exists(f"{outdir}/novel_distribution.tsv"):
    with open(f"{outdir}/novel_distribution.tsv") as f:
        for row in csv.reader(f, delimiter="\t"):
            if row[0] == "family" or len(row) < 7: continue
            fam, parm_bp, qarm_bp = row[0], int(row[2]), int(row[3])
            loc = "parm" if parm_bp > qarm_bp else "qarm"
            if parm_bp > 0 and qarm_bp > 0: loc = "both"
            dist[fam] = loc
            tandem[fam] = bool(int(row[6]))

sd_total = defaultdict(int); sd_hit = defaultdict(int)
if os.path.exists(f"{outdir}/novel.bed"):
    with open(f"{outdir}/novel.bed") as f:
        for row in csv.reader(f, delimiter="\t"):
            if len(row) < 4: continue
            sd_total[row[3]] += 1
if os.path.exists(f"{outdir}/novel_sd_hit_intervals.bed"):
    with open(f"{outdir}/novel_sd_hit_intervals.bed") as f:
        for row in csv.reader(f, delimiter="\t"):
            if len(row) < 4: continue
            sd_hit[row[3]] += 1

families = [l[1:].split()[0] for l in open("${NOVEL_FA}") if l.startswith(">")]

BLAST_THRESH = 1e-5; HMM_THRESH = 1e-5
rows = []
for fam in families:
    fam_short = fam.split("#")[0]
    b = blast_best.get(fam) or blast_best.get(fam_short)
    h = hmm_best.get(fam)  or hmm_best.get(fam_short)

    blast_hit  = b[0] if b else "none"; blast_eval = b[1] if b else 1.0
    hmm_hit    = h[0] if h else "none"; hmm_eval   = h[1] if h else 1.0

    bp_total = family_bp.get(fam_short, 0)
    bp_olap  = family_olap.get(fam_short, 0)
    pct_olap = round(100 * bp_olap / bp_total, 1) if bp_total > 0 else 0.0

    sdt = sd_total.get(fam_short, 0); sdh = sd_hit.get(fam_short, 0)
    pct_sd = round(100 * sdh / sdt, 1) if sdt > 0 else 0.0

    loc = dist.get(fam_short, "unknown")
    is_tandem = tandem.get(fam_short, False)

    known_blast = blast_eval <= BLAST_THRESH
    known_hmm   = hmm_eval   <= HMM_THRESH
    known_olap  = pct_olap   >= 50.0
    known_sd    = pct_sd     >= 50.0

    if known_blast or known_hmm or known_olap:
        call = "CONFIRMED_KNOWN"
    elif known_sd:
        call = "SD_DERIVED"
    elif is_tandem:
        call = "TANDEM_ARRAY"
    elif pct_olap < 10.0:
        call = "CANDIDATE_NOVEL"
    else:
        call = "AMBIGUOUS"

    if loc in ("parm", "both"):
        call += "+PARM_FLAG"

    rows.append([fam, blast_hit, f"{blast_eval:.2e}", hmm_hit, f"{hmm_eval:.2e}",
                 f"{pct_olap:.1f}", f"{pct_sd:.1f}", loc, int(is_tandem), call])

with open(f"{outdir}/family_calls.tsv", "w", newline="") as f:
    w = csv.writer(f, delimiter="\t")
    w.writerow(["family_id","best_dfam_blast","blast_evalue",
                "best_hmm_hit","hmm_evalue","pct_overlap_existing",
                "pct_copies_in_sd","location_class","tandem_array_flag","call"])
    w.writerows(rows)

calls = Counter(r[9].split("+")[0] for r in rows)
print("\n=== Per-family call summary (v1, pre-translation-screen) ===")
for k, v in sorted(calls.items()):
    print(f"  {v:3d}  {k}")
print(f"\nFull table: {outdir}/family_calls.tsv")
PYEOF

# ──────────────────────────────────────────────────────────────────────────────
# Step 11 — Translation screen on CANDIDATE_NOVEL survivors
# Runs on every CANDIDATE_NOVEL survivor; inspect hits for host-gene domain
# signatures before claiming novelty.
# ──────────────────────────────────────────────────────────────────────────────
log "Step 11: Extracting CANDIDATE_NOVEL survivors for translation screen..."
python3 - <<PYEOF
import csv
survivors = set()
with open("validation/family_calls.tsv") as f:
    r = csv.reader(f, delimiter="\t"); next(r)
    for row in r:
        if row[-1].split("+")[0] == "CANDIDATE_NOVEL":
            survivors.add(row[0])

kept = 0
with open("${NOVEL_FA}") as fin, \
     open("${CANDIDATE_FA}", "w") as fout:
    keep = False
    for line in fin:
        if line.startswith(">"):
            name = line[1:].split()[0]; keep = name in survivors
        if keep:
            fout.write(line)
            if line.startswith(">"): kept += 1
print(f"  {kept} CANDIDATE_NOVEL families written to ${CANDIDATE_FA}")
PYEOF

NSURVIVORS=$(grep -c "^>" "$CANDIDATE_FA" || true)
log "Step 11: $NSURVIVORS survivors entering translation screen."

if [[ "$NSURVIVORS" -gt 0 ]]; then
  log "Step 11: 6-frame translation (esl-translate)..."
  esl-translate "$CANDIDATE_FA" > validation/novel_translated.fa

  log "Step 11: blastx vs RepeatPeps.lib (local, TE coding domains)..."
  blastx -query "$CANDIDATE_FA" -db "$REPEATPEPS" \
    -outfmt "6 qseqid sseqid pident length evalue bitscore" -num_threads "$CPUS" \
    > validation/novel_vs_repeatpeps_final.blastx.out

  log "Step 11: blastx -remote -db swissprot (NCBI; may take several minutes)..."
  set +e
  blastx -query "$CANDIDATE_FA" -db swissprot -remote \
    -outfmt '6 qseqid sseqid pident length evalue bitscore stitle' \
    > validation/novel_vs_swissprot.tsv 2> validation/swissprot_remote.log
  SP_STATUS=$?
  set -e
  if [[ $SP_STATUS -ne 0 ]]; then
    log "Step 11: WARNING — remote SwissProt BLAST failed (see validation/swissprot_remote.log). Re-run manually."
    : > validation/novel_vs_swissprot.tsv
  else
    log "Step 11: $(wc -l < validation/novel_vs_swissprot.tsv) SwissProt hits."
  fi

  log "Step 11: Updating decision table (v2/final)..."
  python3 - <<'PYEOF'
import csv
from collections import Counter

sp_best = {}
with open("validation/novel_vs_swissprot.tsv") as f:
    for row in csv.reader(f, delimiter="\t"):
        if len(row) < 6: continue
        qid = row[0]; ev = float(row[4])
        if qid not in sp_best or ev < sp_best[qid][1]:
            sp_best[qid] = (row[1], ev, row[-1] if len(row) > 6 else "")

rp_best = {}
with open("validation/novel_vs_repeatpeps_final.blastx.out") as f:
    for row in csv.reader(f, delimiter="\t"):
        if len(row) < 6: continue
        qid = row[0]; ev = float(row[4])
        if qid not in rp_best or ev < rp_best[qid][1]:
            rp_best[qid] = (row[1], ev)

rows = []
with open("validation/family_calls.tsv") as f:
    r = csv.reader(f, delimiter="\t"); header = next(r)
    for row in r: rows.append(row)

header += ["swissprot_hit","swissprot_evalue","repeatpeps_hit","repeatpeps_evalue","final_call"]
for row in rows:
    fam = row[0]; fam_short = fam.split("#")[0]
    call_base = row[-1].split("+")[0]
    sp = sp_best.get(fam) or sp_best.get(fam_short)
    rp = rp_best.get(fam) or rp_best.get(fam_short)
    sp_hit, sp_ev = (sp[0], sp[1]) if sp else ("none", 1.0)
    rp_hit, rp_ev = (rp[0], rp[1]) if rp else ("none", 1.0)

    if call_base != "CANDIDATE_NOVEL":
        final_call = row[-1]
    elif sp_ev <= 1e-5:
        final_call = "HOST_GENE_DISCARD"
    elif rp_ev <= 1e-5:
        final_call = "CONFIRMED_KNOWN_TE_PROTEIN"
    else:
        final_call = row[-1]

    row += [sp_hit, f"{sp_ev:.2e}", rp_hit, f"{rp_ev:.2e}", final_call]

with open("validation/family_calls.tsv", "w", newline="") as f:
    w = csv.writer(f, delimiter="\t")
    w.writerow(header); w.writerows(rows)

calls = Counter(row[-1].split("+")[0] for row in rows)
print("\n=== Final per-family call summary ===")
for k, v in sorted(calls.items()):
    print(f"  {v:3d}  {k}")
final_survivors = [row[0] for row in rows if row[-1].split("+")[0] == "CANDIDATE_NOVEL"]
print(f"\n{len(final_survivors)} families remain CANDIDATE_NOVEL after Step 11:")
for fm in final_survivors:
    print(f"  {fm}")
PYEOF
else
  log "Step 11: no CANDIDATE_NOVEL survivors — skipping translation screen."
fi

log "========================================"
log "Pipeline complete. Key outputs:"
log "  $UNKNOWN_FA                    Step 3"
log "  $SAT_LIST             Step 3b removals"
log "  $NONSAT_FA       post-Step-3b"
log "  $HAS_HIT_TXT                            Step 4 removals"
log "  $NOVEL_FA                      candidate families"
log "  $HAS_GENOMIC_HIT_TXT                    Step 5 confirmed"
log "  $CANDIDATE_FA            Step 11 input"
log "  validation/novel_vs_dfam.tsv                 9a BLAST"
log "  validation/novel_vs_dfam_sensitive.tsv       9a sensitive"
log "  validation/novel_vs_dfam_hmm.tbl             9b nhmmer"
log "  validation/novel_overlap.tsv                 9c annotation overlap"
log "  validation/novel_distribution.tsv            9d distribution + tandem-array flags"
log "  validation/novel_sd_overlap.tsv              10 SD overlap"
log "  validation/novel_translated.fa               11 6-frame translations"
log "  validation/novel_vs_swissprot.tsv            11 SwissProt blastx"
log "  validation/novel_vs_repeatpeps_final.blastx.out 11 RepeatPeps blastx"
log "  validation/family_calls.tsv                  final per-family decisions"
log ""
log "Next: run Step 12 — append $CHROM section to chromosome_results.md"
log "========================================"
