#!/usr/bin/env bash
# =============================================================================
# pipeline.sh
#
# End-to-end training data generation (single simulation run).
# Sequential — designed for 16 GB RAM, one run at a time.
#
# Steps:
#   1. iss generate --ncbi   downloads genomes + simulates reads in one shot
#   2. megahit_topo          builds SDBG (capped at 10 GB), extracts features
#   3. jellyfish             counts true vs error k-mers against reference
#   4. cleanup               deletes reads + FASTA + jf dbs; keeps JSON + abundance
#
# Usage:
#   ./pipeline.sh <abundance.txt> <output_dir> [n_reads] [model] [kingdom] [--keep-reads]
#
# Arguments:
#   abundance.txt   ISS abundance file (accession<TAB>fraction per line)
#   output_dir      directory for <run_id>.json output
#   n_reads         total paired reads to simulate (default: 1000000)
#   model           HiSeq | MiSeq | NovaSeq (default: HiSeq)
#   kingdom         bacteria | archaea | fungi | viral | protists (default: bacteria)
#   --keep-reads    skip cleanup of reads/FASTA (for debugging)
#
# Output JSON:
#   {
#     "run_id": "...", "n_reads_requested": N, "model": "HiSeq",
#     "abundance": "..._abundance.txt",
#     "node": { 10 graph features },
#     "timing": { "node_ms": ... },
#     "kmers": { "n_true": N, "n_error": N, "error_fraction": 0.044 }
#   }
# =============================================================================

set -euo pipefail

# ---------------------------------------------------------------------------
# Arguments
# ---------------------------------------------------------------------------
if [[ $# -lt 2 ]]; then
    echo "Usage: $0 <abundance.txt> <output_dir> [n_reads] [model] [kingdom] [--keep-reads]" >&2
    exit 1
fi

ABUNDANCE_FILE=$(realpath "$1")
OUTPUT_DIR=$(realpath "$2")
N_READS=${3:-1000000}
MODEL=${4:-HiSeq}
KINGDOM=${5:-bacteria}
KEEP_READS=0
for arg in "$@"; do [[ "$arg" == "--keep-reads" ]] && KEEP_READS=1; done

# ---------------------------------------------------------------------------
# Tool paths
# ---------------------------------------------------------------------------
MEGAHIT_TOPO=~/metaTopoGraph-build/megahit_topo
JELLYFISH=jellyfish
ISS=iss
CPUS=$(nproc)

# RAM budget: leave 4 GB for OS + jellyfish, give 10 GB to megahit_topo
TOPO_MEM_GB=10

# ---------------------------------------------------------------------------
# Work area — wiped on exit regardless of success/failure
# ---------------------------------------------------------------------------
mkdir -p "$OUTPUT_DIR"
RUN_ID=$(basename "$ABUNDANCE_FILE" .txt)
WORK=$(mktemp -d /tmp/pipeline_XXXXXX)
trap 'rm -rf "$WORK"' EXIT

PREFIX="$WORK/$RUN_ID"
FASTA="${PREFIX}_ncbi_genomes.fasta"
R1="${PREFIX}_R1.fastq"
R2="${PREFIX}_R2.fastq"
REF_JF="$WORK/ref.jf"
READS_JF="$WORK/reads.jf"
READS_KMERS_FA="$WORK/reads_kmers.fa"
TOPO_JSON="$WORK/topo.json"
FINAL_JSON="$OUTPUT_DIR/${RUN_ID}.json"

# ---------------------------------------------------------------------------
# Step 1: InSilicoSeq — download genomes from NCBI + simulate reads
# ISS --ncbi fetches each accession from the abundance file automatically
# and writes <prefix>_ncbi_genomes.fasta alongside the reads.
# ---------------------------------------------------------------------------
echo "[1/3] iss generate --ncbi  (n_reads=$N_READS, model=$MODEL, kingdom=$KINGDOM)"
$ISS generate \
    --ncbi "$KINGDOM" \
    --abundance_file "$ABUNDANCE_FILE" \
    --model "$MODEL" \
    --n_reads "$N_READS" \
    --output "$PREFIX" \
    --cpus "$CPUS"

# ISS may write .fastq or .fastq.gz
if [[ -f "${PREFIX}_R1.fastq.gz" ]]; then
    gunzip "${PREFIX}_R1.fastq.gz"
    gunzip "${PREFIX}_R2.fastq.gz"
fi

# ---------------------------------------------------------------------------
# Step 2: megahit_topo — build SDBG (memory capped), extract graph features
# Runs after ISS finishes so full 16 GB is available.
# ---------------------------------------------------------------------------
echo "[2/3] megahit_topo  (mem=${TOPO_MEM_GB}G)"
$MEGAHIT_TOPO \
    --reads  "$R1" \
    --reads2 "$R2" \
    --output "$TOPO_JSON" \
    --threads "$CPUS" \
    --mem "$TOPO_MEM_GB"

# ---------------------------------------------------------------------------
# Step 3: jellyfish — count true vs error k-mers
# megahit_topo is done and its memory is freed before jellyfish starts.
# Reference FASTA is the one ISS downloaded — exact ground truth.
# ---------------------------------------------------------------------------
echo "[3/3] jellyfish kmer counting"

$JELLYFISH count -m 21 -s 200M -t "$CPUS" -C -o "$REF_JF"   "$FASTA"
$JELLYFISH count -m 21 -s 100M -t "$CPUS" -C -o "$READS_JF" "$R1" "$R2"

$JELLYFISH dump "$READS_JF" | awk '{print ">"NR"\n"$1}' > "$READS_KMERS_FA"

N_TRUE=$($JELLYFISH query "$REF_JF" -s "$READS_KMERS_FA" | awk '$2>0{t++} END{print t+0}')
N_ERR=$($JELLYFISH  query "$REF_JF" -s "$READS_KMERS_FA" | awk '$2==0{e++} END{print e+0}')
N_TOTAL=$((N_TRUE + N_ERR))
ERR_FRAC=$(awk "BEGIN{printf \"%.6f\", ($N_TOTAL>0 ? $N_ERR/$N_TOTAL : 0)}")

# ---------------------------------------------------------------------------
# Merge topo JSON + kmer truth + metadata → final output
# ---------------------------------------------------------------------------
cp "$ABUNDANCE_FILE" "$OUTPUT_DIR/${RUN_ID}_abundance.txt"

python3 - <<PYEOF
import json
with open("$TOPO_JSON") as f:
    topo = json.load(f)
topo["run_id"]            = "$RUN_ID"
topo["n_reads_requested"] = $N_READS
topo["model"]             = "$MODEL"
topo["abundance"]         = "${RUN_ID}_abundance.txt"
topo["kmers"] = {
    "n_true":         $N_TRUE,
    "n_error":        $N_ERR,
    "error_fraction": $ERR_FRAC,
}
with open("$FINAL_JSON", "w") as f:
    json.dump(topo, f, indent=2)
print(f"  true k-mers:    {$N_TRUE:,}")
print(f"  error k-mers:   {$N_ERR:,}")
print(f"  error fraction: {float('$ERR_FRAC'):.4%}")
PYEOF

# ---------------------------------------------------------------------------
# Cleanup — reads, FASTA, jellyfish dbs all deleted
# Only JSON + abundance survive in output_dir
# ---------------------------------------------------------------------------
if [[ $KEEP_READS -eq 0 ]]; then
    echo "Cleaning up..."
    rm -f "$R1" "$R2" "$FASTA" "$REF_JF" "$READS_JF" "$READS_KMERS_FA"
    rm -f "$WORK"/*.vcf "$WORK"/*.bam "$WORK"/*.bai
fi

echo "Done -> $FINAL_JSON"
