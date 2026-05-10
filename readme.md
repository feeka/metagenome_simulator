# metagenome-simulator

Simulates Illumina paired-end reads from a synthetic microbial community
using [InSilicoSeq](https://github.com/HadrienG/InSilicoSeq).

## What it does

1. Selects a community of organisms (biome-weighted or random NCBI download)
2. Draws relative abundances from a log-normal distribution
3. Runs ISS with a HiSeq / MiSeq / NovaSeq error profile to simulate reads
4. Optionally: runs jellyfish to count how many unique k-mers in the reads
   are absent from the reference genomes — this is the ISS error injection rate

## Requirements

| Package | Install |
|---|---|
| InSilicoSeq | `pip install InSilicoSeq` |
| numpy | `pip install numpy` |
| pandas | `pip install pandas` (local mode only) |
| biopython | `pip install biopython` (local mode only) |
| jellyfish 2 | `apt/brew/conda install jellyfish` (only with `--jellyfish`) |

## Usage

```bash
# single run, gut biome, random model, local genomes
python run.py --output sims/

# 5 runs, MiSeq error model, with error rate estimation
python run.py --output sims/ --n-runs 5 --model MiSeq --jellyfish

# NCBI mode -- no local genomes needed
python run.py --output sims/ --ncbi-mode

# resume a partial batch, 10 runs
python run.py --output sims/ --n-runs 10 --resume
```

## Flags

| Flag | Default | Description |
|---|---|---|
| `--output DIR` | (required) | Directory for output files |
| `--n-runs N` | `1` | Number of simulation runs |
| `--model` | random | `HiSeq` \| `MiSeq` \| `NovaSeq` |
| `--n-reads MIN MAX` | `300 10000` | Read count range in thousands (log-uniform) |
| `--community-size MIN MAX` | `3 200` | Genomes per community (log-uniform, local mode) |
| `--ncbi-mode` | off | Download genomes from NCBI instead of using local store |
| `--jellyfish` | off | Estimate ISS error rate via k-mer analysis |
| `--seed N` | `42` | Random seed |
| `--cpus N` | all | CPUs per run |
| `--resume` | off | Skip runs whose reads already exist |

## Modes

### LOCAL (default)

Reads genomes from `references/manifest.tsv`.
Run `01_download_references.py` first to populate the local store.
Community composition is sampled from biome-weighted Dirichlet distributions:

| Biome | Dominant kingdoms |
|---|---|
| gut | bacteria 72%, archaea 8%, host 7% |
| soil | bacteria 58%, fungi 15%, protists 8% |
| marine | bacteria 60%, archaea 18%, viral 10% |
| clinical | bacteria 50%, host 35% |
| oral | bacteria 74%, fungi 8%, viral 8% |

### NCBI

ISS fetches genomes directly from NCBI each run.
Randomly picks 1–3 kingdoms from: `bacteria`, `viruses`, `archaea`.
Auto-enabled when `references/manifest.tsv` is missing.

## Output per run

| File | Contents |
|---|---|
| `<run_id>_R1.fastq` | Read 1 (paired-end) |
| `<run_id>_R2.fastq` | Read 2 (paired-end) |
| `<run_id>_abundance.txt` | Genome name -> relative fraction |
| `<run_id>_kmer_stats.txt` | k-mer error stats (`--jellyfish` only) |

## Jellyfish error estimation (`--jellyfish`)

Runs automatically when `--jellyfish` is passed.
When running interactively without the flag, the tool asks before starting.

**What it measures:**

- Counts all unique canonical 21-mers in the simulated reads
- Queries which of those k-mers are absent from the reference genomes
- Missing k-mers = sequencing-error artifacts from the ISS error model

**Output (`<run_id>_kmer_stats.txt`):**

```
total_kmers   1 234 567
error_kmers   89 012
error_rate    7.2124%

error_rate = unique 21-mers in reads NOT found in reference.
One sequencing error corrupts up to 21 adjacent k-mers,
so error_rate > per-base error rate.
```

**Note:** because one substitution error corrupts up to k=21 overlapping k-mers,
`error_rate` is larger than the ISS per-base error rate.
A rough per-base estimate: `error_rate x mean_error_multiplicity / mean_depth`.