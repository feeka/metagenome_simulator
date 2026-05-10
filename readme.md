# metagenome-simulator

Simulates Illumina paired-end reads from a synthetic microbial community
using [InSilicoSeq](https://github.com/HadrienG/InSilicoSeq).

## What it does

1. Picks a set of organisms and assigns them relative abundances
2. Runs ISS with a HiSeq / MiSeq / NovaSeq error profile to simulate paired-end reads
3. Optionally estimates the actual error rate ISS injected using jellyfish k-mer analysis

## Requirements

| Package | Install |
|---|---|
| InSilicoSeq | `pip install InSilicoSeq` |
| numpy | `pip install numpy` |
| pandas + biopython | `pip install pandas biopython` (local mode only) |
| jellyfish 2 | `apt/brew/conda install jellyfish` (only with `--jellyfish`) |

## Usage

```bash
# one run, any model, NCBI genomes
python run.py --output sims/ --ncbi-mode

# 5 runs, MiSeq, with error rate estimation
python run.py --output sims/ --n-runs 5 --model MiSeq --jellyfish --ncbi-mode

# 10 runs, resume if interrupted
python run.py --output sims/ --n-runs 10 --ncbi-mode --resume
```

## Flags

| Flag | Default | Description |
|---|---|---|
| `--output DIR` | (required) | Directory for output files |
| `--n-runs N` | `1` | Number of simulation runs |
| `--model` | random | `HiSeq` \| `MiSeq` \| `NovaSeq` |
| `--n-reads MIN MAX` | `300 10000` | Read count range in thousands (log-uniform) |
| `--community-size MIN MAX` | `3 200` | Genomes per community, log-uniform (local mode only) |
| `--ncbi-mode` | off | Download genomes from NCBI each run (no setup required) |
| `--jellyfish` | off | Estimate ISS error rate via k-mer analysis |
| `--seed N` | `42` | Random seed |
| `--cpus N` | all | CPUs per run |
| `--resume` | off | Skip runs whose reads already exist |

## Modes

### NCBI (recommended)

Pass `--ncbi-mode`. ISS downloads genomes from NCBI each run — no setup needed.
Randomly picks 1–3 kingdoms from `bacteria`, `viruses`, `archaea`.

### LOCAL

Default when a `references/manifest.tsv` file is present.
Expects pre-downloaded FASTA files organised by kingdom.
Community composition varies by biome (`gut`, `soil`, `marine`, `clinical`, `oral`).

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
