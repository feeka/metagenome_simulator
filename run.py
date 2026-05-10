#!/usr/bin/env python3
"""
run.py — metagenome read simulator

Generates Illumina paired-end reads from a synthetic microbial community
using InSilicoSeq (iss).

Modes:
  LOCAL   Genomes already on disk (run 01_download_references.py first).
          Community composition is drawn from biome-weighted Dirichlet
          mixtures — proportions vary realistically run to run.
  NCBI    ISS downloads genomes from NCBI each run (--ncbi-mode).
          No local genome store required; works immediately.

Output per run:
  <run_id>_R1.fastq        read 1 (paired-end)
  <run_id>_R2.fastq        read 2 (paired-end)
  <run_id>_abundance.txt   genome -> relative fraction table

Optional (--jellyfish):
  Counts canonical 21-mers in the simulated reads and in the reference
  genomes, then queries how many unique read k-mers are absent from the
  reference.  These missing k-mers are sequencing-error artifacts
  introduced by the ISS error model.

  Saves <run_id>_kmer_stats.txt and prints a summary:
    total_kmers   unique canonical k-mers in reads
    error_kmers   k-mers absent from reference
    error_rate    error_kmers / total_kmers

  Note: one base error corrupts up to k=21 overlapping k-mers, so
  error_rate exceeds the per-base error rate.  Coverage-weighted
  per-base estimate ~= error_rate x mean_error_multiplicity / mean_depth.

Requirements:
  iss        pip install InSilicoSeq
  numpy      pip install numpy
  jellyfish  apt/brew/conda install jellyfish   (only with --jellyfish)
  pandas     pip install pandas                 (local mode only)
  biopython  pip install biopython              (local mode only)

Usage:
  python run.py --output sims/
  python run.py --output sims/ --n-runs 5 --model MiSeq --jellyfish
  python run.py --output sims/ --ncbi-mode --n-reads 500 2000
  python run.py --output sims/ --n-runs 10 --seed 7 --resume
"""

import argparse
import gzip
import logging
import os
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np

try:
    import pandas as pd
    from Bio import SeqIO
    _LOCAL_DEPS = True
except ImportError:
    _LOCAL_DEPS = False

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

BASE_DIR = Path(__file__).parent
MANIFEST = BASE_DIR / "references" / "manifest.tsv"

ISS       = "iss"
JELLYFISH = "jellyfish"

MODELS = ["HiSeq", "MiSeq", "NovaSeq"]

# Controls Dirichlet spread: lower = more run-to-run variance in community
# composition; higher = proportions stay closer to the biome means.
ALPHA_STRENGTH = 5.0

# Biome-mean kingdom fractions used as Dirichlet concentration parameters.
# Sources: rough order-of-magnitude estimates from HMP, Tara Oceans, etc.
BIOME_MIXES = {
    "gut": {
        "bacteria": 0.72, "archaea": 0.08, "fungi": 0.04,
        "viral": 0.06,    "protists": 0.01, "host": 0.07, "plasmids": 0.02,
    },
    "soil": {
        "bacteria": 0.58, "archaea": 0.07, "fungi": 0.15,
        "viral": 0.05,    "protists": 0.08, "host": 0.02, "plasmids": 0.05,
    },
    "marine": {
        "bacteria": 0.60, "archaea": 0.18, "fungi": 0.02,
        "viral": 0.10,    "protists": 0.08, "host": 0.01, "plasmids": 0.01,
    },
    "clinical": {
        "bacteria": 0.50, "archaea": 0.02, "fungi": 0.05,
        "viral": 0.05,    "protists": 0.01, "host": 0.35, "plasmids": 0.02,
    },
    "oral": {
        "bacteria": 0.74, "archaea": 0.03, "fungi": 0.08,
        "viral": 0.08,    "protists": 0.01, "host": 0.04, "plasmids": 0.02,
    },
}

NCBI_KINGDOMS = ["bacteria", "viruses", "archaea"]

# ---------------------------------------------------------------------------
# Logging
# ---------------------------------------------------------------------------

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s  %(levelname)-7s  %(message)s",
    handlers=[logging.StreamHandler(sys.stdout)],
)
log = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Genome selection -- LOCAL mode
# ---------------------------------------------------------------------------

def load_manifest(path: Path):
    if not path.exists():
        return None
    if not _LOCAL_DEPS:
        log.warning("pandas/biopython not installed -- falling back to --ncbi-mode")
        return None
    df = pd.read_csv(path, sep="\t")
    required = {"path", "kingdom"}
    missing = required - set(df.columns)
    if missing:
        log.error(f"manifest.tsv missing columns: {missing}")
        sys.exit(1)
    df = df[df["path"].apply(lambda p: Path(p).exists())]
    log.info(f"Manifest: {len(df)} genomes across {df['kingdom'].nunique()} kingdoms")
    return df


def _biome_fractions(biome: str, rng: np.random.Generator) -> dict:
    mix = BIOME_MIXES[biome]
    kingdoms = list(mix.keys())
    alpha = np.array([max(mix[k], 0.01) * ALPHA_STRENGTH for k in kingdoms])
    fracs = rng.dirichlet(alpha)
    return dict(zip(kingdoms, fracs))


def select_genomes(community_size, biome, manifest, rng):
    fracs = _biome_fractions(biome, rng)
    kingdoms = list(fracs.keys())

    counts = {}
    remaining = community_size
    for i, k in enumerate(kingdoms):
        if i == len(kingdoms) - 1:
            counts[k] = max(0, remaining)
        else:
            counts[k] = max(0, round(community_size * fracs[k]))
            remaining -= counts[k]

    paths = []
    for k, n in counts.items():
        if n == 0:
            continue
        available = manifest[manifest["kingdom"] == k]["path"].tolist()
        if not available:
            available = manifest[manifest["kingdom"] == "bacteria"]["path"].tolist()
        n = min(n, len(available))
        if n == 0:
            continue
        chosen = rng.choice(available, size=n, replace=False).tolist()
        paths.extend(Path(p) for p in chosen)

    bac = manifest[manifest["kingdom"] == "bacteria"]["path"].tolist()
    while len(paths) < community_size and bac:
        paths.append(Path(rng.choice(bac)))

    paths = paths[:community_size]
    rng.shuffle(paths)

    raw = rng.lognormal(mean=0.0, sigma=1.5, size=len(paths))
    abundances = (raw / raw.sum()).tolist()
    return paths, abundances


def get_first_record_id(fasta_path):
    with open(fasta_path) as fh:
        for rec in SeqIO.parse(fh, "fasta"):
            return rec.id
    raise ValueError(f"Empty FASTA: {fasta_path}")


def write_iss_abundance(genome_paths, abundances, out_path):
    with open(out_path, "w") as fh:
        for gp, ab in zip(genome_paths, abundances):
            seq_id = get_first_record_id(gp)
            fh.write(f"{seq_id}\t{ab:.8f}\n")

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _run(cmd, desc):
    log.info(f"  {desc}")
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        log.error(f"FAILED: {' '.join(cmd)}")
        log.error(result.stderr[-2000:])
        raise RuntimeError(desc)


def _gunzip_if_needed(work, run_id):
    for suf in ["_R1.fastq.gz", "_R2.fastq.gz"]:
        gz = work / f"{run_id}{suf}"
        if gz.exists():
            fq = work / f"{run_id}{suf[:-3]}"
            with gzip.open(gz, "rb") as gi, open(fq, "wb") as fo:
                shutil.copyfileobj(gi, fo)
            gz.unlink()

# ---------------------------------------------------------------------------
# Jellyfish k-mer error analysis
# ---------------------------------------------------------------------------

def run_jellyfish_analysis(r1, r2, ref_fa, work, cpus):
    """
    Count canonical 21-mers in the reference and in the reads.
    Query which read k-mers are absent from the reference -- these are
    error-corrupted k-mers injected by the ISS error model.

    Returns dict:
      total_kmers   unique canonical k-mers found in reads
      error_kmers   k-mers absent from reference
      error_rate    error_kmers / total_kmers
    """
    ref_jf     = work / "ref.jf"
    reads_jf   = work / "reads.jf"
    reads_dump = work / "reads_dump.fa"

    _run(
        [JELLYFISH, "count", "-C", "-m", "21", "-s", "200M",
         "-t", str(cpus), "-o", str(ref_jf), str(ref_fa)],
        "jellyfish count reference",
    )
    _run(
        [JELLYFISH, "count", "-C", "-m", "21", "-s", "100M",
         "-t", str(cpus), "-o", str(reads_jf), str(r1), str(r2)],
        "jellyfish count reads",
    )

    dump = subprocess.run(
        [JELLYFISH, "dump", str(reads_jf)],
        capture_output=True, text=True, check=True,
    )
    with open(reads_dump, "w") as fh:
        idx = 0
        for line in dump.stdout.splitlines():
            line = line.strip()
            if not line.startswith(">"):
                fh.write(f">k{idx}\n{line}\n")
                idx += 1

    query = subprocess.run(
        [JELLYFISH, "query", str(ref_jf), "-s", str(reads_dump)],
        capture_output=True, text=True, check=True,
    )
    n_true = n_err = 0
    for line in query.stdout.splitlines():
        parts = line.split()
        if len(parts) == 2:
            if int(parts[1]) > 0:
                n_true += 1
            else:
                n_err += 1

    total = n_true + n_err
    rate  = n_err / total if total else 0.0
    return {
        "total_kmers": total,
        "error_kmers": n_err,
        "error_rate":  round(rate, 6),
    }


def _save_kmer_stats(stats, path):
    with open(path, "w") as fh:
        fh.write(f"total_kmers   {stats['total_kmers']:,}\n")
        fh.write(f"error_kmers   {stats['error_kmers']:,}\n")
        fh.write(f"error_rate    {stats['error_rate']:.4%}\n")
        fh.write("\n")
        fh.write("error_rate = unique 21-mers in reads NOT found in reference.\n")
        fh.write("One sequencing error corrupts up to 21 adjacent k-mers,\n")
        fh.write("so error_rate > per-base error rate.\n")

# ---------------------------------------------------------------------------
# LOCAL mode -- one run
# ---------------------------------------------------------------------------

def run_one(run_id, genome_paths, abundances, n_reads, model, biome,
            output_dir, cpus, use_jellyfish):
    work = output_dir / ".work" / run_id
    work.mkdir(parents=True, exist_ok=True)
    try:
        prefix  = str(work / run_id)
        r1      = work / f"{run_id}_R1.fastq"
        r2      = work / f"{run_id}_R2.fastq"
        ref_fa  = work / "reference.fasta"
        ab_file = work / "abundance.txt"

        write_iss_abundance(genome_paths, abundances, ab_file)

        log.info(f"  Concatenating {len(genome_paths)} genomes -> reference FASTA")
        with open(ref_fa, "wb") as out:
            for gp in genome_paths:
                with open(gp, "rb") as inp:
                    shutil.copyfileobj(inp, out)

        _run(
            [
                ISS, "generate",
                "--genomes",        *[str(p) for p in genome_paths],
                "--abundance_file",  str(ab_file),
                "--model",           model,
                "--n_reads",         str(n_reads),
                "--output",          prefix,
                "--cpus",            str(cpus),
            ],
            f"iss generate ({n_reads:,} reads, model={model})",
        )
        _gunzip_if_needed(work, run_id)

        out_r1 = output_dir / f"{run_id}_R1.fastq"
        out_r2 = output_dir / f"{run_id}_R2.fastq"
        shutil.move(str(r1), out_r1)
        shutil.move(str(r2), out_r2)

        out_ab = output_dir / f"{run_id}_abundance.txt"
        with open(out_ab, "w") as fh:
            fh.write(f"# biome: {biome}\n")
            for gp, ab in zip(genome_paths, abundances):
                fh.write(f"{gp.name}\t{ab:.8f}\n")

        if use_jellyfish:
            stats = run_jellyfish_analysis(out_r1, out_r2, ref_fa, work, cpus)
            kmer_out = output_dir / f"{run_id}_kmer_stats.txt"
            _save_kmer_stats(stats, kmer_out)
            log.info(
                f"  k-mer stats -> {kmer_out.name}  "
                f"(total: {stats['total_kmers']:,}  "
                f"error: {stats['error_kmers']:,}  "
                f"rate: {stats['error_rate']:.2%})"
            )
    finally:
        shutil.rmtree(work, ignore_errors=True)

# ---------------------------------------------------------------------------
# NCBI mode -- one run
# ---------------------------------------------------------------------------

def run_one_ncbi(run_id, n_reads, model, output_dir, cpus, rng, use_jellyfish):
    n_kingdoms    = int(rng.integers(1, 4))
    kingdoms      = list(rng.choice(NCBI_KINGDOMS, size=n_kingdoms, replace=False))
    total_genomes = max(n_kingdoms, int(np.exp(rng.uniform(np.log(2), np.log(20)))))

    base   = total_genomes // n_kingdoms
    counts = [base] * n_kingdoms
    counts[-1] += total_genomes - sum(counts)
    counts = [max(1, c) for c in counts]

    abundance_dist = str(rng.choice(["uniform", "halfnormal", "exponential"]))

    work = output_dir / ".work" / run_id
    work.mkdir(parents=True, exist_ok=True)
    try:
        prefix = str(work / run_id)
        r1     = work / f"{run_id}_R1.fastq"
        r2     = work / f"{run_id}_R2.fastq"
        ref_fa = work / f"{run_id}_ncbi_genomes.fasta"

        _run(
            [
                ISS, "generate",
                "--ncbi",            *kingdoms,
                "--n_genomes_ncbi",  *[str(n) for n in counts],
                "--abundance",        abundance_dist,
                "--model",            model,
                "--n_reads",          str(n_reads),
                "--output",           prefix,
                "--cpus",             str(cpus),
            ],
            f"iss generate --ncbi {'+'.join(kingdoms)} ({n_reads:,} reads, model={model})",
        )
        _gunzip_if_needed(work, run_id)

        out_r1 = output_dir / f"{run_id}_R1.fastq"
        out_r2 = output_dir / f"{run_id}_R2.fastq"
        shutil.move(str(r1), out_r1)
        shutil.move(str(r2), out_r2)

        out_ab = output_dir / f"{run_id}_abundance.txt"
        with open(out_ab, "w") as fh:
            fh.write("# NCBI mode -- kingdom composition\n")
            for k, n in zip(kingdoms, counts):
                fh.write(f"{k}\t{n} genomes\n")

        if use_jellyfish:
            if ref_fa.exists():
                stats = run_jellyfish_analysis(out_r1, out_r2, ref_fa, work, cpus)
                kmer_out = output_dir / f"{run_id}_kmer_stats.txt"
                _save_kmer_stats(stats, kmer_out)
                log.info(
                    f"  k-mer stats -> {kmer_out.name}  "
                    f"(total: {stats['total_kmers']:,}  "
                    f"error: {stats['error_kmers']:,}  "
                    f"rate: {stats['error_rate']:.2%})"
                )
            else:
                log.warning("  jellyfish skipped -- ISS did not write a reference FASTA")
    finally:
        shutil.rmtree(work, ignore_errors=True)

# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def _ask_jellyfish():
    """Interactive y/N prompt -- only called when stdin is a terminal."""
    try:
        ans = input("Run jellyfish k-mer error analysis on simulated reads? [y/N] ").strip().lower()
        print()
        return ans in ("y", "yes")
    except (EOFError, KeyboardInterrupt):
        return False


def parse_args():
    p = argparse.ArgumentParser(
        description=(
            "Metagenome read simulator.\n"
            "Generates Illumina paired-end reads from a synthetic microbial\n"
            "community using InSilicoSeq (iss)."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
examples:
  python run.py --output sims/
  python run.py --output sims/ --n-runs 5 --model MiSeq --jellyfish
  python run.py --output sims/ --ncbi-mode --n-reads 500 2000
  python run.py --output sims/ --n-runs 10 --seed 7 --resume
""",
    )
    p.add_argument("--output",    required=True, type=Path,
                   help="Output directory for reads and abundance files")
    p.add_argument("--n-runs",    type=int, default=1,
                   help="Number of simulation runs (default: 1)")
    p.add_argument("--seed",      type=int, default=42,
                   help="Random seed (default: 42)")
    p.add_argument("--cpus",      type=int, default=os.cpu_count() or 1,
                   help="CPUs per run (default: all cores)")
    p.add_argument("--model",     choices=MODELS, default=None,
                   help="Sequencing error model: HiSeq | MiSeq | NovaSeq "
                        "(default: random per run)")
    p.add_argument("--resume",    action="store_true",
                   help="Skip runs whose reads already exist in --output")
    p.add_argument("--ncbi-mode", action="store_true",
                   help="Download genomes from NCBI each run -- no local "
                        "genome store required; auto-enabled when no manifest.tsv found")
    p.add_argument("--jellyfish", action="store_true",
                   help="Estimate ISS error rate with jellyfish k-mer analysis "
                        "(jellyfish must be on PATH)")
    p.add_argument(
        "--community-size", nargs=2, type=float, metavar=("MIN", "MAX"),
        default=[3, 200],
        help="[local mode] Community size log-uniform range (default: 3 200)",
    )
    p.add_argument(
        "--n-reads", nargs=2, type=float, metavar=("MIN_K", "MAX_K"),
        default=[300, 10000],
        help="Read count range in thousands, log-uniform (default: 300 10000)",
    )
    return p.parse_args()


def main():
    args = parse_args()
    args.output.mkdir(parents=True, exist_ok=True)

    rng = np.random.default_rng(args.seed)

    manifest = load_manifest(MANIFEST)
    use_ncbi = args.ncbi_mode or (manifest is None)

    if use_ncbi:
        log.info("Mode: NCBI  (ISS downloads genomes per run)")
    else:
        log.info(f"Mode: LOCAL  ({len(manifest)} genomes in manifest)")

    # Jellyfish: use --jellyfish flag, or ask interactively in a terminal
    use_jellyfish = args.jellyfish
    if not use_jellyfish and sys.stdin.isatty():
        use_jellyfish = _ask_jellyfish()

    if use_jellyfish:
        if shutil.which(JELLYFISH) is None:
            log.error("jellyfish not found on PATH -- install it or omit --jellyfish")
            sys.exit(1)
        log.info("jellyfish: enabled  (k=21, canonical k-mer error analysis)")
    else:
        log.info("jellyfish: disabled")

    log.info(f"Starting: {args.n_runs} run(s) -> {args.output}")

    biomes = list(BIOME_MIXES.keys())
    succeeded = failed = skipped = 0

    for i in range(args.n_runs):
        run_id = f"run_{i+1:04d}"
        out_r1 = args.output / f"{run_id}_R1.fastq"

        if args.resume and out_r1.exists():
            log.info(f"[{i+1}/{args.n_runs}] {run_id} -- skipping (already done)")
            skipped += 1
            continue

        n_reads = max(10_000, int(
            np.exp(rng.uniform(
                np.log(args.n_reads[0] * 1000),
                np.log(args.n_reads[1] * 1000),
            ))
        ))
        model = args.model if args.model else str(rng.choice(MODELS))

        try:
            if use_ncbi:
                log.info(f"[{i+1}/{args.n_runs}] {run_id}  reads={n_reads:,}  model={model}")
                run_one_ncbi(
                    run_id=run_id, n_reads=n_reads, model=model,
                    output_dir=args.output, cpus=args.cpus, rng=rng,
                    use_jellyfish=use_jellyfish,
                )
            else:
                biome = str(rng.choice(biomes))
                community_size = max(2, int(np.exp(rng.uniform(
                    np.log(args.community_size[0]),
                    np.log(args.community_size[1]),
                ))))
                genome_paths, abundances = select_genomes(
                    community_size, biome, manifest, rng
                )
                log.info(
                    f"[{i+1}/{args.n_runs}] {run_id}  biome={biome}  "
                    f"size={len(genome_paths)}  reads={n_reads:,}  model={model}"
                )
                run_one(
                    run_id=run_id, genome_paths=genome_paths, abundances=abundances,
                    n_reads=n_reads, model=model, biome=biome,
                    output_dir=args.output, cpus=args.cpus,
                    use_jellyfish=use_jellyfish,
                )
            succeeded += 1
        except Exception as exc:
            log.error(f"  FAILED: {exc}")
            failed += 1

    log.info(f"Done: {succeeded} succeeded, {failed} failed, {skipped} skipped")


if __name__ == "__main__":
    main()