#!/usr/bin/env python3
"""
Generate a synthetic BAM of perfect reads aligned to a local reference FASTA,
with an optional sharply-localised coverage spike.

Purpose: a *smoke test* for panfreebayes core — exercises the realignment +
haplotype-window + genotyping code path without needing the real (5GB+) BAMs.
It is NOT a correctness test: perfect reads should yield ~no variant calls
(a handful from read-end / homopolymer effects is normal).

The spike option targets Milestone 1's stated risk: buildHaplotypeAlleles
window-growth behaviour at a position with coverage far above the regional
average (the real DL238 bubble has a ~23,000x spike over a ~900x baseline).

Requires: samtools on PATH (add `samtools` to the guix shell).

Usage:
  make_synthetic_bam.py --ref REF.fasta --out reads.bam \
      [--depth 30] [--readlen 150] [--sample DL238] \
      [--spike-at 0.5] [--spike-depth 23000] [--spike-width 1] \
      [--seed 1]

  --spike-at      fraction (0..1) OR absolute 0-based coord of the spike centre
                  in the FIRST contig; omit for no spike
  --spike-depth   total reads stacked over the spike window
  --spike-width   width in bp of the spike window (reads fully span it)
"""
import argparse
import os
import subprocess
import sys
import random


def read_fasta(path):
    seqs = []
    name = None
    buf = []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if name is not None:
                    seqs.append((name, "".join(buf)))
                name = line[1:].split()[0]
                buf = []
            else:
                buf.append(line)
    if name is not None:
        seqs.append((name, "".join(buf)))
    return seqs


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ref", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--depth", type=int, default=30)
    ap.add_argument("--readlen", type=int, default=150)
    ap.add_argument("--sample", default="SYN")
    ap.add_argument("--spike-at", default=None)
    ap.add_argument("--spike-depth", type=int, default=23000)
    ap.add_argument("--spike-width", type=int, default=1)
    ap.add_argument("--seed", type=int, default=1)
    args = ap.parse_args()

    random.seed(args.seed)
    seqs = read_fasta(args.ref)
    if not seqs:
        sys.exit("no sequences in " + args.ref)

    rg = "syn1"
    sam_path = args.out + ".sam"
    with open(sam_path, "w") as sam:
        sam.write("@HD\tVN:1.6\tSO:coordinate\n")
        for name, seq in seqs:
            sam.write("@SQ\tSN:%s\tLN:%d\n" % (name, len(seq)))
        sam.write("@RG\tID:%s\tSM:%s\tPL:PACBIO\n" % (rg, args.sample))

        read_no = 0

        def emit(chrom, pos0, s):
            nonlocal read_no
            read_no += 1
            qual = "I" * len(s)
            # 0-based pos0 -> 1-based SAM POS
            sam.write("\t".join([
                "syn_%d" % read_no, "0", chrom, str(pos0 + 1), "60",
                "%dM" % len(s), "*", "0", "0", s, qual, "RG:Z:" + rg,
            ]) + "\n")

        for name, seq in seqs:
            L = len(seq)
            if L < args.readlen:
                continue
            step = max(1, args.readlen // max(1, args.depth))
            for start in range(0, L - args.readlen + 1, step):
                emit(name, start, seq[start:start + args.readlen])

        # coverage spike in the first contig
        if args.spike_at is not None:
            name, seq = seqs[0]
            L = len(seq)
            sa = args.spike_at
            centre = int(float(sa) * L) if ("." in str(sa) or 0.0 <= float(sa) <= 1.0 and float(sa) != int(float(sa))) else int(sa)
            centre = max(args.readlen, min(L - args.readlen - 1, centre))
            win_lo = centre
            win_hi = centre + max(1, args.spike_width)
            rstart = max(0, win_hi - args.readlen)
            rstart = min(rstart, win_lo)  # ensure read spans the whole window
            for _ in range(args.spike_depth):
                jitter = random.randint(-2, 2)
                st = max(0, min(L - args.readlen, rstart + jitter))
                emit(name, st, seq[st:st + args.readlen])
            sys.stderr.write(
                "spike: %s ~%d reads centred at 0-based %d (window %d..%d)\n"
                % (name, args.spike_depth, centre, win_lo, win_hi))

    subprocess.check_call(["samtools", "sort", "-o", args.out, sam_path])
    subprocess.check_call(["samtools", "index", args.out])
    os.remove(sam_path)
    if not os.path.exists(args.ref + ".fai"):
        subprocess.check_call(["samtools", "faidx", args.ref])
    sys.stderr.write("wrote %s (%d reads)\n" % (args.out, read_no))


if __name__ == "__main__":
    main()
