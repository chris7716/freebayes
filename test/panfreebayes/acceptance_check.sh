#!/usr/bin/env bash
#
# panFreebayes acceptance / equivalence check (Milestone 2 harness core).
#
# Runs panfreebayes and stock freebayes on the SAME reference + BAM with the
# SAME parameters, and checks:
#   - both produce the same number of variant calls
#   - their VCF bodies are byte-identical (ignoring ## header lines)
#   - optionally, that the call count equals an expected value (--expect)
#
# It does NOT know where your data lives — pass --ref and --bam. Nothing is
# hard-coded to a machine layout, so this is the script to wire into CI /
# a permanent regression test once real data is available.
#
# Usage:
#   test/panfreebayes/acceptance_check.sh \
#       --ref  regions/chrII/DL238_..._1664999-1748589.fasta \
#       --bam  align/chrII/DL238_DL238_..._1664999-1748589.bam \
#       --expect 1413 \
#       [--outdir .] [--freebayes /path/to/freebayes] [--panfreebayes /path/to/panfreebayes] \
#       [-- <extra args passed verbatim to BOTH tools>]
#
# Default calling parameters (the empirically-needed long-read set) are:
#   --pooled-continuous --min-alternate-count 2 --min-alternate-fraction 0.2 --limit-coverage 200
# Pass `-- <args>` to override/extend them (replaces the defaults entirely).
#
# Must be run inside the guix build environment so the binaries find their libs:
#   guix shell --pure bash grep coreutils diffutils gcc-toolchain \
#     zlib xz bzip2 htslib simde tabixpp fastahack smithwaterman intervaltree \
#     -- bash test/panfreebayes/acceptance_check.sh --ref ... --bam ... --expect 1413
#
set -euo pipefail

here=$(cd "$(dirname "$0")" && pwd)
build_default=$(cd "$here/../.." && pwd)/build

REF="" BAM="" EXPECT="" OUTDIR="."
FB="$build_default/freebayes"
PFB="$build_default/panfreebayes"
EXTRA=()

while [ $# -gt 0 ]; do
  case "$1" in
    --ref)          REF=$2; shift 2 ;;
    --bam)          BAM=$2; shift 2 ;;
    --expect)       EXPECT=$2; shift 2 ;;
    --outdir)       OUTDIR=$2; shift 2 ;;
    --freebayes)    FB=$2; shift 2 ;;
    --panfreebayes) PFB=$2; shift 2 ;;
    --)             shift; EXTRA=("$@"); break ;;
    -h|--help)      sed -n '2,40p' "$0"; exit 0 ;;
    *) echo "unknown arg: $1" >&2; exit 2 ;;
  esac
done

[ -n "$REF" ] && [ -n "$BAM" ] || { echo "need --ref and --bam" >&2; exit 2; }
[ -x "$FB" ]  || { echo "freebayes not executable: $FB" >&2; exit 2; }
[ -x "$PFB" ] || { echo "panfreebayes not executable: $PFB" >&2; exit 2; }
[ -f "$REF" ] || { echo "no such ref: $REF" >&2; exit 2; }
[ -f "$BAM" ] || { echo "no such bam: $BAM" >&2; exit 2; }
mkdir -p "$OUTDIR"

if [ "${#EXTRA[@]}" -eq 0 ]; then
  EXTRA=(--pooled-continuous --min-alternate-count 2 --min-alternate-fraction 0.2 --limit-coverage 200)
fi

base=$(basename "$BAM" .bam)
pfb_vcf="$OUTDIR/${base}.panfreebayes.vcf"
fb_vcf="$OUTDIR/${base}.freebayes.vcf"
diff_out="$OUTDIR/${base}.pfb_vs_fb.diff"

echo "ref:   $REF"
echo "bam:   $BAM"
echo "args:  ${EXTRA[*]}"
echo

echo "running panfreebayes ..."
time "$PFB" --ref "$REF" --bam "$BAM" "${EXTRA[@]}" > "$pfb_vcf"
echo "running freebayes ..."
time "$FB" -f "$REF" "$BAM" "${EXTRA[@]}" > "$fb_vcf"
echo

pfb_n=$(grep -vc '^#' "$pfb_vcf" || true)
fb_n=$(grep -vc  '^#' "$fb_vcf"  || true)
echo "panfreebayes calls: $pfb_n"
echo "freebayes    calls: $fb_n"

rc=0

if diff <(grep -v '^##' "$pfb_vcf") <(grep -v '^##' "$fb_vcf") > "$diff_out"; then
  echo "VCF body vs freebayes: IDENTICAL"
else
  echo "VCF body vs freebayes: DIFFERS  ($(grep -c '^[<>]' "$diff_out") changed lines -> $diff_out)"
  rc=1
fi

if [ -n "$EXPECT" ]; then
  if [ "$pfb_n" = "$EXPECT" ]; then
    echo "expected count ($EXPECT): MATCH"
  else
    echo "expected count ($EXPECT): MISMATCH (got $pfb_n)"
    rc=1
  fi
fi

exit $rc
