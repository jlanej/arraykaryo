"""arraykaryo inspect: what an array export carries, before anything is computed - the samples, the contigs and the assembly
they imply, the records per contig (from the index), the FORMAT and INFO fields, and the conventions seen in a sample of
records (how a no-call, a missing LRR or BAF, an intensity-only probe is written; the BAF range). The gate for the pipeline:
an export without FORMAT/GT, FORMAT/LRR and FORMAT/BAF is refused."""
import json
import math
import os
import shutil
import subprocess

import pysam

from .genome import guess_assembly, normalize_contig

REQUIRED = ("GT", "LRR", "BAF")
INFO_OF_INTEREST = ("INTENSITY_ONLY", "GC", "ALLELE_A", "ALLELE_B", "GenTrain_Score", "Cluster_Sep")


def index_counts(path):
    """Records per contig from the index (`bcftools index -s`), {} without bcftools or an index."""
    exe = shutil.which("bcftools")
    if not exe:
        return {}
    try:
        out = subprocess.run([exe, "index", "-s", path], capture_output=True, text=True, check=True).stdout
    except (subprocess.CalledProcessError, OSError):
        return {}
    counts = {}
    for line in out.splitlines():
        parts = line.split("\t")
        if len(parts) >= 3:
            counts[parts[0]] = int(parts[2])
    return counts


def _isnan(x):
    return x is None or (isinstance(x, float) and math.isnan(x))


def inspect(path, records=500, samples=200):
    """A dict describing the export; `records` records (in file order) over the first `samples` samples are read for the
    conventions (None for all)."""
    vf = pysam.VariantFile(path)
    names = list(vf.header.samples)
    contigs = {c.name: c.length for c in vf.header.contigs.values()}
    formats = {k: dict(number=str(v.number), type=v.type, description=v.description) for k, v in vf.header.formats.items()}
    infos = {k: dict(number=str(v.number), type=v.type) for k, v in vf.header.info.items()}
    missing = [f for f in REQUIRED if f not in formats]
    asm = guess_assembly(contigs)
    chrom_prefixed = all(str(c).startswith("chr") for c in contigs) if contigs else None
    sub = names if samples is None else names[:samples]
    if sub != names:
        vf.subset_samples(sub)
    seen = dict(records=0, by_contig={}, gt_nocall=0, lrr_missing=0, baf_missing=0, baf_min=math.inf, baf_max=-math.inf, lrr_min=math.inf, lrr_max=-math.inf,
                intensity_only=0, alt_missing=0, multiallelic=0, non_snv=0, het=0, calls=0)
    if not missing:
        for rec in vf:
            if seen["records"] >= records:
                break
            seen["records"] += 1
            c = normalize_contig(rec.chrom)
            seen["by_contig"][c] = seen["by_contig"].get(c, 0) + 1
            alts = rec.alts or ()
            if not alts or all(a in (None, ".", "") for a in alts):      # ALT "." (an intensity-only probe): None from a VCF, "." from a BCF
                seen["alt_missing"] += 1
            elif len(alts) > 1:
                seen["multiallelic"] += 1
            elif len(rec.ref) != 1 or len(alts[0]) != 1:
                seen["non_snv"] += 1
            if "INTENSITY_ONLY" in infos and rec.info.get("INTENSITY_ONLY", False):
                seen["intensity_only"] += 1
            for s in sub:
                v = rec.samples[s]
                gt = v.get("GT")
                if gt is None or any(a is None for a in gt):
                    seen["gt_nocall"] += 1
                else:
                    seen["calls"] += 1
                    seen["het"] += int(len(set(gt)) > 1)
                lrr, baf = v.get("LRR"), v.get("BAF")
                if _isnan(lrr):
                    seen["lrr_missing"] += 1
                else:
                    seen["lrr_min"], seen["lrr_max"] = min(seen["lrr_min"], lrr), max(seen["lrr_max"], lrr)
                if _isnan(baf):
                    seen["baf_missing"] += 1
                else:
                    seen["baf_min"], seen["baf_max"] = min(seen["baf_min"], baf), max(seen["baf_max"], baf)
    for k in ("baf_min", "baf_max", "lrr_min", "lrr_max"):
        if math.isinf(seen[k]):
            seen[k] = None
    return dict(file=os.path.basename(path), samples=len(names), sample_names_head=names[:5], contigs=len(contigs), assembly=asm or "unknown",
                chrom_prefixed=chrom_prefixed, contig_lengths={normalize_contig(k): v for k, v in contigs.items()}, records_per_contig=index_counts(path),
                formats=formats, required_missing=missing, info_of_interest={k: infos[k] for k in INFO_OF_INTEREST if k in infos}, info_fields=sorted(infos),
                sampled=dict(seen, samples_read=len(sub)))


def describe(d):
    """The inspection in words."""
    lines = []
    lines.append("%s: %d samples (%s%s), %d contigs, assembly %s%s" % (
        d["file"], d["samples"], ", ".join(d["sample_names_head"]), ", ..." if d["samples"] > 5 else "", d["contigs"], d["assembly"],
        "" if d["chrom_prefixed"] is None else (" (chr-prefixed contigs)" if d["chrom_prefixed"] else " (contigs without the chr prefix)")))
    rpc = d["records_per_contig"]
    if rpc:
        total = sum(rpc.values())
        lines.append("records: %d in the index; per contig: %s" % (total, ", ".join("%s %d" % (normalize_contig(k), v) for k, v in rpc.items())))
    else:
        lines.append("records: no index statistics (bcftools not on PATH, or no .csi/.tbi beside the file)")
    f = d["formats"]
    lines.append("FORMAT: %s" % ", ".join("%s (%s %s)" % (k, v["number"], v["type"]) for k, v in f.items()))
    if d["required_missing"]:
        lines.append("REFUSED: the export lacks FORMAT/%s - arraykaryo needs GT, LRR and BAF" % ", FORMAT/".join(d["required_missing"]))
        return "\n".join(lines)
    lines.append("INFO of interest: %s" % (", ".join(d["info_of_interest"]) or "none of INTENSITY_ONLY, GC, ALLELE_A/B, GenTrain_Score, Cluster_Sep"))
    s = d["sampled"]
    if s["records"]:
        n = s["records"] * s["samples_read"]
        lines.append("conventions over %d records x %d samples: no-call GT %.2f%%, LRR missing %.2f%%, BAF missing %.2f%%; heterozygous %.1f%% of calls; "
                     "BAF range %s..%s, LRR range %s..%s; intensity-only probes %d, records without an ALT %d, multiallelic %d, non-SNV %d; contigs sampled: %s" % (
                         s["records"], s["samples_read"], 100.0 * s["gt_nocall"] / n, 100.0 * s["lrr_missing"] / n, 100.0 * s["baf_missing"] / n,
                         100.0 * s["het"] / max(s["calls"], 1), _f(s["baf_min"]), _f(s["baf_max"]), _f(s["lrr_min"]), _f(s["lrr_max"]),
                         s["intensity_only"], s["alt_missing"], s["multiallelic"], s["non_snv"], ", ".join("%s %d" % kv for kv in s["by_contig"].items())))
        if s["baf_max"] is not None and (s["baf_max"] > 1.0 or s["baf_min"] < 0.0):
            lines.append("note: BAF is not clipped to [0, 1] in this export (gtc2vcf writes it as computed); clip before any band statistic")
    return "\n".join(lines)


def _f(x):
    return "NA" if x is None else ("%.3f" % x)


def write_json(d, path):
    with open(path, "w") as fh:
        json.dump(d, fh, indent=1)
