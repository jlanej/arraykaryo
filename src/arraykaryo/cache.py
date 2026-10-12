"""The cache `arraykaryo extract` writes and every later command reads: per chromosome, memory-mapped arrays of shape
samples x probes - `gt` int8 (0 hom-ref, 1 het, 2 hom-alt, -1 no-call), `lrr` float32, `baf` float16 (BAF noise at a
heterozygous call is about 0.03, so float16's 0.0005 resolution costs nothing; measured in docs/METHODS.md) - plus
`probes.tsv` (and its binary twin `probes.npz`), `samples.txt` in the export's order, and `manifest.json`.

Sample-major layout: one sample's whole chromosome is a contiguous read (the per-sample figures and event calls), and a
block of probes over all samples (QC, the cohort reference) is a strided read of contiguous runs; the extraction writes
blocks of records, which the page cache merges."""
import os

import numpy as np

from .genome import genome_for, load_genome, normalize_contig
from .util import read_json, read_tsv

PROBE_COLS = ("chrom", "pos", "id", "ref", "alt", "snv", "intensity_only", "par", "gc", "gentrain", "cluster_sep", "allele_a", "allele_b")


class Probes:
    """One chromosome's probes: pos (1-based), id, ref, alt, snv (a biallelic single-base change), intensity_only (LRR usable, BAF
    and GT not), par (inside a pseudoautosomal region), gc, gentrain, cluster_sep (NaN when the export lacks them), and the
    sliced index into the whole-cohort probe table (`index`)."""

    def __init__(self, chrom, sl, npz, id_loader):
        self.chrom = chrom
        self.index = sl
        self.pos = npz["pos"][sl]
        self.snv = npz["snv"][sl]
        self.intensity_only = npz["intensity_only"][sl]
        self.par = npz["par"][sl]
        self.gc = npz["gc"][sl]
        self.gentrain = npz["gentrain"][sl]
        self.cluster_sep = npz["cluster_sep"][sl]
        self._ids, self._id_loader = None, id_loader

    @property
    def n(self):
        return int(self.pos.size)

    @property
    def ids(self):
        """The probe identifiers (read from probes.tsv on first use)."""
        if self._ids is None:
            self._ids = self._id_loader(self.chrom)
        return self._ids


class Cache:
    def __init__(self, path):
        self.path = path
        self.manifest = read_json(os.path.join(path, "manifest.json"))
        with open(os.path.join(path, "samples.txt")) as fh:
            self.samples = [l.rstrip("\n") for l in fh if l.strip()]
        self.n_samples = len(self.samples)
        self._index = {s: i for i, s in enumerate(self.samples)}
        self.genome = load_genome(self.manifest["genome_table"]) if self.manifest.get("genome_table") else load_genome(self.manifest["assembly"])
        self.chroms = [c for c in self.genome.chroms if c in self.manifest["chroms"]]
        self.n_probes = {c: int(self.manifest["chroms"][c]) for c in self.chroms}
        self._npz = None
        self._offsets = None
        self._ids = {}

    def sample_index(self, name):
        return self._index[name]

    @property
    def npz(self):
        if self._npz is None:
            z = np.load(os.path.join(self.path, "probes.npz"))
            self._npz = {k: z[k] for k in z.files}
            off, n = {}, 0
            for c in self.chroms:
                off[c] = slice(n, n + self.n_probes[c])
                n += self.n_probes[c]
            self._offsets = off
        return self._npz

    def probes(self, chrom):
        npz = self.npz
        return Probes(chrom, self._offsets[chrom], npz, self.probe_ids)

    def probe_ids(self, chrom):
        """The probe identifiers of one chromosome (read from probes.tsv on first use)."""
        if chrom not in self._ids:
            rows = read_tsv(os.path.join(self.path, "probes.tsv"))
            by = {}
            for r in rows:
                by.setdefault(r["chrom"], []).append(r["id"])
            for c in by:
                self._ids[c] = np.asarray(by[c])
        return self._ids.get(chrom, np.asarray([]))

    def probe_table(self):
        return read_tsv(os.path.join(self.path, "probes.tsv"))

    def _mm(self, kind, chrom, mode="r"):
        return np.load(os.path.join(self.path, "%s.%s.npy" % (kind, chrom)), mmap_mode=mode)

    def gt(self, chrom):
        return self._mm("gt", chrom)

    def lrr(self, chrom):
        return self._mm("lrr", chrom)

    def baf(self, chrom):
        return self._mm("baf", chrom)

    def chrom_done(self, chrom):
        return os.path.exists(os.path.join(self.path, "done", chrom + ".json"))


def array_path(cache_dir, kind, chrom):
    return os.path.join(cache_dir, "%s.%s.npy" % (kind, chrom))


class BlockWriter:
    """A sample-major .npy array written block by block with positioned writes (no memory-mapped dirty pages: the process's
    memory stays at one block). The file is created with numpy's header by open_memmap and then written through os.pwrite."""

    def __init__(self, path, n_samples, n_probes, dtype):
        self.path, self.dtype = path, np.dtype(dtype)
        self.shape = (n_samples, n_probes)
        mm = np.lib.format.open_memmap(path, mode="w+", dtype=self.dtype, shape=self.shape)
        self.offset = mm.offset
        del mm
        self.fd = os.open(path, os.O_WRONLY)
        self.row_bytes = n_probes * self.dtype.itemsize

    def write(self, i0, block):
        """block: (n_samples x k) values for probes i0..i0+k."""
        block = np.ascontiguousarray(block, dtype=self.dtype)
        k = block.shape[1]
        if block.shape[0] != self.shape[0] or i0 + k > self.shape[1]:
            raise ValueError("block does not fit the array")
        start = self.offset + i0 * self.dtype.itemsize
        for j in range(self.shape[0]):
            os.pwrite(self.fd, block[j].tobytes(), start + j * self.row_bytes)

    def close(self):
        os.fsync(self.fd)
        os.close(self.fd)


def open_arrays(cache_dir, chrom, n_samples, n_probes, baf_dtype=np.float16):
    """Fresh block writers for one chromosome (gt int8, lrr float32, baf float16 by default), renamed into place when complete."""
    tmp = lambda k: array_path(cache_dir, k, chrom) + ".part"  # noqa: E731
    return (BlockWriter(tmp("gt"), n_samples, n_probes, np.int8), BlockWriter(tmp("lrr"), n_samples, n_probes, np.float32), BlockWriter(tmp("baf"), n_samples, n_probes, baf_dtype))


def finish_arrays(cache_dir, chrom, arrays):
    for k, a in zip(("gt", "lrr", "baf"), arrays):
        a.close()
        os.replace(array_path(cache_dir, k, chrom) + ".part", array_path(cache_dir, k, chrom))
