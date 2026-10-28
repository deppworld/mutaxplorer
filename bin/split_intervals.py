#!/usr/bin/env python3
"""
MutaXplorer - build calling intervals and split them into N balanced chunks for scatter/gather.

  WGS        : --fai genome.fa.fai  -> primary chromosomes (1-22, X, Y), optional --exclude BED
  WES/panel  : --bed targets.bed    -> cleaned, sorted, merged, padded targets

Writes:  <prefix>.all.bed          all intervals (used for coverage / -L)
         <prefix>.chunk_###.bed    N chunks of ~equal size (never splitting a target region
                                    for WES/panel; WGS contigs are cut into pieces)
"""
import argparse, re, sys

PRIMARY = re.compile(r"^(chr)?([0-9]{1,2}|X|Y)$")


def read_bed(path, pad, lengths):
    iv = []
    with open(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            f = line.split("\t")
            c, s, e = f[0], int(f[1]), int(f[2])
            if lengths and c not in lengths:
                alt = c[3:] if c.startswith("chr") else "chr" + c
                if alt in lengths:
                    c = alt
                else:
                    sys.exit(f"ERROR: BED contig '{f[0]}' is not in the reference")
            iv.append((c, max(0, s - pad), min(lengths.get(c, e + pad), e + pad) if lengths else e + pad))
    return iv


def merge(iv, order):
    iv.sort(key=lambda x: (order.get(x[0], 10**6), x[0], x[1]))
    out = []
    for c, s, e in iv:
        if out and out[-1][0] == c and s <= out[-1][2]:
            out[-1][2] = max(out[-1][2], e)
        else:
            out.append([c, s, e])
    return [tuple(x) for x in out]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fai", required=True)
    ap.add_argument("--bed", help="target BED (WES/panel); omit for WGS")
    ap.add_argument("--exclude", help="BED of regions to drop (WGS), e.g. ENCODE blacklist / N gaps")
    ap.add_argument("--pad", type=int, default=0)
    ap.add_argument("--chunks", type=int, default=1)
    ap.add_argument("--prefix", default="intervals")
    a = ap.parse_args()

    lengths, order = {}, {}
    with open(a.fai) as fh:
        for i, line in enumerate(fh):
            c, l = line.split("\t")[:2]
            lengths[c] = int(l); order[c] = i

    if a.bed:
        iv = merge(read_bed(a.bed, a.pad, lengths), order)
        splittable = False
    else:
        iv = [(c, 0, l) for c, l in lengths.items() if PRIMARY.match(c)]
        if not iv:
            sys.exit("ERROR: no primary chromosomes (1-22,X,Y) found in the .fai")
        if a.exclude:
            ex = merge(read_bed(a.exclude, 0, lengths), order)
            keep = []
            for c, s, e in iv:
                cur = s
                for xc, xs, xe in ex:
                    if xc != c or xe <= cur or xs >= e:
                        continue
                    if xs > cur:
                        keep.append((c, cur, xs))
                    cur = max(cur, xe)
                if cur < e:
                    keep.append((c, cur, e))
            iv = keep
        splittable = True

    with open(f"{a.prefix}.all.bed", "w") as fh:
        for c, s, e in iv:
            fh.write(f"{c}\t{s}\t{e}\n")

    n = max(1, a.chunks)
    total = sum(e - s for _, s, e in iv)
    target = total / n
    chunks, cur, size = [], [], 0
    for c, s, e in iv:
        while splittable and size + (e - s) > target * 1.05 and len(chunks) < n - 1:
            cut = s + int(target - size)
            if cut <= s:
                break
            cur.append((c, s, cut)); chunks.append(cur); cur, size, s = [], 0, cut
        cur.append((c, s, e)); size += e - s
        if not splittable and size >= target and len(chunks) < n - 1:
            chunks.append(cur); cur, size = [], 0
    if cur:
        chunks.append(cur)
    for i, ch in enumerate(chunks, 1):
        with open(f"{a.prefix}.chunk_{i:03d}.bed", "w") as fh:
            for c, s, e in ch:
                fh.write(f"{c}\t{s}\t{e}\n")
    print(f"{len(iv)} intervals, {total/1e6:.2f} Mb -> {len(chunks)} chunks", file=sys.stderr)


if __name__ == "__main__":
    main()
