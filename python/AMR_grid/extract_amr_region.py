#!/usr/bin/env python3
"""
extract_amr_region.py
=====================

Filter a generic AMR file to a spatial region while **preserving the
original BOXLEN**.  Cells whose center is inside the requested box are
kept at their original coordinates; cells outside are dropped.

Why preserve BOXLEN?
--------------------

The LaRT octree builder assumes cell positions lie on the octree's
natural dyadic grid (cx = origin + (2k+1) * BOXLEN/2^(n+1) at level n).
When cells come from a RAMSES simulation, they sit on the *original*
RAMSES boxlen's natural grid.  Shrinking BOXLEN to wrap a sub-region
moves the natural grid off the data — `amr_find_leaf` then fails for
positions that *do* contain cells, and pole traversal returns 0.

This tool keeps BOXLEN equal to the input file's BOXLEN, so the
surviving cells remain on grid.  The box outside the kept region is
sparse; LaRT's gap-skip logic (`amr_gap_exit`) walks through empty
space without crashing.

Streaming
---------

Designed for huge files (e.g. 130M cells, 7.5 GB).  Reads HDF5 in
chunks via h5py.

Usage
-----

    python extract_amr_region.py jellyfish.h5 \\
        --center 50 50 50 --size 4 \\
        -o jellyfish_region.h5
"""

from __future__ import annotations

import argparse
import sys
from typing import Iterable

import numpy as np


def extract_region(
    input_path: str,
    output_path: str,
    bounds,        # (xmin, xmax, ymin, ymax, zmin, zmax)
    recenter: bool = False,
    chunk: int = 5_000_000,
) -> dict:
    """Stream-filter cells in ``bounds`` to ``output_path``.

    The output preserves the input's BOXLEN, ORIGINX/Y/Z, UNITLCGS.
    If ``recenter`` is True, the cell coordinates and ORIGIN are shifted
    by -BOXLEN/2 (only valid when input ORIGIN = 0 and the shift moves
    the data onto the centered natural grid).
    """
    import h5py
    xmin, xmax, ymin, ymax, zmin, zmax = bounds

    with h5py.File(input_path, "r") as fin:
        # Locate the data group: AMRGRID (canonical) or 'section_001' (legacy)
        if "AMRGRID" in fin:
            gin = fin["AMRGRID"]
        elif "section_001" in fin:
            gin = fin["section_001"]
        else:
            gin = fin
        col_names = list(gin.keys())
        n = gin[col_names[0]].shape[0]
        boxlen = float(gin.attrs.get("BOXLEN", 0.0))
        unit_l = float(gin.attrs.get("UNITLCGS", 0.0))
        origin0 = (
            float(gin.attrs.get("ORIGINX", 0.0)),
            float(gin.attrs.get("ORIGINY", 0.0)),
            float(gin.attrs.get("ORIGINZ", 0.0)),
        )
        print(f"  input nleaf = {n}, BOXLEN = {boxlen}, ORIGIN = {origin0}")

        keep_arrays = {name: [] for name in col_names}
        kept = 0
        scanned = 0
        for start in range(0, n, chunk):
            end = min(start + chunk, n)
            x = gin["x"][start:end]
            y = gin["y"][start:end]
            z = gin["z"][start:end]
            mask = (
                (x >= xmin) & (x < xmax) &
                (y >= ymin) & (y < ymax) &
                (z >= zmin) & (z < zmax)
            )
            kept += int(mask.sum())
            scanned = end
            print(f"  scanned {scanned}/{n}, kept so far {kept}", flush=True)
            if mask.any():
                for name in col_names:
                    arr = gin[name][start:end][mask]
                    keep_arrays[name].append(arr)

        if kept == 0:
            raise SystemExit("no cells selected by the requested bounds")

        out_cols = {name: np.concatenate(arrs) for name, arrs in keep_arrays.items() if arrs}

    # Optional recenter (shift by -BOXLEN/2 so the box is [-BOXLEN/2, +BOXLEN/2])
    if recenter:
        shift = -0.5 * boxlen
        out_cols["x"] = out_cols["x"] + (shift - origin0[0])
        out_cols["y"] = out_cols["y"] + (shift - origin0[1])
        out_cols["z"] = out_cols["z"] + (shift - origin0[2])
        new_origin = (shift, shift, shift)
    else:
        new_origin = origin0

    # Write output
    if output_path.endswith(".h5") or output_path.endswith(".hdf5"):
        _write_h5(output_path, out_cols, boxlen, new_origin, unit_l)
    else:
        _write_fits(output_path, out_cols, boxlen, new_origin, unit_l)

    return {
        "nleaf_in":   n,
        "nleaf_out":  kept,
        "boxlen":     boxlen,
        "origin_in":  origin0,
        "origin_out": new_origin,
        "x_range":    (float(out_cols["x"].min()), float(out_cols["x"].max())),
        "y_range":    (float(out_cols["y"].min()), float(out_cols["y"].max())),
        "z_range":    (float(out_cols["z"].min()), float(out_cols["z"].max())),
    }


def _write_h5(path, cols, boxlen, origin, unit_l):
    import h5py
    if "gasDen" in cols and "nH" not in cols:
        cols["nH"] = cols.pop("gasDen")
    mandatory = ["x", "y", "z", "level", "nH", "T", "vx", "vy", "vz"]
    optional  = ["metallicity", "xHI", "n_e", "n_ion", "emissivity", "ndust"]
    for n in mandatory:
        if n not in cols:
            raise SystemExit(f"missing mandatory column '{n}'")
    with h5py.File(path, "w", libver="latest", track_order=True) as f:
        g = f.create_group("AMRGRID", track_order=True)
        for n in mandatory:
            arr = cols[n]
            kwargs = {}
            if arr.size > 4096:
                kwargs = {"chunks": (min(arr.size, 4096),),
                          "compression": "gzip", "compression_opts": 4}
            if n == "level":
                g.create_dataset(n, data=arr.astype(np.int32), **kwargs)
            else:
                g.create_dataset(n, data=arr.astype(np.float64), **kwargs)
        for n in optional:
            if n in cols:
                arr = cols[n]
                kwargs = {}
                if arr.size > 4096:
                    kwargs = {"chunks": (min(arr.size, 4096),),
                              "compression": "gzip", "compression_opts": 4}
                g.create_dataset(n, data=arr.astype(np.float64), **kwargs)
        g.attrs["BOXLEN"]  = float(boxlen)
        g.attrs["ORIGINX"] = float(origin[0])
        g.attrs["ORIGINY"] = float(origin[1])
        g.attrs["ORIGINZ"] = float(origin[2])
        g.attrs["NLEAF"]   = np.int32(len(cols["x"]))
        g.attrs["NAXIS2"]  = np.int32(len(cols["x"]))
        if unit_l > 0:
            g.attrs["UNITLCGS"] = float(unit_l)


def _write_fits(path, cols, boxlen, origin, unit_l):
    from astropy.io import fits
    if "gasDen" in cols and "nH" not in cols:
        cols["nH"] = cols.pop("gasDen")
    mandatory = ["x", "y", "z", "level", "nH", "T", "vx", "vy", "vz"]
    optional  = ["metallicity", "xHI", "n_e", "n_ion", "emissivity", "ndust"]
    fits_cols = []
    for n in mandatory:
        if n not in cols:
            raise SystemExit(f"missing mandatory column '{n}'")
        arr = cols[n]
        if n == "level":
            fits_cols.append(fits.Column(name=n, format="J", array=arr.astype(np.int32)))
        else:
            fits_cols.append(fits.Column(name=n, format="D", array=arr.astype(np.float64)))
    for n in optional:
        if n in cols:
            fits_cols.append(fits.Column(name=n, format="D", array=cols[n].astype(np.float64)))
    hdu = fits.BinTableHDU.from_columns(fits_cols, name="AMRGRID")
    hdu.header["BOXLEN"]  = (float(boxlen),    "AMR box length")
    hdu.header["ORIGINX"] = (float(origin[0]), "Box origin x")
    hdu.header["ORIGINY"] = (float(origin[1]), "Box origin y")
    hdu.header["ORIGINZ"] = (float(origin[2]), "Box origin z")
    hdu.header["NLEAF"]   = (len(cols["x"]),   "Number of leaf cells")
    if unit_l > 0:
        hdu.header["UNITLCGS"] = (float(unit_l), "RAMSES unit_l in cm")
    fits.HDUList([fits.PrimaryHDU(), hdu]).writeto(path, overwrite=True)


def main(argv: Iterable[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("input", help="HDF5 generic AMR file")
    ap.add_argument("-o", "--output", required=True)
    g1 = ap.add_argument_group("explicit bounds")
    g1.add_argument("--xmin", type=float); g1.add_argument("--xmax", type=float)
    g1.add_argument("--ymin", type=float); g1.add_argument("--ymax", type=float)
    g1.add_argument("--zmin", type=float); g1.add_argument("--zmax", type=float)
    g2 = ap.add_argument_group("center + size")
    g2.add_argument("--center", nargs=3, type=float, metavar=("CX", "CY", "CZ"))
    g2.add_argument("--size", type=float)
    ap.add_argument("--recenter", action="store_true",
                    help="also shift cells by -BOXLEN/2 so the output is centered. "
                         "Use only when input ORIGIN = 0 (i.e. cells in [0, BOXLEN]).")
    ap.add_argument("--chunk", type=int, default=5_000_000,
                    help="read chunk size (default 5M)")
    args = ap.parse_args(argv)

    if args.center is not None and args.size is not None:
        cx, cy, cz = args.center; h = 0.5 * args.size
        bounds = (cx-h, cx+h, cy-h, cy+h, cz-h, cz+h)
    else:
        for n in ("xmin","xmax","ymin","ymax","zmin","zmax"):
            if getattr(args, n) is None:
                raise SystemExit("specify either --center+--size or all --xmin..--zmax")
        bounds = (args.xmin, args.xmax, args.ymin, args.ymax, args.zmin, args.zmax)

    info = extract_region(args.input, args.output, bounds, recenter=args.recenter,
                          chunk=args.chunk)
    print(f"\n  input cells   : {info['nleaf_in']}")
    print(f"  output cells  : {info['nleaf_out']}")
    print(f"  BOXLEN (kept) : {info['boxlen']}")
    print(f"  input origin  : {info['origin_in']}")
    print(f"  output origin : {info['origin_out']}")
    print(f"  output x range: {info['x_range']}")
    print(f"  output y range: {info['y_range']}")
    print(f"  output z range: {info['z_range']}")
    print(f"  wrote         : {args.output}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
