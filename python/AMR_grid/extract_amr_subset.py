#!/usr/bin/env python3
"""
extract_amr_subset.py
=====================

Extract a cubic spatial sub-region from a generic-format AMR file
(.fits / .fits.gz / .h5) and write a new generic AMR file with the
selected leaves, recentered into a box centered at the origin.

The output uses the LaRT-preferred "centered" coordinate convention
(``ORIGINX = ORIGINY = ORIGINZ = -boxlen/2``), so that LaRT's box-center
based pole-traversal and ``par%xs_point = 0`` source default land at the
center of the extracted volume.

Behavior
--------

* Selection is by cell center (any cell whose centre falls inside the
  requested cube is kept).
* All optional columns present in the input (xHI, n_e, emissivity, etc.)
  are carried through unchanged.
* Cell *re-indexing*: the output is row-indexed 1..N_subset; the
  original simulation cell ID is **not** preserved (LaRT generic AMR
  files do not carry an explicit cell-ID column).  If you need to round-
  trip, add such a column to ``AMRGrid`` first.
* The new BOXLEN is the side length of the requested cube (must be
  cubic).  Coordinates are shifted from the input frame to a centered
  cube ``[-boxlen_new/2, +boxlen_new/2]``.

Usage
-----

    python extract_amr_subset.py input.fits.gz \\
        --xmin 0 --xmax 30 --ymin 0 --ymax 30 --zmin 0 --zmax 30 \\
        -o output.fits.gz

If ``--center`` and ``--size`` are given instead of explicit min/max,
they define the cube ``[center-size/2, center+size/2]``.
"""

from __future__ import annotations

import argparse
import os
import sys
from typing import Iterable, Tuple

import numpy as np
from astropy.io import fits


# Optional columns (must match read_generic_amr.f90 and AMRGrid.py)
OPTIONAL_COLUMNS = (
    "metallicity",
    "xHI",
    "n_e",
    "n_ion",
    "emissivity",
    "ndust",
)


def _open_table(path: str):
    """Return (table-data, header) for a .fits/.fits.gz file or HDF5 file."""
    lower = path.lower()
    if lower.endswith(".h5") or lower.endswith(".hdf5"):
        try:
            import h5py
        except ImportError as exc:  # pragma: no cover
            raise SystemExit(
                "h5py is required to read HDF5 generic AMR files"
            ) from exc
        with h5py.File(path, "r") as f:
            g = f["amr_grid"] if "amr_grid" in f else f
            cols = {name: np.asarray(g[name]) for name in g.keys()}
            header = {key: g.attrs[key] for key in g.attrs.keys()}
        return cols, header

    # FITS / FITS.gz
    with fits.open(path) as hdul:
        hdu = hdul[1]
        data = hdu.data
        cols = {name: np.array(data[name]) for name in data.dtype.names}
        header = {k: hdu.header[k] for k in hdu.header.keys()}
    return cols, header


def _write_fits(path: str, cols, header_extra) -> None:
    """Write a generic AMR FITS file with one binary-table HDU."""
    names = ["x", "y", "z", "level", "nH", "T", "vx", "vy", "vz"]
    # Map the conventional 'gasDen' column to 'nH' if needed.
    src_names = list(cols.keys())
    if "gasDen" in src_names and "nH" not in src_names:
        cols["nH"] = cols.pop("gasDen")
    # Re-order mandatory then optional columns.
    out_order = []
    for n in names:
        if n in cols:
            out_order.append(n)
        else:
            raise SystemExit(f"input is missing mandatory column '{n}'")
    for n in OPTIONAL_COLUMNS:
        if n in cols:
            out_order.append(n)
    fits_cols = []
    for n in out_order:
        arr = cols[n]
        if n == "level":
            fits_cols.append(fits.Column(name=n, format="J", array=arr.astype(np.int32)))
        else:
            fits_cols.append(fits.Column(name=n, format="D", array=arr.astype(np.float64)))
    hdu = fits.BinTableHDU.from_columns(fits_cols)
    for k, v in header_extra.items():
        hdu.header[k] = v
    hdul = fits.HDUList([fits.PrimaryHDU(), hdu])
    hdul.writeto(path, overwrite=True)


def extract_subset(
    input_path: str,
    output_path: str,
    bounds: Tuple[float, float, float, float, float, float],
) -> dict:
    """Read ``input_path``, keep cells whose centre is inside ``bounds``,
    shift to a centered box, and write ``output_path``.  Returns a dict
    with summary statistics.
    """
    xmin, xmax, ymin, ymax, zmin, zmax = bounds
    dx, dy, dz = xmax - xmin, ymax - ymin, zmax - zmin
    if not (dx > 0 and dy > 0 and dz > 0):
        raise SystemExit("bounds must satisfy xmin < xmax (etc.)")
    if not np.isclose(dx, dy) or not np.isclose(dy, dz):
        raise SystemExit(
            f"Subset must be cubic (got dx={dx}, dy={dy}, dz={dz}). "
            f"LaRT amr_grid only supports cubic boxes."
        )
    boxlen_new = float(dx)

    cols, _ = _open_table(input_path)
    if "x" not in cols or "y" not in cols or "z" not in cols:
        raise SystemExit(f"input {input_path} has no x/y/z columns")

    x = cols["x"]
    y = cols["y"]
    z = cols["z"]
    mask = (
        (x >= xmin)
        & (x < xmax)
        & (y >= ymin)
        & (y < ymax)
        & (z >= zmin)
        & (z < zmax)
    )
    n_in  = len(x)
    n_out = int(mask.sum())
    if n_out == 0:
        raise SystemExit("no cells fall inside the requested subset")

    cube_center = np.array([0.5 * (xmin + xmax),
                            0.5 * (ymin + ymax),
                            0.5 * (zmin + zmax)])

    # Shift to centered frame.  Output positions are in
    # [-boxlen_new/2, +boxlen_new/2].
    out_cols = {}
    for key, arr in cols.items():
        sub = arr[mask]
        if   key == "x": sub = sub - cube_center[0]
        elif key == "y": sub = sub - cube_center[1]
        elif key == "z": sub = sub - cube_center[2]
        out_cols[key] = sub

    header_extra = {
        "BOXLEN":  (boxlen_new, "Side length of the extracted cube"),
        "ORIGINX": (-0.5 * boxlen_new, "Lower corner x (centered)"),
        "ORIGINY": (-0.5 * boxlen_new, "Lower corner y (centered)"),
        "ORIGINZ": (-0.5 * boxlen_new, "Lower corner z (centered)"),
        "NLEAF":   (n_out, "Number of leaf cells"),
        "SUBSET":  (1, "1 = extracted from a larger AMR volume"),
        "SUBCNTX": (cube_center[0], "Subset center x (in input frame)"),
        "SUBCNTY": (cube_center[1], "Subset center y (in input frame)"),
        "SUBCNTZ": (cube_center[2], "Subset center z (in input frame)"),
    }

    if output_path.endswith(".h5") or output_path.endswith(".hdf5"):
        _write_hdf5(output_path, out_cols, header_extra)
    else:
        _write_fits(output_path, out_cols, header_extra)

    return {
        "n_in":         n_in,
        "n_out":        n_out,
        "boxlen_new":   boxlen_new,
        "cube_center":  tuple(float(v) for v in cube_center),
        "output_path":  output_path,
    }


def _write_hdf5(path: str, cols, header_extra) -> None:
    """Write an HDF5 generic AMR file."""
    try:
        import h5py
    except ImportError as exc:  # pragma: no cover
        raise SystemExit("h5py is required to write HDF5 output") from exc

    # Map gasDen -> nH if needed
    if "gasDen" in cols and "nH" not in cols:
        cols["nH"] = cols.pop("gasDen")

    with h5py.File(path, "w") as f:
        g = f.create_group("amr_grid")
        order = ["x", "y", "z", "level", "nH", "T", "vx", "vy", "vz"]
        for n in order:
            if n not in cols:
                raise SystemExit(f"missing mandatory column '{n}'")
            g.create_dataset(n, data=cols[n])
        for n in OPTIONAL_COLUMNS:
            if n in cols:
                g.create_dataset(n, data=cols[n])
        for k, v in header_extra.items():
            value = v[0] if isinstance(v, tuple) else v
            g.attrs[k] = value


def main(argv: Iterable[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("input", help="input generic AMR file (.fits / .fits.gz / .h5)")
    ap.add_argument("-o", "--output", required=True, help="output file path")
    # Two ways to specify the cube: explicit min/max, or center+size.
    g = ap.add_argument_group("explicit bounds")
    g.add_argument("--xmin", type=float)
    g.add_argument("--xmax", type=float)
    g.add_argument("--ymin", type=float)
    g.add_argument("--ymax", type=float)
    g.add_argument("--zmin", type=float)
    g.add_argument("--zmax", type=float)
    c = ap.add_argument_group("center + size")
    c.add_argument("--center", nargs=3, type=float, metavar=("CX", "CY", "CZ"),
                   help="cube center (in input coordinates)")
    c.add_argument("--size", type=float, help="cube side length")
    args = ap.parse_args(argv)

    if args.center is not None and args.size is not None:
        cx, cy, cz = args.center
        h = 0.5 * args.size
        bounds = (cx - h, cx + h, cy - h, cy + h, cz - h, cz + h)
    else:
        for name in ("xmin", "xmax", "ymin", "ymax", "zmin", "zmax"):
            if getattr(args, name) is None:
                raise SystemExit(
                    "specify either all of --xmin/xmax/--ymin/ymax/--zmin/zmax "
                    "or both --center and --size"
                )
        bounds = (args.xmin, args.xmax, args.ymin, args.ymax, args.zmin, args.zmax)

    info = extract_subset(args.input, args.output, bounds)
    print(f"  input cells   : {info['n_in']}")
    print(f"  output cells  : {info['n_out']}")
    print(f"  new BOXLEN    : {info['boxlen_new']}")
    print(f"  subset center : {info['cube_center']}  (in input frame)")
    print(f"  output frame  : centered, [-{info['boxlen_new']/2}, +{info['boxlen_new']/2}]")
    print(f"  wrote         : {info['output_path']}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
