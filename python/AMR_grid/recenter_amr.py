#!/usr/bin/env python3
"""
recenter_amr.py
===============

Re-write a generic AMR file (.fits / .fits.gz / .h5) so that the box
is centered on the origin, i.e.

    ORIGINX = ORIGINY = ORIGINZ = -BOXLEN/2
    cell coordinates are shifted from [origin, origin+BOXLEN]
                                  to  [-BOXLEN/2, +BOXLEN/2]

This is the convention preferred by LaRT (matches the Cartesian grid
and `par%xs_point = 0` placing the source at the box center).  All
optional physics columns (xHI, n_e, emissivity, ...) are carried
through unchanged; the BOXLEN value is preserved.

Cell index handling
-------------------

Cells are written in the same row order as in the input file.  No
re-indexing is done -- there is no notion of an "original cell ID"
in the LaRT generic format, so the (row index) -> (leaf index)
mapping is preserved bit-for-bit.

Usage
-----

    python recenter_amr.py input.fits.gz -o output.fits.gz
    python recenter_amr.py input.h5      -o output.h5
"""

from __future__ import annotations

import argparse
import sys
from typing import Iterable, Tuple

import numpy as np
from astropy.io import fits


OPTIONAL_COLUMNS = (
    "metallicity",
    "xHI",
    "n_e",
    "n_ion",
    "emissivity",
    "ndust",
)


def _read_input(path: str) -> Tuple[dict, dict]:
    lower = path.lower()
    if lower.endswith(".h5") or lower.endswith(".hdf5"):
        try:
            import h5py
        except ImportError as exc:
            raise SystemExit("h5py is required for HDF5 input") from exc
        with h5py.File(path, "r") as f:
            # The file may use either a top-level group 'AMRGRID'
            # (converter output) or a flat layout.
            g = f["AMRGRID"] if "AMRGRID" in f else f
            cols   = {name: np.asarray(g[name]) for name in g.keys()}
            header = {key: g.attrs[key] for key in g.attrs.keys()}
        return cols, header

    with fits.open(path) as hdul:
        hdu = hdul[1]
        data = hdu.data
        cols   = {name: np.array(data[name]) for name in data.dtype.names}
        header = {k: hdu.header[k] for k in hdu.header.keys()}
    return cols, header


def _write_fits(path: str, cols: dict, header_extra: dict) -> None:
    if "gasDen" in cols and "nH" not in cols:
        cols["nH"] = cols.pop("gasDen")

    mandatory = ["x", "y", "z", "level", "nH", "T", "vx", "vy", "vz"]
    for n in mandatory:
        if n not in cols:
            raise SystemExit(f"input is missing mandatory column '{n}'")

    fits_cols = []
    for n in mandatory:
        arr = cols[n]
        if n == "level":
            fits_cols.append(fits.Column(name=n, format="J",
                                         array=arr.astype(np.int32)))
        else:
            fits_cols.append(fits.Column(name=n, format="D",
                                         array=arr.astype(np.float64)))
    for n in OPTIONAL_COLUMNS:
        if n in cols:
            fits_cols.append(fits.Column(name=n, format="D",
                                         array=cols[n].astype(np.float64)))

    hdu = fits.BinTableHDU.from_columns(fits_cols, name="AMRGRID")
    for k, v in header_extra.items():
        hdu.header[k] = v
    fits.HDUList([fits.PrimaryHDU(), hdu]).writeto(path, overwrite=True)


def _write_hdf5(path: str, cols: dict, header_extra: dict) -> None:
    try:
        import h5py
    except ImportError as exc:
        raise SystemExit("h5py is required for HDF5 output") from exc

    if "gasDen" in cols and "nH" not in cols:
        cols["nH"] = cols.pop("gasDen")

    mandatory = ["x", "y", "z", "level", "nH", "T", "vx", "vy", "vz"]
    for n in mandatory:
        if n not in cols:
            raise SystemExit(f"input is missing mandatory column '{n}'")

    with h5py.File(path, "w", libver="latest", track_order=True) as f:
        g = f.create_group("AMRGRID", track_order=True)
        for n in mandatory:
            arr = cols[n]
            if n == "level":
                g.create_dataset(n, data=arr.astype(np.int32))
            else:
                g.create_dataset(n, data=arr.astype(np.float64))
        for n in OPTIONAL_COLUMNS:
            if n in cols:
                g.create_dataset(n, data=cols[n].astype(np.float64))
        for k, v in header_extra.items():
            if isinstance(v, tuple):
                v = v[0]
            g.attrs[k] = v


def recenter(
    input_path: str,
    output_path: str,
    mode: str = "auto",
    origin: Tuple[float, float, float] | None = None,
    boxlen_new: float | None = None,
) -> dict:
    """Re-center a generic AMR file.

    Parameters
    ----------
    mode :
        ``'data-center'`` (default) -- ignore the file's ORIGIN headers.
            Shift the cell coordinates so the **data center**
            ``0.5 * (min + max)`` lands at the origin (0,0,0).  Output
            box is ``[-BOXLEN/2, +BOXLEN/2]``.  Works for both
            corner-based and already-centered inputs.
        ``'from-header'`` -- trust the ORIGINX/Y/Z headers; shift by
            ``-ORIGIN - BOXLEN/2``.  Use only when you are sure the
            headers are correct (the legacy converter wrote ORIGINX=0
            for corner-based, with no automatic data-extent check).
        ``'explicit'`` -- shift by ``-origin`` (caller-supplied
            ``origin`` tuple).
    boxlen_new :
        If given, overrides BOXLEN in the output.  Useful when the
        file's BOXLEN is misleading (e.g., a subset of a larger box).
        Default: keep the input BOXLEN.
    """
    cols, header = _read_input(input_path)

    if "BOXLEN" not in header:
        raise SystemExit("input file does not carry a BOXLEN keyword")
    boxlen = float(header["BOXLEN"])

    cur_ox = float(header.get("ORIGINX", 0.0)) if header.get("ORIGINX") is not None else None
    cur_oy = float(header.get("ORIGINY", 0.0)) if header.get("ORIGINY") is not None else None
    cur_oz = float(header.get("ORIGINZ", 0.0)) if header.get("ORIGINZ") is not None else None

    if mode == "auto":
        # Heuristic:
        #   - if any coordinate is negative, the data is already (close to)
        #     centered  -->  shift = 0.
        #   - if every coordinate is non-negative AND fits within [0, BOXLEN],
        #     the data is corner-based  -->  shift = -BOXLEN/2.
        #   - otherwise, fall back to a data-center shift but warn.
        xmin, ymin, zmin = float(cols["x"].min()), float(cols["y"].min()), float(cols["z"].min())
        xmax, ymax, zmax = float(cols["x"].max()), float(cols["y"].max()), float(cols["z"].max())
        eps = 1e-6 * boxlen
        all_negative_lower_bound = (xmin < -eps) or (ymin < -eps) or (zmin < -eps)
        all_corner_based = (
            xmin >= -eps and xmax <= boxlen + eps and
            ymin >= -eps and ymax <= boxlen + eps and
            zmin >= -eps and zmax <= boxlen + eps and
            not all_negative_lower_bound
        )
        if all_negative_lower_bound:
            print("  auto: data already centered (some coord < 0); no shift applied.")
            shift_x = shift_y = shift_z = 0.0
        elif all_corner_based:
            print("  auto: data is corner-based in [0, BOXLEN]; shifting by -BOXLEN/2.")
            shift_x = shift_y = shift_z = -0.5 * boxlen
        else:
            cx_d = 0.5 * (xmin + xmax)
            cy_d = 0.5 * (ymin + ymax)
            cz_d = 0.5 * (zmin + zmax)
            print(f"  auto: ambiguous data layout, using data-center shift "
                  f"(-{cx_d:.6g}, -{cy_d:.6g}, -{cz_d:.6g}).")
            print( "        WARNING: this may move cells off the octree natural grid.")
            shift_x, shift_y, shift_z = -cx_d, -cy_d, -cz_d
    elif mode == "data-center":
        # Shift the data center to (0, 0, 0).  May move cells off the
        # octree-natural grid -- use only when you know what you are doing.
        cx = 0.5 * (float(cols["x"].min()) + float(cols["x"].max()))
        cy = 0.5 * (float(cols["y"].min()) + float(cols["y"].max()))
        cz = 0.5 * (float(cols["z"].min()) + float(cols["z"].max()))
        shift_x = -cx
        shift_y = -cy
        shift_z = -cz
    elif mode == "from-header":
        ox = cur_ox if cur_ox is not None else 0.0
        oy = cur_oy if cur_oy is not None else 0.0
        oz = cur_oz if cur_oz is not None else 0.0
        new_origin = -0.5 * boxlen
        shift_x = new_origin - ox
        shift_y = new_origin - oy
        shift_z = new_origin - oz
    elif mode == "explicit":
        if origin is None:
            raise SystemExit("mode='explicit' requires --origin OX OY OZ")
        ox, oy, oz = origin
        new_origin = -0.5 * boxlen
        shift_x = new_origin - ox
        shift_y = new_origin - oy
        shift_z = new_origin - oz
    else:
        raise SystemExit(f"unknown mode '{mode}'")

    cols["x"] = cols["x"] + shift_x
    cols["y"] = cols["y"] + shift_y
    cols["z"] = cols["z"] + shift_z

    new_origin = -0.5 * boxlen   # the lower corner in centered convention

    # If the data extent is much smaller than BOXLEN (subset of a larger
    # box), warn the user and offer the tighter cubic enclosing box
    # (sized to the data extent, still centered).
    data_extent = max(
        float(cols["x"].max() - cols["x"].min()),
        float(cols["y"].max() - cols["y"].min()),
        float(cols["z"].max() - cols["z"].min()),
    )
    if boxlen_new is not None:
        boxlen = float(boxlen_new)
        new_origin = -0.5 * boxlen
    elif data_extent < 0.5 * boxlen:
        print(f"  NOTE: data extent ({data_extent:.3g}) is < BOXLEN/2 ({0.5*boxlen:.3g}).")
        print(f"        Consider passing --boxlen to shrink the box around the data;")
        print(f"        otherwise pole traversal starts in empty space and tau falls")
        print(f"        back to tauhomo.")

    # Carry over header but override origin keys.  Drop unknown
    # FITS-internal keys that BinTableHDU writes by itself.
    drop_keys = {"XTENSION", "BITPIX", "NAXIS", "NAXIS1", "NAXIS2",
                 "PCOUNT", "GCOUNT", "TFIELDS", "EXTNAME",
                 "ORIGINX", "ORIGINY", "ORIGINZ"}
    drop_prefixes = ("TTYPE", "TFORM", "TUNIT", "TDISP", "TNULL")
    header_extra = {}
    for k, v in header.items():
        if k in drop_keys:
            continue
        if any(k.startswith(p) for p in drop_prefixes):
            continue
        header_extra[k] = v
    header_extra["BOXLEN"]  = (boxlen,    "Side length of the AMR box")
    header_extra["ORIGINX"] = (new_origin, "Box origin x (centered)")
    header_extra["ORIGINY"] = (new_origin, "Box origin y (centered)")
    header_extra["ORIGINZ"] = (new_origin, "Box origin z (centered)")
    header_extra["RECENTRD"] = (1, "Recentered by recenter_amr.py")
    header_extra["NLEAF"]   = (len(cols["x"]), "Number of leaf cells")

    if output_path.endswith(".h5") or output_path.endswith(".hdf5"):
        _write_hdf5(output_path, cols, header_extra)
    else:
        _write_fits(output_path, cols, header_extra)

    return {
        "boxlen":     boxlen,
        "old_origin": (cur_ox, cur_oy, cur_oz),
        "new_origin": (new_origin,) * 3,
        "shift":      (shift_x, shift_y, shift_z),
        "nleaf":      len(cols["x"]),
        "x_range":    (float(cols["x"].min()), float(cols["x"].max())),
        "y_range":    (float(cols["y"].min()), float(cols["y"].max())),
        "z_range":    (float(cols["z"].min()), float(cols["z"].max())),
    }


def main(argv: Iterable[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("input",  help="input generic AMR file (.fits / .fits.gz / .h5)")
    ap.add_argument("-o", "--output", required=True, help="output file path")
    ap.add_argument("--mode", choices=("auto", "data-center", "from-header", "explicit"),
                    default="auto",
                    help="how to determine the current lower-corner origin (default: auto)")
    ap.add_argument("--origin", nargs=3, type=float, metavar=("OX", "OY", "OZ"),
                    help="explicit current origin (use with --mode explicit)")
    ap.add_argument("--boxlen", type=float,
                    help="override output BOXLEN (default: keep input BOXLEN)")
    args = ap.parse_args(argv)

    info = recenter(args.input, args.output,
                    mode=args.mode,
                    origin=tuple(args.origin) if args.origin else None,
                    boxlen_new=args.boxlen)
    print(f"  mode        : {args.mode}")
    print(f"  nleaf       : {info['nleaf']}")
    print(f"  boxlen      : {info['boxlen']}")
    print(f"  old origin  : {info['old_origin']}")
    print(f"  new origin  : {info['new_origin']}  (centered)")
    print(f"  shift       : {info['shift']}")
    print(f"  data x range: {info['x_range']}")
    print(f"  data y range: {info['y_range']}")
    print(f"  data z range: {info['z_range']}")
    print(f"  wrote       : {args.output}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
