#!/usr/bin/env python3
"""lart_io.py — format-agnostic LaRT I/O.

This module is the Python counterpart of the Fortran ``iofile_mod`` facade.
It reads LaRT output (or input) files in either FITS (.fits, .fits.gz) or
HDF5 (.h5, .hdf5) form and exposes a uniform structure, and it converts
between the two formats by re-emitting the same logical content under the
target schema.

Layout mapping (matches what LaRT writes through ``iofile_mod``):

  FITS                                  HDF5
  -----------------------------------------------------------------
  Empty Primary HDU                     root group, attrs only
  Image Primary HDU (EXTNAME=X)         /X/data dataset + group attrs
  Image extension HDU (EXTNAME=X)       /X/data dataset + group attrs
  BinTable HDU (EXTNAME=X)              /X/<colA>, /X/<colB>, ...
  HDU header keyword K = v              attribute @K on the group

Public API
----------
``load_lart(path)``
    Read either a .fits[.gz] or a .h5/.hdf5 file and return a
    :class:`LartFile` whose ``sections`` list mirrors the ordered HDU /
    group sequence.

``convert(src, dst)``
    Convert between formats by extension.  ``in.fits.gz`` → ``out.h5``
    and ``in.h5`` → ``out.fits.gz`` are both supported.  Round-trip
    preserves data within numerical precision (storage-type rounding only).

Command line
------------
    python lart_io.py info  <file>
    python lart_io.py convert <src> <dst>
"""
from __future__ import annotations

import os
import sys
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, Tuple

import numpy as np


# ---------------------------------------------------------------------------
# Format detection
# ---------------------------------------------------------------------------

_HDF5_EXTS = ('.h5', '.hdf5')
_FITS_EXTS = ('.fits', '.fits.gz', '.fit')


def detect_format(path: str) -> str:
    """Return 'fits' or 'hdf5' based on the filename extension."""
    lower = path.lower()
    if lower.endswith(_HDF5_EXTS):
        return 'hdf5'
    for ext in _FITS_EXTS:
        if lower.endswith(ext):
            return 'fits'
    # Last resort: try a magic-byte sniff
    try:
        with open(path, 'rb') as fh:
            head = fh.read(8)
        if head[:4] == b'\x89HDF':
            return 'hdf5'
        if head[:6] == b'SIMPLE':
            return 'fits'
    except OSError:
        pass
    raise ValueError(f"Cannot determine LaRT format for {path!r}")


# ---------------------------------------------------------------------------
# In-memory representation
# ---------------------------------------------------------------------------

@dataclass
class Section:
    """One logical block — either an image HDU or a BinTable HDU.

    ``name`` is the EXTNAME (or auto-generated ``section_NNN``).  ``data`` is
    a numpy array for image sections; for tables it is ``None`` and the
    columns live in ``columns`` (an ordered dict).  ``attrs`` collects
    metadata.
    """
    name: str
    data: Optional[np.ndarray] = None
    columns: Dict[str, np.ndarray] = field(default_factory=dict)
    attrs: Dict[str, Any] = field(default_factory=dict)

    @property
    def is_table(self) -> bool:
        return self.data is None and len(self.columns) > 0


@dataclass
class LartFile:
    path: str
    fmt: str                        # 'fits' or 'hdf5'
    sections: List[Section] = field(default_factory=list)
    root_attrs: Dict[str, Any] = field(default_factory=dict)

    def section(self, name: str) -> Optional[Section]:
        for s in self.sections:
            if s.name == name:
                return s
        return None

    def __getitem__(self, key: str) -> Section:
        s = self.section(key)
        if s is None:
            raise KeyError(key)
        return s

    def info(self) -> str:
        lines = [f'LartFile: {self.path}  ({self.fmt})']
        if self.root_attrs:
            lines.append('  root attrs: ' + ', '.join(sorted(self.root_attrs)))
        for s in self.sections:
            if s.is_table:
                cols = ', '.join(s.columns.keys())
                lines.append(f'  /{s.name} [table] cols=[{cols}]  attrs={len(s.attrs)}')
            elif s.data is not None:
                lines.append(f'  /{s.name} [image] shape={s.data.shape} dtype={s.data.dtype}  attrs={len(s.attrs)}')
            else:
                lines.append(f'  /{s.name} [empty]  attrs={len(s.attrs)}')
        return '\n'.join(lines)


# ---------------------------------------------------------------------------
# FITS reader
# ---------------------------------------------------------------------------

def _load_fits(path: str) -> LartFile:
    from astropy.io import fits

    out = LartFile(path=path, fmt='fits')
    with fits.open(path) as hdul:
        for i, hdu in enumerate(hdul):
            attrs = _fits_header_to_attrs(hdu.header)
            extname = attrs.pop('EXTNAME', None) or hdu.name or f'section_{i+1:03d}'

            if isinstance(hdu, fits.BinTableHDU):
                cols = {}
                for cname in hdu.columns.names:
                    cols[cname] = np.array(hdu.data[cname])
                sec = Section(name=extname, columns=cols, attrs=attrs)
            elif hdu.data is None or hdu.data.size == 0:
                # Empty Primary HDU; preserve its attrs at root and skip
                if i == 0:
                    out.root_attrs.update(attrs)
                    continue
                sec = Section(name=extname, attrs=attrs)
            else:
                sec = Section(name=extname, data=np.array(hdu.data), attrs=attrs)
            out.sections.append(sec)
    return out


def _fits_header_to_attrs(header) -> Dict[str, Any]:
    """Extract user-relevant FITS header keywords as a plain dict.

    Drops structural keywords (BITPIX, NAXIS*, EXTEND, ...) that the HDF5
    side reconstructs from the dataset shape/dtype itself.
    """
    skip_exact = {
        'SIMPLE', 'BITPIX', 'NAXIS', 'EXTEND', 'BSCALE', 'BZERO',
        'XTENSION', 'PCOUNT', 'GCOUNT', 'TFIELDS', 'CHECKSUM', 'DATASUM',
    }
    skip_prefix = ('NAXIS', 'TFORM', 'TTYPE', 'TUNIT', 'TDIM', 'COMMENT', 'HISTORY')
    out = {}
    for card in header.cards:
        k = card.keyword
        if not k:
            continue
        if k in skip_exact:
            continue
        if any(k.startswith(p) for p in skip_prefix) and k not in ('NAXIS1', 'NAXIS2', 'NAXIS3', 'NAXIS4'):
            # NAXIS itself is structural; NAXIS1.. for tables get exposed via shape too
            if k.startswith(('NAXIS', 'TFORM', 'TTYPE', 'TUNIT', 'TDIM', 'COMMENT', 'HISTORY')):
                continue
        # Keep value; FITS sometimes stores as fits Undefined — skip those.
        v = card.value
        if v is None:
            continue
        out[k] = v
    return out


# ---------------------------------------------------------------------------
# HDF5 reader
# ---------------------------------------------------------------------------

def _load_hdf5(path: str) -> LartFile:
    import h5py

    out = LartFile(path=path, fmt='hdf5')
    with h5py.File(path, 'r') as f:
        # Root-level attributes (if any)
        for k, v in f.attrs.items():
            out.root_attrs[k] = _h5_attr_to_py(v)

        # Walk children in creation order if available, else by name.
        names = _h5_children_in_order(f)
        for nm in names:
            obj = f[nm]
            if not isinstance(obj, h5py.Group):
                # A bare dataset at root: treat as an image section named after the dataset.
                sec = Section(name=nm, data=np.array(obj), attrs={})
                out.sections.append(sec)
                continue
            attrs = {k: _h5_attr_to_py(v) for k, v in obj.attrs.items()}
            children = list(obj.keys())
            if children == ['data']:
                sec = Section(name=nm, data=np.array(obj['data']), attrs=attrs)
            else:
                cols = {}
                for cn in children:
                    if isinstance(obj[cn], h5py.Dataset):
                        cols[cn] = np.array(obj[cn])
                sec = Section(name=nm, columns=cols, attrs=attrs)
            out.sections.append(sec)
    return out


def _h5_children_in_order(group) -> List[str]:
    """Best-effort enumerate of child names in creation order, with name fallback."""
    import h5py  # noqa: F401
    try:
        # h5py exposes get_link_creation_order via low-level API.
        gcpl = group.id.get_create_plist()
        crt_tracked = gcpl.get_link_creation_order() & 0x1   # H5P_CRT_ORDER_TRACKED bit
        if crt_tracked:
            # h5py iterates by name by default; use low-level link iteration.
            n = group.id.get_num_objs()
            order = []
            for i in range(n):
                order.append(group.id.get_objname_by_idx(i).decode()
                             if isinstance(group.id.get_objname_by_idx(i), bytes)
                             else group.id.get_objname_by_idx(i))
            return order
    except Exception:
        pass
    return list(group.keys())


def _h5_attr_to_py(value: Any) -> Any:
    """Convert an HDF5 attribute value to a Python scalar / list / bytes-string."""
    if isinstance(value, np.ndarray):
        if value.shape == (1,):
            return _h5_attr_to_py(value[0])
        return value.tolist()
    if isinstance(value, (bytes, np.bytes_)):
        try:
            return value.decode().rstrip('\x00').rstrip()
        except UnicodeDecodeError:
            return value
    if isinstance(value, np.generic):
        return value.item()
    return value


# ---------------------------------------------------------------------------
# Public reader
# ---------------------------------------------------------------------------

def load_lart(path: str) -> LartFile:
    """Read a LaRT output file in either FITS or HDF5 format.

    The returned :class:`LartFile` has ``.sections`` in HDU / group-creation
    order; each :class:`Section` carries either ``.data`` (image-style) or
    ``.columns`` (table-style) plus ``.attrs``.
    """
    fmt = detect_format(path)
    if fmt == 'fits':
        return _load_fits(path)
    return _load_hdf5(path)


# ---------------------------------------------------------------------------
# Converters
# ---------------------------------------------------------------------------

def _write_fits(lf: LartFile, path: str) -> None:
    from astropy.io import fits

    hdul = fits.HDUList()
    primary_filled = False

    def _attrs_to_header(attrs: Dict[str, Any]) -> fits.Header:
        h = fits.Header()
        for k, v in attrs.items():
            if v is None:
                continue
            try:
                if isinstance(v, (bytes, np.bytes_)):
                    v = v.decode() if isinstance(v, (bytes, np.bytes_)) else v
                h[str(k)[:8].upper()] = v
            except Exception:
                pass
        return h

    # Root attributes go on Primary HDU header.
    primary_header = _attrs_to_header(lf.root_attrs)

    for i, sec in enumerate(lf.sections):
        hdr = _attrs_to_header(sec.attrs)
        if sec.name and 'EXTNAME' not in hdr:
            hdr['EXTNAME'] = sec.name

        if sec.is_table:
            cols = []
            for cn, arr in sec.columns.items():
                cols.append(fits.Column(name=cn, array=arr,
                                        format=_numpy_to_tform(arr)))
            hdu = fits.BinTableHDU.from_columns(cols, header=hdr, name=sec.name)
            if not primary_filled:
                hdul.append(fits.PrimaryHDU(header=primary_header))
                primary_filled = True
            hdul.append(hdu)
        elif sec.data is not None:
            if not primary_filled:
                # Fold the first image-style section into the Primary HDU,
                # matching LaRT's own CFITSIO-driven write pattern.
                merged = fits.Header()
                merged.update(primary_header)
                merged.update(hdr)
                hdul.append(fits.PrimaryHDU(data=sec.data, header=merged))
                primary_filled = True
            else:
                hdul.append(fits.ImageHDU(data=sec.data, header=hdr, name=sec.name))
        # empty sections are dropped

    if not primary_filled:
        hdul.append(fits.PrimaryHDU(header=primary_header))

    if os.path.exists(path):
        os.remove(path)
    hdul.writeto(path)


def _numpy_to_tform(arr: np.ndarray) -> str:
    """Map a 1-D numpy dtype to a FITS BinTable TFORM code."""
    dt = arr.dtype
    if dt == np.float64:
        return 'D'
    if dt == np.float32:
        return 'E'
    if dt in (np.int32,):
        return 'J'
    if dt in (np.int64,):
        return 'K'
    if dt in (np.int16,):
        return 'I'
    if dt.kind in ('U', 'S'):
        n = max(1, arr.dtype.itemsize if dt.kind == 'S' else dt.itemsize // 4)
        return f'{n}A'
    return 'D'


def _write_hdf5(lf: LartFile, path: str) -> None:
    import h5py

    if os.path.exists(path):
        os.remove(path)
    with h5py.File(path, 'w', libver='latest', track_order=True) as f:
        # Root attributes
        for k, v in lf.root_attrs.items():
            _h5_set_attr(f, k, v)

        for sec in lf.sections:
            g = f.create_group(sec.name, track_order=True)
            for k, v in sec.attrs.items():
                _h5_set_attr(g, k, v)
            if sec.is_table:
                for cn, arr in sec.columns.items():
                    _h5_create_dataset(g, cn, arr)
            elif sec.data is not None:
                _h5_create_dataset(g, 'data', sec.data)


def _h5_create_dataset(parent, name: str, arr: np.ndarray) -> None:
    # Chunk + gzip for arrays large enough to benefit.
    kwargs: Dict[str, Any] = {}
    if arr.size > 4096:
        if arr.ndim == 1:
            chunks = (min(arr.shape[0], 4096),)
        else:
            chunks = tuple(min(s, 64) for s in arr.shape)
        kwargs.update(dict(chunks=chunks, compression='gzip', compression_opts=4))
    parent.create_dataset(name, data=arr, **kwargs)


def _h5_set_attr(obj, name: str, value: Any) -> None:
    if isinstance(value, str):
        # Store as fixed-length byte string for FITS round-trip stability.
        obj.attrs.create(name, value)
    else:
        obj.attrs[name] = value


def convert(src: str, dst: str) -> None:
    """Convert ``src`` to ``dst``, picking direction from file extensions."""
    src_fmt = detect_format(src)
    dst_fmt = detect_format(dst)
    lf = load_lart(src)
    if dst_fmt == 'fits':
        _write_fits(lf, dst)
    else:
        _write_hdf5(lf, dst)
    print(f'Converted {src} [{src_fmt}] -> {dst} [{dst_fmt}]  '
          f'({len(lf.sections)} sections)')


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def _usage() -> str:
    return ('Usage:\n'
            '  python lart_io.py info <file>\n'
            '  python lart_io.py convert <src> <dst>\n')


def _main(argv: List[str]) -> int:
    if len(argv) < 2:
        print(_usage(), file=sys.stderr)
        return 2
    cmd = argv[1]
    if cmd == 'info' and len(argv) == 3:
        print(load_lart(argv[2]).info())
        return 0
    if cmd == 'convert' and len(argv) == 4:
        convert(argv[2], argv[3])
        return 0
    print(_usage(), file=sys.stderr)
    return 2


if __name__ == '__main__':
    sys.exit(_main(sys.argv))
