#!/usr/bin/env python3
"""
convert_ramses_to_generic_amr.py
================================

Convert RAMSES AMR + hydro output into the generic AMR format understood by
LaRT v2.00.

This script mirrors the record order and unit conventions documented in:
  - LaRT_v2.00/docs/RAMSES_data_structure.tex
  - LaRT_v2.00/docs/LaRT_AMR_description.tex
  - LaRT_v2.00/read_ramses_amr.f90

Output formats
--------------
1. Plain text (default)
     nleaf  boxlen
     x  y  z  level  nH[cm^-3]  T[K]  vx[km/s]  vy[km/s]  vz[km/s]
2. FITS binary table (.fits or .fits.gz)

The generic file stores total hydrogen number density nH, not neutral
hydrogen density nHI. LaRT computes the neutral fraction later using either
full neutrality or the CIE option.

If ``hydro_file_descriptor.txt`` is present in ``output_XXXXX/``, the script
uses it to infer whether variables 2--4 are velocities or momenta and whether
the thermodynamic variable is stored as total energy density or thermal
pressure. This matters for RAMSES variants such as the jellyfish RMHD runs.
"""

from __future__ import annotations

import argparse
import re
import struct
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np

try:
    from astropy.io import fits
except ImportError:
    fits = None


MASS_H_CGS = 1.6726e-24
KB_CGS = 1.381e-16
KPC_TO_CM = 3.0856775814913673e21
PC_TO_CM = 3.0856775814913673e18
AU_TO_CM = 1.495978707e13


@dataclass
class RamsesInfo:
    repository: Path
    snapnum: int
    ncpu: int
    ndim: int
    nx: int
    ny: int
    nz: int
    nlevelmax: int
    boxlen_code: float
    unit_l: float
    unit_d: float
    unit_t: float
    gamma: float

    @property
    def unit_v(self) -> float:
        return self.unit_l / self.unit_t

    @property
    def boxlen_cm(self) -> float:
        return self.boxlen_code * self.unit_l


@dataclass
class HydroLayout:
    velocity_layout: str = "momentum"
    thermo_mode: str = "energy"
    density_index: int = 1
    velocity_indices: Tuple[int, int, int] = (2, 3, 4)
    thermo_index: int = 5
    descriptor_found: bool = False


class FortranRecordReader:
    """Simple reader for sequential Fortran unformatted files."""

    def __init__(self, filename: Path, endian: Optional[str] = None):
        self.filename = Path(filename)
        self.handle = self.filename.open("rb")
        self.endian = endian or self._detect_endian()

    def close(self) -> None:
        self.handle.close()

    def __enter__(self) -> "FortranRecordReader":
        return self

    def __exit__(self, exc_type, exc, tb) -> None:
        self.close()

    def _detect_endian(self) -> str:
        pos = self.handle.tell()
        lead = self.handle.read(4)
        if len(lead) != 4:
            raise EOFError(f"Cannot read first record marker from {self.filename}")
        payload = self.handle.read(8)
        self.handle.seek(pos)

        for endian in ("<", ">"):
            recl = struct.unpack(f"{endian}i", lead)[0]
            if recl <= 0 or recl > 1024:
                continue
            if len(payload) < recl + 4:
                continue
            trail = payload[recl:recl + 4]
            if len(trail) != 4:
                continue
            if struct.unpack(f"{endian}i", trail)[0] == recl:
                return endian
        raise ValueError(f"Could not determine Fortran record endianness for {self.filename}")

    def read_record(self) -> bytes:
        lead = self.handle.read(4)
        if not lead:
            raise EOFError(f"Unexpected end of file in {self.filename}")
        if len(lead) != 4:
            raise EOFError(f"Truncated record marker in {self.filename}")
        recl = struct.unpack(f"{self.endian}i", lead)[0]
        if recl < 0:
            raise ValueError(f"Invalid record length {recl} in {self.filename}")
        payload = self.handle.read(recl)
        trail = self.handle.read(4)
        if len(payload) != recl or len(trail) != 4:
            raise EOFError(f"Truncated record in {self.filename}")
        tail_recl = struct.unpack(f"{self.endian}i", trail)[0]
        if tail_recl != recl:
            raise ValueError(
                f"Record length mismatch in {self.filename}: {recl} != {tail_recl}"
            )
        return payload

    def skip(self, nrec: int = 1) -> None:
        for _ in range(nrec):
            self.read_record()

    def read_scalar(self, dtype: np.dtype):
        arr = self.read_array(dtype)
        if arr.size != 1:
            raise ValueError(f"Expected one scalar from {self.filename}, got {arr.size}")
        return arr[0]

    def read_array(self, dtype: np.dtype) -> np.ndarray:
        dtype = np.dtype(dtype).newbyteorder(self.endian)
        payload = self.read_record()
        if len(payload) % dtype.itemsize != 0:
            raise ValueError(
                f"Record size {len(payload)} is incompatible with dtype {dtype} "
                f"in {self.filename}"
            )
        return np.frombuffer(payload, dtype=dtype).copy()


def normalize_repository_and_snapnum(path: str, snapnum: Optional[int]) -> Tuple[Path, int]:
    p = Path(path).expanduser().resolve()
    match = re.fullmatch(r"output_(\d{5})", p.name)
    if match:
        inferred = int(match.group(1))
        if snapnum is None:
            snapnum = inferred
        elif snapnum != inferred:
            raise ValueError(
                f"Input path implies snapnum={inferred}, but --snapnum={snapnum} was given."
            )
        p = p.parent
    if snapnum is None:
        raise ValueError("Specify --snapnum, or pass a path ending in output_00042.")
    return p, snapnum


def info_filename(repository: Path, snapnum: int) -> Path:
    return repository / f"output_{snapnum:05d}" / f"info_{snapnum:05d}.txt"


def amr_filename(repository: Path, snapnum: int, icpu: int) -> Path:
    return repository / f"output_{snapnum:05d}" / f"amr_{snapnum:05d}.out{icpu:05d}"


def hydro_filename(repository: Path, snapnum: int, icpu: int) -> Path:
    return repository / f"output_{snapnum:05d}" / f"hydro_{snapnum:05d}.out{icpu:05d}"


def descriptor_filename(repository: Path, snapnum: int) -> Path:
    return repository / f"output_{snapnum:05d}" / "hydro_file_descriptor.txt"


def read_info(repository: Path, snapnum: int) -> Dict[str, float]:
    filename = info_filename(repository, snapnum)
    values: Dict[str, float] = {}
    with filename.open("r", encoding="utf-8") as handle:
        for raw in handle:
            if "=" not in raw:
                continue
            key, value = raw.split("=", 1)
            key = key.strip()
            try:
                values[key] = float(value.strip().split()[0].replace("d", "e").replace("D", "e"))
            except ValueError:
                continue
    return values


def read_hydro_descriptor(repository: Path, snapnum: int) -> Dict[int, str]:
    filename = descriptor_filename(repository, snapnum)
    mapping: Dict[int, str] = {}
    if not filename.exists():
        return mapping
    for raw in filename.read_text(encoding="utf-8").splitlines():
        match = re.match(r"\s*variable\s*#\s*(\d+)\s*:\s*(.+?)\s*$", raw, re.IGNORECASE)
        if match:
            mapping[int(match.group(1))] = match.group(2).strip().lower()
    return mapping


def infer_hydro_layout(repository: Path, snapnum: int) -> HydroLayout:
    names = read_hydro_descriptor(repository, snapnum)
    if not names:
        return HydroLayout()

    layout = HydroLayout(descriptor_found=True)
    if all("velocity_" in names.get(i, "") for i in (2, 3, 4)):
        layout.velocity_layout = "velocity"
        layout.velocity_indices = (2, 3, 4)
    elif all(("momentum_" in names.get(i, "")) or ("rho_v" in names.get(i, "")) for i in (2, 3, 4)):
        layout.velocity_layout = "momentum"
        layout.velocity_indices = (2, 3, 4)

    for idx, name in names.items():
        if name == "density":
            layout.density_index = idx

    for idx, name in names.items():
        if name == "thermal_pressure":
            layout.thermo_mode = "pressure"
            layout.thermo_index = idx
            return layout

    for idx, name in names.items():
        if name in ("pressure", "gas_pressure", "thermal_energy_density", "total_energy", "energy"):
            layout.thermo_mode = "pressure" if "pressure" in name else "energy"
            layout.thermo_index = idx
            return layout
    return layout


def read_basic_amr_header(filename: Path) -> Tuple[str, Tuple[int, int, int, int, int, int, int, float]]:
    with FortranRecordReader(filename) as reader:
        endian = reader.endian
        ncpu = int(reader.read_scalar(np.int32))
        ndim = int(reader.read_scalar(np.int32))
        nx, ny, nz = reader.read_array(np.int32)
        nlevelmax = int(reader.read_scalar(np.int32))
        ngridmax = int(reader.read_scalar(np.int32))
        nboundary = int(reader.read_scalar(np.int32))
        ngrid_current = int(reader.read_scalar(np.int32))
        boxlen = float(reader.read_scalar(np.float64))
    return endian, (ncpu, ndim, int(nx), int(ny), int(nz), nlevelmax, nboundary, boxlen)


def build_ramses_info(repository: Path, snapnum: int) -> RamsesInfo:
    values = read_info(repository, snapnum)
    first_amr = amr_filename(repository, snapnum, 1)
    _, header = read_basic_amr_header(first_amr)
    ncpu_hdr, ndim, nx, ny, nz, nlevelmax, nboundary, boxlen_code = header
    if ndim != 3:
        raise ValueError(f"Only 3D RAMSES outputs are supported, got ndim={ndim}.")
    ncpu = int(values.get("ncpu", ncpu_hdr))
    return RamsesInfo(
        repository=repository,
        snapnum=snapnum,
        ncpu=ncpu,
        ndim=ndim,
        nx=nx,
        ny=ny,
        nz=nz,
        nlevelmax=nlevelmax,
        boxlen_code=float(values.get("boxlen", boxlen_code)),
        unit_l=float(values.get("unit_l", 1.0)),
        unit_d=float(values.get("unit_d", 1.0)),
        unit_t=float(values.get("unit_t", 1.0)),
        gamma=float(values.get("gamma", 5.0 / 3.0)),
    )


def read_amr_layout(reader: FortranRecordReader, ncpu: int, nlevelmax: int, nboundary: int) -> np.ndarray:
    reader.skip(13)
    ngridlevel = reader.read_array(np.int32).reshape((ncpu, nlevelmax), order="F")
    ngridfile = np.zeros((ncpu + nboundary, nlevelmax), dtype=np.int32)
    ngridfile[:ncpu, :] = ngridlevel
    reader.skip(1)
    if nboundary > 0:
        reader.skip(2)
        ngridbound = reader.read_array(np.int32).reshape((nboundary, nlevelmax), order="F")
        ngridfile[ncpu:, :] = ngridbound
    # rd_amr.pro shows six additional header records here:
    #   headf/tailf/numbf/used_mem/used_mem_tot
    #   ordering
    #   bound_key (or five bisection records in older layouts)
    #   son
    #   flag1
    #   cpu_map
    reader.skip(6)
    return ngridfile


def read_hydro_header(reader: FortranRecordReader) -> Tuple[int, float]:
    reader.skip(1)
    nvar = int(reader.read_scalar(np.int32))
    reader.skip(3)
    gamma = float(reader.read_scalar(np.float64))
    return nvar, gamma


def oct_offsets(ndim: int, ilevel: int) -> np.ndarray:
    twotondim = 2 ** ndim
    dx = 0.5 ** ilevel
    xc = np.zeros((twotondim, ndim), dtype=np.float64)
    for ind in range(twotondim):
        if ndim >= 3:
            xc[ind, 2] = (ind // 4)
        if ndim >= 2:
            xc[ind, 1] = ((ind // 2) % 2)
        xc[ind, 0] = (ind % 2)
        xc[ind, :] = (xc[ind, :] - 0.5) * dx
    return xc


def convert_positions_unit(values_cm: np.ndarray, unit: str) -> Tuple[np.ndarray, float]:
    unit = unit.lower()
    if unit == "cm":
        return values_cm, 1.0
    if unit == "kpc":
        return values_cm / KPC_TO_CM, KPC_TO_CM
    if unit == "pc":
        return values_cm / PC_TO_CM, PC_TO_CM
    if unit == "au":
        return values_cm / AU_TO_CM, AU_TO_CM
    raise ValueError(f"Unsupported output unit: {unit}")


def convert_ramses_snapshot(
    repository: Path,
    snapnum: int,
    hydro_precision: int = 8,
    mu: float = 1.22,
    density_index: Optional[int] = None,
    velocity_indices: Optional[Tuple[int, int, int]] = None,
    thermo_index: Optional[int] = None,
    velocity_layout: str = "auto",
    thermo_mode: str = "auto",
) -> Tuple[RamsesInfo, HydroLayout, Dict[str, np.ndarray]]:
    info = build_ramses_info(repository, snapnum)
    inferred = infer_hydro_layout(repository, snapnum)

    layout = HydroLayout(
        velocity_layout=inferred.velocity_layout,
        thermo_mode=inferred.thermo_mode,
        density_index=inferred.density_index,
        velocity_indices=inferred.velocity_indices,
        thermo_index=inferred.thermo_index,
        descriptor_found=inferred.descriptor_found,
    )
    if density_index is not None:
        layout.density_index = density_index
    if velocity_indices is not None:
        layout.velocity_indices = velocity_indices
    if thermo_index is not None:
        layout.thermo_index = thermo_index
    if velocity_layout != "auto":
        layout.velocity_layout = velocity_layout
    if thermo_mode != "auto":
        layout.thermo_mode = thermo_mode

    x_list: List[np.ndarray] = []
    y_list: List[np.ndarray] = []
    z_list: List[np.ndarray] = []
    level_list: List[np.ndarray] = []
    nH_list: List[np.ndarray] = []
    T_list: List[np.ndarray] = []
    vx_list: List[np.ndarray] = []
    vy_list: List[np.ndarray] = []
    vz_list: List[np.ndarray] = []

    xbound = np.array([info.nx // 2, info.ny // 2, info.nz // 2], dtype=np.float64)
    twotondim = 2 ** info.ndim

    for icpu in range(1, info.ncpu + 1):
        amr_path = amr_filename(repository, snapnum, icpu)
        hyd_path = hydro_filename(repository, snapnum, icpu)

        with FortranRecordReader(amr_path) as amr_reader, FortranRecordReader(hyd_path) as hyd_reader:
            ncpu_f = int(amr_reader.read_scalar(np.int32))
            ndim = int(amr_reader.read_scalar(np.int32))
            nx, ny, nz = amr_reader.read_array(np.int32)
            nlevelmax = int(amr_reader.read_scalar(np.int32))
            amr_reader.skip(1)
            nboundary = int(amr_reader.read_scalar(np.int32))
            amr_reader.skip(2)

            if ndim != info.ndim or ncpu_f != info.ncpu:
                raise ValueError(
                    f"Inconsistent AMR header in {amr_path}: ndim={ndim}, ncpu={ncpu_f}"
                )
            if (int(nx), int(ny), int(nz), nlevelmax) != (info.nx, info.ny, info.nz, info.nlevelmax):
                raise ValueError(f"Inconsistent AMR dimensions in {amr_path}")

            ngridfile = read_amr_layout(amr_reader, ncpu_f, nlevelmax, nboundary)
            nvar, gamma_h = read_hydro_header(hyd_reader)
            gamma = gamma_h if gamma_h > 0 else info.gamma

            required_indices = [layout.density_index, *layout.velocity_indices]
            if layout.thermo_mode != "constant":
                required_indices.append(layout.thermo_index)
            max_index = max(required_indices)
            if max_index > nvar:
                raise ValueError(
                    f"{hyd_path} stores nvar={nvar}, but index {max_index} was requested."
                )

            for ilevel in range(1, nlevelmax + 1):
                ngrida = int(ngridfile[icpu - 1, ilevel - 1])
                xc = oct_offsets(info.ndim, ilevel)

                if ngrida > 0:
                    xg = np.empty((ngrida, info.ndim), dtype=np.float64)
                    son = np.empty((ngrida, twotondim), dtype=np.int32)
                    var = np.empty((ngrida, twotondim, nvar), dtype=np.float64)
                else:
                    xg = son = var = None

                for j in range(1, ncpu_f + nboundary + 1):
                    ngrid_dom = int(ngridfile[j - 1, ilevel - 1])
                    if ngrid_dom > 0:
                        amr_reader.skip(3)
                        for idim in range(info.ndim):
                            if j == icpu:
                                xg[:, idim] = amr_reader.read_array(np.float64)
                            else:
                                amr_reader.skip(1)
                        amr_reader.skip(1 + 2 * info.ndim)
                        for ind in range(twotondim):
                            if j == icpu:
                                son[:, ind] = amr_reader.read_array(np.int32)
                            else:
                                amr_reader.skip(1)
                        amr_reader.skip(2 * twotondim)

                    hyd_reader.skip(2)
                    if ngrid_dom > 0:
                        for ind in range(twotondim):
                            for ivar in range(nvar):
                                if j == icpu:
                                    if hydro_precision == 4:
                                        var[:, ind, ivar] = hyd_reader.read_array(np.float32).astype(np.float64)
                                    else:
                                        var[:, ind, ivar] = hyd_reader.read_array(np.float64)
                                else:
                                    hyd_reader.skip(1)

                if ngrida <= 0:
                    continue

                for ind in range(twotondim):
                    mask = son[:, ind] == 0
                    if not np.any(mask):
                        continue

                    x_code = (xg[:, 0] + xc[ind, 0] - xbound[0]) / info.nx + 0.5
                    y_code = (xg[:, 1] + xc[ind, 1] - xbound[1]) / info.ny + 0.5
                    z_code = (xg[:, 2] + xc[ind, 2] - xbound[2]) / info.nz + 0.5

                    dens_code = var[:, ind, layout.density_index - 1]
                    rho_cgs = dens_code * info.unit_d
                    nH = rho_cgs / MASS_H_CGS

                    if layout.velocity_layout == "momentum":
                        with np.errstate(divide="ignore", invalid="ignore"):
                            vx_cms = np.where(
                                dens_code > 0.0,
                                var[:, ind, layout.velocity_indices[0] - 1] / dens_code * info.unit_v,
                                0.0,
                            )
                            vy_cms = np.where(
                                dens_code > 0.0,
                                var[:, ind, layout.velocity_indices[1] - 1] / dens_code * info.unit_v,
                                0.0,
                            )
                            vz_cms = np.where(
                                dens_code > 0.0,
                                var[:, ind, layout.velocity_indices[2] - 1] / dens_code * info.unit_v,
                                0.0,
                            )
                        kinetic_code = 0.5 * dens_code * (
                            np.where(dens_code > 0.0, var[:, ind, layout.velocity_indices[0] - 1] / np.maximum(dens_code, 1.0e-40), 0.0) ** 2
                            + np.where(dens_code > 0.0, var[:, ind, layout.velocity_indices[1] - 1] / np.maximum(dens_code, 1.0e-40), 0.0) ** 2
                            + np.where(dens_code > 0.0, var[:, ind, layout.velocity_indices[2] - 1] / np.maximum(dens_code, 1.0e-40), 0.0) ** 2
                        )
                    elif layout.velocity_layout == "velocity":
                        vx_cms = var[:, ind, layout.velocity_indices[0] - 1] * info.unit_v
                        vy_cms = var[:, ind, layout.velocity_indices[1] - 1] * info.unit_v
                        vz_cms = var[:, ind, layout.velocity_indices[2] - 1] * info.unit_v
                        kinetic_code = 0.5 * dens_code * (
                            var[:, ind, layout.velocity_indices[0] - 1] ** 2
                            + var[:, ind, layout.velocity_indices[1] - 1] ** 2
                            + var[:, ind, layout.velocity_indices[2] - 1] ** 2
                        )
                    else:
                        raise ValueError(f"Unsupported velocity layout: {layout.velocity_layout}")

                    if layout.thermo_mode == "energy":
                        eint_code = np.where(
                            dens_code > 0.0,
                            (var[:, ind, layout.thermo_index - 1] - kinetic_code) / np.maximum(dens_code, 1.0e-40),
                            0.0,
                        )
                        eint_code = np.maximum(eint_code, 0.0)
                        temp = (gamma - 1.0) * eint_code * info.unit_v ** 2 * mu * MASS_H_CGS / KB_CGS
                        temp = np.maximum(temp, 10.0)
                    elif layout.thermo_mode == "pressure":
                        pressure_cgs = var[:, ind, layout.thermo_index - 1] * info.unit_d * info.unit_v ** 2
                        temp = np.where(
                            rho_cgs > 0.0,
                            pressure_cgs * mu * MASS_H_CGS / np.maximum(rho_cgs, 1.0e-40) / KB_CGS,
                            10.0,
                        )
                        temp = np.maximum(temp, 10.0)
                    else:
                        temp = np.full_like(nH, 1.0e4)

                    x_list.append(x_code[mask] * info.boxlen_cm)
                    y_list.append(y_code[mask] * info.boxlen_cm)
                    z_list.append(z_code[mask] * info.boxlen_cm)
                    level_list.append(np.full(np.count_nonzero(mask), ilevel, dtype=np.int32))
                    nH_list.append(nH[mask])
                    T_list.append(temp[mask])
                    vx_list.append(vx_cms[mask] / 1.0e5)
                    vy_list.append(vy_cms[mask] / 1.0e5)
                    vz_list.append(vz_cms[mask] / 1.0e5)

    result = {
        "x_cm": np.concatenate(x_list) if x_list else np.array([], dtype=np.float64),
        "y_cm": np.concatenate(y_list) if y_list else np.array([], dtype=np.float64),
        "z_cm": np.concatenate(z_list) if z_list else np.array([], dtype=np.float64),
        "level": np.concatenate(level_list) if level_list else np.array([], dtype=np.int32),
        "nH": np.concatenate(nH_list) if nH_list else np.array([], dtype=np.float64),
        "T": np.concatenate(T_list) if T_list else np.array([], dtype=np.float64),
        "vx": np.concatenate(vx_list) if vx_list else np.array([], dtype=np.float64),
        "vy": np.concatenate(vy_list) if vy_list else np.array([], dtype=np.float64),
        "vz": np.concatenate(vz_list) if vz_list else np.array([], dtype=np.float64),
    }
    return info, layout, result


def write_generic_text(
    filename: Path,
    x: np.ndarray,
    y: np.ndarray,
    z: np.ndarray,
    level: np.ndarray,
    nH: np.ndarray,
    T: np.ndarray,
    vx: np.ndarray,
    vy: np.ndarray,
    vz: np.ndarray,
    boxlen: float,
) -> None:
    with filename.open("w", encoding="utf-8") as handle:
        handle.write(f"{len(x)}  {boxlen:.12e}\n")
        for values in zip(x, y, z, level, nH, T, vx, vy, vz):
            handle.write(
                "{:.12e}  {:.12e}  {:.12e}  {:d}  {:.12e}  {:.12e}  {:.12e}  {:.12e}  {:.12e}\n".format(
                    values[0], values[1], values[2], int(values[3]),
                    values[4], values[5], values[6], values[7], values[8]
                )
            )


def _auto_bitpix_col(arr: np.ndarray) -> int:
    """Return -32 (float32) or -64 (float64) based on dynamic range of arr.

    Mirrors the Fortran fitsio_mod auto-detect logic:
    max(|data|) / min(|data[nonzero]|) < 1e6  →  float32 sufficient.
    """
    a = np.asarray(arr, dtype=np.float64)
    max_abs = float(np.max(np.abs(a)))
    if max_abs == 0.0:
        return -32
    nz = np.abs(a[a != 0.0])
    return -32 if max_abs / float(np.min(nz)) < 1.0e6 else -64


def write_generic_fits(
    filename: Path,
    x: np.ndarray,
    y: np.ndarray,
    z: np.ndarray,
    level: np.ndarray,
    nH: np.ndarray,
    T: np.ndarray,
    vx: np.ndarray,
    vy: np.ndarray,
    vz: np.ndarray,
    boxlen: float,
    unit_l_cgs: float,
    bitpix: int = 0,
) -> None:
    if fits is None:
        raise ImportError("Writing FITS output requires astropy.")

    float_cols = {
        "x": x, "y": y, "z": z,
        "gasDen": nH, "T": T, "vx": vx, "vy": vy, "vz": vz,
    }
    cols = []
    for name, arr in float_cols.items():
        bp = _auto_bitpix_col(arr) if bitpix == 0 else bitpix
        fmt = "1E" if bp == -32 else "1D"
        col_arr = arr.astype(np.float32 if bp == -32 else np.float64)
        cols.append(fits.Column(name=name, format=fmt, array=col_arr))
        if name == "z":
            cols.append(fits.Column(name="level", format="1J",
                                    array=level.astype(np.int32)))

    primary = fits.PrimaryHDU()
    table   = fits.BinTableHDU.from_columns(cols, name="AMRGRID")
    table.header["CREATOR"]  = ("convert_ramses_to_generic_amr.py", "Created by")
    table.header["DATE"]     = (datetime.now().strftime("%Y-%m-%dT%H:%M"),
                                "Creation date-time (YYYY-MM-DDTHH:MM)")
    table.header["BOXLEN"]   = (float(boxlen),    "Simulation box length")
    table.header["UNITLCGS"] = (float(unit_l_cgs),"RAMSES unit_l in cm")
    table.header["ORIGINX"]  = (0.0, "Box origin x")
    table.header["ORIGINY"]  = (0.0, "Box origin y")
    table.header["ORIGINZ"]  = (0.0, "Box origin z")
    table.header["NLEAF"]    = (len(x), "Number of AMR leaf cells")
    fits.HDUList([primary, table]).writeto(filename, overwrite=True)


def output_is_fits(filename: Path) -> bool:
    name = filename.name.lower()
    return name.endswith(".fits") or name.endswith(".fits.gz")


def suggest_lart_distance(output_unit: str) -> str:
    output_unit = output_unit.lower()
    if output_unit == "cm":
        return "par%distance2cm = 1.0"
    return f"par%distance_unit = '{output_unit}'"


def write_lart_input(
    filename: Path,
    amr_filename: Path,
    output_unit: str,
    boxlen: float,
    nphotons: float = 1.0e6,
    taumax: float = 1.0e4,
) -> None:
    center = 0.5 * boxlen
    content = (
        "&parameters\n"
        " par%use_amr_grid  = .true.\n"
        " par%amr_type      = 'generic'\n"
        f" par%amr_file      = '{amr_filename.name}'\n"
        f" {suggest_lart_distance(output_unit)}\n"
        f" par%no_photons    = {nphotons:.0e}\n"
        f" par%taumax        = {taumax:.1e}\n"
        " par%DGR           = 0.0\n"
        " par%spectral_type = 'voigt'\n"
        " par%source_geometry = 'point'\n"
        f" par%xs_point      = {center:.6f}\n"
        f" par%ys_point      = {center:.6f}\n"
        f" par%zs_point      = {center:.6f}\n"
        " par%nxfreq        = 121\n"
        " par%nxim          = 128\n"
        " par%nyim          = 128\n"
        " par%nprint        = 100000\n"
        " par%out_merge     = .false.\n"
        "/\n"
    )
    filename.write_text(content, encoding="utf-8")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Convert a RAMSES snapshot to LaRT generic AMR format."
    )
    parser.add_argument(
        "repository",
        help="RAMSES repository root, or a specific output_00042 directory.",
    )
    parser.add_argument(
        "-s", "--snapnum",
        type=int,
        help="Snapshot number. Optional if repository ends with output_00042.",
    )
    parser.add_argument(
        "-o", "--output",
        required=True,
        help="Output generic AMR file (.dat, .fits, or .fits.gz). Use .fits.gz for compressed FITS.",
    )
    parser.add_argument(
        "--output-unit",
        choices=("cm", "kpc", "pc", "au"),
        default="kpc",
        help="Position unit written to the generic file (default: kpc).",
    )
    parser.add_argument(
        "--hydro-precision",
        type=int,
        choices=(4, 8),
        default=8,
        help="Precision used in hydro_*.out files (default: 8).",
    )
    parser.add_argument(
        "--mu",
        type=float,
        default=1.22,
        help="Mean molecular weight used for T from var5 (default: 1.22).",
    )
    parser.add_argument(
        "--density-index",
        type=int,
        default=None,
        help="1-based RAMSES hydro variable index for density. Default: auto from descriptor, else 1.",
    )
    parser.add_argument(
        "--velocity-indices",
        type=int,
        nargs=3,
        metavar=("IVX", "IVY", "IVZ"),
        default=None,
        help="1-based hydro indices for velocity or momentum components. Default: auto from descriptor, else 2 3 4.",
    )
    parser.add_argument(
        "--thermo-index",
        type=int,
        default=None,
        help="1-based hydro variable index for thermodynamic variable. Default: auto from descriptor, else 5.",
    )
    parser.add_argument(
        "--velocity-layout",
        choices=("auto", "momentum", "velocity"),
        default="auto",
        help="Interpret hydro velocity slots as rho*v or v. Default: auto from descriptor.",
    )
    parser.add_argument(
        "--thermo-mode",
        choices=("auto", "energy", "pressure", "constant"),
        default="auto",
        help="Thermodynamic variable interpretation. Default: auto from descriptor.",
    )
    parser.add_argument(
        "--write-input",
        action="store_true",
        help="Also write a matching LaRT .in file next to the output file.",
    )
    parser.add_argument(
        "--input-file",
        help="Explicit output path for the LaRT .in file. Implies --write-input.",
    )
    parser.add_argument(
        "--nphotons",
        type=float,
        default=1.0e6,
        help="Photon count written to the optional .in file (default: 1e6).",
    )
    parser.add_argument(
        "--taumax",
        type=float,
        default=1.0e4,
        help="taumax written to the optional .in file (default: 1e4).",
    )
    parser.add_argument(
        "--bitpix",
        type=int,
        choices=(0, -32, -64),
        default=0,
        help="FITS column precision: 0=auto per column (default), -32=float32, -64=float64.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    repository, snapnum = normalize_repository_and_snapnum(args.repository, args.snapnum)
    output = Path(args.output).expanduser().resolve()

    velocity_indices = tuple(args.velocity_indices) if args.velocity_indices is not None else None

    info, layout, data = convert_ramses_snapshot(
        repository=repository,
        snapnum=snapnum,
        hydro_precision=args.hydro_precision,
        mu=args.mu,
        density_index=args.density_index,
        velocity_indices=velocity_indices,
        thermo_index=args.thermo_index,
        velocity_layout=args.velocity_layout,
        thermo_mode=args.thermo_mode,
    )

    x_out, _ = convert_positions_unit(data["x_cm"], args.output_unit)
    y_out, _ = convert_positions_unit(data["y_cm"], args.output_unit)
    z_out, _ = convert_positions_unit(data["z_cm"], args.output_unit)
    boxlen_out, _ = convert_positions_unit(np.array([info.boxlen_cm]), args.output_unit)
    boxlen_out = float(boxlen_out[0])

    if output_is_fits(output):
        write_generic_fits(
            output, x_out, y_out, z_out, data["level"], data["nH"], data["T"],
            data["vx"], data["vy"], data["vz"], boxlen_out, info.unit_l,
            bitpix=args.bitpix,
        )
    else:
        write_generic_text(
            output, x_out, y_out, z_out, data["level"], data["nH"], data["T"],
            data["vx"], data["vy"], data["vz"], boxlen_out
        )

    print(f"Converted snapshot output_{snapnum:05d}")
    print(f"  repository   : {repository}")
    print(f"  output file  : {output}")
    print(f"  nleaf        : {len(data['level'])}")
    print(f"  boxlen       : {boxlen_out:.6e} {args.output_unit}")
    print(f"  velocity     : {layout.velocity_layout} (vars {layout.velocity_indices})")
    print(f"  thermo       : {layout.thermo_mode} (var {layout.thermo_index})")
    print(f"  descriptor   : {'found' if layout.descriptor_found else 'not found'}")
    if len(data["level"]) > 0:
        print(f"  level range  : {data['level'].min()} .. {data['level'].max()}")
        print(f"  nH range     : {data['nH'].min():.6e} .. {data['nH'].max():.6e} cm^-3")
        print(f"  T range      : {data['T'].min():.6e} .. {data['T'].max():.6e} K")
    else:
        print("  level range  : no leaf cells found")
    print("Suggested LaRT settings:")
    print("  par%use_amr_grid = .true.")
    print("  par%amr_type     = 'generic'")
    print(f"  par%amr_file     = '{output}'")
    print(f"  {suggest_lart_distance(args.output_unit)}")

    write_input = args.write_input or (args.input_file is not None)
    if write_input:
        if args.input_file is not None:
            input_path = Path(args.input_file).expanduser().resolve()
        else:
            input_name = output.name
            if input_name.lower().endswith(".fits.gz"):
                input_name = input_name[:-8]
            else:
                input_name = output.stem
            input_path = output.with_name(f"{input_name}.in")
        write_lart_input(
            input_path,
            output,
            args.output_unit,
            boxlen_out,
            nphotons=args.nphotons,
            taumax=args.taumax,
        )
        print(f"  input file   : {input_path}")


if __name__ == "__main__":
    main()
