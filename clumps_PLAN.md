# Plan: Spatially-Varying Clump Properties in LaRT Clumpy-Medium Mode

## Goal

Extend the LaRT clumpy-medium mode so that the following can vary with distance
from the box centre (or, more generally, per clump):

1. Clump radius `r_cl(r)`
2. Clump-internal gas density `n_H(r)` (and therefore `cl_rhokap`, `cl_voigt_a`,
   `cl_Dfreq` if temperature varies)
3. Clump number density `n_cl(r)` (i.e. the spatial volume filling factor)

The current implementation assumes all three are uniform in the spherical box.

---

## Working Mode

The implementing assistant proceeds through the phases below **autonomously,
without asking the user for per-step approval**. The user is consulted only
when a question arises that genuinely cannot be resolved without their
input — for example, an ambiguous physical convention not pinned down here,
a destructive or irreversible action, an unexpected conflict with existing
behaviour, or an open design choice not already decided in this document.
Routine progress, intermediate test results, and the eventual final summary
are reported to the user, but no confirmation is requested before each
file edit, build, or test run.

---

## Current State (verified for v2.10; v2.00 is identical)

- **Per-clump arrays already exist** (shared memory, dimension `N_clumps`):
  `cl_x/y/z`, `cl_vx/y/z`. Velocities are stored dimensionlessly as
  `v / cl_vtherm` and are already a function of position via `velocity_type`.
- **Scalar (uniform) physical properties** in `clump_mod.f90:26–32`:
  `cl_radius`, `cl_radius2`, `cl_rhokap`, `cl_voigt_a`, `cl_Dfreq`,
  `cl_vtherm`, `cl_temperature`.
- **Uniformity assumptions in code**:
  - RSA overlap test uses `(2*cl_radius)²` (`clump_mod.f90:229`).
  - RSA / CSR cell sizes use the single `cl_radius`
    (`clump_mod.f90:191, 382–383`).
  - `clump_cell_range` uses single `cl_radius` (`clump_mod.f90:440–452`).
  - `ray_sphere_isect` hard-codes `cl_radius2` (`clump_mod.f90:469`).
  - Opacity in raytrace: `kap = cl_rhokap * voigt(xfreq, cl_voigt_a)`
    (`raytrace_clump.f90:82`).
  - Column density and peel-off: same scalar lookups.
- **`N_clumps` derivation** (closed-form from uniform assumption):
  `f_vol = N * (r_cl/R)³`, `f_cov = (3/4) * N * (r_cl/R)²`.
- **Velocity field** is already position-dependent (per-clump), so velocity
  variation is orthogonal and unchanged by this plan.

---

## Design Decision

**Convert all scalar clump physical properties to per-clump arrays.**
A `'constant'` profile fills every entry with the same value, exactly
reproducing the current behaviour (modulo RNG noise) for backward compatibility.

Rationale: marginal memory cost (~50 bytes per clump × N_clumps), much simpler
than threading position-dependent function calls through every raytrace
hot-path query, and naturally supports random per-clump scatter as well as
purely radial profiles.

---

## Implementation Phases

### Phase 1 — Convert scalars to per-clump arrays (no behaviour change)

In `clump_mod.f90`:

- Replace these scalars with shared-memory arrays of dimension `N_clumps`:
  `cl_radius(:)`, `cl_radius2(:)`, `cl_rhokap(:)`, `cl_voigt_a(:)`,
  `cl_Dfreq(:)`, `cl_vtherm(:)`, `cl_temperature(:)`.
- Keep a scalar `cl_radius_max` for use as the RSA / CSR grid cell size.
- Update signatures:
  - `ray_sphere_isect(..., icl)` → uses `cl_radius2(icl)`
  - `clump_exit_dist(..., icl)` → uses `cl_radius(icl)`
  - `clump_cell_range(icl, ...)` → uses `cl_radius(icl)`
- RSA overlap test: `d² < (cl_radius(i) + cl_radius_trial)²`.
- RSA grid cell size: `2 * cl_radius_max` (reduces to current behaviour when
  uniform).
- CSR grid cell size: `2 * cl_radius_max / cgx` (same reduction).
- Inside `find_next_clump`, the `t_max` stopping logic is unchanged because
  per-cell `best_te` is still chosen from per-clump intersections.

Fill all per-clump arrays with the existing scalar values so that the input
parameters (`par%clump_radius`, `par%clump_tau0`, `par%temperature`) continue
to work unchanged.

**Validation**: the six reference cases in `examples/clump_sphere/` must agree
with current results to RNG noise level.

### Phase 2 — Add radial-profile input parameters

In `define.f90` `params_type`, add (defaults shown):

```fortran
character(len=32) :: clump_radius_profile  = 'constant'
character(len=32) :: clump_density_profile = 'constant'   ! n_H(r) inside clump
character(len=32) :: clump_number_profile  = 'constant'   ! n_cl(r), spatial
real(wp)          :: clump_radius_alpha    = 0.0_wp
real(wp)          :: clump_radius_r0       = 0.0_wp
real(wp)          :: clump_density_alpha   = 0.0_wp
real(wp)          :: clump_density_r0      = 0.0_wp
real(wp)          :: clump_number_alpha    = 0.0_wp
real(wp)          :: clump_number_r0       = 0.0_wp
character(len=:)  :: clump_profile_file    = ''           ! optional tabulated file
```

Built-in profile shapes (per axis): `'constant'`, `'powerlaw'`, `'gaussian'`,
`'exponential'`, `'file'`.

`'file'` mode reads a plain-text or FITS table with columns
`r, r_cl(r), n_H(r), T(r), n_cl(r)` and interpolates. This is the most
flexible single addition.

When all three profiles are `'constant'`, the existing scalar inputs
(`par%clump_radius`, `par%clump_tau0`, `par%temperature`) are honoured exactly
as today (full backward compatibility).

### Phase 3 — Non-uniform position sampling

`generate_clumps` currently samples uniformly inside the sphere by box
rejection. Replace with:

1. Build a fine 1-D radial grid (~1000 points) of `n_cl(r) * 4π r²`.
2. Compute its CDF and an inverse-CDF lookup table.
3. Per trial: sample `r ~ inverse-CDF`, `(θ, φ) ~` isotropic.
4. Evaluate `r_cl(r)` for the trial radius.
5. Run the RSA overlap test against the linked-list grid.

The RSA grid cell size must be `≥ 2 * cl_radius_max` so the 27-neighbour
overlap search remains complete. If the dynamic range of `r_cl(r)` is large
(say max/min > 10), the simple single-resolution grid becomes inefficient — a
hierarchical RSA grid is a future improvement, not part of the first cut.

### Phase 4 — Generalise `N_clumps` derivation

The closed-form formulas break for non-uniform profiles. Replace with numerical
quadrature in a helper `derive_clump_normalization`:

- **`f_vol` given**:
  `f_vol = ∫₀ᴿ n_cl(r) · (4/3)π r_cl(r)³ · 4π r² dr / ((4/3)π R³)` — solve for
  the normalisation constant `A` of `n_cl(r) = A · f̂(r)`. Then
  `N_clumps = ∫₀ᴿ n_cl(r) · 4π r² dr`.
- **`f_cov` given**: `f_cov ≈ ∫₀ᴿ n_cl(r) · π r_cl(r)² dr` along a radial
  sight line through the centre. (The exact LaRT convention will be preserved
  in the constant-profile limit; this generalisation must be documented.)
- **`N_clumps` given**: invert either of the above to fix `A`.

### Phase 5 — Per-clump property assignment

After positions are placed, for each clump compute
`r = sqrt(x² + y² + z²)` and evaluate:

- `cl_radius(icl) = r_cl(r)`, `cl_radius2(icl) = cl_radius(icl)²`
- `cl_temperature(icl) = T(r)`
- `cl_vtherm(icl)` from `T(r)` (line-dependent)
- `cl_Dfreq(icl)`, `cl_voigt_a(icl)` from `T(r)` via the standard formulas
- `cl_rhokap(icl) = n_H(r) * line%cross0 / cl_Dfreq(icl)`

**Velocity unit consistency**: bulk velocities are stored as
`cl_v* / cl_vtherm`. With temperature varying per clump, this dimensionless
storage must use `cl_vtherm(icl)` consistently in both `assign_clump_velocities`
and in `raytrace_clump.f90` where `u_los` is computed. Verify there is no
implicit assumption that `cl_vtherm` is global.

### Phase 6 — Output / diagnostics

- Extend the `save_clump_info` FITS table with per-clump columns:
  `cl_radius`, `cl_nH` (or `cl_rhokap`), `cl_temperature`.
- Log the realised `f_vol`, `f_cov`, mean `r_cl`, mean `tau_per_clump`, and
  total HI mass on rank 0 so the user can verify the requested normalisation.

### Phase 7 — Docs and validation

- Append a section to `CLUMPS_CHANGES.txt` and update
  `docs/LaRT_clump_description.tex` to describe the new profile inputs,
  normalisation, and any redefinitions of `f_cov`.
- Add 1–2 example inputs under `examples/clump_sphere/`, e.g.
  `clump_powerlaw_density.in`.
- Validation cases:
  1. Constant profile reproduces existing six reference cases (Phase 1
     regression).
  2. Analytic limit checks: e.g. `n_H(r) ∝ r⁻α` with constant `r_cl`, `n_cl`
     should give a known total tau as a function of α; cross-check against an
     equivalent Cartesian or AMR spherically-symmetric run with the same
     volume-averaged opacity.
  3. Spectrum comparison vs. a non-clumpy run with the same effective average
     opacity to confirm physical reasonableness.

---

## Open Design Questions

1. **Profile representation**: support both built-in shapes (powerlaw,
   gaussian, exponential) *and* file-input, or only file-input for the first
   cut? Recommendation: powerlaw + file for v1. => 두 경우 모두 고려해야 합니다.
2. **`f_cov` definition under non-uniformity**: line-of-sight integral
   (default, matches direct physical intuition) vs. a covering-factor weighted
   average. Decide and document explicitly.
3. **`tau0` input mode**: keep `clump_tau0` (re-interpreted as "tau at
   reference radius") vs. require `n_H(r)` directly. Possibly support both via
   a `clump_tau0_ref_r` parameter.
4. **CSR / RSA grid efficiency** for large dynamic range in `r_cl`. First cut
   uses a single resolution; consider a hierarchical grid only if profiling
   shows it matters. => 그렇게 하는 것이 올바른 방법일 듯 하군요.

---

## Phase Summary

| Phase | Scope | Behaviour change |
|-------|-------|------------------|
| 1 | Scalars → per-clump arrays | None (regression-checked) |
| 2 | New input parameters | None when defaults used |
| 3 | Non-uniform position sampling | Active when `clump_number_profile ≠ 'constant'` |
| 4 | Numerical N_clumps normalisation | Active for non-uniform profiles |
| 5 | Per-clump property assignment from profiles | Active for non-constant profiles |
| 6 | Output diagnostics | Additive |
| 7 | Docs + validation | — |

Phase 1 is the most invasive refactor but should be a no-op in behaviour;
once it lands cleanly, Phases 2–5 are localised additions.
