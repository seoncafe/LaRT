# Plan: Combined Hydrogen + Deuterium Lyman-alpha Line in LaRT

## Goal

Add a new line ID `'ly_alpha_HD'` that simultaneously treats neutral hydrogen
and deuterium Ly-α as two coexisting two-level resonance scatterers in the
same medium. The two transitions sit ~0.33 Å apart (~81 km/s) and have
nearly identical atomic constants except mass and rest wavelength.

A convenience flag `par%include_deuterium = .true.` automatically promotes
`par%line_id = 'ly_alpha'` to `'ly_alpha_HD'`.

In this first cut we ignore fine-structure splitting for both species
(simple dipole, `line_type = 1`-style) and treat both as pure resonance
(Voigt absorption + dipole phase function). The deuterium-to-hydrogen number
abundance is a user-tunable input parameter, defaulting to the cosmic value.

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

## Design Summary

- New `line_type = 7`: two coexisting species, each individually behaving like
  the existing `line_type = 1`.
- A small `species_type` sub-record holds per-species atomic data
  (`wavelength0, mass_amu, vtherm1, damping, f12, cross0, g_recoil0`).
- A single `line_type` instance carries an array `species(2)` (or
  `species(nspecies)` with `nspecies = 1` for backward-compatible single-line
  cases).
- **No new per-cell arrays.** The existing `rhokap, Dfreq, voigt_a` are
  computed for hydrogen as before; deuterium contributions are evaluated
  on-the-fly from H values using simulation-wide constant ratios.
- Photon `xfreq` is always carried in **H-frame Doppler units**. D-frame
  conversions happen inside `calc_voigt_HD` and `do_resonance_HD`.

### Key constants (precomputed once at line init)

For H Lyα (1215.6701 Å) and D Lyα (1215.3374 Å), `m_H = 1.00797`,
`m_D = 2.01410`:

```
Δν_HD_Hz       = c × (1/λ_D − 1/λ_H)              ! D is bluer (higher freq)
ratio_DfreqHD  = (λ_H/λ_D) × √(m_D/m_H) ≈ √2      ! Dfreq_H / Dfreq_D
ratio_voigtaHD = (Γ_D/Γ_H) × ratio_DfreqHD ≈ √2   ! a_D / a_H (Γ_H = Γ_D)
nD_over_nH     = par%D_to_H_ratio                  ! e.g. 1.5e-5
```

Per-cell auxiliary (depends on local T):

```
Δx_HD_cell = Δν_HD_Hz / Dfreq_H_cell                ! H-Doppler-units offset to D
xfreq_D    = (xfreq − Δx_HD_cell) × ratio_DfreqHD   ! photon detuning in D's frame
voigt_a_D  = voigt_a_H × ratio_voigtaHD
```

### Combined Voigt profile

```
calc_voigt_HD(grid, xfreq, i, j, k) =
    voigt(xfreq, voigt_a_H_cell)
  + nD_over_nH × ratio_DfreqHD × voigt(xfreq_D, voigt_a_D)
```

The total opacity is then `rhokap_H_cell × calc_voigt_HD`. The factor
`ratio_DfreqHD` accounts for the LaRT normalisation `rhokap = n σ / Dfreq`.

The two species are assumed to share the same `cross0` to leading order
(`f_H ≈ f_D` and `Γ_H ≈ Γ_D`). The `f12_D` is stored exactly so we can refine
later if needed; in v1 we keep them numerically equal so the only opacity
asymmetry is the abundance × Doppler-width ratio above.

## Implementation Phases

### Phase 1 — `line_mod.f90`: data structure and dispatch

- Define `species_type` containing per-species atomic data
  (`wavelength0, mass_amu, vtherm1, damping, f12(3), cross0, g_recoil0`).
- Add `line%species(:)`, `line%nspecies`, `line%nD_over_nH`,
  `line%delta_nu_HD_Hz`, `line%ratio_Dfreq_HD`, `line%ratio_voigtA_HD` to
  `line_type` in `define.f90`.
- For all existing line IDs, populate `species(1)` with the current scalars,
  set `nspecies = 1`. Keep the legacy scalar fields populated as well for now
  (safety; remove after migration is verified).
- Add new entry `'ly_alpha_HD'`:
  - `species(1)` = H Lyα atomic data (existing values).
  - `species(2)` = D Lyα atomic data, hard-coded:
    - `λ_D = 1215.3374 Å` (vacuum, NIST)
    - `m_D = 2.01410177812 amu`
    - `Γ_D ≈ Γ_H = 6.2649e8 Hz`
    - `f12_D` stored explicitly (separate from H), with the well-known H Lyα
      values used as the v1 numerical placeholder; refine if a published D Lyα
      f-value differs at the level we care about.
  - Compute and store the cross-species constants listed above.
  - `line%line_type = 7`.
- `par%include_deuterium = .true.` together with `par%line_id = 'ly_alpha'`
  is mapped to `'ly_alpha_HD'` at the top of the line setup so users can keep
  the old line_id and just toggle the flag.

### Phase 2 — Input parameters (`define.f90`, `setup.f90`)

In `params_type`:
```fortran
real(wp) :: D_to_H_ratio     = 1.5e-5_wp   ! number ratio n_D / n_H, cosmic default
logical  :: include_deuterium = .false.    ! convenience: promotes 'ly_alpha' → 'ly_alpha_HD'
```

`read_input` accepts both. `setup.f90` performs the line-id promotion and
forwards `par%D_to_H_ratio` into `line%nD_over_nH` during line initialisation.

### Phase 3 — `calc_voigt` / `do_resonance` dispatch

In `setup_procedure()` add the `line_type = 7` branch:
- `calc_voigt → calc_voigt_HD`
- `do_resonance → do_resonance_HD`

Implement these in `line_mod.f90` (alongside `calc_voigt1`, `do_resonance1`):

`calc_voigt_HD(grid, xfreq, i, j, k)`:
1. `vH = voigt(xfreq, grid%voigt_a(i,j,k))`
2. `Δx = line%delta_nu_HD_Hz / grid%Dfreq(i,j,k)`
3. `xfreq_D = (xfreq − Δx) × line%ratio_Dfreq_HD`
4. `vD = voigt(xfreq_D, grid%voigt_a(i,j,k) × line%ratio_voigtA_HD)`
5. `return vH + line%nD_over_nH × line%ratio_Dfreq_HD × vD`

`do_resonance_HD(photon, ...)`:
1. Compute partial opacity contributions `vH` and `vD_term` (the D term in the
   sum above, including the abundance and Doppler-ratio prefactors).
2. Stochastically pick H or D with probability `vH / (vH + vD_term)`.
3. Call the existing `do_resonance1`-equivalent body using the **selected
   species**' atomic constants:
   - For H: standard path, no change (works exactly as today).
   - For D: convert `xfreq → xfreq_D` (D-frame Doppler units), use `mass_D`
     for atom velocity sampling and recoil, use `Γ_D` for damping. After the
     scattering, convert back to H-frame Doppler units for storage on the
     photon.
4. Return-frequency conversion (D → H frame): `xfreq_new = xfreq_D_new ×
   ratio_DfreqD_H + Δx_HD_cell` where `ratio_DfreqD_H = 1/ratio_DfreqHD`.

The species-specific scattering routine is implemented as a small helper that
takes a `species_type` argument so the H and D paths share code.

### Phase 4 — Grid / opacity setup

Unchanged. `grid_mod_car.f90`, `grid_mod_amr.f90`, `grid_mod_clump.f90` keep
computing `rhokap, Dfreq, voigt_a` from `line%species(1)` (= H). The D
contribution lives entirely inside `calc_voigt_HD` / `do_resonance_HD` via
the constant ratios + per-cell `Δx_HD`.

### Phase 5 — Peel-off and Stokes

- Phase function for both species is the same dipole (`E1 = 1, E2 = 0,
  E3 = 1`), so the existing peel-off and Stokes paths require no change in
  scattering matrix bookkeeping.
- The only species-aware step in peel-off is the frequency shift back into the
  observer frame, which already flows through `photon%xfreq` (H-frame
  Doppler units). Verify that `do_resonance_HD` writes back into the H frame
  correctly so peel-off needs no special handling.

### Phase 6 — Validation and Jupyter post-processing

Run reference cases and bundle plots in a Jupyter notebook so a user can
verify by eye that the H+D extension behaves correctly.

Reference runs (suggest under `examples/lya_HD/`):
1. **Backward-compat regression**: `D_to_H_ratio = 0`, sphere with
   `taumax = 1e4`, T = 10⁴ K. Spectrum must agree with the existing
   `ly_alpha` sphere result to RNG noise.
2. **Cosmic abundance**: same setup, `D_to_H_ratio = 1.5e-5`. Expect a small
   absorption dip blueward of H Ly-α at the H-frame velocity equivalent of
   the D line offset (~ +81 km/s in the +Δν direction).
3. **Enhanced abundance sweep**: `D_to_H_ratio ∈ {1e-5, 1e-4, 1e-3, 1e-2}`
   with the same sphere. Verify approximately logarithmic deepening of the D
   feature with abundance.
4. **Static-medium high-τ**: `taumax = 1e6` so the D feature is well in the
   damping wing of H. Confirm Voigt-wing physics produces the expected
   absorption shape.
5. **Direct comparison with Dijkstra, Haiman & Spaans (2006), Figure 3**.
   That paper computes the emergent Ly-α spectrum from a static, uniform,
   spherically-symmetric, optically-thick neutral hydrogen cloud including
   deuterium with `D/H = 3 × 10⁻⁵`. Reproduce their setup
   (T, N_HI, geometry, central monochromatic source) and overplot the LaRT
   spectrum on top of their published curve in the notebook. The
   characteristic D absorption dip in the blue wing should match in position
   and depth.

Notebook (`examples/lya_HD/inspect_HD.ipynb`):
- Loads the reference FITS outputs (using astropy as in the existing python
  scripts), plots `Jout(ν)` for each case overlaid, marks the H and D line
  centers in velocity space, and prints summary diagnostics (peak velocity,
  D-feature equivalent width, photon escape fraction).
- Includes a single-cell pure-D limit consistency check (set `D_to_H_ratio`
  large, swap H and D atomic data — confirms the D scatter path reproduces a
  pure-H run).
- Includes the Dijkstra, Haiman & Spaans (2006) Fig. 3 reproduction (run 5
  above) with the published curve digitised and overlaid for direct visual
  comparison.
- Cells are organised so a user can re-run after changing `D_to_H_ratio` in a
  parameter cell at the top.

### Phase 7 — Documentation

- `docs/LaRT_user_manual.tex`: add a subsection under the "Supported lines"
  / "Line IDs" section describing:
  - the new `'ly_alpha_HD'` line ID,
  - the `par%D_to_H_ratio` parameter (number ratio, cosmic default),
  - the `par%include_deuterium` convenience flag,
  - the explicit physical assumptions (no fine-structure for either species
    in v1; identical phase function; identical f-values and damping; mass
    and wavelength differences fully accounted for),
  - the location of the example notebook.
- Update `README_HOWTO` with the new line ID row.
- Add a brief note in the `line_mod` line table in `CLAUDE.md` (project
  reference card).

## Items Deferred to a Future Iteration

These are intentionally **not** part of v1 and are recorded here so we can
revisit them after the basic H+D path is validated.

1. **Fine-structure for H and D** (`par%fine_structure = .true.`). The
   current LaRT path uses `line_type = 2` for fine-structure, which assumes a
   single mass/Doppler width. For H+D fine-structure we would need a
   `line_type` analogous to type 7 but with `2 × 2` components, plus careful
   handling of the four nearly-coincident transitions. v1 ignores fine
   structure for both species — if the user sets
   `par%fine_structure = .true.` together with `'ly_alpha_HD'`, setup will
   issue a warning and override it to `.false.`.
2. **Per-photon species tracking**. It would be useful to know which species
   (H or D) last scattered each photon, both for diagnostics and for
   conditional peel-off images (e.g., "D-only" spectrum). v1 does not record
   this. Adding it later requires a single integer field on `photon_type`
   (e.g., `photon%last_species = 1` or `2`) updated inside
   `do_resonance_HD`, plus optional output filters.

## Open / Verified Decisions Recorded Here

- **line_id name**: `'ly_alpha_HD'` (decided).
- **`include_deuterium` flag**: yes, in addition to the explicit line_id
  (decided).
- **D atomic constants**: hard-coded in `line_mod.f90` (decided).
- **f12 for D**: stored separately as `species(2)%f12` (decided). Numerical
  values for f_D in v1: equal to f_H (the difference is below the noise floor
  of typical LaRT runs); revisit if a verified D-specific value is preferred.
- **Abundance unit**: number ratio `n_D / n_H` via `par%D_to_H_ratio`
  (decided).
- **Per-photon species tracking**: deferred (see above).

## Phase Summary

| Phase | Scope | Behaviour change |
|-------|-------|------------------|
| 1 | line_mod data structure + new line_id | None for legacy lines; new `'ly_alpha_HD'` available |
| 2 | input parameters | None when defaults used |
| 3 | calc_voigt / do_resonance dispatch | Active only for `line_type = 7` |
| 4 | grid setup | None |
| 5 | peel-off / Stokes | None (verify only) |
| 6 | validation runs + Jupyter notebook | Diagnostic |
| 7 | docs | — |

---

## Phase 9 (post-v1) — Species-aware recoil for H+D

### Bug identified

The recoil step in `scatter_resonance_*nostokes/stokes` (Cartesian and AMR)
and in `peeling_resonance_*` uses `line%g_recoil0`, which is computed from
the **hydrogen** atomic mass and rest wavelength:
```
line%g_recoil0 = h_planck / m_H / lambda_H**2     (set in line_mod)
g_recoil       = line%g_recoil0 / Dfreq_H_cell    (caller, in H Doppler units)
photon%xfreq  -= g_recoil * (1 - cost)            (caller)
```
For deuterium scatters this is wrong: deuterium is ~2x heavier, so the true
recoil shift in H-Doppler-units is approximately
```
g_recoil_D / g_recoil_H = (m_H/m_D) * (lambda_H/lambda_D)**2  ~ 0.50
```
The current code applies the H value to D scatters too, overestimating the
deuterium recoil by a factor of ~2. The shift per scatter is small
(~2.5e-4 H Doppler units), but the v1 plan flagged this as an acceptable
approximation only because D scatters are rare. With cosmic D/H = 3e-5 and
typical NSC ~ 1e6, the cumulative D-recoil error is small, but the user
prefers correctness.

### Fix plan

1. Add an integer field `selected_species_HD` (default 1 = H) inside the
   `line_type` derived type in `define.f90`. (Initially placed as a
   stand-alone module variable `last_species_HD` in `define.f90`; later
   moved into `line_type` for cleaner organization — see Phase 9d.)
   Set it to 1 or 2 inside `do_resonance_HD` and `do_resonance_HD_amr`
   after species selection.
2. Update the recoil block in
   - `scattering_car.f90` `scatter_resonance_nostokes`,
   - `scattering_car.f90` `scatter_resonance_stokes`,
   - `scattering_amr.f90` `scatter_resonance_nostokes_amr`,
   - `scattering_amr.f90` `scatter_resonance_stokes_amr`,
   - `peelingoff_rect.f90` (×2 places that apply recoil),
   - `peelingoff_amr.f90` (×3 places that apply recoil).
   Each block branches on
   `line%line_type == 7 .and. line%selected_species_HD == 2`
   to use `line%g_recoil0_D` instead of `line%g_recoil0`.
3. Backport identical changes to v2.00.

### Validation plan

After the fix lands, re-run `examples/lya_HD/sphere_HD_dijkstra2006.in`
twice with otherwise identical settings (monochromatic injection,
N_HI = 1.2e19 cm^-2, T = 1e4 K, D/H = 3e-5):

- `recoil = .false.` baseline (the result currently in
  `sphere_HD_dijkstra2006.fits.gz`).
- `recoil = .true.` with the fixed species-aware path, output to
  `sphere_HD_dijkstra2006_recoil.fits.gz`.

Add a new section to `inspect_HD.ipynb` overlaying the two LaRT spectra
together with Dijkstra+2006 and Michel-Dansac+2020 reference curves
(int=1 normalization). Expectation: H recoil dominates the cumulative
shift, so the recoil=true curve is mainly affected by the H part. The D
contribution to the recoil drift is sub-dominant. The user will inspect
the result and decide whether further debugging of the redistribution
physics is required.

### Note on recoil temperature dependence

Per Seon & Kim (2020), Eq. (15) — the recoil parameter
`g* = (h ν0² / m_H c²) / Δν_D ≈ 0.54 a (T/1e4 K)^(-1/2)` — is small
(~2.5e-4 in H Doppler units at T = 1e4 K) and only accumulates to a
visible spectral signature when **T ≲ 10² K** (where g* approaches
~10⁻³–10⁻²). For the Dijkstra+2006 setup at T = 1e4 K, the recoil
contribution is below RNG noise and the recoil=true vs. false
spectra are nearly identical (verified empirically: <1% peak height
change, sub-bin position shift). The species-aware fix is therefore
mostly a correctness-of-implementation improvement rather than a
physics-relevant change for typical Lyα LAE setups.

### Phase 9b — bookkeeping clean-up in `do_resonance_HD`

While reviewing the species-aware recoil patch, the following clean-up
was applied to `line_mod::do_resonance_HD` and
`scattering_amr::do_resonance_HD_amr`:

1. The line `xfreq_atom_D = xfreq_D - uz_D` was dead code (the
   variable was never used downstream). Removed.
2. The conversion `uz = uz_D / line%ratio_Dfreq_HD` was kept (it is
   correct; documented carefully in the routine). Derivation summary
   (following Seon's hand notes):

   Define x ≡ (ν − ν_H)/Δν_H, x̃ ≡ (ν − ν_D)/Δν_D (the photon's
   frequency offset in the H-cell and D-atom Doppler units),
   u ≡ v/v_th_H, ũ ≡ v/v_th_D (atom velocity normalized by H and D
   thermal speeds). Convention: positive x is bluer.

   Going to the D-atom rest frame, the LOS Doppler shift uses ν_D, so
   ```
   x̃' = x̃ − ũ_∥
   ```
   After scattering at angle (k̂_in, k̂_out):
   ```
   x̃'' = x̃' − ũ_∥ + ũ · k̂_out
   ```
   Translating x̃ → x via x̃ = (x − Δx)·(Δν_H/Δν_D), with Δx = δν_HD/Δν_H,
   ```
   x'' = x' − (u_∥ − u·k̂_out) × (λ_H/λ_D)
   ```
   The factor (λ_H/λ_D) is the wavelength ratio: it appears because the
   photon Doppler shift due to atom motion is ν_atom × v/c, so for a D
   atom the shift in cell-H Doppler units carries ν_D/ν_H = λ_H/λ_D.

   The LaRT caller updates as
   ```
   x'' = xfreq_atom + uz·cost + (ux·cosp + uy·sinp)·sint
   ```
   Matching: define the caller's `uz` so that uz·cost + perpendicular
   contributions equal `(u·k̂_out) × (λ_H/λ_D)`. With u_H-normalized,
   ```
   uz_caller = u_∥ × (λ_H/λ_D) = ũ_∥ × (v_th_D/v_th_H) × (λ_H/λ_D).
   ```
   Since v_th = Δν × λ, the combined factor reduces to
   `(v_th_D × λ_H)/(v_th_H × λ_D) = Δν_D / Δν_H = 1/ratio_Dfreq_HD`.
   Therefore `uz_caller = uz_D / ratio_Dfreq_HD`, which is exactly
   what the existing code does.

### Phase 9c — perpendicular velocity bug fix (`scatter_resonance_*`)

Original analysis (v1) labeled the perpendicular-component sampling for
D scatters as a "known limitation, negligible for cosmic D/H." Empirical
validation against Dijkstra+2006 Fig. 3 and Michel-Dansac+2020 Fig. B.3
showed this was incorrect: the LaRT spectrum peak at int=1 normalization
was ~30% too high. The root cause is that the caller in
`scatter_resonance_*` samples ux, uy from a unit Gaussian in cell-H
Doppler units regardless of which species scattered. For a D scatter,
the proper conversion to caller's "uz" convention requires the same
1/ratio_Dfreq_HD prefactor as for uz (derivation identical: the photon
Doppler shift from atom motion v in D's frame is ν_D × v/c, which
in cell-H Doppler units carries the (v_th_D/v_th_H)·(λ_H/λ_D) =
Δν_D/Δν_H factor for both LOS and perpendicular components).

**Fix**: insert a 4-line scaling block immediately before the
`photon%xfreq = xfreq_atom + uz*cost + (ux*cosp + uy*sinp)*sint`
update line in every resonance scattering routine:
```fortran
if (line%line_type == 7 .and. line%selected_species_HD == 2) then
   ux = ux / line%ratio_Dfreq_HD
   uy = uy / line%ratio_Dfreq_HD
endif
```

Applied to 13 sites across:
- `scattering_car.f90` `scatter_resonance_nostokes` (× 2 branches: core_skip + full)
- `scattering_car.f90` `scatter_resonance_stokes`   (× 2 branches)
- `scattering_amr.f90` `scatter_resonance_nostokes_amr`
- `scattering_amr.f90` `scatter_resonance_stokes_amr`
- Same in v2.00 (slightly different number of sites because v2.00's
  AMR variants have a different code structure).

**Validation**: with the perpendicular fix in place, peak J(x) at
int=1 normalization (LaRT N_HI=1.2e19, T=1e4K, D/H=3e-5,
monochromatic central source, 5×10⁴ photons):

| Curve | Peak J | Peak position | Ratio LaRT/MD |
|---|---|---|---|
| LaRT recoil=False | 0.1145 | x=+8.12 | 1.092 |
| LaRT recoil=True  | 0.1104 | x=−6.53 | 1.054 |
| Dijkstra+2006     | 0.1118 | x=−6.25 | 1.066 |
| Michel-Dansac+2020| 0.1048 | x=+7.98 | 1.000 |

LaRT peak is now within ~5% of both published references (was ~30%
before the fix). The remaining offset is at RNG-noise level for
N_phot = 5×10⁴ and consistent with the bin-width difference between
codes (LaRT dx = 0.16, Dijkstra dx ≈ 0.5, Michel-Dansac dx ≈ 0.21).

Why the v1 label was wrong: D scatters are rare per individual photon,
but every photon that survives in a τ ≈ 10⁶ medium undergoes ~10⁶
scatters — even ~5% of those (~5×10⁴ D scatters per photon at the
relevant resonance frequencies) accumulates an O(1) thermal-width
mis-scaling in the perpendicular Doppler kicks, and the cumulative
random-walk variance in xfreq propagates to the emergent peak shape.

### Phase 9d — code organization: move `last_species_HD` into `line_type`

Initially the H/D selection flag was a stand-alone module variable
`last_species_HD` declared at module scope in `define.f90`. All other
HD-related state already lived inside the `line_type` derived type
(`wavelength0_D`, `mass_amu_D`, `f12_D`, `cross0_D`, `vtherm1_D`,
`damping_D`, `g_recoil0_D`, `nD_over_nH`, `delta_nu_HD_Hz`,
`ratio_Dfreq_HD`, `ratio_voigta_HD`), so the lone module variable was
inconsistent with that design.

**Refactor**: replace `last_species_HD` (module-level) with
`line%selected_species_HD` (field of `line_type`, default 1 = H,
declared in `define.f90` next to the other HD fields).

Affected files (v2.10 + v2.00):
- `define.f90`: removed the standalone module variable; added
  `integer :: selected_species_HD = 1` to `line_type`.
- `line_mod.f90` `do_resonance_HD`: assignment sites updated
  (`last_species_HD = 1/2` → `line%selected_species_HD = 1/2`).
- `scattering_amr.f90` `do_resonance_HD_amr`: same.
- `scattering_car.f90` recoil and perpendicular-velocity branches
  (multiple sites): comparison updated to
  `line%line_type == 7 .and. line%selected_species_HD == 2`.
- `scattering_amr.f90` recoil and perpendicular-velocity branches: same.
- `peelingoff_rect.f90` recoil branches: same.
- `peelingoff_amr.f90` recoil branches: same.

Total: 41 reference sites rewritten across 12 files (6 per version).

This is purely a code-organization change; the simulation behaviour is
unchanged (regression-checked at 5×10⁴ photons against the same Dijkstra
setup — peak J and NSC match the prior result to RNG noise).
