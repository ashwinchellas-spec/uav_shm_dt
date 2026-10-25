# Physics-Informed AI Digital Twin for SHM of Composite UAV Wings

**Option B:** in-house FSDT MATLAB solver → AI digital twin

| Milestone | Workflow steps | Status |
|---|---|---|
| **M1** Problem definition & wing modelling | 1 | ✅ done — `reports/Milestone1_Report.md` |
| **M2** In-house FSDT solver + validation | 2 | ✅ done — `reports/Milestone2_Report.md` (56/56 checks pass) |
| **M3** DOE + high-fidelity simulation database | 3–4 | ✅ done — `reports/Milestone3_Report.md` (3 200 simulations, 0 failures) |
| **M4** Data processing & physics-informed features | 5 | ✅ done — `reports/Milestone4_Report.md` (8 400 records, 117 features) |
| **M5** Physics-informed AI models | 6 | ✅ done — `reports/Milestone5_Report.md` (PI-Net + hybrid physics-AI surrogate) |
| M6 Digital twin, virtual sensing, prognostics, validation | 7–10 | next |

## Quick start (MATLAB R2018b+ or GNU Octave 6+)

```matlab
cd uav_shm_dt
setup_paths            % adds all folders to the path
run_milestone1         % Step 1: wing definition, laminate, V-n, loads   (~5 s)
run_all_validation     % Step 2: 9 benchmark groups, prints PASS/FAIL     (~30 s)
run_milestone2         % Step 2: static / modal / gust / damage on wing   (~1 min)
run_milestone3         % Steps 3-4: DOE + 3200 simulations (~1.6 h CPU, resumable) + analysis
run_milestone4         % Step 5: sensors, virtual measurements, features, evaluation (~1 min)
run_milestone5         % Step 6: physics-informed networks + hybrid surrogate (~20 min Octave)
% or everything at once:
run_all_milestones
```

All figures go to `results/figures/`, and all numbers to `results/data/*.txt`.

No toolboxes are required, only base MATLAB. The code is written to run identically in MATLAB and Octave.

## Baseline wing (Milestone 1)

- Semi-span 1.0 m, chord 0.25 m, NACA 0012 equivalent, clamped at the root.
- Sandwich layup **[0/45/−45/90/core/90/−45/45/0]**: T300/5208 skins (0.125 mm plies) on a 21.75 mm Rohacell 51 WF core.
- 6 kg UAV, V_C = 25 m/s, V_D = 32 m/s. Design limit n = +5.11 (gust-critical), ultimate 7.67.
- FE model: 40 × 10 Q4 FSDT elements (Hellinger–Reissner mixed shear), 2255 DOF.

## Key solver features (Milestone 2)

- FSDT (Reissner–Mindlin) laminated-plate element with **Hellinger–Reissner mixed transverse shear**: locking-free, passes the patch test, no spurious modes.
- **Energy-equivalent shear correction factors** (correct for sandwich; k = 5/6 is 9× wrong).
- Static, modal (mass-normalised, sign-fixed modes) and Newmark-β transient analysis with Rayleigh damping.
- Ply stresses, equilibrium-based interlaminar/core shear stresses, and Tsai–Wu, Hashin, max-stress, core and ILSS indices.
- Physics-based aero loads: Schrenk span-wise, thin-airfoil chord-wise with singular quadrature, inertia relief, V–n envelope with Pratt gust loads.
- Continuum-damage module with 5 damage modes (matrix, fibre, debond, core crush, delamination).
- Closed-form Navier FSDT solution for verification.

## Physics-informed AI (Milestone 5)

- **PI-Net**: multi-task network (damage class, damaged face, damage map, location, extent, residual strength). Physics enters in four ways:
  1. healthy-referenced (Mahalanobis-whitened) physics features as inputs;
  2. an FE damage-evidence map as a prior on the map head;
  3. a **physics-consistency loss**: the predicted damage field must reproduce the measured feature changes through FE sensitivities (eigenvalue perturbation + adjoint statics), which also uses unlabelled data;
  4. hand-derived back-propagation, verified to ~1e-7.
- **Hybrid surrogate**: a reduced-order beam/CLT model (Ritz, Dunkerley, Tsai–Wu) times a neural correction of its model-form error. It is far more accurate than a pure network for the same number of FE runs, and it extrapolates.
- Everything (MLP, Adam, losses, gradients) is written from scratch; there are no toolboxes.

See `reports/` for the full theory and results.
