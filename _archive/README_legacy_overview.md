# `phase_11_secondary/` — secondary confirmation tests on candidate inversions

Two orthogonal secondary confirmation tests that each predict a
specific signature of a real balanced inversion polymorphism. These
are **not** discovery modules — they're confirmation layers applied
after `phase_9_classification/` has produced a classified catalog and
`phase_10_followup/` has generated per-candidate data products.

Each module answers a different prediction:

| Module | Prediction tested |
|---|---|
| `MODULE_5C_Inversion_LD/` | Suppressed recombination → elevated interior LD |
| `MODULE_5D_Inversion_FST/` | Balanced polymorphism → elevated Fst between karyotype groups |

The third prediction — *heterokaryotype carriers → locally elevated
Hobs vs HWE* — used to live here as `MODULE_5E_Inversion_HOBS/`. It is
now handled by the per-group Hobs pair in `phase_qc_shelf/` (Q07b + Q07c);
see "Hobs confirmation moved" below.

A candidate that passes both remaining tests plus the phase_qc_shelf
Hobs check, in addition to the phase 3/4 breakpoint + boundary evidence,
is strongly supported as a real polymorphic inversion in the 226-sample
*C. gariepinus* cohort.

## Layout

```
phase_11_secondary/
├── MODULE_5C_Inversion_LD/           linkage disequilibrium
│   ├── README.md
│   ├── analysis/   compute/   plot/
│   ├── run_ld_candidate_suite_v4.py
│   └── run_ld_candidate_suite_v4.slurm
│
└── MODULE_5D_Inversion_FST/          Hudson Fst scans
    ├── README.md
    ├── STEP19_candidate_FST_scan.R
    └── STEP19_candidate_FST_scan.slurm
```

## Hobs confirmation moved

`MODULE_5E_Inversion_HOBS/` was archived on 2026-04-24 to
`inversion_modules/_archive_superseded/MODULE_5E_Inversion_HOBS_superseded_by_Q07b/`.
The replacement is the per-group pair in `phase_qc_shelf/`:

- **`STEP_Q07b_hobs_per_group.sh`** — runs the patched ANGSD
  (`Isoris/angsd_fixed_HWE`) with `-doHWE 1` separately on each
  karyotype group (Hom1 / Het / Hom2) from `invgt_assignments`. This is
  the scientifically correct Mérot approach — `Hexp` at an inversion SNP
  is computed from within-group allele frequencies rather than a
  cohort-pooled estimate.
- **`STEP_Q07c_hobs_windower.sh`** — runs the `hobs_windower` C binary
  at seven multi-scale sliding windows per group; the inversion
  signature is `HoverE_HET → 2.0` and `HoverE_HOM → 0` at
  arrangement-differentiating SNPs.

Outputs (`hobs_sites.*.tsv.gz`, `hobs_win.*.<SCALE>.tsv.gz`) are
consumed by `phase_10_followup/analysis/04_candidate_overlay.py` and
`phase_10_followup/figures/05_plot_hobs_hwe.R`.

See the archived module's `SUPERSEDED_NOTE.md` for a side-by-side
comparison and rationale.

## Dependencies

Each module reads:

- Candidate table + karyotype / contrast group assignments from phase 4
  (registry-backed; see `phase_7_karyotype_groups/proposal/`
  and `phase_7_karyotype_groups/validation/`).
- Per-candidate dosage / BEAGLE / GL data from phase 5 engines where
  relevant (residual dosage for LD; raw dosage for Fst).

The two remaining phase_6 modules do not depend on each other — they
can run in parallel.

## Naming note

The `MODULE_5*` prefixes persist from the pre-phase-refactor scheme
where these were direct sub-modules of MODULE_5. The folder paths are
referenced from several wrapper scripts
(`run_ld_candidate_suite_v4.py`, `STEP19_candidate_FST_scan.slurm`)
via self-relative paths. Renaming them requires a coordinated sweep
and is tracked in
`inversion_modules/phase_2_discovery/2c_precomp/RENAMING.md`
Category 3 (cross-module path renames, HPC-coordinated).

## Manuscript mapping

Both remaining modules feed into §3.7 (secondary confirmation) of
`MS_Inversions_North_african_catfish`. LD tracks appear as
per-candidate supplementary panels; the Fst shelf contrast is the
headline karyotype-vs-flanking statistic (Fst 0.308 vs flanking
~0.032–0.055 on the LG28 test case). Hobs confirmation comes from
phase_qc_shelf (see above) and plots are generated in
phase_10_followup.
