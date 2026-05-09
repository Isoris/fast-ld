# MODULE_5C — Inversion LD

Linkage disequilibrium analysis module for inversion candidates. Produces whole-chromosome and candidate-zoom LD heatmaps, marker LD contrast scoring with bootstrap confidence, and profile similarity clustering.

## Dependencies

- **Upstream:** MODULE_5A STEP07 (BEAGLE files), STEP08 (dosage), STEP10 (candidate table)
- **Grouping:** MODULE_5B STEP17c (contrast group manifests)
- **External:** ngsLD

## Pipeline

### Compute (`compute/`)

```
STEP_LD_00  : Thin BEAGLE to ~2000 evenly-spaced markers per chromosome
STEP_LD_00b : ngsLD on sparse BEAGLEs (whole-chromosome, unlimited distance)
STEP_LD_00c : ngsLD on candidate zoom regions (regional subset, unlimited distance)
STEP_LD_00d : ngsLD per inversion group (FULL_A, FULL_B, HALF separately)
STEP_LD_01  : Parse ngsLD output → binned numpy matrix cache (pickle)
```

### Analysis (`analysis/`)

```
marker_ld_contrast_v4_fixed_panels.py     — Latest (v4): fixed bootstrap panels,
                                            per-marker contrast, profile similarity,
                                            coherent cluster extraction
marker_ld_contrast_bootstrap_v3_auto_seed.py — Reference (v3): auto-seed with
                                               sample-name remapping
```

### Plot (`plot/`)

```
STEP_LD_02_plot_heatmaps.py                — LD heatmaps from binary cache
                                             (full / upper_half / triangle_down / triangle_up)
plot_ld_contrast_score_vs_pos.py           — Contrast score vs genomic position
plot_marker_ld_contrast_v4.py              — v4 contrast with cluster overlay
plot_direct_ld_from_cluster.py             — Direct LD heatmap of largest profile cluster
plot_marker_profile_similarity_heatmap.py  — Cross-marker LD profile correlation matrix
make_top_marker_ld_heatmap.py              — Focused LD heatmap of top contrast markers
```

### Orchestrator

```
run_ld_candidate_suite_v4.py   — Three-product orchestrator:
                                  whole-chromosome + candidate zoom + negative control
run_ld_candidate_suite_v4.slurm — SLURM launcher
```

## Quick start

```bash
# Run the full LD suite for all candidates
sbatch run_ld_candidate_suite_v4.slurm
```

## Notes

- v1 prototype (`marker_ld_contrast_bootstrap.py`) and old R versions are in `deprecated/`.
- The v4 fixed-panel approach ensures bootstrap reproducibility across markers.
- Group-specific LD (STEP_LD_00d) uses STEP17c sample lists for consistent grouping.
