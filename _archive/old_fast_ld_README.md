# fast_ld — windows-driven LD engine

Single-file C engine + Python wrapper that computes pairwise r² on
imputed dosage data, driven by a per-chromosome **windows JSON** built
from your atlas's chromosome JSON. Replaces ngsLD for inversion
heatmap visualization in the *C. gariepinus* atlas.

This is **not** a general ngsLD replacement. It assumes:

- Your dosages come from BEAGLE-imputed GLs (no NaNs in practice).
- Your goal is plotting LD blocks within a candidate region.
- You can afford the dosage-correlation approximation instead of EM.

## Pipeline

```
  Atlas chromosome JSON   +   sites.tsv.gz              (preprocessing,
  (windows: bp ranges)        (SNP positions)            once per chrom)
                  │            │
                  ▼            ▼
          build_windows_json.py
                  │
                  ▼
       <chrom>.windows.json
       (lookup: window_idx → snp_positions[100])
                  │
                  ▼
  ┌──── per-request input ──────────────────────────────────────────────┐
  │  control file:                                                       │
  │    dosage_path=...                                                   │
  │    windows_json=<chrom>.windows.json                                 │
  │    window_range=2150-2400                                            │
  │    group=HOM_INV:S001,...     group=HOM_REF:S003,...                 │
  │    shelf_start=15200000   shelf_end=18100000   (optional)            │
  │    snp_cap=5000  thin_to=  (optional override)                       │
  └────────────────────────────────────────────────────────────────────────┘
                  │
                  ▼
              fast_ld
                  │
                  ▼
         pairs.<group>.bin (uint8 r²·255, upper triangle)
         sites.bin           (per-SNP per-group MAF/var/n)
         summary.tsv         (medians, shelf-ratio, decay deciles)
```

## Build

```
make            # → src/fast_ld
make test       # → run all tests
make bench      # → realistic-scale benchmark
```

## Performance

On synthetic data resembling the 226-sample *C. gariepinus* hatchery
cohort (5 SNPs/kb, 100-SNP windows step 20):

| use case | unique SNPs | wallclock |
|---|---|---|
| 4-window inversion | 160 | 30 ms |
| 50-window candidate | 1080 | 350 ms |
| 250-window region | 5080 | 8 s |
| whole chrom thinned to 2000 | 2000 | 1.7 s |

Compute is `O(n_snps² × n_samples)` per group. The 8-second case is
near the ceiling we'd want for a live request; for whole-chromosome
views, use `thin_to` to bound compute.

## Preprocessing — `build_windows_json.py`

```
python3 python/build_windows_json.py \
    --atlas-json /path/to/LG28.json \
    --sites      /path/to/LG28.sites.tsv.gz \
    --out        /path/to/LG28.windows.json
```

Reads the atlas JSON's `windows` array (each window has `idx`,
`start_bp`, `end_bp`) and intersects with the SNP positions in
`sites.tsv.gz` (the canonical sites file from `STEP_A01`). Emits a
slim JSON with explicit `snp_positions` per window.

Both formats of `sites.tsv.gz` are auto-detected:

- STEP_A01 5-column: `marker, chrom, pos, allele1, allele2`
- ANGSD 4-column: `chrom, pos, major, minor`

For LG28 (~86K SNPs, 4302 windows), this takes ~1 s and produces a
~4 MB JSON.

## Engine inputs

The C binary takes one argument: a path to a TSV control file. Format:

```
dosage_path=/path/to/<chrom>.dosage.tsv.gz
windows_json=/path/to/<chrom>.windows.json
window_range=2150-2400          # inclusive
snp_cap=5000                     # default; reject if more unique SNPs
thin_to=                         # optional; if set, even-bp-thin to N
shelf_start=15200000             # optional, for shelf-ratio summary
shelf_end=18100000               # optional
out_dir=/path/to/output_dir
threads=4
group=NAME1:id1,id2,id3,...      # one line per group, max 4
group=NAME2:idA,idB,idC,...
```

If `n_unique_snps_in_window_range > snp_cap` and `thin_to` is unset,
the engine fails loudly with a clear message. Three options to
proceed: smaller window range, raise `snp_cap`, or set `thin_to=N`.

## Engine outputs

In `out_dir`:

- **`pairs.<groupname>.bin`** — uint8 r²·255, upper triangle row-major,
  `N*(N-1)/2` bytes. Pair index for `(i, j)`, `i < j`:
  `idx = i * (2N - i - 1) / 2 + (j - i - 1)`.
  N = `n_snps_used` from `summary.tsv`.

- **`sites.bin`** — packed struct, **56 bytes per kept SNP**:
  ```
  int32  idx            (0..N-1)
  int32  pos            (bp)
  float  maf[4]         (per-group, MAX_GROUPS=4)
  float  var[4]         (per-group dosage variance)
  int32  n_complete[4]  (per-group non-NaN count)
  ```
  Little-endian. Trailing slots zero for unused groups.

- **`summary.tsv`** — one `key\tvalue` per line. Fixed top-level keys
  plus per-group keys prefixed by group name. See `python/fast_ld_wrapper.py`
  for the full schema; in particular each group has:
  - `n_samples`
  - `compute_seconds`
  - `median_r2_overall`, `pct_pairs_above_0_8`
  - `median_r2_shelf`, `median_r2_flank`, `shelf_ratio`
  - `decay_decile_0` … `decay_decile_9`

## Python wrapper

```python
from fast_ld_wrapper import compute_ld

result = compute_ld({
    "dosage_path":  "/data/dosage/C_gar_LG28.dosage.tsv.gz",
    "windows_json": "/data/windows/C_gar_LG28.windows.json",
    "window_range": [2150, 2400],
    "groups": {
        "HOM_INV": [...],
        "HOM_REF": [...],
    },
    "shelf_bp":  [15_200_000, 18_100_000],
    "snp_cap":   5000,
    "thin_to":   None,
    "triangle_assign": {"lower": "HOM_REF", "upper": "HOM_INV"},
    "threads":   4,
}, bin_path=Path("/path/to/fast_ld"))

# result["matrices"]["HOM_INV"]["pairs_b64"]   uint8 upper triangle (b64)
# result["matrices"]["HOM_INV"]["r2"]          dense float32 NxN (set return_matrices=False to skip)
# result["sites"]                               per-SNP {idx, pos, maf_*, var_*, n_complete_*}
# result["summary"]                             top-level metrics
# result["triangle_assign"]                     pass-through
```

The wrapper does input validation, writes the control file, calls the
binary, and parses outputs. Return is JSON-safe — the localhost server
serializes `result` and ships to the browser.

## Method

For each group `g`:

1. Build dosage submatrix `D[i, k]` = `dosage[snp_i][group_sample_k]`.
2. Per row, standardize: subtract row mean, divide by row sd
   (population, ddof=0). NaN entries → 0 after standardization.
3. Compute outer product `R[i, j] = (1/k) Σ_k Z[i, k] · Z[j, k]` for
   `i < j`. `R` is the Pearson correlation on the original dosages.
4. Square: `r²[i, j] = R[i, j]²`. Quantize to uint8.
5. Emit upper triangle.

Validation: synthetic data with known structure, fast_ld output matches
`numpy.corrcoef(D)²` to within 1/255 across 1.4M pairs (the quantization
bound). See `tests/test_fast_ld.py::test_known_correlation_matches_numpy`.

## Validation against ngsLD on real data

**Not done in this build** — requires LANTA access plus a candidate
where ngsLD has produced `pairs_r2.tsv`. Procedure:

1. Run fast_ld on one inversion candidate.
2. Run ngsLD on the same candidate, same group, same SNP set.
3. Compare shelf-ratio (median r² in shelf / median r² in flanks)
   between engines. Pass criterion: agreement within ~5%.

If shelf-ratios diverge significantly, dosage-corr is inflating r² in
ways that affect inversion detection — investigate by stratifying on
MAF and INFO score. If they match, fast_ld is ratified for atlas use.

## Limitations

- Dosage input only, no GLs.
- No per-pair complete-case missing handling (NaN dosage → standardized 0).
- No HWE/MAF filter at engine level; filter upstream if needed.
- Max 4 groups (compile-time `MAX_GROUPS`).
- Max 4096 samples (compile-time `MAX_SAMPLES`).
- No automatic fallback if windows JSON is missing — call
  `build_windows_json.py` first.

## Testing

```
make test
```

Runs:
- `test_fast_ld.py` — 8 engine-level tests including JSON parse, snp_cap
  rejection, thin_to override, numpy correctness check, synthetic
  inversion signal, real LG28 windows JSON smoke test
- `test_wrapper.py` — wrapper end-to-end + validation errors
