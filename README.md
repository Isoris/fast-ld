# fast-ld — windows-driven LD engine

Single-file C engine + Python wrappers + FastAPI server endpoint that
compute pairwise r² on imputed dosage data, driven by a per-chromosome
**windows JSON** built from the atlas's chromosome JSON. Replaces ngsLD
for inversion heatmap visualization in the *C. gariepinus* atlas.

This is **not** a general ngsLD replacement. It assumes:

- Your dosages come from BEAGLE-imputed GLs (no NaNs in practice).
- Your goal is plotting LD blocks within a candidate region.
- You can afford the dosage-correlation approximation instead of EM.

## Repo layout

```
fast-ld/
├── README.md
├── LICENSE
├── atlas_ld.js                     Browser renderer (split-heatmap canvas)
│
├── engines/                        C source + Makefile
│   ├── Makefile                    Build with: make -C engines
│   └── fast_ld.c
│
├── wrappers/                       Python wrappers / preprocessors
│   ├── fast_ld_wrapper.py          compute_ld(req) — JSON-shaped
│   └── build_windows_json.py       atlas_json + sites → windows.json
│
├── server/                         FastAPI server endpoint
│   ├── SPLICE_POINTS.md            How to splice into popstats_server
│   ├── fast_ld_endpoint.py
│   └── lazy_windows_json.py
│
├── tests/                          All tests + benchmark
│   ├── test_fast_ld.py             Engine tests (8 cases incl. numpy check)
│   ├── test_wrapper.py             Wrapper end-to-end
│   ├── test_fast_ld_endpoint.py    FastAPI handler tests (11 cases)
│   └── bench.py                    Realistic-scale benchmark
│
└── _archive/                       Legacy / superseded code
    ├── old_fast_ld_README.md
    ├── README_legacy_overview.md
    └── MODULE_5C_Inversion_LD/     Pre-fast_ld inversion-LD pipeline
        ├── README.md
        ├── STEP_LD_00_make_sparse_beagle.sh
        ├── STEP_LD_00b_run_ngsLD_chrsparse.slurm
        ├── STEP_LD_00c_run_ngsLD_invzoom.slurm
        ├── STEP_LD_00d_run_ngsLD_by_group.slurm
        ├── STEP_LD_01_build_binary_cache.py
        ├── marker_ld_contrast_v4_fixed_panels.py
        ├── marker_ld_contrast_bootstrap_v3_auto_seed.py
        ├── run_ld_candidate_suite_v4.py
        ├── run_ld_candidate_suite_v4.slurm
        └── plots/                  All MODULE_5C plot scripts
            ├── STEP_LD_02_plot_heatmaps.py
            ├── plot_ld_panels.R
            └── ...
```

## Build

```
make -C engines              # → engines/fast_ld
make -C engines test         # run engine + wrapper tests
make -C engines bench        # realistic-scale benchmark
make -C engines static       # statically-linked binary (for distribution)
```

## Pipeline

```
  Atlas chromosome JSON   +   sites.tsv.gz              (preprocessing,
  (windows: bp ranges)        (SNP positions)            once per chrom)
                  │            │
                  ▼            ▼
       wrappers/build_windows_json.py
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
              engines/fast_ld
                  │
                  ▼
         pairs.<group>.bin (uint8 r²·255, upper triangle)
         sites.bin           (per-SNP per-group MAF/var/n)
         summary.tsv         (medians, shelf-ratio, decay deciles)
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

## Preprocessing — `wrappers/build_windows_json.py`

```
python3 wrappers/build_windows_json.py \
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
the engine fails loudly. Three options to proceed: smaller window
range, raise `snp_cap`, or set `thin_to=N`.

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
  plus per-group keys prefixed by group name. See
  `wrappers/fast_ld_wrapper.py` for the full schema.

## Python wrapper

```python
import sys
sys.path.insert(0, "wrappers")
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
}, bin_path=Path("engines/fast_ld"))
```

The wrapper does input validation, writes the control file, calls the
binary, and parses outputs. Return is JSON-safe — the localhost server
serializes `result` and ships to the browser.

## Server endpoint

`server/fast_ld_endpoint.py` is a FastAPI handler that wraps the C
engine with caching and named-groups input. See `server/SPLICE_POINTS.md`
for how to splice it into a popstats_server.

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

## Limitations

- Dosage input only, no GLs.
- No per-pair complete-case missing handling (NaN dosage → standardized 0).
- No HWE/MAF filter at engine level; filter upstream if needed.
- Max 4 groups (compile-time `MAX_GROUPS`).
- Max 4096 samples (compile-time `MAX_SAMPLES`).
- No automatic fallback if windows JSON is missing — call
  `wrappers/build_windows_json.py` first.
