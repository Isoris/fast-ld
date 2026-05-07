# engine_fast_ld

Standalone fast LD engine. C binary + Python wrapper + windows-JSON
preprocessor. Imported by `server_turn1/ld_endpoint.py` (the LD endpoint
module of the merged atlas server) to power the page-3 LD panel. Not
specific to the atlas — could be reused by anyone with a BEAGLE-imputed
dosage TSV + an atlas-style windows JSON.

## Files

| file | what |
|---|---|
| `fast_ld.c` | Single-file C engine. ~1100 lines. Reads control TSV, computes Pearson r² on standardized dosages, emits uint8 upper-triangle pairs + summary. |
| `Makefile` | `make` builds `fast_ld`. Needs `gcc` + `zlib` + OpenMP (`libgomp` on Ubuntu, `libomp` on macOS). |
| `fast_ld_wrapper.py` | Python adapter. `compute_ld(req)` → JSON-shaped result. Used by the server endpoint. |
| `build_windows_json.py` | Preprocessor. Reads atlas chromosome JSON + sites file, emits per-chrom windows JSON with explicit SNP membership per window. |
| `test_fast_ld.py` | 8 engine tests. Synthetic data, includes a numpy-correctness check (within 1/255 of `corrcoef²`). |
| `test_wrapper.py` | 2 wrapper tests. End-to-end through the binary. |
| `bench.py` | Realistic-scale benchmark (226 samples × varying SNP counts). |

## Build + test

```
make            # builds fast_ld
make test       # runs both test files
make bench      # realistic-scale benchmark
make static     # statically-linked binary (for distribution)
```

## Method

For each group:
1. Build dosage submatrix `D[snp_i, group_sample_k]` (N × K).
2. Per row, standardize: subtract row mean, divide by row sd
   (population, ddof=0). NaN entries → 0 after standardization.
3. Compute outer product `R[i, j] = (1/K) Σ_k Z[i,k] · Z[j,k]` for
   `i < j`. Square: `r²[i, j] = R[i, j]²`.
4. Quantize to uint8 (`floor(r²·255 + 0.5)`).
5. Emit upper triangle as flat binary.

Output matches `numpy.corrcoef(D)²` to within 1/255 (the quantization
bound) — see `test_fast_ld.py::test_known_correlation_matches_numpy`.

## Inputs

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

The dosage file follows the format from
`Modules/MODULE_x/STEP_A01_beagle_to_dosage_by_chr.py`:
- gzipped TSV
- header: `marker\tsample1\tsample2\t...`
- rows: `<chrom>_<pos>\tdosage1\tdosage2\t...`
- dosage in `[0, 2]`

If `n_unique_snps_in_window_range > snp_cap` and `thin_to` is unset,
the engine fails loudly. Three options: smaller window range, raise
`snp_cap`, or set `thin_to=N`.

## Outputs

In `out_dir`:

- **`pairs.<groupname>.bin`** — uint8 r²·255, upper triangle row-major,
  N*(N-1)/2 bytes. Pair index for `(i, j)`, `i < j`:
  `idx = i * (2N - i - 1) / 2 + (j - i - 1)`.

- **`sites.bin`** — packed struct, **56 bytes per kept SNP**:
  ```
  int32  idx            (0..N-1)
  int32  pos            (bp)
  float  maf[4]         (per-group, MAX_GROUPS=4)
  float  var[4]         (per-group dosage variance)
  int32  n_complete[4]  (per-group non-NaN count)
  ```
  Little-endian. Trailing slots zero for unused groups.

- **`summary.tsv`** — one `key\tvalue` per line. See
  `fast_ld_wrapper.py::_read_summary()` for the full schema.

## Performance reference

On 226 samples × sliding 100-SNP/step-20 windows (matching
`MS_Inversions` LG28):

| use case | unique SNPs | wallclock |
|---|---|---|
| 4-window inversion (smallest case) | ~150 | 30 ms |
| 50-window candidate | ~1080 | 350 ms |
| 250-window large region | ~5080 | 8 s |
| whole chrom thinned to 2000 | 2000 | 1.7 s |

Compute is `O(n_snps² × n_samples)` per group. The 8-second case is
near the live-request ceiling; for whole-chromosome views, use
`thin_to`.

## Validation against ngsLD on real data

**Not done in this build** — requires a candidate where ngsLD has
already produced `pairs_r2.tsv` on real LANTA data. Procedure:

1. Pick one inversion candidate.
2. Run fast_ld on the same SNP set + same group.
3. Compare shelf-ratio between engines. Pass: agreement within ~5%.

If shelf-ratios diverge significantly, investigate by stratifying on
MAF and INFO score — dosage-corr is known to slightly inflate r² on
low-MAF or low-INFO sites.

## Limitations

- Dosage input only, no GLs.
- No per-pair complete-case missing handling (NaN dosage → 0 after
  standardization).
- No HWE/MAF prefilter at engine level — filter upstream if needed.
- Max 4 groups, max 4096 samples (compile-time constants).
