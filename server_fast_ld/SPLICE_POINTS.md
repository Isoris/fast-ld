# turn 11c — server LD endpoint (fast_ld native)

FastAPI handler that wraps the fast_ld C engine, plus a lazy
windows-JSON cache. Replaces `server_turn11a_ld/` (which wrapped ngsLD).

## Files

| file | what |
|---|---|
| `fast_ld_endpoint.py` | FastAPI handler. Pydantic `FastLDReq`, composite caching, named-groups dict input. ~340 lines. |
| `lazy_windows_json.py` | Per-chromosome `<chrom>.windows.json` cache. Builds on first request, rebuilds when atlas JSON or sites file is newer. |
| `test_fast_ld_endpoint.py` | 11 tests including a real end-to-end through the actual fast_ld binary. |

Run tests: `python3 test_fast_ld_endpoint.py` (requires the engine
built — see `../engine_fast_ld/`).

## Drop-in pattern

Copy the two `.py` files into your `popstats_server/` directory
alongside `popstats_server.py`. They're standalone — only `numpy`,
`fastapi`, `pydantic` (already in `server_turn1/requirements.txt`)
plus `fast_ld_wrapper.py` from `engine_fast_ld/`.

```
popstats_server/
├── popstats_server.py
├── fast_ld_endpoint.py     ← THIS
├── lazy_windows_json.py    ← THIS
├── fast_ld_wrapper.py      ← from engine_fast_ld/
├── build_windows_json.py   ← from engine_fast_ld/
├── ld_endpoint.py          ← turn 11a, can stay if you want a ngsLD fallback
└── ...
```

## Patch popstats_server.py

Two edits, both small.

**Edit A — imports near the top:**

```python
from pathlib import Path
from fast_ld_endpoint import FastLDReq, handle_split_heatmap, fast_ld_engine_hash
from lazy_windows_json import WindowsJsonCache
```

**Edit B — bootstrap + route, near other endpoint mounts:**

```python
# At server bootstrap (next to where engines are loaded):
WINDOWS_CACHE = WindowsJsonCache(
    cache_dir=Path(CFG["windows_cache_dir"]),
    atlas_json_dir=Path(CFG["atlas_json_dir"]),
    sites_dir=Path(CFG["sites_dir"]),
    builder_script=Path(CFG["build_windows_json_path"]),
)
FAST_LD_BIN = Path(ENGINES.path("fast_ld"))
FAST_LD_HASH = fast_ld_engine_hash(FAST_LD_BIN)

@app.post("/api/ld/split_heatmap")
async def ld_split_heatmap(req: FastLDReq) -> JSONResponse:
    _ensure_ready()
    def cache_get(k): return CACHE.get_json(k, "ld")
    def cache_put(k, v): CACHE.put_json(k, "ld", v)
    payload = await handle_split_heatmap(
        req,
        fast_ld_bin=FAST_LD_BIN,
        fast_ld_engine_hash=FAST_LD_HASH,
        dosage_dir=Path(CFG["dosage_dir"]),
        windows_cache=WINDOWS_CACHE,
        sample_filter=SAMPLES.filter_known,
        cache_get=cache_get, cache_put=cache_put,
    )
    return JSONResponse(payload)
```

**This replaces the turn-11a route at the same path.** If you want both
running (for validation), rename one — e.g.,
`/api/ld/split_heatmap_ngsld` for the old one, `/api/ld/split_heatmap`
for this.

## Configure popstats_server_config.yaml

Add new keys + the engine path:

```yaml
dosage_dir:           /scratch/lt200308-agbsci/.../dosage_per_chrom
atlas_json_dir:       /scratch/lt200308-agbsci/.../atlas_data
sites_dir:            /scratch/lt200308-agbsci/.../sites_per_chrom
windows_cache_dir:    /scratch/lt200308-agbsci/.../cache/windows_json
build_windows_json_path: /scratch/lt200308-agbsci/.../tools/build_windows_json.py

engines:
  region_popstats: <existing>
  hobs_windower:   <existing>
  angsd_patched:   <existing>
  instant_q:       <existing>
  fast_ld:         /scratch/lt200308-agbsci/.../tools/fast_ld
```

The `build_windows_json_path` should point to the script in the
sibling `engine_fast_ld/build_windows_json.py`. The `fast_ld` engine
path should point to the compiled binary.

## Verify

```bash
curl -s http://127.0.0.1:8765/api/health | jq '.engines.fast_ld'
# → "abcd1234..." (16-char content hash)

curl -X POST http://127.0.0.1:8765/api/ld/split_heatmap \
  -H 'Content-Type: application/json' \
  -d '{
    "chrom": "C_gar_LG28",
    "window_range": [2150, 2200],
    "groups": {
      "ALL": ["..."],
      "HOM_INV": ["..."]
    },
    "shelf_bp": [15200000, 18100000]
  }' | jq '.timing,.cache_state,.matrices.HOM_INV.shelf_ratio'
```

First call builds `LG28.windows.json` (~1 s for LG28-sized chroms),
runs fast_ld, returns. Second call: instant cache hit.

## Request shape (named groups)

This endpoint follows the **named-groups dict pattern** that
`popstats_groupwise` uses, not turn 11a's positional
`all_samples`/`common_samples`.

```json
{
  "chrom": "C_gar_LG28",
  "window_range": [2150, 2200],
  "groups": {
    "ALL": ["S001", "S002", ...],
    "HOM_INV": ["S017", "S041", ...]
  },
  "triangle_assign": {"lower": "ALL", "upper": "HOM_INV"},
  "shelf_bp": [15200000, 18100000],
  "snp_cap": 5000,
  "thin_to": null,
  "threads": 4
}
```

Up to 4 groups. `triangle_assign` is optional — if omitted with exactly
2 groups, defaults to `lower=larger group, upper=other`.

## Response shape

```json
{
  "chrom": "C_gar_LG28",
  "window_range": [2150, 2200],
  "n_snps": 487,
  "n_pairs": 118341,
  "sites": {
    "idx": [0, 1, ...],
    "pos": [15200143, 15200389, ...],
    "maf_ALL": [...], "var_ALL": [...], "n_complete_ALL": [...],
    "maf_HOM_INV": [...], ...
  },
  "matrices": {
    "ALL": {
      "n_samples": 226,
      "n_pairs": 118341,
      "pairs_b64": "...",
      "median_r2_overall": 0.041,
      "median_r2_shelf": 0.387,
      "median_r2_flank": 0.046,
      "shelf_ratio": 8.42,
      "pct_pairs_above_0_8": 0.02,
      "decay_deciles": [...]
    },
    "HOM_INV": { ... }
  },
  "triangle_assign": {"lower": "ALL", "upper": "HOM_INV"},
  "summary": {...},
  "timing": {"compute_seconds": 0.350, "total_wallclock_seconds": 0.420},
  "cache_state": "miss",
  "cache_key": "fast_ld.split.abcd1234..."
}
```

`pairs_b64` is base64-encoded uint8 r²·255, upper triangle row-major.
The atlas decodes with `popgenLD.decodePairsB64()`.

## Cache behavior

Cache key hashes:
- `chrom`
- `window_range`
- `groups` (sample IDs sorted within each group)
- `shelf_bp`
- `snp_cap` / `thin_to`
- engine hash (binary content sha256)

Order of sample IDs within a group doesn't matter (sorted before
hashing). Adding/removing/renaming a group is a different key.
Engine binary changes invalidate all keys.
