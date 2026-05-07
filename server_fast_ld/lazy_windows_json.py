"""
lazy_windows_json.py — manages per-chromosome `<chrom>.windows.json`
files for the fast_ld endpoint.

Built on demand from the atlas chromosome JSON + sites file. Cached on
disk; stale if either source mtime is newer than the cached file.

Used by fast_ld_endpoint.py — but isolated as a separate module so it
can be tested without a running server.
"""

from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path
from typing import Optional


class WindowsJsonCache:
    """Resolves a `<chrom>.windows.json` path, building on demand.

    Parameters
    ----------
    cache_dir : Path
        Directory where built `<chrom>.windows.json` files live.
    atlas_json_dir : Path
        Directory where atlas chromosome JSONs (e.g., LG28.json) live.
    sites_dir : Path
        Directory where per-chrom sites files (`<chrom>.sites.tsv.gz` from
        STEP_A01) live.
    builder_script : Path
        Path to `build_windows_json.py`.
    """

    def __init__(self, *, cache_dir: Path, atlas_json_dir: Path,
                 sites_dir: Path, builder_script: Path):
        self.cache_dir = Path(cache_dir)
        self.atlas_json_dir = Path(atlas_json_dir)
        self.sites_dir = Path(sites_dir)
        self.builder_script = Path(builder_script)
        self.cache_dir.mkdir(parents=True, exist_ok=True)

    def _cached_path(self, chrom: str) -> Path:
        return self.cache_dir / f"{chrom}.windows.json"

    def _atlas_path(self, chrom: str) -> Path:
        return self.atlas_json_dir / f"{chrom}.json"

    def _sites_path(self, chrom: str) -> Path:
        return self.sites_dir / f"{chrom}.sites.tsv.gz"

    def _is_stale(self, cached: Path, atlas: Path, sites: Path) -> bool:
        if not cached.exists():
            return True
        c = cached.stat().st_mtime
        if atlas.exists() and atlas.stat().st_mtime > c:
            return True
        if sites.exists() and sites.stat().st_mtime > c:
            return True
        return False

    def get_or_build(self, chrom: str) -> Path:
        cached = self._cached_path(chrom)
        atlas = self._atlas_path(chrom)
        sites = self._sites_path(chrom)
        if not atlas.exists():
            raise FileNotFoundError(f"atlas JSON missing for {chrom}: {atlas}")
        if not sites.exists():
            raise FileNotFoundError(f"sites file missing for {chrom}: {sites}")
        if not self._is_stale(cached, atlas, sites):
            return cached
        # Build
        cmd = [
            sys.executable, str(self.builder_script),
            "--atlas-json", str(atlas),
            "--sites", str(sites),
            "--out", str(cached),
        ]
        proc = subprocess.run(cmd, capture_output=True, text=True)
        if proc.returncode != 0:
            raise RuntimeError(
                f"build_windows_json.py failed for {chrom} "
                f"(rc={proc.returncode}):\n{proc.stderr}")
        if not cached.exists():
            raise RuntimeError(
                f"build_windows_json.py did not produce {cached}\n"
                f"stderr: {proc.stderr}")
        return cached
