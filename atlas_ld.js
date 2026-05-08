// =============================================================================
// atlas_ld.js — turn 11c (fast_ld native)
// =============================================================================
//
// Linkage-disequilibrium split-heatmap renderer + page-3 panel, driven
// by the fast_ld engine's native uint8-pair-bytes output.
//
// Two triangles, one canvas:
//   - Lower triangle (matrix coords j ≤ i): r² across the "lower" group
//   - Upper triangle (j > i): r² across the "upper" group
// Group→triangle mapping comes from payload.triangle_assign.
//
// Diagonal direction: upper-left to lower-right. Matches sim_mat, the
// live dosage heatmap, and Q09b's correlation matrix. Don't switch only
// LD to Hi-C convention — consistency wins.
//
// Apex colormap (byte-faithful port from STEP_LD_02_plot_heatmaps.py):
//   ["#080808","#1a0a2e","#3d1261","#6a1b8a",
//    "#9c2e7a","#cc4466","#e8734a","#f5a623"]
// Power transform √r² applied by default (matches --power 0.5 in the
// old pipeline).
//
// Payload shape (from fast_ld_endpoint.py):
//   {
//     chrom, window_range, n_snps, n_pairs,
//     sites: { idx[], pos[], maf_<g>[], var_<g>[], n_complete_<g>[], ... },
//     matrices: {
//       <group_name>: {
//         pairs_b64,  // uint8 r²·255, upper triangle row-major
//         n_samples, n_pairs,
//         median_r2_overall, median_r2_shelf, median_r2_flank, shelf_ratio,
//         pct_pairs_above_0_8, decay_deciles[10],
//       },
//       ...
//     },
//     triangle_assign: { lower: <group_name>, upper: <group_name> },
//     summary, timing, cache_state, ...
//   }
//
// API (window.popgenLD):
//   .APEX_COLORS                — the 8 anchor hex strings
//   .apexColormap()             — Uint8ClampedArray length 256*4 (RGBA LUT)
//   .colorAt(lut, r2, opts)     — apply transform + colormap → [r,g,b,a]
//   .decodePairsB64(b64, n)     — Uint8Array of length n*(n-1)/2
//   .pairAt(arr, n, i, j)       — uint8 r²·255 for (i,j); symmetric, NaN diag
//   .drawSplitHeatmap(canvas, payload, opts)
//   .makeLDSplitPanel(opts)     — DOM panel: header, picker, canvas, caption
//   .ldSplitHeatmap(req, opts)  — request wrapper (delegates to popgenLive
//                                 if installed, else direct fetch)
//
// =============================================================================

(function (root, factory) {
  if (typeof module !== 'undefined' && module.exports) {
    module.exports = factory();
  } else {
    root.popgenLD = factory();
  }
}(typeof self !== 'undefined' ? self : this, function () {
  'use strict';

  // ===========================================================================
  // Apex colormap — byte-faithful port
  // ===========================================================================

  const APEX_COLORS = [
    '#080808', '#1a0a2e', '#3d1261', '#6a1b8a',
    '#9c2e7a', '#cc4466', '#e8734a', '#f5a623',
  ];

  function _hexToRGB(hex) {
    const h = hex.replace('#', '');
    return [
      parseInt(h.slice(0, 2), 16),
      parseInt(h.slice(2, 4), 16),
      parseInt(h.slice(4, 6), 16),
    ];
  }

  function apexColormap() {
    const K = APEX_COLORS.length;
    const lut = new Uint8ClampedArray(256 * 4);
    const anchors = APEX_COLORS.map(_hexToRGB);
    const segLen = 1 / (K - 1);
    for (let i = 0; i < 256; i++) {
      const t = i / 255;
      let seg = Math.floor(t / segLen);
      if (seg >= K - 1) seg = K - 2;
      const tLocal = (t - seg * segLen) / segLen;
      const a = anchors[seg];
      const b = anchors[seg + 1];
      lut[i * 4]     = Math.round(a[0] + (b[0] - a[0]) * tLocal);
      lut[i * 4 + 1] = Math.round(a[1] + (b[1] - a[1]) * tLocal);
      lut[i * 4 + 2] = Math.round(a[2] + (b[2] - a[2]) * tLocal);
      lut[i * 4 + 3] = 255;
    }
    return lut;
  }

  function colorAt(lut, r2, opts) {
    opts = opts || {};
    if (r2 == null || (typeof r2 === 'number' && isNaN(r2))) {
      return [255, 255, 255, 0];
    }
    const transform = opts.transform || 'power';
    const power = (opts.power != null) ? opts.power : 0.5;
    const vmax = (opts.vmax != null) ? opts.vmax : 1.0;
    let v;
    if (transform === 'raw') v = r2;
    else if (transform === 'power') v = Math.pow(Math.max(0, r2), power);
    else if (transform === 'log1p') {
      const k = (opts.log1p_k != null) ? opts.log1p_k : 10.0;
      v = Math.log1p(k * Math.max(0, r2));
    } else v = r2;
    let idx = Math.round((v / vmax) * 255);
    if (idx < 0) idx = 0;
    if (idx > 255) idx = 255;
    return [lut[idx * 4], lut[idx * 4 + 1], lut[idx * 4 + 2], lut[idx * 4 + 3]];
  }

  // ===========================================================================
  // base64 → Uint8Array
  // ===========================================================================

  function decodePairsB64(b64, n_expected) {
    let bytes;
    if (typeof Buffer !== 'undefined') {
      // Node
      const b = Buffer.from(b64, 'base64');
      bytes = new Uint8Array(b.buffer, b.byteOffset, b.byteLength);
    } else {
      // Browser
      const bin = atob(b64);
      bytes = new Uint8Array(bin.length);
      for (let i = 0; i < bin.length; i++) bytes[i] = bin.charCodeAt(i);
    }
    if (n_expected != null && bytes.length !== n_expected) {
      throw new Error(
        'decodePairsB64: expected ' + n_expected + ' bytes, got ' + bytes.length);
    }
    return bytes;
  }

  // Pair index for (i, j), i < j: idx = i*(2N-i-1)/2 + (j-i-1)
  function _pairIdx(N, i, j) {
    return (i * (2 * N - i - 1) >> 1) + (j - i - 1);
  }

  // Return uint8 r²·255 for (i, j). Symmetric. Returns NaN-flag (-1) on diagonal.
  function pairAt(arr, N, i, j) {
    if (i === j) return -1;
    if (i > j) { const t = i; i = j; j = t; }
    return arr[_pairIdx(N, i, j)];
  }

  // ===========================================================================
  // Canvas split renderer
  //
  // Cells are (i, j). Diagonal upper-left → bottom-right.
  //   - j > i: upper triangle (above the diagonal)
  //   - j < i: lower triangle (below the diagonal)
  //   - j == i: diagonal, drawn as a thin line in opts.diagColor
  // Each cell is drawn as a (canvas_w / N) × (canvas_h / N) rect. We use a
  // single ImageData write at the end for speed.
  // ===========================================================================

  function drawSplitHeatmap(canvas, payload, opts) {
    opts = opts || {};
    if (!canvas || !payload) return false;
    if (!payload.matrices || !payload.triangle_assign) return false;
    const ta = payload.triangle_assign;
    if (!ta.lower || !ta.upper) return false;
    const lowerMat = payload.matrices[ta.lower];
    const upperMat = payload.matrices[ta.upper];
    if (!lowerMat || !upperMat) return false;

    const N = payload.n_snps;
    if (!N || N < 2) return false;

    const ctx = canvas.getContext('2d');
    const W = canvas.width, H = canvas.height;
    const lut = opts.lut || apexColormap();
    const transformOpts = {
      transform: opts.transform || 'power',
      power:     (opts.power != null) ? opts.power : 0.5,
      vmax:      (opts.vmax  != null) ? opts.vmax  : 1.0,
    };

    const lowerArr = decodePairsB64(lowerMat.pairs_b64, N * (N - 1) / 2);
    const upperArr = decodePairsB64(upperMat.pairs_b64, N * (N - 1) / 2);

    const img = ctx.createImageData(W, H);

    // For each output pixel, map to (i, j) cell coords.
    // i = floor(y * N / H), j = floor(x * N / W).
    // Diagonal: |i - j| <= 0 (drawn flat at 1 cell).
    for (let y = 0; y < H; y++) {
      const i = Math.floor(y * N / H);
      const baseLowerI = (i * (2 * N - i - 1)) >> 1;
      for (let x = 0; x < W; x++) {
        const j = Math.floor(x * N / W);
        let r2;
        if (i === j) {
          // Diagonal: leave as background; we'll overdraw a hairline if requested
          r2 = NaN;
        } else if (j > i) {
          // Upper triangle: read from upperMat at (i, j)
          r2 = upperArr[baseLowerI + (j - i - 1)] / 255;
        } else {
          // Lower triangle: read from lowerMat at (j, i) (symmetric → swap)
          const baseLowerJ = (j * (2 * N - j - 1)) >> 1;
          r2 = lowerArr[baseLowerJ + (i - j - 1)] / 255;
        }
        const rgba = colorAt(lut, isNaN(r2) ? null : r2, transformOpts);
        const off = (y * W + x) << 2;
        img.data[off]     = rgba[0];
        img.data[off + 1] = rgba[1];
        img.data[off + 2] = rgba[2];
        img.data[off + 3] = rgba[3];
      }
    }

    // White background where alpha=0 (NaN) — replace transparent with white
    // (the canvas may sit on a colored panel and we don't want bleed-through).
    if (opts.bgColor !== false) {
      const [br, bg, bb] = opts.bgColor || [255, 255, 255];
      for (let p = 0; p < img.data.length; p += 4) {
        if (img.data[p + 3] === 0) {
          img.data[p]     = br;
          img.data[p + 1] = bg;
          img.data[p + 2] = bb;
          img.data[p + 3] = 255;
        }
      }
    }

    ctx.putImageData(img, 0, 0);

    // Diagonal hairline
    if (opts.diagColor !== null) {
      ctx.strokeStyle = opts.diagColor || 'rgba(0,0,0,0.35)';
      ctx.lineWidth = 1;
      ctx.beginPath();
      ctx.moveTo(0, 0);
      ctx.lineTo(W, H);
      ctx.stroke();
    }

    return true;
  }

  // ===========================================================================
  // Caption rendering
  // ===========================================================================

  function _f2(x) {
    if (x == null || isNaN(x)) return '—';
    return (Math.round(x * 100) / 100).toFixed(2);
  }
  function _f3(x) {
    if (x == null || isNaN(x)) return '—';
    return (Math.round(x * 1000) / 1000).toFixed(3);
  }

  function _renderCaption(parent, payload) {
    const ta = payload.triangle_assign || {};
    const matrices = payload.matrices || {};
    const lines = [];
    const order = [ta.lower, ta.upper].filter(Boolean);
    // Plus any extra groups not assigned to triangles
    Object.keys(matrices).forEach(function (g) {
      if (!order.includes(g)) order.push(g);
    });
    for (const g of order) {
      const m = matrices[g];
      if (!m) continue;
      const role = (g === ta.lower) ? 'lower'
                 : (g === ta.upper) ? 'upper' : '';
      lines.push(
        '<div class="ld-caption-row">' +
        '<span class="ld-caption-tag" data-role="' + role + '">' +
          (role ? '[' + role + '] ' : '') + g +
        '</span>' +
        ' <span class="ld-caption-n">n=' + m.n_samples + '</span>' +
        ' median r²=' + _f3(m.median_r2_overall) +
        ' · pct r²>0.8=' + _f2(m.pct_pairs_above_0_8 * 100) + '%' +
        (m.shelf_ratio != null && !isNaN(m.shelf_ratio)
          ? ' · shelf-ratio=' + _f2(m.shelf_ratio) + '×'
          : '') +
        '</div>');
    }
    const meta = payload.summary || {};
    lines.push(
      '<div class="ld-caption-meta">' +
      payload.n_snps + ' SNPs · ' + payload.n_pairs.toLocaleString() + ' pairs' +
      (meta.thinning_applied ? ' · thinned from ' + meta.n_snps_unique_in_range : '') +
      ' · ' + (payload.timing && payload.timing.compute_seconds != null
              ? _f3(payload.timing.compute_seconds) + 's compute' : '') +
      ' · ' + (payload.cache_state || 'unknown') +
      '</div>');
    parent.innerHTML = lines.join('');
  }

  // ===========================================================================
  // CSS injection
  // ===========================================================================

  let _cssInjected = false;
  function _injectCSS() {
    if (_cssInjected) return;
    _cssInjected = true;
    if (typeof document === 'undefined') return;
    const css = `
      [data-popgen-ld="root"] {
        font-family: -apple-system, BlinkMacSystemFont, "Segoe UI",
                     Roboto, sans-serif;
        background: #1a1a1a;
        color: #ddd;
        border: 1px solid #333;
        border-radius: 4px;
        padding: 8px;
      }
      [data-popgen-ld="header"] {
        display: flex;
        align-items: center;
        gap: 8px;
        margin-bottom: 6px;
        font-size: 13px;
      }
      [data-popgen-ld="title"] { font-weight: 600; flex: 1; }
      [data-popgen-ld="picker"] {
        background: #222;
        color: #ddd;
        border: 1px solid #555;
        padding: 2px 6px;
      }
      [data-popgen-ld="run"] {
        background: #2a4f7a;
        color: #fff;
        border: none;
        padding: 4px 10px;
        cursor: pointer;
        border-radius: 3px;
      }
      [data-popgen-ld="run"]:disabled { opacity: 0.5; cursor: wait; }
      [data-popgen-ld="status"] {
        font-size: 11px;
        color: #999;
        margin: 4px 0;
        font-family: ui-monospace, monospace;
      }
      [data-popgen-ld="canvas"] {
        background: #fff;
        display: block;
        margin: 4px auto;
        border: 1px solid #444;
      }
      [data-popgen-ld="caption"] {
        font-size: 11px;
        margin-top: 6px;
        font-family: ui-monospace, monospace;
        line-height: 1.45;
      }
      .ld-caption-tag[data-role="lower"] { color: #69b6e6; }
      .ld-caption-tag[data-role="upper"] { color: #f0a070; }
      .ld-caption-meta { color: #777; margin-top: 4px; }
    `;
    const style = document.createElement('style');
    style.textContent = css;
    style.setAttribute('data-popgen-ld-css', '1');
    document.head.appendChild(style);
  }

  // ===========================================================================
  // DOM panel factory
  // ===========================================================================

  function makeLDSplitPanel(opts) {
    opts = opts || {};
    _injectCSS();
    const groupOptions = opts.group_options || ['HOM1', 'HOM2', 'HET'];
    const upperDefault = opts.upper_default || groupOptions[0];
    const canvasSize   = opts.canvas_size || 360;
    const inputProvider = opts.input_provider || _defaultInputProvider;
    const requestFn    = opts.requestFn || ldSplitHeatmap;

    const root = document.createElement('div');
    root.setAttribute('data-popgen-ld', 'root');

    const header = document.createElement('div');
    header.setAttribute('data-popgen-ld', 'header');
    root.appendChild(header);

    const title = document.createElement('div');
    title.setAttribute('data-popgen-ld', 'title');
    title.textContent = 'LD heatmap (split)';
    header.appendChild(title);

    const picker = document.createElement('select');
    picker.setAttribute('data-popgen-ld', 'picker');
    for (const g of groupOptions) {
      const o = document.createElement('option');
      o.value = g; o.textContent = g;
      if (g === upperDefault) o.selected = true;
      picker.appendChild(o);
    }
    header.appendChild(picker);

    const runBtn = document.createElement('button');
    runBtn.setAttribute('data-popgen-ld', 'run');
    runBtn.textContent = '▶ run';
    header.appendChild(runBtn);

    const status = document.createElement('div');
    status.setAttribute('data-popgen-ld', 'status');
    status.textContent = 'idle';
    root.appendChild(status);

    const canvas = document.createElement('canvas');
    canvas.setAttribute('data-popgen-ld', 'canvas');
    canvas.width = canvasSize;
    canvas.height = canvasSize;
    root.appendChild(canvas);

    const caption = document.createElement('div');
    caption.setAttribute('data-popgen-ld', 'caption');
    root.appendChild(caption);

    let lastPayload = null;
    let inFlight = false;

    async function run() {
      if (inFlight) return;
      let req;
      try {
        req = inputProvider(picker.value);
      } catch (e) {
        status.textContent = 'cannot build request: ' + (e.message || e);
        return;
      }
      if (!req) {
        status.textContent = 'no candidate selected';
        return;
      }
      inFlight = true;
      runBtn.disabled = true;
      status.textContent = 'computing… (window range: ' +
        req.window_range.join('-') + ')';
      const t0 = Date.now();
      try {
        const payload = await requestFn(req);
        lastPayload = payload;
        const ms = Date.now() - t0;
        const cache = payload.cache_state || 'unknown';
        status.textContent =
          'done in ' + ms + 'ms · ' + cache +
          (payload.cache_key ? ' · key=' + payload.cache_key.slice(0, 24) : '');
        drawSplitHeatmap(canvas, payload, opts.draw || {});
        _renderCaption(caption, payload);
      } catch (e) {
        status.textContent = 'error: ' + (e.message || e);
      } finally {
        inFlight = false;
        runBtn.disabled = false;
      }
    }

    runBtn.addEventListener('click', run);
    picker.addEventListener('change', function () {
      // Mark picker change without auto-running; let user click ▶
      status.textContent = 'picker changed → click ▶ to recompute';
    });

    // Programmatic API on the panel
    root.popgenLDPanel = {
      run: run,
      getPayload: function () { return lastPayload; },
      setUpperGroup: function (g) {
        for (const o of picker.options) o.selected = (o.value === g);
      },
    };

    return root;
  }

  // ===========================================================================
  // Default input provider — looks at window.atlasState if present
  // ===========================================================================

  function _defaultInputProvider(upperGroupName) {
    if (typeof window === 'undefined') return null;
    const state = window.atlasState || {};
    const cand = state.candidate;
    if (!cand) return null;
    const fish = cand.fish_calls || {};
    // Build sample lists from fish_calls regimes:
    //   regime 0 → HOM1, 1 → HET, 2 → HOM2
    const regimeMap = { HOM1: 0, HET: 1, HOM2: 2 };
    const targetRegime = regimeMap[upperGroupName];
    if (targetRegime == null) return null;
    const allSamples = state.allSamples || [];
    const upperSamples = [];
    for (const s of allSamples) {
      const call = fish[s];
      if (call != null && call.regime === targetRegime) upperSamples.push(s);
    }
    if (upperSamples.length < 5) return null;
    return {
      chrom: cand.chrom,
      window_range: [cand.start_w, cand.end_w],
      groups: {
        ALL: allSamples,
        [upperGroupName]: upperSamples,
      },
      triangle_assign: { lower: 'ALL', upper: upperGroupName },
      shelf_bp: cand.shelf_start_bp != null
        ? [cand.shelf_start_bp, cand.shelf_end_bp]
        : null,
      snp_cap: 5000,
    };
  }

  // ===========================================================================
  // Request wrapper — delegates to popgenLive, falls back to direct fetch
  // ===========================================================================

  async function ldSplitHeatmap(reqBody, opts) {
    opts = opts || {};
    const url = opts.url || '/api/ld/split_heatmap';

    // Prefer the request-layer adapter from turn 3 if it exposes the kind
    if (typeof window !== 'undefined' &&
        window.popgenLive &&
        typeof window.popgenLive.ldSplitHeatmap === 'function') {
      return window.popgenLive.ldSplitHeatmap(reqBody, opts);
    }

    // Direct fetch fallback
    const resp = await fetch(url, {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify(reqBody),
      signal: opts.signal,
    });
    if (!resp.ok) {
      const text = await resp.text();
      throw new Error('LD request failed (' + resp.status + '): ' + text);
    }
    return resp.json();
  }

  // ===========================================================================
  // Public surface
  // ===========================================================================

  return {
    APEX_COLORS: APEX_COLORS.slice(),
    apexColormap: apexColormap,
    colorAt: colorAt,
    decodePairsB64: decodePairsB64,
    pairAt: pairAt,
    drawSplitHeatmap: drawSplitHeatmap,
    makeLDSplitPanel: makeLDSplitPanel,
    ldSplitHeatmap: ldSplitHeatmap,
    _renderCaption: _renderCaption,
    _defaultInputProvider: _defaultInputProvider,
  };
}));
