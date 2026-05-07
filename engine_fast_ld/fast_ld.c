/*
 * fast_ld.c — windows-driven LD engine for inversion heatmaps
 * ============================================================
 *
 * Computes pairwise r² (Pearson correlation²) on dosage data, given a
 * pre-built windows JSON that maps window_idx → list of SNP bp positions.
 *
 * Method
 * ------
 *   Input:   one control file specifying paths, window range, groups,
 *            and a SNP cap.
 *   Step 1:  Load windows JSON; collect unique SNP positions across the
 *            requested window range.
 *   Step 2:  If unique_snp_count > snp_cap, error out (unless thin_to is
 *            set, in which case even-bp-thin to thin_to and continue).
 *   Step 3:  Stream the dosage TSV; keep only rows whose pos matches a
 *            requested SNP position.
 *   Step 4:  Per group: build n_snps × n_group dosage submatrix,
 *            standardize each row (mean 0, unit variance), compute
 *            pairwise r² as (Z @ Z') / k. Output upper triangle as uint8
 *            (r²·255 quantized).
 *   Step 5:  Emit summary statistics (median r², shelf-ratio if shelf
 *            bounds set, decay deciles, pct above 0.8).
 *
 * The output heatmap is N × N where N = number of unique SNPs from the
 * requested windows. Every cell is one real SNP-pair r² (no binning).
 *
 * Control file format (key=value, one per line, # for comments):
 *
 *     dosage_path=/path/to/<chrom>.dosage.tsv.gz
 *     windows_json=/path/to/<chrom>.windows.json
 *     window_range=2150-2400         # inclusive on both ends
 *     snp_cap=5000                   # default: 5000
 *     thin_to=                       # optional; if set, thin to N SNPs
 *                                    # (overrides snp_cap)
 *     shelf_start=15200000           # optional; for shelf-ratio summary
 *     shelf_end=18100000             # optional
 *     group=NAME:id1,id2,...         # one line per group, max 4
 *     out_dir=/path/to/output/
 *     threads=8                      # default: 1
 *
 * Outputs in out_dir/:
 *
 *     pairs.<groupname>.bin    uint8 r²·255, upper triangle row-major,
 *                              N*(N-1)/2 bytes
 *     sites.bin                packed struct, 56 bytes per kept SNP
 *     summary.tsv              key\tvalue, one per line
 *
 * Build:
 *     gcc -O3 -march=native -fopenmp -o fast_ld fast_ld.c -lz -lm
 */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
#include <time.h>
#include <ctype.h>
#include <errno.h>
#include <stdarg.h>
#include <zlib.h>
#include <sys/stat.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#define MAX_PATH 4096
#define MAX_LINE (1 << 24)        /* 16 MiB */
#define MAX_SAMPLES 4096
#define MAX_GROUPS 4
#define MAX_GROUP_NAME 64
#define DEFAULT_SNP_CAP 5000

/* =========================================================================
 * Utilities
 * =======================================================================*/

static double now_s(void) {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec + ts.tv_nsec * 1e-9;
}

__attribute__((noreturn))
static void diefmt(const char *fmt, ...) {
    va_list ap;
    va_start(ap, fmt);
    fprintf(stderr, "fast_ld: ");
    vfprintf(stderr, fmt, ap);
    fprintf(stderr, "\n");
    va_end(ap);
    exit(1);
}

static char *strdup_or_die(const char *s) {
    char *out = strdup(s);
    if (!out) diefmt("out of memory in strdup");
    return out;
}

static void *xmalloc(size_t n) {
    void *p = malloc(n);
    if (!p && n > 0) diefmt("out of memory (%zu bytes)", n);
    return p;
}

static void *xrealloc(void *p, size_t n) {
    void *q = realloc(p, n);
    if (!q && n > 0) diefmt("out of memory (realloc %zu bytes)", n);
    return q;
}

static void trim(char *s) {
    char *p = s;
    while (*p && isspace((unsigned char)*p)) p++;
    if (p != s) memmove(s, p, strlen(p) + 1);
    size_t n = strlen(s);
    while (n > 0 && isspace((unsigned char)s[n - 1])) s[--n] = 0;
}

static int split_tsv(char *buf, char **fields, int max_fields) {
    int n = 0;
    fields[n++] = buf;
    for (char *p = buf; *p; p++) {
        if (*p == '\n' || *p == '\r') { *p = 0; break; }
        if (*p == '\t') {
            *p = 0;
            if (n >= max_fields) return -1;
            fields[n++] = p + 1;
        }
    }
    return n;
}

/* =========================================================================
 * Control file
 * =======================================================================*/

typedef struct {
    char dosage_path[MAX_PATH];
    char windows_json[MAX_PATH];
    int  win_lo;
    int  win_hi;
    int  snp_cap;
    int  thin_to;            /* 0 = unset */
    long shelf_start;        /* -1 = unset */
    long shelf_end;
    char out_dir[MAX_PATH];
    int  threads;

    int  n_groups;
    char group_names[MAX_GROUPS][MAX_GROUP_NAME];
    char *group_ids_csv[MAX_GROUPS];
} Control;

static void read_control(const char *path, Control *c) {
    memset(c, 0, sizeof(*c));
    c->win_lo = -1;
    c->win_hi = -1;
    c->shelf_start = -1;
    c->shelf_end = -1;
    c->snp_cap = DEFAULT_SNP_CAP;
    c->threads = 1;

    FILE *fp = fopen(path, "r");
    if (!fp) diefmt("cannot open control file: %s", path);
    char line[8192];
    int line_no = 0;
    while (fgets(line, sizeof(line), fp)) {
        line_no++;
        trim(line);
        if (line[0] == 0 || line[0] == '#') continue;
        char *eq = strchr(line, '=');
        if (!eq) diefmt("control line %d: missing '='", line_no);
        *eq = 0;
        char *key = line, *val = eq + 1;
        trim(key); trim(val);

        if (strcmp(key, "dosage_path") == 0) {
            strncpy(c->dosage_path, val, MAX_PATH - 1);
        } else if (strcmp(key, "windows_json") == 0) {
            strncpy(c->windows_json, val, MAX_PATH - 1);
        } else if (strcmp(key, "window_range") == 0) {
            char *dash = strchr(val, '-');
            if (!dash) diefmt("window_range must be 'lo-hi': %s", val);
            *dash = 0;
            c->win_lo = atoi(val);
            c->win_hi = atoi(dash + 1);
            if (c->win_lo < 0 || c->win_hi < c->win_lo)
                diefmt("invalid window_range: %d-%d", c->win_lo, c->win_hi);
        } else if (strcmp(key, "snp_cap") == 0) {
            c->snp_cap = atoi(val);
        } else if (strcmp(key, "thin_to") == 0) {
            if (val[0]) c->thin_to = atoi(val);
        } else if (strcmp(key, "shelf_start") == 0) {
            c->shelf_start = atol(val);
        } else if (strcmp(key, "shelf_end") == 0) {
            c->shelf_end = atol(val);
        } else if (strcmp(key, "out_dir") == 0) {
            strncpy(c->out_dir, val, MAX_PATH - 1);
        } else if (strcmp(key, "threads") == 0) {
            c->threads = atoi(val);
        } else if (strcmp(key, "group") == 0) {
            if (c->n_groups >= MAX_GROUPS)
                diefmt("too many groups (max %d)", MAX_GROUPS);
            char *colon = strchr(val, ':');
            if (!colon) diefmt("group line missing ':' separator: %s", val);
            *colon = 0;
            strncpy(c->group_names[c->n_groups], val, MAX_GROUP_NAME - 1);
            c->group_ids_csv[c->n_groups] = strdup_or_die(colon + 1);
            c->n_groups++;
        } else {
            fprintf(stderr, "fast_ld: warning: unknown key '%s' (line %d)\n",
                    key, line_no);
        }
    }
    fclose(fp);

    if (!c->dosage_path[0]) diefmt("control: dosage_path required");
    if (!c->windows_json[0]) diefmt("control: windows_json required");
    if (c->win_lo < 0) diefmt("control: window_range required");
    if (!c->out_dir[0]) diefmt("control: out_dir required");
    if (c->n_groups < 1) diefmt("control: at least one group required");
    if (c->snp_cap < 10) diefmt("control: snp_cap too small (min 10)");
}

/* =========================================================================
 * Minimal JSON parser for the windows file
 *
 * We only need to extract:
 *   chrom (string at top level)
 *   windows[] with each containing:
 *     idx (int)
 *     start_bp, end_bp (int)
 *     snp_positions[] (int array)
 *
 * Rather than pull in cJSON, we hand-roll a parser that handles exactly
 * this schema. Tolerates whitespace, requires the schema we wrote.
 * =======================================================================*/

typedef struct {
    const char *p;
    const char *end;
} JParser;

static void j_skip_ws(JParser *j) {
    while (j->p < j->end && isspace((unsigned char)*j->p)) j->p++;
}
static int j_eat(JParser *j, char c) {
    j_skip_ws(j);
    if (j->p < j->end && *j->p == c) { j->p++; return 1; }
    return 0;
}
static long j_int(JParser *j) {
    j_skip_ws(j);
    char *endp;
    long v = strtol(j->p, &endp, 10);
    if (endp == j->p) diefmt("JSON parse: expected integer near '%.20s'", j->p);
    j->p = endp;
    return v;
}
/* Read a string literal "..." into out (may grow). Returns malloc'd string. */
static char *j_str(JParser *j) {
    j_skip_ws(j);
    if (!j_eat(j, '"')) diefmt("JSON parse: expected '\"' near '%.20s'", j->p);
    const char *s = j->p;
    while (j->p < j->end && *j->p != '"') {
        if (*j->p == '\\') {
            if (j->p + 1 < j->end) j->p += 2;
            else diefmt("JSON parse: bad escape");
        } else j->p++;
    }
    if (j->p >= j->end) diefmt("JSON parse: unterminated string");
    size_t len = j->p - s;
    char *out = xmalloc(len + 1);
    memcpy(out, s, len);
    out[len] = 0;
    j->p++;  /* skip closing quote */
    return out;
}
/* Skip a JSON value (object, array, string, number, bool, null). */
static void j_skip_value(JParser *j) {
    j_skip_ws(j);
    if (j->p >= j->end) diefmt("JSON parse: unexpected end");
    if (*j->p == '"') { free(j_str(j)); return; }
    if (*j->p == '{' || *j->p == '[') {
        char open = *j->p, close = (open == '{' ? '}' : ']');
        j->p++;
        int depth = 1;
        while (j->p < j->end && depth > 0) {
            if (*j->p == '"') { free(j_str(j)); continue; }
            if (*j->p == open) depth++;
            else if (*j->p == close) depth--;
            j->p++;
        }
        return;
    }
    /* number / bool / null — skip until comma, } or ] */
    while (j->p < j->end && *j->p != ',' && *j->p != '}' && *j->p != ']' &&
           !isspace((unsigned char)*j->p)) j->p++;
}

/* Window record */
typedef struct {
    int idx;
    long start_bp;
    long end_bp;
    int n_snps;
    long *snp_positions;   /* malloc'd */
} WindowRec;

typedef struct {
    char *chrom;
    int n_windows;
    WindowRec *windows;
} WindowsJSON;

static void parse_windows_json(const char *path, WindowsJSON *out) {
    FILE *fp = fopen(path, "rb");
    if (!fp) diefmt("cannot open windows JSON: %s", path);
    fseek(fp, 0, SEEK_END);
    long sz = ftell(fp);
    fseek(fp, 0, SEEK_SET);
    char *buf = xmalloc(sz + 1);
    if (fread(buf, 1, sz, fp) != (size_t)sz) diefmt("short read of %s", path);
    buf[sz] = 0;
    fclose(fp);

    JParser j = { buf, buf + sz };
    if (!j_eat(&j, '{')) diefmt("JSON parse: expected '{' at top level");

    memset(out, 0, sizeof(*out));
    /* Walk top-level keys; pluck chrom + windows; skip everything else */
    int first_key = 1;
    while (1) {
        j_skip_ws(&j);
        if (j_eat(&j, '}')) break;
        if (!first_key) {
            if (!j_eat(&j, ',')) diefmt("JSON parse: ',' between top-level keys");
            j_skip_ws(&j);
        }
        first_key = 0;
        char *key = j_str(&j);
        if (!j_eat(&j, ':')) diefmt("JSON parse: expected ':' after key '%s'", key);
        if (strcmp(key, "chrom") == 0) {
            out->chrom = j_str(&j);
        } else if (strcmp(key, "windows") == 0) {
            j_skip_ws(&j);
            if (!j_eat(&j, '[')) diefmt("JSON: 'windows' must be array");
            size_t cap = 4096;
            out->windows = xmalloc(sizeof(WindowRec) * cap);
            out->n_windows = 0;
            int first_win = 1;
            while (1) {
                j_skip_ws(&j);
                if (j_eat(&j, ']')) break;
                if (!first_win) {
                    if (!j_eat(&j, ',')) diefmt("JSON: ',' between windows");
                }
                first_win = 0;
                /* Parse one window object */
                j_skip_ws(&j);
                if (!j_eat(&j, '{')) diefmt("JSON: '{' for window");
                if ((size_t)out->n_windows >= cap) {
                    cap *= 2;
                    out->windows = xrealloc(out->windows, sizeof(WindowRec) * cap);
                }
                WindowRec *w = &out->windows[out->n_windows];
                memset(w, 0, sizeof(*w));
                int first_field = 1;
                while (1) {
                    j_skip_ws(&j);
                    if (j_eat(&j, '}')) break;
                    if (!first_field) {
                        if (!j_eat(&j, ',')) diefmt("JSON: ',' between window fields");
                        j_skip_ws(&j);
                    }
                    first_field = 0;
                    char *wk = j_str(&j);
                    if (!j_eat(&j, ':')) diefmt("JSON: ':' after window key");
                    if (strcmp(wk, "idx") == 0) {
                        w->idx = (int)j_int(&j);
                    } else if (strcmp(wk, "start_bp") == 0) {
                        w->start_bp = j_int(&j);
                    } else if (strcmp(wk, "end_bp") == 0) {
                        w->end_bp = j_int(&j);
                    } else if (strcmp(wk, "n_snps") == 0) {
                        w->n_snps = (int)j_int(&j);
                    } else if (strcmp(wk, "snp_positions") == 0) {
                        j_skip_ws(&j);
                        if (!j_eat(&j, '[')) diefmt("JSON: '[' for snp_positions");
                        size_t pcap = w->n_snps > 0 ? (size_t)w->n_snps : 64;
                        w->snp_positions = xmalloc(sizeof(long) * pcap);
                        int pn = 0;
                        int first_p = 1;
                        while (1) {
                            j_skip_ws(&j);
                            if (j_eat(&j, ']')) break;
                            if (!first_p) {
                                if (!j_eat(&j, ',')) diefmt("JSON: ',' in snp_positions");
                            }
                            first_p = 0;
                            if ((size_t)pn >= pcap) {
                                pcap *= 2;
                                w->snp_positions = xrealloc(w->snp_positions,
                                                            sizeof(long) * pcap);
                            }
                            w->snp_positions[pn++] = j_int(&j);
                        }
                        w->n_snps = pn;
                    } else {
                        j_skip_value(&j);
                    }
                    free(wk);
                }
                out->n_windows++;
            }
        } else {
            j_skip_value(&j);
        }
        free(key);
    }
    free(buf);
    if (!out->chrom) diefmt("windows JSON missing 'chrom'");
    if (!out->windows) diefmt("windows JSON missing 'windows'");
}

static void free_windows_json(WindowsJSON *w) {
    if (w->chrom) free(w->chrom);
    for (int i = 0; i < w->n_windows; i++) {
        free(w->windows[i].snp_positions);
    }
    free(w->windows);
    memset(w, 0, sizeof(*w));
}

/* =========================================================================
 * Build the unique SNP-position list for the requested window range
 * =======================================================================*/

static int cmp_long(const void *a, const void *b) {
    long x = *(const long*)a, y = *(const long*)b;
    return (x > y) - (x < y);
}

static long *unique_positions_in_range(const WindowsJSON *wj, int win_lo, int win_hi,
                                       int *n_out)
{
    if (win_hi >= wj->n_windows) win_hi = wj->n_windows - 1;
    if (win_lo > win_hi) diefmt("empty window range %d-%d", win_lo, win_hi);
    /* Total membership count (with duplicates) */
    long total = 0;
    for (int w = win_lo; w <= win_hi; w++) total += wj->windows[w].n_snps;
    long *all = xmalloc(sizeof(long) * total);
    long k = 0;
    for (int w = win_lo; w <= win_hi; w++) {
        WindowRec *wr = &wj->windows[w];
        for (int i = 0; i < wr->n_snps; i++) all[k++] = wr->snp_positions[i];
    }
    qsort(all, total, sizeof(long), cmp_long);
    /* Dedupe in place */
    long n_unique = 0;
    for (long i = 0; i < total; i++) {
        if (i == 0 || all[i] != all[i-1]) all[n_unique++] = all[i];
    }
    *n_out = (int)n_unique;
    return all;  /* length = n_unique */
}

/* Even-bp thinning: choose `target` positions evenly spaced across the
 * sorted input array. */
static long *thin_positions(const long *positions, int n_in, int target,
                            int *n_out)
{
    if (n_in <= target) {
        long *out = xmalloc(sizeof(long) * n_in);
        memcpy(out, positions, sizeof(long) * n_in);
        *n_out = n_in;
        return out;
    }
    long *out = xmalloc(sizeof(long) * target);
    long lo = positions[0], hi = positions[n_in - 1];
    double range = (double)(hi - lo);
    int j = 0, n_kept = 0;
    for (int t = 0; t < target; t++) {
        double ideal = lo + (range * (t + 0.5)) / target;
        while (j + 1 < n_in &&
               fabs((double)positions[j+1] - ideal) <
               fabs((double)positions[j] - ideal)) j++;
        if (n_kept == 0 || out[n_kept - 1] != positions[j]) out[n_kept++] = positions[j];
    }
    *n_out = n_kept;
    return out;
}

/* =========================================================================
 * Dosage TSV reader — but selective: keep only rows whose pos is in
 * `wanted_positions` (sorted, deduped). Uses binary search per row.
 * =======================================================================*/

typedef struct {
    int n_snps;          /* number of kept SNPs */
    int n_samples;
    char **sample_ids;   /* length n_samples */
    long *positions;     /* length n_snps */
    float *dosages;      /* length n_snps * n_samples, row-major */
} DosageMatrix;

static void free_matrix(DosageMatrix *m) {
    if (m->sample_ids) {
        for (int i = 0; i < m->n_samples; i++) free(m->sample_ids[i]);
        free(m->sample_ids);
    }
    free(m->positions);
    free(m->dosages);
    memset(m, 0, sizeof(*m));
}

/* Parse "<chrom>_<pos>" or "<chrom>:<pos>"; returns 0 on success. */
static int parse_marker(const char *marker, const char *expect_chrom,
                        long *pos_out)
{
    size_t n = strlen(marker);
    if (n == 0) return -1;
    ssize_t i = (ssize_t)n - 1;
    while (i >= 0 && isdigit((unsigned char)marker[i])) i--;
    if (i < 0 || i == (ssize_t)n - 1) return -1;
    if (marker[i] != '_' && marker[i] != ':') return -1;
    char *endptr;
    long pos = strtol(marker + i + 1, &endptr, 10);
    if (*endptr != 0) return -1;
    if (strncmp(marker, expect_chrom, (size_t)i) != 0) return -1;
    if (strlen(expect_chrom) != (size_t)i) return -1;
    *pos_out = pos;
    return 0;
}

/* Binary search inside read_dosage_selective uses inline code, no helper. */

static void read_dosage_selective(const Control *c, const char *expect_chrom,
                                  const long *wanted_positions, int n_wanted,
                                  DosageMatrix *m)
{
    gzFile gz = gzopen(c->dosage_path, "rb");
    if (!gz) diefmt("cannot open dosage file: %s", c->dosage_path);
    char *line = xmalloc(MAX_LINE);

    if (gzgets(gz, line, MAX_LINE) == NULL) diefmt("dosage file empty");

    /* Header */
    int n_tabs = 0;
    for (char *p = line; *p; p++) if (*p == '\t') n_tabs++;
    int n_samples = n_tabs;
    if (n_samples < 1 || n_samples > MAX_SAMPLES)
        diefmt("invalid sample count from header: %d", n_samples);

    char **fields = xmalloc(sizeof(char*) * (n_samples + 16));
    int nf = split_tsv(line, fields, n_samples + 16);
    if (nf != n_samples + 1) diefmt("header field count mismatch");

    m->sample_ids = xmalloc(sizeof(char*) * n_samples);
    for (int i = 0; i < n_samples; i++) m->sample_ids[i] = strdup_or_die(fields[i + 1]);
    m->n_samples = n_samples;

    /* Allocate exactly n_wanted rows; we fill them by position-index */
    m->n_snps = n_wanted;
    m->positions = xmalloc(sizeof(long) * n_wanted);
    memcpy(m->positions, wanted_positions, sizeof(long) * n_wanted);
    m->dosages = xmalloc(sizeof(float) * (size_t)n_wanted * n_samples);
    /* Mark all rows as not-yet-filled by setting first column to NaN; we'll
     * verify all were filled at the end. */
    for (long k = 0; k < (long)n_wanted * n_samples; k++) m->dosages[k] = NAN;
    int n_filled = 0;

    long range_lo = wanted_positions[0];
    long range_hi = wanted_positions[n_wanted - 1];

    int line_no = 1;
    int n_chrom_mismatch = 0;
    while (gzgets(gz, line, MAX_LINE)) {
        line_no++;
        nf = split_tsv(line, fields, n_samples + 16);
        if (nf != n_samples + 1)
            diefmt("row %d: %d fields != %d expected", line_no, nf, n_samples + 1);
        long pos;
        if (parse_marker(fields[0], expect_chrom, &pos) != 0) {
            n_chrom_mismatch++;
            continue;
        }
        if (pos < range_lo) continue;
        if (pos > range_hi) break;   /* dosage rows are sorted by pos */
        /* Binary search for `pos` in wanted_positions */
        int lo = 0, hi = n_wanted - 1;
        int found_idx = -1;
        while (lo <= hi) {
            int mid = (lo + hi) >> 1;
            if (wanted_positions[mid] == pos) { found_idx = mid; break; }
            if (wanted_positions[mid] < pos) lo = mid + 1;
            else hi = mid - 1;
        }
        if (found_idx < 0) continue;
        float *row = m->dosages + (size_t)found_idx * n_samples;
        for (int j = 0; j < n_samples; j++) {
            char *endptr;
            float v = strtof(fields[j + 1], &endptr);
            if (*endptr != 0 && !isspace((unsigned char)*endptr)) v = NAN;
            row[j] = v;
        }
        n_filled++;
    }
    gzclose(gz);
    free(fields);
    free(line);

    fprintf(stderr, "fast_ld: dosage selective read: %d/%d wanted SNPs found "
                    "(chrom mismatch: %d)\n",
            n_filled, n_wanted, n_chrom_mismatch);
    if (n_filled < n_wanted) {
        /* Compact the matrix to drop unfilled rows */
        int k = 0;
        for (int i = 0; i < n_wanted; i++) {
            if (!isnan(m->dosages[(size_t)i * n_samples])) {
                if (k != i) {
                    m->positions[k] = m->positions[i];
                    memmove(m->dosages + (size_t)k * n_samples,
                            m->dosages + (size_t)i * n_samples,
                            sizeof(float) * n_samples);
                }
                k++;
            }
        }
        m->n_snps = k;
        if (k < 2) diefmt("only %d SNPs from requested set found in dosage", k);
        fprintf(stderr, "fast_ld: compacted to %d SNPs\n", k);
    }
}

/* =========================================================================
 * Group resolution
 * =======================================================================*/

typedef struct {
    char name[MAX_GROUP_NAME];
    int n;
    int *col_idx;
} GroupResolved;

static int find_sample_col(const DosageMatrix *m, const char *id) {
    for (int i = 0; i < m->n_samples; i++)
        if (strcmp(m->sample_ids[i], id) == 0) return i;
    return -1;
}

static void resolve_group(const DosageMatrix *m, const char *name,
                          const char *ids_csv, GroupResolved *g)
{
    strncpy(g->name, name, MAX_GROUP_NAME - 1);
    int n_commas = 0;
    for (const char *p = ids_csv; *p; p++) if (*p == ',') n_commas++;
    g->col_idx = xmalloc(sizeof(int) * (n_commas + 1));
    g->n = 0;
    char *buf = strdup_or_die(ids_csv);
    char *tok = strtok(buf, ",");
    int n_unknown = 0;
    while (tok) {
        trim(tok);
        if (*tok) {
            int col = find_sample_col(m, tok);
            if (col < 0) {
                n_unknown++;
                if (n_unknown <= 5)
                    fprintf(stderr, "fast_ld: warn: unknown sample in '%s': %s\n",
                            name, tok);
            } else {
                g->col_idx[g->n++] = col;
            }
        }
        tok = strtok(NULL, ",");
    }
    free(buf);
    if (n_unknown > 5)
        fprintf(stderr, "fast_ld: ... %d more unknown samples in '%s'\n",
                n_unknown - 5, name);
    if (g->n < 5) diefmt("group '%s' has only %d valid samples (min 5)", name, g->n);
}

/* =========================================================================
 * r² computation
 * =======================================================================*/

typedef struct {
    int n_snps;
    int n_group_samples;
    long *positions;       /* shared with caller */
    long n_pairs;
    uint8_t *r2_q;
    float *maf;
    float *var;
    int *n_complete;
    /* summary */
    double median_r2_overall;
    double pct_above_0_8;
    double median_r2_shelf;
    double median_r2_flank;
    double shelf_ratio;
    double decay_deciles[10];
    double compute_seconds;
} GroupResult;

static void build_submatrix(const DosageMatrix *m, const GroupResolved *g,
                            float *out)
{
    int k_n = g->n;
    for (int i = 0; i < m->n_snps; i++) {
        const float *src = m->dosages + (size_t)i * m->n_samples;
        float *dst = out + (size_t)i * k_n;
        for (int k = 0; k < k_n; k++) dst[k] = src[g->col_idx[k]];
    }
}

static void standardize_rows(float *mat, int n_rows, int n_cols,
                             float *out_maf, float *out_var, int *out_nc)
{
    #pragma omp parallel for schedule(static)
    for (int i = 0; i < n_rows; i++) {
        float *row = mat + (size_t)i * n_cols;
        double sum = 0.0;
        int n_ok = 0;
        for (int k = 0; k < n_cols; k++) {
            float v = row[k];
            if (!isnan(v)) { sum += v; n_ok++; }
        }
        out_nc[i] = n_ok;
        if (n_ok < 2) {
            out_maf[i] = NAN;
            out_var[i] = 0;
            for (int k = 0; k < n_cols; k++) row[k] = 0.0f;
            continue;
        }
        double mean = sum / n_ok;
        double ss = 0;
        for (int k = 0; k < n_cols; k++) {
            float v = row[k];
            if (!isnan(v)) { double d = v - mean; ss += d * d; }
        }
        double var = ss / n_ok;
        out_maf[i] = (float)(mean * 0.5);
        out_var[i] = (float)var;
        if (var <= 1e-12) {
            for (int k = 0; k < n_cols; k++) row[k] = 0.0f;
            continue;
        }
        double sd = sqrt(var);
        for (int k = 0; k < n_cols; k++) {
            float v = row[k];
            row[k] = isnan(v) ? 0.0f : (float)((v - mean) / sd);
        }
    }
}

static void compute_r2(const float *Z, int n_rows, int n_cols,
                       const float *var, uint8_t *r2_q)
{
    double inv_k = 1.0 / (double)n_cols;
    #pragma omp parallel for schedule(dynamic, 16)
    for (int i = 0; i < n_rows; i++) {
        const float *zi = Z + (size_t)i * n_cols;
        long base = (long)i * (2L * n_rows - i - 1) / 2;
        for (int j = i + 1; j < n_rows; j++) {
            const float *zj = Z + (size_t)j * n_cols;
            float r2;
            if (var[i] <= 1e-12f || var[j] <= 1e-12f) {
                r2 = 0.0f;
            } else {
                double dot = 0.0;
                for (int k = 0; k < n_cols; k++)
                    dot += (double)zi[k] * (double)zj[k];
                double r = dot * inv_k;
                if (r >  1.0) r =  1.0;
                if (r < -1.0) r = -1.0;
                r2 = (float)(r * r);
            }
            int q = (int)(r2 * 255.0f + 0.5f);
            if (q < 0) q = 0;
            if (q > 255) q = 255;
            r2_q[base + (j - i - 1)] = (uint8_t)q;
        }
    }
}

/* =========================================================================
 * Summary stats
 * =======================================================================*/

static int dbl_cmp(const void *a, const void *b) {
    double x = *(const double*)a, y = *(const double*)b;
    return (x > y) - (x < y);
}

static double median_of_q(const uint8_t *q, const int *idx, int n) {
    if (n == 0) return NAN;
    double *buf = xmalloc(sizeof(double) * n);
    for (int i = 0; i < n; i++) buf[i] = q[idx[i]] / 255.0;
    qsort(buf, n, sizeof(double), dbl_cmp);
    double med = (n % 2) ? buf[n / 2] : 0.5 * (buf[n/2 - 1] + buf[n/2]);
    free(buf);
    return med;
}

static void compute_summary(GroupResult *r, const Control *c)
{
    long n_pairs = r->n_pairs;
    int n = r->n_snps;

    /* Median r² and pct above 0.8 across ALL pairs */
    if (n_pairs > 0) {
        double *all = xmalloc(sizeof(double) * n_pairs);
        int n_above = 0;
        uint8_t thr = (uint8_t)(0.8 * 255);
        for (long p = 0; p < n_pairs; p++) {
            all[p] = r->r2_q[p] / 255.0;
            if (r->r2_q[p] > thr) n_above++;
        }
        qsort(all, n_pairs, sizeof(double), dbl_cmp);
        r->median_r2_overall = (n_pairs % 2) ? all[n_pairs/2]
            : 0.5 * (all[n_pairs/2 - 1] + all[n_pairs/2]);
        r->pct_above_0_8 = (double)n_above / n_pairs;
        free(all);
    } else {
        r->median_r2_overall = NAN;
        r->pct_above_0_8 = 0;
    }

    /* Decay deciles */
    double dist_min = 1e18, dist_max = 0;
    for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++) {
        double d = (double)(r->positions[j] - r->positions[i]);
        if (d < dist_min) dist_min = d;
        if (d > dist_max) dist_max = d;
    }
    double dec_sum[10] = {0};
    long dec_n[10] = {0};
    double range = dist_max - dist_min;
    if (range < 1) range = 1;
    long pair_idx = 0;
    for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++) {
        double d = (double)(r->positions[j] - r->positions[i]);
        int bin = (int)(((d - dist_min) / range) * 10.0);
        if (bin < 0) bin = 0;
        if (bin > 9) bin = 9;
        dec_sum[bin] += r->r2_q[pair_idx] / 255.0;
        dec_n[bin]++;
        pair_idx++;
    }
    for (int b = 0; b < 10; b++)
        r->decay_deciles[b] = dec_n[b] > 0 ? dec_sum[b] / dec_n[b] : NAN;

    /* Shelf ratio */
    if (c->shelf_start > 0 && c->shelf_end > 0) {
        int *shelf_p = xmalloc(sizeof(int) * n_pairs);
        int *flank_p = xmalloc(sizeof(int) * n_pairs);
        int n_shelf = 0, n_flank = 0;
        pair_idx = 0;
        for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++) {
            int i_in = (r->positions[i] >= c->shelf_start &&
                        r->positions[i] <= c->shelf_end);
            int j_in = (r->positions[j] >= c->shelf_start &&
                        r->positions[j] <= c->shelf_end);
            if (i_in && j_in) shelf_p[n_shelf++] = pair_idx;
            else if (!i_in && !j_in) flank_p[n_flank++] = pair_idx;
            pair_idx++;
        }
        r->median_r2_shelf = median_of_q(r->r2_q, shelf_p, n_shelf);
        r->median_r2_flank = median_of_q(r->r2_q, flank_p, n_flank);
        r->shelf_ratio = (r->median_r2_flank > 0)
            ? r->median_r2_shelf / r->median_r2_flank : NAN;
        free(shelf_p); free(flank_p);
    } else {
        r->median_r2_shelf = NAN;
        r->median_r2_flank = NAN;
        r->shelf_ratio = NAN;
    }
}

/* =========================================================================
 * Output writers (same schema as before)
 * =======================================================================*/

static void write_pairs_bin(const char *out_dir, const char *gname,
                            const uint8_t *r2_q, long n_pairs)
{
    char path[MAX_PATH];
    snprintf(path, MAX_PATH, "%s/pairs.%s.bin", out_dir, gname);
    FILE *fp = fopen(path, "wb");
    if (!fp) diefmt("cannot write %s", path);
    if (fwrite(r2_q, 1, n_pairs, fp) != (size_t)n_pairs) diefmt("short write");
    fclose(fp);
}

static void write_sites_bin(const char *out_dir, int n_snps,
                            const long *positions,
                            const GroupResult *results, int n_groups)
{
    char path[MAX_PATH];
    snprintf(path, MAX_PATH, "%s/sites.bin", out_dir);
    FILE *fp = fopen(path, "wb");
    if (!fp) diefmt("cannot write %s", path);
    for (int i = 0; i < n_snps; i++) {
        int32_t idx = i;
        int32_t pos = (int32_t)positions[i];
        fwrite(&idx, sizeof(int32_t), 1, fp);
        fwrite(&pos, sizeof(int32_t), 1, fp);
        for (int g = 0; g < MAX_GROUPS; g++) {
            float maf = (g < n_groups) ? results[g].maf[i] : 0.0f;
            fwrite(&maf, sizeof(float), 1, fp);
        }
        for (int g = 0; g < MAX_GROUPS; g++) {
            float var = (g < n_groups) ? results[g].var[i] : 0.0f;
            fwrite(&var, sizeof(float), 1, fp);
        }
        for (int g = 0; g < MAX_GROUPS; g++) {
            int32_t nc = (g < n_groups) ? results[g].n_complete[i] : 0;
            fwrite(&nc, sizeof(int32_t), 1, fp);
        }
    }
    fclose(fp);
}

static void write_summary(const char *out_dir, const Control *c,
                          const WindowsJSON *wj,
                          const GroupResult *results, int n_groups,
                          int n_snps_unique, int n_snps_used,
                          long n_pairs, const long *positions, double t_total,
                          int thinning_applied)
{
    char path[MAX_PATH];
    snprintf(path, MAX_PATH, "%s/summary.tsv", out_dir);
    FILE *fp = fopen(path, "w");
    if (!fp) diefmt("cannot write %s", path);
    fprintf(fp, "engine_version\tfast_ld 0.2.0\n");
    fprintf(fp, "chrom\t%s\n", wj->chrom);
    fprintf(fp, "window_lo\t%d\n", c->win_lo);
    fprintf(fp, "window_hi\t%d\n", c->win_hi);
    fprintf(fp, "shelf_start\t%ld\n", c->shelf_start);
    fprintf(fp, "shelf_end\t%ld\n", c->shelf_end);
    fprintf(fp, "n_snps_unique_in_range\t%d\n", n_snps_unique);
    fprintf(fp, "n_snps_used\t%d\n", n_snps_used);
    fprintf(fp, "thinning_applied\t%d\n", thinning_applied);
    fprintf(fp, "snp_cap\t%d\n", c->snp_cap);
    fprintf(fp, "thin_to\t%d\n", c->thin_to);
    fprintf(fp, "n_pairs\t%ld\n", n_pairs);
    fprintf(fp, "n_groups\t%d\n", n_groups);
    fprintf(fp, "max_groups\t%d\n", MAX_GROUPS);
    fprintf(fp, "sites_record_size_bytes\t%d\n", 8 + 12 * MAX_GROUPS);
    fprintf(fp, "compute_seconds_total\t%.4f\n", t_total);
    if (n_snps_used > 0) {
        fprintf(fp, "first_pos\t%ld\n", positions[0]);
        fprintf(fp, "last_pos\t%ld\n", positions[n_snps_used - 1]);
    }
    for (int g = 0; g < n_groups; g++) {
        const GroupResult *r = &results[g];
        const char *gn = c->group_names[g];
        fprintf(fp, "%s.n_samples\t%d\n", gn, r->n_group_samples);
        fprintf(fp, "%s.compute_seconds\t%.4f\n", gn, r->compute_seconds);
        fprintf(fp, "%s.median_r2_overall\t%.6f\n", gn, r->median_r2_overall);
        fprintf(fp, "%s.pct_pairs_above_0_8\t%.6f\n", gn, r->pct_above_0_8);
        fprintf(fp, "%s.median_r2_shelf\t%.6f\n", gn, r->median_r2_shelf);
        fprintf(fp, "%s.median_r2_flank\t%.6f\n", gn, r->median_r2_flank);
        fprintf(fp, "%s.shelf_ratio\t%.6f\n", gn, r->shelf_ratio);
        for (int b = 0; b < 10; b++)
            fprintf(fp, "%s.decay_decile_%d\t%.6f\n", gn, b, r->decay_deciles[b]);
    }
    fclose(fp);
}

/* =========================================================================
 * main
 * =======================================================================*/

int main(int argc, char **argv) {
    if (argc != 2) {
        fprintf(stderr, "usage: %s <control.tsv>\n", argv[0]);
        return 2;
    }
    Control c;
    read_control(argv[1], &c);

#ifdef _OPENMP
    if (c.threads > 0) omp_set_num_threads(c.threads);
#endif
    mkdir(c.out_dir, 0755);

    double t0 = now_s();

    /* Step 1: parse windows JSON */
    WindowsJSON wj;
    parse_windows_json(c.windows_json, &wj);
    fprintf(stderr, "fast_ld: windows JSON: chrom=%s, n_windows=%d\n",
            wj.chrom, wj.n_windows);
    if (c.win_hi >= wj.n_windows) {
        fprintf(stderr, "fast_ld: warn: window_hi=%d clamped to %d\n",
                c.win_hi, wj.n_windows - 1);
        c.win_hi = wj.n_windows - 1;
    }

    /* Step 2: collect unique SNP positions in window range */
    int n_unique;
    long *unique_pos = unique_positions_in_range(&wj, c.win_lo, c.win_hi, &n_unique);
    fprintf(stderr, "fast_ld: window range [%d, %d] → %d unique SNPs\n",
            c.win_lo, c.win_hi, n_unique);

    /* Step 3: cap or thin */
    long *use_pos;
    int n_use;
    int thinning_applied = 0;
    if (c.thin_to > 0) {
        if (c.thin_to >= n_unique) {
            use_pos = unique_pos;
            n_use = n_unique;
        } else {
            use_pos = thin_positions(unique_pos, n_unique, c.thin_to, &n_use);
            thinning_applied = 1;
            fprintf(stderr, "fast_ld: thin_to=%d applied: %d → %d SNPs\n",
                    c.thin_to, n_unique, n_use);
        }
    } else {
        if (n_unique > c.snp_cap) {
            diefmt(
                "request would compute on %d SNPs, exceeds snp_cap=%d.\n"
                "Either: (a) request a smaller window range, "
                "(b) raise snp_cap, or "
                "(c) set thin_to=N to even-bp-thin to N SNPs.",
                n_unique, c.snp_cap);
        }
        use_pos = unique_pos;
        n_use = n_unique;
    }

    /* Step 4: read dosage selectively */
    DosageMatrix mat = {0};
    read_dosage_selective(&c, wj.chrom, use_pos, n_use, &mat);
    if (mat.n_snps < 2) diefmt("fewer than 2 SNPs read from dosage");

    long n_pairs = (long)mat.n_snps * (mat.n_snps - 1) / 2;

    /* Step 5: resolve groups */
    GroupResolved groups[MAX_GROUPS] = {0};
    for (int g = 0; g < c.n_groups; g++) {
        resolve_group(&mat, c.group_names[g], c.group_ids_csv[g], &groups[g]);
        fprintf(stderr, "fast_ld: group '%s' → %d samples\n",
                groups[g].name, groups[g].n);
    }

    /* Step 6: compute per-group r² */
    GroupResult results[MAX_GROUPS] = {0};
    for (int g = 0; g < c.n_groups; g++) {
        double tg = now_s();
        int k_n = groups[g].n;
        float *Z = xmalloc(sizeof(float) * (size_t)mat.n_snps * k_n);
        build_submatrix(&mat, &groups[g], Z);

        results[g].n_snps = mat.n_snps;
        results[g].n_group_samples = k_n;
        results[g].n_pairs = n_pairs;
        results[g].positions = mat.positions;
        results[g].maf = xmalloc(sizeof(float) * mat.n_snps);
        results[g].var = xmalloc(sizeof(float) * mat.n_snps);
        results[g].n_complete = xmalloc(sizeof(int) * mat.n_snps);
        results[g].r2_q = xmalloc(sizeof(uint8_t) * n_pairs);

        standardize_rows(Z, mat.n_snps, k_n,
                         results[g].maf, results[g].var, results[g].n_complete);
        compute_r2(Z, mat.n_snps, k_n, results[g].var, results[g].r2_q);
        compute_summary(&results[g], &c);
        results[g].compute_seconds = now_s() - tg;
        fprintf(stderr, "fast_ld: group '%s' done in %.3f s "
                        "(median r²=%.4f, shelf_ratio=%.3f)\n",
                c.group_names[g], results[g].compute_seconds,
                results[g].median_r2_overall, results[g].shelf_ratio);
        free(Z);
    }

    /* Step 7: write outputs */
    for (int g = 0; g < c.n_groups; g++)
        write_pairs_bin(c.out_dir, c.group_names[g], results[g].r2_q, n_pairs);
    write_sites_bin(c.out_dir, mat.n_snps, mat.positions, results, c.n_groups);
    double t_total = now_s() - t0;
    write_summary(c.out_dir, &c, &wj, results, c.n_groups,
                  n_unique, mat.n_snps, n_pairs, mat.positions,
                  t_total, thinning_applied);

    fprintf(stderr, "fast_ld: total %.3f s, outputs in %s\n", t_total, c.out_dir);

    /* Cleanup */
    for (int g = 0; g < c.n_groups; g++) {
        free(results[g].maf);
        free(results[g].var);
        free(results[g].n_complete);
        free(results[g].r2_q);
        free(groups[g].col_idx);
        free(c.group_ids_csv[g]);
    }
    if (use_pos != unique_pos) free(use_pos);
    free(unique_pos);
    free_matrix(&mat);
    free_windows_json(&wj);
    return 0;
}
