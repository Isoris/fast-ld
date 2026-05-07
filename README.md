# fast-ld

Small fast regional LD engine.

This repository contains a compact LD engine and wrappers for fast regional linkage-disequilibrium calculations.

## Scope

`fast-ld` owns only LD-like computation, including:

- regional LD summaries
- shelf-LD diagnostics
- fast LD calculations for inversion candidate review
- small wrappers/tests around the LD engine

## What this repo is not

This repository does not own:

- general population-statistics engines such as FST, dXY, theta, or local Q
- inversion discovery logic
- karyotype interpretation
- marker tiering
- Atlas UI pages

## Project boundary

The intended architecture is:

```text
fast-ld
  = small reusable LD engine

unified-ancestry
  = broader population-genetic compute backend

catfish-inversion-analysis
  = project-specific catfish inversion analysis

inversion-atlas
  = browser review and visualization interface
