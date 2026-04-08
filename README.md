This repository contains example MATLAB code accompanying the publication:

“Distance-weighted reflectance for arbitrary source-detector geometries from a single pencil-beam Monte Carlo simulation”

The code is provided as an illustrative implementation of the computational framework described in the paper. Its purpose is to help readers understand the workflow, inspect the geometry-dependent distance distributions, and reproduce the main single-integral approximation (SIA) versus full source-detector (FSD) Monte Carlo comparisons presented in the manuscript.

## Tested environment

This example code was tested with:

- MCX v2025 (“Jumbo Jolt”, Rev 5332f0, build date 2025-02-16)
- MCXLAB
- MATLAB R2025a
- Ubuntu 24.04.2 LTS

## Important note

This repository is not intended to be actively maintained for compatibility with future software updates, package changes, or new MATLAB/MCX releases. It is provided as an example implementation corresponding to the software environment listed above.

Users may therefore need to make minor adjustments when running the code in other environments or with newer software versions.

## Repository structure

The repository is organized around reusable geometry utilities and a single master demonstration script.

- `utils/`
  Contains reusable helper code.

  - `utils/distance_pdfs/`
    Analytical and semi-analytical functions for evaluating the source-detector distance distribution `p(rho)`, or the effective distance distribution `p_eff(rho)` for non-uniform illumination and detection.

  - `utils/sampled_distance_pdfs/`
    Monte Carlo sampling utilities that generate random source and detector point pairs and estimate the corresponding sampled distance distributions. These functions are used for geometry visualization and for comparison against the analytical distance distributions.

- `DEMO_SIA_FSD_all_geometries.m`
  The main demonstration script. This script:
  - compares analytical distance distributions against Monte Carlo-sampled distance distributions for the included geometries,
  - visualizes the source and detector geometries,
  - and compares the single-integral approximation (SIA) against explicit full source-detector (FSD) Monte Carlo simulations.

For the actual SIA computation, the master script uses the analytical distance-distribution functions from `utils/analytical_distance_pdfs/`. The Monte Carlo-sampled distance functions are used only for overview plots and geometry/PDF sanity checks.

## What the code demonstrates

The repository includes examples for overlapping, concentric, non-overlapping, and non-uniformly weighted source-detector geometries. The master script demonstrates how a single localized pencil-beam Monte Carlo simulation can be reused across multiple source-detector configurations by combining the simulated radial reflectance information with the appropriate geometry-dependent distance distribution.

In addition to reflectance comparisons, the code also demonstrates the corresponding absorption-weighted detected pathlength distributions.
