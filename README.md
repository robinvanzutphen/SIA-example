# SIA–FSD Monte Carlo demonstration

This repository contains example MATLAB code accompanying the publication:

> "Distance-weighted reflectance for arbitrary source–detector geometries
> from a single pencil-beam Monte Carlo simulation"

The code is provided as an illustrative implementation of the computational
framework described in the paper. Its purpose is to help readers understand
the workflow, inspect the geometry-dependent distance distributions, and
reproduce the main single-integral approximation (SIA) versus full
source–detector (FSD) Monte Carlo comparisons presented in the manuscript.

---

## Tested environment

- MCX v2025 ("Jumbo Jolt", Rev 5332f0, build date 2025-02-16)
- MCXLAB
- MATLAB R2025a
- Ubuntu 24.04.2 LTS

---

## Important note

This repository is not intended to be actively maintained for compatibility
with future software updates, package changes, or new MATLAB/MCX releases.
It is provided as an example implementation corresponding to the software
environment listed above. Users may therefore need to make minor adjustments
when running the code in other environments or with newer software versions.

---

## Runtime and hardware requirements

Running the full demonstration script with default settings launches a large
number of MCX simulations and can take **15–40 minutes** on modern hardware.
A warning is printed at the top of the script with an estimated runtime.

MCX uses **NVIDIA CUDA** for GPU-accelerated photon transport. An NVIDIA GPU
with a compatible driver is therefore required. AMD and Intel GPUs are not
supported by the standard MCX build.

---

## Path setup

Before running the demonstration script, ensure the following directories are
on your MATLAB path:

```matlab
addpath('/path/to/mcxlab');
addpath('/path/to/utils/core');
addpath('/path/to/utils/helpers');
addpath('/path/to/utils/distance_pdfs');
addpath('/path/to/utils/sampled_distance_pdfs');
```

---

## Repository structure

```
.
├── DEMO_SIA_FSD_ALL_GEOMETRIES.m
└── utils/
    ├── core/
    ├── helpers/
    ├── distance_pdfs/
    └── sampled_distance_pdfs/
```

### `DEMO_SIA_FSD_ALL_GEOMETRIES.m`

The main demonstration script. It:

- visualizes the source and detector geometries for each scenario,
- compares analytical distance distributions against Monte Carlo-sampled distributions,
- and compares SIA against explicit FSD Monte Carlo simulations for both reflectance and absorption-weighted detected pathlength distributions.

For the actual SIA computation the script uses the analytical distance
distributions from `utils/distance_pdfs/`. The Monte Carlo-sampled distance
functions are used only for overview plots and sanity checks.

### `utils/core/`

Contains the two central computational functions that implement the SIA framework described in the paper:

- `sia_reflectance.m` — evaluates the distance-weighted reflectance integral (Eq. 3 of the Letter) from pencil-beam photon data and a precomputed distance distribution.
- `sia_pathlength.m` — evaluates the corresponding absorption-weighted pathlength distribution (Eq. 5 of the Letter).

Both functions are heavily commented and include explicit references to the relevant equations in the manuscript.

### `utils/helpers/`

Contains utility functions that support the main script but are not part of the core SIA computation:

- `build_sia_cases.m` — constructs the cell array of source–detector geometry structs for all demonstration scenarios.
- `build_mcx_cfg_fsd.m` — builds an MCX configuration struct for a given explicit full source–detector simulation, branching off a shared base configuration.
- `plot_sia_fsd_geometries.m` — geometry scatter plots.
- `plot_sia_fsd_prho_overview.m` — sampled vs. analytical p(rho) overview.
- `plot_sia_fsd_reflectance_hists.m` — reflectance distribution histograms.
- `plot_sia_fsd_pathlength_pdfs.m` — detected pathlength PDF curves.

### `utils/distance_pdfs/`

Analytical and semi-analytical functions for evaluating the source–detector
distance distribution `p(rho)`, or the effective distance distribution
`p_eff(rho)` for non-uniform illumination and detection. These are the
distributions that enter the SIA integral directly.

### `utils/sampled_distance_pdfs/`

Monte Carlo sampling utilities that generate random source and detector point
pairs and estimate the corresponding distance distributions. Used only for
geometry visualization and comparison against the analytical distributions.

---

## What the code demonstrates

The repository covers overlapping, concentric, non-overlapping, and
non-uniformly weighted source–detector geometries. The master script
demonstrates how a **single** localized pencil-beam Monte Carlo simulation
can be reused across multiple source–detector configurations by combining the
simulated radial reflectance profile with the appropriate geometry-dependent
distance distribution — the central result of the accompanying publication.
