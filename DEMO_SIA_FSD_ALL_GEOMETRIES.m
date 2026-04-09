clc; clear all; close all;

%% ========================================================================
%  DEMO_SIA_FSD_ALL_GEOMETRIES
%
%  Demonstration script accompanying the Letter:
%    "Distance-weighted reflectance for arbitrary source-detector
%     geometries from a single pencil-beam Monte Carlo simulation"
%
%  PURPOSE
%  -------
%  This script demonstrates how the single-integral
%  approximation (SIA) was evaluated and compared against explicit full
%  source-detector (FSD) Monte Carlo simulations across a representative
%  set of source-detector geometries.
%
%  The central idea is that, for laterally homogeneous media, the lateral
%  source-detector geometry enters only through a distance distribution
%  p_SD(rho) (or p_eff(rho) for non-uniform weighting). Photon transport
%  itself is obtained from ONE localized pencil-beam Monte Carlo run and
%  is then reused across many geometries by post-processing the exiting
%  photon coordinates and absorption-weighted pathlengths.
%
%  WHAT THIS SCRIPT SHOWS
%  ----------------------
%  1) Geometry overview plots:
%       - randomly sampled source and detector points for each case
%       - MC-sampled p(rho) or p_eff(rho)
%       - analytical / semi-analytical overlay used for the actual SIA
%
%  2) Equivalence between SIA and FSD:
%       - reflectance distributions over repeated Monte Carlo runs
%       - detected pathlength distributions
%
%  IMPORTANT IMPLEMENTATION NOTES
%  ------------------------------
%  - For the ACTUAL SIA computation, only the analytical / semi-analytical
%    distance distributions are used.
%  - MC-sampled distance distributions are used only for the geometry/PDF
%    sanity plots.
%  - Within each repetition, a SINGLE pencil-beam MC simulation is run and
%    then reused across all selected geometries. This is the key efficiency
%    gain of the SIA.
%  - Photon pathlength contributions are absorption-weighted by exp(-mu_a L)
%    in both the FSD and SIA paths, consistent with the Letter.
%  - The displayed pathlength curves are formed by pooling the weighted
%    pathlength contributions from ALL repetitions and normalizing only once
%    at the end. The reflectance distributions, in contrast, remain defined
%    at the repetition level.
%  - The rho-integral is discretized here as a simple Riemann sum over a
%    uniform rho-grid. Other quadrature choices (e.g. trapz or adaptive
%    quadrature) would also be possible; convergence-optimization is not
%    investigated in this demonstration.
%
%  RELATION TO THE LETTER
%  ----------------------
%  The Letter validates six scenarios (I-VI). This script additionally
%  includes a concentric core->total case that is not shown in the Letter.
%  It is included here because it is a simple test case, but it was omitted
%  from the manuscript figures because it does not materially broaden the
%  validation beyond the other concentric examples.
%
%  !! RUNTIME WARNING !!
%  ---------------------
%  With default settings (n_repetitions = 10, n_photons = 1e7, 7 cases)
%  this script already launches 80 MCX simulations Reduce n_repetitions or n_photons for a quick test.
%  See Section 1 (USER CONTROLS) below.
%  Moreover, it also executes a stochastic sampling method to approximate
%  and validate the analytical PDFs, which can also take up a minute. 
%
%  GPU REQUIREMENT
%  ---------------
%  MCX (mcxlab) uses NVIDIA CUDA for GPU-accelerated photon transport.
%  An NVIDIA GPU with a recent CUDA-capable driver is therefore required
%  to run this script. The GPU index is set via the gpu_id variable in
%  Section 1. There also exists a CPU-only (MCXLAB-CL) version. 
%
%  PATH SETUP
%  ----------
%  Before running, make sure the following directories are on your MATLAB
%
%    addpath('/path/to/mcxlab');          % MCXlab mex interface
%    addpath('/path/to/sia_scripts');     % analytical p(rho) utilities
%                                         % and plot helper functions
%
%  The specific companion .m files required are listed under
%  "COMPANION FILES REQUIRED" further below.
%
%  COMPANION FILES REQUIRED ON THE MATLAB PATH
%  -------------------------------------------
%  Core SIA functions:
%     sia_reflectance.m
%     sia_pathlength.m
%
%  Geometry builder:
%     build_sia_cases.m
%
%  MCX config builder:
%     build_mcx_cfg_fsd.m
%
%  Analytical distance-distribution utilities:
%     p_rho_uniform_overlap_analytic.m
%     p_rho_core_total_analytic.m
%     p_rho_annulus_annulus_analytic.m
%     p_rho_disk_annulus_analytic.m
%     p_rho_nonoverlap_disks_analytic.m
%     p_rho_square_square_analytic.m
%     p_eff_gaussian_overlap_analytic.m
%
%  MC sampling utilities (overview plots only):
%     sample_p_rho_uniform_overlap_mc.m
%     sample_p_rho_core_total_mc.m
%     sample_p_rho_annulus_annulus_mc.m
%     sample_p_rho_disk_annulus_mc.m
%     sample_p_rho_nonoverlap_disks_mc.m
%     sample_p_rho_square_square_mc.m
%     sample_p_eff_gaussian_overlap_mc.m
%
%  Plot helpers:
%     plot_sia_fsd_geometries.m
%     plot_sia_fsd_prho_overview.m
%     plot_sia_fsd_reflectance_hists.m
%     plot_sia_fsd_pathlength_pdfs.m
%
%  MCX / MCXlab:
%     mcxlab must be installed and available on the MATLAB path.
%     See: https://mcx.space/
%
%  DEFAULTS VS LETTER SETTINGS
%  ---------------------------
%  The defaults below are intentionally modest so that the script remains a
%  practical demo. The Letter itself used substantially larger Monte Carlo
%  budgets (n_repetitions = 5000, n_photons = 1e7). Increase the number of
%  photons and repetitions if manuscript-quality statistics are desired.
% ========================================================================

%% ========================================================================
%  0) USER CONTROLS
% ========================================================================


% Plot switches: set to false to skip the corresponding figure group.
plot_geometry_overview  = true;   % geometry scatter plots and p(rho) overview
plot_equivalence_results = true;  % reflectance histograms and pathlength PDFs

% Case selection (1 = run, 0 = skip).
% The order matches the AllCases array built in build_sia_cases.m:
%   [I, Extra, II, III, IV, V, VI]
%    I     = uniform overlapping disk -> disk
%    Extra = concentric core -> total (not in Letter)
%    II    = annulus -> annulus (concentric)
%    III   = core -> cladding (concentric disk -> annulus)
%    IV    = non-overlapping disks
%    V     = square -> square
%    VI    = Gaussian source + Gaussian-weighted detector
run_case = [1 1 1 1 1 1 1];

%% ========================================================================
%  1) OPTICAL / MONTE CARLO SETTINGS
% ========================================================================

% ----------------------- Optical properties ------------------------------
% These are the medium and launch parameters used in the Letter (Sec. 3).
NA    = 0.22;          % launch and collection numerical aperture [dimensionless]
n_med = 1.35;          % refractive index of the scattering medium [dimensionless]
theta_acc_max = asin(NA / n_med);  % half-angle acceptance cone [rad]

mu_a  = 0.10;          % absorption coefficient [1/mm]
mu_sp = 5.00;          % reduced scattering coefficient mu_s*(1-g) [1/mm]
g     = 0.90;          % Henyey-Greenstein anisotropy factor [dimensionless]
mu_s  = mu_sp / (1 - g); % scattering coefficient [1/mm]

% -------------------- Monte Carlo repetition settings --------------------
% These defaults are intentionally reduced relative to the Letter so that
% the script remains runnable as a demonstration. Increase them for more
% stable distributions at the cost of runtime.
n_repetitions = 2;    % number of independent MC repeat runs per scenario
n_photons     = 5e6;   % number of photons per individual MCX call
gpu_id        = '11';  % MCX GPU selection string (see mcxlab documentation)

% ------------------------ Computational domain ---------------------------
% Here we represent a 100 x 100 x 100 mm homogeneous box using a SINGLE
% voxel. This avoids the overhead of a large 3D fluence grid. Because the
% domain collapses to one voxel, the coordinate axes are defined in voxel
% units where ONE voxel = mm_per_voxel millimetres.
%
% Consequences:
%   - Volume array is 1 x 1 x 1 (one voxel).
%   - All MCX source/detector sizes in srcparam1/2 must be given in voxels
%     (physical mm / mm_per_voxel).
%   - Photon exit coordinates and pathlengths returned by MCX are in voxel
%     units and must be converted back to mm by multiplying by mm_per_voxel.

mm_per_voxel = 100;                      % physical size of one voxel [mm/voxel]
Lxy = 1;                                 % domain extent in x and y [voxels]
Lz  = 1;                                 % domain extent in z [voxels]
volume = uint8(ones(Lxy, Lxy, Lz));      % medium label array (1 = tissue, 0 = background)

% Source is placed at the center of the top face (z = 0) of the single voxel.
% With issrcfrom0 = 1 in MCX, position [0.5, 0.5, 0] maps to the midpoint
% of the top surface in voxel coordinates.
source_position_vox = [0.5 0.5 0];      % source position [voxels]
source_direction    = [0 0 1];          % downward (+z) launch direction [unit vector]

% Physical source center in mm, derived from the voxel position above.
% Used later to compute photon exit distances relative to the source.
source_center_mm    = double(source_position_vox(1:2)) * mm_per_voxel;  % [mm]

% ------------------------- Pathlength histogram --------------------------
% Detected pathlength distributions are accumulated as histograms over
% repetitions. L_max_mm and n_L_bins set the histogram support and resolution.
L_max_mm = 5;                                              % max pathlength shown [mm]
n_L_bins = 1000;                                           % number of histogram bins
L_edges  = linspace(0, L_max_mm, n_L_bins + 1);           % bin edges [mm]
L_mid    = (L_edges(1:end-1) + L_edges(2:end)) / 2;       % bin centres [mm]
dL       = L_edges(2) - L_edges(1);                       % bin width [mm]

% ----------------------- rho discretization for SIA ----------------------
% The SIA reflectance integral (Eq. 3 of the Letter) is evaluated as a
% Riemann sum over a uniform rho-grid. n_rho_bins_mix controls the
% resolution of this grid; 5000 bins provides a good balance between
% accuracy and memory. 
n_rho_bins_mix = 5000;  % number of bins in the rho grid for the SIA integral

% ------------------- MC sampling only for overview plots -----------------
% These parameters control the density of the geometry scatter plots and
% the MC-sampled p(rho) curves shown in the overview figures.
% They do NOT affect the SIA or FSD reflectance/pathlength estimates.
n_rho_bins_plot  = 50;    % histogram bins for the MC-sampled p(rho) overview
n_geometry_plot  = 5000;  % number of source/detector point pairs to scatter-plot
n_pairs_plot_pdf = 2e6;   % number of random pairs drawn to estimate p(rho) for plot

% ---------------------- Gaussian p_eff settings --------------------------
% For Case VI (non-uniform Gaussian illumination/detection), the effective
% distance distribution p_eff(rho) is computed via an analytical integral
% evaluated numerically using a Hankel (Fourier-Bessel) transform.
% The parameters below control the numerical accuracy of that transform.
% See p_eff_gaussian_overlap_analytic.m and the Supplemental Material for
% details on the integration scheme.
NrS_gauss   = 1100;   % radial quadrature points over the Gaussian source cap
NrD_gauss   = 1100;   % radial quadrature points over the detector aperture
NrInt_gauss = 4000;   % integration points for the intermediate real-space integral
Nk_gauss    = 2500;   % number of spatial-frequency (k) points in the Hankel transform
kmax_gauss  = 350;    % maximum spatial frequency [1/mm] in the Hankel transform

%% ========================================================================
%  2) BASE MCX CONFIGURATION
% ========================================================================
% This struct holds all MCX settings that are SHARED between the pencil-beam
% run and every explicit FSD run. Individual runs then branch off this base
% config and set only the source-type-specific fields.
% See mcxlab.m (distributed with MCXlab) for a full description of each field.

mcx_cfg_base = struct();

mcx_cfg_base.nphoton      = n_photons;     % number of photons to launch per call
mcx_cfg_base.vol          = volume;        % medium label volume (1 x 1 x 1 here)
mcx_cfg_base.unitinmm     = mm_per_voxel;  % physical size of one voxel [mm]

mcx_cfg_base.srcpos       = source_position_vox;  % source position [voxels]
mcx_cfg_base.srcdir       = source_direction;      % launch direction [unit vector]

mcx_cfg_base.autopilot    = 1;     % let MCX choose thread/block layout automatically
mcx_cfg_base.tstart       = 0;     % simulation start time [s]
mcx_cfg_base.tend         = 5e-9;  % simulation end time [s]; set large enough to
                                    % capture all photons reaching the surface
mcx_cfg_base.tstep        = 5e-9;  % time-gate width [s]; equals tend for CW mode

mcx_cfg_base.isreflect    = 0;     % 0 = ignore Fresnel reflections at boundaries
mcx_cfg_base.respin       = 1;     % number of sub-runs to accumulate per mcxlab call
mcx_cfg_base.issave2pt    = 0;     % 0 = do not save 3D fluence (saves memory/time)
mcx_cfg_base.savedetflag  = 'dspvx'; % save: photon exit data (d), scattering path (s),
                                      %       partial pathlength (p), exit direction (v),
                                      %       and exit position (x)
mcx_cfg_base.maxdetphoton = n_photons; % pre-allocated buffer for detected photons

mcx_cfg_base.gpuid        = gpu_id;  % GPU index string; change if multiple GPUs present
mcx_cfg_base.issrcfrom0   = 1;       % interpret srcpos relative to the corner of voxel 0

% Random seed: -1 tells MCX to use a time-based seed, giving independent
% draws across repetitions without manual seed management.
mcx_cfg_base.seed         = -1;

% Finite-NA launch: The angular distribution of launched photons is encoded
% as an inverse CDF sampled uniformly from [0, theta_acc_max/pi]. This
% enforces a uniform angular distribution within the acceptance cone of the
% specified NA. Under the SIA assumptions (Letter Sec. 2B), the launch NA
% is embedded in R(rho) and does not change the form of the SIA integral.
mcx_cfg_base.angleinvcdf  = linspace(0, theta_acc_max/pi, 5);

% Detection boundary: 'aa_aaa001000' means the -z face (top surface, where
% photons exit back into air) is the detection plane. All other faces are
% absorbing. See mcxlab documentation for the boundary condition string format.
mcx_cfg_base.bc           = 'aa_aaa001000';

% Optical properties table. Each row is [mu_a, mu_s, g, n] for a medium.
% Row 1 is the background (void); Row 2 is the tissue medium.
mcx_cfg_base.prop         = [0 0 1 1; mu_a mu_s g n_med];

% Total number of photons launched per repetition (used for normalisation).
n_launched_photons = double(mcx_cfg_base.nphoton * mcx_cfg_base.respin);

%% ========================================================================
%  3) GEOMETRY DEFINITIONS
% ========================================================================
% The AllCases cell array defines the source-detector configuration for
% each scenario. It is built by the companion function
% build_sia_cases.m. See that file for documentation on each field.
%
% Manuscript mapping:
%   I   = uniform overlapping disk -> disk
%   II  = annulus -> annulus (concentric)
%   III = core -> cladding (concentric disk -> annulus)
%   IV  = non-overlapping disks
%   V   = square -> square
%   VI  = Gaussian source + Gaussian-weighted detector
%
% Additional demo-only case:
%   Extra = concentric core -> total (not shown in the Letter figures)

% --------------------- Circular geometries (mm) --------------------------
R_total = 0.50;   % outer (total fiber / annulus) radius [mm]
r_core  = 0.25;   % inner (core) radius [mm]

% ------------------- Non-overlapping disk geometry -----------------------
r_source_nonoverlap   = 0.30;   % source disk radius [mm]
r_detector_nonoverlap = 0.45;   % detector disk radius [mm]
% Centre-to-centre separation: disks just touch plus a small gap (0.25*R_total).
sep_nonoverlap = r_source_nonoverlap + r_detector_nonoverlap + 0.25*R_total;  % [mm]

% -------------------------- Square geometry ------------------------------
side_source_square    = 0.50;   % source square side length [mm]
side_detector_square  = 0.50;   % detector square side length [mm]
% Offset places the two squares corner-to-corner (touching diagonally).
detector_offset_square_mm = [side_source_square, side_source_square];  % [mm, mm]

% ------------------- Non-uniform Gaussian geometry -----------------------
gaussian_w_mm     = 0.50;  % 1/e^2 beam waist radius in exp(-2r^2/w^2) [mm]
R_detector_gauss  = 0.50;  % hard aperture radius of the detector [mm]
Nw_cap            = 10;    % source truncation cap: R_src = R_det + Nw_cap*w [dimensionless]
                            % (chosen large enough that the Gaussian tail is negligible)

% Build the cell array of case structs (see build_sia_cases.m for details).
AllCases = build_sia_cases(R_total, r_core, ...
                            r_source_nonoverlap, r_detector_nonoverlap, sep_nonoverlap, ...
                            side_source_square, side_detector_square, detector_offset_square_mm, ...
                            gaussian_w_mm, R_detector_gauss, Nw_cap);

% Select only the cases flagged in run_case.
selected_idx = find(logical(run_case(:).'));
Cases = AllCases(selected_idx);
n_cases = numel(Cases);

if n_cases == 0
    error('No cases selected. Set at least one entry of run_case to 1.');
end

%% ========================================================================
%  4) PRECOMPUTE GEOMETRY OVERVIEW DATA AND ANALYTICAL p(rho)
% ========================================================================
% For each case we prepare two different types of information:
%
%  A) Analytical / semi-analytical distance distribution used in the ACTUAL
%     SIA reflectance and pathlength scaling.  This is stored in
%     p_rho_mix_store and is the quantity that enters Eq. (3) of the Letter.
%
%  B) MC-sampled point clouds and MC-sampled distance histograms used ONLY
%     for the overview figures. They serve as a visual sanity check that
%     the analytical p(rho) matches random point-pair sampling.

% Preallocate storage cells for the precomputed quantities.
overview_source_points    = cell(n_cases,1);  % scatter points for source aperture
overview_detector_points  = cell(n_cases,1);  % scatter points for detector aperture

rho_mid_mix_store         = cell(n_cases,1);  % rho bin centres used in SIA integral
rho_edges_mix_store       = cell(n_cases,1);  % rho bin edges
p_rho_mix_store           = cell(n_cases,1);  % normalised analytical p(rho) for SIA
p_rho_theory_store        = cell(n_cases,1);  % analytical p(rho) for overlay in plot
p_rho_theory_label_store  = cell(n_cases,1);  % legend label for analytical curve

rho_plot_store            = cell(n_cases,1);  % rho axis for MC-sampled overview plot
p_rho_plot_sampled_store  = cell(n_cases,1);  % MC-sampled p(rho) values for overview

detector_area_store       = nan(n_cases,1);   % effective detector area A_d [mm^2]

for c = 1:n_cases
    CaseNow = Cases{c};

    switch CaseNow.type

        % -----------------------------------------------------------------
        % Uniform-illumination aperture cases (disk and annulus shapes)
        % -----------------------------------------------------------------
        case 'uniform_shapes'
            % Determine the maximum source-detector separation that can
            % occur for this geometry. This sets the support of p(rho).
            rho_max_mm = support_rhomax_uniform(CaseNow);

            % Build uniform rho grid for the SIA Riemann sum (Eq. 3).
            rho_edges_mix = linspace(0, rho_max_mm, n_rho_bins_mix + 1);
            rho_mid_mix   = (rho_edges_mix(1:end-1) + rho_edges_mix(2:end)) / 2;
            drho_mix      = rho_edges_mix(2) - rho_edges_mix(1);

            % Choose the correct analytical p(rho) function based on the
            % source and detector shapes and their relative positions.
            if strcmp(CaseNow.source_shape,'disk') && strcmp(CaseNow.detector_shape,'disk') && ...
                    all(CaseNow.detector_offset_mm == 0) && (CaseNow.source_rout == CaseNow.detector_rout)
                % Case I: equal, concentric, fully overlapping disks.
                % p(rho) has a closed-form expression (see Supplemental).
                [rho_plot, pdf_samp_plot, source_pts_plot, detector_pts_plot] = ...
                    sample_p_rho_uniform_overlap_mc(2*CaseNow.detector_rout, n_pairs_plot_pdf, n_geometry_plot, n_rho_bins_plot);
                p_mix    = p_rho_uniform_overlap_analytic(rho_mid_mix, 2*CaseNow.detector_rout);
                p_theory = p_mix;
                p_label  = 'analytic';

            elseif strcmp(CaseNow.source_shape,'disk') && strcmp(CaseNow.detector_shape,'disk') && ...
                    all(CaseNow.detector_offset_mm == 0)
                % Extra case: concentric disks with different radii (core->total).
                % Also fully analytical.
                [rho_plot, pdf_samp_plot, source_pts_plot, detector_pts_plot] = ...
                    sample_p_rho_core_total_mc(CaseNow.source_rout, CaseNow.detector_rout, ...
                                               n_pairs_plot_pdf, n_geometry_plot, n_rho_bins_plot);
                p_mix    = p_rho_core_total_analytic(rho_mid_mix, CaseNow.source_rout, CaseNow.detector_rout);
                p_theory = p_mix;
                p_label  = 'analytic';

            elseif strcmp(CaseNow.source_shape,'annulus') && strcmp(CaseNow.detector_shape,'annulus') && ...
                    all(CaseNow.detector_offset_mm == 0)
                % Case II: concentric annular source and detector (annulus->annulus).
                [rho_plot, pdf_samp_plot, source_pts_plot, detector_pts_plot] = ...
                    sample_p_rho_annulus_annulus_mc(CaseNow.source_rin, CaseNow.source_rout, ...
                                                    CaseNow.detector_rin, CaseNow.detector_rout, ...
                                                    n_pairs_plot_pdf, n_geometry_plot, n_rho_bins_plot);
                p_mix    = p_rho_annulus_annulus_analytic(rho_mid_mix, ...
                                                          CaseNow.source_rin, CaseNow.source_rout, ...
                                                          CaseNow.detector_rin, CaseNow.detector_rout);
                p_theory = p_mix;
                p_label  = 'analytic';

            elseif strcmp(CaseNow.source_shape,'disk') && strcmp(CaseNow.detector_shape,'annulus') && ...
                    all(CaseNow.detector_offset_mm == 0)
                % Case III: central disk source, annular detector (core->cladding).
                [rho_plot, pdf_samp_plot, source_pts_plot, detector_pts_plot] = ...
                    sample_p_rho_disk_annulus_mc(CaseNow.source_rout, CaseNow.detector_rin, CaseNow.detector_rout, ...
                                                 n_pairs_plot_pdf, n_geometry_plot, n_rho_bins_plot);
                p_mix    = p_rho_disk_annulus_analytic(rho_mid_mix, ...
                                                       CaseNow.source_rout, CaseNow.detector_rin, CaseNow.detector_rout);
                p_theory = p_mix;
                p_label  = 'analytic';

            elseif strcmp(CaseNow.source_shape,'disk') && strcmp(CaseNow.detector_shape,'disk') && ...
                    any(CaseNow.detector_offset_mm ~= 0)
                % Case IV: laterally separated (non-overlapping) disks.
                % Rotational symmetry is broken; p(rho) requires a 1D numerical
                % integral over the angular variable theta (Eq. 6 of Letter).
                sep_mm = norm(CaseNow.detector_offset_mm);
                [rho_plot, pdf_samp_plot, source_pts_plot, detector_pts_plot] = ...
                    sample_p_rho_nonoverlap_disks_mc(CaseNow.source_rout, CaseNow.detector_rout, sep_mm, ...
                                                     n_pairs_plot_pdf, n_geometry_plot, n_rho_bins_plot);
                p_mix    = p_rho_nonoverlap_disks_analytic(rho_mid_mix, sep_mm, CaseNow.source_rout, CaseNow.detector_rout);
                p_theory = p_mix;
                p_label  = 'semi-analytic (1D integral)';

            else
                error('Unsupported uniform case encountered during precomputation.');
            end

            % Sanitize and normalise the analytical p(rho) to a proper PDF
            % (unit integral over drho) before use in the SIA integral.
            p_mix = sanitize_pdf(p_mix);
            Zmix  = sum(p_mix) * drho_mix;
            if Zmix > 0
                p_mix = p_mix / Zmix;
            end

            % Physical detector area used for the SIA normalisation (A_d in Eq. 3).
            detector_area_store(c) = detector_area_uniform(CaseNow.detector_shape, ...
                                                           CaseNow.detector_rin, CaseNow.detector_rout);

        % -----------------------------------------------------------------
        % Case V: square source and square detector
        % -----------------------------------------------------------------
        case 'square_square'
            % Maximum possible separation: corner-to-corner diagonal of the
            % combined extent of both squares plus the centre offset.
            rho_max_mm = hypot(CaseNow.detector_offset_mm(1) + (CaseNow.source_side/2 + CaseNow.detector_side/2), ...
                               CaseNow.detector_offset_mm(2) + (CaseNow.source_side/2 + CaseNow.detector_side/2));

            rho_edges_mix = linspace(0, rho_max_mm, n_rho_bins_mix + 1);
            rho_mid_mix   = (rho_edges_mix(1:end-1) + rho_edges_mix(2:end)) / 2;
            drho_mix      = rho_edges_mix(2) - rho_edges_mix(1);

            [rho_plot, pdf_samp_plot, source_pts_plot, detector_pts_plot] = ...
                sample_p_rho_square_square_mc(CaseNow.source_side, CaseNow.detector_side, ...
                                              CaseNow.detector_offset_mm(1), CaseNow.detector_offset_mm(2), ...
                                              n_pairs_plot_pdf, n_geometry_plot, n_rho_bins_plot);

            % p(rho) for square apertures is obtained via the covariogram
            % gSD(u) of the two square regions .
            p_mix = p_rho_square_square_analytic(rho_mid_mix, ...
                                                 CaseNow.source_side, CaseNow.detector_side, ...
                                                 CaseNow.detector_offset_mm(1), CaseNow.detector_offset_mm(2));
            p_mix = sanitize_pdf(p_mix);
            Zmix  = sum(p_mix) * drho_mix;
            if Zmix > 0
                p_mix = p_mix / Zmix;
            end

            p_theory = p_mix;
            p_label  = 'semi-analytic (covariogram)';
            % Detector area for a square aperture is simply side^2 [mm^2].
            detector_area_store(c) = CaseNow.detector_side^2;

        % -----------------------------------------------------------------
        % Case VI: Gaussian-weighted source and detector
        % -----------------------------------------------------------------
        case 'gaussian_overlap'
            w_mm         = CaseNow.w_mm;        % beam waist radius [mm]
            R_det_mm     = CaseNow.R_det_mm;    % detector aperture radius [mm]
            % Source is truncated at R_src_cap_mm to keep p_eff numerically tractable.
            % With Nw_cap = 10, the truncated tail
            R_src_cap_mm = R_det_mm + CaseNow.Nw_cap * w_mm;  % source cap radius [mm]
            rho_max_mm   = R_src_cap_mm + R_det_mm;            % max possible separation [mm]

            rho_edges_mix = linspace(0, rho_max_mm, n_rho_bins_mix + 1);
            rho_mid_mix   = (rho_edges_mix(1:end-1) + rho_edges_mix(2:end)) / 2;

            [rho_plot, pdf_samp_plot, source_pts_plot, detector_pts_plot] = ...
                sample_p_eff_gaussian_overlap_mc(w_mm, R_src_cap_mm, R_det_mm, ...
                                                 n_pairs_plot_pdf, n_geometry_plot, n_rho_bins_plot);

            % Analytical p_eff(rho) computed via the Hankel transform route.
            % NrS, NrD, NrInt, Nk, kmax are the quadrature settings defined above.
            % See p_eff_gaussian_overlap_analytic.m for the algorithm.
            gaussian_curves = p_eff_gaussian_overlap_analytic(rho_edges_mix, w_mm, R_src_cap_mm, R_det_mm, ...
                                                              NrS_gauss, NrD_gauss, NrInt_gauss, Nk_gauss, kmax_gauss);
            p_mix    = gaussian_curves.p_base;  % the normalised p_eff(rho) for the SIA
            p_theory = p_mix;
            p_label  = sprintf('analytic integral (R_src = R_det + %dw)', CaseNow.Nw_cap);

            % Effective detector area A_det,eff = integral_{A_d} w_d(r) dA.
            % For a Gaussian w_d(r) = exp(-2r^2/w^2) clipped to radius R_det_mm,
            % this integral evaluates to (pi*w^2/2) * (1 - exp(-2*R_det^2/w^2)).
            % This replaces A_d in Eq. (3) to account for the non-uniform
            % detection sensitivity (Letter Sec. 2C).
            detector_area_store(c) = (pi*w_mm^2/2) * (1 - exp(-2*R_det_mm^2 / w_mm^2));

        otherwise
            error('Unknown case type encountered during precomputation.');
    end

    % Store the scatter-plot source and detector point clouds for plotting.
    overview_source_points{c}   = source_pts_plot(1:min(n_geometry_plot, size(source_pts_plot,1)), :);
    overview_detector_points{c} = detector_pts_plot(1:min(n_geometry_plot, size(detector_pts_plot,1)), :);

    % Store the precomputed rho grids and distributions.
    rho_mid_mix_store{c}        = rho_mid_mix;
    rho_edges_mix_store{c}      = rho_edges_mix;
    p_rho_mix_store{c}          = p_mix;             % normalised, used in SIA integral
    p_rho_theory_store{c}       = sanitize_pdf(p_theory);
    p_rho_theory_label_store{c} = p_label;

    % Store the MC-sampled p(rho) data for the overview plot.
    rho_plot_store{c}           = rho_plot;
    p_rho_plot_sampled_store{c} = pdf_samp_plot;
end

% Determine the global rho support: the pencil-beam simulation only needs
% to retain photons up to the largest rho required by any selected geometry.
rho_max_global_mm = 0;
for c = 1:n_cases
    rho_max_global_mm = max(rho_max_global_mm, rho_edges_mix_store{c}(end));
end

%% ========================================================================
%  5) OPTIONAL OVERVIEW PLOTS
% ========================================================================
% These figures show the geometry layouts and p(rho) curves for a visual
% sanity check before the main MC simulations start.
if plot_geometry_overview
    plot_sia_fsd_geometries(Cases, overview_source_points, overview_detector_points);
    plot_sia_fsd_prho_overview(Cases, rho_plot_store, p_rho_plot_sampled_store, ...
                               rho_mid_mix_store, p_rho_theory_store, p_rho_theory_label_store);
end

%% ========================================================================
%  6) PREALLOCATE OUTPUTS
% ========================================================================
% Reflectance is stored per repetition so that FSD and SIA distributions
% can be compared as histograms (one value per repetition, per case).
R_full = nan(n_cases, n_repetitions);  % FSD reflectance [dimensionless]
R_sia  = nan(n_cases, n_repetitions);  % SIA reflectance [dimensionless]

% Per-repetition mean pathlength and Wasserstein-1 distance between the
% FSD and SIA pathlength CDFs, used as summary metrics.
meanL_full = nan(n_cases, n_repetitions);  % mean detected pathlength, FSD [mm]
meanL_sia  = nan(n_cases, n_repetitions);  % mean detected pathlength, SIA [mm]
W1_L       = nan(n_cases, n_repetitions);  % W1 distance between FSD and SIA CDFs [mm]

% For the displayed pathlength curves we pool ALL weighted pathlength mass
% across repetitions and normalize only once at the end. This gives smoother
% curves than averaging per-repetition PDFs.
pL_full_mass_total = zeros(n_cases, n_L_bins);  % accumulated FSD pathlength mass
pL_sia_mass_total  = zeros(n_cases, n_L_bins);  % accumulated SIA pathlength mass

%% ========================================================================
%  7) BUILD THE SHARED PENCIL-BEAM CONFIGURATION
% ========================================================================
% This is the configuration for the reusable (SIA) simulation. Within each
% repetition it is run exactly once, and its detected photons are then
% post-processed for every geometry via the corresponding p(rho).

mcx_cfg_pencil = mcx_cfg_base;
mcx_cfg_pencil.srctype = 'pencil';  % point source; localized launch at srcpos
% Remove any source-shape parameters that might be left from a previous run.
mcx_cfg_pencil = rmfield_safe(mcx_cfg_pencil, 'srcparam1');
mcx_cfg_pencil = rmfield_safe(mcx_cfg_pencil, 'srcparam2');

% Total number of simulation steps: one pencil + one FSD per case per repetition.
n_progress_steps = n_repetitions * (n_cases + 1);
progress_step    = 0;

% Waitbar with ETA. The ETA is estimated from elapsed wall-clock time.
wb = waitbar(0, 'Starting...', 'Name', 'MCX SIA vs FSD progress');
cleanupObj = onCleanup(@() close_waitbar_safe(wb)); %#ok<NASGU>
t_loop_start = tic;  % start wall-clock timer for ETA computation

%% ========================================================================
%  8) MAIN LOOP
% ========================================================================
% Loop structure:
%   outer loop  = Monte Carlo repetitions (rep = 1 ... n_repetitions)
%   inner loop  = source-detector geometry cases (c = 1 ... n_cases)
%
% Per repetition:
%   8A) ONE shared pencil-beam run  --> provides photons for ALL SIA cases
%   8B) For each geometry: one FSD run + SIA post-processing of pencil data
%
% This is the core demonstration of the SIA efficiency: (n_cases + 1)
% simulations per repetition versus the 2*n_cases that a naive FSD approach
% would require.

for rep = 1:n_repetitions

    % --------------------------------------------------------------------
    % 8A) Shared pencil-beam run for this repetition
    % --------------------------------------------------------------------
    progress_step = progress_step + 1;
    if ishandle(wb)
        elapsed_s   = toc(t_loop_start);
        eta_s       = elapsed_s / progress_step * (n_progress_steps - progress_step + 1);
        eta_min     = eta_s / 60;
        waitbar(progress_step / n_progress_steps, wb, ...
                sprintf('Rep %d/%d | Pencil run | ETA ~%.1f min', rep, n_repetitions, eta_min));
    end

    % Run the pencil-beam simulation. The output detPencil contains the
    % detected photon data (positions, pathlengths, exit angles).
    [~, detPencil] = mcxlab(mcx_cfg_pencil);

    % Check that the simulation returned detected photon data.
    has_pencil_data = isfield(detPencil, 'p') && ~isempty(detPencil.p) && ...
                      isfield(detPencil, 'data') && ~isempty(detPencil.data);

    if ~has_pencil_data
        % No photons detected in this repetition; initialise empty arrays.
        rho_pencil_all = zeros(0,1);
        L_pencil_all   = zeros(0,1);
        w_abs_pencil   = zeros(0,1);
        dx_pencil_all  = zeros(0,1);
        dy_pencil_all  = zeros(0,1);
        det_data_all   = zeros(4,0);
    else
        % Apply the NA and global rho acceptance cuts.
        % Photons outside the collection cone or beyond the largest geometry
        % support are discarded here; they cannot contribute to any SIA case.
        [dataPencil, posPencil_vox] = mask_by_radius_and_NA_mm(detPencil, mcx_cfg_pencil.unitinmm, ...
                                                               mcx_cfg_pencil.srcpos(1:2), ...
                                                               theta_acc_max, rho_max_global_mm);

        if isempty(posPencil_vox)
            % All photons fell outside the acceptance cuts.
            rho_pencil_all = zeros(0,1);
            L_pencil_all   = zeros(0,1);
            w_abs_pencil   = zeros(0,1);
            dx_pencil_all  = zeros(0,1);
            dy_pencil_all  = zeros(0,1);
            det_data_all   = zeros(4,0);
        else
            % Convert photon exit positions from voxels to mm.
            posPencil_mm  = double(posPencil_vox) * mcx_cfg_pencil.unitinmm;  % [mm]

            % Lateral displacement from the pencil-beam launch point [mm].
            dx_pencil_all = posPencil_mm(:,1) - source_center_mm(1);  % [mm]
            dy_pencil_all = posPencil_mm(:,2) - source_center_mm(2);  % [mm]

            % Radial exit distance from the source [mm].
            rho_pencil_all = hypot(dx_pencil_all, dy_pencil_all);  % [mm]

            % Total photon pathlength, converted from voxels to mm.
            % dataPencil row 3 contains cumulative partial pathlengths in voxel units.
            L_pencil_all = double(dataPencil(3,:)).' * mcx_cfg_pencil.unitinmm;  % [mm]

            % Absorption weighting: exp(-mu_a * L), consistent with Letter Eq. (3).
            w_abs_pencil = exp(-mu_a .* L_pencil_all);

            % Remove any photons with non-finite values or out-of-range rho/L.
            keep = isfinite(rho_pencil_all) & isfinite(L_pencil_all) & isfinite(w_abs_pencil) & ...
                   (rho_pencil_all >= 0) & (rho_pencil_all <= rho_max_global_mm) & (L_pencil_all >= 0);

            rho_pencil_all = rho_pencil_all(keep);
            L_pencil_all   = L_pencil_all(keep);
            w_abs_pencil   = w_abs_pencil(keep);
            dx_pencil_all  = dx_pencil_all(keep);
            dy_pencil_all  = dy_pencil_all(keep);

            % Pack the relevant photon data into a compact 4-row matrix used
            % by sia_pathlength.m. Rows: [unused; L; dx; dy] in mm.
            det_data_all = zeros(4, numel(L_pencil_all));
            det_data_all(2,:) = L_pencil_all.';   % pathlength [mm]
            det_data_all(3,:) = dx_pencil_all.';  % lateral exit x relative to source [mm]
            det_data_all(4,:) = dy_pencil_all.';  % lateral exit y relative to source [mm]
        end
    end

    % --------------------------------------------------------------------
    % 8B) Per-geometry FSD run + SIA post-processing
    % For each geometry: (1) build and run the explicit FSD simulation,
    % (2) apply the corresponding SIA post-processing to the shared pencil
    % photons. Both produce a reflectance estimate and a pathlength histogram
    % for this repetition.
    % --------------------------------------------------------------------
    for c = 1:n_cases
        CaseNow = Cases{c};

        progress_step = progress_step + 1;
        if ishandle(wb)
            elapsed_s = toc(t_loop_start);
            eta_s     = elapsed_s / progress_step * (n_progress_steps - progress_step + 1);
            eta_min   = eta_s / 60;
            waitbar(progress_step / n_progress_steps, wb, ...
                    sprintf('Rep %d/%d | FSD case %d/%d (%s) | ETA ~%.1f min', ...
                            rep, n_repetitions, c, n_cases, CaseNow.letter_group, eta_min));
        end

        % Retrieve the precomputed rho grid and analytical p(rho) for this case.
        rho_mid_mix   = rho_mid_mix_store{c};     % rho bin centres [mm]
        rho_edges_mix = rho_edges_mix_store{c};   % rho bin edges [mm]
        p_rho_mix     = p_rho_mix_store{c};       % normalised p(rho) for SIA integral
        drho_mix      = rho_edges_mix(2) - rho_edges_mix(1);  % bin width [mm]

        % Ring areas for the SIA rho bins: pi*(r_outer^2 - r_inner^2) [mm^2].
        % Dividing the detected photon weight in each ring by its area gives
        % the pencil-beam reflectance per unit area R(rho_k) [1/mm^2].
        annulus_areas = pi * (rho_edges_mix(2:end).^2 - rho_edges_mix(1:end-1).^2);  % [mm^2]
        annulus_areas = double(annulus_areas(:));

        % Combined weight for the SIA Riemann sum: p(rho_k)*drho / A_annulus(rho_k).
        % This factor converts ring-summed photon weights into geometry-weighted
        % reflectance contributions (Eq. 3 of Letter).
        ring_weight   = (double(p_rho_mix(:)) .* double(drho_mix)) ./ max(annulus_areas, eps);

        % Effective detector area for this case [mm^2].
        A_detector_use = detector_area_store(c);

        % ---- Build the explicit FSD MCX configuration for this geometry ----
        % build_mcx_cfg_fsd.m sets the source type and shape parameters;
        % all other fields are inherited from mcx_cfg_base.
        detector_center_mm = source_center_mm;
        if isfield(CaseNow, 'detector_offset_mm')
            detector_center_mm = source_center_mm + double(CaseNow.detector_offset_mm(:)).';
        end
        mcx_cfg_full = build_mcx_cfg_fsd(mcx_cfg_base, CaseNow);

        % ---- Run the explicit (FSD) simulation ----
        [~, detFull] = mcxlab(mcx_cfg_full);

        % Preallocate pathlength histogram accumulators for this repetition.
        pL_mass_full_rep = zeros(1, n_L_bins);  % FSD pathlength mass histogram
        pL_mass_sia_rep  = zeros(1, n_L_bins);  % SIA pathlength mass histogram

        % ================================================================
        % FSD reflectance and pathlength histogram for this repetition.
        % Apply the appropriate detector mask (shape and NA) to the FSD
        % photons and accumulate the absorption-weighted reflectance and
        % pathlength histogram.
        % ================================================================
        switch CaseNow.type
            case 'uniform_shapes'
                % Mask to disk or annular detector aperture plus NA cone.
                [dataFull, posFull_vox] = mask_by_detector_general_mm(detFull, mcx_cfg_full.unitinmm, ...
                                                                      detector_center_mm, theta_acc_max, ...
                                                                      CaseNow.detector_shape, CaseNow.detector_rin, CaseNow.detector_rout);
                if isempty(posFull_vox)
                    R_full(c, rep) = 0;
                else
                    % Pathlength in mm; absorption weight exp(-mu_a*L).
                    L_full = double(dataFull(3,:)).' * mcx_cfg_full.unitinmm;  % [mm]
                    w_full = exp(-mu_a .* L_full);

                    % Reflectance = total detected absorption-weighted power / launched photons.
                    R_full(c, rep)   = sum(w_full) / n_launched_photons;
                    pL_mass_full_rep = weighted_hist_1d_accum(L_full, w_full, L_edges);
                end

            case 'square_square'
                % Mask to square detector aperture plus NA cone.
                [dataFull, posFull_vox] = mask_by_square_and_NA_mm(detFull, mcx_cfg_full.unitinmm, ...
                                                                   detector_center_mm, CaseNow.detector_side, theta_acc_max);
                if isempty(posFull_vox)
                    R_full(c, rep) = 0;
                else
                    L_full = double(dataFull(3,:)).' * mcx_cfg_full.unitinmm;  % [mm]
                    w_full = exp(-mu_a .* L_full);

                    R_full(c, rep)   = sum(w_full) / n_launched_photons;
                    pL_mass_full_rep = weighted_hist_1d_accum(L_full, w_full, L_edges);
                end

            case 'gaussian_overlap'
                % For the Gaussian case the FSD uses a Gaussian source.
                % Detection is post-processed: photons within the detector
                % aperture radius are weighted by the Gaussian sensitivity
                % profile w_d(r) = exp(-2*r^2/w^2) (Letter Sec. 2C).
                pos_xy_vox = detFull.p(:,1:2);
                pos_xy_mm  = double(pos_xy_vox) * mcx_cfg_full.unitinmm;  % [mm]

                % Radial exit distance from detector centre [mm].
                dx = pos_xy_mm(:,1) - source_center_mm(1);
                dy = pos_xy_mm(:,2) - source_center_mm(2);
                r_det = hypot(dx, dy);
                r_det = r_det(:);   % force column vector

                % Exit angle relative to the surface normal [rad].
                ang = acos(-double(detFull.data(9,:)).');
                ang = ang(:);       % force column vector

                % Accept photons within the aperture radius and NA cone.
                mask = (ang <= theta_acc_max) & (r_det <= CaseNow.R_det_mm);
                mask = mask(:);     % force column logical

                if ~any(mask)
                    R_full(c, rep) = 0;
                else
                    L_full = double(detFull.data(3, mask)).' * mcx_cfg_full.unitinmm;  % [mm]
                    L_full = L_full(:);  % force column

                    w_abs = exp(-mu_a .* L_full);   % absorption weight
                    % Gaussian detection sensitivity weight w_d(r) = exp(-2r^2/w^2).
                    w_det = exp(-2 * (r_det(mask).^2) / (CaseNow.w_mm^2));
                    w_det = w_det(:);   % force column

                    w_full = w_abs .* w_det;  % combined weight

                    R_full(c, rep)   = sum(w_full) / n_launched_photons;
                    pL_mass_full_rep = weighted_hist_1d_accum(L_full, w_full, L_edges);
                end

            otherwise
                error('Unknown case type in FSD evaluation.');
        end

        % ================================================================
        % SIA reflectance and pathlength histogram for this repetition.
        % The pencil-beam photons are reweighted by p(rho)*drho/A_ring to
        % implement the distance-weighted integral (Eq. 3 of Letter).
        % See sia_reflectance.m and sia_pathlength.m for the full derivation.
        % ================================================================
        if isempty(rho_pencil_all)
            R_sia(c, rep) = 0;
        else
            % Call the core SIA reflectance function.
            R_sia(c, rep) = sia_reflectance(rho_pencil_all, w_abs_pencil, ...
                                             rho_edges_mix, ring_weight, ...
                                             A_detector_use, n_launched_photons);

            % Build the geometry-specific SIA pathlength mass histogram.
            pL_mass_sia_rep = sia_pathlength(det_data_all, mu_a, n_launched_photons, ...
                                             A_detector_use, rho_edges_mix, L_edges, ...
                                             p_rho_mix, drho_mix);
        end

        % ================================================================
        % Per-repetition pathlength metrics
        % ================================================================
        % Normalise the raw histogram mass to a proper PDF for this repetition.
        pL_full_rep = normalize_hist_mass_to_pdf(pL_mass_full_rep, L_edges);
        pL_sia_rep  = normalize_hist_mass_to_pdf(pL_mass_sia_rep,  L_edges);

        % Mean pathlength <L> = integral L * p(L) dL, evaluated as a discrete sum.
        meanL_full(c, rep) = sum(L_mid .* pL_full_rep) * dL;  % [mm]
        meanL_sia (c, rep) = sum(L_mid .* pL_sia_rep ) * dL;  % [mm]

        % Wasserstein-1 distance between FSD and SIA pathlength CDFs:
        % W1 = integral |CDF_FSD(L) - CDF_SIA(L)| dL [mm].
        % A small W1 indicates close agreement in the pathlength distributions.
        cdf_full    = cumsum(pL_full_rep) * dL;
        cdf_sia     = cumsum(pL_sia_rep ) * dL;
        W1_L(c, rep) = sum(abs(cdf_full - cdf_sia)) * dL;  % [mm]

        % Pool raw pathlength mass across repetitions. Normalisation into a
        % PDF is applied once at the end (Section 9), giving smoother curves.
        pL_full_mass_total(c, :) = pL_full_mass_total(c, :) + pL_mass_full_rep;
        pL_sia_mass_total (c, :) = pL_sia_mass_total (c, :) + pL_mass_sia_rep;

        % Print a one-line progress summary to the Command Window.
        fprintf('rep %2d/%2d | %s | R_full = %.6g | R_sia = %.6g | <L>_F = %.3f | <L>_S = %.3f | W1 = %.4g\n', ...
            rep, n_repetitions, CaseNow.display_name, R_full(c, rep), R_sia(c, rep), ...
            meanL_full(c, rep), meanL_sia(c, rep), W1_L(c, rep));
    end
end

% Close the waitbar cleanly (the onCleanup object handles crashes too).
if ishandle(wb)
    waitbar(1, wb, 'Done');
end

%% ========================================================================
%  9) FINAL DISPLAY PATHLENGTH PDFS (POOLED OVER REPS)
% ========================================================================
% Normalise the pooled pathlength mass (accumulated across all repetitions)
% to obtain smooth PDFs for display. This is done once here rather than
% averaging per-repetition PDFs.
pL_full_display = zeros(n_cases, n_L_bins);
pL_sia_display  = zeros(n_cases, n_L_bins);

for c = 1:n_cases
    pL_full_display(c, :) = normalize_hist_mass_to_pdf(pL_full_mass_total(c, :), L_edges);
    pL_sia_display (c, :) = normalize_hist_mass_to_pdf(pL_sia_mass_total (c, :), L_edges);
end

%% ========================================================================
% 10) TEXT SUMMARY
% ========================================================================
% Print a summary table of mean reflectance, mean pathlength, and the
% Wasserstein-1 distance between FSD and SIA pathlength distributions.
fprintf('\n====================== SUMMARY ======================\n');
for c = 1:n_cases
    CaseNow = Cases{c};

    mean_R_full = mean(R_full(c, :), 'omitnan');
    std_R_full  = std( R_full(c, :), 0, 'omitnan');
    mean_R_sia  = mean(R_sia(c, :),  'omitnan');
    std_R_sia   = std( R_sia(c, :),  0, 'omitnan');

    mean_L_full = mean(meanL_full(c, :), 'omitnan');
    std_L_full  = std( meanL_full(c, :), 0, 'omitnan');
    mean_L_sia  = mean(meanL_sia(c, :),  'omitnan');
    std_L_sia   = std( meanL_sia(c, :),  0, 'omitnan');

    mean_W1   = mean(W1_L(c, :),   'omitnan');
    median_W1 = median(W1_L(c, :), 'omitnan');

    fprintf('\n%s\n', CaseNow.display_name);
    if CaseNow.appears_in_letter
        fprintf('  Included in Letter scenario %s.\n', CaseNow.letter_group);
    else
        fprintf('  Additional demo-only case (not shown in the Letter figures).\n');
    end
    fprintf('  Reflectance: Full = %.6g ± %.3g | SIA = %.6g ± %.3g | Δmean = %.3g\n', ...
        mean_R_full, std_R_full, mean_R_sia, std_R_sia, mean_R_sia - mean_R_full);
    fprintf('  Mean pathlength <L> [mm]: Full = %.6g ± %.3g | SIA = %.6g ± %.3g | Δmean = %.3g\n', ...
        mean_L_full, std_L_full, mean_L_sia, std_L_sia, mean_L_sia - mean_L_full);
    fprintf('  W1(CDF) [mm]: mean = %.4g | median = %.4g\n', mean_W1, median_W1);

    if strcmp(CaseNow.type, 'gaussian_overlap')
        R_src_cap_mm = CaseNow.R_det_mm + CaseNow.Nw_cap * CaseNow.w_mm;
        fprintf('  Gaussian details: w = %.3f mm | R_det = %.3f mm | R_src_cap = %.3f mm | A_det,eff = %.6g mm^2\n', ...
            CaseNow.w_mm, CaseNow.R_det_mm, R_src_cap_mm, detector_area_store(c));
    end
end

%% ========================================================================
% 11) OPTIONAL RESULT PLOTS
% ========================================================================
if plot_equivalence_results
    % Reflectance histograms: FSD vs SIA distributions for each case.
    % Figure uses a compact tiled layout (see plot_sia_fsd_reflectance_hists.m).
    plot_sia_fsd_reflectance_hists(Cases, R_full, R_sia);

    % Detected pathlength PDFs: FSD vs SIA for each case.
    plot_sia_fsd_pathlength_pdfs(Cases, L_mid, pL_full_display, pL_sia_display);
end

%% ========================================================================
% 12) STRUCTS LEFT IN THE WORKSPACE
% ========================================================================
% These are intentionally left behind so that the user can inspect the
% geometry definitions, precomputed p(rho), reflectance distributions, and
% pooled pathlength PDFs without rerunning the script.

Overview = struct();
Overview.Cases                 = Cases;
Overview.source_points         = overview_source_points;
Overview.detector_points       = overview_detector_points;
Overview.rho_plot              = rho_plot_store;
Overview.p_rho_plot_sampled    = p_rho_plot_sampled_store;
Overview.rho_mix_mid           = rho_mid_mix_store;
Overview.rho_mix_edges         = rho_edges_mix_store;
Overview.p_rho_theory          = p_rho_theory_store;
Overview.p_rho_theory_label    = p_rho_theory_label_store;
Overview.p_rho_mix             = p_rho_mix_store;
Overview.detector_area         = detector_area_store;

Results = struct();
Results.Cases                  = Cases;
Results.R_full                 = R_full;
Results.R_sia                  = R_sia;
Results.meanL_full             = meanL_full;
Results.meanL_sia              = meanL_sia;
Results.W1_L                   = W1_L;
Results.pL_full_display        = pL_full_display;
Results.pL_sia_display         = pL_sia_display;
Results.pL_full_mass_total     = pL_full_mass_total;
Results.pL_sia_mass_total      = pL_sia_mass_total;
Results.L_mid                  = L_mid;

%% ========================================================================
%  Local helper functions (non-plotting)
% ========================================================================

function rho_max = support_rhomax_uniform(CaseNow)
    % Compute the maximum source-detector separation [mm] for a uniform
    % aperture case. This determines the upper bound of the rho grid.
    %
    % For two concentric equal disks the maximum separation is 2*r (diameter).
    % For concentric unequal apertures it is r_source + r_detector.
    % For laterally offset apertures the displacement is added.
    if strcmp(CaseNow.source_shape, 'disk') && strcmp(CaseNow.detector_shape, 'disk') && ...
            all(CaseNow.detector_offset_mm == 0) && (CaseNow.source_rout == CaseNow.detector_rout)
        rho_max = 2 * CaseNow.detector_rout;  % diameter of one disk [mm]
        return;
    end

    if all(CaseNow.detector_offset_mm == 0)
        % Concentric but unequal apertures.
        rho_max = CaseNow.source_rout + CaseNow.detector_rout;
    else
        % Laterally offset: add the centre-to-centre distance.
        rho_max = norm(CaseNow.detector_offset_mm) + (CaseNow.source_rout + CaseNow.detector_rout);
    end
end

function A = detector_area_uniform(shape, rin, rout)
    % Compute the geometric area of a uniform disk or annular detector [mm^2].
    switch shape
        case 'disk'
            A = pi * rout^2;
        case 'annulus'
            A = pi * (rout^2 - rin^2);
        otherwise
            error('Unknown uniform detector shape.');
    end
end

function [data_cut, pos_cut_vox] = mask_by_radius_and_NA_mm(det, unitinmm, srcpos_vox_xy, theta_acc_max, r_max_mm)
    % Retain only photons exiting within a radial distance r_max_mm from the
    % source AND within the NA acceptance cone.
    %
    % Inputs:
    %   det            - MCX detector output struct (fields: p, data)
    %   unitinmm       - voxel size [mm/voxel]
    %   srcpos_vox_xy  - source x,y position in voxel units [voxels]
    %   theta_acc_max  - maximum acceptance half-angle [rad]
    %   r_max_mm       - maximum allowed radial exit distance [mm]
    %
    % Outputs:
    %   data_cut       - detected photon data matrix, masked columns
    %   pos_cut_vox    - photon exit positions in voxel coordinates, masked rows
    if ~isfield(det, 'p') || isempty(det.p)
        data_cut    = det.data(:, false(1,0));
        pos_cut_vox = zeros(0,2);
        return;
    end

    % Lateral exit position in voxel units.
    pos_xy_vox = det.p(:,1:2);
    % Displacement from the source in voxels, converted to mm.
    d_vox = double(pos_xy_vox) - double(srcpos_vox_xy);
    r_mm  = hypot(d_vox(:,1), d_vox(:,2)) * unitinmm;  % [mm]

    % Exit angle from the surface normal: data row 9 stores cos(theta_exit).
    ang  = acos(-double(det.data(9,:)));  % [rad]
    mask = (r_mm <= r_max_mm) & (ang.' <= theta_acc_max);

    if any(mask)
        data_cut    = det.data(:, mask);
        pos_cut_vox = pos_xy_vox(mask, :);
    else
        data_cut    = det.data(:, false(1,0));
        pos_cut_vox = zeros(0,2);
    end
end

function [data_cut, pos_cut_vox] = mask_by_detector_general_mm(det, unitinmm, detector_center_mm, theta_acc_max, detector_shape, detector_rin, detector_rout)
    % Retain only photons exiting within a disk or annular detector aperture
    % AND within the NA acceptance cone.
    %
    % Inputs:
    %   det                 - MCX detector output struct
    %   unitinmm            - voxel size [mm/voxel]
    %   detector_center_mm  - detector centre [x, y] in mm
    %   theta_acc_max       - NA acceptance half-angle [rad]
    %   detector_shape      - 'disk' or 'annulus'
    %   detector_rin        - inner radius [mm] (0 for disk)
    %   detector_rout       - outer radius [mm]
    %
    % Outputs:
    %   data_cut       - detected photon data matrix, masked columns
    %   pos_cut_vox    - photon exit positions in voxel coordinates, masked rows
    if ~isfield(det, 'p') || isempty(det.p)
        data_cut    = det.data(:, false(1,0));
        pos_cut_vox = zeros(0,2);
        return;
    end

    % Convert photon exit positions from voxels to mm.
    pos_xy_vox = det.p(:,1:2);
    pos_xy_mm  = double(pos_xy_vox) * unitinmm;  % [mm]

    % Radial distance from the detector centre [mm].
    dx = pos_xy_mm(:,1) - detector_center_mm(1);
    dy = pos_xy_mm(:,2) - detector_center_mm(2);
    rr = hypot(dx, dy);  % [mm]

    % Apply the aperture shape mask.
    switch detector_shape
        case 'disk'
            inshape = (rr <= detector_rout);
        case 'annulus'
            inshape = (rr >= detector_rin) & (rr <= detector_rout);
        otherwise
            error('Unknown detector shape.');
    end

    ang  = acos(-double(det.data(9,:)));  % exit angle [rad]
    mask = inshape & (ang.' <= theta_acc_max);

    if any(mask)
        data_cut    = det.data(:, mask);
        pos_cut_vox = pos_xy_vox(mask, :);
    else
        data_cut    = det.data(:, false(1,0));
        pos_cut_vox = zeros(0,2);
    end
end

function [data_cut, pos_cut_vox] = mask_by_square_and_NA_mm(det, unitinmm, detector_center_mm, detector_side_mm, theta_acc_max)
    % Retain only photons exiting within a square detector aperture AND
    % within the NA acceptance cone.
    %
    % Inputs:
    %   det                 - MCX detector output struct
    %   unitinmm            - voxel size [mm/voxel]
    %   detector_center_mm  - detector centre [x, y] in mm
    %   detector_side_mm    - side length of the square detector [mm]
    %   theta_acc_max       - NA acceptance half-angle [rad]
    if ~isfield(det, 'p') || isempty(det.p)
        data_cut    = det.data(:, false(1,0));
        pos_cut_vox = zeros(0,2);
        return;
    end

    pos_xy_vox = det.p(:,1:2);
    pos_xy_mm  = double(pos_xy_vox) * unitinmm;  % [mm]

    halfside = detector_side_mm / 2;  % half the side length [mm]
    dx = pos_xy_mm(:,1) - detector_center_mm(1);
    dy = pos_xy_mm(:,2) - detector_center_mm(2);

    % A photon is inside the square if |dx| <= halfside AND |dy| <= halfside.
    inshape = (abs(dx) <= halfside) & (abs(dy) <= halfside);

    ang  = acos(-double(det.data(9,:)));  % exit angle [rad]
    mask = inshape & (ang.' <= theta_acc_max);

    if any(mask)
        data_cut    = det.data(:, mask);
        pos_cut_vox = pos_xy_vox(mask, :);
    else
        data_cut    = det.data(:, false(1,0));
        pos_cut_vox = zeros(0,2);
    end
end

function hist_mass = weighted_hist_1d_accum(x, w, edges)
    % Accumulate a weighted 1D histogram of values x with weights w.
    % Returns the total weight mass per bin (not normalised to a PDF).
    %
    % Inputs:
    %   x      - data values (column or row vector)
    %   w      - corresponding weights (same size as x)
    %   edges  - histogram bin edges (length n_bins + 1)
    %
    % Output:
    %   hist_mass - (1 x n_bins) array of accumulated weighted mass per bin
    hist_mass = zeros(1, numel(edges)-1);
    if isempty(x) || isempty(w)
        return;
    end

    x = double(x(:));
    w = double(w(:));

    % Discard non-finite entries and entries with non-positive weight.
    ok = isfinite(x) & isfinite(w) & (w > 0);
    x  = x(ok);
    w  = w(ok);
    if isempty(x)
        return;
    end

    % Assign each value to a bin index.
    bin = discretize(x, edges);

    % Keep only values that fall within the histogram range.
    ok_bin = ~isnan(bin);
    if ~any(ok_bin)
        return;
    end

    % Sum weights within each bin using accumarray.
    hist_mass = accumarray(double(bin(ok_bin)), double(w(ok_bin)), [numel(edges)-1 1], @sum, 0).';
end

function pdf = normalize_hist_mass_to_pdf(hist_mass, edges)
    % Normalise a raw histogram mass array to a proper probability density
    % function (PDF): pdf integrates to 1 over the bin edges.
    %
    % Inputs:
    %   hist_mass - (1 x n_bins) raw weighted mass per bin
    %   edges     - histogram bin edges (length n_bins + 1)
    %
    % Output:
    %   pdf       - (1 x n_bins) normalised PDF values [1/unit]
    pdf = zeros(1, numel(edges)-1);
    if isempty(hist_mass)
        return;
    end

    hist_mass = double(hist_mass(:)).';
    hist_mass(~isfinite(hist_mass)) = 0;
    hist_mass(hist_mass < 0)        = 0;

    de  = diff(edges(:)).';  % bin widths
    pdf = hist_mass ./ max(de, eps);

    % Normalise so that sum(pdf .* de) = 1.
    Z = sum(pdf .* de);
    if Z > 0
        pdf = pdf / Z;
    end
end



function P = sanitize_pdf(P)
    % Set non-finite and negative values in a PDF array to zero.
    if isempty(P)
        return;
    end
    P = double(P(:)).';
    P(~isfinite(P)) = 0;
    P(P < 0)        = 0;
end

function S = rmfield_safe(S, fieldname)
    % Remove a field from a struct only if it exists (avoids errors).
    if isfield(S, fieldname)
        S = rmfield(S, fieldname);
    end
end

function close_waitbar_safe(wb)
    % Close the waitbar handle if it is still valid (used by onCleanup).
    if ~isempty(wb) && ishandle(wb)
        close(wb);
    end
end
