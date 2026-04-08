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
%  This script demonstrates, in a single place, how the single-integral
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
%  It is included here because it is a simple test
%  case, but it was omitted from the manuscript figures because it does not
%  materially broaden the validation beyond the other concentric examples.
%
%  COMPANION FILES REQUIRED ON THE MATLAB PATH
%  -------------------------------------------
%  Analytical utilities:
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
%
%  DEFAULTS VS LETTER SETTINGS
%  ---------------------------
%  The defaults below are intentionally modest so that the script remains a
%  practical demo. The Letter itself used substantially larger Monte Carlo
%  budgets. Increase the number of photons and repetitions if manuscript-
%  quality statistics are desired.
% ========================================================================

%% ========================================================================
%  0) USER CONTROLS
% ========================================================================
% Plot switches
plot_geometry_overview = true;
plot_equivalence_results = true;

% Case selection (1 = run, 0 = skip)
% Order is defined in Section 3.
run_case = [1 1 1 1 1 1 1];

%% ========================================================================
%  1) OPTICAL / MONTE CARLO SETTINGS
% ========================================================================

% ----------------------- Optical properties ------------------------------
NA    = 0.22;
n_med = 1.35;
theta_acc_max = asin(NA / n_med);

mu_a  = 0.10;      % absorption coefficient [1/mm]
mu_sp = 5.00;      % reduced scattering coefficient [1/mm]
g     = 0.90;      % anisotropy factor used in the demo script
mu_s  = mu_sp / (1 - g);

% -------------------- Monte Carlo repetition settings --------------------
% These defaults are intentionally reduced relative to the Letter so that
% the script remains runnable as a demonstration. Increase them for more
% stable distributions.
n_repetitions = 10;
n_photons     = 1e7;
gpu_id        = '11'; % MCX GPU selection

% ------------------------ Computational domain ---------------------------
% Here we represent a 100 x 100 x 100 mm homogeneous box using a SINGLE voxel.
% Therefore:
%   - the volume is 1 x 1 x 1 voxels
%   - each voxel corresponds to 100 mm
%
% All MCX source sizes are therefore interpreted in voxel units and must be
% obtained by dividing physical dimensions [mm] by unitinmm.
%
% Likewise, detected photon pathlengths and exit coordinates returned by MCX
% are in grid units and must be converted back to mm by multiplying by
% unitinmm.

mm_per_voxel = 100;      % [mm / voxel]
Lxy = 1;
Lz  = 1;
volume = uint8(ones(Lxy, Lxy, Lz));

% Put the source at the center of the top surface of the single voxel.
% With issrcfrom0 = 1, [0.5 0.5 0] corresponds to the top-center location.
source_position_vox = [0.5 0.5 0];
source_direction    = [0 0 1];

% Physical source center in mm
source_center_mm    = double(source_position_vox(1:2)) * mm_per_voxel;

% ------------------------- Pathlength histogram --------------------------
% Detected pathlength distributions are shown after absorption weighting.
L_max_mm = 5;
n_L_bins = 1000;
L_edges  = linspace(0, L_max_mm, n_L_bins + 1);
L_mid    = (L_edges(1:end-1) + L_edges(2:end)) / 2;
dL       = L_edges(2) - L_edges(1);

% ----------------------- rho discretization for SIA ----------------------
% The SIA reflectance integral is evaluated here as a Riemann sum on a
% uniform rho-grid.
n_rho_bins_mix = 5000;

% ------------------- MC sampling only for overview plots -----------------
n_rho_bins_plot  = 50;
n_geometry_plot  = 5000;
n_pairs_plot_pdf = 2e6;

% ---------------------- Gaussian p_eff settings --------------------------
% For the non-uniform Gaussian case, the analytical utility evaluates the
% effective distance distribution numerically. These are the corresponding
% quadrature / transform settings.
NrS_gauss   = 1100;
NrD_gauss   = 1100;
NrInt_gauss = 4000;
Nk_gauss    = 2500;
kmax_gauss  = 350;

%% ========================================================================
%  2) BASE MCX CONFIGURATION
% ========================================================================
% Shared MCX configuration used as the starting point for both the pencil
% run and all explicit FSD runs.

mcx_cfg_base = struct();
mcx_cfg_base.nphoton      = n_photons;
mcx_cfg_base.vol          = volume;
mcx_cfg_base.unitinmm     = mm_per_voxel;

mcx_cfg_base.srcpos       = source_position_vox;
mcx_cfg_base.srcdir       = source_direction;

mcx_cfg_base.autopilot    = 1;
mcx_cfg_base.tstart       = 0;
mcx_cfg_base.tend         = 5e-9;
mcx_cfg_base.tstep        = 5e-9;

mcx_cfg_base.isreflect    = 0;
mcx_cfg_base.respin       = 1;
mcx_cfg_base.issave2pt    = 0;
mcx_cfg_base.savedetflag  = 'dspvx';
mcx_cfg_base.maxdetphoton = n_photons;

mcx_cfg_base.gpuid        = gpu_id;
mcx_cfg_base.issrcfrom0   = 1;
mcx_cfg_base.seed         = -1;

% Finite-NA launch is encoded directly in the photon transport. Under the
% assumptions of the Letter, launch angular effects therefore remain
% embedded in R(rho). Here we set the launch NA
mcx_cfg_base.angleinvcdf  = linspace(0, theta_acc_max/pi, 5);

% Use the top (-z) boundary as the detection plane.
mcx_cfg_base.bc           = 'aa_aaa001000';

% Medium optical properties.
mcx_cfg_base.prop         = [0 0 1 1; mu_a mu_s g n_med];

n_launched_photons = double(mcx_cfg_base.nphoton * mcx_cfg_base.respin);

%% ========================================================================
%  3) GEOMETRY DEFINITIONS
% ========================================================================
% Manuscript mapping:
%   I   = uniform overlapping disk -> disk
%   II  = annulus -> annulus (concentric)
%   III = core -> cladding (concentric)
%   IV  = non-overlapping disks
%   V   = square -> square
%   VI  = Gaussian source + Gaussian-weighted detector
%
% Additional demo-only case:
%   extra = core -> total (concentric)

% --------------------- Circular geometries (mm) --------------------------
R_total = 0.50;   % outer radius [mm]
r_core  = 0.25;   % inner/core radius [mm]

% ------------------- Non-overlapping disk geometry -----------------------
r_source_nonoverlap = 0.30;
r_detector_nonoverlap = 0.45;
sep_nonoverlap = r_source_nonoverlap + r_detector_nonoverlap + 0.25*R_total;

% -------------------------- Square geometry ------------------------------
side_source_square   = 0.50;
side_detector_square = 0.50;
detector_offset_square_mm = [side_source_square, side_source_square];

% ------------------- Non-uniform Gaussian geometry -----------------------
gaussian_w_mm     = 0.50;  % 1/e^2 radius in exp(-2 r^2 / w^2) [mm]
R_detector_gauss  = 0.50;  % detector aperture radius [mm]
Nw_cap            = 10;    % source cap: R_src = R_det + Nw_cap * w

AllCases = { ...
    struct('display_name','I) Disk -> disk (uniform overlap)', ...
           'letter_group','I', ...
           'appears_in_letter',true, ...
           'type','uniform_shapes', ...
           'source_shape','disk',   'source_rin',0,      'source_rout',R_total, 'source_side',[], ...
           'detector_shape','disk', 'detector_rin',0,    'detector_rout',R_total, 'detector_side',[], ...
           'detector_offset_mm',[0 0]), ...

    struct('display_name','Extra) Core -> total (concentric disk -> disk)', ...
           'letter_group','extra', ...
           'appears_in_letter',false, ...
           'type','uniform_shapes', ...
           'source_shape','disk',   'source_rin',0,      'source_rout',r_core, 'source_side',[], ...
           'detector_shape','disk', 'detector_rin',0,    'detector_rout',R_total, 'detector_side',[], ...
           'detector_offset_mm',[0 0]), ...

    struct('display_name','II) Annulus -> annulus (concentric)', ...
           'letter_group','II', ...
           'appears_in_letter',true, ...
           'type','uniform_shapes', ...
           'source_shape','annulus',   'source_rin',r_core, 'source_rout',R_total, 'source_side',[], ...
           'detector_shape','annulus', 'detector_rin',r_core, 'detector_rout',R_total, 'detector_side',[], ...
           'detector_offset_mm',[0 0]), ...

    struct('display_name','III) Core -> cladding (concentric disk -> annulus)', ...
           'letter_group','III', ...
           'appears_in_letter',true, ...
           'type','uniform_shapes', ...
           'source_shape','disk',      'source_rin',0,      'source_rout',r_core, 'source_side',[], ...
           'detector_shape','annulus', 'detector_rin',r_core, 'detector_rout',R_total, 'detector_side',[], ...
           'detector_offset_mm',[0 0]), ...

    struct('display_name','IV) Non-overlapping disks', ...
           'letter_group','IV', ...
           'appears_in_letter',true, ...
           'type','uniform_shapes', ...
           'source_shape','disk',   'source_rin',0, 'source_rout',r_source_nonoverlap, 'source_side',[], ...
           'detector_shape','disk', 'detector_rin',0, 'detector_rout',r_detector_nonoverlap, 'detector_side',[], ...
           'detector_offset_mm',[sep_nonoverlap 0]), ...

    struct('display_name','V) Square -> square (touching corners)', ...
           'letter_group','V', ...
           'appears_in_letter',true, ...
           'type','square_square', ...
           'source_shape','square',   'source_rin',[], 'source_rout',[], 'source_side',side_source_square, ...
           'detector_shape','square', 'detector_rin',[], 'detector_rout',[], 'detector_side',side_detector_square, ...
           'detector_offset_mm',detector_offset_square_mm), ...

    struct('display_name','VI) Gaussian source + Gaussian-weighted detector', ...
           'letter_group','VI', ...
           'appears_in_letter',true, ...
           'type','gaussian_overlap', ...
           'w_mm',gaussian_w_mm, 'R_det_mm',R_detector_gauss, 'Nw_cap',Nw_cap, ...
           'detector_offset_mm',[0 0]) ...
    };

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
%     SIA reflectance and pathlength scaling.
%
%  B) MC-sampled point clouds and MC-sampled distance histograms used ONLY
%     for the overview figures.

overview_source_points    = cell(n_cases,1);
overview_detector_points  = cell(n_cases,1);

rho_mid_mix_store         = cell(n_cases,1);
rho_edges_mix_store       = cell(n_cases,1);
p_rho_mix_store           = cell(n_cases,1);
p_rho_theory_store        = cell(n_cases,1);
p_rho_theory_label_store  = cell(n_cases,1);

rho_plot_store            = cell(n_cases,1);
p_rho_plot_sampled_store  = cell(n_cases,1);

detector_area_store       = nan(n_cases,1);

for c = 1:n_cases
    CaseNow = Cases{c};

    switch CaseNow.type
        case 'uniform_shapes'
            rho_max_mm = support_rhomax_uniform(CaseNow);
            rho_edges_mix = linspace(0, rho_max_mm, n_rho_bins_mix + 1);
            rho_mid_mix   = (rho_edges_mix(1:end-1) + rho_edges_mix(2:end)) / 2;
            drho_mix      = rho_edges_mix(2) - rho_edges_mix(1);

            if strcmp(CaseNow.source_shape,'disk') && strcmp(CaseNow.detector_shape,'disk') && ...
                    all(CaseNow.detector_offset_mm == 0) && (CaseNow.source_rout == CaseNow.detector_rout)

                [rho_plot, pdf_samp_plot, source_pts_plot, detector_pts_plot] = ...
                    sample_p_rho_uniform_overlap_mc(2*CaseNow.detector_rout, n_pairs_plot_pdf, n_geometry_plot, n_rho_bins_plot);
                p_mix = p_rho_uniform_overlap_analytic(rho_mid_mix, 2*CaseNow.detector_rout);
                p_theory = p_mix;
                p_label = 'analytic';

            elseif strcmp(CaseNow.source_shape,'disk') && strcmp(CaseNow.detector_shape,'disk') && ...
                    all(CaseNow.detector_offset_mm == 0)

                [rho_plot, pdf_samp_plot, source_pts_plot, detector_pts_plot] = ...
                    sample_p_rho_core_total_mc(CaseNow.source_rout, CaseNow.detector_rout, ...
                                               n_pairs_plot_pdf, n_geometry_plot, n_rho_bins_plot);
                p_mix = p_rho_core_total_analytic(rho_mid_mix, CaseNow.source_rout, CaseNow.detector_rout);
                p_theory = p_mix;
                p_label = 'analytic';

            elseif strcmp(CaseNow.source_shape,'annulus') && strcmp(CaseNow.detector_shape,'annulus') && ...
                    all(CaseNow.detector_offset_mm == 0)

                [rho_plot, pdf_samp_plot, source_pts_plot, detector_pts_plot] = ...
                    sample_p_rho_annulus_annulus_mc(CaseNow.source_rin, CaseNow.source_rout, ...
                                                    CaseNow.detector_rin, CaseNow.detector_rout, ...
                                                    n_pairs_plot_pdf, n_geometry_plot, n_rho_bins_plot);
                p_mix = p_rho_annulus_annulus_analytic(rho_mid_mix, ...
                                                       CaseNow.source_rin, CaseNow.source_rout, ...
                                                       CaseNow.detector_rin, CaseNow.detector_rout);
                p_theory = p_mix;
                p_label = 'analytic';

            elseif strcmp(CaseNow.source_shape,'disk') && strcmp(CaseNow.detector_shape,'annulus') && ...
                    all(CaseNow.detector_offset_mm == 0)

                [rho_plot, pdf_samp_plot, source_pts_plot, detector_pts_plot] = ...
                    sample_p_rho_disk_annulus_mc(CaseNow.source_rout, CaseNow.detector_rin, CaseNow.detector_rout, ...
                                                 n_pairs_plot_pdf, n_geometry_plot, n_rho_bins_plot);
                p_mix = p_rho_disk_annulus_analytic(rho_mid_mix, ...
                                                    CaseNow.source_rout, CaseNow.detector_rin, CaseNow.detector_rout);
                p_theory = p_mix;
                p_label = 'analytic';

            elseif strcmp(CaseNow.source_shape,'disk') && strcmp(CaseNow.detector_shape,'disk') && ...
                    any(CaseNow.detector_offset_mm ~= 0)

                sep_mm = norm(CaseNow.detector_offset_mm);
                [rho_plot, pdf_samp_plot, source_pts_plot, detector_pts_plot] = ...
                    sample_p_rho_nonoverlap_disks_mc(CaseNow.source_rout, CaseNow.detector_rout, sep_mm, ...
                                                     n_pairs_plot_pdf, n_geometry_plot, n_rho_bins_plot);
                p_mix = p_rho_nonoverlap_disks_analytic(rho_mid_mix, sep_mm, CaseNow.source_rout, CaseNow.detector_rout);
                p_theory = p_mix;
                p_label = 'semi-analytic (1D integral)';

            else
                error('Unsupported uniform case encountered during precomputation.');
            end

            p_mix = sanitize_pdf(p_mix);
            Zmix = sum(p_mix) * drho_mix;
            if Zmix > 0
                p_mix = p_mix / Zmix;
            end

            detector_area_store(c) = detector_area_uniform(CaseNow.detector_shape, ...
                                                           CaseNow.detector_rin, CaseNow.detector_rout);

        case 'square_square'
            rho_max_mm = hypot(CaseNow.detector_offset_mm(1) + (CaseNow.source_side/2 + CaseNow.detector_side/2), ...
                               CaseNow.detector_offset_mm(2) + (CaseNow.source_side/2 + CaseNow.detector_side/2));

            rho_edges_mix = linspace(0, rho_max_mm, n_rho_bins_mix + 1);
            rho_mid_mix   = (rho_edges_mix(1:end-1) + rho_edges_mix(2:end)) / 2;
            drho_mix      = rho_edges_mix(2) - rho_edges_mix(1);

            [rho_plot, pdf_samp_plot, source_pts_plot, detector_pts_plot] = ...
                sample_p_rho_square_square_mc(CaseNow.source_side, CaseNow.detector_side, ...
                                              CaseNow.detector_offset_mm(1), CaseNow.detector_offset_mm(2), ...
                                              n_pairs_plot_pdf, n_geometry_plot, n_rho_bins_plot);

            p_mix = p_rho_square_square_analytic(rho_mid_mix, ...
                                                 CaseNow.source_side, CaseNow.detector_side, ...
                                                 CaseNow.detector_offset_mm(1), CaseNow.detector_offset_mm(2));
            p_mix = sanitize_pdf(p_mix);
            Zmix = sum(p_mix) * drho_mix;
            if Zmix > 0
                p_mix = p_mix / Zmix;
            end

            p_theory = p_mix;
            p_label  = 'semi-analytic (covariogram)';
            detector_area_store(c) = CaseNow.detector_side^2;

        case 'gaussian_overlap'
            w_mm = CaseNow.w_mm;
            R_det_mm = CaseNow.R_det_mm;
            R_src_cap_mm = R_det_mm + CaseNow.Nw_cap * w_mm;
            rho_max_mm = R_src_cap_mm + R_det_mm;

            rho_edges_mix = linspace(0, rho_max_mm, n_rho_bins_mix + 1);
            rho_mid_mix   = (rho_edges_mix(1:end-1) + rho_edges_mix(2:end)) / 2;

            [rho_plot, pdf_samp_plot, source_pts_plot, detector_pts_plot] = ...
                sample_p_eff_gaussian_overlap_mc(w_mm, R_src_cap_mm, R_det_mm, ...
                                                 n_pairs_plot_pdf, n_geometry_plot, n_rho_bins_plot);

            gaussian_curves = p_eff_gaussian_overlap_analytic(rho_edges_mix, w_mm, R_src_cap_mm, R_det_mm, ...
                                                              NrS_gauss, NrD_gauss, NrInt_gauss, Nk_gauss, kmax_gauss);
            p_mix = gaussian_curves.p_base;
            p_theory = p_mix;
            p_label = sprintf('analytic integral (R_src = R_det + %dw)', CaseNow.Nw_cap);

            % For non-uniform detection, p_eff(rho) is normalized and the
            % overall detection strength is carried separately through the
            % effective detector area A_det,eff = integral w_d(r) dA.
            detector_area_store(c) = (pi*w_mm^2/2) * (1 - exp(-2*R_det_mm^2 / w_mm^2));

        otherwise
            error('Unknown case type encountered during precomputation.');
    end

    overview_source_points{c}   = source_pts_plot(1:min(n_geometry_plot, size(source_pts_plot,1)), :);
    overview_detector_points{c} = detector_pts_plot(1:min(n_geometry_plot, size(detector_pts_plot,1)), :);

    rho_mid_mix_store{c}        = rho_mid_mix;
    rho_edges_mix_store{c}      = rho_edges_mix;
    p_rho_mix_store{c}          = p_mix;
    p_rho_theory_store{c}       = sanitize_pdf(p_theory);
    p_rho_theory_label_store{c} = p_label;

    rho_plot_store{c}           = rho_plot;
    p_rho_plot_sampled_store{c} = pdf_samp_plot;
end

% The shared pencil-beam run only needs to retain photons up to the largest
% rho required by any selected geometry.
rho_max_global_mm = 0;
for c = 1:n_cases
    rho_max_global_mm = max(rho_max_global_mm, rho_edges_mix_store{c}(end));
end

%% ========================================================================
%  5) OPTIONAL OVERVIEW PLOTS
% ========================================================================
if plot_geometry_overview
    plot_sia_fsd_geometries(Cases, overview_source_points, overview_detector_points);
    plot_sia_fsd_prho_overview(Cases, rho_plot_store, p_rho_plot_sampled_store, ...
                               rho_mid_mix_store, p_rho_theory_store, p_rho_theory_label_store);
end

%% ========================================================================
%  6) PREALLOCATE OUTPUTS
% ========================================================================
% Reflectance is stored per repetition so that FSD and SIA distributions can
% be compared directly.
R_full = nan(n_cases, n_repetitions);
R_sia  = nan(n_cases, n_repetitions);

% Per-repetition pathlength metrics.
meanL_full = nan(n_cases, n_repetitions);
meanL_sia  = nan(n_cases, n_repetitions);
W1_L       = nan(n_cases, n_repetitions);

% For the displayed pathlength curves we pool ALL weighted pathlength mass
% across repetitions and normalize only once at the end.
pL_full_mass_total = zeros(n_cases, n_L_bins);
pL_sia_mass_total  = zeros(n_cases, n_L_bins);

%% ========================================================================
%  7) BUILD THE SHARED PENCIL-BEAM CONFIGURATION
% ========================================================================
% This is the reusable MC simulation from which all SIA estimates are
% derived. Within a given repetition it is run exactly once.

mcx_cfg_pencil = mcx_cfg_base;
mcx_cfg_pencil.srctype = 'pencil';
mcx_cfg_pencil = rmfield_safe(mcx_cfg_pencil, 'srcparam1');
mcx_cfg_pencil = rmfield_safe(mcx_cfg_pencil, 'srcparam2');

n_progress_steps = n_repetitions * (n_cases + 1);  % 1 pencil + n_cases FSD runs per repetition
progress_step = 0;

wb = waitbar(0, 'Starting...', 'Name', 'MCX progress');
cleanupObj = onCleanup(@() close_waitbar_safe(wb)); %#ok<NASGU>

%% ========================================================================
%  8) MAIN LOOP
% ========================================================================
% Loop structure:
%   outer loop  = Monte Carlo repetitions
%   per repeat  = ONE pencil-beam run + one explicit FSD run per geometry
%
% This structure is essential: the SIA reuses the SAME pencil-beam photons
% across all geometries in a repetition, whereas the FSD must explicitly run
% each source-detector configuration separately.

for rep = 1:n_repetitions

    % --------------------------------------------------------------------
    % 8A) Shared pencil-beam run for this repetition
    % --------------------------------------------------------------------
    progress_step = progress_step + 1;
    if ishandle(wb)
        waitbar(progress_step / n_progress_steps, wb, ...
                sprintf('Rep %d/%d | Shared pencil run', rep, n_repetitions));
    end

    [~, detPencil] = mcxlab(mcx_cfg_pencil);

    % Retain only photons that:
    %   - exit within the collection cone corresponding to the selected NA
    %   - fall within the largest rho support needed by any selected geometry
    % All further geometry dependence is then introduced only in post-
    % processing through p_SD(rho) or p_eff(rho).
    has_pencil_data = isfield(detPencil, 'p') && ~isempty(detPencil.p) && ...
                      isfield(detPencil, 'data') && ~isempty(detPencil.data);

    if ~has_pencil_data
        rho_pencil_all = zeros(0,1);
        L_pencil_all   = zeros(0,1);
        w_abs_pencil   = zeros(0,1);
        dx_pencil_all  = zeros(0,1);
        dy_pencil_all  = zeros(0,1);
        det_data_all   = zeros(4,0);
    else
        [dataPencil, posPencil_vox] = mask_by_radius_and_NA_mm(detPencil, mcx_cfg_pencil.unitinmm, ...
                                                               mcx_cfg_pencil.srcpos(1:2), ...
                                                               theta_acc_max, rho_max_global_mm);

        if isempty(posPencil_vox)
            rho_pencil_all = zeros(0,1);
            L_pencil_all   = zeros(0,1);
            w_abs_pencil   = zeros(0,1);
            dx_pencil_all  = zeros(0,1);
            dy_pencil_all  = zeros(0,1);
            det_data_all   = zeros(4,0);
        else
            posPencil_mm = double(posPencil_vox) * mcx_cfg_pencil.unitinmm;
            dx_pencil_all = posPencil_mm(:,1) - source_center_mm(1);
            dy_pencil_all = posPencil_mm(:,2) - source_center_mm(2);

            rho_pencil_all = hypot(dx_pencil_all, dy_pencil_all);
            L_pencil_all   = double(dataPencil(3,:)).' * mcx_cfg_pencil.unitinmm;

            % Absorption weighting used throughout the Letter.
            w_abs_pencil = exp(-mu_a .* L_pencil_all);

            keep = isfinite(rho_pencil_all) & isfinite(L_pencil_all) & isfinite(w_abs_pencil) & ...
                   (rho_pencil_all >= 0) & (rho_pencil_all <= rho_max_global_mm) & (L_pencil_all >= 0);

            rho_pencil_all = rho_pencil_all(keep);
            L_pencil_all   = L_pencil_all(keep);
            w_abs_pencil   = w_abs_pencil(keep);
            dx_pencil_all  = dx_pencil_all(keep);
            dy_pencil_all  = dy_pencil_all(keep);

            % Compact representation used later for SIA pathlength mixing.
            % Row 2 = pathlength L, rows 3-4 = exit coordinates relative to
            % the pencil-beam launch point.
            det_data_all = zeros(4, numel(L_pencil_all));
            det_data_all(2,:) = L_pencil_all.';
            det_data_all(3,:) = dx_pencil_all.';
            det_data_all(4,:) = dy_pencil_all.';
        end
    end

    % --------------------------------------------------------------------
    % 8B) Per-geometry FSD run + SIA post-processing
    % --------------------------------------------------------------------
    for c = 1:n_cases
        CaseNow = Cases{c};

        progress_step = progress_step + 1;
        if ishandle(wb)
            waitbar(progress_step / n_progress_steps, wb, ...
                    sprintf('Rep %d/%d | FSD case %d/%d', rep, n_repetitions, c, n_cases));
        end

        rho_mid_mix   = rho_mid_mix_store{c};
        rho_edges_mix = rho_edges_mix_store{c};
        p_rho_mix     = p_rho_mix_store{c};
        drho_mix      = rho_edges_mix(2) - rho_edges_mix(1);

        % Ring areas associated with the rho discretization. The pencil-beam
        % run provides detected power per ring; dividing by the ring area
        % yields the per-unit-area reflectance needed in the SIA integral.
        annulus_areas = pi * (rho_edges_mix(2:end).^2 - rho_edges_mix(1:end-1).^2);
        annulus_areas = double(annulus_areas(:));
        ring_weight   = (double(p_rho_mix(:)) .* double(drho_mix)) ./ max(annulus_areas, eps);

        A_detector_use = detector_area_store(c);

        % ---------------- Build explicit FSD MCX source ------------------
        mcx_cfg_full = mcx_cfg_base;
        detector_center_mm = source_center_mm;
        if isfield(CaseNow, 'detector_offset_mm')
            detector_center_mm = source_center_mm + double(CaseNow.detector_offset_mm(:)).';
        end

        switch CaseNow.type
            case 'uniform_shapes'
                mcx_cfg_full.srctype = 'disk';
                mcx_cfg_full = rmfield_safe(mcx_cfg_full, 'srcparam1');
                mcx_cfg_full = rmfield_safe(mcx_cfg_full, 'srcparam2');

                if strcmp(CaseNow.source_shape, 'disk')
                    mcx_cfg_full.srcparam1 = [CaseNow.source_rout, 0] / mcx_cfg_full.unitinmm;
                elseif strcmp(CaseNow.source_shape, 'annulus')
                    mcx_cfg_full.srcparam1 = [CaseNow.source_rout, CaseNow.source_rin] / mcx_cfg_full.unitinmm;
                else
                    error('Unsupported source shape in uniform FSD build.');
                end

            case 'square_square'
                mcx_cfg_full.srctype = 'planar';
                mcx_cfg_full = rmfield_safe(mcx_cfg_full, 'srcparam1');
                mcx_cfg_full = rmfield_safe(mcx_cfg_full, 'srcparam2');

                side_vox = CaseNow.source_side / mcx_cfg_full.unitinmm;
                corner_vox = [double(mcx_cfg_base.srcpos(1)) - side_vox/2, ...
                              double(mcx_cfg_base.srcpos(2)) - side_vox/2, ...
                              double(mcx_cfg_base.srcpos(3))];
                mcx_cfg_full.srcpos    = corner_vox;
                mcx_cfg_full.srcparam1 = [side_vox 0 0 0];
                mcx_cfg_full.srcparam2 = [0 side_vox 0 0];
                mcx_cfg_full.srcdir    = [0 0 1];

            case 'gaussian_overlap'
                mcx_cfg_full.srctype   = 'gaussian';
                mcx_cfg_full.srcdir    = [0 0 1 0];
                mcx_cfg_full.srcparam1 = [CaseNow.w_mm / mcx_cfg_full.unitinmm 0 0 0];
                mcx_cfg_full.srcparam2 = [0 0 0 0];

            otherwise
                error('Unknown case type in FSD source construction.');
        end

        % ---------------------- Run explicit FSD MC ----------------------
        [~, detFull] = mcxlab(mcx_cfg_full);

        pL_mass_full_rep = zeros(1, n_L_bins);
        pL_mass_sia_rep  = zeros(1, n_L_bins);

        % ================================================================
        % FSD reflectance and pathlength histogram for this repetition
        % ================================================================
        switch CaseNow.type
            case 'uniform_shapes'
                [dataFull, posFull_vox] = mask_by_detector_general_mm(detFull, mcx_cfg_full.unitinmm, ...
                                                                      detector_center_mm, theta_acc_max, ...
                                                                      CaseNow.detector_shape, CaseNow.detector_rin, CaseNow.detector_rout);
                if isempty(posFull_vox)
                    R_full(c, rep) = 0;
                else
                    L_full = double(dataFull(3,:)).' * mcx_cfg_full.unitinmm;
                    w_full = exp(-mu_a .* L_full);

                    R_full(c, rep) = sum(w_full) / n_launched_photons;
                    pL_mass_full_rep = weighted_hist_1d_accum(L_full, w_full, L_edges);
                end

            case 'square_square'
                [dataFull, posFull_vox] = mask_by_square_and_NA_mm(detFull, mcx_cfg_full.unitinmm, ...
                                                                   detector_center_mm, CaseNow.detector_side, theta_acc_max);
                if isempty(posFull_vox)
                    R_full(c, rep) = 0;
                else
                    L_full = double(dataFull(3,:)).' * mcx_cfg_full.unitinmm;
                    w_full = exp(-mu_a .* L_full);

                    R_full(c, rep) = sum(w_full) / n_launched_photons;
                    pL_mass_full_rep = weighted_hist_1d_accum(L_full, w_full, L_edges);
                end
            case 'gaussian_overlap'
                pos_xy_vox = detFull.p(:,1:2);
                pos_xy_mm  = double(pos_xy_vox) * mcx_cfg_full.unitinmm;
                
                dx = pos_xy_mm(:,1) - source_center_mm(1);
                dy = pos_xy_mm(:,2) - source_center_mm(2);
                r_det = hypot(dx, dy);
                r_det = r_det(:);   % force column
                
                ang = acos(-double(detFull.data(9,:)).');
                ang = ang(:);       % force column
                
                mask = (ang <= theta_acc_max) & (r_det <= CaseNow.R_det_mm);
                mask = mask(:);     % force column logical
                
                if ~any(mask)
                    R_full(c, rep) = 0;
                else
                    L_full = double(detFull.data(3, mask)).' * mcx_cfg_full.unitinmm;
                    L_full = L_full(:);  % force column
                
                    w_abs = exp(-mu_a .* L_full);
                    w_det = exp(-2 * (r_det(mask).^2) / (CaseNow.w_mm^2));
                    w_det = w_det(:);    % force column
                
                    w_full = w_abs .* w_det;
                
                    R_full(c, rep) = sum(w_full) / n_launched_photons;
                    pL_mass_full_rep = weighted_hist_1d_accum(L_full, w_full, L_edges);
                end
        
            otherwise
                error('Unknown case type in FSD evaluation.');
        end

        % ================================================================
        % SIA reflectance and pathlength histogram for this repetition
        % ================================================================
        if isempty(rho_pencil_all)
            R_sia(c, rep) = 0;
        else
            % Reflectance via discretized single-integral approximation:
            %
            %   R ~= A_d * sum_k [ R(rho_k) * p(rho_k) * drho ]
            %
            % where the pencil-beam simulation provides an estimate of
            % ring-averaged R(rho_k) through the detected absorption-weighted
            % photon sums divided by the ring area.
            bin_rho = discretize(rho_pencil_all, rho_edges_mix);
            ok_rho  = ~isnan(bin_rho);

            if ~any(ok_rho)
                R_sia(c, rep) = 0;
            else
                ring_sums = accumarray(double(bin_rho(ok_rho)), double(w_abs_pencil(ok_rho)), ...
                                       [numel(rho_mid_mix), 1], @sum, 0);
                R_sia(c, rep) = (A_detector_use / n_launched_photons) * sum(double(ring_sums(:)) .* double(ring_weight(:)));
            end

            % Build the absorption-weighted geometry-specific pathlength mass
            % by reweighting all pencil-beam photons with the same SIA ring
            % factors used for the reflectance integral.
            pL_mass_sia_rep = weighted_hist_L_sia_accum(det_data_all, mu_a, n_launched_photons, ...
                                                        A_detector_use, rho_edges_mix, L_edges, p_rho_mix, drho_mix);
        end

        % ================================================================
        % Per-repetition pathlength metrics
        % ================================================================
        pL_full_rep = normalize_hist_mass_to_pdf(pL_mass_full_rep, L_edges);
        pL_sia_rep  = normalize_hist_mass_to_pdf(pL_mass_sia_rep,  L_edges);

        meanL_full(c, rep) = sum(L_mid .* pL_full_rep) * dL;
        meanL_sia (c, rep) = sum(L_mid .* pL_sia_rep ) * dL;

        cdf_full = cumsum(pL_full_rep) * dL;
        cdf_sia  = cumsum(pL_sia_rep ) * dL;
        W1_L(c, rep) = sum(abs(cdf_full - cdf_sia)) * dL;

        % Pool mass across repetitions for the final displayed pathlength
        % curves. This uses all individual photon contributions together.
        pL_full_mass_total(c, :) = pL_full_mass_total(c, :) + pL_mass_full_rep;
        pL_sia_mass_total (c, :) = pL_sia_mass_total (c, :) + pL_mass_sia_rep;

        fprintf('rep %2d/%2d | %s | R_full = %.6g | R_sia = %.6g | <L>_F = %.3f | <L>_S = %.3f | W1 = %.4g\n', ...
            rep, n_repetitions, CaseNow.display_name, R_full(c, rep), R_sia(c, rep), ...
            meanL_full(c, rep), meanL_sia(c, rep), W1_L(c, rep));
    end
end

if ishandle(wb)
    waitbar(1, wb, 'Done');
end

%% ========================================================================
%  9) FINAL DISPLAY PATHLENGTH PDFS (POOLED OVER REPS)
% ========================================================================
pL_full_display = zeros(n_cases, n_L_bins);
pL_sia_display  = zeros(n_cases, n_L_bins);

for c = 1:n_cases
    pL_full_display(c, :) = normalize_hist_mass_to_pdf(pL_full_mass_total(c, :), L_edges);
    pL_sia_display (c, :) = normalize_hist_mass_to_pdf(pL_sia_mass_total (c, :), L_edges);
end

%% ========================================================================
% 10) TEXT SUMMARY
% ========================================================================
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

    mean_W1     = mean(W1_L(c, :),   'omitnan');
    median_W1   = median(W1_L(c, :), 'omitnan');

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
    plot_sia_fsd_reflectance_hists(Cases, R_full, R_sia);
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
Results.L_mid                  = L_mid;
Results.pL_full_display        = pL_full_display;
Results.pL_sia_display         = pL_sia_display;
Results.pL_full_mass_total     = pL_full_mass_total;
Results.pL_sia_mass_total      = pL_sia_mass_total;

%% ========================================================================
%  Local helper functions (non-plotting)
% ========================================================================

function rho_max = support_rhomax_uniform(CaseNow)
    if strcmp(CaseNow.source_shape, 'disk') && strcmp(CaseNow.detector_shape, 'disk') && ...
            all(CaseNow.detector_offset_mm == 0) && (CaseNow.source_rout == CaseNow.detector_rout)
        rho_max = 2 * CaseNow.detector_rout;
        return;
    end

    if all(CaseNow.detector_offset_mm == 0)
        rho_max = CaseNow.source_rout + CaseNow.detector_rout;
    else
        rho_max = norm(CaseNow.detector_offset_mm) + (CaseNow.source_rout + CaseNow.detector_rout);
    end
end

function A = detector_area_uniform(shape, rin, rout)
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
    if ~isfield(det, 'p') || isempty(det.p)
        data_cut = det.data(:, false(1,0));
        pos_cut_vox = zeros(0,2);
        return;
    end

    pos_xy_vox = det.p(:,1:2);
    d_vox = double(pos_xy_vox) - double(srcpos_vox_xy);
    r_mm = hypot(d_vox(:,1), d_vox(:,2)) * unitinmm;

    ang = acos(-double(det.data(9,:)));
    mask = (r_mm <= r_max_mm) & (ang.' <= theta_acc_max);

    if any(mask)
        data_cut = det.data(:, mask);
        pos_cut_vox = pos_xy_vox(mask, :);
    else
        data_cut = det.data(:, false(1,0));
        pos_cut_vox = zeros(0,2);
    end
end

function [data_cut, pos_cut_vox] = mask_by_detector_general_mm(det, unitinmm, detector_center_mm, theta_acc_max, detector_shape, detector_rin, detector_rout)
    if ~isfield(det, 'p') || isempty(det.p)
        data_cut = det.data(:, false(1,0));
        pos_cut_vox = zeros(0,2);
        return;
    end

    pos_xy_vox = det.p(:,1:2);
    pos_xy_mm  = double(pos_xy_vox) * unitinmm;

    dx = pos_xy_mm(:,1) - detector_center_mm(1);
    dy = pos_xy_mm(:,2) - detector_center_mm(2);
    rr = hypot(dx, dy);

    switch detector_shape
        case 'disk'
            inshape = (rr <= detector_rout);
        case 'annulus'
            inshape = (rr >= detector_rin) & (rr <= detector_rout);
        otherwise
            error('Unknown detector shape.');
    end

    ang = acos(-double(det.data(9,:)));
    mask = inshape & (ang.' <= theta_acc_max);

    if any(mask)
        data_cut = det.data(:, mask);
        pos_cut_vox = pos_xy_vox(mask, :);
    else
        data_cut = det.data(:, false(1,0));
        pos_cut_vox = zeros(0,2);
    end
end

function [data_cut, pos_cut_vox] = mask_by_square_and_NA_mm(det, unitinmm, detector_center_mm, detector_side_mm, theta_acc_max)
    if ~isfield(det, 'p') || isempty(det.p)
        data_cut = det.data(:, false(1,0));
        pos_cut_vox = zeros(0,2);
        return;
    end

    pos_xy_vox = det.p(:,1:2);
    pos_xy_mm  = double(pos_xy_vox) * unitinmm;

    halfside = detector_side_mm / 2;
    dx = pos_xy_mm(:,1) - detector_center_mm(1);
    dy = pos_xy_mm(:,2) - detector_center_mm(2);

    inshape = (abs(dx) <= halfside) & (abs(dy) <= halfside);

    ang = acos(-double(det.data(9,:)));
    mask = inshape & (ang.' <= theta_acc_max);

    if any(mask)
        data_cut = det.data(:, mask);
        pos_cut_vox = pos_xy_vox(mask, :);
    else
        data_cut = det.data(:, false(1,0));
        pos_cut_vox = zeros(0,2);
    end
end

function hist_mass = weighted_hist_1d_accum(x, w, edges)
    hist_mass = zeros(1, numel(edges)-1);
    if isempty(x) || isempty(w)
        return;
    end

    x = double(x(:));
    w = double(w(:));

    ok = isfinite(x) & isfinite(w) & (w > 0);
    x = x(ok);
    w = w(ok);
    if isempty(x)
        return;
    end

    bin = discretize(x, edges);
    ok_bin = ~isnan(bin);
    if ~any(ok_bin)
        return;
    end

    hist_mass = accumarray(double(bin(ok_bin)), double(w(ok_bin)), [numel(edges)-1 1], @sum, 0).';
end

function pdf = normalize_hist_mass_to_pdf(hist_mass, edges)
    pdf = zeros(1, numel(edges)-1);
    if isempty(hist_mass)
        return;
    end

    hist_mass = double(hist_mass(:)).';
    hist_mass(~isfinite(hist_mass)) = 0;
    hist_mass(hist_mass < 0) = 0;

    de = diff(edges(:)).';
    pdf = hist_mass ./ max(de, eps);

    Z = sum(pdf .* de);
    if Z > 0
        pdf = pdf / Z;
    end
end

function hist_mass = weighted_hist_L_sia_accum(det_data, mu_a, launched_photons, A_detector_use, edges_rho, edges_L, p_rho_mix, drho_mix)
    %#ok<INUSD> launched_photons is carried for interface clarity.
    hist_mass = zeros(1, numel(edges_L)-1);

    if isempty(det_data) || isempty(p_rho_mix)
        return;
    end

    L  = double(det_data(2,:)).';
    dx = double(det_data(3,:)).';
    dy = double(det_data(4,:)).';
    rho = hypot(dx, dy);

    ok = isfinite(L) & isfinite(rho) & (L >= 0) & (rho >= 0);
    L   = L(ok);
    rho = rho(ok);

    if isempty(L)
        return;
    end

    w_abs = exp(-mu_a .* L);

    annulus_areas = pi * (edges_rho(2:end).^2 - edges_rho(1:end-1).^2);
    annulus_areas = double(annulus_areas(:));
    ring_weight   = (double(p_rho_mix(:)) .* double(drho_mix)) ./ max(annulus_areas, eps);

    bin_rho = discretize(rho, edges_rho);
    ok_rho  = ~isnan(bin_rho);
    if ~any(ok_rho)
        return;
    end

    w_sia = A_detector_use * w_abs(ok_rho) .* ring_weight(bin_rho(ok_rho));

    bin_L = discretize(L(ok_rho), edges_L);
    ok_L  = ~isnan(bin_L);
    if ~any(ok_L)
        return;
    end

    hist_mass = accumarray(double(bin_L(ok_L)), double(w_sia(ok_L)), [numel(edges_L)-1 1], @sum, 0).';
end

function P = sanitize_pdf(P)
    if isempty(P)
        return;
    end
    P = double(P(:)).';
    P(~isfinite(P)) = 0;
    P(P < 0) = 0;
end

function S = rmfield_safe(S, fieldname)
    if isfield(S, fieldname)
        S = rmfield(S, fieldname);
    end
end

function close_waitbar_safe(wb)
    if ~isempty(wb) && ishandle(wb)
        close(wb);
    end
end
