function hist_mass = sia_pathlength(det_data, mu_a, n_launched, A_detector, ...
                                     rho_edges, L_edges, p_rho, drho)
% SIA_PATHLENGTH  Compute the SIA-weighted pathlength mass histogram for a
%                 given source-detector geometry.
%
%  hist_mass = sia_pathlength(det_data, mu_a, n_launched, A_detector, ...
%                              rho_edges, L_edges, p_rho, drho)
%
%  This function implements the pathlength analogue of the SIA reflectance
%  integral, corresponding to Eq. (5) of the Letter:
%
%      P(L) = integral_0^rho_max  f(L | rho) * p_SD(rho) d_rho
%
%  where f(L | rho) is the conditional (absorption-weighted) pathlength
%  distribution for photons exiting at radius rho, and p_SD(rho) is the
%  source-detector distance distribution.
%
%  DISCRETISATION
%  --------------
%  The integral is evaluated as a double Riemann sum: first the pencil-beam
%  photons are sorted into rho bins (same grid as for reflectance), and then
%  within each rho bin their absorption-weighted pathlengths are binned into
%  a 1D histogram. The histogram bins are then combined across all rho bins
%  using the same ring_weight factors as in sia_reflectance:
%
%      hist_mass(L_bin) ~= A_d * sum_k  ring_weight(k)
%                               * sum_{photons in rho bin k, L bin j} w_abs_i
%
%  where ring_weight(k) = p_SD(rho_k) * drho / A_ring(rho_k).
%
%  The result is a raw (unnormalised) pathlength mass histogram. Normalise
%  it with normalize_hist_mass_to_pdf (in the main script) to obtain a PDF.
%
%  INPUTS
%  ------
%  det_data    - (4 x N) matrix of detected pencil-beam photon data:
%                  row 1: (unused placeholder for compatibility)
%                  row 2: pathlength L [mm]
%                  row 3: lateral exit coordinate dx relative to source [mm]
%                  row 4: lateral exit coordinate dy relative to source [mm]
%                Only photons passing the NA acceptance cut should be included.
%
%  mu_a        - absorption coefficient of the medium [1/mm].
%                Used to compute the per-photon absorption weight
%                w_abs = exp(-mu_a * L).
%
%  n_launched  - total number of photons launched in the pencil-beam
%                simulation (carried for interface symmetry with
%                sia_reflectance; not used in the pathlength weighting
%                itself, which is relative).
%
%  A_detector  - effective detector area [mm^2]. Applied as an overall
%                scaling factor consistent with Eq. (3)/(5) of the Letter.
%
%  rho_edges   - (1 x n_rho_bins+1) rho bin edges [mm]. Must match those
%                used to build p_rho and to run sia_reflectance.
%
%  L_edges     - (1 x n_L_bins+1) pathlength bin edges [mm].
%
%  p_rho       - (n_rho_bins x 1) normalised source-detector distance PDF
%                at the rho bin centres, p_SD(rho_k). Must integrate to 1
%                over drho (i.e. sum(p_rho)*drho = 1).
%
%  drho        - scalar rho bin width [mm].
%
%  OUTPUT
%  ------
%  hist_mass   - (1 x n_L_bins) raw pathlength mass histogram [mm^2 / mm].
%                Not normalised; call normalize_hist_mass_to_pdf to convert
%                to a proper PDF. Pool across repetitions before normalising
%                for smoother curves (see main demo script Section 9).
%
%  NOTES
%  -----
%  - The function returns a zero histogram silently if no photons or no
%    p_rho values are provided, matching the behaviour of sia_reflectance.
%  - The pathlength mass histogram shares the same rho-weighting structure
%    as the reflectance integral, so the two quantities are internally
%    consistent.
%
%  RELATIONSHIP TO THE LETTER
%  --------------------------
%  Implements Eq. (5):
%      P(L) = integral f(L|rho) p_SD(rho) drho
%  discretised as described above, with absorption weighting exp(-mu_a*L)
%  applied to each photon weight as per the convention in the Letter.
%
%  SEE ALSO
%  --------
%  sia_reflectance.m       - reflectance integral (Eq. 3)
%  normalize_hist_mass_to_pdf (in DEMO_SIA_FSD_ALL_GEOMETRIES.m)

% --- Initialise output ---
hist_mass = zeros(1, numel(L_edges) - 1);

% Guard: return zeros silently if inputs are empty.
if isempty(det_data) || isempty(p_rho)
    return;
end

% --- Extract photon properties from the packed data matrix ---
L   = double(det_data(2,:)).';   % pathlength for each photon [mm]
dx  = double(det_data(3,:)).';   % x exit position relative to source [mm]
dy  = double(det_data(4,:)).';   % y exit position relative to source [mm]
rho = hypot(dx, dy);             % radial exit distance [mm]

% Remove photons with non-finite or negative pathlengths / exit distances.
ok  = isfinite(L) & isfinite(rho) & (L >= 0) & (rho >= 0);
L   = L(ok);
rho = rho(ok);

if isempty(L)
    return;
end

% --- Compute absorption weight for each photon ---
% w_abs = exp(-mu_a * L) 
w_abs = exp(-mu_a .* L);

% --- Precompute the SIA ring weights ---
% ring_weight(k) = p_SD(rho_k) * drho / A_ring(rho_k)
% where A_ring(k) = pi * (rho_outer_k^2 - rho_inner_k^2) [mm^2].
% This converts a ring-summed photon weight into a geometry-weighted
% reflectance contribution (same factor as in sia_reflectance).
annulus_areas = pi * (rho_edges(2:end).^2 - rho_edges(1:end-1).^2);  % [mm^2]
annulus_areas = double(annulus_areas(:));
ring_weight   = (double(p_rho(:)) .* double(drho)) ./ max(annulus_areas, eps);

% --- Step 1: Assign each photon to its rho bin ---
bin_rho = discretize(rho, rho_edges);  % integer bin index or NaN

% Discard photons outside the rho grid.
ok_rho = ~isnan(bin_rho);
if ~any(ok_rho)
    return;
end

% --- Step 2: Compute the combined SIA weight for each accepted photon ---
% w_sia = A_detector * w_abs * ring_weight(rho_bin)
% This weight encodes both the absorption attenuation and the geometry
% contribution of the photon's rho bin.
w_sia = A_detector * w_abs(ok_rho) .* ring_weight(bin_rho(ok_rho));

% --- Step 3: Assign each accepted photon to its pathlength bin ---
bin_L = discretize(L(ok_rho), L_edges);  % integer L bin index or NaN

% Discard photons outside the pathlength histogram range.
ok_L = ~isnan(bin_L);
if ~any(ok_L)
    return;
end

% --- Step 4: Accumulate weighted pathlength mass into the histogram ---
% For each L bin, sum the SIA weights of all photons falling in that bin.
hist_mass = accumarray(double(bin_L(ok_L)), ...    % L bin indices
                        double(w_sia(ok_L)), ...    % SIA weights
                        [numel(L_edges)-1, 1], ...  % output size
                        @sum, 0).';                 % row vector output

end % function sia_pathlength
