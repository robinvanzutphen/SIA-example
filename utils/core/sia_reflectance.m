function R = sia_reflectance(rho_all, w_abs_all, rho_edges, ring_weight, A_detector, n_launched)
% SIA_REFLECTANCE  Estimate reflectance for a finite source-detector geometry
%                  from a single pencil-beam Monte Carlo simulation.
%
%  R = sia_reflectance(rho_all, w_abs_all, rho_edges, ring_weight, ...
%                      A_detector, n_launched)
%
%  This function implements the discretised form of the single-integral
%  approximation (SIA) described in Eq. (3) of the Letter:
%
%      R = A_d * integral_0^rho_max  R(rho) * p_SD(rho) d_rho
%
%  where R(rho) is the pencil-beam radial reflectance per unit area [1/mm^2]
%  and p_SD(rho) is the source-detector distance distribution.
%
%  The integral is evaluated as a Riemann sum over a uniform rho grid.
%  For each rho bin k:
%
%      R(rho_k) ~ (1/n_launched) * sum_{photons in ring k} w_abs_i
%                 / A_ring(rho_k)
%
%  where A_ring = pi*(rho_outer^2 - rho_inner^2) is the area of annular bin k.
%  Substituting into the Riemann sum and collecting terms:
%
%      R ~= (A_d / n_launched) * sum_k ring_sums(k) * ring_weight(k)
%
%  with ring_weight(k) = p_SD(rho_k) * drho / A_ring(rho_k)
%  precomputed externally (see the main demo script).
%
%  INPUTS
%  ------
%  rho_all       - (N x 1) radial exit distances of detected pencil-beam
%                  photons, measured from the source position [mm].
%                  Only photons passing the NA acceptance cut should be
%                  included (pre-filtered in the main script).
%
%  w_abs_all     - (N x 1) absorption weights for each detected photon,
%                  w_abs = exp(-mu_a * L), where L is the photon pathlength
%                  [mm] and mu_a is the absorption coefficient [1/mm].
%
%  rho_edges     - (1 x n_bins+1) bin edges of the uniform rho grid [mm].
%                  Must span [0, rho_max] with n_bins equally spaced bins.
%
%  ring_weight   - (n_bins x 1) combined weight for each rho bin:
%                    ring_weight(k) = p_SD(rho_k) * drho / A_ring(rho_k)
%                  where p_SD is the normalised source-detector distance PDF,
%                  drho is the bin width [mm], and A_ring is the bin area [mm^2].
%                  Precomputed in the main loop for efficiency.
%
%  A_detector    - scalar effective detector area [mm^2]. For uniform
%                  detection this is the geometric aperture area; for
%                  non-uniform detection it is A_det,eff = integral w_d(r) dA
%                  (Letter Sec. 2C).
%
%  n_launched    - scalar total number of photons launched in the pencil-beam
%                  simulation. Used to normalise the ring sums to a per-photon
%                  reflectance.
%
%  OUTPUT
%  ------
%  R             - scalar: estimated reflectance (collected power fraction,
%                  dimensionless). Returns 0 if no photons are available.
%
%  RELATIONSHIP TO THE LETTER
%  --------------------------
%  This implements Eq. (3):
%      P_R / P_0 = A_d * integral R(rho) p_SD(rho) drho
%
%  with the pencil-beam R(rho) estimated from the annular photon sums and
%  p_SD(rho) provided analytically or semi-analytically (see companion
%  functions p_rho_*_analytic.m).
%


% NOTE: We use the functions discretize and accumarray to build R(rho) in a
% vectorized manner. A loop over rho would be the simpler, but costlier
% alternative


% --- Guard: return 0 if no photons are available ---
if isempty(rho_all) || isempty(w_abs_all) || n_launched <= 0
    R = 0;
    return;
end

% Ensure column vectors for consistent indexing.
rho_all   = double(rho_all(:));
w_abs_all = double(w_abs_all(:));

% --- Step 1: Assign each photon to its rho bin ---
% discretize returns NaN for photons outside the rho_edges range.
bin_rho = discretize(rho_all, rho_edges);  % integer bin index or NaN

% Retain only photons that fall within the defined rho grid.
ok_rho = ~isnan(bin_rho);

if ~any(ok_rho)
    % No photons fall within the rho support of this geometry.
    R = 0;
    return;
end

n_bins = numel(rho_edges) - 1;  % total number of rho bins

% --- Step 2: Sum absorption weights within each rho bin ---
% ring_sums(k) = sum of w_abs for all photons exiting in annular bin k.
% Dividing by n_launched gives the (unnormalised) pencil-beam reflectance
% per unit area for that ring (before the area normalisation in ring_weight).
ring_sums = accumarray(double(bin_rho(ok_rho)), ...   % bin indices
                        double(w_abs_all(ok_rho)), ... % weights
                        [n_bins, 1], ...               % output size
                        @sum, 0);                      % sum within each bin, default 0

% --- Step 3: Evaluate the SIA Riemann sum ---
% Multiply the ring sums by the precomputed ring_weight (which contains
% p_SD(rho_k)*drho/A_ring) and sum over all bins. Scaling by A_detector
% and dividing by n_launched completes the Riemann approximation of Eq. (3).
R = (A_detector / n_launched) * sum(double(ring_sums(:)) .* double(ring_weight(:)));

% Guard against non-finite results (e.g. from empty simulations).
if ~isfinite(R)
    R = 0;
end

end % function sia_reflectance
