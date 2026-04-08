function p = p_rho_square_square_analytic(rho, s_src, s_det, dx, dy, Ntheta)
%P_RHO_SQUARE_SQUARE_ANALYTIC  Semi-analytic p(rho) via covariogram + 1D θ integral.

if nargin < 6 || isempty(Ntheta)
    Ntheta = 4096;
end

rho_in = rho;
rho    = double(rho(:));

hs = double(s_src)/2;
hd = double(s_det)/2;

As = double(s_src)^2;
Ad = double(s_det)^2;

theta = linspace(0, 2*pi, Ntheta);
ct = cos(theta);
st = sin(theta);

tx = double(dx) + rho * ct;
ty = double(dy) + rho * st;

ovx = max(0, (hs + hd) - abs(tx));
ovy = max(0, (hs + hd) - abs(ty));

Acap = ovx .* ovy;
Itheta = trapz(theta, Acap, 2);

p = (rho ./ (As*Ad)) .* Itheta;

p(~isfinite(p)) = 0;
p(p < 0) = 0;

p = reshape(p, size(rho_in));
end
