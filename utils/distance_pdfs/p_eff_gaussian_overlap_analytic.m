function curves = p_eff_gaussian_overlap_analytic(rho_edges, w_mm, R_s, R_d, NrS, NrD, Nr_int, Nk, kmax)
%P_EFF_GAUSSIAN_OVERLAP_ANALYTIC  Non-MC p_eff(rho) curves for overlapping,
% truncated Gaussian source/detector weights.
%
% Output struct fields:
%   curves.rho_mid
%   curves.drho
%   curves.p_base
%   curves.p_hankel
%   curves.p_naive1
%   curves.p_naive2

if any([w_mm, R_s, R_d] <= 0)
    error('w_mm, R_s, R_d must be > 0');
end

rho_edges = rho_edges(:).';
drho      = rho_edges(2) - rho_edges(1);
rho_mid   = ((rho_edges(1:end-1) + rho_edges(2:end))/2).';

wfun = @(r) exp(-2*(r.^2)/(w_mm^2));
Z    = @(R) (pi*w_mm^2/2) * (1 - exp(-2*R.^2/w_mm^2));

Zs = Z(R_s);
Zd = Z(R_d);

qs = @(r) wfun(r) / Zs;
qd = @(r) wfun(r) / Zd;

pr_s = @(r) 2*pi*r.*qs(r);
pr_d = @(r) 2*pi*r.*qd(r);

%% 1) Baseline real-space quadrature
rS = linspace(0, R_s, NrS).';   drS = rS(2)-rS(1);
rD = linspace(0, R_d, NrD).';   drD = rD(2)-rD(1);

pRS = pr_s(rS);  pRS = pRS ./ max(trapz(rS, pRS), eps);
pRD = pr_d(rD);  pRD = pRD ./ max(trapz(rD, pRD), eps);

Wrt = (pRS * pRD.') * (drS * drD);

rS_safe = rS; rS_safe(1) = eps;
rD_safe = rD; rD_safe(1) = eps;
[RR, TT] = ndgrid(rS_safe, rD_safe);
RT = 2*RR.*TT;

p_bin = zeros(numel(rho_edges)-1, 1);

for b = 1:numel(p_bin)
    rho1 = rho_edges(b);
    rho2 = rho_edges(b+1);

    rmin = abs(RR - TT);
    rmax = RR + TT;

    a = max(rho1, rmin);
    c = min(rho2, rmax);

    ok = (c > a) & (RT > 0);
    if ~any(ok(:))
        p_bin(b) = 0;
        continue;
    end

    cos_a = (RR.^2 + TT.^2 - a.^2) ./ RT;
    cos_c = (RR.^2 + TT.^2 - c.^2) ./ RT;

    cos_a = min(max(cos_a, -1), 1);
    cos_c = min(max(cos_c, -1), 1);

    dphi = 2 * (acos(cos_c) - acos(cos_a));
    dphi(~ok) = 0;

    prob = dphi / (2*pi);
    p_bin(b) = sum(sum(prob .* Wrt));
end

p_base = p_bin ./ (diff(rho_edges(:)) + eps);
p_base(~isfinite(p_base)) = 0;
p_base(p_base < 0) = 0;
p_base = p_base ./ (sum(p_base)*drho + eps);

%% 2) Hankel / Fourier-Bessel method
k = linspace(0, kmax, Nk);

rS_h = linspace(0, R_s, Nr_int);
rD_h = linspace(0, R_d, Nr_int);

qs_r = qs(rS_h) .* rS_h;
qd_r = qd(rD_h) .* rD_h;

Q_s = zeros(size(k));
Q_d = zeros(size(k));

for ii = 1:numel(k)
    ki = k(ii);
    Q_s(ii) = 2*pi * trapz(rS_h, qs_r .* besselj(0, ki*rS_h));
    Q_d(ii) = 2*pi * trapz(rD_h, qd_r .* besselj(0, ki*rD_h));
end

p_hankel = zeros(numel(rho_mid),1);
for j = 1:numel(rho_mid)
    r = rho_mid(j);
    integrand = (Q_s .* Q_d) .* besselj(0, k*r) .* k;
    p_hankel(j) = r * trapz(k, integrand);
end

p_hankel(~isfinite(p_hankel)) = 0;
p_hankel(p_hankel < 0) = 0;
p_hankel = p_hankel ./ (sum(p_hankel)*drho + eps);

%% 3) Naive overlays
Rref = min(R_s, R_d);
df   = 2*Rref;

p_unif = p_rho_circle_local(rho_mid, df);
p_unif = p_unif ./ (sum(p_unif)*drho + eps);

g1 = exp(-2*(rho_mid.^2)/(w_mm^2));
g2 = exp(-4*(rho_mid.^2)/(w_mm^2));

p_naive1 = p_unif .* g1;  p_naive1 = p_naive1 ./ (sum(p_naive1)*drho + eps);
p_naive2 = p_unif .* g2;  p_naive2 = p_naive2 ./ (sum(p_naive2)*drho + eps);

curves.rho_mid   = rho_mid;
curves.drho      = drho;
curves.p_base    = p_base;
curves.p_hankel  = p_hankel;
curves.p_naive1  = p_naive1;
curves.p_naive2  = p_naive2;
end

function p = p_rho_circle_local(rho, df)
rho = double(rho(:));
df  = double(df);

p = zeros(size(rho));
in = (rho >= 0) & (rho <= df);
if ~any(in), return; end

x = rho(in)./df;
x = min(max(x, 0), 1);

p(in) = (16.*rho(in))./(pi*df.^2).*acos(x) ...
      - (16.*rho(in).^2)./(pi*df.^3).*sqrt(max(0, 1 - x.^2));

p(~isfinite(p)) = 0;
p(p < 0) = 0;
end
