function p = p_rho_nonoverlap_disks_analytic(rho, sep, rA, rB, Ntheta)
%P_RHO_NONOVERLAP_DISKS_ANALYTIC  Semi-analytic p(rho) for two uniform
% disks separated by sep.

if nargin < 5 || isempty(Ntheta)
    Ntheta = 4096;
end

rho_in = rho;
rho    = double(rho(:));
sep    = double(sep);
rA     = double(rA);
rB     = double(rB);

AreaA = pi*rA^2;
AreaB = pi*rB^2;

theta = linspace(0, 2*pi, Ntheta);
cth   = cos(theta);

p = zeros(size(rho));

rho_min = abs(sep - (rA + rB));
rho_max = sep + (rA + rB);

in = (rho >= rho_min) & (rho <= rho_max);
if ~any(in)
    p = reshape(p, size(rho_in));
    return;
end

for i = find(in).'
    r = rho(i);
    s = sqrt(sep^2 + r^2 + 2*sep*r*cth);
    Aov = overlap_area_two_circles_local(s, rA, rB);
    I = trapz(theta, Aov);
    p(i) = (r / (AreaA*AreaB)) * I;
end

p(~isfinite(p)) = 0;
p(p < 0) = 0;

p = reshape(p, size(rho_in));
end

function A = overlap_area_two_circles_local(s, a, b)
s = double(s(:));
a = double(a);
b = double(b);

A = zeros(size(s));

no = (s >= (a + b));
A(no) = 0;

co = (s <= abs(a - b));
A(co) = pi * min(a,b)^2;

pa = ~(no | co);
sp = s(pa);

c1 = (sp.^2 + a^2 - b^2) ./ (2*sp*a);
c2 = (sp.^2 + b^2 - a^2) ./ (2*sp*b);
c1 = min(max(c1, -1), 1);
c2 = min(max(c2, -1), 1);

term = (-sp + a + b) .* (sp + a - b) .* (sp - a + b) .* (sp + a + b);
term = max(term, 0);

A(pa) = a^2 .* acos(c1) + b^2 .* acos(c2) - 0.5 .* sqrt(term);
end
