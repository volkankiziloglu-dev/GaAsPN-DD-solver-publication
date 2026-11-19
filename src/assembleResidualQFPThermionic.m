function [F, n, p, JnSeg, JpSeg] = assembleResidualQFPThermionic(dev, psi, phin, phip, V, Gprofile)
%ASSEMBLERESIDUALQFPTHERMIONIC Same as assembleResidualQFP.m (NOT
% modified, kept as a separate working reference), except that the
% Scharfetter-Gummel segment flux is REPLACED by thermionicFluxQFP.m's
% saturating-weight flux at any segment(s) flagged via
% dev.thermionicJunction (a junctionNodes(k) value, scalar or vector for
% multiple junctions) / dev.vn_thermionic / dev.vp_thermionic (scalar or
% per-junction vector, matching dev.thermionicJunction's length).
%
% dev.thermionicJunction empty/absent -> behaves IDENTICALLY to
% assembleResidualQFP.m (see tests/test_thermionicJunctionQFP_regression.m).
%
% See docs/THERMIONIC_INTERFACE_PLAN.md (Faz 2) and thermionicFluxQFP.m
% for the physics/derivation. This file is a Newton-solver RESIDUAL, not
% a diagnostic -- Gprofile/interface-recombination/boundary handling are
% UNCHANGED from assembleResidualQFP.m; only the flux LAW at the flagged
% segment(s) differs.
c = constants();
N = dev.N;
h = dev.h;

[n, p] = carrierDensitiesQFP(dev, psi, phin, phip);

% ---- Poisson (interior nodes) ---- (identical to assembleResidualQFP.m)
idxP = (2:N-1)';
epsL = c.EPS0*0.5*(dev.eps(idxP)+dev.eps(idxP-1));
epsR = c.EPS0*0.5*(dev.eps(idxP)+dev.eps(idxP+1));
hL = h(idxP-1); hR = h(idxP); hAvgP = 0.5*(hL+hR);
rho = c.Q*(p(idxP)-n(idxP)+dev.Nd(idxP)-dev.Na(idxP));

Fpsi = zeros(N,1);
Fpsi(idxP) = (epsR.*(psi(idxP+1)-psi(idxP))./hR - epsL.*(psi(idxP)-psi(idxP-1))./hL)./hAvgP + rho;
Fpsi(1) = psi(1) - dev.psiEq(1);
Fpsi(N) = psi(N) - (dev.psiEq(N) + V);

% ---- Segment (Scharfetter-Gummel) currents, i=1..N-1 ----
idxJ = (1:N-1)';
psiEffN = psi + dev.psiOffsetN;
psiEffP = psi + dev.psiOffsetP;
dPsiN = (psiEffN(idxJ+1)-psiEffN(idxJ))/c.VT;
dPsiP = (psiEffP(idxJ+1)-psiEffP(idxJ))/c.VT;
DnSeg = c.VT*0.5*(dev.muN(idxJ)+dev.muN(idxJ+1));
DpSeg = c.VT*0.5*(dev.muP(idxJ)+dev.muP(idxJ+1));
JnSeg = c.Q*DnSeg./h(idxJ) .* (n(idxJ+1).*sgBernoulli(dPsiN) - n(idxJ).*sgBernoulli(-dPsiN));
JpSeg = c.Q*DpSeg./h(idxJ) .* (p(idxJ).*sgBernoulli(dPsiP) - p(idxJ+1).*sgBernoulli(-dPsiP));

% ---- Thermionic-limited override at flagged junction segment(s) ----
if isfield(dev,'thermionicJunction') && ~isempty(dev.thermionicJunction)
  jNodes = dev.thermionicJunction(:)';
  vnAll = dev.vn_thermionic; if isscalar(vnAll), vnAll = repmat(vnAll, size(jNodes)); end
  vpAll = dev.vp_thermionic; if isscalar(vpAll), vpAll = repmat(vpAll, size(jNodes)); end
  for kj = 1:numel(jNodes)
    jSeg = jNodes(kj); % segment jSeg: node jSeg -> jSeg+1
    if jSeg < 1 || jSeg > N-1, continue; end
    [JnSeg(jSeg), JpSeg(jSeg)] = thermionicFluxQFP(c, vnAll(kj), vpAll(kj), ...
      dPsiN(jSeg), dPsiP(jSeg), n(jSeg), n(jSeg+1), p(jSeg), p(jSeg+1));
  end
end

% ---- Continuity (interior nodes) ---- (identical to assembleResidualQFP.m)
idxC = (2:N-1)';
hLc = h(idxC-1); hRc = h(idxC); hAvgC = 0.5*(hLc+hRc);
npExcess = n(idxC).*p(idxC)-dev.ni(idxC).^2;
Rbulk = bulkSRHRate(dev, n, p, idxC) + dev.Brad(idxC).*npExcess + ...
  (dev.CnAuger(idxC).*n(idxC)+dev.CpAuger(idxC).*p(idxC)).*npExcess;   % dd_solver v1.1.1: Auger term was missing in this variant
assertIfaceModelSupported(dev, mfilename);   % 'four'/'cross' only in assembleResidualQFP (dd_solver v1.1.1)
Rint = interfaceRecombination(dev, n, p, idxC, hAvgC);
Rtot = Rbulk + Rint;

Fphin = zeros(N,1);
Fphip = zeros(N,1);
Fphin(idxC) = (JnSeg(idxC) - JnSeg(idxC-1)) - c.Q*hAvgC.*(Rtot - Gprofile(idxC));
Fphip(idxC) = (JpSeg(idxC) - JpSeg(idxC-1)) + c.Q*hAvgC.*(Rtot - Gprofile(idxC));

% ---- Front boundary (i=1): always Dirichlet phi=0 ----
Fphin(1) = phin(1);
Fphip(1) = phip(1);

% ---- Back boundary (i=N) ----
Sback = dev.Sback;
if dev.Na(N) > 0, pRight = dev.Na(N); else, pRight = dev.ni(N)^2/dev.Nd(N); end
nRight = dev.ni(N)^2/pRight;

if dev.Nd(N) > 0
  Fphin(N) = phin(N) - V;
else
  Fphin(N) = JnSeg(N-1) + c.Q*Sback*(n(N)-nRight);
end
if dev.Na(N) > 0
  Fphip(N) = phip(N) - V;
else
  Fphip(N) = JpSeg(N-1) - c.Q*Sback*(p(N)-pRight);
end

F = zeros(3*N,1);
F(1:3:end) = Fpsi;
F(2:3:end) = Fphin;
F(3:3:end) = Fphip;
end
