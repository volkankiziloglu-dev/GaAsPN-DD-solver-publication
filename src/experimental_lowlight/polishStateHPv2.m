function [out, info] = polishStateHPv2(dev, r, V, Gprofile, kind, opts)
%POLISHSTATEHPV2 EXPERIMENTAL (v2: convergence flag from the polishing's own Newton step; opts.convFlag) mixed-precision Newton polishing (iterative refinement).
% Starts from a converged gummelSolveQFPAnalytic(X) output (r); keeps the state in two parts hi+lo,
% computes the residual with assembleResidualQFPHP(kind) (extended precision), solves the correction with the production
% analytic block Jacobians (double) and adds it to the lo part. Block Gauss-Seidel
% (psi, phin, phip) sweeps; stops when the residual norm no longer decreases.
% opts.maxSweeps (6), opts.nDigits (40). maxSweeps=0 -> HP evaluation only (no polishing).
if nargin < 6, opts = struct(); end
if ~isfield(opts,'maxSweeps'), opts.maxSweeps = 6; end
if ~isfield(opts,'nDigits'), opts.nDigits = 40; end
if ~isfield(opts,'convFlag'), opts.convFlag = 'polish'; end   % 'polish' | 'inherit'
if ~isfield(opts,'stepTol'), opts.stepTol = 1e-12; end
c = constants(); N = dev.N;
S = struct('psi',r.psi,'psiLo',zeros(N,1),'phin',r.phin,'phinLo',zeros(N,1),'phip',r.phip,'phipLo',zeros(N,1));
res = @(SS) assembleResidualQFPHP(dev, SS, V, Gprofile, kind, opts.nDigits);
F = res(S);
nrm = @(FF) [norm(FF(2:3:end)), norm(FF(3:3:end))];
hist = nrm(F); best = S; bestN = sum(hist);
names = {'psi','phin','phip'};
stepHist = [];
for sweep = 1:opts.maxSweeps
  sweepStep = 0;
  for blk = 1:3
    Sd = S; % Jacobian at rounded state
    psiC = S.psi + S.psiLo; phinC = S.phin + S.phinLo; phipC = S.phip + S.phipLo;
    switch blk
      case 1, J = assemblePsiJacobianAnalyticQFP(dev, psiC, phinC, phipC);
      case 2, J = assemblePhinJacobianAnalyticQFP(dev, psiC, phinC, phipC, V, zeros(N,1));
      case 3, J = assemblePhipJacobianAnalyticQFP(dev, psiC, phinC, phipC, V, zeros(N,1));
    end
    Fb = F(blk:3:end);
    rowNorm = full(max(abs(J),[],2)); good = isfinite(rowNorm) & rowNorm > 0;
    floorNorm = max(min(rowNorm(good))*1e-12, realmin); sc = 1./max(rowNorm, floorNorm);
    delta = (spdiags(sc,0,N,N)*J) \ (-(sc.*Fb));
    if ~all(isfinite(delta)), sweepStep = Inf; continue; end
    sweepStep = max(sweepStep, max(abs(delta)));
    hi = Sd.(names{blk}); lo = Sd.([names{blk} 'Lo']) + delta;
    t = hi + lo; lo = lo - (t - hi); hi = t;          % renormalize (fast two-sum)
    S.(names{blk}) = hi; S.([names{blk} 'Lo']) = lo;
    F = res(S);
  end
  hn = nrm(F); hist(end+1,:) = hn; stepHist(end+1) = sweepStep; %#ok<AGROW>
  if all(isfinite(hn)) && sum(hn) < bestN
    best = S; bestN = sum(hn);
  elseif sweep > 1
    break;
  end
end
S = best;
[F, n, p, JnSeg, JpSeg] = res(S);
psi = S.psi + S.psiLo; phin = S.phin + S.phinLo; phip = S.phip + S.phipLo;
JtotSeg = ifaceCorrectTotalCurrent(dev, JnSeg + JpSeg, n, p);
idxP = 2:N-1;
chargeScale = c.Q * max([max(dev.Nd(idxP)), max(dev.Na(idxP)), max(n(idxP)), max(p(idxP)), 1e10]);
out = r;
out.psi = psi; out.phin = phin; out.phip = phip; out.n = n; out.p = p;
out.JnSeg = JnSeg; out.JpSeg = JpSeg; out.Jseg = JtotSeg; out.J = median(JtotSeg);
out.isFinite = r.isFinite && all(isfinite(F)) && all(isfinite(n)) && all(isfinite(p));
out.poissonResidual = max(abs(F(3*(idxP-1)+1))) / chargeScale;
out.poissonConverged = out.poissonResidual < 1e-6;
out.currentConservation = checkCurrentConservation(JtotSeg, out.J, 1e-7, 1e-3);
out.continuityResiduals = checkContinuityResiduals(dev, Gprofile, n, p, JnSeg, JpSeg, psi);
polishConverged = ~isempty(stepHist) && min(stepHist) < opts.stepTol;
if strcmp(opts.convFlag,'inherit'), out.numericalConverged = r.numericalConverged && out.isFinite;
else, out.numericalConverged = polishConverged && out.isFinite; end
info = struct('normHistory',hist,'stepHistory',stepHist,'polishConverged',polishConverged,'inheritedConverged',r.numericalConverged,'sweeps',size(hist,1)-1, ...
  'lo',struct('psi',S.psiLo,'phin',S.phinLo,'phip',S.phipLo));
end
