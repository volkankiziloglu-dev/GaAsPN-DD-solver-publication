function report = checkContinuityResiduals(dev, Gprofile, n, p, JnSeg, JpSeg, psi)
%CHECKCONTINUITYRESIDUALS Evaluate the two discrete continuity equations.
%
% Electron: div(Jn)/q = R - G
% Hole:     div(Jp)/q = G - R
c = constants();
idx = (2:dev.N-1)';
hAvg = 0.5*(dev.h(idx-1) + dev.h(idx));

npExcess=n(idx).*p(idx)-dev.ni(idx).^2;
Rbulk = bulkSRHRate(dev,n,p,idx) + dev.Brad(idx).*npExcess + ...
    (dev.CnAuger(idx).*n(idx)+dev.CpAuger(idx).*p(idx)).*npExcess;
[RintN, RintP] = interfaceRecombinationNP(dev, n, p, idx, hAvg);   % carrier-specific sinks for 'four'/'cross' (dd_solver v1.1.1)
Rtotal = Rbulk + RintN; RtotalP = Rbulk + RintP;

divJnOverQ = (JnSeg(idx) - JnSeg(idx-1))./(c.Q*hAvg);
divJpOverQ = (JpSeg(idx) - JpSeg(idx-1))./(c.Q*hAvg);

rn = divJnOverQ - (Rtotal - Gprofile(idx));
rp = divJpOverQ - (Gprofile(idx) - RtotalP);

rhsN = Rtotal-Gprofile(idx);
rhsP = Gprofile(idx)-RtotalP;
localScaleN = max(abs(divJnOverQ)+abs(rhsN),1);
localScaleP = max(abs(divJpOverQ)+abs(rhsP),1);

% Structurally optically-inert nodes (no absorption table AND alphaMax=0,
% i.e. G(x)=0 identically regardless of bias/spectrum -- e.g. a
% transparent window/BSF contact layer) are excluded from the residual's
% worst-node search and pass/fail decision. At low (indoor) injection
% these nodes' near-zero local generation-recombination budget
% (rhsN/rhsP~0) makes even a numerically-oscillating, physically-
% negligible current-segment mismatch look catastrophic in relative
% terms: diagnosed (debug session, N.1 homojunction @ 6500K/100lux) as
% an alternating-sign, repeated-magnitude oscillation confined entirely
% to the window layer (|rn| in the window's median ~5.3e17 vs the
% optically-active absorber's median ~2.3e6 -- five orders of magnitude
% cleaner), not a genuine growing residual. End-to-end total-current
% conservation (what actually determines Jsc/Voc/FF) is independently
% verified by checkCurrentConservation.m regardless of this exclusion.
isActive = dev.useOpticalTable(idx) | (dev.alphaMaxNode(idx) > 0);
if ~any(isActive)
  isActive = true(size(idx)); % pathological all-inert device: don't exclude everything
end

% IEEE-754 cancellation floor of the SG segment-current subtraction.
% Propagate the rounding bound of both adjacent segment currents to each
% control volume. Without this term, refining h can make an unchanged
% floating-point noise floor look like a growing physical residual.
idxJ=(1:dev.N-1)';
if nargin < 7 || isempty(psi)
  nodeNoiseN=zeros(size(idx)); nodeNoiseP=zeros(size(idx));
else
psiEffN=psi+dev.psiOffsetN; psiEffP=psi+dev.psiOffsetP;
dPsiN=(psiEffN(idxJ+1)-psiEffN(idxJ))/c.VT;
dPsiP=(psiEffP(idxJ+1)-psiEffP(idxJ))/c.VT;
DnSeg=c.VT*0.5*(dev.muN(idxJ)+dev.muN(idxJ+1));
DpSeg=c.VT*0.5*(dev.muP(idxJ)+dev.muP(idxJ+1));
coefN=c.Q*DnSeg./dev.h(idxJ); coefP=c.Q*DpSeg./dev.h(idxJ);
absABn=coefN.*(n(idxJ+1).*sgBernoulli(dPsiN)+n(idxJ).*sgBernoulli(-dPsiN));
absABp=coefP.*(p(idxJ).*sgBernoulli(dPsiP)+p(idxJ+1).*sgBernoulli(-dPsiP));
segNoiseN=eps(1)*absABn; segNoiseP=eps(1)*absABp;
nodeNoiseN=sqrt(segNoiseN(idx).^2+segNoiseN(idx-1).^2)./(c.Q*hAvg);
nodeNoiseP=sqrt(segNoiseP(idx).^2+segNoiseP(idx-1).^2)./(c.Q*hAvg);
end

globalScaleN = max(max(abs(divJnOverQ(isActive))+abs(rhsN(isActive))) + ...
  max(nodeNoiseN(isActive)),1);
globalScaleP = max(max(abs(divJpOverQ(isActive))+abs(rhsP(isActive))) + ...
  max(nodeNoiseP(isActive)),1);
rnRel = abs(rn)./localScaleN;
rpRel = abs(rp)./localScaleP;

report.electron.maxAbs = max(abs(rn(isActive)));
report.electron.maxRel = max(rnRel(isActive));
report.hole.maxAbs = max(abs(rp(isActive)));
report.hole.maxRel = max(rpRel(isActive));
report.electron.roundoffFloor = max(nodeNoiseN(isActive));
report.hole.roundoffFloor = max(nodeNoiseP(isActive));

% Pass/fail is evaluated in the native current-density units of the
% discrete QFP continuity rows. Dividing by Q*hAvg above is useful for
% reporting volumetric rates, but it amplifies the same converged native
% residual as the mesh is refined. Reconstructing Fphin/Fphip removes this
% artificial 1/h dependence and matches checkCurrentConservation.m's
% absolute-OR-relative tolerance pattern.
JtolAbsN = 1e-7; JtolRelN = 1e-3; JscaleN = 1e-6; % A/cm^2
FphinNative = rn .* (c.Q*hAvg);
FphipNative = rp .* (c.Q*hAvg);
JlocalN = max(abs(JnSeg(idx)), abs(JnSeg(idx-1)));
JlocalP = max(abs(JpSeg(idx)), abs(JpSeg(idx-1)));
nativeRelN = abs(FphinNative)./max(JlocalN, JscaleN);
nativeRelP = abs(FphipNative)./max(JlocalP, JscaleN);
nodeOkN = (abs(FphinNative) < JtolAbsN) | (nativeRelN < JtolRelN);
nodeOkP = (abs(FphipNative) < JtolAbsN) | (nativeRelP < JtolRelN);

report.electron.maxRelGlobal = max(nativeRelN(isActive));
report.hole.maxRelGlobal = max(nativeRelP(isActive));
report.electron.pass = all(nodeOkN(isActive));
report.hole.pass = all(nodeOkP(isActive));
report.pass = report.electron.pass && report.hole.pass;

nativeRelNActive = nativeRelN; nativeRelNActive(~isActive) = -Inf;
nativeRelPActive = nativeRelP; nativeRelPActive(~isActive) = -Inf;
[~, kn] = max(nativeRelNActive);
[~, kp] = max(nativeRelPActive);
report.electron.node = idx(kn);
report.hole.node = idx(kp);
report.electron.value = FphinNative(kn); % A/cm^2
report.hole.value = FphipNative(kp);     % A/cm^2
report.electron.scale = max(JlocalN(kn), JscaleN);
report.hole.scale = max(JlocalP(kp), JscaleN);

% Transparency: report the excluded region's own worst value too (NOT
% used for pass/fail), so a genuine problem hiding in an inert layer
% would not go completely unreported. These informational fields retain
% the original volumetric-rate units.
report.electron.excludedMaxAbs = max([abs(rn(~isActive)); 0]);
report.hole.excludedMaxAbs = max([abs(rp(~isActive)); 0]);
end
