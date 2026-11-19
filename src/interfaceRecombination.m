function R_int_perNode = interfaceRecombination(dev, nPrev, pPrev, idx, hAvg, carrier)
% dd_solver v1.1.0: optional 6th argument carrier ('n'/'p') and opt-in models dev.ifaceModel = 'onesided' | 'four' | 'cross'
% (set by dd_enableInterfaceModel.m; ported from Sb2S3_Sb2Se3_Paper). Without dev.ifaceModel the production path is unchanged.
% Interface (2D, areal) recombination at the ETL/Abs heterojunction, added
% on top of bulk SRH (D.1) and radiative (D.2) -- see MATHEMATICAL_MODEL.md
% Section J. R_int [cm^-2 s^-1], midgap-trap SRH form:
%   R_int = (n_int*p_int - ni_int^2) / [(n_int+ni_int)/Sp_int + (p_int+ni_int)/Sn_int]
% n_int, p_int: averaged from the two nodes straddling the junction (frozen
% at the previous Gummel iterate, fully explicit -- consistent with how
% bulk SRH's "other carrier" is already frozen). ni_int: geometric mean of
% the two sides' ni (simple, symmetric choice for an assumed-midgap trap;
% NOT tied to a specific interface trap energy Et).
%
% R_int is converted to a volumetric-equivalent sink and split (half each)
% onto the two nodes immediately adjacent to the junction (dev.N_ETL,
% dev.N_ETL+1), each normalized by ITS OWN control-volume hAvg. This IS the
% exactly-conservative split: since each node's contribution to the
% discretized system ultimately gets weighted by that SAME node's own
% hAvg elsewhere in the equation, integrating [hAvg(posE)*R_int_perNode(posE)]
% + [hAvg(posA)*R_int_perNode(posA)] = 0.5*Rint + 0.5*Rint = Rint EXACTLY,
% for any hAvg(posE), hAvg(posA) -- regardless of how differently the two
% regions happen to be meshed.
%
% (Prior history, corrected by external review: an earlier version divided
% by the single shared interface spacing h(N_ETL) instead, believing the
% per-node-hAvg split was mesh-dependent/the cause of an observed
% instability. That diagnosis was wrong -- the h(N_ETL) split is NOT
% exactly conservative (the two nodes' hAvg generally differ substantially
% from h(N_ETL) itself, so the integrated total deviates from R_int by a
% mesh-dependent factor), while the per-node-hAvg split IS exactly
% conservative by construction, as shown above. The real cause of the
% earlier instability was a sign error in solveHoleContinuity.m's use of
% this function's output (see that file's note) -- fixed independently,
% and re-confirmed here by reverting to the conservative per-node-hAvg
% split and re-running test_interface_recomb_sanity.m.)
%
% Returns a zero vector (sized like idx) when Sn_int=0 or Sp_int=0
% (default), exactly recovering the "ideal interface" behavior used
% throughout this model until now.
%
% MULTI-JUNCTION GENERALIZATION (RearHJT/buildDeviceML.m): if dev carries
% dev.junctionNodes (array of junction locations) and dev.Sn_int_j/Sp_int_j
% (per-junction velocities), this loops over ALL of them -- each junction
% contributes independently and additively, same formula/split as below.
% For a 2-region device built via the original common/buildDevice.m
% (dev.N_ETL/Sn_int/Sp_int, no junctionNodes field), falls back to the
% original single-junction path, UNCHANGED, for exact backward
% compatibility.
R_int_perNode = zeros(numel(idx), 1); % column, matching h/ni/Gprofile(idx) shape (NOT size(idx), which is a row)
ni = dev.ni;
if nargin < 6, carrier = ''; end

if isfield(dev,'junctionNodes') && ~isempty(dev.junctionNodes)
  for jk = 1:numel(dev.junctionNodes)
    Sn = dev.Sn_int_j(jk); Sp = dev.Sp_int_j(jk);
    if Sn <= 0 || Sp <= 0, continue; end
    iE = dev.junctionNodes(jk); iA = iE + 1;
    model = ''; if isfield(dev,'ifaceModel'), model = dev.ifaceModel; end
    switch model
      case 'onesided'
        R_int_perNode = addJunctionContributionOneSided(R_int_perNode, nPrev, pPrev, dev.ifaceN0A(jk), dev.ifaceP0A(jk), dev.ifaceNodeA(jk), Sn, Sp, idx, hAvg, ifaceEtOf(dev));
      case 'four'
        R_int_perNode = addJunctionContributionFour(R_int_perNode, nPrev, pPrev, dev.ifaceN0(jk), dev.ifaceP0(jk), iE, iA, Sn, Sp, idx, hAvg, carrier, ifaceEtOf(dev));
      case 'cross'
        R_int_perNode = addJunctionContributionCross(R_int_perNode, nPrev, pPrev, dev.ifaceNP0(jk), iE, iA, Sn, Sp, idx, hAvg, carrier);
      otherwise
        [wE, wA] = intRecWeights(dev, iE, iA);
        R_int_perNode = addJunctionContribution(R_int_perNode, nPrev, pPrev, ni, iE, iA, Sn, Sp, idx, hAvg, wE, wA);
    end
  end
  return;
end

if ~isfield(dev,'Sn_int') || ~isfield(dev,'Sp_int') || dev.Sn_int <= 0 || dev.Sp_int <= 0
  return;
end
iE = dev.N_ETL; iA = dev.N_ETL + 1;
[wE, wA] = intRecWeights(dev, iE, iA);
R_int_perNode = addJunctionContribution(R_int_perNode, nPrev, pPrev, ni, iE, iA, dev.Sn_int, dev.Sp_int, idx, hAvg, wE, wA);
end

function R_int_perNode = addJunctionContribution(R_int_perNode, nPrev, pPrev, ni, iE, iA, Sn, Sp, idx, hAvg, wE, wA)
% wE = wA = 0.5: geometric mean (production); see intRecWeights.m for the lower-gap-side alternative.
% N.19 fix: n_int/p_int MUST use the geometric mean, matching ni_int's
% own sqrt(ni(iE)*ni(iA)) below -- an arithmetic mean here is WRONG and
% was the actual root cause of N.16.16's "interface-recombination
% stagnation" (previously misdiagnosed as a Gummel block-decomposition
% limitation). Proof: at a real p-n junction n(iE) and p(iA) differ from
% n(iA) and p(iE) by many orders of magnitude (majority vs minority on
% each side). Even when BOTH nodes are individually at their own exact
% local equilibrium (n(iE)*p(iE)=ni(iE)^2 AND n(iA)*p(iA)=ni(iA)^2
% exactly), the arithmetic-mean product 0.5*(n(iE)+n(iA)) *
% 0.5*(p(iE)+p(iA)) is dominated by the CROSS terms n(iE)*p(iA) and
% n(iA)*p(iE) (majority*majority), NOT n(iE)*p(iE)/n(iA)*p(iA) -- giving
% a spurious nonzero R_int even at true equilibrium (verified directly:
% 39x violation on a minimal test device, vs 1.000000000000003 with the
% geometric mean). The geometric mean is exact here specifically because
% sqrt(n(iE)*p(iE)) * sqrt(n(iA)*p(iA)) = sqrt(ni(iE)^2) * sqrt(ni(iA)^2)
% = ni(iE)*ni(iA) = ni_int^2 whenever each node individually satisfies
% mass action -- an algebraic identity, not an approximation. See
% MATHEMATICAL_MODEL.md N.19.
if wE == 0.5   % production path, bit-identical to the original geometric mean
  n_int = sqrt(nPrev(iE) * nPrev(iA));
  p_int = sqrt(pPrev(iE) * pPrev(iA));
  ni_int = sqrt(ni(iE) * ni(iA));
else
  n_int = nPrev(iE)^wE * nPrev(iA)^wA;
  p_int = pPrev(iE)^wE * pPrev(iA)^wA;
  ni_int = ni(iE)^wE * ni(iA)^wA;
end
denom_int = (n_int + ni_int)/Sp + (p_int + ni_int)/Sn;
Rint = (n_int*p_int - ni_int^2) / denom_int;

posE = iE - idx(1) + 1; posA = iA - idx(1) + 1; % position within idx (1-indexed)
R_int_perNode(posE) = R_int_perNode(posE) + 0.5*Rint/hAvg(posE);
R_int_perNode(posA) = R_int_perNode(posA) + 0.5*Rint/hAvg(posA);
end

function Et = ifaceEtOf(dev)
% Optional interface trap level relative to the effective midgap (eV); 0 = midgap.
if isfield(dev,'ifaceEt') && ~isempty(dev.ifaceEt), Et = dev.ifaceEt; else, Et = 0; end
end

function R_int_perNode = addJunctionContributionOneSided(R_int_perNode, nPrev, pPrev, n0A, p0A, k, Sn, Sp, idx, hAvg, Et)
% One-sided planar defect: the trap sees only the densities of node k (absorber side, chosen in dd_enableInterfaceModel):
%   R_int = (nA*pA - n0A*p0A) / [(nA + n1)/Sp + (pA + p1)/Sn],  n1 = sqrt(n0A*p0A)*exp(Et/VT), p1 = sqrt(n0A*p0A)*exp(-Et/VT),
% exactly zero in equilibrium; the whole sink of both carriers is placed on node k.
nA = nPrev(k); pA = pPrev(k); [n1, p1] = ifaceN1P1(n0A, p0A, Et);
Rint = (nA*pA - n0A*p0A) / ((nA + n1)/Sp + (pA + p1)/Sn);
pos = k - idx(1) + 1;
R_int_perNode(pos) = R_int_perNode(pos) + Rint/hAvg(pos);
end

function R_int_perNode = addJunctionContributionCross(R_int_perNode, nPrev, pPrev, np0, iE, iA, Sn, Sp, idx, hAvg, carrier)
% Cross channel: electrons of the E side (n-type transport layer) with holes of the A side (absorber),
%   R_int = (nE*pA - np0) / [(nE + n1)/Sp + (pA + n1)/Sn], n1 = sqrt(np0); electron sink on iE, hole sink on iA
% (half/half without a carrier argument, legacy callers).
nE = nPrev(iE); pA = pPrev(iA); n1 = sqrt(np0);
Rint = (nE*pA - np0) / ((nE + n1)/Sp + (pA + n1)/Sn);
posE = iE - idx(1) + 1; posA = iA - idx(1) + 1;
if isempty(carrier)
  R_int_perNode(posE) = R_int_perNode(posE) + 0.5*Rint/hAvg(posE);
  R_int_perNode(posA) = R_int_perNode(posA) + 0.5*Rint/hAvg(posA);
elseif carrier == 'n'
  R_int_perNode(posE) = R_int_perNode(posE) + Rint/hAvg(posE);
else
  R_int_perNode(posA) = R_int_perNode(posA) + Rint/hAvg(posA);
end
end

function R_int_perNode = addJunctionContributionFour(R_int_perNode, nPrev, pPrev, N0, P0, iE, iA, Sn, Sp, idx, hAvg, carrier, Et)
% Four-channel interface SRH: Nn = nE + nA, Np = pE + pA, R_int = (Nn*Np - N0*P0) / [(Nn + n1)/Sp + (Np + p1)/Sn];
% each carrier sink is split between the nodes by that node's share of the carrier density (electron split without carrier).
nE = nPrev(iE); nA = nPrev(iA); pE = pPrev(iE); pA = pPrev(iA);
Nn = nE + nA; Np = pE + pA; [n1, p1] = ifaceN1P1(N0, P0, Et);
Rint = (Nn*Np - N0*P0) / ((Nn + n1)/Sp + (Np + p1)/Sn);
posE = iE - idx(1) + 1; posA = iA - idx(1) + 1;
if isempty(carrier) || carrier == 'n', wE = nE/Nn; wA = nA/Nn; else, wE = pE/Np; wA = pA/Np; end
R_int_perNode(posE) = R_int_perNode(posE) + wE*Rint/hAvg(posE);
R_int_perNode(posA) = R_int_perNode(posA) + wA*Rint/hAvg(posA);
end
