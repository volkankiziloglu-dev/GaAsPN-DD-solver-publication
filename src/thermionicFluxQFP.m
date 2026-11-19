function [JnSeg, JpSeg] = thermionicFluxQFP(c, vn, vp, dPsiN, dPsiP, nL, nR, pL, pR)
%THERMIONICFLUXQFP Thermionic-emission-limited segment current, replacing
% the standard Scharfetter-Gummel flux at a single designated abrupt
% junction segment.
%
% Ported from the already-validated (this project's own prior work,
% GaAsNP_RHJ_final/solveElectronContinuity_thermionic.m /
% solveHoleContinuity_thermionic.m, see docs/THERMIONIC_INTERFACE_PLAN.md
% Faz 0) v2 saturating-weight formula. Unified substitution rule (derived
% directly from assembleResidualQFP.m's own JnSeg/JpSeg formulas, verified
% to reproduce the old code's per-carrier role assignment exactly):
% wherever sgBernoulli(x) appears in the standard SG segment-current
% formula, replace it with thermionicWeight(x) = exp(min(0,-x)) -- same
% substitution for BOTH carriers, no special-casing needed.
%
% Detailed balance (hand-verified, see plan doc): at equilibrium
% (n(R)/n(L) = exp(dPsiN) exactly, from carrierDensitiesQFP.m's own
% definitions), JnSeg reduces to Q*vn*n(L)*(exp(dPsiN)*W(dPsiN)-W(-dPsiN))
% which is IDENTICALLY ZERO for both dPsiN>=0 and dPsiN<0 -- confirms the
% same equilibrium-consistency property the old code's own comments
% claimed (W(x)/W(-x) = exp(-x) = B(x)/B(-x), same ratio as the Bernoulli
% weights it replaces), re-derived independently here for the QFP
% formulation's psiOffsetN/psiOffsetP convention.
%
% c: constants() struct. vn,vp: thermionic velocities (cm/s), NOT divided
% by any mesh spacing (this IS the whole coefficient, unlike DnSeg/h).
% dPsiN,dPsiP: (psiEffN/P(R)-psiEffN/P(L))/VT, SAME sign convention as
% assembleResidualQFP.m's own dPsiN/dPsiP (L=segment's lower-index node,
% R=lower-index+1). nL,nR,pL,pR: carrier densities at the segment's two
% nodes.
WnDPsiN  = thermionicWeight(dPsiN);
WnMDPsiN = thermionicWeight(-dPsiN);
WpDPsiP  = thermionicWeight(dPsiP);
WpMDPsiP = thermionicWeight(-dPsiP);
JnSeg = c.Q*vn .* (nR.*WnDPsiN - nL.*WnMDPsiN);
JpSeg = c.Q*vp .* (pL.*WpDPsiP - pR.*WpMDPsiP);
end

function w = thermionicWeight(x)
w = exp(min(0,-x));
end
