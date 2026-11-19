function T = run_recombination_decomposition(codes, tauAbs, tag)
%RUN_RECOMBINATION_DECOMPOSITION Current balance of HJ-0, RBSF-GaP and RHJ under AM1.5G (production devices, recombination,
% ARC and QFP solver as run_paper_spatial_profiles.m). At each voltage the terminal current is split exactly as
%   J = q*int(G) - q*int(R_SRH + R_rad + R_Aug) - q*int(R_int) - J_rear,minority - J_front,minority
% with the integrals taken over the control volumes of the interior nodes, i.e. the discrete balance of the solver's own
% continuity rows (assembleResidualQFP.m). All devices have the n-type contact at the front and the p-type contact at the
% rear, so the minority-carrier losses are holes at the front and electrons at the rear contact.
% Voltages: 0 V (J_sc), the V_oc of HJ-0 (common voltage for the V_oc comparison) and each device's own V_oc.
% Output: results/paper_recombination_decomposition/decomposition_AM15G.csv (+ per-layer csv).
% Optional: codes (string array; the first one sets the common comparison voltage), tauAbs (absorber SRH lifetime, s),
% tag (file suffix). Defaults reproduce the production run: ["HJ-0","RBSF-GaP","RHJ"], 0.6e-6, "".
if nargin < 1 || isempty(codes), codes = ["HJ-0","RBSF-GaP","RHJ"]; end
if nargin < 2 || isempty(tauAbs), tauAbs = 0.6e-6; end
if nargin < 3, tag = ""; end
assertCanonicalPath();
warning('off','MATLAB:nearlySingularMatrix'); warning('off','MATLAB:singularMatrix');
root = fileparts(fileparts(mfilename('fullpath')));
outDir = fullfile(root,'results','paper_recombination_decomposition'); if ~isfolder(outDir), mkdir(outDir); end
arc = load(fullfile(root,'results','paper_stageA','stageA_comparison.mat'),'d1nm','d2nm');
arcFn = @(lam) tmmReflectance(lam,[1.0,1.38,2.35,gapRefractiveIndexReal(lam)],[arc.d1nm,arc.d2nm]);
Bmat = paperRadiativeCoefficients(300);
spec = loadRealAM15G(); c = constants();
VocHJ0 = NaN; rows = {}; lrows = {};
for ic = 1:numel(codes)
  [dev,meta] = buildPaperArchitecture(codes(ic), BradAbs=Bmat.GaAsPN_cm3ps, BradGaP=Bmat.GaP_cm3ps, BradGaNP=Bmat.GaPN_Ege18_cm3ps, ...
    CnAugerAbs=3e-30,CpAugerAbs=3e-30,CnAugerGaP=3e-30,CpAugerGaP=3e-30,CnAugerGaNP=3e-30,CpAugerGaNP=3e-30,AbsorberTau_s=tauAbs);
  G = deviceGenAndTransmitML(dev,spec,arcFn); G = G(:);
  [s0, sOc, Voc, sweepStates, Vg] = solveStates(dev,G);
  if ic == 1, VocHJ0 = Voc; end
  sets = {s0, 0, 'Jsc'; sOc, Voc, 'Voc_own'};
  if ic > 1
    sC = solveAt(dev,G,VocHJ0,sweepStates,Vg); sets(end+1,:) = {sC, VocHJ0, 'V_HJ0oc'};
  end
  for k = 1:size(sets,1)
    [b, L] = balance(dev,meta,G,sets{k,1},c);
    rows(end+1,:) = {char(codes(ic)), sets{k,3}, sets{k,2}, b.J, b.Gen, b.SRH, b.Rad, b.Aug, b.Int, b.Rear, b.Front, b.Closure}; %#ok<AGROW>
    for j = 1:numel(L.name)
      lrows(end+1,:) = {char(codes(ic)), sets{k,3}, sets{k,2}, L.name{j}, L.Gen(j), L.SRH(j), L.Rad(j), L.Aug(j), L.Int(j)}; %#ok<AGROW>
    end
    fprintf('%-9s %-8s V=%.4f  J=%.4f | Gen %.4f SRH %.4f Rad %.4f Aug %.4f Int %.4f Rear %.4f Front %.4f | closure %.2e (mA/cm2)\n', ...
      codes(ic), sets{k,3}, sets{k,2}, b.J, b.Gen, b.SRH, b.Rad, b.Aug, b.Int, b.Rear, b.Front, b.Closure);
  end
end
T = cell2table(rows,'VariableNames',{'Architecture','Point','V_V','J_mAcm2','Gen_mAcm2','SRH_mAcm2','Rad_mAcm2','Aug_mAcm2','Interface_mAcm2','RearContact_mAcm2','FrontContact_mAcm2','Closure_mAcm2'});
writetable(T,fullfile(outDir,'decomposition_AM15G' + string(tag) + '.csv'));
TL = cell2table(lrows,'VariableNames',{'Architecture','Point','V_V','Layer','Gen_mAcm2','SRH_mAcm2','Rad_mAcm2','Aug_mAcm2','Interface_mAcm2'});
writetable(TL,fullfile(outDir,'decomposition_AM15G_layers' + string(tag) + '.csv')); disp(T);
end

function [b, L] = balance(dev,meta,G,s,c)
N = dev.N; x = dev.x(:); h = diff(x); idxC = (2:N-1)'; hAvg = 0.5*(h(idxC-1)+h(idxC));
n = s.n(:); p = s.p(:); npx = n(idxC).*p(idxC) - dev.ni(idxC).^2;
Rsrh = bulkSRHRate(dev,n,p,idxC); Rrad = dev.Brad(idxC).*npx; Raug = (dev.CnAuger(idxC).*n(idxC)+dev.CpAuger(idxC).*p(idxC)).*npx;
Rint = interfaceRecombination(dev,n,p,idxC,hAvg);
w = c.Q*hAvg*1e3;                                  % A/cm2 per (cm-3 s-1) -> mA/cm2
b.Gen = sum(w.*G(idxC)); b.SRH = sum(w.*Rsrh); b.Rad = sum(w.*Rrad); b.Aug = sum(w.*Raug); b.Int = sum(w.*Rint);
b.Rear = -s.JnSeg(N-1)*1e3; b.Front = -s.JpSeg(1)*1e3; b.J = s.J*1e3;
b.Closure = b.Gen - b.SRH - b.Rad - b.Aug - b.Int - b.Rear - b.Front - b.J;
lsn = [dev.layerStartNode(:); N+1]; L.name = cellstr(string(meta.layerNames));
for j = 1:numel(L.name)
  in = idxC >= lsn(j) & idxC < lsn(j+1);
  L.Gen(j) = sum(w(in).*G(idxC(in))); L.SRH(j) = sum(w(in).*Rsrh(in)); L.Rad(j) = sum(w(in).*Rrad(in));
  L.Aug(j) = sum(w(in).*Raug(in)); L.Int(j) = sum(w(in).*Rint(in));
end
end

function [s0, sOc, Voc, states, Vg] = solveStates(dev,G)
N = dev.N; Vg = (0.02:0.02:1.56)'; ws = @(s) struct('psiInit',s.psi,'phinInit',s.phin,'phipInit',s.phip);
r = gummelSolveQFPAnalytic(dev,0.02,G,struct('psiInit',dev.psiEq,'phinInit',zeros(N,1),'phipInit',zeros(N,1)));
prev = r; for v = [0.01 0.005 0], s0 = gummelSolveQFPAnalytic(dev,v,G,ws(prev)); prev = s0; end
states = cell(numel(Vg),1); states{1} = r; prev = r; J = nan(numel(Vg),1); J(1) = r.J;
for k = 2:numel(Vg)
  rk = gummelSolveQFPAnalytic(dev,Vg(k),G,ws(prev)); states{k} = rk;
  if rk.isFinite, prev = rk; J(k) = rk.J; end
  if rk.isFinite && rk.J < 0, break; end
end
kpos = find(J > 0,1,'last'); V1 = Vg(kpos); J1 = J(kpos); V2 = Vg(kpos+1); J2 = J(kpos+1); prev = states{kpos};
for it = 1:10
  Voc = V1 - J1*(V2-V1)/(J2-J1); sOc = gummelSolveQFPAnalytic(dev,Voc,G,ws(prev));
  if abs(sOc.J) < 1e-5*abs(s0.J), break; end
  if sOc.J > 0, V1 = Voc; J1 = sOc.J; prev = sOc; else, V2 = Voc; J2 = sOc.J; end
end
end

function s = solveAt(dev,G,V,states,Vg)
ws = @(s) struct('psiInit',s.psi,'phinInit',s.phin,'phipInit',s.phip);
k = find(Vg <= V,1,'last'); prev = states{k};
s = gummelSolveQFPAnalytic(dev,V,G,ws(prev));
end
