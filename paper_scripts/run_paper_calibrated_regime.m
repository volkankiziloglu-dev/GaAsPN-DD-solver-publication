function T = run_paper_calibrated_regime(code)
%RUN_PAPER_CALIBRATED_REGIME Architecture comparison in the recombination regime calibrated to the measured Geisz MD421 GaNPAs cell.
% The staged calibration of run_geisz2002_staged_calibration.m (results/geisz2002_staged_calibration.csv, rank 1) matched J_sc = 5.8 mA/cm2, V_oc = 1.12 V, FF = 60 % and 3.9 % with
%   two SRH centres: midgap (tau_n = tau_p = 0.5 ns) and E_i + 0.25 eV (tau_n = 15 ns, tau_p = 1.5 ns), an effective collection scale of 0.46 on the optical generation, and a terminal shunt R_sh = 631 Ohm cm2 (R_s = 0).
% Here the same centres are placed in the GaAsPN absorber layers of one architecture (code), the centres have fixed capture times, so the minority-carrier lifetimes of the n-type RHJ absorber and of the p-type absorbers follow from the SRH expression;
% all other parameters are nominal. Collection scale 0.46 applied to the generation; AM1.5G also with the terminal shunt. The multi-centre SRH is supported by gummelSolveQFP (finite-difference Jacobian) only, hence the coarse voltage grids.
% Output: results/paper_calibrated_regime/calibrated_regime_<code>.csv and curves.
assertCanonicalPath(); warning('off','all');
root = fileparts(fileparts(mfilename('fullpath'))); outDir = fullfile(root,'results','paper_calibrated_regime'); if ~isfolder(outDir), mkdir(outDir); end
arc = load(fullfile(root,'results','paper_stageA','stageA_comparison.mat'),'d1nm','d2nm');
arcFn = @(lam) tmmReflectance(lam,[1.0,1.38,2.35,gapRefractiveIndexReal(lam)],[arc.d1nm,arc.d2nm]);
Bmat = paperRadiativeCoefficients(300);
ctr(1) = struct('tauN',5e-10,'tauP',5e-10,'EtRelEi_eV',0);
ctr(2) = struct('tauN',1.5e-8,'tauP',1.5e-9,'EtRelEi_eV',0.25);
conds = {'AM1.5G'; '6500K_1000lx'};
specs = {loadRealAM15G(); buildCIEIndoorSpectrum('6500K',1000)};
Vgrids = {[0 0.02 0.1:0.1:1.0 1.02:0.02:1.40]'; [0 0.02 0.1:0.1:0.8 0.82:0.02:1.30]'};
Pins = cellfun(@spectrumPower,specs);
cs = 0.46;
rows = {};
dev = buildPaperArchitecture(code,'AbsorberSRHCenters',ctr, ...
  'BradAbs',Bmat.GaAsPN_cm3ps,'BradGaP',Bmat.GaP_cm3ps,'BradGaNP',Bmat.GaPN_Ege18_cm3ps, ...
  'CnAugerAbs',3e-30,'CpAugerAbs',3e-30,'CnAugerGaP',3e-30,'CpAugerGaP',3e-30,'CnAugerGaNP',3e-30,'CpAugerGaNP',3e-30);
for icd = 1:2
  t0 = tic; G = cs*deviceGenAndTransmitML(dev,specs{icd},arcFn);
  try
    [Vv,Jv,nfin] = sweepLocal(dev,G,Vgrids{icd});
    pts = [Vv,Jv]; m = pvMetrics(pts,Pins(icd));
    writematrix(pts,fullfile(outDir,sprintf('%s_%s.csv',code,conds{icd})));
    rows(end+1,:) = {code, conds{icd}, cs, Inf, m.Jsc*1e3, m.Voc*1e3, m.FF*100, m.PCE, nfin, 'ok'}; %#ok<AGROW>
    if icd == 1
      mt = pvMetrics(applyTerminalResistance(pts,0,631),Pins(icd));
      rows(end+1,:) = {code, conds{icd}, cs, 631, mt.Jsc*1e3, mt.Voc*1e3, mt.FF*100, mt.PCE, nfin, 'ok'}; %#ok<AGROW>
    end
  catch ME
    rows(end+1,:) = {code, conds{icd}, cs, Inf, NaN, NaN, NaN, NaN, 0, ['ERR ' strrep(ME.message,',',';')]}; %#ok<AGROW>
  end
  T = cell2table(rows,'VariableNames',{'Architecture','Condition','CollectionScale','Rsh','Jsc','Voc_mV','FF_pct','PCE_pct','nFinite','Status'});
  writetable(T,fullfile(outDir,sprintf('calibrated_regime_%s.csv',code)));
  disp([code ' ' conds{icd} ' PCE=' num2str(rows{end,8},'%.3f') ' ' rows{end,10} ' (' num2str(toc(t0),'%.0f') ' s)']);
end
end

function [Vv,Jv,nfin] = sweepLocal(dev,G,Vgrid)
N = dev.N; solver = @gummelSolveQFP;
r = solver(dev,0.02,G,struct('psiInit',dev.psiEq,'phinInit',zeros(N,1),'phipInit',zeros(N,1))); prev = r;
fin = false(numel(Vgrid),1); J = nan(numel(Vgrid),1); idx = find(abs(Vgrid-0.02)<1e-12,1);
J(idx) = mean(r.Jseg(N-20:N-1)); fin(idx) = r.isFinite;
for v = [0.01 0.005 0]
  r0 = solver(dev,v,G,struct('psiInit',prev.psi,'phinInit',prev.phin,'phipInit',prev.phip)); prev = r0;
  k = find(abs(Vgrid-v)<1e-12,1); if ~isempty(k), J(k) = mean(r0.Jseg(N-20:N-1)); fin(k) = r0.isFinite; end
end
J(1) = mean(r0.Jseg(N-20:N-1)); fin(1) = r0.isFinite; prev = r;
for k = idx+1:numel(Vgrid)
  if ~prev.isFinite, break; end
  rk = solver(dev,Vgrid(k),G,struct('psiInit',prev.psi,'phinInit',prev.phin,'phipInit',prev.phip));
  fin(k) = rk.isFinite; if rk.isFinite, J(k) = mean(rk.Jseg(N-20:N-1)); prev = rk; end
end
Vv = Vgrid(fin); Jv = J(fin); nfin = numel(Vv);
end
