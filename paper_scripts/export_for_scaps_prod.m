function export_for_scaps_prod()
%EXPORT_FOR_SCAPS_PROD Export of the PRODUCTION configuration for the SCAPS-1D cross-check (simulated referee R1-1, docs/HAKEM_SIMULASYONU_2026-11-06.md).
% Unlike export_for_sesame.m (SRH only, ideal heterointerfaces), the devices keep the production recombination: bulk SRH, radiative
% (B_rad) and Auger (3e-30 cm6/s) terms, and the GaP/GaAsPN interface recombination velocity of 8910 cm/s at every heterointerface
% (geometric-mean interface densities in the MATLAB model). Per node: x, layer, Eg, chi, eps, Nc, Nv, Nd, Na, muN, muP, tau, Brad, Cn, Cp, G;
% meta: Sback, interface velocities between consecutive layers, MATLAB metrics (unified extractor: all finite points, J_sc at 0 V,
% linear V_oc, parabolic P_max) and J-V.
% Cases (<= 6 layers for the SCAPS builder): RBSF-GaP, RHJ, F-HJ at AM1.5G (0.6 us and 10 ns) and 6500 K at 1000, 200 and 100 lx (0.6 us);
% the 100 and 200 lx cases check the low-illuminance curves that the residual diagnostic classifies as SCREEN (referee R3-1).
% Output: external_validation/scaps_prod_input/<tag>_nodes.csv, <tag>_mlJV.csv, <tag>_meta.json.
assertCanonicalPath(); warning('off','all');
root = fileparts(fileparts(mfilename('fullpath'))); outDir = fullfile(root,'external_validation','scaps_prod_input'); if ~isfolder(outDir), mkdir(outDir); end
arc = load(fullfile(root,'results','paper_stageA','stageA_comparison.mat'),'d1nm','d2nm');
arcFn = @(lam) tmmReflectance(lam,[1.0,1.38,2.35,gapRefractiveIndexReal(lam)],[arc.d1nm,arc.d2nm]);
Bmat = paperRadiativeCoefficients(300);
common = {'BradAbs',Bmat.GaAsPN_cm3ps,'BradGaP',Bmat.GaP_cm3ps,'BradGaNP',Bmat.GaPN_Ege18_cm3ps, ...
  'CnAugerAbs',3e-30,'CpAugerAbs',3e-30,'CnAugerGaP',3e-30,'CpAugerGaP',3e-30,'CnAugerGaNP',3e-30,'CpAugerGaNP',3e-30};
specA = loadRealAM15G(); specI = buildCIEIndoorSpectrum('6500K',1000);
cases = {'RBSF-GaP', 0.6e-6, 'P_RBSF-GaP', 'A'; ...
         'RHJ',      0.6e-6, 'P_RHJ',      'A'; ...
         'F-HJ',     0.6e-6, 'P_F-HJ',     'A'; ...
         'RBSF-GaP', 10e-9,  'P_RBSF-GaP_tau10ns', 'A'; ...
         'RHJ',      10e-9,  'P_RHJ_tau10ns',      'A'; ...
         'RBSF-GaP', 0.6e-6, 'P_RBSF-GaP_indoor', 'I'; ...
         'RHJ',      0.6e-6, 'P_RHJ_indoor',      'I'; ...
         'F-HJ',     0.6e-6, 'P_F-HJ_indoor',     'I'; ...
         'RBSF-GaP', 0.6e-6, 'P_RBSF-GaP_indoor200', '2'; ...
         'RHJ',      0.6e-6, 'P_RHJ_indoor200',      '2'; ...
         'F-HJ',     0.6e-6, 'P_F-HJ_indoor200',     '2'; ...
         'RBSF-GaP', 0.6e-6, 'P_RBSF-GaP_indoor100', '1'; ...
         'RHJ',      0.6e-6, 'P_RHJ_indoor100',      '1'; ...
         'F-HJ',     0.6e-6, 'P_F-HJ_indoor100',     '1'};
for ic = 1:size(cases,1)
  [code,tau,tag,sp] = cases{ic,:};
  switch sp
    case 'A', spec = specA; Vg = (0:0.02:1.70)';
    case 'I', spec = specI; Vg = (0:0.02:1.56)';
    case '2', spec = buildCIEIndoorSpectrum('6500K',200); Vg = (0:0.02:1.56)';
    case '1', spec = buildCIEIndoorSpectrum('6500K',100); Vg = (0:0.02:1.56)';
  end
  Pin = spectrumPower(spec);
  [dev,meta0] = buildPaperArchitecture(code,'AbsorberTau_s',tau,common{:});
  G = deviceGenAndTransmitML(dev,spec,arcFn);
  T = table(dev.x(:), dev.layerIdxOfNode(:), dev.Eg(:), dev.chi(:), dev.eps(:), dev.Nc(:), dev.Nv(:), dev.Nd(:), dev.Na(:), dev.muN(:), dev.muP(:), dev.tau(:), ...
    dev.Brad(:), dev.CnAuger(:), dev.CpAuger(:), G(:), ...
    'VariableNames',{'x_cm','layer','Eg','chi','eps','Nc','Nv','Nd','Na','muN','muP','tau','Brad','Cn','Cp','G'});
  writetable(T,fullfile(outDir,sprintf('%s_nodes.csv',tag)));
  [cv,m] = sweepAll(dev,G,Vg,Pin);
  writetable(cv,fullfile(outDir,sprintf('%s_mlJV.csv',tag)));
  Sint = arrayfun(@(j) j.Sn_int, meta0.junctions);
  meta = struct('code',code,'tag',tag,'N',dev.N,'Sback',dev.Sback,'Sint_cmps',Sint,'Pin_Wm2',Pin, ...
    'Jsc_mAcm2',m.Jsc*1e3,'Voc_mV',m.Voc*1e3,'FF_pct',m.FF*100,'PCE_pct',m.PCE);
  fid = fopen(fullfile(outDir,sprintf('%s_meta.json',tag)),'w'); fprintf(fid,'%s',jsonencode(meta)); fclose(fid);
  fprintf('%s: Jsc=%.5g mA/cm2 Voc=%.2f FF=%.3f PCE=%.4f (N=%d, Sint=%s)\n',tag,m.Jsc*1e3,m.Voc*1e3,m.FF*100,m.PCE,dev.N,mat2str(Sint));
end
end

function [cv,m] = sweepAll(dev,G,Vg,Pin)
% Production continuation; terminal current of every finite point (unified extractor).
N = dev.N; idx = find(abs(Vg-0.02)<1e-12,1);
r02 = gummelSolveQFPAnalytic(dev,0.02,G,struct('psiInit',dev.psiEq,'phinInit',zeros(N,1),'phipInit',zeros(N,1))); prev = r02;
for v = [0.01 0.005 0], r0 = gummelSolveQFPAnalytic(dev,v,G,struct('psiInit',prev.psi,'phinInit',prev.phin,'phipInit',prev.phip)); prev = r0; end
st = repmat(r02,numel(Vg),1); st(1) = r0; st(idx) = r02; prev = r02;
for k = idx+1:numel(Vg)
  rk = gummelSolveQFPAnalytic(dev,Vg(k),G,struct('psiInit',prev.psi,'phinInit',prev.phin,'phipInit',prev.phip)); st(k) = rk; if rk.isFinite, prev = rk; end
end
J = nan(numel(Vg),1); fin = false(numel(Vg),1);
for k = 1:numel(Vg), if st(k).isFinite, J(k) = st(k).J; fin(k) = true; end, end
cv = table(Vg,J,double(fin),'VariableNames',{'V_V','J_Acm2','Finite'});
m = pvMetrics([Vg(fin),J(fin)],Pin);
end
