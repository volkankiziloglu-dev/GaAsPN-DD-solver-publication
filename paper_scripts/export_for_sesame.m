function export_for_sesame()
%EXPORT_FOR_SESAME Export paper architectures for an independent Sesame cross-check (external_validation/sesame_hetero_crosscheck.py).
% SRH only (Brad = Auger = 0, the defaults of buildPaperArchitecture) and frontIRV = rearIRV = 0 so that both codes solve the same equations; the MATLAB
% generation profile G(x) (cm^-3 s^-1 at the mesh nodes) is exported so that the optical input is identical. Per node: x, layer index, Eg, chi, eps, Nc, Nv,
% Nd, Na, muN, muP, tau, G. Also the MATLAB J-V (accepted sweep). Cases: nominal HJ-0/RBSF-GaP/RHJ (AM1.5G), short-lifetime RBSF-GaP and RHJ (tau = 10 ns,
% RHJ also 30 ns, i.e. lifetime ratio 3), and RBSF-GaP/RHJ under CIE 6500 K / 1000 lx.
% Output: external_validation/sesame_input/<tag>_nodes.csv, <tag>_mlJV.csv, <tag>_meta.json.
assertCanonicalPath(); warning('off','all');
root = fileparts(fileparts(mfilename('fullpath'))); outDir = fullfile(root,'external_validation','sesame_input'); if ~isfolder(outDir), mkdir(outDir); end
arc = load(fullfile(root,'results','paper_stageA','stageA_comparison.mat'),'d1nm','d2nm');
arcFn = @(lam) tmmReflectance(lam,[1.0,1.38,2.35,gapRefractiveIndexReal(lam)],[arc.d1nm,arc.d2nm]);
specA = loadRealAM15G(); specI = buildCIEIndoorSpectrum('6500K',1000);
%        code        tau_abs   tag            spectrum
cases = {'HJ-0',     0.6e-6,  'HJ-0',        'A'; ...
         'RBSF-GaP', 0.6e-6,  'RBSF-GaP',    'A'; ...
         'RHJ',      0.6e-6,  'RHJ',         'A'; ...
         'RBSF-GaP', 10e-9,   'RBSF-GaP_tau10ns', 'A'; ...
         'RHJ',      10e-9,   'RHJ_tau10ns', 'A'; ...
         'RHJ',      30e-9,   'RHJ_tau30ns', 'A'; ...
         'RBSF-GaP', 0.6e-6,  'RBSF-GaP_indoor', 'I'; ...
         'RHJ',      0.6e-6,  'RHJ_indoor',  'I'};
for ic = 1:size(cases,1)
  code = cases{ic,1}; tau = cases{ic,2}; tag = cases{ic,3};
  if cases{ic,4} == 'A', spec = specA; Vg = (0:0.02:1.70)'; else, spec = specI; Vg = (0:0.02:1.56)'; end
  Pin = spectrumPower(spec);
  dev = buildPaperArchitecture(code,'frontIRV',0,'rearIRV',0,'AbsorberTau_s',tau);
  G = deviceGenAndTransmitML(dev,spec,arcFn);
  N = dev.N;
  T = table(dev.x(:), dev.layerIdxOfNode(:), dev.Eg(:), dev.chi(:), dev.eps(:), dev.Nc(:), dev.Nv(:), dev.Nd(:), dev.Na(:), dev.muN(:), dev.muP(:), dev.tau(:), G(:), ...
    'VariableNames',{'x_cm','layer','Eg','chi','eps','Nc','Nv','Nd','Na','muN','muP','tau','G'});
  writetable(T,fullfile(outDir,sprintf('%s_nodes.csv',tag)));
  if cases{ic,4} == 'A'
    [curve,~,m] = runAcceptedJVSweep(dev,G,Vg,Pin);
    cv = curve(:,{'V_V','J_Acm2','Accepted'});
  else
    [cv,m] = indoorCurve(dev,G,Vg,Pin);
  end
  writetable(cv,fullfile(outDir,sprintf('%s_mlJV.csv',tag)));
  meta = struct('code',code,'tag',tag,'N',N,'Sback',dev.Sback,'Pin_Wm2',Pin,'Jsc_mAcm2',m.Jsc*1e3,'Voc_mV',m.Voc*1e3,'FF_pct',m.FF*100,'PCE_pct',m.PCE);
  fid = fopen(fullfile(outDir,sprintf('%s_meta.json',tag)),'w'); fprintf(fid,'%s',jsonencode(meta)); fclose(fid);
  fprintf('%s: Jsc=%.4g mA/cm2 Voc=%.1f FF=%.2f PCE=%.3f (N=%d)\n',tag,m.Jsc*1e3,m.Voc*1e3,m.FF*100,m.PCE,N);
end
end

function [cv,m] = indoorCurve(dev,G,Vg,Pin)
N = dev.N; idx = find(abs(Vg-0.02)<1e-12,1);
r02 = gummelSolveQFPAnalytic(dev,0.02,G,struct('psiInit',dev.psiEq,'phinInit',zeros(N,1),'phipInit',zeros(N,1))); prev = r02;
for v = [0.01 0.005 0], r0 = gummelSolveQFPAnalytic(dev,v,G,struct('psiInit',prev.psi,'phinInit',prev.phin,'phipInit',prev.phip)); prev = r0; end
st = repmat(r02,numel(Vg),1); st(1) = r0; st(idx) = r02; prev = r02;
for k = idx+1:numel(Vg)
  rk = gummelSolveQFPAnalytic(dev,Vg(k),G,struct('psiInit',prev.psi,'phinInit',prev.phin,'phipInit',prev.phip)); st(k) = rk; if rk.isFinite, prev = rk; end
end
J = nan(numel(Vg),1); acc = false(numel(Vg),1);
for k = 1:numel(Vg), if st(k).isFinite, J(k) = mean(st(k).Jseg(N-20:N-1)); acc(k) = qfpResultAccepted(st(k)); end, end
cv = table(Vg,J,double(acc),'VariableNames',{'V_V','J_Acm2','Accepted'});
m = pvMetrics([Vg(acc),J(acc)],Pin);
end
