function T = run_paper_voc_lux(code, tauList_ns)
%RUN_PAPER_VOC_LUX Light-intensity dependence of V_oc (slope dV_oc/dlog10(lux), ideality) of one architecture under the CIE 6500 K spectrum (referee check).
% lux = 50, 100, 200, 500, 1000, 2000, 5000; absorber lifetime tauList_ns (default 600 and 10 ns; the RHJ absorber has the same lifetime, r = 1); production settings,
% all finite points of every curve are saved (cases/analyze_voc_lux.py extracts the metrics with the common extractor).
% Output: results/paper_voc_lux/<code>_tau<ns>_<lux>lx.csv and voc_lux_<code>.csv (class, timing).
if nargin < 2, tauList_ns = [600 10]; end
assertCanonicalPath(); warning('off','all');
root = fileparts(fileparts(mfilename('fullpath'))); outDir = fullfile(root,'results','paper_voc_lux'); if ~isfolder(outDir), mkdir(outDir); end
arc = load(fullfile(root,'results','paper_stageA','stageA_comparison.mat'),'d1nm','d2nm');
arcFn = @(lam) tmmReflectance(lam,[1.0,1.38,2.35,gapRefractiveIndexReal(lam)],[arc.d1nm,arc.d2nm]);
Bmat = paperRadiativeCoefficients(300);
Vgrid = (0:0.02:1.56)';
luxList = [50 100 200 500 1000 2000 5000];
rows = {};
for tau = tauList_ns(:)'
  dev = buildPaperArchitecture(code,'AbsorberTau_s',tau*1e-9, ...
    'BradAbs',Bmat.GaAsPN_cm3ps,'BradGaP',Bmat.GaP_cm3ps,'BradGaNP',Bmat.GaPN_Ege18_cm3ps, ...
    'CnAugerAbs',3e-30,'CpAugerAbs',3e-30,'CnAugerGaP',3e-30,'CpAugerGaP',3e-30,'CnAugerGaNP',3e-30,'CpAugerGaNP',3e-30);
  for lux = luxList
    t0 = tic;
    spec = buildCIEIndoorSpectrum('6500K',lux); Pin = spectrumPower(spec);
    G = deviceGenAndTransmitML(dev,spec,arcFn);
    try
      gate = indoorSweepLocal(dev,G,Vgrid,Pin);
      f = gate.isFinite(:);
      out = table(gate.Vs(f),gate.Jterm(f),double(gate.pointPass(f)),'VariableNames',{'V_V','J_Acm2','Accepted'});
      writetable(out,fullfile(outDir,sprintf('%s_tau%g_%glx.csv',code,tau,lux)));
      cls = char(gate.curveClass);
    catch ME
      cls = ['ERR ' strrep(ME.message,',',';')]; Pin = NaN;
    end
    rows(end+1,:) = {code, tau, lux, Pin, cls, toc(t0)}; %#ok<AGROW>
    T = cell2table(rows,'VariableNames',{'Architecture','tauP_ns','lux','Pin_Wm2','Class','Seconds'});
    writetable(T,fullfile(outDir,sprintf('voc_lux_%s.csv',code)));
    disp([code ' tau=' num2str(tau) ' lux=' num2str(lux) ' ' cls ' (' num2str(toc(t0),'%.0f') ' s)']);
  end
end
end

function gate = indoorSweepLocal(dev,G,Vgrid,Pin)
N = dev.N; idxStart = find(abs(Vgrid-0.02)<1e-12,1);
r002 = gummelSolveQFPAnalytic(dev,0.02,G,struct('psiInit',dev.psiEq,'phinInit',zeros(N,1),'phipInit',zeros(N,1)));
prev = r002;
for v = [0.01,0.005,0], r0 = gummelSolveQFPAnalytic(dev,v,G,struct('psiInit',prev.psi,'phinInit',prev.phin,'phipInit',prev.phip)); prev = r0; end
states = repmat(r002,numel(Vgrid),1); states(1) = r0; states(idxStart) = r002; prev = r002;
for k = idxStart+1:numel(Vgrid)
  rk = gummelSolveQFPAnalytic(dev,Vgrid(k),G,struct('psiInit',prev.psi,'phinInit',prev.phin,'phipInit',prev.phip));
  states(k) = rk; if rk.isFinite, prev = rk; end
end
gate = indoorAcceptanceGate(dev,G,Vgrid,states,Pin);
end
