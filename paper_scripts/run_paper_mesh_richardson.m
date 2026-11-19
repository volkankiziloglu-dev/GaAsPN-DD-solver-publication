function T = run_paper_mesh_richardson(code, scales, tauNs, outName)
%RUN_PAPER_MESH_RICHARDSON Three-level mesh study (1x, 2x, 4x nodes) for one architecture, AM1.5G and CIE 6500 K / 1000 lx.
% Production setup (Stage-A ARC, SRH + radiative + Auger, 0.02 V grid). The full J-V curves are saved (all finite points and the gate flag), so that the terminal
% metrics can be extracted offline with one common extractor (cases/analyze_mesh_richardson.py).
% Output: results/paper_mesh_richardson/<code>_x<scale>_<cond>.csv and mesh_richardson_<code>.csv (class and timing).
if nargin < 2, scales = [1 2 4]; end
if nargin < 3, tauNs = 600; end   % absorber lifetime (ns; RHJ r = 1); 1 and 10 ns in addition to the nominal 600 ns
sfx = ''; if tauNs ~= 600, sfx = sprintf('_tau%g',tauNs); end
assertCanonicalPath(); warning('off','all');
root = fileparts(fileparts(mfilename('fullpath'))); if nargin < 4 || isempty(outName), outName = 'paper_mesh_richardson'; end   % e.g. paper_mesh_richardson_expm1 (re-run after the sgBernoulli fix)
outDir = fullfile(root,'results',outName); if ~isfolder(outDir), mkdir(outDir); end
arc = load(fullfile(root,'results','paper_stageA','stageA_comparison.mat'),'d1nm','d2nm');
arcFn = @(lam) tmmReflectance(lam,[1.0,1.38,2.35,gapRefractiveIndexReal(lam)],[arc.d1nm,arc.d2nm]);
Bmat = paperRadiativeCoefficients(300);
conds = {'AM1.5G'; '6500K_1000lx'};
specs = {loadRealAM15G(); buildCIEIndoorSpectrum('6500K',1000)};
Vgrids = {(0:0.02:1.70)'; (0:0.02:1.56)'};
Pins = cellfun(@spectrumPower,specs);
rows = {};
for sc = scales(:)'
  dev = buildPaperArchitecture(code,'meshScale',sc,'AbsorberTau_s',tauNs*1e-9, ...
    'BradAbs',Bmat.GaAsPN_cm3ps,'BradGaP',Bmat.GaP_cm3ps,'BradGaNP',Bmat.GaPN_Ege18_cm3ps, ...
    'CnAugerAbs',3e-30,'CpAugerAbs',3e-30,'CnAugerGaP',3e-30,'CpAugerGaP',3e-30,'CnAugerGaNP',3e-30,'CpAugerGaNP',3e-30);
  for icd = 1:2
    t0 = tic; G = deviceGenAndTransmitML(dev,specs{icd},arcFn);
    try
      if icd == 1
        [curve,diag] = runAcceptedJVSweep(dev,G,Vgrids{icd},Pins(icd));
        out = table(curve.V_V,curve.J_Acm2,double(curve.Accepted),'VariableNames',{'V_V','J_Acm2','Accepted'});
        out = out(isfinite(out.J_Acm2),:); cls = sprintf('AM15G %d/%d',diag.nAccepted,diag.nAttempted);
      else
        gate = indoorSweepLocal(dev,G,Vgrids{icd},Pins(icd));
        f = gate.isFinite(:);
        out = table(gate.Vs(f),gate.Jterm(f),double(gate.pointPass(f)),'VariableNames',{'V_V','J_Acm2','Accepted'});
        cls = char(gate.curveClass);
      end
      writetable(out,fullfile(outDir,sprintf('%s%s_x%g_%s.csv',code,sfx,sc,conds{icd})));
    catch ME
      cls = ['ERR ' strrep(ME.message,',',';')];
    end
    rows(end+1,:) = {code, sc, dev.N, conds{icd}, cls, toc(t0)}; %#ok<AGROW>
    T = cell2table(rows,'VariableNames',{'Architecture','MeshScale','Nodes','Condition','Class','Seconds'});
    writetable(T,fullfile(outDir,sprintf('mesh_richardson_%s%s.csv',code,sfx)));
    disp([code ' x' num2str(sc) ' ' conds{icd} ' ' cls ' (' num2str(toc(t0),'%.0f') ' s)']);
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
