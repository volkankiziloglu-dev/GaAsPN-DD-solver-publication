function T = run_paper_mesh_cvgen(code, scales)
%RUN_PAPER_MESH_CVGEN Mesh study with control-volume-integrated generation.
% The production generation is sampled at the nodes, G_i = alpha_i*phi(x_i), and multiplied by the control-volume width; its error is
% O(alpha*h), i.e. first order in the mesh. Here G is the exact Beer-Lambert absorption integrated over each control volume
% (deviceGenAndTransmitML(...,true)), everything else as run_paper_mesh_richardson.m (production setup, 0.02 V grid, unified extractor:
% all finite points, J_sc at 0 V, linear V_oc, parabolic P_max). Also records the generated current q*int(G) for both generation forms.
% Output: results/paper_mesh_cvgen/mesh_cvgen_<code>.csv.
if nargin < 2, scales = [0.5 1 2 4]; end
assertCanonicalPath(); warning('off','all');
root = fileparts(fileparts(mfilename('fullpath'))); outDir = fullfile(root,'results','paper_mesh_cvgen'); if ~isfolder(outDir), mkdir(outDir); end
arc = load(fullfile(root,'results','paper_stageA','stageA_comparison.mat'),'d1nm','d2nm');
arcFn = @(lam) tmmReflectance(lam,[1.0,1.38,2.35,gapRefractiveIndexReal(lam)],[arc.d1nm,arc.d2nm]);
Bmat = paperRadiativeCoefficients(300); c = constants();
conds = {'AM1.5G'; '6500K_1000lx'};
specs = {loadRealAM15G(); buildCIEIndoorSpectrum('6500K',1000)};
Vgrids = {(0:0.02:1.70)'; (0:0.02:1.56)'};
Pins = cellfun(@spectrumPower,specs);
rows = {};
for sc = scales(:)'
  dev = buildPaperArchitecture(code,'meshScale',sc, ...
    'BradAbs',Bmat.GaAsPN_cm3ps,'BradGaP',Bmat.GaP_cm3ps,'BradGaNP',Bmat.GaPN_Ege18_cm3ps, ...
    'CnAugerAbs',3e-30,'CpAugerAbs',3e-30,'CnAugerGaP',3e-30,'CpAugerGaP',3e-30,'CnAugerGaNP',3e-30,'CpAugerGaNP',3e-30);
  x = dev.x(:); h = diff(x); w = zeros(dev.N,1); w(1:end-1) = w(1:end-1)+h/2; w(2:end) = w(2:end)+h/2;
  for icd = 1:2
    t0 = tic;
    Gn = deviceGenAndTransmitML(dev,specs{icd},arcFn,false);
    G  = deviceGenAndTransmitML(dev,specs{icd},arcFn,true);
    JgenNodal = c.Q*sum(Gn(2:end-1).*w(2:end-1)); JgenCV = c.Q*sum(G(2:end-1).*w(2:end-1));   % interior control volumes (A/cm2)
    try
      if icd == 1
        [~,~,m] = runAcceptedJVSweep(dev,G,Vgrids{icd},Pins(icd),true); cls = 'AM15G';
      else
        gate = indoorSweepLocal(dev,G,Vgrids{icd},Pins(icd)); cls = char(gate.curveClass);
        f = gate.isFinite(:); m = pvMetrics([gate.Vs(f), gate.Jterm(f)],Pins(icd));
      end
      Jsc = m.Jsc; Voc = m.Voc*1e3; FF = m.FF*100; PCE = m.PCE;
    catch ME
      Jsc = NaN; Voc = NaN; FF = NaN; PCE = NaN; cls = ['ERR ' strrep(ME.message,',',';')];
    end
    rows(end+1,:) = {code, sc, dev.N, conds{icd}, JgenNodal, JgenCV, Jsc, Voc, FF, PCE, cls, toc(t0)}; %#ok<AGROW>
    T = cell2table(rows,'VariableNames',{'Architecture','MeshScale','Nodes','Condition','JgenNodal_Acm2','JgenCV_Acm2','Jsc_Acm2','Voc_mV','FF_pct','PCE_pct','Class','Seconds'});
    writetable(T,fullfile(outDir,sprintf('mesh_cvgen_%s.csv',code)));
    fprintf('%s x%g %s: Jsc=%.6e PCE=%.4f %s (%.0f s)\n',code,sc,conds{icd},Jsc,PCE,cls,toc(t0));
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
