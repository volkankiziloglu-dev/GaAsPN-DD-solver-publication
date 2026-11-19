function T = run_paper_interface_mesh(slice)
%RUN_PAPER_INTERFACE_MESH Local mesh refinement at the heterointerfaces.
% The abrupt band step lies inside one Scharfetter-Gummel segment whose width is the node spacing of the left layer (1 nm at the
% n+-GaP/GaAsPN front interface, 6.25 nm at the GaAsPN/p+-GaP rear interface). buildPaperArchitecture(...,'InterfaceSlabNodes',ns) splits
% a 5 nm slab of the same material off each side of every heterointerface and meshes it with ns nodes (25 -> 0.2 nm, 50 -> 0.1 nm).
% As the segment shrinks, D/h grows (about 2.6e7 cm/s at 5 nm, 6e8 cm/s at 0.2 nm), i.e. the SG flux approaches the limit of fast
% interface transport; a change of the metrics with refinement would show interface-transport (mesh) dependence.
% Jobs: HJ-0, F-HJ, RBSF-GaP, RHJ x chi_GaP in {3.8 (dEc = dEv = 0.20 eV), 3.65 (dEc = 0.35, dEv = 0.05 eV)} x slab nodes {0, 25, 50}
%       x tau_p in {600, 10} ns x {AM1.5G, 6500 K / 1000 lx}; production sweeps, unified extractor.
% Output: results/paper_interface_mesh/interface_mesh_<k>of<n>.csv
if nargin < 1 || isempty(slice), slice = [1 1]; end
assertCanonicalPath(); warning('off','all');
root = fileparts(fileparts(mfilename('fullpath'))); outDir = fullfile(root,'results','paper_interface_mesh'); if ~isfolder(outDir), mkdir(outDir); end
arc = load(fullfile(root,'results','paper_stageA','stageA_comparison.mat'),'d1nm','d2nm');
arcFn = @(lam) tmmReflectance(lam,[1.0,1.38,2.35,gapRefractiveIndexReal(lam)],[arc.d1nm,arc.d2nm]);
Bmat = paperRadiativeCoefficients(300);
conds = {'AM1.5G'; '6500K_1000lx'};
specs = {loadRealAM15G(); buildCIEIndoorSpectrum('6500K',1000)};
Vgrids = {(0:0.02:1.70)'; (0:0.02:1.56)'};
Pins = cellfun(@spectrumPower,specs);
common = {'BradAbs',Bmat.GaAsPN_cm3ps,'BradGaP',Bmat.GaP_cm3ps,'BradGaNP',Bmat.GaPN_Ege18_cm3ps, ...
  'CnAugerAbs',3e-30,'CpAugerAbs',3e-30,'CnAugerGaP',3e-30,'CpAugerGaP',3e-30,'CnAugerGaNP',3e-30,'CpAugerGaNP',3e-30};
jobs = {};
for tau = [600 10]*1e-9
  for chi = [3.8 3.65]
    for ns = [0 25 50]
      for code = ["HJ-0","F-HJ","RBSF-GaP","RHJ"], jobs(end+1,:) = {char(code), tau, chi, ns}; end %#ok<AGROW>
    end
  end
end
jobs = jobs(slice(1):slice(2):end,:);
rows = {}; ck = fullfile(outDir,sprintf('interface_mesh_%dof%d.csv',slice(1),slice(2)));
for j = 1:size(jobs,1)
  [code,tauP,chi,ns] = jobs{j,:};
  for icd = 1:2
    t0 = tic; nNodes = NaN;
    try
      dev = buildPaperArchitecture(code,'AbsorberTau_s',tauP,'GaPChi_eV',chi,'InterfaceSlabNodes',ns,common{:}); nNodes = dev.N;
      G = deviceGenAndTransmitML(dev,specs{icd},arcFn);
      if icd == 1
        [~,~,m] = runAcceptedJVSweep(dev,G,Vgrids{icd},Pins(icd),true); cls = 'AM15G';
        Jsc = m.Jsc*1e3; Voc = m.Voc*1e3; FF = m.FF*100; PCE = m.PCE;
      else
        gate = indoorSweepLocal(dev,G,Vgrids{icd},Pins(icd)); cls = char(gate.curveClass);
        Jsc = gate.Jsc*1e6; Voc = gate.Voc*1e3; FF = gate.FF*100; PCE = gate.PCE;
      end
    catch ME
      Jsc = NaN; Voc = NaN; FF = NaN; PCE = NaN; cls = ['ERR ' strrep(ME.message,',',';')];
    end
    rows(end+1,:) = {code, tauP, chi, ns, nNodes, conds{icd}, Jsc, Voc, FF, PCE, cls}; %#ok<AGROW>
    T = cell2table(rows,'VariableNames',{'Architecture','tauP_s','chiGaP','SlabNodes','Nodes','Condition','Jsc','Voc_mV','FF_pct','PCE_pct','Class'});
    writetable(T,ck);
    fprintf('%s tau=%.1e chi=%.2f ns=%d %s: PCE=%.4f %s (%.0f s)\n',code,tauP,chi,ns,conds{icd},PCE,cls,toc(t0));
  end
end
end

function gate = indoorSweepLocal(dev,G,Vgrid,Pin)
% Identical to run_paper_referee_checks>indoorSweepLocal (unified extractor).
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
fU = gate.isFinite(:); mU = pvMetrics([gate.Vs(fU), gate.Jterm(fU)],Pin);
gate.Jsc = mU.Jsc; gate.Voc = mU.Voc; gate.FF = mU.FF; gate.PCE = mU.PCE;
end
