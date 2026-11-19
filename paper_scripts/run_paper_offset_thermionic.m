function T = run_paper_offset_thermionic(job, tag)
%RUN_PAPER_OFFSET_THERMIONIC Robustness check: conduction-band-dominated band offset and thermionic-emission-limited heterointerface transport.
% job.chi        list of GaP electron affinities (eV); absorber chi = 4.0 eV, so chi = 3.8 gives dEc = 0.20 / dEv = 0.20 eV (nominal) and
%                chi = 3.65 gives dEc = 0.35 / dEv = 0.05 eV (conduction-band dominated, as expected from the band-anticrossing picture).
% job.codes      architectures; job.v  list of thermionic velocities (cm/s, vn = vp), NaN = standard Scharfetter-Gummel flux (reference).
% job.tauP_ns, job.r   p-type lifetime and RHJ n/p lifetime ratio (RHJ absorber lifetime r*tauP).
% Thermionic flux (thermionicFluxQFP.m) is applied at the segments of the heteromaterial junctions only (chi or Eg jump across the junction).
% All curves are evaluated with the same extractor (all finite points, Jsc = J(0), linear Voc, grid maximum of V*J), so that the thermionic and the
% reference curves of one device are directly comparable; the AM1.5G and indoor production values differ slightly (0.02 V grid, extractor).
% Output: results/paper_offset_thermionic/offset_thermionic_<tag>.csv
assertCanonicalPath();
warning('off','MATLAB:nearlySingularMatrix'); warning('off','MATLAB:singularMatrix');
root = fileparts(fileparts(mfilename('fullpath'))); outDir = fullfile(root,'results','paper_offset_thermionic'); if ~isfolder(outDir), mkdir(outDir); end
arc = load(fullfile(root,'results','paper_stageA','stageA_comparison.mat'),'d1nm','d2nm');
arcFn = @(lam) tmmReflectance(lam,[1.0,1.38,2.35,gapRefractiveIndexReal(lam)],[arc.d1nm,arc.d2nm]);
Bmat = paperRadiativeCoefficients(300);
conds = {'AM1.5G'; '6500K_1000lx'};
specs = {loadRealAM15G(); buildCIEIndoorSpectrum('6500K',1000)};
% coarse grids (the thermionic solver uses a finite-difference Jacobian and is slow): 0.1 V steps up to the MPP region, 0.02 V steps around the MPP and Voc
Vgrids = {[0 0.02 0.1:0.1:1.2 1.22:0.02:1.52]'; [0 0.02 0.1:0.1:1.1 1.12:0.02:1.40]'};
Pins = cellfun(@spectrumPower,specs);
if ~isfield(job,'conds'), job.conds = [1 2]; end
ck = fullfile(outDir,sprintf('offset_thermionic_%s.csv',tag)); rows = {};
for chi = job.chi(:)'
  for ic = 1:numel(job.codes)
    code = job.codes{ic}; r = 1; if strcmp(code,'RHJ') && isfield(job,'r'), r = job.r; end
    tauAbs = job.tauP_ns*1e-9*r;
    dev0 = buildPaperArchitecture(code,'GaPChi_eV',chi,'AbsorberTau_s',tauAbs, ...
      'BradAbs',Bmat.GaAsPN_cm3ps,'BradGaP',Bmat.GaP_cm3ps,'BradGaNP',Bmat.GaPN_Ege18_cm3ps, ...
      'CnAugerAbs',3e-30,'CpAugerAbs',3e-30,'CnAugerGaP',3e-30,'CpAugerGaP',3e-30,'CnAugerGaNP',3e-30,'CpAugerGaNP',3e-30);
    jn = dev0.junctionNodes; hetero = false(size(jn));
    for k = 1:numel(jn)
      a = jn(k); b = jn(k)+1;
      hetero(k) = abs(dev0.chi(a)-dev0.chi(b)) > 1e-6 || abs(dev0.Eg(a)-dev0.Eg(b)) > 1e-6;
    end
    for icd = job.conds
      G = deviceGenAndTransmitML(dev0,specs{icd},arcFn);
      for v = job.v(:)'
        dev = dev0; t0 = tic;
        if isnan(v), solver = @gummelSolveQFPAnalytic; dev.thermionicJunction = [];
        else, solver = @gummelSolveQFPThermionic; dev.thermionicJunction = jn(hetero); dev.vn_thermionic = v; dev.vp_thermionic = v; end
        try
          [Vv,Jv,nfin] = sweepLocal(dev,G,Vgrids{icd},solver);
          m = metricsLocal(Vv,Jv,Pins(icd));
          Jsc = m(1)*1e3; Voc = m(2)*1e3; FF = m(3)*100; PCE = m(4); cls = 'ok';
        catch ME
          Jsc = NaN; Voc = NaN; FF = NaN; PCE = NaN; nfin = 0; cls = ['ERR ' strrep(ME.message,',',';')];
        end
        rows(end+1,:) = {chi, code, job.tauP_ns, r, conds{icd}, v, nnz(hetero), Jsc, Voc, FF, PCE, nfin, toc(t0), cls}; %#ok<AGROW>
        Tt = cell2table(rows,'VariableNames',{'chiGaP','Architecture','tauP_ns','ratio','Condition','v_cm_s','nHeteroJunctions','Jsc','Voc_mV','FF_pct','PCE_pct','nFinite','Seconds','Status'});
        writetable(Tt,ck);
        fprintf('chi=%.2f %s %s v=%g: Voc=%.1f PCE=%.3f (%.0f s)\n',chi,code,conds{icd},v,Voc,PCE,toc(t0));
      end
    end
  end
end
T = Tt;
end

function [Vv,Jv,nfin] = sweepLocal(dev,G,Vgrid,solver)
N = dev.N; idxStart = find(abs(Vgrid-0.02)<1e-12,1);
seed = struct('psiInit',dev.psiEq,'phinInit',zeros(N,1),'phipInit',zeros(N,1));
r = solver(dev,0.02,G,seed); prev = r;
st = repmat(r,numel(Vgrid),1); fin = false(numel(Vgrid),1); st(idxStart) = r; fin(idxStart) = r.isFinite;
for v = [0.01 0.005 0]
  r0 = solver(dev,v,G,struct('psiInit',prev.psi,'phinInit',prev.phin,'phipInit',prev.phip)); prev = r0;
  k = find(abs(Vgrid-v)<1e-12,1); if ~isempty(k), st(k) = r0; fin(k) = r0.isFinite; end
end
st(1) = r0; fin(1) = r0.isFinite; prev = r;
for k = idxStart+1:numel(Vgrid)
  if ~prev.isFinite, break; end
  rk = solver(dev,Vgrid(k),G,struct('psiInit',prev.psi,'phinInit',prev.phin,'phipInit',prev.phip));
  st(k) = rk; fin(k) = rk.isFinite; if rk.isFinite, prev = rk; end
end
J = nan(numel(Vgrid),1);
for k = 1:numel(Vgrid), if fin(k), J(k) = mean(st(k).Jseg(N-20:N-1)); end, end
Vv = Vgrid(fin); Jv = J(fin); nfin = numel(Vv);
end

function m = metricsLocal(V,J,Pin)
[V,o] = sort(V); J = J(o);
Jsc = J(1); Voc = NaN;
for k = 1:numel(V)-1
  if V(k) >= 0 && J(k) >= 0 && J(k+1) < 0, Voc = V(k)+(0-J(k))/(J(k+1)-J(k))*(V(k+1)-V(k)); break; end
end
if isnan(Voc) || Jsc <= 0, m = [Jsc Voc NaN NaN]; return; end
P = 0; for k = 1:numel(V), if V(k) >= 0 && V(k) <= Voc && J(k) >= 0, P = max(P,V(k)*J(k)); end, end
m = [Jsc Voc P/(Voc*Jsc) P*1e4/Pin*100];
end
