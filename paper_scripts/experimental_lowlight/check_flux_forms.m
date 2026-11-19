function check_flux_forms(code, lux, Vpick)
%CHECK_FLUX_FORMS EXPERIMENTAL diagnostic: compares production difference, expm1, excess-carrier and vpa fluxes on the same (double) state.
if nargin < 1, code = "RHJ"; end
if nargin < 2, lux = 100; end
if nargin < 3, Vpick = 0.5; end
root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(fullfile(root,'src','experimental_lowlight'));
Bmat = paperRadiativeCoefficients(300);
dev = buildPaperArchitecture(code, BradAbs=Bmat.GaAsPN_cm3ps, BradGaP=Bmat.GaP_cm3ps, BradGaNP=Bmat.GaPN_Ege18_cm3ps, ...
  CnAugerAbs=3e-30,CpAugerAbs=3e-30,CnAugerGaP=3e-30,CpAugerGaP=3e-30,CnAugerGaNP=3e-30,CpAugerGaNP=3e-30);
L = load(fullfile(root,'results','experimental_lowlight','states',sprintf('%s_%glx_prod.mat',code,lux)));
[~,k] = min(abs(L.Vgrid-Vpick)); r = L.states(k); V = L.Vgrid(k); N = dev.N; G = zeros(N,1);
modes = {'prod','expm1','excess'}; J = struct();
for m = modes
  d = dev; d.fluxModeX = m{1}; [~,~,~,Jn,Jp] = assembleResidualQFPX(d, r.psi, r.phin, r.phip, V, G); J.(m{1}) = [Jn Jp];
end
S = struct('psi',r.psi,'psiLo',zeros(N,1),'phin',r.phin,'phinLo',zeros(N,1),'phip',r.phip,'phipLo',zeros(N,1));
tic; [~,~,~,Jn,Jp] = assembleResidualQFPHP(dev, S, V, G, 'vpadiff', 40); tv = toc; J.vpa = [Jn Jp];
[~,~,~,Jn,Jp] = assembleResidualQFPHP(dev, S, V, G, 'ddexcess'); J.dd = [Jn Jp];
fprintf('%s %g lx, V=%.2f, N=%d, vpa evaluation %.2f s, |J| ~ %.3g A/cm2\n', code, lux, V, N, tv, abs(r.J));
ref = J.vpa;
for m = {'prod','expm1','excess','dd'}
  dJ = J.(m{1}) - ref;
  [mx,iw] = max(abs(dJ(:,1))); [mxp,iwp] = max(abs(dJ(:,2)));
  fprintf('  %-7s vs vpa: max|dJn| %.3e (seg %d)  max|dJp| %.3e (seg %d)\n', m{1}, mx, iw, mxp, iwp);
end
dpp = J.prod - J.excess; fprintf('  prod vs excess: max|dJn| %.3e max|dJp| %.3e\n', max(abs(dpp(:,1))), max(abs(dpp(:,2))));
end
