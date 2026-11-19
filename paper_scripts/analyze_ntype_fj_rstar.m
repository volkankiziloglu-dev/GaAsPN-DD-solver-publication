function T = analyze_ntype_fj_rstar()
%ANALYZE_NTYPE_FJ_RSTAR r* of the RHJ and of the n-type front-junction cell (RBSF-GaP-n) against RBSF-GaP, per tau_p, with the
% interpolation of the design map (linear in log r between the grid ratios). RHJ and RBSF-GaP from the design-map runs
% (results/paper_design_map_unified{,_am}, S_back = 1e5), RBSF-GaP-n from run_ntype_fj_rstar.m.
% Output: results/paper_ntype_front_junction/rstar_ntype_fj.csv
root = fileparts(fileparts(mfilename('fullpath'))); out = fullfile(root,'results','paper_ntype_front_junction');
DM = [readDM(fullfile(root,'results','paper_design_map_unified_am')); readDM(fullfile(root,'results','paper_design_map_unified'))];
DM = DM(DM.Sback==1e5,:);
F = []; d = dir(fullfile(out,'ntype_fj_taun_part*.csv')); for k = 1:numel(d), F = [F; readtable(fullfile(out,d(k).name))]; end %#ok<AGROW>
rList = [1 1.5 2 3 5 10 20 30 60]; conds = ["AM1.5G","6500K_1000lx"]; rows = {};
for c = conds
  for tp = [1 3 10 30 100 300 600]
    D = DM(DM.Condition==c & DM.tauP_ns==tp,:);
    pR = D.PCE_pct(D.Architecture=="RBSF-GaP"); pR = pR(1);
    pH = arrayfun(@(r) D.PCE_pct(D.Architecture=="RHJ" & abs(D.ratio-r)<1e-9), rList);
    Fc = F(string(F.Condition)==c,:);
    pN = arrayfun(@(r) Fc.PCE_pct(abs(Fc.tauN_ns - round(tp*r,6)) < 1e-6*tp*r), rList);
    rows(end+1,:) = {char(c), tp, pR, pH(1), pN(1), crossing(rList,pH-pR), crossing(rList,pN-pR)}; %#ok<AGROW>
  end
end
T = cell2table(rows,'VariableNames',{'Condition','tauP_ns','PCE_RBSF','PCE_RHJ_r1','PCE_nFJ_r1','rstar_RHJ','rstar_nFJ'});
writetable(T,fullfile(out,'rstar_ntype_fj.csv')); disp(T)
end

function r = crossing(rList,d)
k = find(d(1:end-1) < 0 & d(2:end) >= 0,1);
if d(1) >= 0, r = rList(1); return, end
if isempty(k), r = NaN; return, end
t = -d(k)/(d(k+1)-d(k)); r = exp(log(rList(k)) + t*(log(rList(k+1))-log(rList(k))));
end

function D = readDM(folder)
d = dir(fullfile(folder,'design_map_*.csv')); D = [];
for k = 1:numel(d), t = readtable(fullfile(folder,d(k).name)); t.Condition = string(t.Condition); t.Architecture = string(t.Architecture); D = [D; t]; end %#ok<AGROW>
end
