function analyze_uncertainty_sampling_ns()
%ANALYZE_UNCERTAINTY_SAMPLING_NS Ranking frequencies and Spearman correlations of the ns-regime paired sampling.
root = fileparts(fileparts(mfilename('fullpath')));
d = fullfile(root,'results','paper_uncertainty_ns_merged');
T = readtable(fullfile(d,'ns_samples_all.csv'));   % merged file: AM1.5G rows of the earlier run + indoor rows with the unified extractor
T.Condition = string(T.Condition); T.Architecture = string(T.Architecture); T.Class = string(T.Class);
fprintf('rows %d, samples %d\n', height(T), numel(unique(T.Sample)));
archs = ["HJ-0","F-HJ","RBSF-GaP","RHJ"];
for c = ["AM1.5G","6500K_1000lx"]
  fprintf('\n===== %s =====\n', c);
  S = unique(T.Sample)';
  M = nan(numel(S),4); Cl = strings(numel(S),4);
  for i = 1:numel(S)
    for a = 1:4
      r = T(T.Sample==S(i) & T.Condition==c & T.Architecture==archs(a),:);
      if ~isempty(r), M(i,a) = r.PCE_pct; Cl(i,a) = r.Class; end
    end
  end
  % valid = finite PCE and (AM1.5G or class PASS-A/B); SCREEN counted separately
  ok = ~isnan(M);
  % unified extractor: indoor uses all finite points (class is a label only)
  fprintf('valid samples per arch: '); for a = 1:4, fprintf('%s=%d ', archs(a), sum(ok(:,a))); end; fprintf('\n');
  pairs = {[3 4],[3 1],[4 1],[2 1],[3 2],[4 2]};
  for p = 1:numel(pairs)
    a = pairs{p}(1); b = pairs{p}(2); v = ok(:,a) & ok(:,b);
    dlt = M(v,a) - M(v,b);
    fprintf('%s > %s: %d/%d  (diff P5/P50/P95 = %.2f/%.2f/%.2f)\n', archs(a), archs(b), sum(dlt>0), sum(v), prctile(dlt,5), prctile(dlt,50), prctile(dlt,95));
  end
  % best architecture per sample among those valid
  allv = all(ok,2);
  [~,ib] = max(M(allv,:),[],2);
  fprintf('best architecture (samples with all four valid, n=%d): ', sum(allv));
  for a = 1:4, fprintf('%s=%d ', archs(a), sum(ib==a)); end; fprintf('\n');
  [~,iw] = min(M(allv,:),[],2);
  fprintf('worst architecture: '); for a = 1:4, fprintf('%s=%d ', archs(a), sum(iw==a)); end; fprintf('\n');
  % Spearman with parameters, for RBSF-RHJ difference and PCE of each arch
  Pn = {'tauP','ratio','Sback','Sint','Nabs','chiGaP','tauGaP','muN'};
  prm = nan(numel(S),numel(Pn));
  for i = 1:numel(S), r = T(T.Sample==S(i),:); for j = 1:numel(Pn), prm(i,j) = r.(Pn{j})(1); end; end
  v = ok(:,3) & ok(:,4);
  fprintf('Spearman rho, RBSF-RHJ difference (n=%d): ', sum(v));
  for j = 1:numel(Pn), fprintf('%s=%+.2f ', Pn{j}, corr(prm(v,j), M(v,3)-M(v,4), 'Type','Spearman')); end; fprintf('\n');
  for a = 1:4
    v = ok(:,a);
    fprintf('Spearman rho, PCE %s (n=%d): ', archs(a), sum(v));
    for j = 1:numel(Pn), fprintf('%s=%+.2f ', Pn{j}, corr(prm(v,j), M(v,a), 'Type','Spearman')); end; fprintf('\n');
  end
  fprintf('PCE P5/P50/P95: ');
  for a = 1:4, q = prctile(M(ok(:,a),a),[5 50 95]); fprintf('%s %.2f/%.2f/%.2f (n=%d)  ', archs(a), q, sum(ok(:,a))); end; fprintf('\n');
end
end
