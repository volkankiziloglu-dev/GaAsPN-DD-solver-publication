function J = assembleJacobianFDColored(resFun, x0, N, varsPerNode, stencil)
%ASSEMBLEJACOBIANFDCOLORED Central-difference Jacobian of a node-local residual by column colouring (2026-11-07).
% For a residual whose rows at node i depend only on nodes i-stencil..i+stencil, the columns of nodes i, i+(2*stencil+1), ... do not
% share rows, so they are perturbed together: 2*varsPerNode*(2*stencil+1) residual evaluations instead of 2*varsPerNode*N
% (assembleJacobianFD.m). Step as in assembleJacobianFD.m (max(1e-7, 1e-6*max(1,|x|))). Returns a sparse matrix.
if nargin < 5 || isempty(stencil), stencil = 1; end
nTotal = N*varsPerNode; P = 2*stencil + 1;
rowsAll = []; colsAll = []; valsAll = [];
for grp = 1:P
  nodes = grp:P:N;
  for v = 1:varsPerNode
    cols = varsPerNode*(nodes-1) + v;
    epsFD = max(1e-7, 1e-6*max(1, abs(x0(cols))));
    xp = x0; xp(cols) = xp(cols) + epsFD;
    xm = x0; xm(cols) = xm(cols) - epsFD;
    dF = resFun(xp) - resFun(xm);
    for k = 1:numel(nodes)
      i = nodes(k);
      jnLo = max(1, i-stencil); jnHi = min(N, i+stencil);
      rows = (varsPerNode*(jnLo-1)+1 : varsPerNode*jnHi)';
      rowsAll = [rowsAll; rows]; %#ok<AGROW>
      colsAll = [colsAll; repmat(cols(k), numel(rows), 1)]; %#ok<AGROW>
      valsAll = [valsAll; dF(rows)/(2*epsFD(k))]; %#ok<AGROW>
    end
  end
end
J = sparse(rowsAll, colsAll, valsAll, nTotal, nTotal);
end
