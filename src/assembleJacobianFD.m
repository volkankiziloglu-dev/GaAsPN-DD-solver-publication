function J = assembleJacobianFD(resFun, x0, F0, N, varsPerNode, stencil)
%ASSEMBLEJACOBIANFD Bandwidth-restricted sparse finite-difference Jacobian.
%
% resFun: function handle F = resFun(x) for a packed state vector x of
%   length N*varsPerNode (block-interleaved: x(varsPerNode*(i-1)+k), k=1..varsPerNode).
% x0, F0: baseline state/residual (F0 = resFun(x0), passed in so the
%   caller's already-computed baseline residual is reused, not recomputed).
% stencil: each unknown at node i is assumed to affect only residual rows
%   at nodes max(1,i-stencil):min(N,i+stencil) (default 1 -- nearest
%   neighbor only, matching this project's Scharfetter-Gummel discretization
%   and interfaceRecombination.m's junction-adjacent-pair coupling, neither
%   of which ever reaches beyond immediate neighbors).
%
% Central differences are used throughout (not one-sided): cost is not the
% bottleneck at this problem's scale (N~60-600), and central differencing
% removes an entire class of "sign bug vs truncation error" ambiguity
% during debugging -- see docs/MATHEMATICAL_MODEL.md N.16.6/N.17 for why
% an analytic Jacobian was deliberately NOT used for the initial
% implementation of this solver.
if nargin < 6 || isempty(stencil), stencil = 1; end
nTotal = N*varsPerNode;
rowsPerCol = varsPerNode*(2*stencil+1);
nnzEst = nTotal*rowsPerCol;
rowIdx = zeros(nnzEst,1); colIdx = zeros(nnzEst,1); val = zeros(nnzEst,1);
k = 0;

for i = 1:N
  jnLo = max(1,i-stencil); jnHi = min(N,i+stencil);
  rowList = (varsPerNode*(jnLo-1)+1 : varsPerNode*jnHi)';
  for v = 1:varsPerNode
    col = varsPerNode*(i-1)+v;
    epsFD = max(1e-7, 1e-6*max(1, abs(x0(col))));

    xp = x0; xp(col) = xp(col)+epsFD;
    xm = x0; xm(col) = xm(col)-epsFD;
    Fp = resFun(xp);
    Fm = resFun(xm);
    dF = (Fp(rowList) - Fm(rowList)) / (2*epsFD);

    nAdd = numel(rowList);
    rowIdx(k+1:k+nAdd) = rowList;
    colIdx(k+1:k+nAdd) = col;
    val(k+1:k+nAdd) = dF;
    k = k+nAdd;
  end
end
J = sparse(rowIdx(1:k), colIdx(1:k), val(1:k), nTotal, nTotal); %#ok<*NASGU>
end
