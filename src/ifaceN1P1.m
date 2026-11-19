function [n1, p1] = ifaceN1P1(N0, P0, Et)
%IFACEN1P1 SRH emission densities of the opt-in interface models for a trap at Et (eV) above the effective interface midgap:
% n1 = sqrt(N0*P0)*exp(Et/VT), p1 = sqrt(N0*P0)*exp(-Et/VT); Et = 0 gives n1 = p1 = sqrt(N0*P0).
% (Ported from Sb2S3_Sb2Se3_Paper/src, dd_solver v1.1.0.)
if nargin < 3 || isempty(Et), Et = 0; end
c = constants(); m = sqrt(N0*P0);
n1 = m*exp(Et/c.VT); p1 = m*exp(-Et/c.VT);
end
