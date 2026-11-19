function ptsTerminal = applyTerminalResistance(ptsDevice, Rs_ohmcm2, Rsh_ohmcm2)
%APPLYTERMINALRESISTANCE Map intrinsic device J-V to terminal J-V.
%
% Sign convention follows pvMetrics: generated current is positive.
%   Vdevice = Vterminal + Jterminal*Rs
%   Jterminal = Jdevice(Vdevice) - Vterminal/Rsh
%
% ptsDevice is an Nx2 [Vdevice, Jdevice(A/cm2)] curve. Rs and Rsh are
% area-normalized resistances in ohm*cm2. The transformation is evaluated
% parametrically at every intrinsic curve point and sorted by Vterminal.
arguments
  ptsDevice (:,2) double
  Rs_ohmcm2 (1,1) double {mustBeNonnegative}
  Rsh_ohmcm2 (1,1) double {mustBePositive}
end
ratio = Rs_ohmcm2/Rsh_ohmcm2;
if ratio >= 1
  error('applyTerminalResistance:InvalidNetwork', ...
    'Rs must be smaller than Rsh for a passive, single-valued terminal map.');
end
Vd = ptsDevice(:,1);
Jd = ptsDevice(:,2);
Jt = (Jd - Vd/Rsh_ohmcm2) ./ (1-ratio);
Vt = Vd - Jt*Rs_ohmcm2;
ptsTerminal = sortrows([Vt Jt],1);
if any(~isfinite(ptsTerminal),'all')
  error('applyTerminalResistance:Nonfinite','Terminal map produced NaN/Inf.');
end
end
