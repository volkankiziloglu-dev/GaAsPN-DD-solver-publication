function assertIfaceModelSupported(dev, caller)
%ASSERTIFACEMODELSUPPORTED Stop when a residual/solver variant without carrier-specific interface sinks is used with 'four'/'cross'.
% Only assembleResidualQFP (gummelSolveQFPAnalytic, gummelSolveQFP, solveCoupledQFP) implements these two models; the other
% variants would silently use a different sink placement than the production residual (dd_solver v1.1.1).
if ifaceCarrierSpecific(dev)
  error('dd_solver:ifaceModelUnsupported', ['%s does not support the interface model ''%s'' (carrier-specific sinks). Use ' ...
    'gummelSolveQFPAnalytic, gummelSolveQFP or solveCoupledQFP, or the models ''geometric''/''onesided''.'], caller, dev.ifaceModel);
end
end
