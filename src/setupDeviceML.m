function dev = setupDeviceML(layers, absorberIdx, junctions, opts)
% Copy of RearHJT/setupDeviceML.m -- it calls gummelSolveML_v2
% (eps0-corrected).
c = constants();
dev = buildDeviceML(layers, absorberIdx, junctions, opts);

psiRough = solveEquilibrium(dev); % common/ -- ONLY an initial guess
chiRel = dev.chi-dev.chiRef;
nRough = dev.Nc .* exp((psiRough+chiRel)/c.VT);
pRough = dev.Nv .* exp(-(psiRough+chiRel+dev.Eg)/c.VT);

dev.psiEq = psiRough;
r = gummelSolveML(dev, 0.0, zeros(dev.N,1), struct('psiInit',psiRough,'nInit',nRough,'pInit',pRough));
dev.psiEq = r.psi;
dev.nEq = r.n;
dev.pEq = r.p;
end
