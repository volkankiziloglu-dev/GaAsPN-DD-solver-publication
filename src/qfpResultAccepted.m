function accepted = qfpResultAccepted(r)
%QFPRESULTACCEPTED Single canonical accept/reject gate for a gummelSolveQFP
% output struct. Any script that turns r.Jseg into a Jsc/Voc/FF/EQE number
% MUST filter through this first -- checking r.isFinite alone (the old
% pattern in every case/test script) accepts points where the electron or
% hole continuity block never actually converged (see N.16.13 in
% MATHEMATICAL_MODEL.md): r.isFinite only means "no NaN/Inf appeared", not
% "the physics is self-consistent". This is what MD421's V=0 point was
% doing -- isFinite=true, physicsPass=false, silently averaged into Jsc
% anyway via the back-region estimator.
accepted = r.isFinite && r.converged && r.physicsPass ...
    && r.poissonConverged ...
    && r.currentConservation.pass ...
    && r.continuityResiduals.electron.pass ...
    && r.continuityResiduals.hole.pass;
end
