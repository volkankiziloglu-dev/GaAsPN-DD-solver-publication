function p = gapParams()
% GaP (gallium phosphide) material parameters -- for the n+ window and p+ rear
% emitter layers of the 4-layer target structure (see Task #9).
%
% IMPORTANT: the lesson learned in the a-Si:H/c-Si calibration (Si/run_asi_csi_hjt_validate.m,
% MATHEMATICAL_MODEL.md Sec. I.2) applies here as well --
% literature values (especially the electron affinity chi and the mobility under heavy
% doping) can vary SIGNIFICANTLY depending on the reference source. Below, each
% field is EXPLICITLY labeled as LITERATURE (direct measurement/table value, with source) or
% ESTIMATED/CALIBRATION-REQUIRED (no clear source, trial-and-error/empirical tuning
% may be needed). The first run in Task #9
% should therefore be treated not as a VALIDATION but as a CONSISTENCY/EXPLORATION run
% (no comparison data exist in the literature for this EXACT 4-layer
% structure).

% --- LITERATURE (Ioffe NSM Archive, ioffe.ru/SVA/NSM/Semicond/GaP/,
%     compilation based on Levinshtein/Rumyantsev/Shur "Handbook Series on Semiconductor
%     Parameters"; 300 K) ---
p.Eg = 2.26;              % eV, fundamental (indirect, X-valley) band gap, 300K.
                           % NOTE: the direct (Gamma-Gamma) transition is at 2.78 eV --
                           % for the optical absorption edge (since it is a thin/transparent window
                           % layer, absorption is neglected in this project by assuming
                           % alphaMax=0, see below).
p.EgT_a = 2.34;            % Eg(T) = EgT_a - EgT_b*T^2/(T+EgT_c) [eV], 0-1200K
p.EgT_b = 6.0e-4;
p.EgT_c = 460;             % K
p.eps = 11.1;              % static dielectric constant (9-11.1 in some
                           % sources; not to be confused with the high-frequency
                           % eps_inf=9.11 -- the static/DC value is used here,
                           % which is the correct one for the Poisson equation).
p.NcNv = 1.85e19;          % cm^-3, sqrt(Nc*Nv), Nc=1.8e19 Nv=1.9e19 (300K,
                           % Ioffe). Used as a single value per this project's convention (see buildDevice.m
                           % NcNv_ETL/NcNv_Abs).
p.muN_max = 250;           % cm^2/Vs, UNDOPED/lightly doped upper limit (300K).
p.muP_max = 150;           % cm^2/Vs, likewise undoped upper limit.

% --- LITERATURE (mobility under heavy doping, UPDATED -- Task "parameter
%     hunt": Ma & Lundstrom-type ab-initio ionized-impurity scattering
%     calculation, arXiv:2301.02323, cross-validated against measured values) ---
p.muN_doped = 100;         % cm^2/Vs, for Nd~1e19 cm^-3. Source: at 1e18
                           % calc=157, meas~100 cm2/Vs (calculation systematically
                           % higher than measurement); at 1e19 no explicit electron
                           % value is given, but for holes a ~%24 drop is seen
                           % in the 1e18->1e19 transition (124->94) --
                           % assuming a similar relative drop for electrons + the meas-calc
                           % offset, 100 cm2/Vs was CHOSEN
                           % (old value was 20 cm2/Vs -- UNDERESTIMATED BY
                           % 5x).
p.muP_doped = 80;          % cm^2/Vs, for Na~1e19 cm^-3. The source gives it
                           % DIRECTLY: calc=94 cm2/Vs @ 1e19 (at 1e18
                           % calc=124, meas=90 -- again calc>meas trend);
                           % taking this offset into account, 80 cm2/Vs was CHOSEN
                           % (old value was 10 cm2/Vs -- UNDERESTIMATED BY
                           % 8-9x). This has NOT YET been applied in the RearHJT/run_4layer_*.m
                           % scripts -- the calling code reads
                           % gap.muN_doped/gap.muP_doped automatically,
                           % so it will be updated AUTOMATICALLY in the next
                           % run.

% --- ESTIMATED / CALIBRATION-REQUIRED (NO clear, single source value) ---
p.chi = 3.8;               % eV, electron affinity -- LARGE UNCERTAINTY: different
                           % sources give values in the 3.8-4.3 eV range
                           % (classic Sze-type table ~4.3 eV; some recent
                           % compilations report ~3.8 eV). Since Anderson's rule
                           % is approximate anyway (see buildDevice.m /
                           % buildDeviceML.m comments, MATHEMATICAL_MODEL.md
                           % Sec. A.2), as in the a-Si:H/c-Si calibration,
                           % this value should be treated as a starting point that will likely be
                           % TUNED to the experimental Voc/band-offset behavior,
                           % not as an exact physical constant.
                           %
                           % TASK #11 FINDING (see RearHJT/gummelSolveML.m
                           % investigation): with 4.3 (deltaEc=+0.3eV, a BLOCKING
                           % spike for electrons at the n+GaP window/
                           % n-GaAsNP absorber interface), under this model's NON-
                           % DEGENERATE Boltzmann statistics (n=ni*exp(psiEff/VT),
                           % unbounded from above) the local electron pile-up under forward bias
                           % grew WITHOUT BOUND and destabilized the Gummel loop
                           % (confirmed: n~5.8e23 cm^-3,
                           % ~30000x NcNv=1.85e19). With 3.8 (deltaEc=-0.2eV,
                           % NOT BLOCKING) the same structure gives the FULL J(V) curve
                           % without problems. Therefore 3.8 was CHOSEN here as the
                           % default -- both cited in the literature and the
                           % physically/electrically PREFERRED direction for electron transport
                           % in this architecture (n-n isotype front interface).
                           % If a truly blocking interface is needed (some
                           % real devices may have deliberate spike/notch
                           % designs), this model's non-degenerate statistics
                           % will be UNRELIABLE in that regime -- see Task #11.
                           %
                           % NOTE: this "blocking" issue is critical ONLY at majority-
                           % carrier (isotype, n-n or p-p) interfaces.
                           % At a MINORITY-carrier interface (e.g.
                           % at the absorber/rear-emitter boundary, which is a real p-n
                           % junction) a barrier in the same direction is a BSF-like
                           % BENEFICIAL effect and does not trigger the same
                           % instability (confirmed, see RearHJT/
                           % run_4layer_explore.m: while this p.chi=3.8 is used for the window,
                           % DELIBERATELY using chi=4.3 for the rear emitter
                           % gives the full J(V) curve).
                           % So p.chi is NOT a SINGLE "GaP value", but a
                           % starting point that should be CHOSEN according to
                           % position/role.
p.tau = 1e-9;              % s, SRH lifetime -- ESTIMATE; in the thin (20-150nm) window/emitter
                           % layer most carriers cross the layer anyway before being
                           % collected, so this value is not critical (low
                           % sensitivity expected, but not verified).

% --- OPTICAL (window/emitter assumed TRANSPARENT in this project) ---
p.alphaMax = 0;            % cm^-1, default: the thin GaP window/emitter layer
                           % is taken as OPTICALLY TRANSPARENT (alphaMax=0),
                           % similar to the alphaETLmax=0 default in common/buildDevice.m
                           % -- real GaP absorbs weakly anyway because of its indirect
                           % band transition (see the 2.78 eV direct transition in the Eg
                           % note), but without a full n(E),k(E)
                           % table (planned for RearHJT/data/) this
                           % simplification is used. If a clear need arises in Task #9,
                           % a real optical table should be added under
                           % RearHJT/data/.
end
