function Out = main_step3_radius_domain_sweep()
%MAIN_STEP3_RADIUS_DOMAIN_SWEEP
% Third executable gate of the SIF asymmetric-mesh audit.
%
% The FEM field is solved ONCE, then both SIF extractors are swept over
% matched outer radii. EDI is additionally swept over inner-domain radius
% at fixed outer radius.
%
% No mesh-asymmetry perturbation is introduced in this step.

    here = fileparts(mfilename('fullpath'));
    repoRoot = fileparts(fileparts(here));
    addpath(genpath(repoRoot));

    fprintf('\n============================================================\n');
    fprintf('SIF AUDIT STEP 3: RADIUS / DOMAIN SENSITIVITY\n');
    fprintf('============================================================\n');

    C = cfg_crack_path_two_leg_control(2.0, ...
        'plotGeom',false,'plotMesh',false);

    Out = sweep_same_field_radii(C,C.sigma0, ...
        'radiusFractions',[0.2 0.3 0.4 0.5 0.6], ...
        'innerFactor',0.1, ...
        'innerFactorsSweep',[0.05 0.1 0.2 0.3], ...
        'fixedOuterFraction',0.5, ...
        'nthet',100, ...
        'Verbose',true);

    fprintf('\nSTEP 3 completed.\n');
    fprintf(['Interpret radius trends before changing the mesh: ', ...
        'this gate separates extraction-domain sensitivity from ', ...
        'mesh-asymmetry sensitivity.\n']);
end
