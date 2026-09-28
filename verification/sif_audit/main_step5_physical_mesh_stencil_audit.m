function Out=main_step5_physical_mesh_stencil_audit(varargin)
%MAIN_STEP5_PHYSICAL_MESH_STENCIL_AUDIT
% Bridge the synthetic topology mechanism to the actual two-leg FEM control.
%
% One physical FEM field is solved once. At several circular old/J contour
% radii we report:
%   - old KI,KII,
%   - FE-nodal EDI KI,KII on matched outer domains,
%   - old-vs-EDI vector difference,
%   - direct reflected-T3 stencil mismatch statistics,
%   - fraction of exactly mirrored P/Q parent triangles.
%
% This does not assume either extractor is exact truth; it asks whether the
% same discrete-stencil mechanism identified in Steps 4F--4G is present in
% the real Crack-Path-style mesh.

ip=inputParser;
addParameter(ip,'radiusFractions',[0.2 0.3 0.4 0.5 0.6], ...
    @(x)isnumeric(x)&&isvector(x)&&all(isfinite(x))&&all(x>0));
addParameter(ip,'innerFactor',0.1, ...
    @(x)isnumeric(x)&&isscalar(x)&&x>0&&x<1);
addParameter(ip,'nthet',240, ...
    @(x)isnumeric(x)&&isscalar(x)&&x>=40);
parse(ip,varargin{:}); opt=ip.Results;

here=fileparts(mfilename('fullpath'));
addpath(genpath(fileparts(fileparts(here))));

fprintf('\n============================================================\n');
fprintf('SIF AUDIT STEP 5: PHYSICAL-MESH STENCIL AUDIT\n');
fprintf('============================================================\n');

C=cfg_crack_path_two_leg_control(2.0,'plotGeom',false,'plotMesh',false);
Sol=solve_crack_path_polyline_field(C,C.sigma0);

V=Sol.V;
if size(V,1)<2
    error('step5:CrackPolyline','Need at least two crack-polyline points.');
end
lastLeg=norm(V(end,:)-V(end-1,:));

rf=opt.radiusFractions(:);
n=numel(rf);
R=nan(n,18);
DebugOld=cell(n,1);
DebugEDI=cell(n,1);

for i=1:n
    rI=rf(i)*lastLeg;

    [KIo,KIIo,DO]=SIF_LEFM_circle2_debug( ...
        Sol.mesh,Sol.U,V,Sol.mat,rI, ...
        'nthet',opt.nthet,'plot',false,'verbose',false);

    domain=struct('r_inner',opt.innerFactor*rI,'r_outer',rI);
    [KIe,KIIe,AE]=SIF_LEFM_interaction_EDI( ...
        Sol.mesh,Sol.U,V,Sol.mat,domain, ...
        'UsePlaneStrain',Sol.mat.ps==1,'Verbose',false, ...
        'WeightFunction','fe_nodal');

    detDen=0.5*(abs(DO.detJP)+abs(DO.detJQ));
    good=detDen>0;
    detMismatch=NaN;
    if any(good)
        detMismatch=median(abs(DO.detJP(good)-DO.detJQ(good))./detDen(good));
    end
    baryMismatch=median(abs(DO.baryMinP-DO.baryMinQ));

    S=DO.diagnostics;
    oldVec=hypot(KIo,KIIo);
    dVec=hypot(KIo-KIe,KIIo-KIIe)/max(oldVec,eps);

    R(i,:)=[ ...
        rf(i),rI,KIo,KIIo,KIe,KIIe, ...
        KIo-KIe,KIIo-KIIe,dVec, ...
        detMismatch,baryMismatch, ...
        S.mirrorT3Mismatch_median,S.mirrorT3Mismatch_p95, ...
        S.mirrorT3Mismatch_max,S.frac_exact_mirror_T3, ...
        S.min_baryP,S.min_baryQ,S.frac_same_elem_PQ];

    DebugOld{i}=DO;
    DebugEDI{i}=AE;
end

T=array2table(R,'VariableNames',{ ...
    'r_over_lastLeg','rI','KI_old','KII_old','KI_EDI','KII_EDI', ...
    'dKI_old_minus_EDI','dKII_old_minus_EDI','vector_difference_rel', ...
    'median_detJ_PQ_mismatch','median_bary_PQ_mismatch', ...
    'mirror_T3_mismatch_median','mirror_T3_mismatch_p95', ...
    'mirror_T3_mismatch_max','fraction_exact_mirror_T3', ...
    'min_baryP','min_baryQ','frac_same_elem_PQ'});

fprintf('\nPHYSICAL-FIELD STENCIL TABLE\n');
disp(T);

% Descriptive association across contour radii only; do not interpret as a
% universal convergence law because changing radius also changes the sampled
% field and EDI domain.
broken=1-T.fraction_exact_mirror_T3;
if numel(broken)>=3 && std(broken)>0
    Ccorr=corrcoef(broken,abs(T.dKII_old_minus_EDI));
    rho=Ccorr(1,2);
else
    rho=NaN;
end

fprintf('Descriptive corr(1-f_exact, |dKII old-EDI|) across radii = %.6f\n',rho);

Out=struct('table',T,'solution',Sol,'caseConfig',C, ...
    'debugOld',{DebugOld},'debugEDI',{DebugEDI}, ...
    'correlationBrokenPairVsAbsDKII',rho, ...
    'settings',opt);

fprintf('\nSTEP 5 completed.\n');
fprintf(['Interpret this only as a bridge to the physical mesh: exact truth is ', ...
    'not known here. The synthetic Steps 4D--4G remain the causal proof.\n']);
end
