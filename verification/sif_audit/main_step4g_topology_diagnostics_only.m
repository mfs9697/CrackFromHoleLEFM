function Out=main_step4g_topology_diagnostics_only(varargin)
%MAIN_STEP4G_TOPOLOGY_DIAGNOSTICS_ONLY
% Lightweight follow-up to Step 4F. Reuses the same A_conn family but runs
% only one exact pure-I old/J extraction at rI=0.08 per mesh to populate the
% new direct mirrored-T3 stencil diagnostics.
%
% No EDI sweep and no multi-radius sweep are repeated.

ip=inputParser;
addParameter(ip,'Fractions',[0 0.10 0.25 0.50 1.00], ...
    @(x)isnumeric(x)&&isvector(x)&&all(isfinite(x))&&all(x>=0)&&all(x<=1));
addParameter(ip,'OldNtheta',240,@(x)isnumeric(x)&&isscalar(x)&&x>=40);
parse(ip,varargin{:}); opt=ip.Results;

here=fileparts(mfilename('fullpath'));
addpath(genpath(fileparts(fileparts(here))));

[baseMesh,baseInfo]=build_literal_ring_crack_cut_mesh();

E=4e3; nu=0.30; ps=1;
coef=E/((1+nu)*(1-2*nu));
D=coef*[1-nu,nu,0;nu,1-nu,0;0,0,(1-2*nu)/2];
mat=struct('E',E,'nu',nu,'ps',ps,'D',D,'Dmat',D);
V=[-1,0;0,0];

frac=opt.Fractions(:);
R=nan(numel(frac),11);

fprintf('\n============================================================\n');
fprintf('SIF AUDIT STEP 4G: DIRECT MIRRORED-STENCIL DIAGNOSTICS\n');
fprintf('============================================================\n');

for i=1:numel(frac)
    [mesh,meta]=flip_literal_cut_mesh_lower_diagonals( ...
        baseMesh,baseInfo,frac(i));

    U=exact_williams_displacement_audit( ...
        mesh.coord,1,0,E,nu,ps, ...
        'UpperFaceIDs',mesh.crackUpperT6IDs, ...
        'LowerFaceIDs',mesh.crackLowerT6IDs);

    [KI,KII,Dg]=SIF_LEFM_circle2_debug( ...
        mesh,U,V,mat,0.08, ...
        'nthet',opt.OldNtheta,'plot',false,'verbose',false);

    detDen=0.5*(abs(Dg.detJP)+abs(Dg.detJQ));
    good=detDen>0;
    detMismatch=median(abs(Dg.detJP(good)-Dg.detJQ(good))./detDen(good));
    baryMismatch=median(abs(Dg.baryMinP-Dg.baryMinQ));

    S=Dg.diagnostics;

    R(i,:)=[ ...
        frac(i),meta.actualFraction,meta.nFlipped, ...
        KI,KII,detMismatch,baryMismatch, ...
        S.mirrorT3Mismatch_median,S.mirrorT3Mismatch_p95, ...
        S.mirrorT3Mismatch_max,S.frac_exact_mirror_T3];
end

T=array2table(R,'VariableNames',{ ...
    'requested_flip_fraction','actual_flip_fraction','n_flipped', ...
    'KI_old','false_KII_old', ...
    'median_detJ_PQ_mismatch','median_bary_PQ_mismatch', ...
    'mirror_T3_mismatch_median','mirror_T3_mismatch_p95', ...
    'mirror_T3_mismatch_max','fraction_exact_mirror_T3'});

disp(T);

Out=struct('table',T,'settings',opt);
end
