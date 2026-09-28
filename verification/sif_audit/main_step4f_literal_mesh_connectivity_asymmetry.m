function Out=main_step4f_literal_mesh_connectivity_asymmetry(varargin)
%MAIN_STEP4F_LITERAL_MESH_CONNECTIVITY_ASYMMETRY
% A_conn experiment: break only lower-half triangulation while preserving
% the approved literal mesh T3 vertex coordinates exactly.
%
% This isolates topology/interpolation-stencil asymmetry from geometric
% asymmetry. Natural annular-cell diagonals are flipped in the strict lower
% half; upper connectivity and the crack neighborhood are frozen.
%
% Exact fields:
%   pure I      (1,0)
%   pure II     (0,1)
%   mixed 5%    (1,0.05)
%   mixed 1%    (1,0.01)
%
% Both the historical mirror/J method and FE-nodal EDI are compared against
% the prescribed exact SIFs.
%
% Usage:
%   O4F = main_step4f_literal_mesh_connectivity_asymmetry();

ip=inputParser;
addParameter(ip,'Fractions',[0 0.10 0.25 0.50 1.00], ...
    @(x)isnumeric(x)&&isvector(x)&&all(isfinite(x))&&all(x>=0)&&all(x<=1));
addParameter(ip,'OldNtheta',240,@(x)isnumeric(x)&&isscalar(x)&&x>=40);
addParameter(ip,'Plot',true,@(x)islogical(x)||isnumeric(x));
parse(ip,varargin{:}); opt=ip.Results;

here=fileparts(mfilename('fullpath'));
addpath(genpath(fileparts(fileparts(here))));

fprintf('\n============================================================\n');
fprintf('SIF AUDIT STEP 4F: CONNECTIVITY-ONLY ASYMMETRY (A_conn)\n');
fprintf('============================================================\n');

[baseMesh,baseInfo,parent]=build_literal_ring_crack_cut_mesh();
assert(baseInfo.auditT3.passed&&baseInfo.auditT6.passed, ...
    'step4f:GeometryGate','Approved S0 geometry gate must pass first.');

E=4e3; nu=0.30; ps=1;
coef=E/((1+nu)*(1-2*nu));
D=coef*[1-nu,nu,0;nu,1-nu,0;0,0,(1-2*nu)/2];
mat=struct('E',E,'nu',nu,'ps',ps,'D',D,'Dmat',D);
V=[-1,0;0,0];

Kcases=[1,0;0,1;1,0.05;1,0.01];
caseNames=["pure_I";"pure_II";"mixed_5pct";"mixed_1pct"];

oldR=[0.02 0.04 0.06 0.08 0.10 0.12 0.14];
rOldRef=0.08;

ediDomains=[0.02,0.10;0.02,0.12;0.03,0.14];
ediRef=2;

frac=opt.Fractions(:);
nF=numel(frac); nK=size(Kcases,1);

Meshes=cell(nF,1); Meta=cell(nF,1);
Old=cell(nF,nK); EDI=cell(nF,nK);

meshRows=nan(nF,11);
rows=nan(nF*nK,21);
rr=0;

for jf=1:nF
    [mesh,meta]=flip_literal_cut_mesh_lower_diagonals( ...
        baseMesh,baseInfo,frac(jf));

    Meshes{jf}=mesh; Meta{jf}=meta;

    meshRows(jf,:)=[ ...
        jf,frac(jf),meta.actualFraction,meta.nCandidateDiagonals, ...
        meta.nEligible,meta.nFlipped,meta.nChangedTriangles, ...
        meta.qualityMin,meta.qualityP05,meta.minAngleDeg, ...
        double(meta.geometryChanged)];

    fprintf('\n------------------------------------------------------------\n');
    fprintf('requested flip fraction = %.3f; actual = %.6f\n', ...
        frac(jf),meta.actualFraction);
    fprintf('candidate/disjoint-eligible/flipped lower diagonals = %d / %d / %d\n', ...
        meta.nCandidateDiagonals,meta.nEligible,meta.nFlipped);
    fprintf('changed T3 triangles = %d / %d\n', ...
        meta.nChangedTriangles,size(mesh.connect3,1));
    fprintf('quality min/p05 = %.5f / %.5f; min angle = %.3f deg\n', ...
        meta.qualityMin,meta.qualityP05,meta.minAngleDeg);

    % The defining isolation check for A_conn.
    assert(isequal(mesh.coord3,baseMesh.coord3), ...
        'step4f:GeometryChanged','A_conn must preserve every T3 vertex coordinate.');

    for ic=1:nK
        KItrue=Kcases(ic,1); KIItrue=Kcases(ic,2);

        U=exact_williams_displacement_audit( ...
            mesh.coord,KItrue,KIItrue,E,nu,ps, ...
            'UpperFaceIDs',mesh.crackUpperT6IDs, ...
            'LowerFaceIDs',mesh.crackLowerT6IDs);

        % Historical mirror/J contour sweep.
        A=nan(numel(oldR),10);
        Dold=cell(numel(oldR),1);
        for ir=1:numel(oldR)
            [KIo,KIIo,DO]=SIF_LEFM_circle2_debug( ...
                mesh,U,V,mat,oldR(ir), ...
                'nthet',opt.OldNtheta,'plot',false,'verbose',false);

            detDen=0.5*(abs(DO.detJP)+abs(DO.detJQ));
            good=detDen>0;
            detMismatch=NaN;
            if any(good)
                detMismatch=median(abs(DO.detJP(good)-DO.detJQ(good))./detDen(good));
            end
            baryMismatch=median(abs(DO.baryMinP-DO.baryMinQ));

            A(ir,:)=[ ...
                oldR(ir),KIo,KIIo, ...
                KIo-KItrue,KIIo-KIItrue, ...
                hypot(KIo-KItrue,KIIo-KIItrue), ...
                DO.JI_over_absint,DO.JII_over_absint, ...
                detMismatch,baryMismatch];
            Dold{ir}=DO;
        end

        Told=array2table(A,'VariableNames',{ ...
            'rI','KI_old','KII_old','KI_error','KII_error', ...
            'vector_error_abs','JI_over_absint','JII_over_absint', ...
            'median_detJ_PQ_mismatch','median_baryMin_PQ_mismatch'});
        Old{jf,ic}=struct('table',Told,'debug',{Dold});

        % FE-nodal EDI domain check.
        B=nan(size(ediDomains,1),7);
        Dedi=cell(size(ediDomains,1),1);
        for id=1:size(ediDomains,1)
            dom=struct('r_inner',ediDomains(id,1), ...
                       'r_outer',ediDomains(id,2));
            [KIe,KIIe,AE]=SIF_LEFM_interaction_EDI( ...
                mesh,U,V,mat,dom, ...
                'UsePlaneStrain',true,'Verbose',false, ...
                'WeightFunction','fe_nodal');

            B(id,:)=[ ...
                dom.r_inner,dom.r_outer,KIe,KIIe, ...
                KIe-KItrue,KIIe-KIItrue, ...
                hypot(KIe-KItrue,KIIe-KIItrue)];
            Dedi{id}=AE;
        end

        Tedi=array2table(B,'VariableNames',{ ...
            'r_inner','r_outer','KI_EDI','KII_EDI', ...
            'KI_error','KII_error','vector_error_abs'});
        EDI{jf,ic}=struct('table',Tedi,'debug',{Dedi});

        [~,io]=min(abs(oldR-rOldRef));
        KIo=Told.KI_old(io); KIIo=Told.KII_old(io);
        KIe=Tedi.KI_EDI(ediRef); KIIe=Tedi.KII_EDI(ediRef);

        trueScale=max(hypot(KItrue,KIItrue),1);
        if abs(KIItrue)>0
            relKIIold=(KIIo-KIItrue)/abs(KIItrue);
            relKIIedi=(KIIe-KIItrue)/abs(KIItrue);
        else
            relKIIold=NaN; relKIIedi=NaN;
        end

        rr=rr+1;
        rows(rr,:)=[ ...
            jf,ic,meta.actualFraction,KItrue,KIItrue, ...
            KIo,KIIo,KIe,KIIe, ...
            KIo-KItrue,KIIo-KIItrue,KIe-KItrue,KIIe-KIItrue, ...
            hypot(KIo-KItrue,KIIo-KIItrue)/trueScale, ...
            hypot(KIe-KItrue,KIIe-KIItrue)/trueScale, ...
            relKIIold,relKIIedi, ...
            max(Told.KI_old)-min(Told.KI_old), ...
            max(Told.KII_old)-min(Told.KII_old), ...
            max(Tedi.KI_EDI)-min(Tedi.KI_EDI), ...
            max(Tedi.KII_EDI)-min(Tedi.KII_EDI)];
    end
end

Tmesh=array2table(meshRows,'VariableNames',{ ...
    'meshID','requested_flip_fraction','actual_flip_fraction', ...
    'n_candidate_diagonals','n_eligible_diagonals','n_flipped_diagonals', ...
    'n_changed_triangles','quality_min','quality_p05','min_angle_deg', ...
    'geometry_changed'});

T=array2table(rows,'VariableNames',{ ...
    'meshID','caseID','flip_fraction','KI_true','KII_true', ...
    'KI_old','KII_old','KI_EDI','KII_EDI', ...
    'KI_old_error','KII_old_error','KI_EDI_error','KII_EDI_error', ...
    'old_vector_error_rel','EDI_vector_error_rel', ...
    'relative_KII_error_old','relative_KII_error_EDI', ...
    'old_KI_radius_range','old_KII_radius_range', ...
    'EDI_KI_domain_range','EDI_KII_domain_range'});
T.caseName=caseNames(T.caseID);
T=movevars(T,'caseName','After','caseID');

detective=nan(nF,13);
for jf=1:nF
    Told=Old{jf,1}.table;
    [~,io]=min(abs(Told.rI-rOldRef));

    pureI=T(T.meshID==jf & T.caseID==1,:);
    pureII=T(T.meshID==jf & T.caseID==2,:);
    mix1=T(T.meshID==jf & T.caseID==4,:);

    detective(jf,:)=[ ...
        jf,Tmesh.actual_flip_fraction(jf),Tmesh.n_flipped_diagonals(jf), ...
        Told.median_detJ_PQ_mismatch(io), ...
        Told.median_baryMin_PQ_mismatch(io), ...
        pureI.KII_old,pureII.KI_old,mix1.relative_KII_error_old, ...
        pureI.KII_EDI,pureII.KI_EDI,mix1.relative_KII_error_EDI, ...
        Tmesh.quality_min(jf),Tmesh.min_angle_deg(jf)];
end

Tdetective=array2table(detective,'VariableNames',{ ...
    'meshID','flip_fraction','n_flipped_diagonals', ...
    'old_detJ_PQ_mismatch','old_bary_PQ_mismatch', ...
    'false_KII_old_from_pureI','false_KI_old_from_pureII', ...
    'mixed1pct_relative_KII_error_old', ...
    'false_KII_EDI_from_pureI','false_KI_EDI_from_pureII', ...
    'mixed1pct_relative_KII_error_EDI', ...
    'quality_min','min_angle_deg'});

fprintf('\n============================================================\n');
fprintf('CONNECTIVITY-ONLY MESH SUMMARY\n');
fprintf('All T3 coordinates are IDENTICAL to S0.\n');
fprintf('============================================================\n');
disp(Tmesh);

fprintf('\n============================================================\n');
fprintf('DETECTIVE TABLE: CONNECTIVITY ASYMMETRY -> MODAL CONTAMINATION\n');
fprintf('old/J reference rI=%.4f; EDI reference annulus=[%.4f,%.4f]\n', ...
    rOldRef,ediDomains(ediRef,1),ediDomains(ediRef,2));
fprintf('============================================================\n');
disp(Tdetective);

fprintf('\nFull exact-field summary\n');
disp(T);

if logical(opt.Plot)
    plot_detective(Tdetective);
end

Out=struct();
Out.meshSummary=Tmesh;
Out.detective=Tdetective;
Out.results=T;
Out.meshes=Meshes;
Out.meta=Meta;
Out.old=Old;
Out.edi=EDI;
Out.baseMesh=baseMesh;
Out.baseInfo=baseInfo;
Out.parent=parent;
Out.settings=struct('fractions',frac,'oldR',oldR,'rOldRef',rOldRef, ...
    'ediDomains',ediDomains,'ediRef',ediRef, ...
    'Kcases',Kcases,'caseNames',caseNames,'OldNtheta',opt.OldNtheta);

fprintf('\nSTEP 4F completed.\n');
fprintf(['This is the topology-only interrogation: geometry is exactly mirror ', ...
    'symmetric. Any new old-method cross-mode leakage is therefore caused ', ...
    'by loss of mirrored interpolation connectivity, not node positions.\n']);
end


function plot_detective(T)

x=T.flip_fraction;

figure('Name','Step 4F: topology-only false cross-mode SIFs','NumberTitle','off');
plot(x,abs(T.false_KII_old_from_pureI),'o-'); hold on;
plot(x,abs(T.false_KI_old_from_pureII),'s-');
plot(x,abs(T.false_KII_EDI_from_pureI),'o--');
plot(x,abs(T.false_KI_EDI_from_pureII),'s--');
grid on;
xlabel('fraction of eligible lower diagonals flipped');
ylabel('absolute false cross-mode SIF');
legend('|K_{II}| old from pure I','|K_I| old from pure II', ...
    '|K_{II}| EDI from pure I','|K_I| EDI from pure II', ...
    'Location','best');
title('Connectivity asymmetry versus modal contamination');

figure('Name','Step 4F: P-Q interpolation diagnostics','NumberTitle','off');
semilogy(x,max(T.old_detJ_PQ_mismatch,eps),'o-'); hold on;
semilogy(x,max(T.old_bary_PQ_mismatch,eps),'s-');
grid on;
xlabel('fraction of eligible lower diagonals flipped');
ylabel('P/Q mismatch diagnostic');
legend('median detJ mismatch','median barycentric mismatch','Location','best');
title('Topology-only loss of mirrored interpolation correspondence');

figure('Name','Step 4F: 1-percent KII error','NumberTitle','off');
semilogy(x,max(abs(T.mixed1pct_relative_KII_error_old),eps),'o-'); hold on;
semilogy(x,max(abs(T.mixed1pct_relative_KII_error_EDI),eps),'s-');
grid on;
xlabel('fraction of eligible lower diagonals flipped');
ylabel('relative error in K_{II}');
legend('old mirror/J','FE-nodal EDI','Location','best');
title('K_{II}/K_I=0.01 under connectivity-only asymmetry');
end
