function Out=main_step4d_literal_mesh_sif_validation(varargin)
%MAIN_STEP4D_LITERAL_MESH_SIF_VALIDATION
% Exact-SIF validation on the approved literal ring lattice after insertion
% of the validated negative-x crack cut.
%
% This is the first SIF calculation accepted on the approved mesh geometry.
% The parent lattice and crack cut are not altered here.
%
% Exact fields:
%   pure_I        (KI,KII) = (1,0)
%   pure_II       (0,1)
%   mixed_5pct    (1,0.05)
%   mixed_1pct    (1,0.01)   representative of the small-mode-II regime
%
% Historical mirror/J extraction:
%   circular contours rI = [0.02 0.04 0.06 0.08 0.10 0.12 0.14].
%
% Canonical interaction EDI:
%   FE-nodal q; several annuli are checked for domain stability.
%
% Both methods are compared directly with the prescribed exact SIFs.
% Neither extractor is treated as the reference truth.
%
% Usage:
%   O4D = main_step4d_literal_mesh_sif_validation();

ip=inputParser;
addParameter(ip,'OldNtheta',240,@(x)isnumeric(x)&&isscalar(x)&&x>=40);
addParameter(ip,'Plot',true,@(x)islogical(x)||isnumeric(x));
parse(ip,varargin{:}); opt=ip.Results;

here=fileparts(mfilename('fullpath'));
addpath(genpath(fileparts(fileparts(here))));

fprintf('\n============================================================\n');
fprintf('SIF AUDIT STEP 4D: APPROVED LITERAL-MESH SIF VALIDATION\n');
fprintf('============================================================\n');

[mesh,ginfo,parent]=build_literal_ring_crack_cut_mesh();
assert(ginfo.auditT3.passed && ginfo.auditT6.passed, ...
    'step4d:GeometryGate','Literal crack-cut geometry gate must pass.');

E=4e3; nu=0.30; ps=1;
coef=E/((1+nu)*(1-2*nu));
D=coef*[1-nu,nu,0;nu,1-nu,0;0,0,(1-2*nu)/2];
mat=struct('E',E,'nu',nu,'ps',ps,'D',D,'Dmat',D);

% Mathematical crack tip at the origin; the mesh begins at r0=0.005.
V=[-1,0;0,0];

Kcases=[ ...
    1,0; ...
    0,1; ...
    1,0.05; ...
    1,0.01];
caseNames=["pure_I";"pure_II";"mixed_5pct";"mixed_1pct"];

oldR=[0.02 0.04 0.06 0.08 0.10 0.12 0.14];
rOldRef=0.08;

ediDomains=[ ...
    0.015,0.08; ...
    0.020,0.10; ...
    0.020,0.12; ...
    0.030,0.14];
ediRef=3; % [0.02,0.12]

nK=size(Kcases,1);
Old=cell(nK,1);
EDI=cell(nK,1);
rows=nan(nK,20);

for ic=1:nK
    KItrue=Kcases(ic,1);
    KIItrue=Kcases(ic,2);

    U=exact_williams_displacement_audit( ...
        mesh.coord,KItrue,KIItrue,E,nu,ps, ...
        'UpperFaceIDs',mesh.crackUpperT6IDs, ...
        'LowerFaceIDs',mesh.crackLowerT6IDs);

    % ------------------------------------------------------------
    % Historical circular mirror/J sweep
    % ------------------------------------------------------------
    A=nan(numel(oldR),10);
    DebugOld=cell(numel(oldR),1);

    for ir=1:numel(oldR)
        [KIold,KIIold,DO]=SIF_LEFM_circle2_debug( ...
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
            oldR(ir),KIold,KIIold, ...
            KIold-KItrue,KIIold-KIItrue, ...
            hypot(KIold-KItrue,KIIold-KIItrue), ...
            DO.JI_over_absint,DO.JII_over_absint, ...
            detMismatch,baryMismatch];

        DebugOld{ir}=DO;
    end

    Told=array2table(A,'VariableNames',{ ...
        'rI','KI_old','KII_old','KI_error','KII_error', ...
        'vector_error_abs','JI_over_absint','JII_over_absint', ...
        'median_detJ_PQ_mismatch','median_baryMin_PQ_mismatch'});
    Old{ic}=struct('table',Told,'debug',{DebugOld});

    % ------------------------------------------------------------
    % Canonical FE-nodal EDI annulus sweep
    % ------------------------------------------------------------
    B=nan(size(ediDomains,1),9);
    DebugEDI=cell(size(ediDomains,1),1);

    for id=1:size(ediDomains,1)
        domain=struct('r_inner',ediDomains(id,1), ...
                      'r_outer',ediDomains(id,2));

        [KIe,KIIe,AE]=SIF_LEFM_interaction_EDI( ...
            mesh,U,V,mat,domain, ...
            'UsePlaneStrain',true,'Verbose',false, ...
            'WeightFunction','fe_nodal');

        B(id,:)=[ ...
            domain.r_inner,domain.r_outer, ...
            KIe,KIIe,KIe-KItrue,KIIe-KIItrue, ...
            hypot(KIe-KItrue,KIIe-KIItrue), ...
            AE.nElem_used,AE.nGP_used];

        DebugEDI{id}=AE;
    end

    Tedi=array2table(B,'VariableNames',{ ...
        'r_inner','r_outer','KI_EDI','KII_EDI', ...
        'KI_error','KII_error','vector_error_abs', ...
        'nElem_used','nGP_used'});
    EDI{ic}=struct('table',Tedi,'debug',{DebugEDI});

    [~,io]=min(abs(oldR-rOldRef));
    KIo=Told.KI_old(io); KIIo=Told.KII_old(io);
    KIe=Tedi.KI_EDI(ediRef); KIIe=Tedi.KII_EDI(ediRef);

    trueScale=max(hypot(KItrue,KIItrue),1);

    oldKIRange=max(Told.KI_old)-min(Told.KI_old);
    oldKIIRange=max(Told.KII_old)-min(Told.KII_old);
    ediKIRange=max(Tedi.KI_EDI)-min(Tedi.KI_EDI);
    ediKIIRange=max(Tedi.KII_EDI)-min(Tedi.KII_EDI);

    rows(ic,:)=[ ...
        ic,KItrue,KIItrue, ...
        KIo,KIIo,KIe,KIIe, ...
        KIo-KItrue,KIIo-KIItrue, ...
        KIe-KItrue,KIIe-KIItrue, ...
        hypot(KIo-KItrue,KIIo-KIItrue)/trueScale, ...
        hypot(KIe-KItrue,KIIe-KIItrue)/trueScale, ...
        oldKIRange,oldKIIRange,ediKIRange,ediKIIRange, ...
        Told.median_detJ_PQ_mismatch(io), ...
        Told.median_baryMin_PQ_mismatch(io), ...
        ginfo.auditT6.qualityCutMin];
end

T=array2table(rows,'VariableNames',{ ...
    'caseID','KI_true','KII_true', ...
    'KI_old_ref','KII_old_ref','KI_EDI_ref','KII_EDI_ref', ...
    'KI_old_error','KII_old_error','KI_EDI_error','KII_EDI_error', ...
    'old_vector_error_rel','EDI_vector_error_rel', ...
    'old_KI_radius_range','old_KII_radius_range', ...
    'EDI_KI_domain_range','EDI_KII_domain_range', ...
    'old_detJ_PQ_mismatch_ref','old_bary_PQ_mismatch_ref', ...
    'cut_quality_min'});

T.caseName=caseNames(T.caseID);
T=movevars(T,'caseName','After','caseID');

fprintf('\nGeometry gate: PASS\n');
fprintf('Parent: Ntheta=%d, Nr=%d, T3=%d\n', ...
    ginfo.parent.Ntheta,ginfo.parent.Nr,size(parent.connect3,1));
fprintf('Cut mesh: T3 nodes/elements=%d/%d; T6 nodes/elements=%d/%d\n', ...
    size(mesh.coord3,1),size(mesh.connect3,1), ...
    size(mesh.coord,1),size(mesh.connect,1));
fprintf('Crack-face pairs: T3=%d, T6=%d\n', ...
    numel(ginfo.cut.crackUpperIDs),numel(mesh.crackUpperT6IDs));
fprintf('Near-cut quality min=%.6f, angle min=%.6f deg\n', ...
    ginfo.auditT6.qualityCutMin,ginfo.auditT6.minAngleCutDeg);

fprintf('\n============================================================\n');
fprintf('REFERENCE EXTRACTION SUMMARY\n');
fprintf('old/J reference contour: rI=%.4f\n',rOldRef);
fprintf('EDI reference annulus: [%.4f, %.4f], FE-nodal q\n', ...
    ediDomains(ediRef,1),ediDomains(ediRef,2));
fprintf('============================================================\n');
disp(T);

for ic=1:nK
    fprintf('\n------------------------------------------------------------\n');
    fprintf('%s: old mirror/J contour sweep\n',caseNames(ic));
    fprintf('------------------------------------------------------------\n');
    disp(Old{ic}.table(:,{ ...
        'rI','KI_old','KII_old','KI_error','KII_error', ...
        'vector_error_abs','median_detJ_PQ_mismatch', ...
        'median_baryMin_PQ_mismatch'}));

    fprintf('\n%s: FE-nodal EDI annulus sweep\n',caseNames(ic));
    disp(EDI{ic}.table);
end

if logical(opt.Plot)
    plot_results(caseNames,Kcases,Old,EDI);
end

Out=struct();
Out.summary=T;
Out.old=Old;
Out.edi=EDI;
Out.mesh=mesh;
Out.geometryInfo=ginfo;
Out.parent=parent;
Out.settings=struct( ...
    'oldR',oldR,'rOldRef',rOldRef, ...
    'ediDomains',ediDomains,'ediRef',ediRef, ...
    'Kcases',Kcases,'caseNames',caseNames, ...
    'OldNtheta',opt.OldNtheta);

fprintf('\nSTEP 4D completed.\n');
fprintf(['Gate interpretation: S0 should show negligible pure-mode cross ', ...
    'leakage, small error against exact KI/KII, modest old-contour ', ...
    'variation, and very small FE-nodal EDI annulus variation.\n']);
end


function plot_results(caseNames,Kcases,Old,EDI)

nK=numel(caseNames);

figure('Name','Step 4D: old mirror-J contour stability','NumberTitle','off');
tiledlayout(2,2,'Padding','compact','TileSpacing','compact');
for ic=1:nK
    nexttile;
    T=Old{ic}.table;
    plot(T.rI,T.KI_old,'o-'); hold on;
    plot(T.rI,T.KII_old,'s-');
    yline(Kcases(ic,1),'--');
    yline(Kcases(ic,2),'--');
    grid on;
    xlabel('r_I');
    ylabel('Recovered SIF');
    title(caseNames(ic),'Interpreter','none');
    legend('K_I old','K_{II} old','K_I exact','K_{II} exact', ...
        'Location','best');
end

figure('Name','Step 4D: FE-nodal EDI domain stability','NumberTitle','off');
tiledlayout(2,2,'Padding','compact','TileSpacing','compact');
for ic=1:nK
    nexttile;
    T=EDI{ic}.table;
    x=1:height(T);
    plot(x,T.KI_EDI,'o-'); hold on;
    plot(x,T.KII_EDI,'s-');
    yline(Kcases(ic,1),'--');
    yline(Kcases(ic,2),'--');
    grid on;
    xticks(x);
    labels=strings(height(T),1);
    for k=1:height(T)
        labels(k)=sprintf('[%.3f,%.3f]',T.r_inner(k),T.r_outer(k));
    end
    xticklabels(labels);
    xtickangle(25);
    ylabel('Recovered SIF');
    title(caseNames(ic),'Interpreter','none');
    legend('K_I EDI','K_{II} EDI','K_I exact','K_{II} exact', ...
        'Location','best');
end
end
