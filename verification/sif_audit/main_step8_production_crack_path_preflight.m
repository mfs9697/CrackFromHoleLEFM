function Out=main_step8_production_crack_path_preflight(varargin)
%MAIN_STEP8_PRODUCTION_CRACK_PATH_PREFLIGHT
% Exact-field gate on the ACTUAL hole+short-crack production mesh.
%
% This is the first production-workflow comparison after the synthetic audit.
% It intentionally runs only one representative crack direction first
% (theta=0 by default), generates the normal production Stage-II mesh, and
% replaces the numerical displacement field by exact local Williams fields.
%
% The old mirror/J extractor and FE-nodal interaction EDI are then compared
% against known KI/KII on exactly the mesh that the hole-crack workflow uses.
%
% After this gate is accepted, run the full physical theta sweep with both
% extractors on the same solved FEM fields.

ip=inputParser;
addParameter(ip,'thetaDeg',0,@(x)isnumeric(x)&&isscalar(x)&&isfinite(x));
addParameter(ip,'OldNtheta',240,@(x)isnumeric(x)&&isscalar(x)&&x>=40);
addParameter(ip,'PlotMesh',true,@(x)islogical(x)||isnumeric(x));
parse(ip,varargin{:});
O=ip.Results;

addpath(genpath(pwd));

fprintf('\n============================================================\n');
fprintf('SIF AUDIT STEP 8: PRODUCTION HOLE-CRACK PREFLIGHT\n');
fprintf('============================================================\n');

C=cfg_hole_initiation();

% Stage I is the unchanged production initiation calculation.
G=geom_hole_only(C);
S1=solve_hole_only(C,G,'lambda',1.0);
B=sample_hole_boundary_stress(C,G,S1);
I=find_hole_initiation_point(C,B);

theta=deg2rad(O.thetaDeg);
[G2,D,M,Mc]=build_stage2_cracked_mesh_for_theta( ...
    C,I,theta, ...
    'PlotGeom',false,'PlotMesh',false,'PlotCollapsed',false);

S2=solve_cracked_LEFM(C,Mc);

if logical(O.PlotMesh)
    local_plot_production_mesh(S2.mesh,Mc,G2,O.thetaDeg);
end

% Local final-leg frame.
V=Mc.crack.Pmid;
tip=V(end,:);
prev=V(end-1,:);
e1=(tip-prev)/norm(tip-prev);
e2=[-e1(2),e1(1)];
Llast=norm(tip-prev);

% Recover complete T6 crack-face node lists from the collapsed T3 face sets.
up3=Mc.crack.upperNodes(:);
lo3=Mc.crack.lowerNodes(:);
up6=quadratic_face(S2.mesh,up3);
lo6=quadratic_face(S2.mesh,lo3);

% Use the PRE-COLLAPSE offsets to assign the +e2 and -e2 branches.
du=mean((Mc.p0(up3,:)-Mc.p(up3,:))*e2.');
dl=mean((Mc.p0(lo3,:)-Mc.p(lo3,:))*e2.');
if du<dl
    tmp=up6; up6=lo6; lo6=tmp;
    faceStatus='stored face labels swapped so upper=+e2';
else
    faceStatus='stored upper face corresponds to +e2';
end

shared=intersect(up6,lo6);
if ~isempty(shared)
    dtip=sqrt(sum((S2.mesh.coord(shared,:)-tip).^2,2));
    tolTip=max(1e-12,1e-9*Llast);
    if any(dtip>tolTip)
        error('step8:SharedNonTipNode', ...
            'Upper/lower T6 face sets share a non-tip node.');
    end
end

fprintf('theta = %.6f deg\n',O.thetaDeg);
fprintf('short-crack length Llast = %.8g\n',Llast);
fprintf('T3 face nodes upper/lower = %d / %d\n',numel(up3),numel(lo3));
fprintf('T6 face nodes upper/lower = %d / %d\n',numel(up6),numel(lo6));
fprintf('shared T6 face nodes = %d (tip-only expected)\n',numel(shared));
fprintf('face orientation = %s\n',faceStatus);

% Exact Williams displacement field in local coordinates.
X=S2.mesh.coord;
Xrel=X-tip;
Xloc=[Xrel*e1.', Xrel*e2.'];

Kcases=[1,0;0,1;1,0.01];
caseNames=["pure_I";"pure_II";"mixed_1pct"];

rows=nan(size(Kcases,1),25);
Cmp=cell(size(Kcases,1),1);

for ic=1:size(Kcases,1)
    KItrue=Kcases(ic,1);
    KIItrue=Kcases(ic,2);

    Uloc=exact_williams_displacement_audit( ...
        Xloc,KItrue,KIItrue,C.E,C.nu,C.ps, ...
        'UpperFaceIDs',up6,'LowerFaceIDs',lo6);

    ul=[Uloc(1:2:end),Uloc(2:2:end)];
    ug=ul(:,1)*e1 + ul(:,2)*e2;

    U=zeros(2*size(X,1),1);
    U(1:2:end)=ug(:,1);
    U(2:2:end)=ug(:,2);

    Sx=S2;
    Sx.U=U;

    R=compute_SIF_for_stage2_compare( ...
        C,G2,Mc,Sx, ...
        'OldNtheta',O.OldNtheta);

    Cmp{ic}=R;

    if abs(KIItrue)>0
        relKIIold=(R.KII_old-KIItrue)/abs(KIItrue);
        relKIIedi=(R.KII_EDI-KIItrue)/abs(KIItrue);
    else
        relKIIold=NaN;
        relKIIedi=NaN;
    end

    trueScale=max(hypot(KItrue,KIItrue),1);

    rows(ic,:)=[ ...
        ic,KItrue,KIItrue, ...
        R.KI_old,R.KII_old,R.KI_EDI,R.KII_EDI, ...
        R.KI_old-KItrue,R.KII_old-KIItrue, ...
        R.KI_EDI-KItrue,R.KII_EDI-KIItrue, ...
        hypot(R.KI_old-KItrue,R.KII_old-KIItrue)/trueScale, ...
        hypot(R.KI_EDI-KItrue,R.KII_EDI-KIItrue)/trueScale, ...
        relKIIold,relKIIedi, ...
        R.r_old,R.domain_EDI.r_inner,R.domain_EDI.r_outer, ...
        R.tipMeshScale.min,R.tipMeshScale.median,R.tipMeshScale.max, ...
        R.stencil.mirror_T3_mismatch_median, ...
        R.stencil.mirror_T3_mismatch_p95, ...
        R.stencil.mirror_T3_mismatch_max, ...
        R.stencil.fraction_exact_mirror_T3];
end

T=array2table(rows,'VariableNames',{ ...
    'caseID','KI_true','KII_true', ...
    'KI_old','KII_old','KI_EDI','KII_EDI', ...
    'KI_old_error','KII_old_error','KI_EDI_error','KII_EDI_error', ...
    'old_vector_error_rel','EDI_vector_error_rel', ...
    'relative_KII_error_old','relative_KII_error_EDI', ...
    'old_radius','EDI_r_inner','EDI_r_outer', ...
    'h_tip_min','h_tip_median','h_tip_max', ...
    'mirror_T3_mismatch_median','mirror_T3_mismatch_p95', ...
    'mirror_T3_mismatch_max','fraction_exact_mirror_T3'});

T.caseName=caseNames(T.caseID);
T=movevars(T,'caseName','After','caseID');

fprintf('\nPRODUCTION-MESH EXACT-FIELD PREFLIGHT\n');
disp(T);

fprintf('\n1%%-MIXED PREFLIGHT VIEW\n');
Tm=T(T.caseID==3,:);
disp(Tm(:,{ ...
    'KI_old','KII_old','KI_EDI','KII_EDI', ...
    'relative_KII_error_old','relative_KII_error_EDI', ...
    'old_radius','EDI_r_inner','EDI_r_outer','h_tip_median', ...
    'mirror_T3_mismatch_p95','fraction_exact_mirror_T3'}));

Out=struct();
Out.config=C;
Out.initiation=I;
Out.G2=G2;
Out.domainDescription=D;
Out.meshOriginal=M;
Out.meshCollapsed=Mc;
Out.solve=S2;
Out.results=T;
Out.mixed1pct=Tm;
Out.compare=Cmp;
Out.thetaDeg=O.thetaDeg;
Out.faceStatus=faceStatus;
Out.upperFaceT6=up6;
Out.lowerFaceT6=lo6;

fprintf('\nSTEP 8 completed.\n');
fprintf(['This is the gate before the full production theta sweep. ', ...
    'Both extractors have now been tested against exact truth on the ', ...
    'actual hole-crack Stage-II mesh.\n']);
end


function ids=quadratic_face(mesh,corners)
edges=[mesh.connect(:,[1 2]);mesh.connect(:,[2 3]);mesh.connect(:,[3 1])];
mids=[mesh.connect(:,4);mesh.connect(:,5);mesh.connect(:,6)];
onFace=all(ismember(edges,corners),2);
ids=unique([corners(:);mids(onFace)]);
end


function local_plot_production_mesh(mesh,Mc,G2,thetaDeg)
figure('Name','Step 8: production hole-crack mesh','Color','w');
clf;
triplot(mesh.connect3,mesh.coord3(:,1),mesh.coord3(:,2), ...
    'Color',[0.72 0.72 0.72]);
hold on; axis equal; box on;
plot(Mc.crack.Pmid(:,1),Mc.crack.Pmid(:,2),'k-','LineWidth',2);
plot(Mc.crack.xtip(1),Mc.crack.xtip(2),'ko','MarkerSize',6,'LineWidth',1.4);
xlabel('x'); ylabel('y');
title(sprintf('Production Stage-II mesh, \\theta=%.2f^\\circ',thetaDeg));

% Crack-tip zoom in a separate figure for thesis/audit export.
figure('Name','Step 8: production crack-tip mesh zoom','Color','w');
clf;
triplot(mesh.connect3,mesh.coord3(:,1),mesh.coord3(:,2), ...
    'Color',[0.65 0.65 0.65]);
hold on; axis equal; box on;
plot(Mc.crack.Pmid(:,1),Mc.crack.Pmid(:,2),'k-','LineWidth',2);
plot(Mc.crack.xtip(1),Mc.crack.xtip(2),'ko','MarkerSize',6,'LineWidth',1.4);
R=max(G2.tip.radiusJ,0.6*G2.crack.a0);
xlim(Mc.crack.xtip(1)+[-R,R]);
ylim(Mc.crack.xtip(2)+[-R,R]);
xlabel('x'); ylabel('y');
title('Production Stage-II crack-tip mesh zoom');
end
