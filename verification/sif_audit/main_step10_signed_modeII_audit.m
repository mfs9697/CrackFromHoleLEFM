function Out=main_step10_signed_modeII_audit(varargin)
%MAIN_STEP10_SIGNED_MODEII_AUDIT
% Explicitly audit the physical sign of KII.
%
% Part A: exact symmetric S0 crack-cut mesh.
% Part B: actual production Stage-II hole-crack mesh at theta=0 deg.
%
% Exact cases include positive and negative Mode II:
%   (0,+1), (0,-1), (1,+0.01), (1,-0.01)
%
% The purpose is not to re-test magnitude accuracy. It is to demonstrate that
% the modal JII contribution is quadratic in KII and therefore cannot, by
% itself, determine the physical sign of KII. The interaction EDI should
% recover both signs.

ip=inputParser;
addParameter(ip,'OldNtheta',240,@(x)isnumeric(x)&&isscalar(x)&&x>=40);
parse(ip,varargin{:});
O=ip.Results;

addpath(genpath(pwd));

fprintf('\n============================================================\n');
fprintf('SIF AUDIT STEP 10: SIGNED MODE-II AUDIT\n');
fprintf('============================================================\n');

Kcases=[0,1;0,-1;1,0.01;1,-0.01];
caseNames=["pure_II_plus";"pure_II_minus";"mixed_plus_1pct";"mixed_minus_1pct"];

%% ------------------------------------------------------------
% Part A: exact symmetric S0 mesh
%% ------------------------------------------------------------
[meshS0,ginfo]=build_literal_ring_crack_cut_mesh(); %#ok<ASGLU>
E=4e3; nu=0.30; ps=1;
coef=E/((1+nu)*(1-2*nu));
D=coef*[1-nu,nu,0;nu,1-nu,0;0,0,(1-2*nu)/2];
mat=struct('E',E,'nu',nu,'ps',ps,'D',D,'Dmat',D);
V=[-1,0;0,0];
oldR=0.08;
dom=struct('r_inner',0.02,'r_outer',0.12);

A=nan(size(Kcases,1),11);
for ic=1:size(Kcases,1)
    KItrue=Kcases(ic,1); KIItrue=Kcases(ic,2);
    U=exact_williams_displacement_audit( ...
        meshS0.coord,KItrue,KIItrue,E,nu,ps, ...
        'UpperFaceIDs',meshS0.crackUpperT6IDs, ...
        'LowerFaceIDs',meshS0.crackLowerT6IDs);

    [KIo,KIIo,DO]=SIF_LEFM_circle2_debug( ...
        meshS0,U,V,mat,oldR, ...
        'nthet',O.OldNtheta,'plot',false,'verbose',false);

    [KIe,KIIe]=SIF_LEFM_interaction_EDI( ...
        meshS0,U,V,mat,dom, ...
        'UsePlaneStrain',true,'Verbose',false, ...
        'WeightFunction','fe_nodal');

    A(ic,:)=[ic,KItrue,KIItrue,KIo,KIIo,abs(KIIo),DO.JII, ...
             KIe,KIIe,KIIo-KIItrue,KIIe-KIItrue];
end

TS0=array2table(A,'VariableNames',{ ...
    'caseID','KI_true','KII_true','KI_old','KII_old_legacy', ...
    'abs_KII_old','JII_old','KI_EDI','KII_EDI', ...
    'KII_old_error','KII_EDI_error'});
TS0.caseName=caseNames(TS0.caseID);
TS0=movevars(TS0,'caseName','After','caseID');

fprintf('\nPART A: EXACT SYMMETRIC S0 MESH\n');
disp(TS0);

%% ------------------------------------------------------------
% Part B: actual production mesh at theta=0
%% ------------------------------------------------------------
C=cfg_hole_initiation();
G=geom_hole_only(C);
S1=solve_hole_only(C,G,'lambda',1.0);
B=sample_hole_boundary_stress(C,G,S1);
I=find_hole_initiation_point(C,B);

[G2,Ddesc,M,Mc]=build_stage2_cracked_mesh_for_theta( ...
    C,I,0,'PlotGeom',false,'PlotMesh',false,'PlotCollapsed',false); %#ok<ASGLU>
S2=solve_cracked_LEFM(C,Mc);

V2=Mc.crack.Pmid;
tip=V2(end,:); prev=V2(end-1,:);
e1=(tip-prev)/norm(tip-prev);
e2=[-e1(2),e1(1)];

up3=Mc.crack.upperNodes(:);
lo3=Mc.crack.lowerNodes(:);
up6=quadratic_face(S2.mesh,up3);
lo6=quadratic_face(S2.mesh,lo3);

du=mean((Mc.p0(up3,:)-Mc.p(up3,:))*e2.');
dl=mean((Mc.p0(lo3,:)-Mc.p(lo3,:))*e2.');
if du<dl
    tmp=up6; up6=lo6; lo6=tmp;
end

X=S2.mesh.coord;
Xrel=X-tip;
Xloc=[Xrel*e1.',Xrel*e2.'];

Btab=nan(size(Kcases,1),14);
for ic=1:size(Kcases,1)
    KItrue=Kcases(ic,1); KIItrue=Kcases(ic,2);

    Uloc=exact_williams_displacement_audit( ...
        Xloc,KItrue,KIItrue,C.E,C.nu,C.ps, ...
        'UpperFaceIDs',up6,'LowerFaceIDs',lo6);

    ul=[Uloc(1:2:end),Uloc(2:2:end)];
    ug=ul(:,1)*e1+ul(:,2)*e2;
    U=zeros(2*size(X,1),1);
    U(1:2:end)=ug(:,1);
    U(2:2:end)=ug(:,2);

    Sx=S2; Sx.U=U;
    R=compute_SIF_for_stage2_compare(C,G2,Mc,Sx,'OldNtheta',O.OldNtheta);

    Btab(ic,:)=[ic,KItrue,KIItrue,R.KI_old,R.KII_old,abs(R.KII_old), ...
        R.old.JII,R.KI_EDI,R.KII_EDI, ...
        R.KII_old-KIItrue,R.KII_EDI-KIItrue, ...
        R.domain_EDI.r_inner,R.domain_EDI.r_outer, ...
        R.stencil.fraction_exact_mirror_T3];
end

Tprod=array2table(Btab,'VariableNames',{ ...
    'caseID','KI_true','KII_true','KI_old','KII_old_legacy', ...
    'abs_KII_old','JII_old','KI_EDI','KII_EDI', ...
    'KII_old_error','KII_EDI_error', ...
    'EDI_r_inner','EDI_r_outer','fraction_exact_mirror_T3'});
Tprod.caseName=caseNames(Tprod.caseID);
Tprod=movevars(Tprod,'caseName','After','caseID');

fprintf('\nPART B: ACTUAL PRODUCTION STAGE-II MESH, THETA=0\n');
disp(Tprod);

fprintf('\nSIGN CHECK SUMMARY\n');
fprintf(['For exact +/-KII pairs, a modal-energy JII should remain nearly ', ...
    'unchanged because JII is quadratic in KII. Therefore the sign of ', ...
    'KII cannot be reconstructed from sign(JII).\n']);

Out=struct();
Out.S0=TS0;
Out.production=Tprod;
Out.Kcases=Kcases;
Out.caseNames=caseNames;
Out.settings=O;

fprintf('\nSTEP 10 completed.\n');
end


function ids=quadratic_face(mesh,corners)
edges=[mesh.connect(:,[1 2]);mesh.connect(:,[2 3]);mesh.connect(:,[3 1])];
mids=[mesh.connect(:,4);mesh.connect(:,5);mesh.connect(:,6)];
onFace=all(ismember(edges,corners),2);
ids=unique([corners(:);mids(onFace)]);
end
