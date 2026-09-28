function Out=main_step6_exact_williams_on_physical_mesh(varargin)
%MAIN_STEP6_EXACT_WILLIAMS_ON_PHYSICAL_MESH
% Exact local Williams fields sampled on the actual two-leg Crack-Path mesh.
%
% This is the strongest bridge from the synthetic mechanism to the physical
% mesh: exact KI/KII are known, but the node layout/connectivity are exactly
% those of the actual two-leg FEM control.
%
% Only contours/domains with r <= 0.6*Llast are used, so the extraction region
% lies entirely inside the straight final crack leg. The earlier kink is
% outside the local Williams audit domain.
%
% Exact cases:
%   pure I      (1,0)
%   pure II     (0,1)
%   mixed 1%    (1,0.01)
%
% The historical mirror/J and FE-nodal EDI methods are compared directly
% against the prescribed SIFs.

ip=inputParser;
addParameter(ip,'radiusFractions',[0.2 0.3 0.4 0.5 0.6], ...
    @(x)isnumeric(x)&&isvector(x)&&all(isfinite(x))&&all(x>0)&&all(x<=0.6));
addParameter(ip,'innerFactor',0.1, ...
    @(x)isnumeric(x)&&isscalar(x)&&x>0&&x<1);
addParameter(ip,'nthet',240, ...
    @(x)isnumeric(x)&&isscalar(x)&&x>=40);
parse(ip,varargin{:}); opt=ip.Results;

here=fileparts(mfilename('fullpath'));
addpath(genpath(fileparts(fileparts(here))));

fprintf('\n============================================================\n');
fprintf('SIF AUDIT STEP 6: EXACT WILLIAMS FIELD ON PHYSICAL MESH\n');
fprintf('============================================================\n');

C=cfg_crack_path_two_leg_control(2.0,'plotGeom',false,'plotMesh',false);
M=build_crack_path_polyline_LEFM_mesh(C);

mesh=struct('coord3',M.coord3,'connect3',M.connect3, ...
            'coord',M.coord,'connect',M.connect);

% Material used by the audit extractors.
E=C.E2; nu=C.nu; ps=C.ps;
coef=E/((1+nu)*(1-2*nu));
D=coef*[1-nu,nu,0;nu,1-nu,0;0,0,(1-2*nu)/2];
mat=struct('E',E,'nu',nu,'ps',ps,'D',D,'Dmat',D);

V=C.Pmid;
P1=V(end-1,:); tip=V(end,:);
d=tip-P1; Llast=norm(d);
e1=d/Llast;
e2=[-e1(2),e1(1)];

% Complete T6 node lists on upper/lower faces of the FINAL straight leg.
up3=last_leg_face_corners(M.crack.upperNodesT3,mesh.coord3,P1,tip);
lo3=last_leg_face_corners(M.crack.lowerNodesT3,mesh.coord3,P1,tip);

up6=quadratic_face(mesh,up3);
lo6=quadratic_face(mesh,lo3);

% Determine whether the semantic upper/lower chains agree with +e2/-e2.
% The pre-collapse pencil chains retain their offset and give an independent
% orientation check near the final leg.
[up6,lo6,faceOrientation]=orient_face_labels_from_pencil( ...
    up6,lo6,M.Ggeom,P1,tip,e2);

mesh.crackUpperT6IDs=up6;
mesh.crackLowerT6IDs=lo6;

fprintf('Final leg length = %.8g\n',Llast);
fprintf('Final-leg T3 face corners upper/lower = %d / %d\n',numel(up3),numel(lo3));
fprintf('Final-leg T6 face nodes upper/lower = %d / %d\n',numel(up6),numel(lo6));
fprintf('Face-label orientation check = %s\n',faceOrientation);

% Local coordinates relative to the current crack tip and final-leg frame.
X=mesh.coord;
Xrel=X-tip;
Xloc=[Xrel*e1.', Xrel*e2.'];

Kcases=[1,0;0,1;1,0.01];
caseNames=["pure_I";"pure_II";"mixed_1pct"];

rf=opt.radiusFractions(:);
nR=numel(rf); nK=size(Kcases,1);

rows=nan(nR*nK,22);
DebugOld=cell(nR,nK);
DebugEDI=cell(nR,nK);
rr=0;

for ic=1:nK
    KItrue=Kcases(ic,1); KIItrue=Kcases(ic,2);

    Uloc=exact_williams_displacement_audit( ...
        Xloc,KItrue,KIItrue,E,nu,ps, ...
        'UpperFaceIDs',up6,'LowerFaceIDs',lo6);

    ul=[Uloc(1:2:end),Uloc(2:2:end)];
    ug=ul(:,1)*e1 + ul(:,2)*e2;
    U=zeros(2*size(X,1),1);
    U(1:2:end)=ug(:,1);
    U(2:2:end)=ug(:,2);

    for ir=1:nR
        rI=rf(ir)*Llast;

        [KIo,KIIo,DO]=SIF_LEFM_circle2_debug( ...
            mesh,U,V,mat,rI, ...
            'nthet',opt.nthet,'plot',false,'verbose',false);

        dom=struct('r_inner',opt.innerFactor*rI,'r_outer',rI);
        [KIe,KIIe,AE]=SIF_LEFM_interaction_EDI( ...
            mesh,U,V,mat,dom, ...
            'UsePlaneStrain',true,'Verbose',false, ...
            'WeightFunction','fe_nodal');

        detDen=0.5*(abs(DO.detJP)+abs(DO.detJQ));
        good=detDen>0;
        detMis=NaN;
        if any(good)
            detMis=median(abs(DO.detJP(good)-DO.detJQ(good))./detDen(good));
        end
        baryMis=median(abs(DO.baryMinP-DO.baryMinQ));

        S=DO.diagnostics;
        trueScale=max(hypot(KItrue,KIItrue),1);

        if abs(KIItrue)>0
            relKIIold=(KIIo-KIItrue)/abs(KIItrue);
            relKIIedi=(KIIe-KIItrue)/abs(KIItrue);
        else
            relKIIold=NaN; relKIIedi=NaN;
        end

        rr=rr+1;
        rows(rr,:)=[ ...
            ic,rf(ir),rI,KItrue,KIItrue, ...
            KIo,KIIo,KIe,KIIe, ...
            KIo-KItrue,KIIo-KIItrue,KIe-KItrue,KIIe-KIItrue, ...
            hypot(KIo-KItrue,KIIo-KIItrue)/trueScale, ...
            hypot(KIe-KItrue,KIIe-KIItrue)/trueScale, ...
            relKIIold,relKIIedi, ...
            detMis,baryMis, ...
            S.mirrorT3Mismatch_median,S.mirrorT3Mismatch_p95, ...
            S.frac_exact_mirror_T3];

        DebugOld{ir,ic}=DO;
        DebugEDI{ir,ic}=AE;
    end
end

T=array2table(rows,'VariableNames',{ ...
    'caseID','r_over_lastLeg','rI','KI_true','KII_true', ...
    'KI_old','KII_old','KI_EDI','KII_EDI', ...
    'KI_old_error','KII_old_error','KI_EDI_error','KII_EDI_error', ...
    'old_vector_error_rel','EDI_vector_error_rel', ...
    'relative_KII_error_old','relative_KII_error_EDI', ...
    'median_detJ_PQ_mismatch','median_bary_PQ_mismatch', ...
    'mirror_T3_mismatch_median','mirror_T3_mismatch_p95', ...
    'fraction_exact_mirror_T3'});

T.caseName=caseNames(T.caseID);
T=movevars(T,'caseName','After','caseID');

fprintf('\nEXACT-FIELD RESULTS ON THE PHYSICAL MESH\n');
disp(T);

% Compact detective view for the most relevant 1%-mixed case.
Tm=T(T.caseID==3,:);
fprintf('\n1%%-MIXED DETECTIVE VIEW\n');
disp(Tm(:,{ ...
    'r_over_lastLeg','KI_old','KII_old','KI_EDI','KII_EDI', ...
    'relative_KII_error_old','relative_KII_error_EDI', ...
    'mirror_T3_mismatch_median','mirror_T3_mismatch_p95', ...
    'fraction_exact_mirror_T3'}));

Out=struct();
Out.results=T;
Out.mixed1pct=Tm;
Out.mesh=mesh;
Out.geometry=M;
Out.caseConfig=C;
Out.debugOld=DebugOld;
Out.debugEDI=DebugEDI;
Out.settings=opt;
Out.localFrame=struct('tip',tip,'e1',e1,'e2',e2, ...
    'Llast',Llast,'faceOrientation',faceOrientation);

fprintf('\nSTEP 6 completed.\n');
fprintf(['Because exact KI/KII are prescribed here, any extraction error is ', ...
    'attributable to the extractor plus this physical mesh/interpolation ', ...
    'layout, not to the FEM solution field.\n']);
end


function ids=last_leg_face_corners(faceIDs,X,P1,P2)
faceIDs=faceIDs(:);
d=P2-P1; L2=dot(d,d);
Q=X(faceIDs,:);
t=((Q-P1)*d.')/L2;
proj=P1+t.*d;
dist=sqrt(sum((Q-proj).^2,2));
tol=max(1e-10,1e-7*sqrt(L2));
keep=dist<=tol & t>=-1e-8 & t<=1+1e-8;
ids=faceIDs(keep);
[~,ord]=sort(t(keep));
ids=ids(ord);
if numel(ids)<2
    error('step6:FinalLegFace','Failed to identify final-leg face corners.');
end
end


function ids=quadratic_face(mesh,corners)
edges=[mesh.connect(:,[1 2]);mesh.connect(:,[2 3]);mesh.connect(:,[3 1])];
mids=[mesh.connect(:,4);mesh.connect(:,5);mesh.connect(:,6)];
onFace=all(ismember(edges,corners),2);
ids=unique([corners(:);mids(onFace)]);
end


function [upper,lower,status]=orient_face_labels_from_pencil( ...
    upper,lower,Ggeom,P1,P2,e2)

% The final portions of the pre-collapse offset chains indicate which stored
% chain lies on +e2. Use the last noncoincident offset sample available.

u=Ggeom.up_chain;
d=Ggeom.dn_chain;

nu=size(u,1); nd=size(d,1);
iu=max(1,nu-2):nu;
id=max(1,nd-2):nd;
su=mean((u(iu,:)-P2)*e2.');
sd=mean((d(id,:)-P2)*e2.');

if su>sd
    status='stored upper chain is +e2';
else
    tmp=upper; upper=lower; lower=tmp;
    status='stored face labels swapped to make upper=+e2';
end
end
