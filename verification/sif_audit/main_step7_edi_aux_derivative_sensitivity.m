function Out=main_step7_edi_aux_derivative_sensitivity(varargin)
%MAIN_STEP7_EDI_AUX_DERIVATIVE_SENSITIVITY
% Final implementation-specific gate for the EDI prototype.
%
% Uses the ACTUAL two-leg Crack-Path mesh and exact Williams nodal fields,
% but fixes the extraction domain at the well-resolved Step-6 choice
% r_outer=0.5*Llast, r_inner=0.05*Llast (=2*htip for ncoh=40).
%
% Sweeps the multiplier applied to the centered finite-difference step used
% for auxiliary displacement gradients. The exact SIFs are known.
%
% Usage:
%   O7 = main_step7_edi_aux_derivative_sensitivity();

ip=inputParser;
addParameter(ip,'Scales',[0.01 0.03 0.1 0.3 1 3 10 30], ...
    @(x)isnumeric(x)&&isvector(x)&&all(isfinite(x))&&all(x>0));
parse(ip,varargin{:}); opt=ip.Results;

here=fileparts(mfilename('fullpath'));
addpath(genpath(fileparts(fileparts(here))));

fprintf('\n============================================================\n');
fprintf('SIF AUDIT STEP 7: EDI AUXILIARY-DERIVATIVE SENSITIVITY\n');
fprintf('============================================================\n');

C=cfg_crack_path_two_leg_control(2.0,'plotGeom',false,'plotMesh',false);
M=build_crack_path_polyline_LEFM_mesh(C);

mesh=struct('coord3',M.coord3,'connect3',M.connect3, ...
            'coord',M.coord,'connect',M.connect);

E=C.E2; nu=C.nu; ps=C.ps;
coef=E/((1+nu)*(1-2*nu));
D=coef*[1-nu,nu,0;nu,1-nu,0;0,0,(1-2*nu)/2];
mat=struct('E',E,'nu',nu,'ps',ps,'D',D,'Dmat',D);

V=C.Pmid;
P1=V(end-1,:); tip=V(end,:);
d=tip-P1; Llast=norm(d);
e1=d/Llast;
e2=[-e1(2),e1(1)];

up3=last_leg_face_corners(M.crack.upperNodesT3,mesh.coord3,P1,tip);
lo3=last_leg_face_corners(M.crack.lowerNodesT3,mesh.coord3,P1,tip);
up6=quadratic_face(mesh,up3);
lo6=quadratic_face(mesh,lo3);
[up6,lo6]=orient_face_labels_from_pencil(up6,lo6,M.Ggeom,tip,e2);

X=mesh.coord;
Xrel=X-tip;
Xloc=[Xrel*e1.', Xrel*e2.'];

rOuter=0.5*Llast;
rInner=0.05*Llast;
domain=struct('r_inner',rInner,'r_outer',rOuter);

fprintf('Llast = %.8g, htip = %.8g\n',Llast,M.htip);
fprintf('EDI annulus = [%.8g, %.8g] = [%.3f, %.3f] Llast\n', ...
    rInner,rOuter,rInner/Llast,rOuter/Llast);
fprintf('r_inner / htip = %.3f\n',rInner/M.htip);

Kcases=[1,0;0,1;1,0.01];
caseNames=["pure_I";"pure_II";"mixed_1pct"];

sc=opt.Scales(:);
nS=numel(sc); nK=size(Kcases,1);
rows=nan(nS*nK,17);
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

    for is=1:nS
        [KIe,KIIe,AE]=SIF_LEFM_interaction_EDI( ...
            mesh,U,V,mat,domain, ...
            'UsePlaneStrain',true,'Verbose',false, ...
            'WeightFunction','fe_nodal', ...
            'AuxDerivativeScale',sc(is));

        trueScale=max(hypot(KItrue,KIItrue),1);
        if abs(KIItrue)>0
            relKII=(KIIe-KIItrue)/abs(KIItrue);
        else
            relKII=NaN;
        end

        rr=rr+1;
        rows(rr,:)=[ ...
            ic,sc(is),KItrue,KIItrue, ...
            KIe,KIIe,KIe-KItrue,KIIe-KIItrue, ...
            hypot(KIe-KItrue,KIIe-KIItrue)/trueScale, ...
            relKII, ...
            AE.auxEpsMismatchI_median,AE.auxEpsMismatchI_max, ...
            AE.auxEpsMismatchII_median,AE.auxEpsMismatchII_max, ...
            AE.nElem_used,AE.nGP_used,rInner/M.htip];
    end
end

T=array2table(rows,'VariableNames',{ ...
    'caseID','derivative_scale','KI_true','KII_true', ...
    'KI_EDI','KII_EDI','KI_error','KII_error', ...
    'vector_error_rel','relative_KII_error', ...
    'aux_eps_mismatch_I_median','aux_eps_mismatch_I_max', ...
    'aux_eps_mismatch_II_median','aux_eps_mismatch_II_max', ...
    'nElem_used','nGP_used','r_inner_over_htip'});

T.caseName=caseNames(T.caseID);
T=movevars(T,'caseName','After','caseID');

fprintf('\nAUXILIARY-DERIVATIVE SENSITIVITY TABLE\n');
disp(T);

Tm=T(T.caseID==3,:);
fprintf('\n1%%-MIXED VIEW\n');
disp(Tm(:,{ ...
    'derivative_scale','KI_EDI','KII_EDI', ...
    'KI_error','KII_error','relative_KII_error', ...
    'aux_eps_mismatch_I_median','aux_eps_mismatch_I_max', ...
    'aux_eps_mismatch_II_median','aux_eps_mismatch_II_max'}));

Out=struct('results',T,'mixed1pct',Tm,'mesh',mesh,'geometry',M, ...
    'caseConfig',C,'domain',domain,'settings',opt);

fprintf('\nSTEP 7 completed.\n');
fprintf(['A broad flat region versus derivative_scale would show that the ', ...
    'remaining EDI error is not controlled by the finite-difference step.\n']);
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
    error('step7:FinalLegFace','Failed to identify final-leg face corners.');
end
end


function ids=quadratic_face(mesh,corners)
edges=[mesh.connect(:,[1 2]);mesh.connect(:,[2 3]);mesh.connect(:,[3 1])];
mids=[mesh.connect(:,4);mesh.connect(:,5);mesh.connect(:,6)];
onFace=all(ismember(edges,corners),2);
ids=unique([corners(:);mids(onFace)]);
end


function [upper,lower]=orient_face_labels_from_pencil(upper,lower,Ggeom,P2,e2)
u=Ggeom.up_chain; d=Ggeom.dn_chain;
nu=size(u,1); nd=size(d,1);
iu=max(1,nu-2):nu; id=max(1,nd-2):nd;
su=mean((u(iu,:)-P2)*e2.');
sd=mean((d(id,:)-P2)*e2.');
if su<=sd
    tmp=upper; upper=lower; lower=tmp;
end
end
