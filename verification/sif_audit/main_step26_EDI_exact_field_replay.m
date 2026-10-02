function O26=main_step26_EDI_exact_field_replay(O25,varargin)
%MAIN_STEP26_EDI_EXACT_FIELD_REPLAY
% Replay exact leading-order Williams displacement fields on the ACTUAL
% Step-25 cracked T6 mesh, with all quadrature settings unchanged.
%
% Purpose: separate the EDI extraction error of this same mesh from
% errors in the actual solved FEM displacement field.
%
% The crack-face duplicated T3 node IDs are supplied by O25.crack.
% Original mesh T6 mid-edge nodes are assigned to a face if both edge
% endpoints belong to that face. Exact displacements on upper and lower
% crack faces are evaluated with theta=+pi and -pi, respectively.
% If an unclassified crack-face node is found inside r_outer, abort.
%
% For each EDI annulus, compute recovery matrix M for:
%   U_exactI  -> [KI,KII] for imposed [1,0]
%   U_exactII -> [KI,KII] for imposed [0,1]
% If M is near identity and path-independent, same-mesh EDI can recover
% an exact singular field: poor path independence of real U is then
% attributable to the computed displacement field or other non-Williams
% effects, not a blanket normalization failure.
%
% NOTE: M^-1*[KI_real,KII_real] is a diagnostic only. The exact field
% differs from the actual boundary-value problem, so same-mesh
% calibration CANNOT be accepted as a physical correction on its own.
%
% Usage:
%   O26=main_step26_EDI_exact_field_replay(O25);

p=inputParser;
addParameter(p,'ROuterOverA0',O25.rOuterOverA0(:).', ...
    @(x)isnumeric(x)&&isvector(x)&&all(isfinite(x))&&all(x>0)&&all(x<1));
addParameter(p,'CommonInner',O25.settings.CommonInner, ...
    @(x)isnumeric(x)&&isscalar(x)&&isfinite(x)&&x>0);
addParameter(p,'MaxMatrixError',0.05, ...
    @(x)isnumeric(x)&&isscalar(x)&&isfinite(x)&&x>0);
parse(p,varargin{:});
O=p.Results;

must(O25,'mesh'); must(O25,'mat'); must(O25,'crack'); must(O25,'U');
mesh=O25.mesh;mat=O25.mat;crack=O25.crack;
V=crack.Pmid;
tip=V(end,:);axis=V(end,:)-V(end-1,:);
axis=axis/norm(axis);
eperp=[-axis(2),axis(1)];
x=mesh.coord;
xl=(x-tip)*[axis(:),eperp(:)];
r=hypot(xl(:,1),xl(:,2));
n=size(x,1);
nt3=size(mesh.coord3,1);
if nt3>n
    error('step26:BadMesh','T6 node count smaller than T3 node count.');
end

% Assign duplicate face-side signs using the authoritative T3 topology.
face=zeros(n,1);
upper=unique(crack.upperNodes(:));
lower=unique(crack.lowerNodes(:));
tipID=crack.tipNode;
upper=upper(upper~=tipID);
lower=lower(lower~=tipID);
uOnly=setdiff(upper,lower);
lOnly=setdiff(lower,upper);
face(uOnly)=+1;
face(lOnly)=-1;
if any(ismember(uOnly,lOnly))
    error('step26:FaceOverlap','Upper and lower face classifications overlap.');
end

% Classify all midside nodes on boundary edges by their T3 endpoints.
T6=mesh.connect;
triEdges=[1 2 4;2 3 5;3 1 6];
for j=1:3
    edges=T6(:,triEdges(j,:));
    markU=face(edges(:,1))==+1 & face(edges(:,2))==+1;
    markL=face(edges(:,1))==-1 & face(edges(:,2))==-1;
    midsU=unique(edges(markU,3));
    midsL=unique(edges(markL,3));
    if any(face(midsU)==-1)||any(face(midsL)==+1)
        error('step26:MidfaceConflict', ...
            'A midside node belongs to both upper/lower crack faces.');
    end
    face(midsU)=+1;
    face(midsL)=-1;
end

a0=norm(V(end,:)-V(1,:));
rMax=max(O.ROuterOverA0)*a0;
tolLine=max(5e-12,1e-8*a0);
onNegAxis=xl(:,1)<-tolLine & abs(xl(:,2))<tolLine & r<rMax+tolLine;
unclassified=find(onNegAxis & face==0);
if ~isempty(unclassified)
    error('step26:UnclassifiedFace', ...
        ['%d T6 nodes on the crack-face negative axis inside EDI ', ...
         'support are not classified as upper/lower. ', ...
         'Exact Williams replay would be discontinuous there.'], ...
        numel(unclassified));
end

th=atan2(xl(:,2),xl(:,1));
nearFace=xl(:,1)<-tolLine & abs(xl(:,2))<tolLine;
th(nearFace & face==+1)=+pi;
th(nearFace & face==-1)=-pi;
th(r<tolLine)=0;

fprintf('\n============================================================\n');
fprintf('STEP 26: EXACT WILLIAMS FIELD REPLAY ON ACTUAL T6 MESH\n');
fprintf('============================================================\n');
fprintf('  tip-edge median=%.8e; a0=%.8e; T6 nodes=%d\n', ...
    O25.hTip,a0,n);
fprintf('  face T6 nodes assigned: upper=%d lower=%d; unclassified=%d\n', ...
    nnz(face==1),nnz(face==-1),numel(unclassified));

UI=make_exact_displacements(xl,r,th,mat,1,0);
UII=make_exact_displacements(xl,r,th,mat,0,1);

if ~isfield(mat,'Dmat'),mat.Dmat=mat.D;end
rat=O.ROuterOverA0(:).';
nR=numel(rat);
M=nan(2,2,nR);
Kactual=nan(2,nR);
calibrated=nan(2,nR);
rows=nan(nR,15);

for ir=1:nR
    ro=rat(ir)*a0;
    ri=O.CommonInner;
    if ~(ro>ri && ri>=2*O25.hTip)
        error('step26:BadAnnulus', ...
            'Need 2*h_tip<=r_inner<r_outer: %.6g, %.6g, %.6g.', ...
            O25.hTip,ri,ro);
    end
    dom=struct('r_inner',ri,'r_outer',ro);

    [a,b,A1]=SIF_LEFM_interaction_EDI( ...
        mesh,UI,V,mat,dom,'UsePlaneStrain',mat.ps==1, ...
        'WeightFunction','fe_nodal','Verbose',false);
    [c,d,A2]=SIF_LEFM_interaction_EDI( ...
        mesh,UII,V,mat,dom,'UsePlaneStrain',mat.ps==1, ...
        'WeightFunction','fe_nodal','Verbose',false);
    [e,f,AR]=SIF_LEFM_interaction_EDI( ...
        mesh,O25.U,V,mat,dom,'UsePlaneStrain',mat.ps==1, ...
        'WeightFunction','fe_nodal','Verbose',false);

    M(:,:,ir)=[a,c;b,d];
    Kactual(:,ir)=[e;f];
    misfit=norm(M(:,:,ir)-eye(2),'fro');
    condM=cond(M(:,:,ir));
    if condM<10 && isfinite(condM)
        calibrated(:,ir)=M(:,:,ir)\Kactual(:,ir);
    end
    rows(ir,:)=[rat(ir),ri/ro,a,b,c,d, ...
        misfit,condM,e,f,f/e,calibrated(1,ir), ...
        calibrated(2,ir),A1.nGP_used,AR.nGP_used];

    fprintf([' r_o/a0=%.2f | exact I=[%+.7e %+.7e], ', ...
        'exact II=[%+.7e %+.7e], ||M-I||F=%.3e, cond(M)=%.3g\n'], ...
        rat(ir),a,b,c,d,misfit,condM);
    fprintf('           actual=[%.8e,%+.8e], ratio=%+.6e, ', ...
        'M-inv diagnostic=[%.8e,%+.8e]\n'], ...
        e,f,f/e,calibrated(1,ir),calibrated(2,ir));
end
T=array2table(rows,'VariableNames',{ ...
    'r_outer_over_a0','r_inner_over_outer', ...
    'exactI_KI','exactI_KII_leakage', ...
    'exactII_KI_leakage','exactII_KII', ...
    'recovery_matrix_error','matrix_condition', ...
    'actual_KI','actual_KII','actual_KII_over_KI', ...
    'M_inverse_diagnostic_KI','M_inverse_diagnostic_KII', ...
    'exactI_GP_used','actual_GP_used'});
fprintf('\nSAME-MESH WILLIAMS RECOVERY TEST\n');disp(T);
fprintf(['Acceptance: small off-diagonal leakage and M approximately I ', ...
    'for every domain. Do not apply M-inverse diagnostic as a ', ...
    'production SIF correction without separate validation.\n']);

O26=struct('settings',O,'table',T,'recoveryMatrix',M, ...
    'K_actual',Kactual,'K_diagnostic',calibrated, ...
    'UI_exact',UI,'UII_exact',UII,'crack_face_side',face, ...
    'crack_axis',axis,'crack_tip',tip,'rOuterOverA0',rat, ...
    'maxMatrixError',max(squeeze(sum(sum((M-repmat(eye(2),[1,1,nR])).^2,1),2)).^0.5));
if O26.maxMatrixError>O.MaxMatrixError
    fprintf(['WARNING: max matrix error %.3e exceeds test threshold %.3e; ', ...
        'do not interpret physical KII yet.\n'], ...
        O26.maxMatrixError,O.MaxMatrixError);
end
fprintf('STEP 26 completed without any new FEM solve.\n');
end

function U=make_exact_displacements(xl,r,th,mat,KI,KII)
mu=mat.E/(2*(1+mat.nu));
if mat.ps==1
    kappa=3-4*mat.nu;
else
    kappa=(3-mat.nu)/(1+mat.nu);
end
c=cos(th/2);s=sin(th/2);
fac=sqrt(r/(2*pi))/(2*mu);
u1=fac.*(KI.*c.*(kappa-1+2*s.^2) ...
        +KII.*s.*(kappa+1+2*c.^2));
u2=fac.*(KI.*s.*(kappa+1-2*c.^2) ...
        -KII.*c.*(kappa-1-2*s.^2));
% Return global two-component interleaved U.
% xl is in the crack frame and the crack axis was fixed before this call.
% Caller rotates the local displacement by [axis; perpendicular].
% Store rotation in a separate immutable input instead of inferring from xl.
U=[u1,u2];  %#ok<NASGU>
error('step26:MissingRotation','Internal rotation argument is required.');
end

function must(S,f)
if ~isstruct(S)||~isfield(S,f)||isempty(S.(f))
    error('step26:MissingField','Required O25 field %s missing.',f);
end
end
