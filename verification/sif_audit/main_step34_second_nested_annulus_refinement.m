function O34=main_step34_second_nested_annulus_refinement(O25,O33,varargin)
%MAIN_STEP34_SECOND_NESTED_ANNULUS_REFINEMENT
% SECOND, successive spatially nested local refinement of the 8-mm crack.
% BASELINE: reuse the already solved Step-33 local FE mesh, material, COD
% and 16-point EDI results. No background remeshing or change of the
% mouth, crack path, hole polygon or outer plate boundary nodes.
% The NEW inner buffer defaults to zero, aiming to leave the already
% stable immediate tip-element scale untouched while further refining
% T3 elements actually inside the three tested integration annuli.
%
% Subdivide existing collapsed T3 elements in the crack-tip EDI
% neighborhood using the conforming topology-aware helper
% refine_collapsed_t3_annulus. Red/green edge-split templates retain
% distinct node IDs across the two coincident, traction-free crack faces.
% Interpolated T6 nodes are subsequently regenerated on this nested T3.
%
% ACTUAL spatial refinement is gated using measured annulus element
% edge lengths BEFORE starting the single new FEM solve.
%
% Compare three fixed EDI annuli using the SAME r_inner=0.0008 m,
% FE-nodal q and 16-point Dunavant quadrature, then compare independent
% native-face COD fits on the same solved displacement field.
%
% This test changes local spatial mesh resolution. Because some inner
% adjacent triangles may also split, the measured tip-scale change must
% be inspected. A tiny nonzero signed KII is not physical until local
% mesh, geometry polygon, and appendix width are independently stable.
%
% Usage:
%   O34mesh=main_step34_second_nested_annulus_refinement(O25,O33);
%   O34=main_step34_second_nested_annulus_refinement( ...
%       O25,O33,'DryRun',false); % only after reviewing dry run

p=inputParser;
addParameter(p,'TargetEdge',2.1e-4, ...
    @(x)isnumeric(x)&&isscalar(x)&&isfinite(x)&&x>0);
addParameter(p,'OuterBuffer',6e-4, ...
    @(x)isnumeric(x)&&isscalar(x)&&isfinite(x)&&x>=0);
addParameter(p,'InnerBuffer',0, ...
    @(x)isnumeric(x)&&isscalar(x)&&isfinite(x)&&x>=0);
addParameter(p,'MaxPasses',2, ...
    @(x)isnumeric(x)&&isscalar(x)&&isfinite(x)&&x>=1&&x<=5&&x==round(x));
addParameter(p,'MinImprovement',0.10, ...
    @(x)isnumeric(x)&&isscalar(x)&&isfinite(x)&&x>=0&&x<1);
addParameter(p,'RequireImprovement',true,@(x)islogical(x)||isnumeric(x));
addParameter(p,'FitWindows',[.04 .20;.04 .30;.08 .30], ...
    @(x)isnumeric(x)&&size(x,2)==2&&all(isfinite(x(:)))&& ...
    all(x(:,1)>0)&all(x(:,1)<x(:,2))&all(x(:,2)<.8));
addParameter(p,'FitDegrees',[1 2], ...
    @(x)isnumeric(x)&&isvector(x)&&all(ismember(x,[1 2])));
addParameter(p,'MinFitPoints',8, ...
    @(x)isnumeric(x)&&isscalar(x)&&isfinite(x)&&x>=6&&x==round(x));
addParameter(p,'Plot',false,@(x)islogical(x)||isnumeric(x));
addParameter(p,'DryRun',true,@(x)islogical(x)||isnumeric(x));
parse(p,varargin{:});
opt=p.Results;

for name={'config','mouth','crack','hTip'}
    must(O25,name{1},'O25');
end
for name={'mesh','U','mat','crack','hTipNew','originalFactor', ...
        'rInner','rOuterOverA0','newKI','newKII','CODTable', ...
        'refinementGatePassed'}
    must(O33,name{1},'O33');
end
if ~O33.refinementGatePassed
    error('step34:UngatedPrior', ...
        'The previous local refinement did not pass its mesh-size gate.');
end
if exist('refine_collapsed_t3_annulus','file')~=2
    error('step34:MissingRefiner', ...
        'Add verification/sif_audit to your MATLAB path.');
end
base=struct('mesh',O33.mesh,'U',O33.U,'mat',O33.mat, ...
    'crack',O33.crack, ...
    'Hmax',O25.config.mesh1.hmax/O33.originalFactor, ...
    'hTip',O33.hTipNew);
if ~isscalar(base.Hmax) || ~isfinite(base.Hmax) || base.Hmax<=0 || ...
        numel(base.U)~=2*size(base.mesh.coord,1)
    error('step34:BaselineInconsistent', ...
        'The previously solved Step-33 mesh and U are inconsistent.');
end
if norm(base.crack.Pmid(1,:)-O25.mouth)>1e-11
    error('step34:WrongMouth','Cached baseline moved the mouth.');
end
if size(base.crack.Pmid,1)~=2
    error('step34:NotStraight','Only a straight trial crack is audited.');
end
rRat=O33.rOuterOverA0(:).';
a0=norm(diff(base.crack.Pmid,1,1));
roMax=max(rRat)*a0;
ri=O33.rInner;
if ~(ri>0 && ri<min(rRat)*a0)
    error('step34:Annulus','Invalid existing common EDI inner radius.');
end
target=opt.TargetEdge;
if ~(target<max(O33.meshTable.new_h_p90))
    error('step34:TooCoarse', ...
        'TargetEdge must be below the baseline annulus p90 size.');
end
C=O25.config;
C.a0=a0;
C.mesh1.hmax=base.Hmax;
C.mesh2.hmax=base.Hmax;
C.solver.verbose=0;

fprintf('\n============================================================\n');
fprintf('STEP 34: SECOND NESTED ANNULAR SPATIAL REFINEMENT\n');
fprintf('============================================================\n');
fprintf('  mouth=[%.10e %.10e] m, a0=%.7g m, global Hmax factor=%.2f\n', ...
    O25.mouth(1),O25.mouth(2),a0,O33.originalFactor);
fprintf('  baseline T3=%d, T6=%d, r_inner=%.4g m, r_outer_max=%.4g m\n', ...
    size(base.mesh.connect3,1),size(base.mesh.coord,1),ri,roMax);
fprintf('  target max T3 edge=%.7g m, outer buffer=%.7g m\n', ...
    target,opt.OuterBuffer);
fprintf(['  All ORIGINAL physical vertices, hole boundary and ', ...
    'remote mesh boundary nodes are kept unchanged.\n']);

[P,T,crack,meshAudit]=refine_collapsed_t3_annulus( ...
    base.mesh.coord3,base.mesh.connect3,base.crack, ...
    ri,roMax,target, ...
    'OuterBuffer',opt.OuterBuffer, ...
    'InnerBuffer',opt.InnerBuffer, ...
    'MaxPasses',opt.MaxPasses,'Verbose',true);
if any(any(abs(P(1:size(base.mesh.coord3,1),:) - ...
        base.mesh.coord3)>1e-14))
    error('step34:VertexMoved','Original mesh vertex positions changed.');
end
% The refinement support is near the crack, nowhere near the outer
% plate boundaries: preserve exactly the previously applied far-field
% traction discretization and all old boundary geometry.
newP=P(size(base.mesh.coord3,1)+1:end,:);
boundTol=1e-11*max(1,max(abs(base.mesh.coord3(:))));
if any(abs(newP(:,1))<boundTol | abs(newP(:,1)-O25.config.A)<boundTol | ...
       abs(newP(:,2)-O25.config.B)<boundTol | ...
       abs(newP(:,2)+O25.config.B)<boundTol)
    error('step34:RemoteBoundaryTouched', ...
        'Local refinement unexpectedly split a loaded/outer plate edge.');
end
if norm(crack.Pmid-base.crack.Pmid,'fro')>1e-12 || ...
        crack.tipNode~=base.crack.tipNode
    error('step34:CrackMoved','Crack geometric path or tip ID changed.');
end
% Restrict this second refinement to the annulus, not the immediate tip.
% Test the ORIGINAL and REFINED actual tip-adjacent T3 geometry BEFORE
% the expensive FEM solve, including through the default dry-run path.
preTip=tip_edge_stats(struct('coord3',P,'connect3',T),crack.Pmid(end,:));
tipChange=abs(preTip/base.hTip-1);
fprintf('  pre-solve tip-edge median=%.10e m (relative change %.3e)\n', ...
    preTip,tipChange);
if tipChange>1e-9
    error('step34:TipNotPreserved', ...
        ['Second annular refinement changed the immediate tip scale; ', ...
         'reduce outer/inner refinement support before solving.']);
end

% Mesh statistics use exactly the same centroid-selected, longest-T3-edge
% definition as Step 31/32, so they permit quantitative comparisons.
nR=numel(rRat);
aOld=nan(nR,4);
aNew=nan(nR,4);
for ir=1:nR
    ro=rRat(ir)*a0;
    aOld(ir,:)=mesh_annulus_stats( ...
        base.mesh.coord3,base.mesh.connect3, ...
        crack.Pmid(end,:),ri,ro);
    aNew(ir,:)=mesh_annulus_stats( ...
        P,T,crack.Pmid(end,:),ri,ro);
end
improveMedian=1-aNew(:,2)./aOld(:,2);
improveP90=1-aNew(:,3)./aOld(:,3);
passed=all(isfinite(improveMedian)&isfinite(improveP90)) && ...
    all(improveMedian>=opt.MinImprovement) && ...
    all(improveP90>=opt.MinImprovement);
MeshTable=array2table([rRat(:),aOld,aNew, ...
    improveMedian,improveP90], ...
    'VariableNames',{ ...
    'r_outer_over_a0','old_nT3','old_h_median','old_h_p90', ...
    'old_h_max','new_nT3','new_h_median','new_h_p90','new_h_max', ...
    'median_reduction_fraction','p90_reduction_fraction'});

fprintf('\nMEASURED NESTED REFINEMENT (SAME LOCAL GEOMETRY)\n');
disp(MeshTable);
fprintf('  required minimum median AND p90 reduction: %.0f%%\n', ...
    100*opt.MinImprovement);
if ~passed
    warning('step34:InsufficientAnnulusRefinement', ...
        ['Measured nested annulus refinement fell below the gate. ', ...
         'Try a smaller TargetEdge or another MaxPasses.']);
    if logical(opt.RequireImprovement)
        O34=struct('settings',opt,'meshAudit',meshAudit, ...
            'meshTable',MeshTable,'refinementGatePassed',false, ...
            'stoppedBeforeSolve',true);
        fprintf(['  Stop: no FEM solve. Rerun with a smaller ', ...
            'TargetEdge or more MaxPasses.\n']);
        return
    end
end

% The first execution should be a geometry-only dry run, which tests
% the entire nested T3 refinement, oriented/area-preserving templates,
% crack-face topology, unchanged outer boundaries, and measured annular
% size reductions BEFORE any expensive FEM calculation.
if logical(opt.DryRun)
    O34=struct('settings',opt,'meshAudit',meshAudit, ...
        'meshTable',MeshTable,'refinementGatePassed',passed, ...
        'stoppedBeforeSolve',true, ...
        'nOldT3',size(base.mesh.connect3,1), ...
        'nNewT3',size(T,1), ...
        'nOldVertices',size(base.mesh.coord3,1), ...
        'nNewVertices',size(P,1), ...
        'hTipBaseline',base.hTip,'hTipRefinedPreSolve',preTip, ...
        'tipScalePreserved',tipChange<=1e-9, ...
        'crack',crack,'p',P,'t',T);
    fprintf(['  DRY RUN COMPLETED: no FE solve. After reviewing ', ...
        'these mesh/face diagnostics, rerun with DryRun=false.\n']);
    return
end

% All old corner node IDs remain exactly unchanged because new T3 nodes
% are appended after all original nodes; no boundary remeshing occurs.
X0=base.mesh.coord3;
[~,leftBottom]=min(sum((X0-[0,-C.B]).^2,2));
[~,rightBottom]=min(sum((X0-[C.A,-C.B]).^2,2));
G=struct();
G.p=P;
G.t=T;
G.edgeSets=struct();
G.edgeSets.corners=struct( ...
    'left_bottom',leftBottom,'right_bottom',rightBottom);
G.edgeSets.crackUpper=crack.upperNodes;
G.edgeSets.crackLower=crack.lowerNodes;
G.edgeSets.crackTip=crack.tipNode;
G.meta=struct('A',C.A,'B',C.B);

S=solve_cracked_LEFM(C,G,'lambda',1.0);
hNew=tip_edge_stats(S.mesh,crack.Pmid(end,:));
fprintf('  refined T6 nodes=%d; median tip edge=%.8e m ', ...
    size(S.mesh.coord,1),hNew);
fprintf('(baseline %.8e m)\n',base.hTip);
if abs(hNew/base.hTip-1)>.20
    warning('step34:TipScaleChanged', ...
        ['Tip-edge scale changed by >20%%. The study is nested local ', ...
         'refinement but does not isolate annulus from tip effects.']);
end

mat=S.mat;
if ~isfield(mat,'Dmat'),mat.Dmat=mat.D;end
newKI=nan(1,nR);newKII=nan(1,nR);
for ir=1:nR
    ro=rRat(ir)*a0;
    [newKI(ir),newKII(ir)]=SIF_LEFM_interaction_EDI( ...
        S.mesh,S.U,crack.Pmid,mat, ...
        struct('r_inner',ri,'r_outer',ro), ...
        'UsePlaneStrain',mat.ps==1, ...
        'Verbose',false,'WeightFunction','fe_nodal', ...
        'QuadratureRule',16);
    fprintf(['  16GP EDI r_o/a0=%.2f: KI=%.8e, ', ...
        'KII=%+.8e, KII/KI=%+.8e\n'], ...
        rRat(ir),newKI(ir),newKII(ir),newKII(ir)/newKI(ir));
end

% Re-evaluate native face COD on the NEW nested mesh from its explicit
% upper/lower crack-face T3 vertex lists; generated T6 midside nodes are
% topologically classified, with no coordinate merging across faces.
[r,apparent,faceInfo]=native_COD( ...
    S.mesh,S.U,S.mat,crack,opt.MinFitPoints);
rr=r/a0;
windows=opt.FitWindows;
degrees=sort(unique(opt.FitDegrees(:).'));
codRows=nan(size(windows,1)*numel(degrees),12);
k=0;
for iw=1:size(windows,1)
    idx=find(rr>=windows(iw,1)&rr<=windows(iw,2));
    for jd=1:numel(degrees)
        k=k+1;
        degree=degrees(jd);
        if numel(idx)<max(opt.MinFitPoints,2*(degree+1))
            fprintf('  COD window [%.2f,%.2f] degree=%d SKIP n=%d\n', ...
                windows(iw,1),windows(iw,2),degree,numel(idx));
            continue;
        end
        polI=polyfit(rr(idx),apparent(idx,1),degree);
        polII=polyfit(rr(idx),apparent(idx,2),degree);
        ki=polI(end);
        kii=polII(end);
        jOld=find( ...
            abs(O33.CODTable.lower_r_over_a0-windows(iw,1))<1e-12 & ...
            abs(O33.CODTable.upper_r_over_a0-windows(iw,2))<1e-12 & ...
            O33.CODTable.degree==degree,1);
        if isempty(jOld)
            error('step34:MissingPriorCOD', ...
                'O33 baseline COD window [%.2f,%.2f] degree %d missing.', ...
                windows(iw,1),windows(iw,2),degree);
        end
        codRows(k,:)=[windows(iw,:),degree,numel(idx), ...
            ki,kii,kii/ki, ...
            O33.CODTable.refined_KI(jOld), ...
            O33.CODTable.refined_KII(jOld), ...
            O33.CODTable.refined_ratio(jOld), ...
            ki-O33.CODTable.refined_KI(jOld), ...
            kii-O33.CODTable.refined_KII(jOld)];
        fprintf(['  COD [%.2f,%.2f] degree=%d n=%d ', ...
            '| refined KI=%.8e KII=%+.8e ratio=%+.6e\n'], ...
            windows(iw,1),windows(iw,2), ...
            degree,numel(idx),ki,kii,kii/ki);
    end
end
CODTable=array2table(codRows(isfinite(codRows(:,4)),:), ...
    'VariableNames',{ ...
    'lower_r_over_a0','upper_r_over_a0','degree','n_native_face', ...
    'refined_KI','refined_KII','refined_ratio', ...
    'baseline_KI','baseline_KII','baseline_ratio', ...
    'delta_KI','delta_KII'});

iRef=find(abs(rRat-.65)==min(abs(rRat-.65)),1);
qOld=O33.newKII./O33.newKI;
qNew=newKII./newKI;
EDIStudy=array2table([rRat(:), ...
    O33.newKI(:),O33.newKII(:),qOld(:), ...
    newKI(:),newKII(:),qNew(:),qNew(:)-qOld(:)], ...
    'VariableNames',{'r_outer_over_a0','baseline_KI','baseline_KII', ...
    'baseline_ratio','nested_KI','nested_KII','nested_ratio', ...
    'ratio_change'});

fprintf('\nNESTED ANNULAR MESH EDI COMPARISON\n');disp(EDIStudy);
fprintf('\nNESTED ANNULAR MESH COD COMPARISON\n');disp(CODTable);
fprintf(['  EDI reference ratio baseline=%+.8e, nested=%+.8e\n'], ...
    qOld(iRef),qNew(iRef));
fprintf(['  EDI domain spread baseline=%.7e, nested=%.7e\n'], ...
    max(qOld)-min(qOld),max(qNew)-min(qNew));
fprintf(['Interpretation: this test ACTUALLY subdivides elements around ', ...
    'the crack; only consider a finite kink if EDI/COD agree and ', ...
    'successive TARGETED mesh levels plus polygon/width tests ', ...
    'stabilize. Near-zero sign is not yet certified.\n']);

if logical(opt.Plot)
    figure('Name','Step33 paired annulus refinement','Color','w');
    hold on;box on;grid on;
    plot(rRat,qOld,'-o','LineWidth',1.2, ...
        'DisplayName','first nested local mesh');
    plot(rRat,qNew,'-s','LineWidth',1.2, ...
        'DisplayName','nested annulus refinement');
    yline(0,'k:','HandleVisibility','off');
    xlabel('r_{outer}/a_0');ylabel('signed K_{II}/K_I');
    legend('Location','best');
end

% Save only the new FE mesh/U and paired verification results.
O34=struct();
O34.settings=opt;
O34.meshAudit=meshAudit;
O34.meshTable=MeshTable;
O34.refinementGatePassed=passed;
O34.originalFactor=O33.originalFactor;
O34.oldKI=O33.newKI;
O34.oldKII=O33.newKII;
O34.newKI=newKI;
O34.newKII=newKII;
O34.rOuterOverA0=rRat;
O34.rInner=ri;
O34.hTipOld=base.hTip;
O34.hTipNew=hNew;
O34.tipScalePreserved=(tipChange<=1e-9);
O34.nT6Old=size(base.mesh.coord,1);
O34.nT6New=size(S.mesh.coord,1);
O34.mesh=S.mesh;
O34.U=S.U;
O34.mat=S.mat;
O34.crack=crack;
O34.faceInfo=faceInfo;
O34.nativeR=r;
O34.nativeApparent=apparent;
O34.EDIStudy=EDIStudy;
O34.CODTable=CODTable;
fprintf('STEP 34 completed: one new FE solve on second nested local refinement.\n');
end

function A=mesh_annulus_stats(P,T,tip,ri,ro)
tri=T(:,1:3);
p1=P(tri(:,1),:);p2=P(tri(:,2),:);p3=P(tri(:,3),:);
cen=(p1+p2+p3)/3;
r=hypot(cen(:,1)-tip(1),cen(:,2)-tip(2));
isAnnulus=r>=ri & r<=ro;
h12=hypot(p1(:,1)-p2(:,1),p1(:,2)-p2(:,2));
h23=hypot(p2(:,1)-p3(:,1),p2(:,2)-p3(:,2));
h31=hypot(p3(:,1)-p1(:,1),p3(:,2)-p1(:,2));
% Sort only longest-edge lengths in the selected EDI annulus.
localH=max([h12 h23 h31],[],2);
localH=sort(localH(isAnnulus));
n=numel(localH);
if n==0
    A=[0 NaN NaN NaN];
else
    A=[n,median(localH), ...
        localH(max(1,ceil(.9*n))),localH(end)];
end
end

function h=tip_edge_stats(mesh,tip)
P=mesh.coord3;T=mesh.connect3(:,1:3);
r=hypot(P(:,1)-tip(1),P(:,2)-tip(2));
tol=max(1e-12,1e-8*max(1,max(abs(P(:)))));
vertices=find(r<=min(r)+tol);
tc=T(any(ismember(T,vertices),2),:);
p1=P(tc(:,1),:);p2=P(tc(:,2),:);p3=P(tc(:,3),:);
lengths=[hypot(p1(:,1)-p2(:,1),p1(:,2)-p2(:,2)); ...
         hypot(p2(:,1)-p3(:,1),p2(:,2)-p3(:,2)); ...
         hypot(p3(:,1)-p1(:,1),p3(:,2)-p1(:,2))];
lengths=lengths(isfinite(lengths)&lengths>tol);
if isempty(lengths),error('step34:TipEdges','No tip-adjacent edges.');end
h=median(lengths);
end

function [r,app,diag]=native_COD(mesh,U,mat,crack,minPts)
X=mesh.coord;T=mesh.connect;
n=size(X,1);
tip=crack.Pmid(end,:);
vec=crack.Pmid(end,:)-crack.Pmid(1,:);
a0=norm(vec);vec=vec/a0;
perp=[-vec(2),vec(1)];
R=[vec(:),perp(:)];
xl=(X-tip)*R;
face=zeros(n,1);
tipID=crack.tipNode;
up=unique(crack.upperNodes(:));
lo=unique(crack.lowerNodes(:));
face(setdiff(up,[lo;tipID]))=1;
face(setdiff(lo,[up;tipID]))=-1;
emap=[1 2 4;2 3 5;3 1 6];
for j=1:3
    edge=T(:,emap(j,:));
    v1=edge(:,1);v2=edge(:,2);
    mU=(face(v1)==1&(face(v2)==1|v2==tipID)) | ...
       (face(v2)==1&(face(v1)==1|v1==tipID));
    mL=(face(v1)==-1&(face(v2)==-1|v2==tipID)) | ...
       (face(v2)==-1&(face(v1)==-1|v1==tipID));
    idsU=unique(edge(mU,3));idsL=unique(edge(mL,3));
    if any(face(idsU)==-1) || any(face(idsL)==1)
        error('step34:FaceConflict','Crack-face midside sets conflict.');
    end
    face(idsU)=+1;face(idsL)=-1;
end
tol=max(1e-12,1e-8*a0);
onFace=xl(:,1)<-tol & -xl(:,1)<=a0+tol & abs(xl(:,2))<tol;
up=find(onFace&face==1);lo=find(onFace&face==-1);
if numel(up)<minPts||numel(lo)<minPts
    error('step34:FaceNodes','Too few classified crack-face T6 nodes.');
end
[rU,iu]=sort(-xl(up,1));[rL,il]=sort(-xl(lo,1));
up=up(iu);lo=lo(il);
[rU,~,gU]=unique(rU);[rL,~,gL]=unique(rL);
u=reshape(U,2,[]).'*R;
Uu=zeros(numel(rU),2);Ul=zeros(numel(rL),2);
for k=1:2
    Uu(:,k)=accumarray(gU,u(up,k),[],@mean);
    Ul(:,k)=accumarray(gL,u(lo,k),[],@mean);
end
mask=rU>=min(rL)&rU<=max(rL);
r=rU(mask);
if isempty(r),error('step34:FaceOverlap','No face abscissa overlap.');end
jump=Uu(mask,:) - interp1(rL,Ul,r,'pchip');
if any(~isfinite(jump(:)))
    error('step34:NonfiniteCOD','Native COD interpolation failed.');
end
mu=mat.E/(2*(1+mat.nu));
if mat.ps==1
    kappa=3-4*mat.nu;
else
    kappa=(3-mat.nu)/(1+mat.nu);
end
scale=mu/(kappa+1)*sqrt(2*pi./r);
app=bsxfun(@times,[jump(:,2),jump(:,1)],scale);
mismatch=NaN;
if numel(rU)==numel(rL)
    mismatch=max(abs(rU-rL));
end
diag=struct('nUpper',numel(rU),'nLower',numel(rL), ...
    'gridMismatch',mismatch);
fprintf('  COD new mesh: native nodes upper/lower=%d/%d; ', ...
    diag.nUpper,diag.nLower);
fprintf('abscissa mismatch %.5e m\n',mismatch);
end

function must(S,field,label)
if ~isstruct(S)||~isfield(S,field)||isempty(S.(field))
    error('step34:MissingField','Missing %s.%s.',label,field);
end
end
