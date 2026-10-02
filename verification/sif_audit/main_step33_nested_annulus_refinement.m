function O33=main_step33_nested_annulus_refinement(O25,O32,varargin)
%MAIN_STEP33_NESTED_ANNULUS_REFINEMENT
% Spatially targeted, nested local refinement of the 8-mm cracked mesh.
% BASELINE: reuse the already solved Step-32 Hmax factor-2 case, including
% existing FE-nodal, 16-point EDI results. No background remeshing, no
% change of the mouth, crack path, polygon or remote boundary nodes.
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
%   O33=main_step33_nested_annulus_refinement(O25,O32);
%   O33=main_step33_nested_annulus_refinement(O25,O32, ...
%       'TargetEdge',3.5e-4,'MaxPasses',3);

p=inputParser;
addParameter(p,'BaselineHmaxFactor',2, ...
    @(x)isnumeric(x)&&isscalar(x)&&isfinite(x)&&x>=1);
addParameter(p,'TargetEdge',3.8e-4, ...
    @(x)isnumeric(x)&&isscalar(x)&&isfinite(x)&&x>0);
addParameter(p,'OuterBuffer',6e-4, ...
    @(x)isnumeric(x)&&isscalar(x)&&isfinite(x)&&x>=0);
addParameter(p,'InnerBuffer',8e-4, ...
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
parse(p,varargin{:});
opt=p.Results;

for name={'config','mouth','crack','hTip'}
    must(O25,name{1},'O25');
end
for name={'Cases','factors','rInner','rOuterOverA0','KI','KII','CODTable'}
    must(O32,name{1},'O32');
end
if exist('refine_collapsed_t3_annulus','file')~=2
    error('step33:MissingRefiner', ...
        'Add verification/sif_audit to your MATLAB path.');
end
iBase=find(abs(O32.factors(:)-opt.BaselineHmaxFactor)<1e-12,1);
if isempty(iBase)
    error('step33:MissingBaseline', ...
        'Baseline factor %.3f is unavailable in O32.', ...
        opt.BaselineHmaxFactor);
end
base=O32.Cases{iBase};
for name={'mesh','U','mat','crack','Hmax','hTip'}
    must(base,name{1},'O32.Cases{baseline}');
end
if norm(base.crack.Pmid(1,:)-O25.mouth)>1e-11
    error('step33:WrongMouth','Cached baseline moved the mouth.');
end
if size(base.crack.Pmid,1)~=2
    error('step33:NotStraight','Only a straight trial crack is audited.');
end
rRat=O32.rOuterOverA0(:).';
a0=norm(diff(base.crack.Pmid,1,1));
roMax=max(rRat)*a0;
ri=O32.rInner;
if ~(ri>0 && ri<min(rRat)*a0)
    error('step33:Annulus','Invalid existing common EDI inner radius.');
end
target=opt.TargetEdge;
if ~(target<max(O32.annulusP90Edge(iBase,:)))
    error('step33:TooCoarse', ...
        'TargetEdge must be below the baseline annulus p90 size.');
end
C=O25.config;
C.a0=a0;
C.mesh1.hmax=base.Hmax;
C.mesh2.hmax=base.Hmax;
C.solver.verbose=0;

fprintf('\n============================================================\n');
fprintf('STEP 33: NESTED LOCAL REFINEMENT OF THE SAME COLLAPSED MESH\n');
fprintf('============================================================\n');
fprintf('  mouth=[%.10e %.10e] m, a0=%.7g m, Hmax factor=%.2f\n', ...
    O25.mouth(1),O25.mouth(2),a0,opt.BaselineHmaxFactor);
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
    error('step33:VertexMoved','Original mesh vertex positions changed.');
end
if norm(crack.Pmid-base.crack.Pmid,'fro')>1e-12 || ...
        crack.tipNode~=base.crack.tipNode
    error('step33:CrackMoved','Crack geometric path or tip ID changed.');
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
    warning('step33:InsufficientAnnulusRefinement', ...
        ['Measured nested annulus refinement fell below the gate. ', ...
         'Try a smaller TargetEdge or another MaxPasses.']);
    if logical(opt.RequireImprovement)
        error('step33:StopBeforeSolve', ...
            'Local spatial refinement gate failed; FEM solve omitted.');
    end
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
    warning('step33:TipScaleChanged', ...
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
        jOld=find(abs(O32.CODTable.global_Hmax_refine_factor- ...
                opt.BaselineHmaxFactor)<1e-12 & ...
            abs(O32.CODTable.fit_lower_r_over_a0-windows(iw,1))<1e-12 & ...
            abs(O32.CODTable.fit_upper_r_over_a0-windows(iw,2))<1e-12 & ...
            O32.CODTable.degree==degree,1);
        if isempty(jOld)
            error('step33:MissingPriorCOD', ...
                'O32 baseline COD window [%.2f,%.2f] degree %d missing.', ...
                windows(iw,1),windows(iw,2),degree);
        end
        codRows(k,:)=[windows(iw,:),degree,numel(idx), ...
            ki,kii,kii/ki, ...
            O32.CODTable.KI_COD(jOld), ...
            O32.CODTable.KII_COD(jOld), ...
            O32.CODTable.ratio_COD(jOld), ...
            ki-O32.CODTable.KI_COD(jOld), ...
            kii-O32.CODTable.KII_COD(jOld)];
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
qOld=O32.KII(iBase,:)./O32.KI(iBase,:);
qNew=newKII./newKI;
EDIStudy=array2table([rRat(:), ...
    O32.KI(iBase,:).',O32.KII(iBase,:).',qOld(:), ...
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
        'DisplayName','original factor-2 mesh');
    plot(rRat,qNew,'-s','LineWidth',1.2, ...
        'DisplayName','nested annulus refinement');
    yline(0,'k:','HandleVisibility','off');
    xlabel('r_{outer}/a_0');ylabel('signed K_{II}/K_I');
    legend('Location','best');
end

% Save ONLY the new FE mesh, U and crack topology; do not retain S.K or
% the old O32.Cases (which may contain several large solved FEM fields).
O33=struct();
O33.settings=opt;
O33.meshAudit=meshAudit;
O33.meshTable=MeshTable;
O33.refinementGatePassed=passed;
O33.originalFactor=opt.BaselineHmaxFactor;
O33.oldKI=O32.KI(iBase,:);
O33.oldKII=O32.KII(iBase,:);
O33.newKI=newKI;
O33.newKII=newKII;
O33.rOuterOverA0=rRat;
O33.rInner=ri;
O33.hTipOld=base.hTip;
O33.hTipNew=hNew;
O33.nT6Old=size(base.mesh.coord,1);
O33.nT6New=size(S.mesh.coord,1);
O33.mesh=S.mesh;
O33.U=S.U;
O33.mat=S.mat;
O33.crack=crack;
O33.faceInfo=faceInfo;
O33.nativeR=r;
O33.nativeApparent=apparent;
O33.EDIStudy=EDIStudy;
O33.CODTable=CODTable;
fprintf('STEP 33 completed: one new FE solve on nested local refinement.\n');
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
h=sort(max([h12 h23 h31],[],2));
% Recompute sorted edge lengths within the selected annulus ONLY.
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
if isempty(lengths),error('step33:TipEdges','No tip-adjacent edges.');end
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
        error('step33:FaceConflict','Crack-face midside sets conflict.');
    end
    face(idsU)=+1;face(idsL)=-1;
end
tol=max(1e-12,1e-8*a0);
onFace=xl(:,1)<-tol & -xl(:,1)<=a0+tol & abs(xl(:,2))<tol;
up=find(onFace&face==1);lo=find(onFace&face==-1);
if numel(up)<minPts||numel(lo)<minPts
    error('step33:FaceNodes','Too few classified crack-face T6 nodes.');
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
if isempty(r),error('step33:FaceOverlap','No face abscissa overlap.');end
jump=Uu(mask,:) - interp1(rL,Ul,r,'pchip');
if any(~isfinite(jump(:)))
    error('step33:NonfiniteCOD','Native COD interpolation failed.');
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
    error('step33:MissingField','Missing %s.%s.',label,field);
end
end
