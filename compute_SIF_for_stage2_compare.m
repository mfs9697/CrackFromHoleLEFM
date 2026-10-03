function Cmp = compute_SIF_for_stage2_compare(C,G2,Mc,S2,varargin)
%COMPUTE_SIF_FOR_STAGE2_COMPARE
% Side-by-side SIF extraction on ONE Stage-II FEM field.
%
% Historical extractor:
%   circular J-integral + Ishikawa--Kitagawa--Okamura mirror separation.
%
% Mesh-general extractor:
%   interaction equivalent-domain integral with FE-nodal q.
%
% No displacement solve is repeated here. Both methods receive exactly
% S2.mesh and S2.U.
%
% Options
%   OldRadius              default G2.tip.radiusJ
%   EDIOuterRadius         default OldRadius
%   EDIInnerRadius         explicit value; default automatic
%   EDIInnerRatio          default 0.10
%   EDIMinTipSizes         default 2.0
%   OldNtheta              default 240
%   AuxDerivativeScale     default 1.0
%
% Automatic EDI inner radius:
%   max(EDIInnerRatio*r_outer, EDIMinTipSizes*h_tip)
% where h_tip is estimated from T3 edge lengths of elements incident on the
% crack-tip corner node.

ip=inputParser;
addParameter(ip,'OldRadius',[],@(x)isempty(x)||(isnumeric(x)&&isscalar(x)&&x>0));
addParameter(ip,'EDIOuterRadius',[],@(x)isempty(x)||(isnumeric(x)&&isscalar(x)&&x>0));
addParameter(ip,'EDIInnerRadius',[],@(x)isempty(x)||(isnumeric(x)&&isscalar(x)&&x>=0));
addParameter(ip,'EDIInnerRatio',0.10,@(x)isnumeric(x)&&isscalar(x)&&x>0&&x<1);
addParameter(ip,'EDIMinTipSizes',2.0,@(x)isnumeric(x)&&isscalar(x)&&x>0);
addParameter(ip,'OldNtheta',240,@(x)isnumeric(x)&&isscalar(x)&&x>=40);
addParameter(ip,'AuxDerivativeScale',1.0,@(x)isnumeric(x)&&isscalar(x)&&x>0);
parse(ip,varargin{:});
O=ip.Results;

mesh=S2.mesh;
mat=S2.mat;
if ~isfield(mat,'Dmat')
    if isfield(mat,'D')
        mat.Dmat=mat.D;
    else
        error('compute_SIF_for_stage2_compare:MissingD','Material needs D or Dmat.');
    end
end

if isfield(Mc,'crack') && isfield(Mc.crack,'Pmid') && ~isempty(Mc.crack.Pmid)
    V=Mc.crack.Pmid;
else
    V=G2.crack.polyline;
end

if size(V,1)<2
    error('compute_SIF_for_stage2_compare:BadCrackPath','Need at least two crack points.');
end

Llast=norm(V(end,:)-V(end-1,:));
if Llast<=eps
    error('compute_SIF_for_stage2_compare:DegenerateLastLeg','Last crack leg is degenerate.');
end

rOld=O.OldRadius;
if isempty(rOld), rOld=G2.tip.radiusJ; end

rOut=O.EDIOuterRadius;
if isempty(rOut), rOut=rOld; end

H=local_tip_mesh_scale(mesh,V(end,:));

rIn=O.EDIInnerRadius;
if isempty(rIn)
    rIn=max(O.EDIInnerRatio*rOut,O.EDIMinTipSizes*H.median);
end

if ~(rIn>=0 && rIn<rOut)
    error('compute_SIF_for_stage2_compare:BadEDIDomain', ...
        ['Automatic/selected EDI domain is invalid: r_inner=%.6g, ', ...
         'r_outer=%.6g, h_tip(median)=%.6g. Refine the tip mesh or ', ...
         'increase the outer radius.'],rIn,rOut,H.median);
end

[KIold,KIIold,Dbg]=SIF_LEFM_circle2_debug( ...
    mesh,S2.U,V,mat,rOld, ...
    'nthet',O.OldNtheta,'plot',false,'verbose',false);

dom=struct('r_inner',rIn,'r_outer',rOut);
[KIedi,KIIedi,AE]=SIF_LEFM_interaction_EDI( ...
    mesh,S2.U,V,mat,dom, ...
    'UsePlaneStrain',mat.ps==1, ...
    'Verbose',false, ...
    'WeightFunction','fe_nodal', ...
    'AuxDerivativeScale',O.AuxDerivativeScale);

Knorm=max(hypot(KIedi,KIIedi),eps);

Cmp=struct();
Cmp.KI_old=KIold;
Cmp.KII_old=KIIold;
Cmp.KI_EDI=KIedi;
Cmp.KII_EDI=KIIedi;
Cmp.dKI_old_minus_EDI=KIold-KIedi;
Cmp.dKII_old_minus_EDI=KIIold-KIIedi;
Cmp.vector_difference_rel=hypot(KIold-KIedi,KIIold-KIIedi)/Knorm;

Cmp.old=Dbg;
Cmp.edi=AE;
Cmp.V=V;
Cmp.Llast=Llast;
Cmp.r_old=rOld;
Cmp.domain_EDI=dom;
Cmp.tipMeshScale=H;
Cmp.settings=O;

Cmp.stencil=struct();
Cmp.stencil.mirror_T3_mismatch_median=Dbg.diagnostics.mirrorT3Mismatch_median;
Cmp.stencil.mirror_T3_mismatch_p95=Dbg.diagnostics.mirrorT3Mismatch_p95;
Cmp.stencil.mirror_T3_mismatch_max=Dbg.diagnostics.mirrorT3Mismatch_max;
Cmp.stencil.fraction_exact_mirror_T3=Dbg.diagnostics.frac_exact_mirror_T3;
end


function H=local_tip_mesh_scale(mesh,tip)
if ~isfield(mesh,'coord3') || ~isfield(mesh,'connect3')
    error('compute_SIF_for_stage2_compare:MissingT3','Need coord3/connect3.');
end

X=mesh.coord3;
T=mesh.connect3;
d=sqrt(sum((X-tip).^2,2));
dmin=min(d);
tol=max(1e-12,1e-8*max(1,max(abs(X(:)))));
tipNodes=find(d<=dmin+tol);

hit=any(ismember(T,tipNodes),2);
Te=T(hit,:);
if isempty(Te)
    error('compute_SIF_for_stage2_compare:NoTipElements','No T3 tip elements found.');
end

L=[];
for k=1:size(Te,1)
    P=X(Te(k,:),:);
    lk=[norm(P(2,:)-P(1,:)),norm(P(3,:)-P(2,:)),norm(P(1,:)-P(3,:))];
    L=[L,lk]; %#ok<AGROW>
end
L=L(isfinite(L)&L>tol);
if isempty(L)
    error('compute_SIF_for_stage2_compare:NoTipEdges','No positive tip-edge lengths.');
end

H=struct();
H.min=min(L);
H.median=median(L);
H.max=max(L);
H.nEdges=numel(L);
H.nTipElements=size(Te,1);
H.tipNodeDistance=dmin;
end
