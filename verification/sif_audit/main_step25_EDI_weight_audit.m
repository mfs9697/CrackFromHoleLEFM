function O25=main_step25_EDI_weight_audit(O22,O24,varargin)
%MAIN_STEP25_EDI_WEIGHT_AUDIT
% Re-solve exactly ONE Step-24 factor-2 FE field and investigate why EDI
% radius sensitivity deteriorates despite refinement at the crack tip.
%
% The physical geometry and global Hmax are identical for every extractor
% call. Reuse the SAME S2.mesh and S2.U in all EDI evaluations. Compare:
%   A) FE-nodal q interpolation versus analytic radial grad(q);
%   B) legacy r_inner=max(0.1*r_outer,2*h_tip) versus COMMON r_inner;
%   C) integration-point counts and aux-field strain/stress consistency.
%
% The analytic q option is a diagnostic, NOT assumed more accurate: abrupt
% radial cutoffs require adequate Gauss-point resolution.
% An inner-radius effect or q-method disagreement implicates integration
%/weight construction; agreement of both methods with continued outer
% radius dependence leaves finite-element field/equilibrium errors possible.
%
% Usage:
%   O25=main_step25_EDI_weight_audit(O22,O24);

p=inputParser;
addParameter(p,'a0',0.008,@(v)isnumeric(v)&&isscalar(v)&&v>0);
addParameter(p,'NArc',480,@(v)isnumeric(v)&&isscalar(v)&&v>=32);
addParameter(p,'RefineFactor',2,@(v)isnumeric(v)&&isscalar(v)&&v>=1);
addParameter(p,'CommonInner',0.0008, ...
    @(v)isnumeric(v)&&isscalar(v)&&isfinite(v)&&v>0);
addParameter(p,'Methods',{'fe_nodal','analytic_radial'}, ...
    @(v)iscell(v)&&~isempty(v)&&all(cellfun(@ischar,v)));
parse(p,varargin{:});
opt=p.Results;
if nargin<2,error('step25:NeedResults','Provide both O22 and O24.');end

C=O22.Stage1.C;
I=O22.Stage1.I;
C.a0=opt.a0;
baseH=C.mesh1.hmin;
C.mesh1.hmin=baseH/opt.RefineFactor;
C.mesh1.hhole=C.mesh1.hmin;
C.mesh2.hcrack=C.mesh1.hmin;
C.mesh2.hhole=C.mesh1.hmin;
% IMPORTANT: global Hmax is unchanged, matching Step 24.
C.solver.verbose=0;
rat=O22.rOuterOverA0(:).';
nR=numel(rat);

fprintf('\n============================================================\n');
fprintf('STEP 25: SAME-FIELD SIGNED EDI WEIGHT AND INNER-RADIUS AUDIT\n');
fprintf('============================================================\n');
fprintf('  phi=%+.9f deg; a0=%.6g; NArc=%d; refinement factor=%.2f\n', ...
    rad2deg(local_wrap(I.phi_star)),C.a0,opt.NArc,opt.RefineFactor);
fprintf('  fixed common inner radius=%.6g m\n',opt.CommonInner);

[~,D,~,Mc]=build_stage2_cracked_mesh_for_theta(C,I,0, ...
    'NArc',opt.NArc,'PlotGeom',false,'PlotMesh',false, ...
    'PlotCollapsed',false);
if norm(D.Pmid(1,:)-I.x_star)>1e-11 || ...
        abs(norm(D.Pmid(end,:)-D.Pmid(1,:))-C.a0)>1e-10
    error('step25:GeometryMoved','Fixed mouth or a0 moved unexpectedly.');
end
S2=solve_cracked_LEFM(C,Mc,'lambda',1.0);
H=local_tip_mesh_scale(S2.mesh,Mc.crack.Pmid(end,:));
mat=S2.mat;
if ~isfield(mat,'Dmat'),mat.Dmat=mat.D;end
fprintf('  nT6nodes=%d, nT3=%d, median tip edge=%.8e m\n', ...
    size(S2.mesh.coord,1),size(S2.mesh.connect3,1),H.median);

regimes={'adaptive_inner','common_inner'};
methods=opt.Methods;
nM=numel(methods);
nA=numel(regimes);
rows=nan(nM*nA*nR,13);
K1=nan(nM,nA,nR);
K2=nan(nM,nA,nR);
GP=nan(nM,nA,nR);
k=0;

for jm=1:nM
    for jr=1:nA
        for ir=1:nR
            ro=rat(ir)*C.a0;
            if jr==1
                ri=max(0.1*ro,2*H.median);
            else
                ri=opt.CommonInner;
            end
            if ~(ri>=2*H.median && ri<ro)
                error('step25:BadAnnulus', ...
                    'Bad annulus with ro=%.6g and ri=%.6g, htip=%.6g.', ...
                    ro,ri,H.median);
            end

            [k1,k2,A]=SIF_LEFM_interaction_EDI( ...
                S2.mesh,S2.U,Mc.crack.Pmid,mat, ...
                struct('r_inner',ri,'r_outer',ro), ...
                'UsePlaneStrain',mat.ps==1, ...
                'Verbose',false,'WeightFunction',methods{jm});

            k=k+1;
            K1(jm,jr,ir)=k1;
            K2(jm,jr,ir)=k2;
            GP(jm,jr,ir)=A.nGP_used;
            rows(k,:)=[jm,jr,ro/C.a0,ri/ro, ...
                k1,k2,k2/k1,A.nGP_used,A.nElem_used, ...
                A.auxEpsMismatchI_median,A.auxEpsMismatchI_max, ...
                A.auxEpsMismatchII_median,A.auxEpsMismatchII_max];

            fprintf(['  %-15s | %-14s | ro/a0=%.2f ri/ro=%.3f | ', ...
                'KI=%.8e KII=%+.8e ratio=%+.6e | nGP=%d | ', ...
                'aux mismatch med I/II=[%.2e %.2e]\n'], ...
                methods{jm},regimes{jr},rat(ir),ri/ro,k1,k2,k2/k1, ...
                A.nGP_used,A.auxEpsMismatchI_median, ...
                A.auxEpsMismatchII_median);
        end
    end
end

T=array2table(rows,'VariableNames',{ ...
    'method_code','inner_regime_code','r_outer_over_a0','r_inner_over_outer', ...
    'KI','KII','KII_over_KI','nGP_used','nElem_used', ...
    'aux_mismatch_I_median','aux_mismatch_I_max', ...
    'aux_mismatch_II_median','aux_mismatch_II_max'});

srows=nan(nM*nA,6);
k=0;
[~,refidx]=min(abs(rat-0.65));
for jm=1:nM
    for jr=1:nA
        k=k+1;
        k1=squeeze(K1(jm,jr,:));
        q=squeeze(K2(jm,jr,:)./K1(jm,jr,:));
        srows(k,:)=[jm,jr,k1(refidx),q(refidx), ...
            max(q)-min(q),(max(k1)-min(k1))/abs(mean(k1))];
    end
end
Summary=array2table(srows,'VariableNames',{ ...
    'method_code','inner_regime_code','KI_reference', ...
    'KII_over_KI_reference','KII_over_KI_domain_spread','KI_relative_domain_spread'});
fprintf('\nSAME-DISPLACEMENT-FIELD EDI COMPARISON\n');disp(T);
fprintf('\nMETHOD/INNER-RADIUS DOMAIN SPREAD\n');disp(Summary);
fprintf('Method codes: ');
for jm=1:nM,fprintf('%d=%s ',jm,methods{jm});end
fprintf('\nInner-radius codes: 1=adaptive, 2=common.\n');

% Compare to Step 24 factor-2, same nominal field generation.
match=find(abs(O24.factors-opt.RefineFactor)<1e-10,1);
if isempty(match)
    warning('step25:NoStep24Match','No matching Step-24 refine factor.');
else
    [~,jFE]=find_method(methods,'fe_nodal');
    if ~isempty(jFE)
        oldRatio=O24.KII(match,:)./O24.KI(match,:);
        newRatio=squeeze(K2(jFE,1,:)./K1(jFE,1,:)).';
        fprintf('  Step24 versus Step25 SAME-CONFIG ratio max diff: %.6e\n', ...
            max(abs(oldRatio-newRatio)));
    end
end

O25=struct();
O25.config=C;
O25.mouth=I.x_star;
O25.methods=methods;
O25.innerRegimes=regimes;
O25.rOuterOverA0=rat;
O25.KI=K1;
O25.KII=K2;
O25.table=T;
O25.summary=Summary;
O25.mesh=S2.mesh;
O25.U=S2.U;
O25.mat=S2.mat;
O25.crack=Mc.crack;
O25.hTip=H.median;
O25.settings=opt;
fprintf('\nSTEP 25 completed; all methods used the same solved FEM field.\n');
end

function [found,idx]=find_method(c,needle)
idx=find(strcmpi(c,needle),1);
found=~isempty(idx);
end

function H=local_tip_mesh_scale(mesh,tip)
X=mesh.coord3;T=mesh.connect3;
d=hypot(X(:,1)-tip(1),X(:,2)-tip(2));
tol=max(1e-12,1e-8*max(1,max(abs(X(:)))));
hit=find(d<=min(d)+tol);
Te=T(any(ismember(T,hit),2),:);
L=[];
for j=1:size(Te,1)
    p=X(Te(j,:),:);
    L=[L,norm(p(2,:)-p(1,:)),norm(p(3,:)-p(2,:)), ...
       norm(p(1,:)-p(3,:))]; %#ok<AGROW>
end
L=L(isfinite(L)&L>tol);
if isempty(L),error('step25:NoTipEdges','No tip-edge lengths found.');end
H=struct('median',median(L));
end

function a=local_wrap(a)
a=mod(a+pi,2*pi)-pi;
end
