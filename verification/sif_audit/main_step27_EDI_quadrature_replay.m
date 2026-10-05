function O27=main_step27_EDI_quadrature_replay(O25,O26,varargin)
%MAIN_STEP27_EDI_QUADRATURE_REPLAY
% Same-mesh, same-field quadrature audit of signed FE-nodal EDI.
%
% Reuses BOTH original FEM displacement field (O25.U) and the two
% manufactured Williams fields (O26.UI_exact, O26.UII_exact).
% Also adds an EXACT AFFINE non-singular T-stress field, with local stress
% [sigma11,sigma22,sigma12]=[1,0,0].
%
% NO new FEM solve. Compare EDI 7 (historical degree-5) against 16
% (degree-8) Dunavant integration on the same T6 mesh, same radial q and
% fixed r_inner. Degree-6 12-point rule is available by passing [7 12 16].
% Quadrature-point counts and exact modal leakage are explicitly reported.
%
% WARNING: a Williams test using the same auxiliary field formulas as the
% EDI extractor is a consistency test, not an independent proof of accuracy.
% Also, T-stress=1 is synthetic, not an estimate of the actual BVP T-stress.
%
% Usage:
%   O27=main_step27_EDI_quadrature_replay(O25,O26);
%   O27=main_step27_EDI_quadrature_replay(O25,O26, ...
%         'Rules',[7 12 16]);

p=inputParser;
addParameter(p,'Rules',[7 16], ...
    @(x)isnumeric(x)&&isvector(x)&&all(ismember(x,[7 12 16])));
addParameter(p,'ROuterOverA0',O26.rOuterOverA0(:).', ...
    @(x)isnumeric(x)&&isvector(x)&&all(x>0)&&all(x<1));
addParameter(p,'CommonInner',O26.settings.CommonInner, ...
    @(x)isnumeric(x)&&isscalar(x)&&x>0);
addParameter(p,'UseCachedSeven',true,@(x)islogical(x)&&isscalar(x));
parse(p,varargin{:});
O=p.Results;
rules=sort(unique(O.Rules(:).'));
rr=O.ROuterOverA0(:).';
nRule=numel(rules);nR=numel(rr);

must(O25,'mesh');must(O25,'U');must(O25,'mat');must(O25,'crack');
must(O26,'UI_exact');must(O26,'UII_exact');
mesh=O25.mesh;
mat=O25.mat;
V=O25.crack.Pmid;
a0=norm(V(end,:)-V(1,:));
if ~isfield(mat,'Dmat'),mat.Dmat=mat.D;end

if numel(O26.UI_exact)~=numel(O25.U) || ...
        numel(O26.UII_exact)~=numel(O25.U)
    error('step27:BadExactField', ...
        'Exact Williams fields must use the exact same nodal DOFs as O25.U.');
end

% Manufactured nonsingular crack-tip T-stress: sigma_11=1, other
% traction components zero. For isotropic plane strain/stress this is a
% physical elastic affine displacement field that is exactly representable
% by T6 interpolation. The crack faces have zero sigma_22 and sigma_12.
axis=V(end,:)-V(end-1,:);
axis=axis/norm(axis);
Rgl=[axis(:),[-axis(2);axis(1)]];
xl=(mesh.coord-V(end,:))*Rgl;
epsT=mat.Dmat\[1;0;0];
if abs(epsT(3))>1e-12*norm(epsT)
    error('step27:UnexpectedShear','Manufactured T-stress has shear strain.');
end
ut=[epsT(1)*xl(:,1),epsT(2)*xl(:,2)];
UT=reshape((ut*Rgl.').',[],1);

nfields=4;
fieldNames={'exact_I','exact_II','actual_FEM','affine_T'};
M=nan(2,2,nRule,nR);
KActual=nan(2,nRule,nR);
KT=nan(2,nRule,nR);
GP=nan(nfields,nRule,nR);
Rows=nan(nRule*nR,15);
irow=0;

fprintf('\n============================================================\n');
fprintf('STEP 27: SAME-FIELD EDI QUADRATURE AND AFFINE T-STRESS AUDIT\n');
fprintf('============================================================\n');
fprintf('  a0=%.6g, fixed r_inner=%.7e m\n',a0,O.CommonInner);
fprintf('  T6 nodes=%d, triangles=%d, rules=%s\n', ...
    size(mesh.coord,1),size(mesh.connect,1),mat2str(rules));
fprintf('  Synthetic unit T-stress is NOT actual BVP T-stress.\n');

for iq=1:nRule
    nq=rules(iq);
    fprintf('\n--- %d-point Dunavant quadrature ---\n',nq);
    for ir=1:nR
        ro=rr(ir)*a0;
        ri=O.CommonInner;
        if ~(ri<ro && ri>=2*O25.hTip)
            error('step27:BadDomain', ...
                'For r_outer/a0=%.2f, r_inner=%.6g is invalid.', ...
                rr(ir),ri);
        end

        dom=struct('r_inner',ri,'r_outer',ro);
        canReuse=logical(O.UseCachedSeven) && nq==7 && ...
            abs(O26.settings.CommonInner-ri)<1e-12 && ...
            any(abs(O26.rOuterOverA0(:)-rr(ir))<1e-12);
        if canReuse
            iOld=find(abs(O26.rOuterOverA0(:)-rr(ir))<1e-12,1);
            M(:,:,iq,ir)=O26.recoveryMatrix(:,:,iOld);
            KActual(:,iq,ir)=O26.K_actual(:,iOld);
            GP(1,iq,ir)=O26.table.exactI_GP_used(iOld);
            GP(2,iq,ir)=GP(1,iq,ir);
            GP(3,iq,ir)=O26.table.actual_GP_used(iOld);
        else
            [a,b,A1]=SIF_LEFM_interaction_EDI( ...
                mesh,O26.UI_exact,V,mat,dom, ...
                'UsePlaneStrain',mat.ps==1,'Verbose',false, ...
                'WeightFunction','fe_nodal','QuadratureRule',nq);
            [c,d,A2]=SIF_LEFM_interaction_EDI( ...
                mesh,O26.UII_exact,V,mat,dom, ...
                'UsePlaneStrain',mat.ps==1,'Verbose',false, ...
                'WeightFunction','fe_nodal','QuadratureRule',nq);
            [e,f,A3]=SIF_LEFM_interaction_EDI( ...
                mesh,O25.U,V,mat,dom, ...
                'UsePlaneStrain',mat.ps==1,'Verbose',false, ...
                'WeightFunction','fe_nodal','QuadratureRule',nq);
            M(:,:,iq,ir)=[a c;b d];
            KActual(:,iq,ir)=[e;f];
            GP(1,iq,ir)=A1.nGP_used;
            GP(2,iq,ir)=A2.nGP_used;
            GP(3,iq,ir)=A3.nGP_used;
        end

        [kt1,kt2,AT]=SIF_LEFM_interaction_EDI( ...
            mesh,UT,V,mat,dom, ...
            'UsePlaneStrain',mat.ps==1,'Verbose',false, ...
            'WeightFunction','fe_nodal','QuadratureRule',nq);
        KT(:,iq,ir)=[kt1;kt2];
        GP(4,iq,ir)=AT.nGP_used;

        Mm=M(:,:,iq,ir);
        Ka=KActual(:,iq,ir);
        mErr=norm(Mm-eye(2),'fro');
        leak=Mm(2,1)/Mm(1,1);
        realRatio=Ka(2)/Ka(1);
        diffDiag=Ka(2)-Ka(1)*leak;

        irow=irow+1;
        Rows(irow,:)=[nq,rr(ir),ri/ro, ...
            Mm(1,1),Mm(2,1),Mm(1,2),Mm(2,2),mErr, ...
            Ka(1),Ka(2),realRatio, ...
            diffDiag,KT(1,iq,ir),KT(2,iq,ir),AT.nGP_used];

        fprintf([' r_o/a0=%.2f | ||M-I||F=%.4e ', ...
            'leak(II<-I)=%+.6e | actual ratio=%+.6e | ', ...
            'T=[%+.3e %+.3e] | GP=%d\n'], ...
            rr(ir),mErr,leak,realRatio,kt1,kt2,AT.nGP_used);
    end
end

T=array2table(Rows(1:irow,:),'VariableNames',{ ...
    'quadrature_rule','r_outer_over_a0','r_inner_over_outer', ...
    'exactI_KI','exactI_KII_leak','exactII_KI_leak','exactII_KII', ...
    'matrix_error','actual_KI','actual_KII','actual_KII_over_KI', ...
    'KI_leakage_adjusted_KII_diagnostic', ...
    'affine_T_KI','affine_T_KII','nGP_affine_T'});

S=nan(nRule,9);
for iq=1:nRule
    q=squeeze(KActual(2,iq,:)./KActual(1,iq,:));
    leak=squeeze(M(2,1,iq,:)./M(1,1,iq,:));
    qT=squeeze(KT(2,iq,:));
    tmp=zeros(nR,1);
    for ir=1:nR,tmp(ir)=norm(M(:,:,iq,ir)-eye(2),'fro');end
    [~,j]=min(abs(rr-0.65));
    S(iq,:)=[rules(iq),q(j),max(q)-min(q), ...
        max(tmp),max(abs(leak)),max(leak)-min(leak), ...
        max(abs(qT)),max(qT)-min(qT),...
        max(abs(squeeze(KT(1,iq,:))))];
end
Summary=array2table(S,'VariableNames',{ ...
    'quadrature_rule','actual_ratio_ref','actual_ratio_domain_spread', ...
    'max_recovery_matrix_error','max_abs_modeI_leakage', ...
    'modeI_leakage_domain_spread','max_abs_T_KII_per_unit_T', ...
    'T_KII_domain_spread_per_unit_T','max_abs_T_KI_per_unit_T'});

fprintf('\nCOMPLETE QUADRATURE AND MANUFACTURED FIELD RESULTS\n');
disp(T);
fprintf('\nQUADRATURE SENSITIVITY SUMMARY\n');
disp(Summary);
fprintf(['Interpretation: higher quadrature order should reduce exact-field ', ...
    'leakage AND improve domain consistency of actual FEM data if ', ...
    'Gauss integration was the leading source of error. ', ...
    'Residual T-field leakage tests nonsingular-field cancellation, ', ...
    'not the unknown T amplitude in the actual BVP.\n']);

O27=struct();
O27.settings=O;
O27.rules=rules;
O27.rOuterOverA0=rr;
O27.recoveryMatrix=M;
O27.Kactual=KActual;
O27.KaffineT=KT;
O27.gaussPoints=GP;
O27.affineT_U=UT;
O27.table=T;
O27.summary=Summary;
O27.inputO25=O25;
fprintf('STEP 27 completed: no new FEM solve.\n');
end

function must(S,f)
if ~isstruct(S)||~isfield(S,f)||isempty(S.(f))
    error('step27:MissingField','Missing required field %s.',f);
end
end
