function O28=main_step28_EDI_interpolation_audit(O25,O26,O27,varargin)
%MAIN_STEP28_EDI_INTERPOLATION_AUDIT
% Compare synthetic Williams fields as:
%   A) exactly sampled displacements interpolated by the existing T6 mesh;
%   B) exact Williams stress, strain and displacement gradient evaluated
%      DIRECTLY at each Gauss point (diagnostic AnalyticActualK option).
% Geometry, mesh, FE-nodal q, annulus and Dunavant rule are identical.
% Actual FEM-field results from O27 are included for context, but never
% corrected by this artificial-field recovery calibration.
%
% Requires existing O25, O26, O27. No new FEM solves.
% Example:
%    O28=main_step28_EDI_interpolation_audit(O25,O26,O27);

p=inputParser;
addParameter(p,'Rules',[7 16], ...
    @(x)isnumeric(x)&&isvector(x)&&all(isfinite(x))&&all(ismember(x,[7 12 16])));
addParameter(p,'CommonInner',O27.settings.CommonInner, ...
    @(x)isnumeric(x)&&isscalar(x)&&isfinite(x)&&x>0);
parse(p,varargin{:});
O=p.Results;

requiredO25={'mesh','U','mat','crack','hTip'};
requiredO26={'UI_exact','UII_exact'};
requiredO27={'rules','rOuterOverA0','recoveryMatrix','Kactual'};
for j=1:numel(requiredO25),must(O25,requiredO25{j},'O25');end
for j=1:numel(requiredO26),must(O26,requiredO26{j},'O26');end
for j=1:numel(requiredO27),must(O27,requiredO27{j},'O27');end
if abs(O.CommonInner-O27.settings.CommonInner)>1e-12
    error('step28:InnerMismatch', ...
        'CommonInner must equal O27.settings.CommonInner when using O27 caches.');
end

mesh=O25.mesh;
mat=O25.mat;
if ~isfield(mat,'Dmat'),mat.Dmat=mat.D;end
V=O25.crack.Pmid;
a0=norm(V(end,:)-V(1,:));
rat=O27.rOuterOverA0(:).';
rules=sort(unique(O.Rules(:).'));
nR=numel(rat);nQ=numel(rules);

MNodal=nan(2,2,nQ,nR);
MGauss=nan(2,2,nQ,nR);
KActual=nan(2,nQ,nR);
rows=nan(nQ*nR,15);
k=0;

fprintf('\n============================================================\n');
fprintf('STEP 28: WILLIAMS INTERPOLATION VS GAUSS-POINT RECOVERY\n');
fprintf('============================================================\n');
fprintf('  a0=%.8e | T6 nodes=%d | r_inner=%.7e m\n', ...
    a0,size(mesh.coord,1),O.CommonInner);
fprintf('  All comparisons: SAME mesh, FE-nodal q, annulus and quadrature.\n');

for iq=1:nQ
    npts=rules(iq);
    iCache=find(O27.rules==npts,1);
    for ir=1:nR
        ro=rat(ir)*a0;
        ri=O.CommonInner;
        if ri<2*O25.hTip || ri>=ro
            error('step28:InvalidAnnulus', ...
                'Need 2*h_tip<=r_inner<r_outer at r_o/a0=%.3f.',rat(ir));
        end
        dom=struct('r_inner',ri,'r_outer',ro);
        args={'UsePlaneStrain',mat.ps==1, ...
            'Verbose',false,'WeightFunction','fe_nodal', ...
            'QuadratureRule',npts};

        % Use O27 cached nodal recovery and actual-field extraction when
        % the rule and annuli are exactly the same.
        cacheOK=~isempty(iCache) && ...
            abs(O27.settings.CommonInner-ri)<1e-12;
        if cacheOK
            MNodal(:,:,iq,ir)=O27.recoveryMatrix(:,:,iCache,ir);
            KActual(:,iq,ir)=O27.Kactual(:,iCache,ir);
        else
            [k1,k2]=SIF_LEFM_interaction_EDI( ...
                mesh,O26.UI_exact,V,mat,dom,args{:});
            [k3,k4]=SIF_LEFM_interaction_EDI( ...
                mesh,O26.UII_exact,V,mat,dom,args{:});
            [ka,kb]=SIF_LEFM_interaction_EDI( ...
                mesh,O25.U,V,mat,dom,args{:});
            MNodal(:,:,iq,ir)=[k1,k3;k2,k4];
            KActual(:,iq,ir)=[ka;kb];
        end

        % Only this block is NEW: actual fields are evaluated at the
        % quadrature point using the same auxiliary Williams formula.
        % The recovered field no longer incurs T6 nodal interpolation
        % error but STILL uses the exact same FE-nodal q gradient.
        [g11,g21]=SIF_LEFM_interaction_EDI( ...
            mesh,O25.U,V,mat,dom,args{:},'AnalyticActualK',[1 0]);
        [g12,g22]=SIF_LEFM_interaction_EDI( ...
            mesh,O25.U,V,mat,dom,args{:},'AnalyticActualK',[0 1]);
        MGauss(:,:,iq,ir)=[g11,g12;g21,g22];

        A=MNodal(:,:,iq,ir);B=MGauss(:,:,iq,ir);
        act=KActual(:,iq,ir);
        errA=norm(A-eye(2),'fro');
        errB=norm(B-eye(2),'fro');
        leakA=A(2,1)/A(1,1);
        leakB=B(2,1)/B(1,1);
        k=k+1;
        rows(k,:)=[npts,rat(ir),ri/ro, ...
            errA,errB,leakA,leakB, ...
            A(1,1),B(1,1),A(2,2),B(2,2), ...
            act(1),act(2),act(2)/act(1), ...
            errA/max(errB,eps)];

        fprintf([' %2dGP r_o/a0=%.2f | nodal ||M-I||=%.4e, ', ...
            'exact-GP ||M-I||=%.4e | nodal leakage=%+.4e, ', ...
            'exact-GP leakage=%+.4e | actual KII/KI=%+.5e\n'], ...
            npts,rat(ir),errA,errB,leakA,leakB,act(2)/act(1));
    end
end

T=array2table(rows(1:k,:),'VariableNames',{ ...
    'quadrature_rule','r_outer_over_a0','r_inner_over_outer', ...
    'nodal_matrix_error','gauss_matrix_error', ...
    'nodal_modeI_leakage','gauss_modeI_leakage', ...
    'nodal_I_recovery','gauss_I_recovery', ...
    'nodal_II_recovery','gauss_II_recovery', ...
    'actual_KI','actual_KII','actual_KII_over_KI', ...
    'nodal_over_gauss_matrix_error_ratio'});

sr=nan(nQ,8);
for iq=1:nQ
    nerr=nan(1,nR);gerr=nan(1,nR);
    nl=nan(1,nR);gl=nan(1,nR);
    aq=squeeze(KActual(2,iq,:)./KActual(1,iq,:));
    for ir=1:nR
        nerr(ir)=norm(MNodal(:,:,iq,ir)-eye(2),'fro');
        gerr(ir)=norm(MGauss(:,:,iq,ir)-eye(2),'fro');
        nl(ir)=MNodal(2,1,iq,ir)/MNodal(1,1,iq,ir);
        gl(ir)=MGauss(2,1,iq,ir)/MGauss(1,1,iq,ir);
    end
    sr(iq,:)=[rules(iq),max(nerr),max(gerr),max(abs(nl)), ...
        max(abs(gl)),max(nl)-min(nl),max(gl)-min(gl),max(aq)-min(aq)];
end
Summary=array2table(sr,'VariableNames',{ ...
    'quadrature_rule','max_nodal_matrix_error','max_gauss_matrix_error', ...
    'max_abs_nodal_modeI_leakage','max_abs_gauss_modeI_leakage', ...
    'nodal_leakage_domain_spread','gauss_leakage_domain_spread', ...
    'actual_ratio_domain_spread'});

fprintf('\nWILLIAMS INTERPOLATION AUDIT RESULTS\n');disp(T);
fprintf('\nCOMPARISON SUMMARY\n');disp(Summary);
fprintf(['Interpretation: exact-GP vs nodal-Williams DIFFERENCE isolates ', ...
    'the synthetic field''s T6 nodal interpolation effect under ', ...
    'unchanged FE q and quadrature. Any exact-GP residual reflects ', ...
    'q interpolation/integration or other extractor errors. ', ...
    'It does NOT by itself establish the source of the actual FEM ', ...
    'displacement error.\n']);

O28=struct('settings',O,'rules',rules,'rOuterOverA0',rat, ...
    'MNodal',MNodal,'MGauss',MGauss,'Kactual',KActual, ...
    'table',T,'summary',Summary);
fprintf('STEP 28 completed; zero additional FEM solves.\n');
end

function must(S,f,label)
if ~isstruct(S)||~isfield(S,f)||isempty(S.(f))
    error('step28:MissingInput','Required %s field %s is missing.',label,f);
end
end
