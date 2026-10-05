function O31=main_step31_full_annulus_mesh_audit(O22,O25,O27,O30,varargin)
%MAIN_STEP31_FULL_ANNULUS_MESH_AUDIT
% Hold Stage-I mouth, crack direction/length, hole polygon (NArc=480),
% temporary mouth half-shift and NOMINAL LOCAL HMIN fixed. Vary the global
% HMAX that controls the mesh outside the immediate crack-tip region.
%
% The existing Step-25 mesh/solution (global HMAX factor 1) and Step-27
% 16-point EDI values are reused. Default factor 2 creates exactly ONE
% additional cracked FE mesh/solution; optional factor 4 is a later step.
% For every field, report actual mesh scales ACROSS each EDI annulus
% (centroid-selected T3 triangles), use the same fixed FE-nodal-q EDI
% radii and 16-point quadrature, and extract native-face COD intercepts
% from the same displacement field using linear/quadratic fits.
%
% This is a combined global mesh/remeshing convergence experiment, not
% exact nested-mesh refinement. Actual tip scales and annulus mesh
% statistics must be checked. An EDI/COD disagreement at tiny signed KII
% is NOT proof of a physical finite-length kink.
%
% Examples:
%   O31=main_step31_full_annulus_mesh_audit(O22,O25,O27,O30);
%   O31=main_step31_full_annulus_mesh_audit(O22,O25,O27,O30,...
%          'HmaxFactors',[1 2 4],'Prior',O31); % reuse factor=2

p=inputParser;
addParameter(p,'HmaxFactors',[1 2], ...
    @(x)isnumeric(x)&&isvector(x)&&all(isfinite(x))&& ...
        all(x>=1)&&any(abs(x-1)<1e-12));
addParameter(p,'FitWindows',[.04 .20;.04 .30;.08 .30], ...
    @(x)isnumeric(x)&&size(x,2)==2&&all(isfinite(x(:)))&& ...
        all(x(:,1)>0)&all(x(:,2)<0.8)&all(x(:,1)<x(:,2)));
addParameter(p,'FitDegrees',[1 2], ...
    @(x)isnumeric(x)&&isvector(x)&&all(ismember(x,[1 2])));
addParameter(p,'MinFitPoints',8, ...
    @(x)isnumeric(x)&&isscalar(x)&&x>=5&&x==round(x));
addParameter(p,'Plot',false,@(x)islogical(x)||isnumeric(x));
addParameter(p,'Prior',[],@(x)isempty(x)||(isstruct(x)&&isscalar(x)));
parse(p,varargin{:});
opt=p.Results;
factors=sort(unique(opt.HmaxFactors(:).'));
windows=opt.FitWindows;
degrees=sort(unique(opt.FitDegrees(:).'));
nF=numel(factors);
rRat=O27.rOuterOverA0(:).';
nR=numel(rRat);

required25={'config','mesh','U','mat','crack','hTip','settings','mouth'};
for j=1:numel(required25),must(O25,required25{j},'O25');end
must(O22,'Stage1','O22');
must(O22.Stage1,'I','O22.Stage1');
must(O27,'Kactual','O27');
must(O27,'rules','O27');
must(O27,'settings','O27');
must(O30,'table','O30');
q16=find(O27.rules==16,1);
if isempty(q16)
    error('step31:NeedEDI16','O27 must include the 16-point EDI results.');
end
if ~isfield(O25.settings,'NArc')
    error('step31:NeedNArc','O25 must record the Stage-II NArc value.');
end
Cbase=O25.config;
I=O22.Stage1.I;
if norm(O25.mouth-I.x_star)>1e-11
    error('step31:MouthMismatch','O25 and O22 have different crack mouths.');
end
if ~isfield(O27.settings,'CommonInner')
    error('step31:NeedCommonInner','O27 must record the fixed inner radius.');
end
ri=O27.settings.CommonInner;
a0=norm(O25.crack.Pmid(end,:)-O25.crack.Pmid(1,:));
if ~(ri>=2*O25.hTip && ri<min(rRat)*a0)
    error('step31:BadInnerRadius','Existing common EDI inner radius is invalid.');
end
NArc=O25.settings.NArc;
baseHmax=Cbase.mesh1.hmax;
baseHmin=Cbase.mesh1.hmin;
mouthWidth=Cbase.mesh2.chw;
prior=opt.Prior;
if ~isempty(prior)
    fields={'factors','Cases','KI','KII','rOuterOverA0','rInner'};
    for j=1:numel(fields)
        must(prior,fields{j},'Prior');
    end
    if numel(prior.Cases)~=numel(prior.factors) || ...
            size(prior.KI,1)~=numel(prior.factors) || ...
            size(prior.KII,1)~=numel(prior.factors) || ...
            size(prior.KI,2)~=nR || size(prior.KII,2)~=nR || ...
            numel(prior.rOuterOverA0)~=nR || ...
            any(abs(prior.rOuterOverA0(:)-rRat(:))>1e-12) || ...
            abs(prior.rInner-ri)>1e-12
        error('step31:PriorIncompatible', ...
            'Prior result has incompatible factor count or EDI settings.');
    end
    for jj=1:numel(prior.Cases)
        Pcase=prior.Cases{jj};
        required={'mesh','U','mat','crack','Hmax','hTip'};
        for jk=1:numel(required)
            must(Pcase,required{jk}, ...
                sprintf('Prior.Cases{%d}',jj));
        end
        if abs(Pcase.Hmax-baseHmax/prior.factors(jj))>1e-12 || ...
                norm(Pcase.crack.Pmid(1,:)-I.x_star)>1e-11 || ...
                size(Pcase.crack.Pmid,1)~=2 || ...
                abs(norm(Pcase.crack.Pmid(end,:)- ...
                    Pcase.crack.Pmid(1,:))-a0)>1e-10 || ...
                numel(Pcase.U)~=2*size(Pcase.mesh.coord,1)
            error('step31:PriorGeometry', ...
                'Prior case %d is inconsistent with fixed mouth/mesh.',jj);
        end
    end
    fprintf('  PRIOR case cache validated (%d factors).\n', ...
        numel(prior.factors));
end
if ~(baseHmax>baseHmin)
    error('step31:BadH','Hmax must be greater than Hmin in baseline.');
end

KI=nan(nF,nR);
KII=nan(nF,nR);
hTip=nan(nF,1);
nodeCount=nan(nF,1);
annulusN=nan(nF,nR);
annulusMedianEdge=nan(nF,nR);
annulusP90Edge=nan(nF,nR);
annulusMaxEdge=nan(nF,nR);
Cases=cell(nF,1);
fitRows=nan(nF*size(windows,1)*numel(degrees),11);
row=0;

fprintf('\n============================================================\n');
fprintf('STEP 31: FIXED-MOUTH FULL-ANNULUS SPATIAL MESH AUDIT\n');
fprintf('============================================================\n');
fprintf('  fixed phi=%+.9f deg, mouth=[%.10e,%.10e] m\n', ...
    mod(rad2deg(I.phi_star)+180,360)-180,I.x_star(1),I.x_star(2));
fprintf('  a0=%.7g m, NArc=%d, appendix mouth shift=%.7g m\n', ...
    a0,NArc,mouthWidth);
fprintf('  fixed nominal Hmin=%.7g m, Hgrad=%.3f, r_inner=%.7g m\n', ...
    baseHmin,Cbase.mesh1.hgrad,ri);
fprintf('  global Hmax base=%.7g m, factors=%s\n', ...
    baseHmax,mat2str(factors));

for jf=1:nF
    facH=factors(jf);
    hmax=baseHmax/facH;
    if hmax<=1.02*baseHmin
        error('step31:TooSmallHmax', ...
            'Factor %.3f would make Hmax<=1.02*fixed Hmin.',facH);
    end
    fprintf('\n--- Hmax factor %.3f: Hmax=%.8e m ---\n',facH,hmax);

    if abs(facH-1)<1e-12
        % Exactly the same displacement field as the earlier audited run.
        Smesh=O25.mesh;
        U=O25.U;
        mat=O25.mat;
        crack=O25.crack;
        H.median=O25.hTip;
        if ~isfield(O27,'rOuterOverA0') || ...
                any(abs(O27.rOuterOverA0(:)-rRat(:))>1e-12)
            error('step31:CacheMismatch','O27 EDI radius data are inconsistent.');
        end
        for ir=1:nR
            KI(jf,ir)=O27.Kactual(1,q16,ir);
            KII(jf,ir)=O27.Kactual(2,q16,ir);
        end
        source='O25 FEM / O27 EDI16 cache';
    elseif ~isempty(prior) && any(abs(prior.factors-facH)<1e-12)
        % Reuse an ALREADY SOLVED full-domain mesh and EDI outputs.
        % This allows the factor-4 study without repeating factor 2.
        kPrior=find(abs(prior.factors-facH)<1e-12,1);
        old=prior.Cases{kPrior};
        Smesh=old.mesh;
        U=old.U;
        mat=old.mat;
        crack=old.crack;
        H=tip_edge_stats(Smesh,crack.Pmid(end,:));
        if abs(H.median-old.hTip)>1e-10
            error('step31:PriorTipMeshChanged', ...
                'Cached median tip-edge scale is not reproducible.');
        end
        KI(jf,:)=prior.KI(kPrior,:);
        KII(jf,:)=prior.KII(kPrior,:);
        source=sprintf('cached full FEM + EDI16, prior factor=%g',facH);
        fprintf('  REUSING prior factor %.3f (no new FEM solve).\n', ...
            facH);
    else
        C=Cbase;
        C.a0=a0;
        C.mesh1.hmax=hmax;
        C.mesh2.hmax=hmax;
        % Preserve the Step-25 fine tip scale and all other configuration.
        C.mesh1.hmin=baseHmin;
        C.mesh1.hhole=Cbase.mesh1.hhole;
        C.mesh1.hgrad=Cbase.mesh1.hgrad;
        C.mesh2.hcrack=Cbase.mesh2.hcrack;
        C.mesh2.hhole=Cbase.mesh2.hhole;
        C.mesh2.chw=mouthWidth;
        C.solver.verbose=0;

        [~,D,~,Mc]=build_stage2_cracked_mesh_for_theta(C,I,0, ...
            'NArc',NArc,'PlotGeom',false, ...
            'PlotMesh',false,'PlotCollapsed',logical(opt.Plot)&&jf==nF);
        if norm(D.Pmid(1,:)-I.x_star)>1e-11 || ...
                abs(norm(D.Pmid(end,:)-D.Pmid(1,:))-a0)>1e-10
            error('step31:GeometryMismatch', ...
                'Fixed Stage-I mouth or short crack length changed.');
        end
        S=solve_cracked_LEFM(C,Mc,'lambda',1.0);
        Smesh=S.mesh;
        U=S.U;
        mat=S.mat;
        crack=Mc.crack;
        H=tip_edge_stats(Smesh,crack.Pmid(end,:));
        if ri<2*H.median
            warning('step31:InnerNearTip', ...
                'Common inner radius is under 2*new tip h: ri/htip=%.2f.', ...
                ri/H.median);
        end
        for ir=1:nR
            [KI(jf,ir),KII(jf,ir)]=SIF_LEFM_interaction_EDI( ...
                Smesh,U,crack.Pmid,mat, ...
                struct('r_inner',ri,'r_outer',rRat(ir)*a0), ...
                'UsePlaneStrain',mat.ps==1, ...
                'Verbose',false,'WeightFunction','fe_nodal', ...
                'QuadratureRule',16);
        end
        source='new FEM / fresh EDI16';
    end
    hTip(jf)=H.median;
    nodeCount(jf)=size(Smesh.coord,1);
    fprintf('  T6 nodes=%d; actual tip-edge median=%.8e m\n', ...
        nodeCount(jf),H.median);
    Cases{jf}=struct('mesh',Smesh,'U',U,'mat',mat, ...
        'crack',crack,'source',source,'Hmax',hmax, ...
        'hTip',H.median);
    for ir=1:nR
        [ne,hm,hp,hx]=mesh_annulus_stats( ...
            Smesh,crack.Pmid(end,:),ri,rRat(ir)*a0);
        annulusN(jf,ir)=ne;
        annulusMedianEdge(jf,ir)=hm;
        annulusP90Edge(jf,ir)=hp;
        annulusMaxEdge(jf,ir)=hx;
        fprintf([' r_o/a0=%.2f | annulus nT3=%d, ', ...
            'h_median=%.6e, h_p90=%.6e | ', ...
            'KI=%.8e KII=%+.8e ratio=%+.7e\n'], ...
            rRat(ir),ne,hm,hp,KI(jf,ir),KII(jf,ir), ...
            KII(jf,ir)/KI(jf,ir));
    end

    % Evaluate COD ONLY at original upper-face T6 node distances. If
    % opposite-face distances differ, interpolate only that opposite face.
    [rNative,apparent,diagFace]=native_face_COD( ...
        Smesh,U,mat,crack,8);
    rr=rNative/a0;
    fprintf(['  COD native faces: up=%d, lo=%d, ', ...
        'grid mismatch=%.4e m\n'], ...
        diagFace.nUpper,diagFace.nLower,diagFace.gridMismatch);
    for iw=1:size(windows,1)
        idx=find(rr>=windows(iw,1) & rr<=windows(iw,2));
        for id=1:numel(degrees)
            row=row+1;
            degree=degrees(id);
            if numel(idx)<max(opt.MinFitPoints,2*(degree+1))
                fprintf('  COD window [%.2f,%.2f] degree=%d SKIP (only %d native nodes)\n', ...
                    windows(iw,1),windows(iw,2),degree,numel(idx));
                continue
            end
            polyI=polyfit(rr(idx),apparent(idx,1),degree);
            polyII=polyfit(rr(idx),apparent(idx,2),degree);
            k1=polyI(end);
            k2=polyII(end);
            % The factor-1 native-COD implementation must reproduce the
            % already completed Step-30 upper-native fits BEFORE we
            % compare any newly solved FEM fields.
            if abs(facH-1)<1e-12
                old=O30.table;
                jOld=find(old.side_1upper_2lower==1 & ...
                    abs(old.lower_r_over_a0-windows(iw,1))<1e-12 & ...
                    abs(old.upper_r_over_a0-windows(iw,2))<1e-12 & ...
                    old.polynomial_degree==degree,1);
                if isempty(jOld)
                    error('step31:MissingCODBaseline', ...
                        'Step-30 baseline missing window [%.2f,%.2f] degree %d.', ...
                        windows(iw,1),windows(iw,2),degree);
                end
                if abs(k1-old.KI_COD(jOld))>1e-9 || ...
                        abs(k2-old.KII_COD(jOld))>1e-10
                    error('step31:CODBaselineMismatch', ...
                        ['Baseline native COD mismatch at [%.2f,%.2f], ', ...
                         'degree %d: ΔKI=%g, ΔKII=%g.'], ...
                         windows(iw,1),windows(iw,2),degree, ...
                         k1-old.KI_COD(jOld),k2-old.KII_COD(jOld));
                end
            end
            fitRows(row,:)=[facH,hmax,H.median, ...
                windows(iw,1),windows(iw,2), ...
                degree,numel(idx),k1,k2,k2/k1, ...
                diagFace.gridMismatch];
            fprintf(['  COD [%4.2f,%4.2f] degree=%d, n=%d | ', ...
                'KI=%.8e KII=%+.8e ratio=%+.7e\n'], ...
                windows(iw,1),windows(iw,2), ...
                degree,numel(idx),k1,k2,k2/k1);
        end
    end
end

nRows=nF*nR;
rows=nan(nRows,11);
k=0;
for jf=1:nF
    for ir=1:nR
        k=k+1;
        rows(k,:)=[factors(jf),baseHmax/factors(jf), ...
            hTip(jf),nodeCount(jf),rRat(ir),annulusN(jf,ir), ...
            annulusMedianEdge(jf,ir),annulusP90Edge(jf,ir), ...
            KI(jf,ir),KII(jf,ir),KII(jf,ir)/KI(jf,ir)];
    end
end
T=array2table(rows,'VariableNames',{ ...
    'global_Hmax_refine_factor','Hmax','tip_h_actual','n_T6', ...
    'r_outer_over_a0','n_annulus_T3','annulus_h_median', ...
    'annulus_h_p90','KI_EDI16','KII_EDI16','ratio_EDI16'});
TCOD=array2table(fitRows(isfinite(fitRows(:,7)),:), ...
    'VariableNames',{'global_Hmax_refine_factor','Hmax','tip_h_actual', ...
    'fit_lower_r_over_a0','fit_upper_r_over_a0','degree', ...
    'n_native_face_nodes','KI_COD','KII_COD','ratio_COD', ...
    'face_grid_mismatch_m'});
S=nan(nF,8);
[~,refIdx]=min(abs(rRat-0.65));
for jf=1:nF
    q=KII(jf,:)./KI(jf,:);
    S(jf,:)=[factors(jf),baseHmax/factors(jf),hTip(jf), ...
        nodeCount(jf),KI(jf,refIdx),q(refIdx), ...
        max(q)-min(q),max(annulusP90Edge(jf,:))];
end
Summary=array2table(S,'VariableNames',{ ...
    'global_Hmax_refine_factor','Hmax','tip_h_actual','n_T6', ...
    'KI_EDI_ref','ratio_EDI_ref','EDI_ratio_domain_spread', ...
    'worst_annulus_h_p90'});

fprintf('\nFULL-ANNULUS EDI MESH STUDY\n');disp(T);
fprintf('\nNATIVE-FACE COD FIT STUDY\n');disp(TCOD);
fprintf('\nGLOBAL-HMAX CONVERGENCE SUMMARY\n');disp(Summary);
% Crucial spatial gate: do not confuse a lower prescribed global Hmax
% with an actual reduction of element sizes in the measured EDI annuli.
% Baseline is always factor 1, which is required by the input parser.
iBase=find(abs(factors-1)<1e-12,1);
refRows=nan(nF,7);
for jf=1:nF
    medRatio=annulusMedianEdge(jf,:)./annulusMedianEdge(iBase,:);
    p90Ratio=annulusP90Edge(jf,:)./annulusP90Edge(iBase,:);
    refRows(jf,:)=[factors(jf), ...
        medRatio(end),p90Ratio(end), ...
        max(medRatio),max(p90Ratio), ...
        min(annulusN(jf,:)./annulusN(iBase,:)), ...
        double(all(medRatio<0.90) && all(p90Ratio<0.90))];
end
Refinement=array2table(refRows,'VariableNames',{ ...
    'Hmax_factor','outer_annulus_median_vs_base', ...
    'outer_annulus_p90_vs_base','worst_median_ratio_vs_base', ...
    'worst_p90_ratio_vs_base','min_annulus_count_ratio_vs_base', ...
    'all_annuli_10pct_smaller'});
fprintf('\nMEASURED ANNULUS REFINEMENT CHECK\n');disp(Refinement);
fprintf(['A smaller input Hmax is NOT proof of local refinement; ', ...
    'inspect the measured ratios and native COD agreement.\n']);
fprintf(['Gate: a physical signed KII requires actual mesh resolution ', ...
    'across the EDI annuli, path-independent EDI, native COD ', ...
    'extrapolation stability, and subsequent polygon/mouth-width ', ...
    'convergence. Do not choose a kink from current residuals.\n']);

if logical(opt.Plot)
    figure('Name','Step31 fixed-geometry spatial convergence','Color','w');
    hold on;grid on;box on;
    for ir=1:nR
        plot(baseHmax./factors, KII(:,ir)./KI(:,ir),'-o', ...
            'DisplayName',sprintf('r_o/a_0=%.2f',rRat(ir)));
    end
    yline(0,'k:','HandleVisibility','off');
    xlabel('global Hmax [m]');
    ylabel('signed K_{II}/K_I');
    legend('Location','best');
end

% Do not embed the previous huge mesh/solution cache inside settings.
% The reusable fields are already returned once through Cases below.
optOut=opt;
optOut=rmfield(optOut,'Prior');
O31=struct('settings',optOut,'factors',factors,'rOuterOverA0',rRat, ...
    'rInner',ri,'KI',KI,'KII',KII,'hTip',hTip, ...
    'nNodes',nodeCount,'annulusN',annulusN, ...
    'annulusMedianEdge',annulusMedianEdge, ...
    'annulusP90Edge',annulusP90Edge, ...
    'annulusMaxEdge',annulusMaxEdge, ...
    'Cases',{Cases},'table',T,'CODTable',TCOD,'summary',Summary, ...
    'refinement',Refinement);
fprintf('STEP 31 completed.\n');
end

function [n,med,p90,mx]=mesh_annulus_stats(mesh,tip,ri,ro)
X=mesh.coord3;
T=mesh.connect3(:,1:3);
P1=X(T(:,1),:);P2=X(T(:,2),:);P3=X(T(:,3),:);
C=(P1+P2+P3)/3;
rc=hypot(C(:,1)-tip(1),C(:,2)-tip(2));
sel=(rc>=ri & rc<=ro);
h12=hypot(P1(:,1)-P2(:,1),P1(:,2)-P2(:,2));
h23=hypot(P2(:,1)-P3(:,1),P2(:,2)-P3(:,2));
h31=hypot(P3(:,1)-P1(:,1),P3(:,2)-P1(:,2));
h=max([h12 h23 h31],[],2);
h=sort(h(sel));
n=numel(h);
if n<1,med=NaN;p90=NaN;mx=NaN;return;end
med=median(h);
p90=h(max(1,ceil(.9*n)));
mx=h(end);
end

function H=tip_edge_stats(mesh,tip)
X=mesh.coord3;T=mesh.connect3(:,1:3);
d=hypot(X(:,1)-tip(1),X(:,2)-tip(2));
tol=max(1e-12,1e-8*max(1,max(abs(X(:)))));
ids=find(d<=min(d)+tol);
te=T(any(ismember(T,ids),2),:);
L=[];
for j=1:size(te,1)
    P=X(te(j,:),:);
    L=[L,norm(P(2,:)-P(1,:)), ...
        norm(P(3,:)-P(2,:)),norm(P(1,:)-P(3,:))]; %#ok<AGROW>
end
L=L(isfinite(L)&L>tol);
if isempty(L),error('step31:NoTipEdges','No tip-adjacent edges.');end
H=struct('median',median(L));
end

function [r,apparent,diag]=native_face_COD(mesh,U,mat,crack,minNative)
X=mesh.coord;T=mesh.connect;
n=size(X,1);
tip=crack.Pmid(end,:);
vec=crack.Pmid(end,:)-crack.Pmid(1,:);
a0=norm(vec);
vec=vec/a0;
perp=[-vec(2),vec(1)];
R=[vec(:),perp(:)];
xl=(X-tip)*R;
face=zeros(n,1);
tipID=crack.tipNode;
up=unique(crack.upperNodes(:));
lo=unique(crack.lowerNodes(:));
face(setdiff(up,[lo;tipID]))=+1;
face(setdiff(lo,[up;tipID]))=-1;
edgesMap=[1 2 4;2 3 5;3 1 6];
for j=1:3
    edge=T(:,edgesMap(j,:));
    v1=edge(:,1);v2=edge(:,2);
    mU=(face(v1)==1 & (face(v2)==1 | v2==tipID)) | ...
       (face(v2)==1 & (face(v1)==1 | v1==tipID));
    mL=(face(v1)==-1 & (face(v2)==-1 | v2==tipID)) | ...
       (face(v2)==-1 & (face(v1)==-1 | v1==tipID));
    mu=unique(edge(mU,3));
    ml=unique(edge(mL,3));
    if any(face(mu)==-1)||any(face(ml)==1)
        error('step31:ConflictingFaces', ...
            'Duplicate crack-face midside classification conflicted.');
    end
    face(mu)=1;
    face(ml)=-1;
end
tol=max(1e-12,1e-8*a0);
onFace=xl(:,1)<-tol & -xl(:,1)<=a0+tol & abs(xl(:,2))<tol;
up=find(onFace & face==1);
lo=find(onFace & face==-1);
if numel(up)<minNative||numel(lo)<minNative
    error('step31:TooFewCOD','Not enough classified native crack-face nodes.');
end
[rUp,iu]=sort(-xl(up,1));
[rLo,il]=sort(-xl(lo,1));
up=up(iu);lo=lo(il);
[rUp,~,gU]=unique(rUp);
[rLo,~,gL]=unique(rLo);
u=reshape(U,2,[]).'*R;
uU=zeros(numel(rUp),2);
uL=zeros(numel(rLo),2);
for k=1:2
    uU(:,k)=accumarray(gU,u(up,k),[],@mean);
    uL(:,k)=accumarray(gL,u(lo,k),[],@mean);
end
valid=(rUp>=min(rLo) & rUp<=max(rLo));
r=rUp(valid);
jump=uU(valid,:)-interp1(rLo,uL,r,'pchip');
if isempty(r)||any(~isfinite(jump(:)))
    error('step31:BadCOD','COD native-face interpolation failed.');
end
muMat=mat.E/(2*(1+mat.nu));
if mat.ps==1
    kappa=3-4*mat.nu;
else
    kappa=(3-mat.nu)/(1+mat.nu);
end
scale=muMat/(kappa+1)*sqrt(2*pi./r);
apparent=bsxfun(@times,[jump(:,2),jump(:,1)],scale);
mismatch=NaN;
if numel(rUp)==numel(rLo)
    mismatch=max(abs(rUp-rLo));
end
diag=struct('nUpper',numel(rUp),'nLower',numel(rLo), ...
    'gridMismatch',mismatch);
end

function must(S,key,label)
if ~isstruct(S)||~isfield(S,key)||isempty(S.(key))
    error('step31:MissingInput','Missing required %s.%s.',label,key);
end
end
