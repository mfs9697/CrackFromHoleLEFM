function O24=main_step24_fixed_geometry_FE_refinement(O22,varargin)
%MAIN_STEP24_FIXED_GEOMETRY_FE_REFINEMENT
% Focused verification of the 8-mm residual revealed in Step 23.
% Exactly the same Stage-I mouth, normal direction, short-crack length,
% retained-hole polygon (NArc=480) and appendix mouth width are used.
% Only nominal FE Hmin/Hhole/Hcrack are changed. Global Hmax and Hgrad
% remain fixed to separate local resolution from far-field refinement.
% A regenerated unstructured mesh is expected for every refinement;
% actual tip-edge scale is printed, not presumed to track nominal Hmin.
%
% The previous Step-23 MATLAB diagnostic error was after all solves,
% leaving O23 unassigned. This independent driver needs only O22.
%
% Examples:
% O24=main_step24_fixed_geometry_FE_refinement(O22);
% O24=main_step24_fixed_geometry_FE_refinement(O22,'RefineFactors',[1 2 4]);

p=inputParser;
addParameter(p,'a0',0.008,@(x)isnumeric(x)&&isscalar(x)&&x>0);
addParameter(p,'NArc',480,@(x)isnumeric(x)&&isscalar(x)&&x>=32&&x==round(x));
addParameter(p,'RefineFactors',[1 2], ...
    @(x)isnumeric(x)&&isvector(x)&&all(isfinite(x))&&all(x>=1));
addParameter(p,'Plot',false,@(x)islogical(x)||isnumeric(x));
parse(p,varargin{:});
O=p.Results;

if nargin<1||~isstruct(O22)||~isfield(O22,'Stage1')|| ...
        ~isfield(O22.Stage1,'I')
    error('step24:NeedO22','Pass existing Step-22 output O22.');
end
Cbase=O22.Stage1.C;
I=O22.Stage1.I;
rat=O22.rOuterOverA0(:).';
factors=sort(unique(O.RefineFactors(:).'));
nF=numel(factors);nR=numel(rat);
hBase=Cbase.mesh1.hmin;
hmaxFixed=Cbase.mesh1.hmax;
KI=nan(nF,nR);KII=nan(nF,nR);rIn=nan(nF,nR);
hTip=nan(nF,1);nNodes=nan(nF,1);nElem=nan(nF,1);
nSig=nan(nF,1);nTri=nan(nF,1);

fprintf('\n============================================================\n');
fprintf('STEP 24: FIXED POLYGON, MOUTH AND GLOBAL HMAX; LOCAL FE REFINEMENT\n');
fprintf('============================================================\n');
fprintf('  fixed phi=%+.9f deg; mouth=[%.10e %.10e]\n', ...
    rad2deg(local_wrap(I.phi_star)),I.x_star(1),I.x_star(2));
fprintf('  a0=%.6g; NArc=%d; mouth half-shift=%.7g\n', ...
    O.a0,O.NArc,Cbase.mesh2.chw);
fprintf('  global Hmax FIXED at %.7g m; nominal base Hmin=%.7g m\n', ...
    hmaxFixed,hBase);

for jf=1:nF
    C=Cbase;
    C.a0=O.a0;
    C.mesh1.hmin=hBase/factors(jf);
    C.mesh1.hhole=C.mesh1.hmin;
    C.mesh1.hmax=hmaxFixed;
    C.mesh2.hmax=hmaxFixed;
    C.mesh2.hhole=C.mesh1.hmin;
    C.mesh2.hcrack=C.mesh1.hmin;
    C.solver.verbose=0;
    fprintf('\n--- factor=%.3f | nominal Hmin=%.7e ---\n', ...
        factors(jf),C.mesh1.hmin);

    [~,D,~,Mc]=build_stage2_cracked_mesh_for_theta(C,I,0, ...
        'NArc',O.NArc,'PlotGeom',false,'PlotMesh',false, ...
        'PlotCollapsed',logical(O.Plot)&&jf==1);
    if norm(D.Pmid(1,:)-I.x_star)>1e-11 || ...
            abs(norm(D.Pmid(end,:)-D.Pmid(1,:))-C.a0)>1e-10
        error('step24:GeometryChanged','Fixed mouth or a0 unexpectedly moved.');
    end

    S2=solve_cracked_LEFM(C,Mc,'lambda',1.0);
    H=tip_edges(S2.mesh,Mc.crack.Pmid(end,:));
    hTip(jf)=H.median;
    nNodes(jf)=size(S2.mesh.coord,1);
    nElem(jf)=size(S2.mesh.connect,1);
    nTri(jf)=size(S2.mesh.connect3,1);
    nSig(jf)=numel(S2.U);
    fprintf('  Actual median tip edge=%.7e; nT6nodes=%d; nT3=%d\n', ...
        H.median,nNodes(jf),nTri(jf));

    mat=S2.mat;
    if ~isfield(mat,'Dmat'),mat.Dmat=mat.D;end
    for ir=1:nR
        ro=rat(ir)*C.a0;
        ri=max(0.10*ro,2*H.median);
        if ri>=ro
            error('step24:BadAnnulus', ...
                'At factor %.3f, r_i/r_o=%.3f >=1.',factors(jf),ri/ro);
        end
        [KI(jf,ir),KII(jf,ir)]=SIF_LEFM_interaction_EDI( ...
            S2.mesh,S2.U,Mc.crack.Pmid,mat, ...
            struct('r_inner',ri,'r_outer',ro), ...
            'UsePlaneStrain',mat.ps==1, ...
            'Verbose',false,'WeightFunction','fe_nodal');
        rIn(jf,ir)=ri;
        fprintf('  ro/a0=%.2f | ri/ro=%.3f | KI=%.8e | KII=%+.8e | ratio=%+.8e\n', ...
            rat(ir),ri/ro,KI(jf,ir),KII(jf,ir),KII(jf,ir)/KI(jf,ir));
    end
end

rows=nan(nF*nR,9);k=0;
for jf=1:nF
    for ir=1:nR
        k=k+1;
        rows(k,:)=[factors(jf),hBase/factors(jf), ...
            hTip(jf),nNodes(jf),rat(ir),KI(jf,ir),KII(jf,ir), ...
            KII(jf,ir)/KI(jf,ir),rIn(jf,ir)/(rat(ir)*O.a0)];
    end
end
T=array2table(rows,'VariableNames',{ ...
    'refinement_factor','Hmin_nominal','h_tip_actual','n_T6nodes', ...
    'r_outer_over_a0','KI','KII','KII_over_KI','r_inner_over_outer'});

S=nan(nF,8);
[~,iR]=min(abs(rat-0.65));
for jf=1:nF
    q=KII(jf,:)./KI(jf,:);
    S(jf,:)=[factors(jf),hBase/factors(jf),hTip(jf), ...
        KI(jf,iR),q(iR),min(q),max(q),max(q)-min(q)];
end
Summary=array2table(S,'VariableNames',{ ...
    'refinement_factor','Hmin_nominal','h_tip_actual', ...
    'KI_ref','KII_over_KI_ref','min_KII_over_KI', ...
    'max_KII_over_KI','EDI_domain_spread'});

fprintf('\nALL EDI RESULTS\n');disp(T);
fprintf('\nFIXED GEOMETRY LOCAL MESH SENSITIVITY\n');disp(Summary);

if logical(O.Plot)
    figure('Name','Step 24: fixed-polygon mesh convergence','Color','w');
    clf;hold on;box on;grid on;
    for ir=1:nR
        plot(hTip,KII(:,ir)./KI(:,ir),'-o','LineWidth',1.2, ...
            'DisplayName',sprintf('ro/a0=%.2f',rat(ir)));
    end
    yline(0,'k:','HandleVisibility','off');
    xlabel('actual median tip edge [m]');
    ylabel('K_{II}/K_I');
    title('Fixed hole polygon and fixed crack mouth, a0=8 mm');
    legend('Location','best');
end

O24=struct('Stage22',O22,'settings',O,'fixedMouth',I.x_star, ...
    'fixedPhiDeg',rad2deg(local_wrap(I.phi_star)), ...
    'factors',factors,'rOuterOverA0',rat,'KI',KI,'KII',KII, ...
    'rInner',rIn,'hTip',hTip,'nNodes',nNodes,'nElements',nElem, ...
    'nT3',nTri,'table',T,'summary',Summary);
fprintf('\nSTEP 24 completed.\n');
fprintf(['Gate: the signed ratio must stabilize as ACTUAL local tip edges ', ...
    'decrease with the same polygon, mouth, width and global Hmax. ', ...
    'If it stabilizes, test polygon/appendix-width dependence separately.\n']);
end

function H=tip_edges(mesh,tip)
X=mesh.coord3;T=mesh.connect3;
d=hypot(X(:,1)-tip(1),X(:,2)-tip(2));
tol=max(1e-12,1e-8*max(1,max(abs(X(:)))));
v=find(d<=min(d)+tol);
te=T(any(ismember(T,v),2),:);
L=[];
for j=1:size(te,1)
    P=X(te(j,:),:);
    L=[L,norm(P(2,:)-P(1,:)),norm(P(3,:)-P(2,:)), ...
       norm(P(1,:)-P(3,:))]; %#ok<AGROW>
end
L=L(isfinite(L)&L>tol);
if isempty(L),error('step24:NoTipEdges','No positive tip-adjacent edges.');end
H=struct('median',median(L),'min',min(L),'max',max(L));
end

function a=local_wrap(a)
a=mod(a+pi,2*pi)-pi;
end
