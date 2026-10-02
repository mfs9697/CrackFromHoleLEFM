function Out=main_step23_asymmetric_fixed_mouth_convergence(O22,varargin)
%MAIN_STEP23_ASYMMETRIC_FIXED_MOUTH_CONVERGENCE
% Disentangle Stage-II retained-hole polygon resolution and FE mesh size.
%
% Reuse Step-22 fine Stage-I mouth (no new initiation optimization).
% At a0=4,8 mm, compare:
%   fine h (nominal Npoly=480) at NArc=160,320,480;
%   coarse h (nominal Npoly=240) at NArc=480 only.
%
% Existing O22 NArc=160 fine-mesh solutions and EDI outputs are reused
% (not recomputed). Polygon refinement also changes generated mesh
% topology; the NArc sweep is therefore a combined geometry/remeshing
% convergence test, not a mathematically pure geometric perturbation.
% A fixed NArc=480 coarse/fine comparison isolates nominal FE refinement,
% although enforced polygon edges limit the actual local mesh scale.
%
% This is a fixed-mouth NORMAL-crack test (theta=0). Small nonzero KII is
% NOT automatically interpreted as a physical crack-direction correction.
%
% Example:
%   O23=main_step23_asymmetric_fixed_mouth_convergence(O22);

ip=inputParser;
addParameter(ip,'Lengths',[0.004 0.008], ...
    @(v)isnumeric(v)&&isvector(v)&&all(isfinite(v))&&all(v>0));
addParameter(ip,'NArcList',[160 320 480], ...
    @(v)isnumeric(v)&&isvector(v)&&all(v>=32)&&all(v==round(v)));
addParameter(ip,'CoarseMeshNpoly',240, ...
    @(v)isnumeric(v)&&isscalar(v)&&v>=120&&v==round(v));
addParameter(ip,'Plot',true,@(v)islogical(v)||isnumeric(v));
parse(ip,varargin{:});
O=ip.Results;

if nargin<1 || ~isstruct(O22) || ~isfield(O22,'Stage1') || ...
        ~isfield(O22,'Cases') || ~isfield(O22,'KII') || ...
        ~isfield(O22,'rOuterOverA0')
    error('step23:NeedStep22', ...
        'Pass O22 returned by main_step22_asymmetric_length_sensitivity.');
end
if ~ismember(160,O.NArcList)
    error('step23:NeedBaseline', ...
        'NArcList must contain 160 to reuse the existing Step-22 baseline.');
end

Cfixed=O22.Stage1.C;
I=O22.Stage1.I;
aList=O.Lengths(:).';
arcList=sort(unique(O.NArcList(:).'));
rRat=O22.rOuterOverA0(:).';
nA=numel(aList);nG=numel(arcList);nR=numel(rRat);
fineH=Cfixed.mesh1.hmin;
coarseH=2*pi*Cfixed.hole.r/O.CoarseMeshNpoly;
if ~(coarseH>fineH)
    error('step23:BadMeshOrdering', ...
        'CoarseMeshNpoly must be less than the fine Stage-I Npoly.');
end

% NArc=160 fine-mesh values come from O22. Verify they correspond to
% exactly the same fixed Stage-I point and requested length.
for ia=1:nA
    old=find(abs(O22.lengthList-aList(ia))<=1e-11,1);
    if isempty(old)
        error('step23:MissingCachedLength', ...
            'Length %.6g is not in O22.lengthList.',aList(ia));
    end
end

KI=nan(nA,nG,2,nR);
KII=nan(nA,nG,2,nR);
rInner=nan(nA,nG,2,nR);
tipH=nan(nA,nG,2);
source=cell(nA,nG,2);
% Third index: 1=coarse mesh (only max NArc), 2=fine mesh.
Gmax=max(arcList);

fprintf('\n============================================================\n');
fprintf('SIF AUDIT STEP 23: FIXED-MOUTH GEOMETRY AND MESH CONVERGENCE\n');
fprintf('============================================================\n');
fprintf('  fixed phi_*= %+.10f deg\n',rad2deg(local_wrap(I.phi_star)));
fprintf('  fixed mouth = [%.10e %.10e] m\n',I.x_star(1),I.x_star(2));
fprintf('  fine nominal h=%.7e, coarse nominal h=%.7e m\n',fineH,coarseH);
fprintf('  NArc=[%s] (points in retained hole arc)\n',num2str(arcList));

for ia=1:nA
    a0=aList(ia);
    idxCache=find(abs(O22.lengthList-a0)<=1e-11,1);
    for ig=1:nG
        nArc=arcList(ig);
        for imesh=1:2
            if imesh==1 && nArc~=Gmax
                continue
            end

            if imesh==2 && nArc==160
                KI(ia,ig,imesh,:)=reshape(O22.KI(idxCache,:),[1 1 1 nR]);
                KII(ia,ig,imesh,:)=reshape(O22.KII(idxCache,:),[1 1 1 nR]);
                rInner(ia,ig,imesh,:)=reshape(O22.rInner(idxCache,:),[1 1 1 nR]);
                tipH(ia,ig,imesh)=O22.hTip(idxCache);
                source{ia,ig,imesh}='O22 cached';
                fprintf('\n  a0=%.4f m | NArc=%d | fine h: reused Step 22\n', ...
                    a0,nArc);
                continue
            end

            C=Cfixed;
            C.a0=a0;
            if imesh==1
                h=coarseH;
                meshLabel='coarse';
            else
                h=fineH;
                meshLabel='fine';
            end
            C.mesh1.hmin=h;
            C.mesh1.hhole=h;
            C.mesh1.hmax=20*h;
            C.mesh2.hmax=C.mesh1.hmax;
            C.mesh2.hhole=h;
            C.mesh2.hcrack=h;
            C.solver.verbose=0;

            fprintf('\n--- a0=%.4f m | NArc=%d | mesh=%s ---\n', ...
                a0,nArc,meshLabel);

            [~,D,~,Mc]=build_stage2_cracked_mesh_for_theta( ...
                C,I,0,'NArc',nArc, ...
                'PlotGeom',false,'PlotMesh',false, ...
                'PlotCollapsed',false);
            % Ensure the intended control variables really did not move.
            if norm(D.Pmid(1,:)-I.x_star)>1e-11 || ...
                    abs(norm(D.Pmid(end,:)-D.Pmid(1,:))-a0)>1e-10
                error('step23:ChangedMouth', ...
                    'The fixed mouth or nominal crack length changed.');
            end

            S2=solve_cracked_LEFM(C,Mc,'lambda',1.0);
            H=local_tip_mesh_scale(S2.mesh,Mc.crack.Pmid(end,:));
            tipH(ia,ig,imesh)=H.median;

            mat=S2.mat;
            if ~isfield(mat,'Dmat'),mat.Dmat=mat.D;end

            for ir=1:nR
                rOut=rRat(ir)*a0;
                rin=max(0.1*rOut,2*H.median);
                if rin>=rOut
                    error('step23:BadAnnulus', ...
                        'Invalid EDI annulus: a0=%.4g,NArc=%d,%s,r/a0=%.2f.', ...
                        a0,nArc,meshLabel,rRat(ir));
                end
                [k1,k2]=SIF_LEFM_interaction_EDI( ...
                    S2.mesh,S2.U,Mc.crack.Pmid,mat, ...
                    struct('r_inner',rin,'r_outer',rOut), ...
                    'UsePlaneStrain',mat.ps==1, ...
                    'Verbose',false,'WeightFunction','fe_nodal');
                KI(ia,ig,imesh,ir)=k1;
                KII(ia,ig,imesh,ir)=k2;
                rInner(ia,ig,imesh,ir)=rin;
                fprintf(['  r_o/a0=%.2f | KI=%.8e | KII=%+.8e | ', ...
                    'KII/KI=%+.6e | ri/ro=%.3f\n'], ...
                    rRat(ir),k1,k2,k2/k1,rin/rOut);
            end
            source{ia,ig,imesh}='fresh FEM';
        end
    end
end

% Comprehensive per-EDI-domain table (all available configurations).
rows=nan(0,11);
for ia=1:nA
    for ig=1:nG
        for imesh=1:2
            if ~isfinite(KI(ia,ig,imesh,1)),continue;end
            for ir=1:nR
                k1=KI(ia,ig,imesh,ir);
                k2=KII(ia,ig,imesh,ir);
                h=(imesh==1)*coarseH+(imesh==2)*fineH;
                rows(end+1,:)=[aList(ia),arcList(ig), ...
                    imesh,h,tipH(ia,ig,imesh),rRat(ir), ...
                    k1,k2,k2/k1,rInner(ia,ig,imesh,ir)/(rRat(ir)*aList(ia)), ...
                    double(strcmp(source{ia,ig,imesh},'O22 cached'))]; %#ok<AGROW>
            end
        end
    end
end
T=array2table(rows,'VariableNames',{ ...
    'a0','NArc','mesh_code_1coarse_2fine','h_nominal','h_tip', ...
    'r_outer_over_a0','KI','KII','KII_over_KI', ...
    'r_inner_over_outer','reused_step22'});

% Compare geometry and mesh effects on the reference EDI domain.
[~,irRef]=min(abs(rRat-0.65));
summary=nan(nA,10);
for ia=1:nA
    qFine=squeeze(KII(ia,:,2,irRef)./KI(ia,:,2,irRef));
    kFine=squeeze(KI(ia,:,2,irRef));
    [~,igFine]=max(arcList);
    [~,igBase]=min(abs(arcList-160));
    qCoarse=KII(ia,igFine,1,irRef)/KI(ia,igFine,1,irRef);

    summary(ia,:)=[ ...
        aList(ia), ...
        qFine(igBase),qFine(igFine), ...
        qFine(igFine)-qFine(igBase), ...
        max(qFine)-min(qFine), ...
        qCoarse,qFine(igFine)-qCoarse, ...
        kFine(igBase),kFine(igFine), ...
        tipH(ia,igFine,2)];
end
S=array2table(summary,'VariableNames',{ ...
    'a0','ratio_arc160_fine_h','ratio_max_arc_fine_h', ...
    'ratio_polygon_endpoint_change','ratio_polygon_range', ...
    'ratio_max_arc_coarse_h','ratio_FE_change_fixed_max_arc', ...
    'KI_arc160_fine_h','KI_max_arc_fine_h','h_tip_max_arc_fine'});

fprintf('\nFULL GEOMETRY/MESH/EDI RESULTS\n');disp(T);
fprintf('\nREFERENCE EDI DOMAIN SENSITIVITY (r_o/a0=%.2f)\n',rRat(irRef));
disp(S);
fprintf(['NOTE: polygon refinement regenerates the FEM mesh; ',
    'the polygon series measures combined geometry/mesh changes. ',
    'The fixed-NArc coarse/fine comparison isolates nominal FE size.\n']);

if logical(O.Plot)
    figure('Name','Step 23: fixed-mouth polygon convergence','Color','w');
    clf;hold on;box on;grid on;
    for ia=1:nA
        q=squeeze(KII(ia,:,2,irRef)./KI(ia,:,2,irRef));
        plot(arcList,q,'-o','LineWidth',1.1, ...
            'DisplayName',sprintf('a_0=%.0f mm',1000*aList(ia)));
    end
    yline(0,'k:','HandleVisibility','off');
    xlabel('NArc (retained hole boundary points)');
    ylabel('K_{II}/K_I at fixed mouth, fine mesh');
    title('Geometry/remeshing sensitivity at theta=0');
    legend('Location','best');
end

Out=struct();
Out.settings=O;
Out.fixedPhiDeg=rad2deg(local_wrap(I.phi_star));
Out.fixedMouth=I.x_star;
Out.arcList=arcList;
Out.lengthList=aList;
Out.rOuterOverA0=rRat;
Out.KI=KI;
Out.KII=KII;
Out.rInner=rInner;
Out.hTip=tipH;
Out.source=source;
Out.table=T;
Out.summary=S;
Out.referenceDomainRatio=rRat(irRef);
fprintf('\nSTEP 23 completed.\n');
fprintf(['Gate: do not interpret the finite-length residual as physical ', ...
    'unless it is stable under retained-hole polygon refinement AND ', ...
    'under fixed-geometry FE refinement, both at the identical mouth.\n']);
end

function a=local_wrap(a)
a=mod(a+pi,2*pi)-pi;
end

function H=local_tip_mesh_scale(mesh,tip)
X=mesh.coord3;T=mesh.connect3;
d=sqrt(sum((X-tip).^2,2));
tol=max(1e-12,1e-8*max(1,max(abs(X(:)))));
tipNodes=find(d<=min(d)+tol);
Te=T(any(ismember(T,tipNodes),2),:);
L=[];
for k=1:size(Te,1)
    P=X(Te(k,:),:);
    L=[L,norm(P(2,:)-P(1,:)),norm(P(3,:)-P(2,:)), ...
        norm(P(1,:)-P(3,:))]; %#ok<AGROW>
end
L=L(isfinite(L)&L>tol);
if isempty(L),error('step23:NoTipEdges','Cannot determine tip mesh scale.');end
H=struct('median',median(L),'min',min(L),'max',max(L));
end
