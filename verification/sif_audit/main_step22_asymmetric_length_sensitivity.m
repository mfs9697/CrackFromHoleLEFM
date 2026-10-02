function Out=main_step22_asymmetric_length_sensitivity(O20,varargin)
%MAIN_STEP22_ASYMMETRIC_LENGTH_SENSITIVITY
% Refine the asymmetric full-domain Stage-II result and investigate the
% effect of the initial short-crack length.
%
% Scientific question:
%   At the Stage-I maximum of boundary tangential stress, does a normal
%   short crack retain KII approximately zero as a0 is varied?
%
% Scope:
%   - one independent Npoly=480 Stage-I solve with c=3 peak fitting;
%   - one NEW cracked FEM solve per short-crack length;
%   - signed FE-nodal interaction EDI on three domains per solved mesh;
%   - Npoly=240,a0=0.004 Step-20 result as reference only.
% Stage-II's existing full-domain geometry builder uses nArc=160 for its
% retained-hole polygon, independent of Stage-I Npoly. This experiment
% therefore refines the FE mesh, not the Stage-II hole polygon. Keep the
% polygon fixed across a0 to isolate the crack-length response.
%
% IMPORTANT: this is a length-sensitivity experiment at theta=0, not a
% local-symmetry root search. For any resolved nonzero KII, a separate
% theta sweep is required to establish the true correction angle.
%
% Usage:
%   O22=main_step22_asymmetric_length_sensitivity(O20);
%   O22=main_step22_asymmetric_length_sensitivity(O20,'Step21',O21);

ip=inputParser;
addParameter(ip,'Npoly',480,@(x)isnumeric(x)&&isscalar(x)&&isfinite(x)&&x>=120);
addParameter(ip,'LengthList',[0.002 0.004 0.006 0.008], ...
    @(x)isnumeric(x)&&isvector(x)&&all(isfinite(x))&&all(x>0));
addParameter(ip,'ROuterOverA0',[0.50 0.65 0.80], ...
    @(x)isnumeric(x)&&isvector(x)&&all(isfinite(x))&&all(x>0)&&all(x<1));
addParameter(ip,'Step21',[],@(x)isempty(x)||isstruct(x));
addParameter(ip,'Plot',true,@(x)islogical(x)||isnumeric(x));
parse(ip,varargin{:});
O=ip.Results;

if nargin<1 || ~isstruct(O20) || ~isfield(O20,'Stage1') || ...
        ~isfield(O20,'KII') || ~isfield(O20,'KI')
    error('step22:NeedStep20', ...
        'Pass the existing Step-20 output O20 for the coarse reference.');
end
addpath(genpath(pwd));

C=cfg_asymmetric_full_domain();
C.hole.npoly=round(O.Npoly);
C.holes={C.hole};
hArc=2*pi*C.hole.r/C.hole.npoly;
C.mesh1.hmin=hArc;
C.mesh1.hhole=hArc;
C.mesh1.hmax=20*hArc;
C.mesh2.hmax=C.mesh1.hmax;
C.mesh2.hhole=C.mesh1.hmin;
C.mesh2.hcrack=C.mesh1.hmin;
C.stage1.angular_fit_halfwidth_factor=3.0;
C.solver.verbose=0;
C.plot.show_mesh1=false;

R1=run_stage1_hole_initiation(C);
I=R1.I;
phiFine=local_wrap(I.phi_star);
phiCoarse=local_wrap(O20.Stage1.I.phi_star);
dPhiDeg=rad2deg(local_wrap(phiFine-phiCoarse));

aList=O.LengthList(:).';
ratios=O.ROuterOverA0(:).';
nA=numel(aList);
nR=numel(ratios);
KI=nan(nA,nR);
KII=nan(nA,nR);
rInner=nan(nA,nR);
hTip=nan(nA,1);
Cases=cell(nA,1);

fprintf('\n============================================================\n');
fprintf('SIF AUDIT STEP 22: ASYMMETRIC LENGTH SENSITIVITY\n');
fprintf('============================================================\n');
fprintf('  Npoly=%d | h_hole/R=%.6f deg\n',C.hole.npoly,rad2deg(hArc/C.hole.r));
fprintf('  NOTE: Stage-II retained-hole arc = 160 points at every a0 (fixed geometry).\n');
fprintf('  Stage-I phi refined = %+.8f deg\n',rad2deg(phiFine));
fprintf('  Step-20 coarse phi  = %+.8f deg\n',rad2deg(phiCoarse));
fprintf('  angular shift       = %+.8f deg\n',dPhiDeg);
fprintf('  refined max stress  = %.10e\n',I.sig_tt_pos_unit);
fprintf('  refined lambda_ini  = %.10e\n',I.lambda_ini);

if ~isempty(O.Step21) && isfield(O.Step21,'phiOffsetsDeg')
    [~,ia]=min(abs(O.Step21.phiOffsetsDeg+1.5));
    [~,ib]=min(abs(O.Step21.phiOffsetsDeg-1.5));
    if abs(O.Step21.phiOffsetsDeg(ia)+1.5)<1e-8 && ...
       abs(O.Step21.phiOffsetsDeg(ib)-1.5)<1e-8
        [~,oldR]=min(abs(O20.rOuterOverA0-0.65));
        dQdPhi=(O.Step21.KII(ib,oldR)-O.Step21.KII(ia,oldR))/3.0;
        fprintf(['  Step-21 coarse dKII/dphi ~ %.8e per degree; ', ...
                 'angular-shift scale |dKII| ~ %.8e\n'], ...
            dQdPhi,abs(dQdPhi*dPhiDeg));
    end
end

for ia=1:nA
    Ck=C;
    Ck.a0=aList(ia);
    if Ck.a0<=8*Ck.mesh2.chw
        error('step22:ShortPencil', ...
            'a0=%.6g is too close to the fixed mouth-shift scale.',Ck.a0);
    end

    fprintf('\n--- a0=%.6g m | a0/R=%.4f ---\n', ...
        Ck.a0,Ck.a0/Ck.hole.r);

    % The fitted hole point and normal are kept identical across a0.
    [G2,D,M,Mc]=build_stage2_cracked_mesh_for_theta( ...
        Ck,I,0, ...
        'PlotGeom',false,'PlotMesh',false, ...
        'PlotCollapsed',logical(O.Plot) && ia==1);

    S2=solve_cracked_LEFM(Ck,Mc,'lambda',1.0);
    H=local_tip_mesh_scale(S2.mesh,Mc.crack.Pmid(end,:));
    hTip(ia)=H.median;

    mat=S2.mat;
    if ~isfield(mat,'Dmat'),mat.Dmat=mat.D;end

    for ir=1:nR
        rOuter=ratios(ir)*Ck.a0;
        rin=max(0.1*rOuter,2*H.median);
        if rin>=rOuter
            error('step22:BadDomain', ...
                'EDI domain invalid: a0=%.6g,r/a0=%.3f,rIn/rOut=%.3f', ...
                Ck.a0,ratios(ir),rin/rOuter);
        end
        if rin/rOuter>0.7
            warning('step22:ThinAnnulus', ...
                'Narrow EDI annulus at a0=%.6g,r/a0=%.3f,rIn/rOut=%.3f.', ...
                Ck.a0,ratios(ir),rin/rOuter);
        end

        [KI(ia,ir),KII(ia,ir)]=SIF_LEFM_interaction_EDI( ...
            S2.mesh,S2.U,Mc.crack.Pmid,mat, ...
            struct('r_inner',rin,'r_outer',rOuter), ...
            'UsePlaneStrain',mat.ps==1,'Verbose',false, ...
            'WeightFunction','fe_nodal');

        rInner(ia,ir)=rin;

        fprintf(['  r_o/a0=%.2f | r_i/r_o=%.3f | KI=%.8e | ', ...
            'KII=%+.8e | KII/KI=%+.6e\n'], ...
            ratios(ir),rin/rOuter,KI(ia,ir), ...
            KII(ia,ir),KII(ia,ir)/KI(ia,ir));
    end

    Cases{ia}=struct('G2',G2,'D',D,'M',M,'Mc',Mc,'S2',S2,'H',H);
end

rows=nan(nA*nR,8);
k=0;
for ia=1:nA
    for ir=1:nR
        k=k+1;
        rows(k,:)=[ ...
            aList(ia),aList(ia)/C.hole.r,ratios(ir), ...
            KI(ia,ir),KII(ia,ir),KII(ia,ir)/KI(ia,ir), ...
            hTip(ia)/aList(ia),rInner(ia,ir)/(ratios(ir)*aList(ia))];
    end
end
T=array2table(rows,'VariableNames',{ ...
    'a0','a0_over_R','r_outer_over_a0','KI','KII','KII_over_KI', ...
    'h_tip_over_a0','r_inner_over_r_outer'});

summaryRows=nan(nA,7);
for ia=1:nA
    [~,iRef]=min(abs(ratios-0.65));
    q=KII(ia,:)./KI(ia,:);
    summaryRows(ia,:)=[aList(ia),aList(ia)/C.hole.r, ...
        KI(ia,iRef),q(iRef),max(abs(q)), ...
        max(q)-min(q),max(KI(ia,:))-min(KI(ia,:))];
end
Summary=array2table(summaryRows,'VariableNames',{ ...
    'a0','a0_over_R','KI_reference','KII_over_KI_reference', ...
    'max_abs_KII_over_KI','KII_over_KI_domain_spread','KI_domain_spread'});

fprintf('\nSTAGE-II LENGTH-SENSITIVITY TABLE\n');
disp(T);
fprintf('\nLENGTH AND EDI-DOMAIN SUMMARY\n');
disp(Summary);

% Matched-length coarse comparison, using nearest matching EDI radii.
coarse=nan(nR,5);
[~,i0]=min(abs(O20.thetaDeg));
[~,iFine]=min(abs(aList-0.004));
if abs(aList(iFine)-0.004)<1e-12
    for ir=1:nR
        [d,j]=min(abs(O20.rOuterOverA0-ratios(ir)));
        if d>1e-9,continue;end
        coarse(ir,:)=[ratios(ir), ...
            O20.KII(i0,j)/O20.KI(i0,j), ...
            KII(iFine,ir)/KI(iFine,ir), ...
            rad2deg(phiCoarse),rad2deg(phiFine)];
    end
end
CoarseComparison=array2table(coarse,'VariableNames',{ ...
    'r_outer_over_a0','ratio_N240','ratio_refined', ...
    'phi_N240_deg','phi_refined_deg'});
fprintf('\nMATCHED LENGTH a0=0.004: MESH AND INITIATION-ANGLE CHANGE\n');
disp(CoarseComparison);
fprintf(['Note: this comparison changes BOTH Npoly and the fitted Stage-I ', ...
    'angle. It is not a fixed-mouth mesh convergence test.\n']);

if logical(O.Plot)
    figure('Name','Step 22: KII versus short-crack length','Color','w');
    clf;hold on;box on;grid on
    for ir=1:nR
        plot(aList/C.hole.r,KII(:,ir)./KI(:,ir),'-o','LineWidth',1.1, ...
            'DisplayName',sprintf('r_o/a_0=%.2f',ratios(ir)));
    end
    yline(0,'k--','HandleVisibility','off');
    xlabel('a_0/R'); ylabel('K_{II}/K_I at theta=0');
    title('Asymmetric crack-length sensitivity at Stage-I peak');
    legend('Location','best');
end

Out=struct();
Out.settings=O;
Out.Stage1=R1;
Out.CoarseStage20=O20;
Out.lengthList=aList;
Out.rOuterOverA0=ratios;
Out.KI=KI;
Out.KII=KII;
Out.rInner=rInner;
Out.hTip=hTip;
Out.Cases=Cases;
Out.table=T;
Out.summary=Summary;
Out.coarseComparison=CoarseComparison;
Out.phiFineDeg=rad2deg(phiFine);
Out.phiCoarseDeg=rad2deg(phiCoarse);

fprintf('\nSTEP 22 completed.\n');
fprintf(['Gate: test whether KII/KI remains below the combined EDI-domain ', ...
    'and Stage-I angular uncertainty or becomes resolved and grows ', ...
    'systematically with a0. A nonzero ratio alone does not establish ', ...
    'a solved nonzero propagation direction.\n']);
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
if isempty(L),error('step22:NoTipEdges','Cannot determine tip mesh scale.');end
H=struct('median',median(L),'min',min(L),'max',max(L), ...
    'nTipElements',size(Te,1));
end
