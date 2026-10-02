function Out=main_step21_asymmetric_initiation_perturbation(O20,varargin)
%MAIN_STEP21_ASYMMETRIC_INITIATION_PERTURBATION
% Diagnostic test: does signed Mode II respond when the crack mouth is
% deliberately displaced from the Stage-I stress maximum?
%
% Uses the existing O20 Stage-I FEM solution and the signed KII(theta)
% slope from Step 20. For each offset in hole polar angle:
%   - rebuild the exact hole mouth and material-normal/tangent frame;
%   - insert a normal (theta=0) short crack at that perturbed point;
%   - solve a NEW full-domain cracked FEM mesh;
%   - evaluate KI, signed KII at several FE-nodal interaction-EDI domains.
%
% The inferred theta_correction=-KII(0)/[dKII/dtheta]_Step20 is only
% a linear-response DIAGNOSTIC. It is not an independently solved
% Stage-II root at the perturbed initiation point.
%
% Example:
%   O21=main_step21_asymmetric_initiation_perturbation(O20);

ip=inputParser;
addParameter(ip,'PhiOffsetsDeg',[-1,-0.5,0,0.5,1], ...
    @(v)isnumeric(v)&&isvector(v)&&all(isfinite(v)));
addParameter(ip,'Plot',true,@(v)islogical(v)||isnumeric(v));
parse(ip,varargin{:});
O=ip.Results;

if nargin<1 || ~isstruct(O20) || ~isfield(O20,'Stage1') || ...
        ~isfield(O20,'KII') || ~isfield(O20,'thetaDeg')
    error('step21:NeedStep20','Pass the existing output O20 from Step 20.');
end

C=O20.C;
I0=O20.Stage1.I;
B=O20.Stage1.B;
phi0=atan2(sin(I0.phi_star),cos(I0.phi_star));
offsetDeg=O.PhiOffsetsDeg(:).';
rRat=O20.rOuterOverA0(:).';
thetaDeg=O20.thetaDeg(:);
thetaRad=deg2rad(thetaDeg);

[~,iPlus]=min(abs(thetaDeg-1));
[~,iMinus]=min(abs(thetaDeg+1));
[~,iZero]=min(abs(thetaDeg));
if any(abs([thetaDeg(iPlus)-1,thetaDeg(iMinus)+1,thetaDeg(iZero)])>1e-9)
    error('step21:NeedProbeAngles', ...
        'Step 20 must include theta=-1,0,+1 deg for local sensitivity.');
end
if ~any(abs(offsetDeg)<1e-10)
    error('step21:NeedBaseline','PhiOffsetsDeg must contain zero.');
end

nP=numel(offsetDeg);
nR=numel(rRat);
KI=nan(nP,nR);
KII=nan(nP,nR);
rIn=nan(nP,nR);
hTip=nan(nP,1);
phiPertDeg=nan(nP,1);
radialStress=nan(nP,1);

slope=zeros(1,nR);
for ir=1:nR
    slope(ir)=(O20.KII(iPlus,ir)-O20.KII(iMinus,ir))/ ...
        (thetaRad(iPlus)-thetaRad(iMinus));
    if ~(isfinite(slope(ir)) && abs(slope(ir))>eps)
        error('step21:BadSlope','Signed KII slope is invalid for domain %d.',ir);
    end
end

fprintf('\n============================================================\n');
fprintf('SIF AUDIT STEP 21: ASYMMETRIC INITIATION PERTURBATION\n');
fprintf('============================================================\n');
fprintf('  baseline phi = %+.8f deg; Npoly=%d; a0=%.6g m\n', ...
    rad2deg(phi0),C.hole.npoly,C.a0);
fprintf('  independent Step-20 slopes dKII/dtheta [per rad]: %s\n', ...
    mat2str(slope,8));

Cases=cell(nP,1);
for iphi=1:nP
    delta=deg2rad(offsetDeg(iphi));
    phi=phi0+delta;

    n=[cos(phi),sin(phi)];
    t=[-sin(phi),cos(phi)];

    I=I0;
    I.phi_star=mod(phi,2*pi);
    I.x_star=C.hole.center(:).'+C.hole.r*n;
    I.n_mat_star=n;
    I.n_hole_star=-n;
    I.t_hat_star=t;

    phiPertDeg(iphi)=rad2deg(atan2(sin(phi),cos(phi)));
    radialStress(iphi)=local_boundary_interp(B,phi);

    fprintf('\n--- delta_phi=%+6.3f deg | phi=%+10.6f deg ---\n', ...
        offsetDeg(iphi),phiPertDeg(iphi));

    plotCase=logical(O.Plot)&&abs(offsetDeg(iphi))<1e-10;
    [G2,D,M,Mc]=build_stage2_cracked_mesh_for_theta(C,I,0, ...
        'PlotGeom',false,'PlotMesh',false, ...
        'PlotCollapsed',plotCase);

    S2=solve_cracked_LEFM(C,Mc,'lambda',1.0);
    H=local_tip_mesh_scale(S2.mesh,Mc.crack.Pmid(end,:));
    hTip(iphi)=H.median;

    mat=S2.mat;
    if ~isfield(mat,'Dmat'),mat.Dmat=mat.D;end

    for ir=1:nR
        rOut=rRat(ir)*C.a0;
        rin=max(0.10*rOut,2*H.median);
        if rin>=rOut
            error('step21:BadDomain','EDI annulus invalid at delta_phi=%g.', ...
                offsetDeg(iphi));
        end
        [KI(iphi,ir),KII(iphi,ir)]=SIF_LEFM_interaction_EDI( ...
            S2.mesh,S2.U,Mc.crack.Pmid,mat, ...
            struct('r_inner',rin,'r_outer',rOut), ...
            'UsePlaneStrain',mat.ps==1,'Verbose',false, ...
            'WeightFunction','fe_nodal');
        rIn(iphi,ir)=rin;

        thetaProxy=-KII(iphi,ir)/slope(ir);
        fprintf([' r/a0=%.2f | KI=%.8e | KII=%+.8e | ', ...
            'KII/KI=%+.5e | proxy correction=%+.5f deg\n'], ...
            rRat(ir),KI(iphi,ir),KII(iphi,ir), ...
            KII(iphi,ir)/KI(iphi,ir),rad2deg(thetaProxy));
    end

    Cases{iphi}=struct('I',I,'G2',G2,'D',D,'M',M, ...
        'Mc',Mc,'S2',S2,'H',H);
end

thetaProxyDeg=-rad2deg(KII./repmat(slope,nP,1));
globalProxyDeg=phiPertDeg+thetaProxyDeg;

rows=nan(nP*nR,9);
k=0;
for iphi=1:nP
    for ir=1:nR
        k=k+1;
        rows(k,:)=[ ...
            offsetDeg(iphi),phiPertDeg(iphi),rRat(ir), ...
            KI(iphi,ir),KII(iphi,ir),KII(iphi,ir)/KI(iphi,ir), ...
            thetaProxyDeg(iphi,ir),globalProxyDeg(iphi,ir), ...
            radialStress(iphi)/I0.sig_tt_pos_unit];
    end
end

T=array2table(rows,'VariableNames',{ ...
    'phi_offset_deg','mouth_phi_deg','r_outer_over_a0', ...
    'KI','KII','KII_over_KI','theta_proxy_deg', ...
    'global_heading_proxy_deg','boundary_stress_over_peak'});

[~,iBase]=min(abs(offsetDeg));
baselineRows=nan(nR,5);
for ir=1:nR
    qBase=KII(iBase,ir);
    qPrior=O20.KII(iZero,ir);
    baselineRows(ir,:)=[rRat(ir),qBase,qPrior, ...
        abs(qBase-qPrior)/max(abs(KI(iBase,ir)),eps), ...
        KII(iBase,ir)/KI(iBase,ir)];
end
Baseline=array2table(baselineRows,'VariableNames',{ ...
    'r_outer_over_a0','KII_rerun','KII_previous', ...
    'KII_difference_over_KI','KII_rerun_over_KI'});

fprintf('\nMOUTH-PERTURBATION SENSITIVITY\n');
disp(T);
fprintf('\nBASELINE REPRODUCIBILITY CHECK (phi offset = 0)\n');
disp(Baseline);

if logical(O.Plot)
    figure('Name','Step 21: KII versus mouth perturbation','Color','w');clf
    hold on;grid on;box on
    for ir=1:nR
        plot(offsetDeg,KII(:,ir)./KI(:,ir),'-o','LineWidth',1.1, ...
            'DisplayName',sprintf('r_o/a_0=%.2f',rRat(ir)));
    end
    xline(0,'k--','HandleVisibility','off');
    yline(0,'k:','HandleVisibility','off');
    xlabel('mouth offset from Stage-I peak [deg]');
    ylabel('K_{II}/K_I at normal appendix');
    title('Does signed Mode II detect displaced initiation?');
    legend('Location','best');

    figure('Name','Step 21: linearized direction proxy','Color','w');clf
    hold on;grid on;box on
    for ir=1:nR
        plot(offsetDeg,thetaProxyDeg(:,ir),'-o','LineWidth',1.1, ...
            'DisplayName',sprintf('r_o/a_0=%.2f',rRat(ir)));
    end
    xlabel('mouth offset [deg]');
    ylabel('linear-response angle proxy [deg]');
    title('Diagnostic only: -KII(0)/(dKII/dtheta)');
    legend('Location','best');
end

Out=struct();
Out.settings=O;
Out.Stage20=O20;
Out.phiOffsetsDeg=offsetDeg;
Out.phiPertDeg=phiPertDeg;
Out.KI=KI;
Out.KII=KII;
Out.rInner=rIn;
Out.hTip=hTip;
Out.slopeFromStep20=slope;
Out.thetaProxyDeg=thetaProxyDeg;
Out.globalHeadingProxyDeg=globalProxyDeg;
Out.table=T;
Out.baseline=Baseline;
Out.Cases=Cases;

fprintf('\nSTEP 21 completed.\n');
fprintf(['Gate: off-peak mouth perturbations should generate resolved signed ', ...
    'KII while the zero-offset baseline reproduces Step 20. ', ...
    'Linear-response angles here are diagnostics, not solved roots.\n']);
end

function sig=local_boundary_interp(B,phi)
ph=B.phi(:);
y=B.sig_tt_eff(:);
q=mod(phi,2*pi);
sig=interp1([ph-2*pi;ph;ph+2*pi],[y;y;y],q,'linear');
end

function H=local_tip_mesh_scale(mesh,tip)
X=mesh.coord3;
T=mesh.connect3;
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
if isempty(L)
    error('step21:NoTipEdges','Cannot estimate the tip mesh scale.');
end
H=struct('median',median(L),'min',min(L),'max',max(L), ...
    'nTipElements',size(Te,1));
end
