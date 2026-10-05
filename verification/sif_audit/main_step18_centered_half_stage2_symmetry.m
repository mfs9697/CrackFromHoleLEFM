function Out=main_step18_centered_half_stage2_symmetry(varargin)
%MAIN_STEP18_CENTERED_HALF_STAGE2_SYMMETRY
% Stage-II symmetry verification on the centered right-half benchmark.
%
% The crack mouth is placed at the exact rightmost hole point.  For each
% trial appendix angle a fresh half-domain cracked mesh is generated and
% solved. Signed KI,KII are extracted with the FE-nodal interaction EDI.
%
% Expected symmetry limit:
%   KII(theta=0) -> 0,
%   KII(-theta) -> -KII(+theta),
%   KI(-theta)  ->  KI(+theta),
%   theta_root  -> 0.

ip=inputParser;
addParameter(ip,'NpolyList',[240 480],@(x)isnumeric(x)&&isvector(x)&&all(x>=32));
addParameter(ip,'ThetaDeg',[-2 -1 -0.5 0 0.5 1 2],@(x)isnumeric(x)&&isvector(x));
addParameter(ip,'ROuterOverA0',[0.50 0.65 0.80],@(x)isnumeric(x)&&isvector(x)&&all(x>0));
addParameter(ip,'Plot',true,@(x)islogical(x)||isnumeric(x));
parse(ip,varargin{:});
O=ip.Results;

addpath(genpath(pwd));

fprintf('\n============================================================\n');
fprintf('SIF AUDIT STEP 18: CENTERED HALF-DOMAIN STAGE-II SYMMETRY\n');
fprintf('============================================================\n');

Nlist=round(O.NpolyList(:));
thetaDeg=O.ThetaDeg(:).';
thetaRad=deg2rad(thetaDeg);
rRat=O.ROuterOverA0(:).';

nL=numel(Nlist);
nT=numel(thetaDeg);
nR=numel(rRat);

KI=nan(nL,nT,nR);
KII=nan(nL,nT,nR);
rin=nan(nL,nT,nR);
htip=nan(nL,nT);
uxsym=nan(nL,nT);
Stage1=cell(nL,1);
Cases=cell(nL,nT);

for il=1:nL
    C=cfg_centered_half_domain();
    C.hole.npoly=Nlist(il);
    C.holes={C.hole};
    hArc=2*pi*C.hole.r/C.hole.npoly;
    C.mesh1.hmin=hArc;
    C.mesh1.hhole=hArc;
    C.mesh1.hmax=20*hArc;
    C.mesh2.hmax=C.mesh1.hmax;
    C.mesh2.hhole=C.mesh1.hmin;
    C.mesh2.hcrack=C.mesh1.hmin;
    C.solver.verbose=0;
    C.plot.show_mesh1=false;

    R1=run_stage1_centered_half_domain(C);
    Stage1{il}=R1;

    fprintf('\n--- Npoly=%d | h/R=%.8g rad = %.6f deg ---\n', ...
        Nlist(il),hArc/C.hole.r,rad2deg(hArc/C.hole.r));
    fprintf('  Stage-I fitted phi = %+.8f deg; exact Stage-II mouth uses 0 deg\n', ...
        rad2deg(R1.I.phi_star));

    for it=1:nT
        th=thetaRad(it);
        fprintf('  theta=%+7.3f deg ... ',thetaDeg(it));

        doPlotCase=logical(O.Plot) && il==1 && abs(thetaDeg(it))<1e-12;

        [G2,D,M,Mc]=build_stage2_centered_half_cracked_mesh_for_theta( ...
            C,R1.I,th, ...
            'PlotGeom',false,'PlotMesh',false,'PlotCollapsed',doPlotCase);

        S2=solve_cracked_LEFM(C,Mc,'lambda',1.0);

        H=local_tip_mesh_scale(S2.mesh,Mc.crack.Pmid(end,:));
        htip(il,it)=H.median;

        if isfield(S2.bc,'symmetry_nodes')
            sn=S2.bc.symmetry_nodes;
            uxsym(il,it)=max(abs(S2.U(2*sn-1)));
        end

        mat=S2.mat;
        if ~isfield(mat,'Dmat')
            mat.Dmat=mat.D;
        end

        V=Mc.crack.Pmid;
        for ir=1:nR
            rOut=rRat(ir)*C.a0;
            rIn=max(0.10*rOut,2.0*H.median);
            if ~(rIn>=0 && rIn<rOut)
                error('step18:BadEDIDomain', ...
                    'r_inner >= r_outer at Npoly=%d theta=%g r/a0=%g.', ...
                    Nlist(il),thetaDeg(it),rRat(ir));
            end

            dom=struct('r_inner',rIn,'r_outer',rOut);
            [KI(il,it,ir),KII(il,it,ir)]=SIF_LEFM_interaction_EDI( ...
                S2.mesh,S2.U,V,mat,dom, ...
                'UsePlaneStrain',mat.ps==1, ...
                'Verbose',false, ...
                'WeightFunction','fe_nodal');

            rin(il,it,ir)=rIn;
        end

        Cases{il,it}=struct('C',C,'G2',G2,'D',D,'M',M,'Mc',Mc,'S2',S2,'H',H);

        [~,iRef]=min(abs(rRat-0.65));
        fprintf('KI=%.6e KII=%+.6e KII/KI=%+.3e | htip=%.3e\n', ...
            KI(il,it,iRef),KII(il,it,iRef), ...
            KII(il,it,iRef)/KI(il,it,iRef),H.median);
    end
end

% -------------------------------------------------------------------------
% Symmetry/root diagnostics for each mesh and EDI outer radius
% -------------------------------------------------------------------------
rows=nan(nL*nR,16);
rr=0;

for il=1:nL
    for ir=1:nR
        rr=rr+1;

        ki=squeeze(KI(il,:,ir));
        kii=squeeze(KII(il,:,ir));

        [~,i0]=min(abs(thetaDeg));
        KI0=ki(i0);
        KII0=kii(i0);
        ratio0=KII0/KI0;

        p=polyfit(thetaRad,kii,1);
        thetaRoot=-p(2)/p(1);

        thetaBracketRoot=local_sign_change_root(thetaRad,kii);

        oddDef=0;
        evenDef=0;
        for it=1:nT
            th=thetaDeg(it);
            if th<=0, continue; end
            [dm,j]=min(abs(thetaDeg+th));
            if dm>1e-12, continue; end
            oddDef=max(oddDef,abs(kii(it)+kii(j)));
            evenDef=max(evenDef,abs(ki(it)-ki(j)));
        end
        oddDef=oddDef/max(abs(KI0),eps);
        evenDef=evenDef/max(abs(KI0),eps);

        yfit=polyval(p,thetaRad);
        fitRel=norm(kii-yfit)/max(norm(kii),eps);

        rows(rr,:)=[ ...
            Nlist(il),rRat(ir), ...
            KI0,KII0,ratio0, ...
            rad2deg(thetaRoot),rad2deg(thetaBracketRoot), ...
            oddDef,evenDef,fitRel, ...
            min(rin(il,:,ir)),max(rin(il,:,ir)), ...
            min(htip(il,:)),max(htip(il,:)), ...
            max(uxsym(il,:)), ...
            rad2deg(Stage1{il}.I.phi_star)];
    end
end

T=array2table(rows,'VariableNames',{ ...
    'Npoly','r_outer_over_a0', ...
    'KI_theta0','KII_theta0','KII_over_KI_theta0', ...
    'theta_root_linear_deg','theta_root_bracket_deg', ...
    'KII_odd_defect_rel','KI_even_defect_rel','linear_fit_rel_residual', ...
    'r_inner_min','r_inner_max','h_tip_min','h_tip_max', ...
    'max_ux_symmetry','stage1_phi_fit_deg'});

fprintf('\nSTAGE-II SYMMETRY SUMMARY\n');
disp(T);

% Domain spread at theta=0 for each mesh.
Drows=nan(nL,7);
for il=1:nL
    [~,i0]=min(abs(thetaDeg));
    k0=squeeze(KI(il,i0,:));
    q0=squeeze(KII(il,i0,:));
    roots=T.theta_root_linear_deg(T.Npoly==Nlist(il));

    Drows(il,:)=[ ...
        Nlist(il), ...
        max(k0)-min(k0), ...
        (max(k0)-min(k0))/max(abs(k0)), ...
        max(abs(q0)), ...
        max(abs(q0./k0)), ...
        max(roots)-min(roots), ...
        max(abs(roots))];
end

Td=array2table(Drows,'VariableNames',{ ...
    'Npoly','KI_domain_abs_spread','KI_domain_rel_spread', ...
    'max_abs_KII_theta0','max_abs_KII_over_KI_theta0', ...
    'theta_root_domain_spread_deg','max_abs_theta_root_deg'});

fprintf('\nEDI DOMAIN-SPREAD SUMMARY\n');
disp(Td);

if logical(O.Plot)
    local_plot_KII(thetaDeg,KII,KI,Nlist,rRat);
    local_plot_theta0(T);
    local_plot_root(T);
end

Out=struct();
Out.settings=O;
Out.NpolyList=Nlist;
Out.thetaDeg=thetaDeg;
Out.rOuterOverA0=rRat;
Out.KI=KI;
Out.KII=KII;
Out.rInner=rin;
Out.hTip=htip;
Out.uxSymmetry=uxsym;
Out.Stage1=Stage1;
Out.Cases=Cases;
Out.summary=T;
Out.domainSummary=Td;

fprintf('\nSTEP 18 completed.\n');
fprintf(['Gate: KII(0)/KI should approach zero, KII should be odd and KI even ', ...
    'in theta, and the signed-EDI root should converge to theta=0 with ', ...
    'small EDI-domain spread.\n']);
end


function root=local_sign_change_root(theta,y)
root=NaN;
best=inf;
for k=1:numel(theta)-1
    if y(k)==0
        cand=theta(k);
    elseif y(k)*y(k+1)<0
        cand=theta(k)-y(k)*(theta(k+1)-theta(k))/(y(k+1)-y(k));
    else
        continue;
    end
    if abs(cand)<best
        best=abs(cand);
        root=cand;
    end
end
end


function H=local_tip_mesh_scale(mesh,tip)
X=mesh.coord3;
T3=mesh.connect3;
d=sqrt(sum((X-tip).^2,2));
dmin=min(d);
tol=max(1e-12,1e-8*max(1,max(abs(X(:)))));
tipNodes=find(d<=dmin+tol);

hit=any(ismember(T3,tipNodes),2);
Te=T3(hit,:);
L=[];
for k=1:size(Te,1)
    P=X(Te(k,:),:);
    L=[L,norm(P(2,:)-P(1,:)),norm(P(3,:)-P(2,:)),norm(P(1,:)-P(3,:))]; %#ok<AGROW>
end
L=L(isfinite(L)&L>tol);

H=struct();
H.min=min(L);
H.median=median(L);
H.max=max(L);
H.nTipElements=size(Te,1);
H.tipNodeDistance=dmin;
end


function local_plot_KII(thetaDeg,KII,KI,Nlist,rRat)
for il=1:numel(Nlist)
    figure('Name',sprintf('Step 18: signed KII, Npoly=%d',Nlist(il)),'Color','w');
    clf; hold on; box on; grid on
    for ir=1:numel(rRat)
        q=squeeze(KII(il,:,ir));
        k=squeeze(KI(il,:,ir));
        plot(thetaDeg,q./k,'-o','LineWidth',1.1, ...
            'DisplayName',sprintf('r_o/a_0=%.2f',rRat(ir)));
    end
    xline(0,'k--','HandleVisibility','off');
    yline(0,'k:','HandleVisibility','off');
    xlabel('\theta [deg]');
    ylabel('K_{II}/K_I');
    title(sprintf('Centered half-domain signed EDI, Npoly=%d',Nlist(il)));
    legend('Location','best');
end
end


function local_plot_theta0(T)
figure('Name','Step 18: KII at theta=0','Color','w');
clf; hold on; box on; grid on
N=unique(T.Npoly);
for k=1:numel(N)
    Q=T(T.Npoly==N(k),:);
    plot(Q.r_outer_over_a0,abs(Q.KII_over_KI_theta0),'-o','LineWidth',1.1, ...
        'DisplayName',sprintf('Npoly=%d',N(k)));
end
set(gca,'YScale','log');
xlabel('r_{outer}/a_0');
ylabel('|K_{II}(0)/K_I(0)|');
title('Symmetry residual at normal extension');
legend('Location','best');
end


function local_plot_root(T)
figure('Name','Step 18: local-symmetry root','Color','w');
clf; hold on; box on; grid on
N=unique(T.Npoly);
for k=1:numel(N)
    Q=T(T.Npoly==N(k),:);
    plot(Q.r_outer_over_a0,Q.theta_root_linear_deg,'-o','LineWidth',1.1, ...
        'DisplayName',sprintf('Npoly=%d',N(k)));
end
yline(0,'k--','HandleVisibility','off');
xlabel('r_{outer}/a_0');
ylabel('\theta_* [deg]');
title('Signed-EDI local-symmetry root');
legend('Location','best');
end
