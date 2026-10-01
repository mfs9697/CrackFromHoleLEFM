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
