function Results = main_check_centered_hole_EDI_parity(varargin)
%MAIN_CHECK_CENTERED_HOLE_EDI_PARITY  Independent sign/parity gate for EDI SIFs.
%
% Checks, on the centered-hole benchmark,
%   KII(theta=0) ~ 0
%   KII(-theta) ~ -KII(+theta)
% and compares the EDI result with the historical mirror-J extractor.
%
% This is a verification gate only; it does not advance the crack.

    ip=inputParser;
    addParameter(ip,'ThetaDeg',[-10 -5 0 5 10],@isnumeric);
    addParameter(ip,'rOverA0',[0.50 0.60 0.70],@isnumeric);
    addParameter(ip,'InnerFraction',0.20,@(x)isnumeric(x)&&isscalar(x)&&x>=0&&x<1);
    addParameter(ip,'MeshScale',0.40,@(x)isnumeric(x)&&isscalar(x)&&x>0);
    addParameter(ip,'a0Factor',4.0,@(x)isnumeric(x)&&isscalar(x)&&x>0);
    addParameter(ip,'SaveOutputs',true,@(x)islogical(x)||isnumeric(x));
    addParameter(ip,'OutputDir',fullfile('verification','sif','outputs'),@(x)ischar(x)||isstring(x));
    addParameter(ip,'Verbose',true,@(x)islogical(x)||isnumeric(x));
    parse(ip,varargin{:});

    addpath(genpath(pwd));

    C0=cfg_hole_initiation();
    C0.hole.center=[0.5*C0.A,0];
    C0.holes={C0.hole};

    G=geom_hole_only(C0);
    S1=solve_hole_only(C0,G,'lambda',1.0);
    B=sample_hole_boundary_stress(C0,G,S1);
    % For a parity benchmark, do NOT use the numerically selected first
    % maximizer on the sampled hole boundary: the discrete sampler may return
    % phi = 2*pi-dphi instead of exactly zero and thereby break the intended
    % reflection symmetry before the SIF extractor is even called.
    %
    % Freeze the geometrically exact rightmost point of the centered circle.
    I=struct();
    I.idx_star=NaN;
    I.phi_star=0;
    I.x_star=C0.hole.center+[C0.hole.r,0];
    I.n_mat_star=[1,0];
    I.n_hole_star=[-1,0];
    I.t_hat_star=[0,1];
    I.sig_tt_unit=NaN;
    I.sig_tt_pos_unit=NaN;
    I.lambda_ini=NaN;
    I.sig_applied_ini=NaN;
    I.all_max_idx=[];
    I.selection_rule='exact_centered_hole_parity_point';

    C=local_scale_mesh_config(C0,ip.Results.MeshScale);
    C.a0=ip.Results.a0Factor*C0.a0;

    thetaDeg=ip.Results.ThetaDeg(:).';
    rfac=ip.Results.rOverA0(:).';
    innerFrac=ip.Results.InnerFraction;

    rows=[];
    for td=thetaDeg
        th=deg2rad(td);

        [G2,~,~,Mc]=build_stage2_cracked_mesh_for_theta( ...
            C,I,th,'PlotGeom',false,'PlotMesh',false,'PlotCollapsed',false);
        S2=solve_cracked_LEFM(C,Mc);

        if isfield(Mc,'crack') && isfield(Mc.crack,'Pmid') && ~isempty(Mc.crack.Pmid)
            V=Mc.crack.Pmid;
        else
            V=G2.crack.polyline;
        end

        mat=S2.mat;
        if ~isfield(mat,'Dmat') && isfield(mat,'D'),mat.Dmat=mat.D;end

        for rr=rfac
            r=rr*G2.crack.a0;

            [KIold,KIIold,Dbg]=SIF_LEFM_circle2_debug( ...
                S2.mesh,S2.U,V,mat,r,'nthet',240,'plot',false,'verbose',false);

            M=SIF_mesh_mirror_metrics(S2.mesh,V,r,'theta',Dbg.theta);

            dom=struct('r_inner',innerFrac*r,'r_outer',r);
            [KIedi,KIIedi,A]=SIF_LEFM_interaction_EDI( ...
                S2.mesh,S2.U,V,mat,dom,'AuxDerivativeMode','analytic');

            rows=[rows; ...
                td,rr,r,M.A_h_median, ...
                KIold,KIIold,KIIold/KIold, ...
                KIedi,KIIedi,KIIedi/KIedi, ...
                KIold-KIedi,KIIold-KIIedi, ...
                A.auxConsistency.maxI,A.auxConsistency.maxII]; %#ok<AGROW>

            if logical(ip.Results.Verbose)
                fprintf(['theta=%+6.1f deg r/a0=%.2f Ah=%.3e | ', ...
                    'old KII/KI=%+.4e EDI KII/KI=%+.4e\n'], ...
                    td,rr,M.A_h_median,KIIold/KIold,KIIedi/KIedi);
            end
        end
    end

    T=array2table(rows,'VariableNames',{ ...
        'thetaDeg','rOuter_over_a0','rOuter','A_h_median', ...
        'KI_old','KII_old','ratio_old','KI_EDI','KII_EDI','ratio_EDI', ...
        'dKI_old_minus_EDI','dKII_old_minus_EDI', ...
        'auxMismatchMaxI','auxMismatchMaxII'});

    S=local_theta_summary(T);
    P=local_parity(S);

    Results=struct();
    Results.C0=C0;
    Results.C=C;
    Results.Stage1=struct('G',G,'S1',S1,'B',B,'I',I);
    Results.Table=T;
    Results.Summary=S;
    Results.Parity=P;

    if logical(ip.Results.SaveOutputs)
        outdir=char(ip.Results.OutputDir);
        if ~exist(outdir,'dir'),mkdir(outdir);end
        writetable(T,fullfile(outdir,'centered_hole_old_vs_edi_parity.csv'));
        writetable(S,fullfile(outdir,'centered_hole_old_vs_edi_parity_summary.csv'));
        writetable(P,fullfile(outdir,'centered_hole_old_vs_edi_parity_pairs.csv'));
        save(fullfile(outdir,'centered_hole_old_vs_edi_parity.mat'),'Results','-v7.3');
    end
end


function S=local_theta_summary(T)
    th=unique(T.thetaDeg);
    rows=[];
    for k=1:numel(th)
        idx=abs(T.thetaDeg-th(k))<1e-12;
        rows=[rows;th(k), ...
            median(T.KI_old(idx)),median(T.KII_old(idx)),median(T.ratio_old(idx)), ...
            median(T.KI_EDI(idx)),median(T.KII_EDI(idx)),median(T.ratio_EDI(idx)), ...
            max(T.A_h_median(idx))]; %#ok<AGROW>
    end
    S=array2table(rows,'VariableNames',{ ...
        'thetaDeg','KI_old','KII_old','ratio_old', ...
        'KI_EDI','KII_EDI','ratio_EDI','A_h_max'});
end


function P=local_parity(S)
    pos=sort(unique(abs(S.thetaDeg(S.thetaDeg>0))));
    rows=[];
    for a=pos(:).'
        im=find(abs(S.thetaDeg+a)<1e-12,1);
        ip=find(abs(S.thetaDeg-a)<1e-12,1);
        if isempty(im)||isempty(ip),continue;end

        om=S.KII_old(im); op=S.KII_old(ip);
        em=S.KII_EDI(im); ep=S.KII_EDI(ip);

        rows=[rows;a, ...
            om,op,0.5*(om+op),0.5*(op-om), ...
            em,ep,0.5*(em+ep),0.5*(ep-em)]; %#ok<AGROW>
    end

    P=array2table(rows,'VariableNames',{ ...
        'absThetaDeg', ...
        'KIIoldMinus','KIIoldPlus','KIIoldEven','KIIoldOdd', ...
        'KIIEDIMinus','KIIEDIPlus','KIIEDIEven','KIIEDIOdd'});

    i0=find(abs(S.thetaDeg)<1e-12,1);
    if ~isempty(i0)
        P.Properties.UserData.KII0_old=S.KII_old(i0);
        P.Properties.UserData.KII0_EDI=S.KII_EDI(i0);
        P.Properties.UserData.ratio0_old=S.ratio_old(i0);
        P.Properties.UserData.ratio0_EDI=S.ratio_EDI(i0);
    end
end


function C=local_scale_mesh_config(C0,s)
    C=C0;
    if isfield(C,'mesh1')
        if isfield(C.mesh1,'hmin'),C.mesh1.hmin=s*C0.mesh1.hmin;end
        if isfield(C.mesh1,'hmax'),C.mesh1.hmax=s*C0.mesh1.hmax;end
    end
    if isfield(C,'mesh2')
        if isfield(C.mesh2,'hmin'),C.mesh2.hmin=s*C0.mesh2.hmin;end
        if isfield(C.mesh2,'hmax'),C.mesh2.hmax=s*C0.mesh2.hmax;end
        if isfield(C.mesh2,'hcrack'),C.mesh2.hcrack=s*C0.mesh2.hcrack;end
        if isfield(C.mesh2,'hhole'),C.mesh2.hhole=s*C0.mesh2.hhole;end
    end
end
