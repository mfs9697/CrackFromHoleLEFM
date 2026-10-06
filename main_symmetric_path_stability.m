function Results = main_symmetric_path_stability(varargin)
%MAIN_SYMMETRIC_PATH_STABILITY Full-domain symmetry-breaking experiment.
%
% Scientific question
% -------------------
% A centered circular hole under symmetric remote-y tension admits the exact
% straight crack path on y=0. This driver asks whether that fixed path is
% locally stable to a small angular perturbation of the SECOND 4-mm segment.
%
% The first segment is always prescribed by symmetry:
%
%   P0 = [A/2+R,0]
%   P1 = P0 + Delta a*[1,0]
%   theta_1 = 0
%
% For each prescribed probe theta_2 in a symmetric set (default
% [-0.1 -0.05 0 +0.05 +0.1] deg), the driver constructs
%
%   P2(theta_2) = P1 + Delta a*[cos(theta_2),sin(theta_2)]
%
% and performs exactly one qualified physical solve at P2. MTS then predicts
% theta_3. No angle is selected and no third segment is generated.
%
% Symmetry expectations
% ---------------------
%   KI(+theta)       =  KI(-theta)
%   KII(+theta)      = -KII(-theta)
%   DeltaTheta3(+t)  = -DeltaTheta3(-t)
%   theta3(+theta)   = -theta3(-theta)
%   KII(0) = 0, theta3(0) = 0
%
% Local stability
% ---------------
% Define F(theta_2)=theta_3 predicted by one MTS update. The symmetric path
% is a fixed point F(0)=0. A finite-difference/least-squares estimate of
%
%       m = dF/dtheta |_(theta=0)
%
% is reported. |m|<1 is locally restoring, |m|>1 amplifying, and m<0 with
% |m|<1 corresponds to damped sign-alternating correction.
%
% IMPORTANT: this test uses the FULL plate. A half-domain symmetry model
% cannot represent the +/- perturbations without constraining the crack.
%
% Usage
% -----
%   R = main_symmetric_path_stability('AllowPhysicalSolves',true);
%
% Options
% -------
%   ThetaProbeDeg       default [-0.1 -0.05 0 0.05 0.1]
%   AllowPhysicalSolves default false
%   RunSynthetic        default true
%   FastEDI             default false
%   ReuseCandidates     default true
%   ExteriorVerbose     default false
%   Plot                default true
%   OutputDir           default verification/crack_path/symmetric_stability
%
% No production crack-path rule is changed by this experiment.

    ip=inputParser;
    addParameter(ip,'ThetaProbeDeg',[-0.1 -0.05 0 0.05 0.1], ...
        @(x)isnumeric(x)&&isvector(x)&&all(isfinite(x)));
    addParameter(ip,'AllowPhysicalSolves',false,@(x)islogical(x)&&isscalar(x));
    addParameter(ip,'RunSynthetic',true,@(x)islogical(x)&&isscalar(x));
    addParameter(ip,'FastEDI',false,@(x)islogical(x)&&isscalar(x));
    addParameter(ip,'ReuseCandidates',true,@(x)islogical(x)&&isscalar(x));
    addParameter(ip,'ExteriorVerbose',false,@(x)islogical(x)&&isscalar(x));
    addParameter(ip,'Plot',true,@(x)islogical(x)&&isscalar(x));
    addParameter(ip,'OutputDir','',@(x)ischar(x)||isstring(x));
    parse(ip,varargin{:});
    opt=ip.Results;

    root=fileparts(mfilename('fullpath'));
    addpath(genpath(root));

    thetaProbe=sort(opt.ThetaProbeDeg(:));
    local_validate_probe(thetaProbe);

    R0=make_symmetric_crack_path_state();
    C=R0.C;
    S0=R0.summary;
    da=S0.a0_reserved_m;
    P0=[S0.x_star_m,S0.y_star_m];
    nMat=[S0.nmat_x,S0.nmat_y];
    tHat=[S0.that_x,S0.that_y];

    P1=P0+da*nMat;

    outDir=char(opt.OutputDir);
    if isempty(strtrim(outDir))
        outDir=fullfile(root,'verification','crack_path','symmetric_stability');
    elseif ~local_is_absolute_path(outDir)
        outDir=fullfile(root,outDir);
    end
    if exist(outDir,'dir')~=7,mkdir(outDir);end

    fprintf('\n============================================================\n');
    fprintf('SYMMETRIC CRACK-PATH LOCAL-STABILITY EXPERIMENT\n');
    fprintf('============================================================\n');
    fprintf('  Full plate        : A=%.3f m, B=%.3f m\n',C.A,C.B);
    fprintf('  centered hole     : [%.3f, %.3f] m, R=%.3f m\n', ...
        C.hole.center(1),C.hole.center(2),C.hole.r);
    fprintf('  P0                : [%.12g, %.12g] m\n',P0(1),P0(2));
    fprintf('  P1                : [%.12g, %.12g] m\n',P1(1),P1(2));
    fprintf('  increment         : %.9f mm\n',1e3*da);
    fprintf('  theta_2 probes    : %s deg\n',mat2str(thetaProbe.',8));
    fprintf('  physical solves   : %d\n',logical(opt.AllowPhysicalSolves));
    fprintf('  output directory  : %s\n',outDir);
    fprintf(['  NOTE: no Stage-I solve is needed; P0 and the local frame ', ...
        'are exact by symmetry.\n']);

    n=numel(thetaProbe);
    rows=nan(n,17);
    stepResults=cell(n,1);
    qualification=cell(n,1);
    paths=cell(n,1);

    for i=1:n
        th=thetaProbe(i);
        e2=cosd(th)*nMat+sind(th)*tHat;
        e2=e2/norm(e2);
        P2=P1+da*e2;
        path=[P0;P1;P2];
        paths{i}=path;

        tag=local_angle_tag(th);
        candidateFile=fullfile(outDir,['theta_' tag '_candidate.mat']);
        qualFile=fullfile(outDir,['theta_' tag '_qualification_small.mat']);
        checkpointFile=fullfile(outDir,['theta_' tag '_physical_solved.mat']);
        physicalFile=fullfile(outDir,['theta_' tag '_physical_small.mat']);

        fprintf('\n------------------------------------------------------------\n');
        fprintf('PROBE theta_2 = %+.9f deg\n',th);
        fprintf('  P2 = [%.12g, %.12g] m\n',P2(1),P2(2));

        candidate=[];
        if opt.ReuseCandidates && exist(candidateFile,'file')==2
            d=load(candidateFile,'candidate');
            if isfield(d,'candidate')&&isstruct(d.candidate)&& ...
                    isfield(d.candidate,'path')&& ...
                    size(d.candidate.path,1)==3&& ...
                    norm(d.candidate.path-path,'fro')<=2e-12&& ...
                    isfield(d.candidate,'scientificallyReadyForIncrementalPhysicalSolve')&& ...
                    logical(d.candidate.scientificallyReadyForIncrementalPhysicalSolve)
                candidate=d.candidate;
                fprintf('  Reusing qualified candidate: %s\n',candidateFile);
            end
        end

        if isempty(candidate)
            Q=qualify_incremental_crack_candidate(path, ...
                'FrozenState',R0, ...
                'RunSynthetic',opt.RunSynthetic, ...
                'FastEDI',opt.FastEDI, ...
                'AllowStraightPath',abs(th)<=1e-14, ...
                'ExteriorVerbose',opt.ExteriorVerbose, ...
                'SaveCandidate',true, ...
                'CandidateFile',candidateFile, ...
                'SaveCompact',true, ...
                'CompactFile',qualFile, ...
                'Plot',false);
            assert(Q.pass,'symstab:Qualification', ...
                'Qualification failed for theta_2=%+.9g deg.',th);
            candidate=Q.candidate;
            qualification{i}=Q.summary;
        elseif exist(qualFile,'file')==2
            qd=load(qualFile,'Small');
            if isfield(qd,'Small')&&isfield(qd.Small,'summary')
                qualification{i}=qd.Small.summary;
            end
        end

        R=solve_incremental_crack_tip(candidate, ...
            'FrozenState',R0, ...
            'AllowSolve',opt.AllowPhysicalSolves, ...
            'FastEDI',opt.FastEDI, ...
            'CheckpointFile',checkpointFile, ...
            'SaveFile',physicalFile);

        assert(R.pass,'symstab:Physical', ...
            'Physical solve/postprocessing failed for theta_2=%+.9g deg.',th);
        stepResults{i}=R;

        qsum=qualification{i};
        clearance=NaN; minAngle=NaN; maxNeighbor=NaN;
        if istable(qsum)&&height(qsum)==1
            if ismember('physical_clearance_m',qsum.Properties.VariableNames)
                clearance=qsum.physical_clearance_m;
            end
            if ismember('min_angle_deg',qsum.Properties.VariableNames)
                minAngle=qsum.min_angle_deg;
            end
            if ismember('max_neighbor_ratio',qsum.Properties.VariableNames)
                maxNeighbor=qsum.max_neighbor_ratio;
            end
        end

        rows(i,:)=[ ...
            th,P2(1),P2(2), ...
            R.EDI.KI_unit,R.EDI.KII_unit,R.EDI.KII_over_KI, ...
            R.deltaThetaNextDeg,R.thetaNextDeg, ...
            R.solverInfo.iter,R.solverInfo.relres,R.solverInfo.trueRelResidual, ...
            R.EDI.EDI_elements,clearance,minAngle,maxNeighbor, ...
            double(candidate.pathIsNonStraight),double(R.newSolve)];
    end

    T=array2table(rows,'VariableNames',{ ...
        'theta2_deg','tip_x_m','tip_y_m', ...
        'KI_unit','KII_unit','KII_over_KI', ...
        'delta_theta3_MTS_deg','theta3_pred_deg', ...
        'PCG_iterations','PCG_relres','true_rel_residual', ...
        'EDI_elements','physical_clearance_m','min_angle_deg', ...
        'max_neighbor_ratio','path_is_nonstraight','newSolve'});
    T.pass=true(height(T),1);

    [Parity,Linearization,Control]=local_analyze(T);

    fprintf('\n============================================================\n');
    fprintf('SYMMETRIC STABILITY RESULTS\n');
    fprintf('============================================================\n');
    disp(T);

    fprintf('\nPAIRED SYMMETRY DIAGNOSTICS\n');
    disp(Parity);

    fprintf('\nSTRAIGHT CONTROL theta_2=0\n');
    disp(Control);

    fprintf('\nLINEARIZED ONE-STEP MAP\n');
    disp(struct2table(Linearization,'AsArray',true));
    fprintf('  F(theta_2)=theta_3, fitted multiplier m = %+.8g\n', ...
        Linearization.theta3_over_theta_slope);
    fprintf('  provisional classification: %s\n',Linearization.classification);

    writetable(T,fullfile(outDir,'symmetric_stability_probes.csv'));
    writetable(Parity,fullfile(outDir,'symmetric_stability_parity.csv'));
    writetable(struct2table(Linearization,'AsArray',true), ...
        fullfile(outDir,'symmetric_stability_linearization.csv'));
    writetable(Control,fullfile(outDir,'symmetric_stability_straight_control.csv'));

    Figures=struct();
    if opt.Plot
        Figures=local_plot(T,Linearization,outDir);
    end

    Results=struct();
    Results.R0=R0;
    Results.config=C;
    Results.P0=P0;
    Results.P1=P1;
    Results.paths=paths;
    Results.Table=T;
    Results.Parity=Parity;
    Results.Control=Control;
    Results.Linearization=Linearization;
    Results.stepResults=stepResults;
    Results.qualification=qualification;
    Results.figures=Figures;
    Results.outputDir=outDir;
    Results.experiment=[ ...
        'Full-domain centered-hole perturbation test of the straight ', ...
        'symmetric MTS fixed path after the first 4-mm segment.'];

    save(fullfile(outDir,'symmetric_stability_result.mat'),'Results','-v7');

    fprintf('\nSYMMETRIC STABILITY EXPERIMENT COMPLETE.\n');
    fprintf('  Result MAT: %s\n', ...
        fullfile(outDir,'symmetric_stability_result.mat'));
end


function local_validate_probe(theta)
    tol=1e-12;
    if ~any(abs(theta)<=tol)
        error('symstab:NeedZero','ThetaProbeDeg must include theta=0.');
    end
    pos=theta(theta>tol);
    neg=theta(theta<-tol);
    if isempty(pos)||isempty(neg)
        error('symstab:NeedPairs','ThetaProbeDeg must contain +/- perturbations.');
    end
    for k=1:numel(pos)
        if ~any(abs(neg+pos(k))<=tol)
            error('symstab:UnpairedProbe', ...
                'Positive probe %+g deg has no negative mirror.',pos(k));
        end
    end
    for k=1:numel(neg)
        if ~any(abs(pos+neg(k))<=tol)
            error('symstab:UnpairedProbe', ...
                'Negative probe %+g deg has no positive mirror.',neg(k));
        end
    end
end


function [Parity,L,C0]=local_analyze(T)
    tol=1e-12;
    pos=unique(T.theta2_deg(T.theta2_deg>tol));
    prow=nan(numel(pos),10);

    for k=1:numel(pos)
        a=pos(k);
        ip=find(abs(T.theta2_deg-a)<=tol,1);
        im=find(abs(T.theta2_deg+a)<=tol,1);

        KIp=T.KI_unit(ip); KIm=T.KI_unit(im);
        qp=T.KII_over_KI(ip); qm=T.KII_over_KI(im);
        dp=T.delta_theta3_MTS_deg(ip); dm=T.delta_theta3_MTS_deg(im);
        fp=T.theta3_pred_deg(ip); fm=T.theta3_pred_deg(im);

        qEven=0.5*(qp+qm);
        qOdd=0.5*(qp-qm);
        dEven=0.5*(dp+dm);
        fEven=0.5*(fp+fm);

        prow(k,:)=[ ...
            a,KIm,KIp, ...
            abs(KIp-KIm)/max(0.5*(abs(KIp)+abs(KIm)),eps), ...
            qEven,qOdd,abs(qEven)/max(abs(qOdd),eps), ...
            dEven,fEven,(fp-fm)/(2*a)];
    end

    Parity=array2table(prow,'VariableNames',{ ...
        'theta_abs_deg','KI_minus','KI_plus','KI_pair_rel_difference', ...
        'ratio_even','ratio_odd','ratio_even_over_odd', ...
        'delta_theta_even_deg','theta3_even_deg','map_multiplier_pair'});

    iz=find(abs(T.theta2_deg)<=tol,1);
    C0=T(iz,{'theta2_deg','KI_unit','KII_unit','KII_over_KI', ...
        'delta_theta3_MTS_deg','theta3_pred_deg','PCG_iterations', ...
        'true_rel_residual','EDI_elements'});

    use=abs(T.theta2_deg)>tol & ...
        isfinite(T.theta2_deg) & isfinite(T.KII_over_KI) & ...
        isfinite(T.delta_theta3_MTS_deg) & isfinite(T.theta3_pred_deg);
    x=T.theta2_deg(use);

    L=struct();
    L.nNonzeroProbes=nnz(use);
    L.ratio_per_degree_slope=(x.'*T.KII_over_KI(use))/(x.'*x);
    L.deltaTheta3_over_theta_slope= ...
        (x.'*T.delta_theta3_MTS_deg(use))/(x.'*x);
    L.theta3_over_theta_slope= ...
        (x.'*T.theta3_pred_deg(use))/(x.'*x);
    L.straight_ratio=T.KII_over_KI(iz);
    L.straight_deltaTheta3_deg=T.delta_theta3_MTS_deg(iz);
    L.straight_theta3_deg=T.theta3_pred_deg(iz);

    m=L.theta3_over_theta_slope;
    if abs(m)<1
        if m>=0
            cls='locally restoring (same-sign decay; finite-probe estimate)';
        else
            cls='locally restoring (sign-alternating decay; finite-probe estimate)';
        end
    elseif abs(m)>1
        cls='locally amplifying / symmetry-breaking (finite-probe estimate)';
    else
        cls='neutral at current finite-probe resolution';
    end
    L.classification=cls;
end


function F=local_plot(T,L,outDir)
    F=struct();

    F.modeMixity=figure('Color','w','Name','Symmetric stability: mode mixity');
    plot(T.theta2_deg,T.KII_over_KI,'-o','LineWidth',1.4,'MarkerSize',5);
    hold on; xline(0,'--'); yline(0,'--');
    xlabel('\theta_2 [deg]'); ylabel('K_{II}/K_I');
    title('Mode mixity under symmetric angular perturbations');
    grid on; box on; hold off;

    F.turn=figure('Color','w','Name','Symmetric stability: MTS correction');
    plot(T.theta2_deg,T.delta_theta3_MTS_deg,'-o', ...
        'LineWidth',1.4,'MarkerSize',5);
    hold on; xline(0,'--'); yline(0,'--');
    xlabel('\theta_2 [deg]'); ylabel('\Delta\theta_3^{MTS} [deg]');
    title('MTS correction after the perturbed second segment');
    grid on; box on; hold off;

    F.map=figure('Color','w','Name','Symmetric stability: one-step map');
    plot(T.theta2_deg,T.theta3_pred_deg,'-o', ...
        'LineWidth',1.4,'MarkerSize',5,'DisplayName','F(\theta_2)=\theta_3');
    hold on;
    xx=[min(T.theta2_deg),max(T.theta2_deg)];
    plot(xx,xx,'--','DisplayName','identity');
    plot(xx,L.theta3_over_theta_slope*xx,':', ...
        'LineWidth',1.2,'DisplayName',sprintf('origin fit, m=%+.4f', ...
        L.theta3_over_theta_slope));
    xline(0,'--','HandleVisibility','off'); yline(0,'--','HandleVisibility','off');
    xlabel('\theta_2 [deg]'); ylabel('\theta_3 predicted [deg]');
    title('One-step crack-direction map near the symmetric path');
    grid on; box on; legend('Location','best'); hold off;

    F.KI=figure('Color','w','Name','Symmetric stability: KI parity');
    plot(T.theta2_deg,T.KI_unit,'-o','LineWidth',1.4,'MarkerSize',5);
    hold on; xline(0,'--');
    xlabel('\theta_2 [deg]'); ylabel('K_I at unit traction [MPa sqrt(m)]');
    title('Mode-I response should be even in the perturbation angle');
    grid on; box on; hold off;

    local_export(F.modeMixity,outDir,'symmetric_01_mode_mixity');
    local_export(F.turn,outDir,'symmetric_02_MTS_turn');
    local_export(F.map,outDir,'symmetric_03_one_step_map');
    local_export(F.KI,outDir,'symmetric_04_KI_parity');
end


function local_export(fig,outDir,stem)
    print(fig,fullfile(outDir,[stem '.eps']),'-depsc','-painters');
    exportgraphics(fig,fullfile(outDir,[stem '.png']),'Resolution',300);
end


function tag=local_angle_tag(theta)
    if abs(theta)<=1e-14
        prefix='z';
    elseif theta<0
        prefix='m';
    else
        prefix='p';
    end
    tag=[prefix strrep(sprintf('%.3f',abs(theta)),'.','p')];
end


function tf=local_is_absolute_path(p)
    p=char(p);
    if isempty(p),tf=false;return,end
    tf=startsWith(p,filesep) || startsWith(p,'\\') || ...
        ~isempty(regexp(p,'^[A-Za-z]:[\\/]','once'));
end
