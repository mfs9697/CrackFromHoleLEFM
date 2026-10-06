function Results = main_symmetric_stability_geometry_check(varargin)
%MAIN_SYMMETRIC_STABILITY_GEOMETRY_CHECK
% Qualification-only preflight for the centered-hole stability experiment.
%
% No physical FEM solve is performed. The five default two-segment paths are
% built, meshed, and passed through all structural and prescribed-Williams
% qualification gates. Candidate files are saved for reuse by
% main_symmetric_path_stability.

    ip=inputParser;
    addParameter(ip,'ThetaProbeDeg',[-0.1 -0.05 0 0.05 0.1], ...
        @(x)isnumeric(x)&&isvector(x)&&all(isfinite(x)));
    addParameter(ip,'RunSynthetic',true,@(x)islogical(x)&&isscalar(x));
    addParameter(ip,'FastEDI',false,@(x)islogical(x)&&isscalar(x));
    addParameter(ip,'ExteriorVerbose',false,@(x)islogical(x)&&isscalar(x));
    addParameter(ip,'OutputDir','',@(x)ischar(x)||isstring(x));
    parse(ip,varargin{:});
    opt=ip.Results;

    root=fileparts(mfilename('fullpath'));
    addpath(genpath(root));

    thetaProbe=sort(opt.ThetaProbeDeg(:));
    local_validate_probe(thetaProbe);

    R0=make_symmetric_crack_path_state();
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
    fprintf('SYMMETRIC STABILITY: GEOMETRY / EDI QUALIFICATION ONLY\n');
    fprintf('============================================================\n');
    fprintf('  NO physical FEM solve will be performed.\n');
    fprintf('  probes: %s deg\n',mat2str(thetaProbe.',8));
    fprintf('  P0 = [%.12g, %.12g] m\n',P0(1),P0(2));
    fprintf('  P1 = [%.12g, %.12g] m\n',P1(1),P1(2));

    n=numel(thetaProbe);
    rows=nan(n,13);

    for i=1:n
        th=thetaProbe(i);
        e2=cosd(th)*nMat+sind(th)*tHat;
        e2=e2/norm(e2);
        P2=P1+da*e2;
        path=[P0;P1;P2];

        tag=local_angle_tag(th);
        candidateFile=fullfile(outDir,['theta_' tag '_candidate.mat']);
        qualFile=fullfile(outDir,['theta_' tag '_qualification_small.mat']);

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

        assert(Q.pass,'symstabgeom:Qualification', ...
            'Qualification failed for theta_2=%+.9g deg.',th);

        S=Q.summary;
        rows(i,:)=[ ...
            th,P2(1),P2(2), ...
            S.T3_nodes,S.T3_elements,S.T6_nodes, ...
            S.physical_clearance_m,S.EDI_elements, ...
            S.min_angle_deg,S.max_neighbor_ratio, ...
            S.recovery_matrix_error,S.tiny_mixed_KII_rel_error, ...
            double(Q.candidate.pathIsNonStraight)];
    end

    T=array2table(rows,'VariableNames',{ ...
        'theta2_deg','tip_x_m','tip_y_m', ...
        'T3_nodes','T3_elements','T6_nodes', ...
        'physical_clearance_m','EDI_elements', ...
        'min_angle_deg','max_neighbor_ratio', ...
        'recovery_matrix_error','tiny_mixed_KII_rel_error', ...
        'path_is_nonstraight'});
    T.pass=true(height(T),1);

    Mirror=local_mirror_table(T);

    fprintf('\nQUALIFICATION SUMMARY\n');
    disp(T);
    fprintf('\nMIRROR-GEOMETRY CHECK\n');
    disp(Mirror);

    writetable(T,fullfile(outDir,'symmetric_stability_geometry.csv'));
    writetable(Mirror,fullfile(outDir,'symmetric_stability_geometry_mirror.csv'));

    Results=struct();
    Results.R0=R0;
    Results.Table=T;
    Results.Mirror=Mirror;
    Results.outputDir=outDir;
    Results.pass=all(T.pass) && all(Mirror.pass);

    save(fullfile(outDir,'symmetric_stability_geometry_result.mat'), ...
        'Results','-v7');

    assert(Results.pass,'symstabgeom:Failed', ...
        'At least one symmetric geometry preflight gate failed.');

    fprintf('\nSYMMETRIC STABILITY GEOMETRY PREFLIGHT PASS.\n');
    fprintf('  Qualified candidates are ready for guarded physical solves.\n');
end


function M=local_mirror_table(T)
    tol=1e-12;
    pos=unique(T.theta2_deg(T.theta2_deg>tol));
    rows=nan(numel(pos),6);

    for k=1:numel(pos)
        a=pos(k);
        ip=find(abs(T.theta2_deg-a)<=tol,1);
        im=find(abs(T.theta2_deg+a)<=tol,1);

        dx=abs(T.tip_x_m(ip)-T.tip_x_m(im));
        ysum=abs(T.tip_y_m(ip)+T.tip_y_m(im));
        dn=abs(T.T3_nodes(ip)-T.T3_nodes(im));
        de=abs(T.T3_elements(ip)-T.T3_elements(im));

        % Equal mesh counts are NOT required: the unstructured full-domain
        % mesher is not forced to be a reflected mesh. Only the prescribed
        % path geometry must be exact mirrors.
        pass=dx<=2e-12 && ysum<=2e-12;

        rows(k,:)=[a,dx,ysum,dn,de,double(pass)];
    end

    M=array2table(rows,'VariableNames',{ ...
        'theta_abs_deg','tip_x_mismatch_m','tip_y_mirror_sum_m', ...
        'T3_node_count_difference','T3_element_count_difference','pass'});
    M.pass=logical(M.pass);
end


function local_validate_probe(theta)
    tol=1e-12;
    if ~any(abs(theta)<=tol)
        error('symstabgeom:NeedZero','ThetaProbeDeg must include theta=0.');
    end
    pos=theta(theta>tol);
    neg=theta(theta<-tol);
    if isempty(pos)||isempty(neg)
        error('symstabgeom:NeedPairs','ThetaProbeDeg must contain +/- perturbations.');
    end
    for k=1:numel(pos)
        if ~any(abs(neg+pos(k))<=tol)
            error('symstabgeom:UnpairedProbe', ...
                'Positive probe %+g deg has no negative mirror.',pos(k));
        end
    end
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
