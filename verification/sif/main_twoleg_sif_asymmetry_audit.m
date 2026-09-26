function Results = main_twoleg_sif_asymmetry_audit(varargin)
%MAIN_TWOLEG_SIF_ASYMMETRY_AUDIT  Old mirror-J versus interaction EDI.
%
% This is the main verification driver for the historical two-leg geometry.
% It does NOT advance the crack.  It compares SIF extraction methods over
% mesh refinement, deliberate local mesh asymmetry, and integration size.
%
% Name-value options
%   'MeshScaleList'       default [1.20 1.00 0.80]
%   'AsymmetryLevels'     default from cfg_twoleg_historical
%   'rOverDelta'          default [0.2 0.3 0.4 0.5 0.6]
%   'EDIInnerFraction'    default 0.20
%   'SaveOutputs'         true
%   'OutputDir'           verification/sif/outputs
%   'Verbose'             true

    ip=inputParser;
    addParameter(ip,'MeshScaleList',[1.20 1.00 0.80],@(x)isnumeric(x)&&all(x>0));
    addParameter(ip,'AsymmetryLevels',[],@isnumeric);
    addParameter(ip,'rOverDelta',[],@isnumeric);
    addParameter(ip,'EDIInnerFraction',0.20,@(x)isnumeric(x)&&isscalar(x)&&x>=0&&x<1);
    addParameter(ip,'SaveOutputs',true,@(x)islogical(x)||isnumeric(x));
    addParameter(ip,'OutputDir',fullfile('verification','sif','outputs'),@(x)ischar(x)||isstring(x));
    addParameter(ip,'Verbose',true,@(x)islogical(x)||isnumeric(x));
    parse(ip,varargin{:});

    C=cfg_twoleg_historical();
    meshScaleList=ip.Results.MeshScaleList(:).';
    asym=ip.Results.AsymmetryLevels;
    if isempty(asym), asym=C.asymmetry_levels; end
    asym=asym(:).';
    rr=ip.Results.rOverDelta;
    if isempty(rr), rr=C.rI_over_delta; end
    rr=rr(:).';
    innerFrac=ip.Results.EDIInnerFraction;
    verbose=logical(ip.Results.Verbose);

    if verbose
        fprintf('\n============================================================\n');
        fprintf('TWO-LEG SIF ASYMMETRY AUDIT\n');
        fprintf('theta1 = %.3f deg, theta2 = %.3f deg, delta/a = %.5f\n', ...
            C.theta1_deg,C.theta2_deg,C.delta/C.a);
        fprintf('IMPORTANT: second-tip radii are normalized by delta, not a.\n');
        fprintf('============================================================\n');
    end

    rows=[];
    debugStore={};
    irow=0;

    for ims=1:numel(meshScaleList)
        ms=meshScaleList(ims);

        for ia=1:numel(asym)
            alev=asym(ia);

            if verbose
                fprintf('\n--- meshScale %.3f | asymmetry request %.3f ---\n',ms,alev);
            end

            S=build_twoleg_sharp_benchmark(C, ...
                'MeshScale',ms,'AsymmetryLevel',alev,'Plot',false);

            for ir=1:numel(rr)
                rI=rr(ir)*C.delta;

                [KIold,KIIold,Dbg]=SIF_LEFM_circle2_debug( ...
                    S.mesh,S.U,S.V,S.mat,rI, ...
                    'nthet',C.nthet,'eps_th',C.eps_th, ...
                    'plot',false,'verbose',false);

                MM=SIF_mesh_mirror_metrics( ...
                    S.mesh,S.V,rI,'theta',Dbg.theta,'edge_tol',1e-4);

                domain=struct('r_inner',innerFrac*rI,'r_outer',rI);
                [KIedi,KIIedi,AE]=SIF_LEFM_interaction_EDI( ...
                    S.mesh,S.U,S.V,S.mat,domain, ...
                    'AuxDerivativeMode','analytic','Verbose',false);

                dKI=KIold-KIedi;
                dKII=KIIold-KIIedi;
                relVec=hypot(dKI,dKII)/max(hypot(KIedi,KIIedi),eps);

                irow=irow+1;
                rows(irow,:)=[ ...
                    ms,alev,S.meshInfo.asymmetryLevelUsed, ...
                    rr(ir),rI, ...
                    MM.A_h_median,MM.A_h_rms,MM.A_centroid_median, ...
                    MM.A_bary_median,MM.A_quality_median, ...
                    MM.frac_near_edge_P,MM.frac_near_edge_Q, ...
                    Dbg.diagnostics.min_baryP,Dbg.diagnostics.min_baryQ, ...
                    Dbg.diagnostics.J1_sign_changes,Dbg.diagnostics.J2_sign_changes, ...
                    Dbg.diagnostics.JI_over_absint,Dbg.diagnostics.JII_over_absint, ...
                    KIold,KIIold,KIedi,KIIedi,dKI,dKII,relVec, ...
                    AE.auxConsistency.maxI,AE.auxConsistency.maxII, ...
                    size(S.mesh.coord,1),size(S.mesh.connect,1)];

                debugStore{irow,1}=struct('old',Dbg,'meshMetric',MM,'edi',AE); %#ok<AGROW>

                if verbose
                    fprintf([' r/d=%.2f  Ah=%.3e | old=(%.6e,%+.6e) ', ...
                        'EDI=(%.6e,%+.6e) | relVec=%.3e\n'], ...
                        rr(ir),MM.A_h_median,KIold,KIIold,KIedi,KIIedi,relVec);
                end
            end
        end
    end

    T=array2table(rows,'VariableNames',{ ...
        'meshScale','asymmetryRequested','asymmetryUsed', ...
        'rI_over_delta','rI', ...
        'A_h_median','A_h_rms','A_centroid_median', ...
        'A_bary_median','A_quality_median', ...
        'fracNearEdgeP','fracNearEdgeQ','minBaryP','minBaryQ', ...
        'J1SignChanges','J2SignChanges','JI_over_absint','JII_over_absint', ...
        'KI_old','KII_old','KI_EDI','KII_EDI', ...
        'dKI_old_minus_EDI','dKII_old_minus_EDI','relVectorDifference', ...
        'auxMismatchMaxI','auxMismatchMaxII','nNode6','nElem6'});

    % Per-mesh summary over contour radii.
    Summary=local_summary(T);

    % Correlation is descriptive only; no causal interpretation is built in.
    good=isfinite(T.A_h_median)&isfinite(T.relVectorDifference);
    if nnz(good)>=3
        Corr=table( ...
            corr(T.A_h_median(good),T.relVectorDifference(good),'Type','Pearson'), ...
            corr(T.A_h_median(good),T.relVectorDifference(good),'Type','Spearman'), ...
            'VariableNames',{'Pearson_Ah_vs_error','Spearman_Ah_vs_error'});
    else
        Corr=table(NaN,NaN,'VariableNames', ...
            {'Pearson_Ah_vs_error','Spearman_Ah_vs_error'});
    end

    Results=struct();
    Results.C=C;
    Results.Table=T;
    Results.Summary=Summary;
    Results.Correlation=Corr;
    Results.DebugStore=debugStore;

    if logical(ip.Results.SaveOutputs)
        outdir=char(ip.Results.OutputDir);
        if ~exist(outdir,'dir'), mkdir(outdir); end
        writetable(T,fullfile(outdir,'twoleg_old_vs_edi.csv'));
        writetable(Summary,fullfile(outdir,'twoleg_old_vs_edi_summary.csv'));
        writetable(Corr,fullfile(outdir,'twoleg_asymmetry_error_correlation.csv'));
        save(fullfile(outdir,'twoleg_sif_asymmetry_audit.mat'),'Results','-v7.3');
    end
end


function S=local_summary(T)
    ms=unique(T.meshScale);
    aa=unique(T.asymmetryRequested);
    rows=[];

    for i=1:numel(ms)
        for j=1:numel(aa)
            idx=abs(T.meshScale-ms(i))<1e-12 & ...
                abs(T.asymmetryRequested-aa(j))<1e-12;
            if ~any(idx), continue; end

            KIo=T.KI_old(idx); KIIo=T.KII_old(idx);
            KIe=T.KI_EDI(idx); KIIe=T.KII_EDI(idx);

            rows=[rows; ...
                ms(i),aa(j),median(T.asymmetryUsed(idx)), ...
                median(T.A_h_median(idx)), ...
                mean(KIo),std(KIo),mean(KIIo),std(KIIo), ...
                mean(KIe),std(KIe),mean(KIIe),std(KIIe), ...
                max(T.relVectorDifference(idx)), ...
                max(T.fracNearEdgeP(idx)),max(T.fracNearEdgeQ(idx))]; %#ok<AGROW>
        end
    end

    S=array2table(rows,'VariableNames',{ ...
        'meshScale','asymmetryRequested','asymmetryUsedMedian', ...
        'A_h_median_over_radii', ...
        'KI_old_mean','KI_old_std','KII_old_mean','KII_old_std', ...
        'KI_EDI_mean','KI_EDI_std','KII_EDI_mean','KII_EDI_std', ...
        'maxRelVectorDifference','maxFracNearEdgeP','maxFracNearEdgeQ'});
end
