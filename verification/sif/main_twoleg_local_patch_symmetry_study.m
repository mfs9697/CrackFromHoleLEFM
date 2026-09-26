function Results = main_twoleg_local_patch_symmetry_study(varargin)
%MAIN_TWOLEG_LOCAL_PATCH_SYMMETRY_STUDY
% Controlled two-leg SIF test with an exactly mirror-paired local tip mesh.
%
% Workflow:
%   1. Solve the full sharp two-leg benchmark once.
%   2. Around V2, cut a circular submodel with Rpatch < delta.
%   3. Prescribe on the outer circle the displacement field interpolated
%      from the global two-leg solution.
%   4. Solve the submodel on an exactly symmetric polar crack mesh.
%   5. Perturb only upper-side interior nodes and re-solve with identical
%      physical boundary data.
%   6. Compare old mirror-J and interaction EDI.
%
% This isolates mesh asymmetry from changes in loading/geometry.
%
% Name-value options
%   'Refinement'       [Nr Ntheta] rows, default [8 48;12 72;16 96]
%   'AsymmetryLevels'  default [0 0.08 0.16 0.24]
%   'PatchOverDelta'   default 0.80
%   'rOverPatch'       default [0.30 0.45 0.60]
%   'EDIInnerFraction' default 0.20
%   'SaveOutputs'      true

    ip=inputParser;
    addParameter(ip,'Refinement',[8 48;12 72;16 96],@(x)isnumeric(x)&&size(x,2)==2);
    addParameter(ip,'AsymmetryLevels',[0 0.08 0.16 0.24],@isnumeric);
    addParameter(ip,'PatchOverDelta',0.80,@(x)isnumeric(x)&&isscalar(x)&&x>0&&x<0.95);
    addParameter(ip,'rOverPatch',[0.30 0.45 0.60],@isnumeric);
    addParameter(ip,'EDIInnerFraction',0.20,@(x)isnumeric(x)&&isscalar(x)&&x>=0&&x<1);
    addParameter(ip,'SaveOutputs',true,@(x)islogical(x)||isnumeric(x));
    addParameter(ip,'OutputDir',fullfile('verification','sif','outputs'),@(x)ischar(x)||isstring(x));
    addParameter(ip,'Verbose',true,@(x)islogical(x)||isnumeric(x));
    parse(ip,varargin{:});

    C=cfg_twoleg_historical();

    % A single global physical solution supplies identical submodel BC data
    % for every local mesh family.
    G=build_twoleg_sharp_benchmark(C,'MeshScale',1.0,'AsymmetryLevel',0,'Plot',false);

    Rpatch=ip.Results.PatchOverDelta*C.delta;
    ref=ip.Results.Refinement;
    aa=ip.Results.AsymmetryLevels(:).';
    rr=ip.Results.rOverPatch(:).';
    innerFrac=ip.Results.EDIInnerFraction;

    rows=[];
    dbg={};
    nrow=0;

    for im=1:size(ref,1)
        Nr=ref(im,1);
        Nt=ref(im,2);

        for ia=1:numel(aa)
            P=build_tip_submodel_from_global(G,Rpatch, ...
                'Nr',Nr,'Ntheta',Nt,'AsymmetryLevel',aa(ia));

            for ir=1:numel(rr)
                rI=rr(ir)*Rpatch;

                [KIold,KIIold,Dold]=SIF_LEFM_circle2_debug( ...
                    P.mesh,P.U,P.V,P.mat,rI, ...
                    'nthet',240,'eps_th',1e-3,'plot',false,'verbose',false);

                M=SIF_mesh_mirror_metrics(P.mesh,P.V,rI,'theta',Dold.theta);

                D=struct('r_inner',innerFrac*rI,'r_outer',rI);
                [KIe,KIIe,Ae]=SIF_LEFM_interaction_EDI( ...
                    P.mesh,P.U,P.V,P.mat,D, ...
                    'AuxDerivativeMode','analytic','Verbose',false);

                rel=hypot(KIold-KIe,KIIold-KIIe)/max(hypot(KIe,KIIe),eps);

                nrow=nrow+1;
                rows(nrow,:)=[ ...
                    Nr,Nt,aa(ia),P.meshInfo.asymmetryUsed, ...
                    rr(ir),rI/C.delta,rI, ...
                    M.A_h_median,M.A_h_rms,M.A_centroid_median,M.A_bary_median, ...
                    KIold,KIIold,KIe,KIIe,rel, ...
                    Dold.diagnostics.min_baryP,Dold.diagnostics.min_baryQ, ...
                    Dold.diagnostics.JI_over_absint,Dold.diagnostics.JII_over_absint, ...
                    Ae.auxConsistency.maxI,Ae.auxConsistency.maxII, ...
                    size(P.mesh.coord,1),size(P.mesh.connect,1)];

                dbg{nrow,1}=struct('old',Dold,'metric',M,'edi',Ae); %#ok<AGROW>

                if logical(ip.Results.Verbose)
                    fprintf(['patch Nr=%d Nt=%d asym=%.3f r/R=%.2f ', ...
                        'Ah=%.3e | old=(%.7e,%+.7e) EDI=(%.7e,%+.7e) ', ...
                        'diff=%.3e\n'], ...
                        Nr,Nt,P.meshInfo.asymmetryUsed,rr(ir),M.A_h_median, ...
                        KIold,KIIold,KIe,KIIe,rel);
                end
            end
        end
    end

    T=array2table(rows,'VariableNames',{ ...
        'Nr','Ntheta','asymmetryRequested','asymmetryUsed', ...
        'rI_over_Rpatch','rI_over_delta','rI', ...
        'A_h_median','A_h_rms','A_centroid_median','A_bary_median', ...
        'KI_old','KII_old','KI_EDI','KII_EDI','relVectorDifference', ...
        'minBaryP','minBaryQ','JI_over_absint','JII_over_absint', ...
        'auxMismatchMaxI','auxMismatchMaxII','nNode6','nElem6'});

    Summary=local_summary(T);

    Results=struct('C',C,'Global',G,'Rpatch',Rpatch, ...
        'Table',T,'Summary',Summary,'DebugStore',{dbg});

    if logical(ip.Results.SaveOutputs)
        outdir=char(ip.Results.OutputDir);
        if ~exist(outdir,'dir'),mkdir(outdir);end
        writetable(T,fullfile(outdir,'twoleg_local_patch_symmetry_study.csv'));
        writetable(Summary,fullfile(outdir,'twoleg_local_patch_symmetry_summary.csv'));
        save(fullfile(outdir,'twoleg_local_patch_symmetry_study.mat'),'Results','-v7.3');
    end
end


function S=local_summary(T)
    nr=unique(T.Nr);
    aa=unique(T.asymmetryRequested);
    rows=[];

    for i=1:numel(nr)
        for j=1:numel(aa)
            idx=T.Nr==nr(i) & abs(T.asymmetryRequested-aa(j))<1e-12;
            if ~any(idx),continue;end

            rows=[rows; ...
                nr(i),median(T.Ntheta(idx)),aa(j),median(T.asymmetryUsed(idx)), ...
                median(T.A_h_median(idx)),max(T.A_h_median(idx)), ...
                mean(T.KI_old(idx)),std(T.KI_old(idx)), ...
                mean(T.KII_old(idx)),std(T.KII_old(idx)), ...
                mean(T.KI_EDI(idx)),std(T.KI_EDI(idx)), ...
                mean(T.KII_EDI(idx)),std(T.KII_EDI(idx)), ...
                max(T.relVectorDifference(idx))]; %#ok<AGROW>
        end
    end

    S=array2table(rows,'VariableNames',{ ...
        'Nr','Ntheta','asymmetryRequested','asymmetryUsed', ...
        'A_h_median_over_radii','A_h_max', ...
        'KI_old_mean','KI_old_std','KII_old_mean','KII_old_std', ...
        'KI_EDI_mean','KI_EDI_std','KII_EDI_mean','KII_EDI_std', ...
        'maxRelVectorDifference'});
end
