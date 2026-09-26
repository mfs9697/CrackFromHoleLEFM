function Results = main_twoleg_edi_domain_sensitivity(varargin)
%MAIN_TWOLEG_EDI_DOMAIN_SENSITIVITY  Sweep EDI inner/outer radii.
%
% The outer radius is normalized by the second-leg length delta because an
% integration domain around V2 must not engulf the knee V1.

    ip=inputParser;
    addParameter(ip,'MeshScaleList',[1.20 1.00 0.80],@(x)isnumeric(x)&&all(x>0));
    addParameter(ip,'AsymmetryLevels',[0 0.20],@isnumeric);
    addParameter(ip,'OuterOverDelta',[0.30 0.40 0.50 0.60],@isnumeric);
    addParameter(ip,'InnerFractions',[0.00 0.20 0.35],@isnumeric);
    addParameter(ip,'SaveOutputs',true,@(x)islogical(x)||isnumeric(x));
    addParameter(ip,'OutputDir',fullfile('verification','sif','outputs'),@(x)ischar(x)||isstring(x));
    parse(ip,varargin{:});

    C=cfg_twoleg_historical();
    rows=[];

    for ms=ip.Results.MeshScaleList(:).'
        for aa=ip.Results.AsymmetryLevels(:).'
            S=build_twoleg_sharp_benchmark(C,'MeshScale',ms,'AsymmetryLevel',aa,'Plot',false);

            for ro=ip.Results.OuterOverDelta(:).'
                rout=ro*C.delta;

                MM=SIF_mesh_mirror_metrics(S.mesh,S.V,rout,'nthet',C.nthet, ...
                    'eps_th',C.eps_th);

                for fi=ip.Results.InnerFractions(:).'
                    rin=fi*rout;
                    if rin>=rout, continue; end

                    D=struct('r_inner',rin,'r_outer',rout);
                    [KI,KII,A]=SIF_LEFM_interaction_EDI(S.mesh,S.U,S.V,S.mat,D, ...
                        'AuxDerivativeMode','analytic','Verbose',false);

                    rows=[rows; ...
                        ms,aa,S.meshInfo.asymmetryLevelUsed, ...
                        ro,fi,rin,rout,MM.A_h_median, ...
                        KI,KII,KII/KI,A.I_modeI,A.I_modeII, ...
                        A.nGP_used,A.nElem_used, ...
                        A.auxConsistency.maxI,A.auxConsistency.maxII]; %#ok<AGROW>
                end
            end
        end
    end

    T=array2table(rows,'VariableNames',{ ...
        'meshScale','asymmetryRequested','asymmetryUsed', ...
        'rOuter_over_delta','rInner_over_rOuter','rInner','rOuter', ...
        'A_h_median','KI_EDI','KII_EDI','KII_over_KI', ...
        'I_modeI','I_modeII','nGPused','nElemUsed', ...
        'auxMismatchMaxI','auxMismatchMaxII'});

    Summary=local_summary(T);

    Results=struct('C',C,'Table',T,'Summary',Summary);

    if logical(ip.Results.SaveOutputs)
        outdir=char(ip.Results.OutputDir);
        if ~exist(outdir,'dir'),mkdir(outdir);end
        writetable(T,fullfile(outdir,'twoleg_edi_domain_sensitivity.csv'));
        writetable(Summary,fullfile(outdir,'twoleg_edi_domain_summary.csv'));
        save(fullfile(outdir,'twoleg_edi_domain_sensitivity.mat'),'Results','-v7.3');
    end
end


function S=local_summary(T)
    ms=unique(T.meshScale);
    aa=unique(T.asymmetryRequested);
    fi=unique(T.rInner_over_rOuter);
    rows=[];

    for i=1:numel(ms)
        for j=1:numel(aa)
            for k=1:numel(fi)
                idx=abs(T.meshScale-ms(i))<1e-12 & ...
                    abs(T.asymmetryRequested-aa(j))<1e-12 & ...
                    abs(T.rInner_over_rOuter-fi(k))<1e-12;
                if ~any(idx),continue;end
                ki=T.KI_EDI(idx);
                kii=T.KII_EDI(idx);
                rows=[rows;ms(i),aa(j),fi(k), ...
                    mean(ki),std(ki),min(ki),max(ki), ...
                    mean(kii),std(kii),min(kii),max(kii)]; %#ok<AGROW>
            end
        end
    end

    S=array2table(rows,'VariableNames',{ ...
        'meshScale','asymmetryRequested','rInner_over_rOuter', ...
        'KImean','KIstd','KImin','KImax', ...
        'KIImean','KIIstd','KIImin','KIImax'});
end
