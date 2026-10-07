function R = main_tip_core_mesh_level_comparison(varargin)
%MAIN_TIP_CORE_MESH_LEVEL_COMPARISON
% Build and plot the deterministic near-tip core at three coarse-to-fine
% levels for publication/diagnostic comparison:
%
%   H2 : CoreScale = 4, hTip = 0.1080493016 mm
%   H1 : CoreScale = 2, hTip = 0.0540246508 mm
%   H0 : CoreScale = 1, hTip = 0.0270123254 mm
%
% This is a MESH-ONLY diagnostic. It performs no full-domain qualification,
% no stiffness assembly, no physical solve, and no EDI evaluation.
%
% Important: H2 is intentionally NOT admitted by
% qualify_incremental_crack_candidate / solve_incremental_crack_tip.
% Its native COD counts in the four established windows are
% 10/15/12/9, so it does not satisfy the existing all-windows >= 12
% sampling gate. The gate is not weakened here.
%
% Output:
%   - one vector PDF per level, with no panel letter/title;
%   - one PNG per level;
%   - one combined PNG preview;
%   - table of mesh counts and native COD sampling.
%
% LaTeX should provide panel letters/captions.

    ip=inputParser;
    addParameter(ip,'Increment_mm',4,@(x)isnumeric(x)&&isscalar(x)&&isfinite(x)&&x>0);
    addParameter(ip,'OutputDir','',@(x)ischar(x)||isstring(x));
    addParameter(ip,'XLim_mm',[-3.15 0.85],@(x)isnumeric(x)&&numel(x)==2&&all(isfinite(x))&&x(2)>x(1));
    addParameter(ip,'YLim_mm',[-1.65 1.65],@(x)isnumeric(x)&&numel(x)==2&&all(isfinite(x))&&x(2)>x(1));
    addParameter(ip,'MeshLineWidth',0.25,@(x)isnumeric(x)&&isscalar(x)&&isfinite(x)&&x>0);
    parse(ip,varargin{:});
    opt=ip.Results;

    root=fileparts(fileparts(fileparts(mfilename('fullpath'))));
    addpath(genpath(root));

    outDir=char(opt.OutputDir);
    if isempty(outDir)
        outDir=fullfile(root,'paper','figures','mesh_levels');
    elseif ~local_is_absolute_path(outDir)
        outDir=fullfile(root,outDir);
    end
    if exist(outDir,'dir')~=7,mkdir(outDir);end

    a0=1e-3*opt.Increment_mm;
    scales=[4 2 1];
    labels=["H2","H1","H0"];
    expectedT3=[864 3318 12678];
    expectedNative=[10 15 12 9;19 28 23 18;38 55 44 34];
    windows=[.04 .20;.04 .30;.08 .30;.12 .30];

    rows=nan(3,8);
    cores=cell(3,1);

    fprintf('\n============================================================\n');
    fprintf('THREE-LEVEL NEAR-TIP CORE MESH COMPARISON\n');
    fprintf('============================================================\n');
    fprintf('  levels = H2 / H1 / H0 = scale 4 / 2 / 1\n');
    fprintf('  Delta a = %.9f mm\n',opt.Increment_mm);
    fprintf('  common rCore = %.9f mm\n',0.75*opt.Increment_mm);
    fprintf('  NO physical solve or EDI evaluation is performed.\n');

    for i=1:3
        C=build_stage2_scaled_audited_core([0 0],[1 0],a0,'Scale',scales(i));
        cores{i}=C;

        nT3=size(C.local.connect3,1);
        nT6=size(C.local.coord,1);
        assert(nT3==expectedT3(i),'meshlevels:CoreFingerprint', ...
            '%s T3 fingerprint changed: got %d expected %d.',labels(i),nT3,expectedT3(i));

        up=C.crack.upperT6(:);
        rr=vecnorm(C.local.coord(up,:),2,2)/a0;
        rr=unique(round(rr,14));
        native=zeros(1,4);
        for k=1:4
            native(k)=nnz(rr>=windows(k,1)-1e-13 & rr<=windows(k,2)+1e-13);
        end
        assert(isequal(native,expectedNative(i,:)), ...
            'meshlevels:NativeFingerprint', ...
            '%s native COD fingerprint changed.',labels(i));

        hTip_mm=1e3*C.hTip;
        rows(i,:)=[scales(i),hTip_mm,nT3,nT6,native];

        fprintf('  %s: scale=%g, hTip=%.10f mm, T3=%d, T6 nodes=%d, native=%s\n', ...
            labels(i),scales(i),hTip_mm,nT3,nT6,mat2str(native));

        fig=figure('Color','w','Units','centimeters','Position',[2 2 6.4 5.2]);
        ax=axes(fig); hold(ax,'on');

        P=1e3*C.local.coord3;
        T=C.local.connect3;
        triplot(T,P(:,1),P(:,2),'Color',[0.25 0.25 0.25], ...
            'LineWidth',opt.MeshLineWidth);

        % Emphasize the crack line without adding labels or panel letters.
        plot(ax,[-1e3*C.rCore 0],[0 0],'k-','LineWidth',0.75);

        axis(ax,'equal');
        xlim(ax,opt.XLim_mm);
        ylim(ax,opt.YLim_mm);
        box(ax,'on');
        grid(ax,'off');
        set(ax,'FontName','Times New Roman','FontSize',8.5, ...
            'TickDir','out','LineWidth',0.75,'Layer','top');
        xlabel(ax,'x_1 [mm]','Interpreter','tex');
        ylabel(ax,'x_2 [mm]','Interpreter','tex');

        pdfFile=fullfile(outDir,sprintf('tip_mesh_%s.pdf',labels(i)));
        pngFile=fullfile(outDir,sprintf('tip_mesh_%s.png',labels(i)));
        exportgraphics(ax,pdfFile,'ContentType','vector');
        exportgraphics(ax,pngFile,'Resolution',600);
        close(fig);
    end

    Summary=array2table(rows,'VariableNames',{ ...
        'core_scale','hTip_mm','core_T3','T6_nodes', ...
        'native_w1','native_w2','native_w3','native_w4'});

    Summary.native_sampling_gate_12= ...
        all(Summary{:,{'native_w1','native_w2','native_w3','native_w4'}}>=12,2);

    fprintf('\nMESH-LEVEL SUMMARY\n');
    disp(Summary);
    fprintf('  H2 native-sampling gate (>=12 in every window): %d\n', ...
        Summary.native_sampling_gate_12(1));
    fprintf('  H2 is mesh-visualization only unless the scientific qualification policy is changed explicitly.\n');

    % Combined preview only; final manuscript composition should use the
    % separate vector panels so LaTeX controls panel letters and caption.
    fig=figure('Color','w','Units','centimeters','Position',[2 2 19.2 5.2]);
    tl=tiledlayout(fig,1,3,'TileSpacing','compact','Padding','compact');
    for i=1:3
        ax=nexttile(tl); hold(ax,'on');
        C=cores{i};
        P=1e3*C.local.coord3; T=C.local.connect3;
        triplot(T,P(:,1),P(:,2),'Color',[0.25 0.25 0.25], ...
            'LineWidth',opt.MeshLineWidth);
        plot(ax,[-1e3*C.rCore 0],[0 0],'k-','LineWidth',0.75);
        axis(ax,'equal');
        xlim(ax,opt.XLim_mm); ylim(ax,opt.YLim_mm);
        box(ax,'on'); grid(ax,'off');
        set(ax,'FontName','Times New Roman','FontSize',8.5, ...
            'TickDir','out','LineWidth',0.75,'Layer','top');
        xlabel(ax,'x_1 [mm]','Interpreter','tex');
        if i==1
            ylabel(ax,'x_2 [mm]','Interpreter','tex');
        else
            ax.YTickLabel=[];
        end
    end
    previewFile=fullfile(outDir,'tip_mesh_H2_H1_H0_preview.png');
    exportgraphics(fig,previewFile,'Resolution',220);
    close(fig);

    writetable(Summary,fullfile(outDir,'tip_mesh_levels_summary.csv'));

    R=struct();
    R.summary=Summary;
    R.outputDir=outDir;
    R.levels=labels;
    R.scales=scales;
    R.h2PhysicalQualificationAdmitted=false;
    R.h2NativeSamplingGatePass=logical(Summary.native_sampling_gate_12(1));
    R.source='deterministic structured core only; no physical solve';

    save(fullfile(outDir,'tip_mesh_levels_summary.mat'),'R','-v7');
end

function tf=local_is_absolute_path(p)
    p=char(p);
    tf=startsWith(p,filesep)|| ...
        ~isempty(regexp(p,'^[A-Za-z]:[\\/]','once'))||startsWith(p,'\\');
end
