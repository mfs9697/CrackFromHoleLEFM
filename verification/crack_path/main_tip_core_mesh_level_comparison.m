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
%   - three fixed-physical-width vector PDFs for 0.315\linewidth subfigures;
%   - one PNG per level and a compact combined preview;
%   - unchanged mesh/COD fingerprints and metadata.
%
% The mesh is reconstructed by the audited deterministic core builder.
% No elastic solve is performed. The default MATLAB LaTeX interpreter
% matches the approved Figure 1 math labels. Euclid is optional.
% LaTeX provides panel letters/captions; only external graphical PDFs are made.

    ip=inputParser;
    addParameter(ip,'Increment_mm',4,@(x)isnumeric(x)&&isscalar(x)&&isfinite(x)&&x>0);
    addParameter(ip,'OutputDir','',@(x)ischar(x)||isstring(x));
    addParameter(ip,'XLim_mm',[-3.15 0.85],@(x)isnumeric(x)&&numel(x)==2&&all(isfinite(x))&&x(2)>x(1));
    addParameter(ip,'YLim_mm',[-1.65 1.65],@(x)isnumeric(x)&&numel(x)==2&&all(isfinite(x))&&x(2)>x(1));
    addParameter(ip,'MeshLineWidth',0.25,@(x)isnumeric(x)&&isscalar(x)&&isfinite(x)&&x>0);
    addParameter(ip,'MeshGrayLevels',[0.48 0.64 0.77], ...
        @(x)isnumeric(x)&&numel(x)==3&&all(isfinite(x))&&all(x>=0)&&all(x<=1));
    addParameter(ip,'FontMode','latex',@(x)ischar(x)||isstring(x));
    parse(ip,varargin{:});
    opt=ip.Results;
    opt.FontMode=lower(char(opt.FontMode));
    assert(ismember(opt.FontMode,{'euclid','latex'}), ...
        'meshlevels:FontMode','FontMode must be ''euclid'' or ''latex''.');

    % Default MATLAB LaTeX math matches approved Figure 1. Euclid is
    % optional; MATLAB's LaTeX interpreter cannot select an external font.
    if strcmp(opt.FontMode,'euclid')
        fonts=listfonts;
        assert(any(strcmpi(fonts,'Euclid'))&&any(strcmpi(fonts,'Euclid Symbol')), ...
            'meshlevels:EuclidFontsMissing', ...
            ['Install licensed Euclid and Euclid Symbol fonts before ' ...
             'using optional Euclid mode. The default ''FontMode'',''latex'' ' ...
             'matches the approved Figure 1 typography.']);
    end
    Pub=publication_style();
    % Exactly the final physical width of each 0.315\linewidth panel.
    panelWidth_cm=Pub.thirdWidth_cm;
    panelHeight_cm=4.50;

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

        % View the accepted topology at the same limits and spatial scale.
        % Lightening the ELEMENT EDGES does not alter mesh connectivity.
        fig=figure('Color','w','Visible','off','Units','centimeters', ...
            'Position',[2 2 panelWidth_cm panelHeight_cm]);
        ax=axes(fig);hold(ax,'on');
        local_draw_core(ax,C,opt,i);
        local_style_axes(ax,opt,Pub,true);

        pdfFile=fullfile(outDir,sprintf('tip_mesh_%s.pdf',labels(i)));
        pngFile=fullfile(outDir,sprintf('tip_mesh_%s.png',labels(i)));
        % A fixed physical PDF page prevents later LaTeX scaling from
        % silently reducing embedded lettering below the manuscript size.
        set(fig,'PaperUnits','centimeters', ...
            'PaperSize',[panelWidth_cm panelHeight_cm], ...
            'PaperPosition',[0 0 panelWidth_cm panelHeight_cm], ...
            'PaperPositionMode','manual','Renderer','painters');
        print(fig,pdfFile,'-dpdf','-painters');
        print(fig,pngFile,'-dpng','-r350');
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

    % Combined preview only; final manuscript remains three separate PDFs.
    fig=figure('Color','w','Visible','off','Units','centimeters', ...
        'Position',[2 2 Pub.textWidth_cm panelHeight_cm+0.35]);
    tl=tiledlayout(fig,1,3,'TileSpacing','compact','Padding','compact');
    for i=1:3
        ax=nexttile(tl);hold(ax,'on');
        local_draw_core(ax,cores{i},opt,i);
        local_style_axes(ax,opt,Pub,i==1);
        if i>1
            % Show one vertical scale in the preview to avoid repetition.
            ax.YTickLabel=[];
        end
    end
    previewFile=fullfile(outDir,'tip_mesh_H2_H1_H0_preview.png');
    exportgraphics(fig,previewFile,'Resolution',300);
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
    R.style=struct('fontMode',opt.FontMode,'label_pt',Pub.label_pt, ...
        'tick_pt',Pub.tick_pt,'panelWidth_cm',panelWidth_cm, ...
        'panelHeight_cm',panelHeight_cm,'meshGrayLevels',opt.MeshGrayLevels, ...
        'meshLineWidth_pt',opt.MeshLineWidth, ...
        'vectorPdf',true,'meshesPhysicallyUnchanged',true);

    save(fullfile(outDir,'tip_mesh_levels_summary.mat'),'R','-v7');
end

function local_draw_core(ax,C,opt,level)
% Render the precise audited mesh; no manipulation of nodes/connectivity.
    xy=double(1e3*C.local.coord3);
    triangles=double(C.local.connect3);
    gray=opt.MeshGrayLevels(level);
    widths=opt.MeshLineWidth*[1.15 0.90 0.70];
    % triplot(ax,...) is not supported in some MATLAB versions: ax is
    % interpreted as triangle connectivity and triangulation() fails.
    % Use an explicitly parented, edge-only patch on the EXACT topology.
    % This changes only rendering, not nodes, triangles or geometry.
    patch('Parent',ax,'Faces',triangles,'Vertices',xy, ...
        'FaceColor','none','EdgeColor',[gray gray gray], ...
        'LineWidth',widths(level));
    % Keep the crack faces conspicuously black at every resolution.
    plot(ax,[-1e3*C.rCore 0],[0 0],'k-','LineWidth',1.05);
end

function local_style_axes(ax,opt,Pub,showY)
% Axes placed identically in all three exported pages.
    axis(ax,'equal');
    xlim(ax,opt.XLim_mm);ylim(ax,opt.YLim_mm);
    box(ax,'on');grid(ax,'off');
    set(ax,'Units','normalized','Position',[.24 .25 .72 .68], ...
        'FontSize',Pub.tick_pt,'TickDir','out', ...
        'LineWidth',.75,'Layer','top');
    if strcmp(opt.FontMode,'euclid')
        set(ax,'FontName','Euclid','TickLabelInterpreter','none');
        xlabel(ax,'\fontname{Euclid}x_{1} [mm]', ...
            'Interpreter','tex','FontName','Euclid','FontSize',Pub.label_pt);
        if showY
            ylabel(ax,'\fontname{Euclid}x_{2} [mm]', ...
                'Interpreter','tex','FontName','Euclid','FontSize',Pub.label_pt);
        end
    else
        set(ax,'FontName','Times New Roman','TickLabelInterpreter','latex');
        xlabel(ax,'$x_1$ [mm]','Interpreter','latex','FontSize',Pub.label_pt);
        if showY
            ylabel(ax,'$x_2$ [mm]','Interpreter','latex','FontSize',Pub.label_pt);
        end
    end
    if ~showY
        % Keep the same viewport and y tick values; omit only repeated label.
        ax.YLabel.String='';
    end
end

function tf=local_is_absolute_path(p)
    p=char(p);
    tf=startsWith(p,filesep)|| ...
        ~isempty(regexp(p,'^[A-Za-z]:[\\/]','once'))||startsWith(p,'\\');
end
