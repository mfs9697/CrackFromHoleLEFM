function Review=review_figures_7_9(varargin)
%REVIEW_FIGURES_7_9 Prepare the common four-panel Figures 7--9 for review.
%
% All three figures use one wide trajectory panel above three companion
% panels, as already defined in paper/main.tex. This helper does NOT
% change published data, solve elasticity, replay EDI, fit COD, or
% fabricate unavailable histories.
%
% Figure 7: redraw four panels from the audited isolated-tip dataset.
% Figure 8: WRAP the four immutable M1 vector PDFs (history unavailable).
% Figure 9: RETAIN the exact native-trajectory PDF; regenerate only three
%           companion panels from the saved 4/2/1-mm plotting records.
%
% Review outputs default to a folder OUTSIDE the Git working tree, to
% avoid GitHub Desktop repeatedly flagging generated review files:
%
%   addpath(genpath(pwd));
%   R=review_figures_7_9();
%   winopen(R.outputRoot);                  % Windows
%
% To choose a destination:
%   R=review_figures_7_9('OutputRoot','D:\FigureReview');
%
% Also produces a three-page assembled review PDF and a geometry manifest.

    ip=inputParser;
    addParameter(ip,'OutputRoot','',@(x)ischar(x)||isstring(x));
    parse(ip,varargin{:});

    paperDir=fileparts(mfilename('fullpath'));
    root=fileparts(paperDir);
    addpath(genpath(root));

    target=char(ip.Results.OutputRoot);
    if isempty(target)
        target=fullfile(tempdir,'CrackFromHoleLEFM_Figures_7_9_Review');
    elseif ~local_absolute_path(target)
        target=fullfile(root,target);
    end
    if ~isfolder(target),mkdir(target);end

    dir7=fullfile(target,'figure7_tip_resolution');
    dir8=fullfile(target,'figure8_M1_exterior');
    dir9=fullfile(target,'figure9_increment');

    names7={'trajectory','vertical_deviation', ...
            'direction_deviation','mode_mixity'};
    names8=names7;
    names9={'trajectory','vertical_deviation', ...
            'direction_deviation','turn_density'};

    % Preflight all RETAINED numerical evidence before drawing anything.
    incData=fullfile(paperDir,'data','increment_plot_records.mat');
    assert(isfile(incData),'fig789:MissingSavedIncrementRecords', ...
        'The audited 4/2/1-mm plotting records are missing: %s',incData);
    retained8=fullfile(paperDir,'data','retained_figure_7_9');
    for k=1:numel(names8)
        local_require_pdf(fullfile(retained8,['m1_' names8{k} '.pdf']));
    end
    local_require_pdf(fullfile(retained8,'increment_trajectory.pdf'));

    % Figure 7: audited independent core-size sensitivity.
    if ~isfolder(dir7),mkdir(dir7);end
    result7=plot_tip2h0_vs_reference_publication('OutputDir',dir7);
    for k=1:numel(names7)
        local_require_pdf(fullfile(dir7,[names7{k} '.pdf']));
    end

    % Figure 8: preserve accepted vector paths, reframe axes and relocate the
    % original legend glyphs into the reserved strip. No numerical history.
    if ~isfolder(dir8),mkdir(dir8);end
    for k=1:numel(names8)
        publication_normalize_retained_pdf(['m1_' names8{k}], ...
            fullfile(dir8,[names8{k} '.pdf']));
    end

    % Figure 9: use ONLY the saved exact common-length and native records.
    if ~isfolder(dir9),mkdir(dir9);end
    publication_normalize_retained_pdf('increment_trajectory',fullfile(dir9,'trajectory.pdf'));
    result9=plot_saved_increment_records(incData,dir9);
    for k=1:numel(names9)
        local_require_pdf(fullfile(dir9,[names9{k} '.pdf']));
    end

    Review=struct();
    Review.outputRoot=target;
    Review.figure7=struct('directory',dir7, ...
        'panels',{names7},'status','redrawn from audited isolated-core evidence');
    Review.figure8=struct('directory',dir8, ...
        'panels',{names8},'status','wrapped original M1 vector content; history unavailable');
    Review.figure9=struct('directory',dir9, ...
        'panels',{names9}, ...
        'status','original trajectory retained; three saved-data panels redrawn');
    Review.plotOnly=true;
    Review.physicalSolves=0;
    Review.CODRefits=0;
    Review.EDIReplays=0;
    Review.latexStructure='one 0.94 linewidth top panel; three 0.315 linewidth lower panels';
    Review.figure7Source=result7;
    Review.figure9Source=result9;
    Review.layout=jsondecode(fileread(fullfile(paperDir,'figure_7_9_layout.json')));
    fid=fopen(fullfile(target,'review_manifest.json'),'w');assert(fid>=0);
    fprintf(fid,'%s\n',jsonencode(Review));fclose(fid);
    publication_pdf_tool({'--review-root',target});
    Review.reviewPdf=fullfile(target,'figures_7_9_review.pdf');

    fprintf('\nFIGURES 7--9 READY FOR GRAPHICAL REVIEW\n');
    fprintf('  %s\n',target);
    fprintf('  Fig. 7: four panels from audited isolated-core data\n');
    fprintf('  Fig. 8: four VERIFIED vector layers wrapped; no history reconstruction\n');
    fprintf('  Fig. 9: original trajectory plus three saved-record panels\n');
    fprintf('  All 12 PDFs checked for existence and PDF file signatures.\n');
    fprintf('  Assembled preview and physical geometry checks: %s\n',Review.reviewPdf);
end

function local_require_pdf(file)
    assert(isfile(file),'fig789:MissingPanel','Missing panel: %s',file);
    d=dir(file);
    assert(d.bytes>1000,'fig789:InvalidPanel','PDF is too small: %s',file);
    f=fopen(file,'rb');
    assert(f>=0,'fig789:UnreadablePanel','Cannot open: %s',file);
    header=fread(f,5,'*uint8')';
    fclose(f);
    assert(isequal(header,uint8('%PDF-')), ...
        'fig789:BadPDF','Invalid PDF signature: %s',file);
end

function tf=local_absolute_path(path)
    tf=startsWith(path,filesep) || ...
        ~isempty(regexp(path,'^[A-Za-z]:[\\/]','once')) || ...
        startsWith(path,'\\');
end
