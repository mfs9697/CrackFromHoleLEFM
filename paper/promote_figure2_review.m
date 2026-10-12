function Report = promote_figure2_review(varargin)
%PROMOTE_FIGURE2_REVIEW Install locally approved Figure 2 panels.
%
% From the repository root, after pulling paper-figure2-mesh-typography:
%   addpath(genpath(pwd));
%   Report = promote_figure2_review();
%
% Copies the three APPROVED review PDFs and PNGs into the manuscript's
% canonical figures/mesh_levels directory. This operation does not build
% a mesh, solve elasticity, or alter SIF/COD calculations.
% Source data are checked against the established exact mesh fingerprints.
%
% To preview files before promotion, inspect:
%   paper/figures/mesh_levels_review/tip_mesh_H2_H1_H0_preview.png
%
% The publication PDFs remain ordinary Git files. Commit the six copied
% files with GitHub Desktop after local validation to record their version.

    ip=inputParser;
    addParameter(ip,'ReviewDir','',@(x)ischar(x)||isstring(x));
    parse(ip,varargin{:});

    paperDir=fileparts(mfilename('fullpath'));
    sourceDir=char(ip.Results.ReviewDir);
    if isempty(sourceDir)
        sourceDir=fullfile(paperDir,'figures','mesh_levels_review');
    end
    destinationDir=fullfile(paperDir,'figures','mesh_levels');
    assert(isfolder(sourceDir),'fig2:MissingReviewDir', ...
        'Cannot find the reviewed mesh panels: %s',sourceDir);
    assert(~strcmpi(sourceDir,destinationDir),'fig2:SameDirectory', ...
        'Review and publication folders must differ.');

    summaryPath=fullfile(sourceDir,'tip_mesh_levels_summary.mat');
    assert(isfile(summaryPath),'fig2:MissingSummary', ...
        'Mesh fingerprint summary missing: %s',summaryPath);
    saved=load(summaryPath,'R');
    assert(isfield(saved,'R')&&isfield(saved.R,'summary'), ...
        'fig2:InvalidSummary','Figure 2 review summary has unexpected contents.');
    T=saved.R.summary;

    expectedScale=[4;2;1];
    expectedT3=[864;3318;12678];
    expectedT6=[1807;6787;25649];
    expectedNative=[10 15 12 9;19 28 23 18;38 55 44 34];
    expectedTip=[0.1080493016;0.0540246508;0.0270123254];

    assert(height(T)==3 ...
        && isequal(double(T.core_scale(:)),expectedScale) ...
        && isequal(double(T.core_T3(:)),expectedT3) ...
        && isequal(double(T.T6_nodes(:)),expectedT6), ...
        'fig2:MeshFingerprint','Figure 2 element/node fingerprints changed.');
    native=double([T.native_w1,T.native_w2,T.native_w3,T.native_w4]);
    assert(isequal(native,expectedNative), ...
        'fig2:CODFingerprint','Figure 2 COD sampling fingerprints changed.');
    assert(all(abs(double(T.hTip_mm(:))-expectedTip)<1e-7), ...
        'fig2:TipFingerprint','Figure 2 tip sizes changed.');
    assert(isequal(all(native>=12,2),[false;true;true]), ...
        'fig2:SamplingGate','H2 must remain illustration only; H1/H0 qualified.');

    names={'tip_mesh_H2.pdf','tip_mesh_H1.pdf','tip_mesh_H0.pdf', ...
           'tip_mesh_H2.png','tip_mesh_H1.png','tip_mesh_H0.png'};
    % Preflight EVERY input before overwriting any canonical graphics.
    for i=1:numel(names)
        file=fullfile(sourceDir,names{i});
        assert(isfile(file),'fig2:MissingGraphic', ...
            'Expected reviewed graphics file is missing: %s',file);
        d=dir(file);
        assert(d.bytes>1000,'fig2:InvalidGraphic', ...
            'Reviewed graphics file is unexpectedly small: %s',file);
        fid=fopen(file,'rb');
        assert(fid~=-1,'fig2:UnreadableGraphic','Cannot read %s',file);
        header=fread(fid,8,'*uint8')';
        fclose(fid);
        if endsWith(names{i},'.pdf')
            assert(numel(header)>=5&&isequal(header(1:5),uint8('%PDF-')), ...
                'fig2:InvalidPDF','File is not a PDF: %s',file);
        else
            assert(numel(header)==8 ...
                && isequal(header,uint8([137 80 78 71 13 10 26 10])), ...
                'fig2:InvalidPNG','File is not a PNG: %s',file);
        end
    end

    if ~isfolder(destinationDir),mkdir(destinationDir);end
    for i=1:numel(names)
        src=fullfile(sourceDir,names{i});
        dst=fullfile(destinationDir,names{i});
        [ok,msg]=copyfile(src,dst,'f');
        assert(ok,'fig2:CopyFailed','Cannot install %s: %s',names{i},msg);
        a=dir(src);b=dir(dst);
        assert(a.bytes==b.bytes,'fig2:IncompleteCopy', ...
            'File size differs after copying %s.',names{i});
    end

    fprintf('\nFIGURE 2 APPROVED FILES INSTALLED\n');
    fprintf('  Source:      %s\n',sourceDir);
    fprintf('  Publication: %s\n',destinationDir);
    fprintf('  H2: geometric illustration only; native gate FAIL (expected)\n');
    fprintf('  H1/H0: native sampling gates PASS\n');
    fprintf('  Commit the six updated PDFs/PNGs in GitHub Desktop when ready.\n');
    Report=struct('sourceDir',sourceDir,'publicationDir',destinationDir, ...
        'files',{names},'summary',T,'physicalSolvePerformed',false);
end
