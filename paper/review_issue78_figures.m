function R=review_issue78_figures(varargin)
% Immutable vector drawings -> internal legends -> common tight crops.
% No history reconstruction, FEM/meshing/EDI/COD/MTS computation.
ip=inputParser;addParameter(ip,'OutputRoot','',@(x)ischar(x)||isstring(x));parse(ip,varargin{:});
paper=fileparts(mfilename('fullpath'));root=char(ip.Results.OutputRoot);
if isempty(root),root=fullfile(paper,'reviews','issue78');end
source=fullfile(paper,'reviews','issue78');
if ~strcmp(root,source)
    if ~isfolder(root),mkdir(root);end
    copyfile(fullfile(source,'baseline_manifest.json'),fullfile(root,'baseline_manifest.json'),'f');
    copyfile(fullfile(source,'before'),fullfile(root,'before'),'f');
end
publication_m0_label_stamps(fullfile(root,'stamps'));
publication_pdf_tool({'--asset-review-root',root});
R=struct('outputRoot',root,'figureRoot',fullfile(root,'after','figures'), ...
    'geometryReport',fullfile(root,'geometry_report.json'),'plotOnly',true, ...
    'physicalSolves',0,'numericalHistoriesRead',false,'promotion','requires author approval under Issue #78');
end
