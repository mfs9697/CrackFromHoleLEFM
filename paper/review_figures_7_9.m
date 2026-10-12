function R=review_figures_7_9(varargin)
% Pure-vector, signed-source review. No numerical history reconstruction.
ip=inputParser;addParameter(ip,'OutputRoot','',@(x)ischar(x)||isstring(x));parse(ip,varargin{:});
paper=fileparts(mfilename('fullpath'));target=char(ip.Results.OutputRoot);
if isempty(target),target=fullfile(tempdir,'CrackFromHoleLEFM_Figures_7_9_Review');end
V=review_issue78_figures();
src={'tip_resolution_sensitivity','m1_mesh_sensitivity','increment_sensitivity'};
dst={'figure7_tip_resolution','figure8_M1_exterior','figure9_increment'};
for k=1:3
    out=fullfile(target,dst{k});if ~isfolder(out),mkdir(out);end
    names={'trajectory','vertical_deviation','direction_deviation','mode_mixity'};
    if k==3,names{4}='turn_density';end
    for j=1:4
        for ext={'.pdf','.png','.layout.json'}
            name=[names{j} ext{1}];copyfile(fullfile(V.figureRoot,src{k},name),fullfile(out,name),'f');
        end
    end
end
publication_pdf_tool({'--review-root',target});
R=struct('outputRoot',target,'reviewPdf',fullfile(target,'figures_7_9_review.pdf'), ...
    'plotOnly',true,'physicalSolves',0,'figure7Source','signed accepted isolated-core vector evidence', ...
    'figure8Source','signed retained M0/M1 vector evidence; no history reconstruction', ...
    'figure9Source','signed accepted increment vector evidence; no native history reconstruction');
end
