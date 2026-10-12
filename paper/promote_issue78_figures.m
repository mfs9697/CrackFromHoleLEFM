function P=promote_issue78_figures(Review,varargin)
% Explicit author action only; draft review generation never promotes assets.
ip=inputParser;addParameter(ip,'AuthorApproved',false,@(x)islogical(x)&&isscalar(x));parse(ip,varargin{:});
assert(ip.Results.AuthorApproved,'publication:AuthorApprovalRequired', ...
    'Issue #78: "author must inspect the new publication PDFs and approve before promotion/merge". https://github.com/mfs9697/CrackFromHoleLEFM/issues/78');
paper=fileparts(mfilename('fullpath'));root=Review.outputRoot;
verified=jsondecode(fileread(fullfile(root,'verification.json')));assert(verified.passed&&verified.panelCount==31);
report=jsondecode(fileread(fullfile(root,'geometry_report.json')));files=cell(0,2);
for k=1:numel(report.panels)
    row=report.panels(k);source=fullfile(root,'after',row.asset);
    assert(strcmp(publication_sha256(source),row.output_sha256),'publication:ChangedCandidate','Changed reviewed PDF');
    for ext={'.pdf','.png','.layout.json'}
        relative=strrep(row.asset,'.pdf',ext{1});input=fullfile(root,'after',relative);
        assert(isfile(input),'publication:IncompleteCandidate','Missing %s',input);
        files(end+1,:)={input,fullfile(paper,relative)}; %#ok<AGROW>
    end
end
for k=1:size(files,1),copyfile(files{k,1},files{k,2},'f');end
P=struct('assetCount',31,'authorApproved',true,'scientificDataChanged',false);
end
