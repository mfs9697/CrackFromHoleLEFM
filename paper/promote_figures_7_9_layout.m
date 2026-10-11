function P=promote_figures_7_9_layout(Review)
% Copy only the twelve visually reviewed panels and their layout metadata.
% Numerical data, metrics, captions and all other figures are untouched.
root=Review.outputRoot;
audit=jsondecode(fileread(fullfile(root,'layout_verification.json')));
assert(audit.passed&&numel(audit.panels)==12,'fig789:UnverifiedReview','Layout checks must pass before promotion.');
paper=fileparts(mfilename('fullpath'));
src={'figure7_tip_resolution','figure8_M1_exterior','figure9_increment'};
dst={'tip_resolution_sensitivity','m1_mesh_sensitivity','increment_sensitivity'};
files=cell(0,2);
for k=1:3
    names={'trajectory','vertical_deviation','direction_deviation','mode_mixity'};
    if k==3,names{4}='turn_density';end
    for j=1:4
        file=fullfile(root,src{k},[names{j} '.pdf']);
        row=audit.panels((k-1)*4+j);
        assert(strcmp(publication_sha256(file),row.sha256),'fig789:ChangedReview', ...
            'Review PDF changed since its geometry checks: %s',file);
        for ext={'.pdf','.png','.layout.json'}
            source=fullfile(root,src{k},[names{j} ext{1}]);
            assert(isfile(source),'fig789:IncompleteReview','Missing review output: %s',source);
            files(end+1,:)={source,fullfile(paper,'figures',dst{k},[names{j} ext{1}])}; %#ok<AGROW>
        end
    end
end
% All preflight checks complete before writing any publication graphic.
for k=1:size(files,1),copyfile(files{k,1},files{k,2},'f');end
copyfile(fullfile(root,'layout_verification.json'),fullfile(paper,'figure_7_9_layout_verification.json'),'f');
P=struct('copiedFiles',size(files,1),'figureNumbers',[7,8,9], ...
    'dataChanged',false,'otherFiguresChanged',false,'reviewRoot',root);
end
