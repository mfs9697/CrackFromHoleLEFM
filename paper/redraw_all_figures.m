function Report=redraw_all_figures(varargin)
% Complete plot-only workflow. Never invokes a solver, mesher, EDI or COD fit.
ip=inputParser;addParameter(ip,'OutputRoot','',@(x)ischar(x)||isstring(x));parse(ip,varargin{:});
paper=fileparts(mfilename('fullpath'));root=fileparts(paper);addpath(genpath(root));
output=char(ip.Results.OutputRoot);if isempty(output),output=fullfile(paper,'figures');end
if exist(output,'dir')~=7,mkdir(output);end
Report=struct('plotOnly',true,'physicalSolves',0,'meshStudies',0,'EDIReplays',0,'CODRefits',0, ...
    'style',publication_style(),'arrangement','one wide trajectory and three companion panels', ...
    'retainedGraphics',struct([]));
Report.specimen=plot_specimen_geometry('OutputDir',fullfile(output,'geometry_loading'));
Report.reference=plot_existing_manuscript_figures('OutputRoot',output);close all
% One layout path for all sensitivity figures. Immutable original PDFs are
% wrapped when histories are unavailable; no fallback to historical data.
R=review_figures_7_9();close all
Report.sensitivity=R;
Report.isolated=R.figure7Source;Report.exterior=R.figure8;Report.increment=R.figure9Source;
sources={'figure7_tip_resolution','figure8_M1_exterior','figure9_increment'};
targets={'tip_resolution_sensitivity','m1_mesh_sensitivity','increment_sensitivity'};
for k=1:3
    dest=fullfile(output,targets{k});if ~isfolder(dest),mkdir(dest);end
    names={'trajectory','vertical_deviation','direction_deviation','mode_mixity'};
    if k==3,names{4}='turn_density';end
    for j=1:numel(names)
        for ext={'.pdf','.png','.layout.json'}
            name=[names{j} ext{1}];
            copyfile(fullfile(R.outputRoot,sources{k},name),fullfile(dest,name),'f');
        end
    end
    for j=1:numel(names)
        info=jsondecode(fileread(fullfile(dest,[names{j} '.layout.json'])));
        if isfield(info,'source_sha256')
            row=struct('file',fullfile(dest,[names{j} '.pdf']), ...
                'source',info.source,'sha256',info.source_sha256, ...
                'status','retained vector content, layout normalized', ...
                'reason','Exact numerical history unavailable; original plotting commands preserved', ...
                'fontTargetVerified',false);
            if isempty(Report.retainedGraphics),Report.retainedGraphics=row;else,Report.retainedGraphics(end+1)=row;end
        end
    end
end
retain('mesh_levels',{'tip_mesh_H2','tip_mesh_H1','tip_mesh_H0'}, ...
    'Original coordinates/connectivity not saved; core generator deliberately not called');
fid=fopen(fullfile(output,'redraw_manifest.json'),'wb');assert(fid>=0);
portable=local_portable_manifest(Report,root,output);
fwrite(fid,unicode2native(jsonencode(portable),'UTF-8'),'uint8');fclose(fid);
fprintf('Plot-only workflow complete. Retained-source limitations are in redraw_manifest.json.\n');

    function retain(folder,names,reason)
        dest=fullfile(output,folder);if exist(dest,'dir')~=7,mkdir(dest);end
        for j=1:numel(names)
            input=fullfile(paper,'figures',folder,[names{j} '.pdf']);
            assert(exist(input,'file')==2,'paperfig:MissingRetainedGraphic','Missing verified graphic: %s',input);
            target=fullfile(dest,[names{j} '.pdf']);
            if ~strcmp(input,target),copyfile(input,target);end
            row=struct('file',target,'source',input,'sha256',publication_sha256(input), ...
                'status','retained, not regenerated','reason',reason,'fontTargetVerified',false);
            if isempty(Report.retainedGraphics),Report.retainedGraphics=row;else,Report.retainedGraphics(end+1)=row;end
        end
        fprintf('Retained %s: %s\n',folder,reason);
    end
end

function value=local_portable_manifest(value,repositoryRoot,outputRoot)
% Serialize stable, machine-independent paths without changing the caller's
% returned MATLAB Report struct or the original source provenance.
    if isstruct(value)
        for j=1:numel(value)
            keys=fieldnames(value(j));
            for k=1:numel(keys)
                key=keys{k};
                value(j).(key)=local_portable_manifest( ...
                    value(j).(key),repositoryRoot,outputRoot);
            end
        end
    elseif iscell(value)
        for j=1:numel(value)
            value{j}=local_portable_manifest(value{j},repositoryRoot,outputRoot);
        end
    elseif isstring(value)
        for j=1:numel(value)
            value(j)=string(local_portable_manifest(char(value(j)),repositoryRoot,outputRoot));
        end
    elseif ischar(value) && isrow(value)
        if strcmpi(value,repositoryRoot)
            value='repo:.';
        elseif startsWith(value,[repositoryRoot filesep],'IgnoreCase',true)
            value=['repo:' strrep(value(numel(repositoryRoot)+2:end),filesep,'/')];
        elseif strcmpi(value,outputRoot)
            value='output:.';
        elseif startsWith(value,[outputRoot filesep],'IgnoreCase',true)
            value=['output:' strrep(value(numel(outputRoot)+2:end),filesep,'/')];
        end
    end
end
