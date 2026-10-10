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
Report.isolated=plot_tip2h0_vs_reference_publication('OutputDir',fullfile(output,'tip_resolution_sensitivity'));close all
% Preserve known verified graphics if original raw data are absent.
m1=fullfile(root,'verification','crack_path','m1_independent_run','trajectory','path_run_state.mat');
if exist(m1,'file')==2
    Report.exterior=plot_m1_vs_reference_publication('OutputDir',fullfile(output,'m1_mesh_sensitivity'));close all
else
    retain('m1_mesh_sensitivity',{'trajectory','vertical_deviation','direction_deviation','mode_mixity'}, ...
        'M1 raw history unavailable; no numerical reconstruction; native font sizes retained');
end
base=fullfile(root,'verification','crack_path','increment_sensitivity_coarse');
if all(cellfun(@(d)exist(fullfile(base,d,'increment_study_summary.mat'),'file')==2,{'da_4mm','da_2mm','da_1mm'}))
    Report.increment=plot_increment_sensitivity_coarse_publication('OutputDir',fullfile(output,'increment_sensitivity'));close all
else
    retain('increment_sensitivity',{'trajectory'},'Full native vertex histories unavailable; retain actual trajectory graphic');
    records=fullfile(paper,'data','increment_plot_records.mat');
    if exist(records,'file')==2
        Report.increment=plot_saved_increment_records(records,fullfile(output,'increment_sensitivity'));
    else
        retain('increment_sensitivity',{'vertical_deviation','direction_deviation','turn_density'}, ...
            'Increment histories/plot records unavailable; no incomplete or two-level substitution');
    end
end
retain('mesh_levels',{'tip_mesh_H2','tip_mesh_H1','tip_mesh_H0'}, ...
    'Original coordinates/connectivity not saved; core generator deliberately not called');
fid=fopen(fullfile(output,'redraw_manifest.json'),'wb');assert(fid>=0);
fwrite(fid,unicode2native(jsonencode(Report),'UTF-8'),'uint8');fclose(fid);
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
