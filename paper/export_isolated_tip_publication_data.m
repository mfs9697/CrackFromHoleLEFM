function D=export_isolated_tip_publication_data(varargin)
% Read accepted archives and export portable publication records only.
% No FEM, mesh study, EDI replay, COD refit, or new physical interpretation.
ip=inputParser;
addParameter(ip,'ArchiveDir','',@(x)ischar(x)||isstring(x));
addParameter(ip,'CorroboratingArchiveDir','',@(x)ischar(x)||isstring(x));
addParameter(ip,'OutputDir','',@(x)ischar(x)||isstring(x));parse(ip,varargin{:});o=ip.Results;
root=fileparts(fileparts(mfilename('fullpath')));addpath(genpath(root));
base=fullfile(root,'verification','crack_path');
authoritative=resolve(root,o.ArchiveDir,fullfile(base,'isolated_tip_resolution_study_20261009T200424169'));
corroborating=resolve(root,o.CorroboratingArchiveDir,fullfile(base,'isolated_tip_resolution_study_20261009T163513088'));
output=resolve(root,o.OutputDir,fullfile(root,'paper','data'));if exist(output,'dir')~=7,mkdir(output);end
z=load(fullfile(root,'paper','data','accepted_stage1_source.mat'),'R0');R0=z.R0;
z=load(fullfile(root,'paper','data','evidence_exact.mat'),'E','T');E=z.E;Ref=z.T;
[fixed,independent,inputs]=read_study(authoritative);
[fixedCheck,independentCheck,checkInputs]=read_study(corroborating);
for j=1:8
    assert(isequal(numeric_row(fixed(j)),numeric_row(fixedCheck(j))), ...
        'isolatedpub:StudyDisagreement','Fixed mechanical results differ between archives.');
end
for j=1:23
    assert(isequal(numeric_row(independent(j)),numeric_row(independentCheck(j))), ...
        'isolatedpub:StudyDisagreement','Independent mechanical results differ between archives.');
end
[~,authorTag]=fileparts(authoritative);[~,checkTag]=fileparts(corroborating);
D=struct('schemaVersion',1,'datasetId','isolated_core_tip_resolution', ...
    'authoritativeStudy',authorTag,'corroboratingStudy',checkTag, ...
    'coreScales',[2,1,.5],'exteriorScale',1,'increment_m',R0.summary.a0_reserved_m, ...
    'fixedSegments',[17,21,22,23],'mouth_m',E.vertices_m(1,:), ...
    'referencePhysicalStates',23,'independentPhysicalStates',23, ...
    'requestedExteriorLawFixed',true,'identicalExteriorConnectivityRequired',false, ...
    'fixed',fixed,'independent',independent);
dataset=fullfile(output,'isolated_tip_resolution.json');write_json(dataset,D);
fixedCsv=fullfile(output,'isolated_tip_fixed.csv');indCsv=fullfile(output,'isolated_tip_states.csv');
writetable(struct2table(rmfield(fixed,{'physicalGates','nativePoints'})),fixedCsv);
writetable(struct2table(rmfield(independent,{'physicalGates','nativePoints'})),indCsv);
normalize_lines(fixedCsv);normalize_lines(indCsv);
manifest=struct('schemaVersion',1,'datasetId',D.datasetId,'authoritativeStudy',authorTag, ...
    'corroboratingStudy',checkTag,'mechanicalAgreement','bitwise identical checked saved mechanical records', ...
    'validatedPhysicalRecordsPerStudy',31,'physicalSolvesPerformed',0, ...
    'sourceCommit',strtrim(git_commit(root)),'sourceInputs',[inputs,checkInputs], ...
    'referenceInputs',[entry(fullfile(root,'paper','data','evidence_exact.mat')), ...
        entry(fullfile(root,'paper','data','evidence.json')),entry(fullfile(root,'paper','data','accepted_stage1_source.mat'))], ...
    'portableFiles',[entry(dataset),entry(fixedCsv),entry(indCsv)]);
write_json(fullfile(output,'isolated_tip_resolution_manifest.json'),manifest);
fprintf('Exported validated isolated-core publication data: 8 fixed + 23 independent; no physical solves.\n');

    function [F,I,sourceInputs]=read_study(dir)
        z=load(fullfile(dir,'isolated_tip_resolution_study.mat'),'D');study=z.D;
        assert(study.exteriorScale==1&&isequal(study.coreScales,[2,.5])&& ...
            height(study.fixedPhysical)==8&&height(study.independentPath.stepTable)==23&& ...
            strcmp(study.independentPath.stopReason,'max_segments_reached'), ...
            'isolatedpub:IncompleteStudy','Wrong or incomplete study.');
        trajectory=fullfile(dir,'independent_core_2_exterior_1','trajectory');
        z=load(fullfile(trajectory,'path_run_state.mat'),'State');S=z.State;
        assert(S.completedPhysicalSegments==23&&S.exteriorMeshControls.coreScale==2&& ...
            S.exteriorMeshControls.exteriorScale==1,'isolatedpub:Family','Wrong path family.');
        F=struct([]);I=struct([]);sourceInputs=entry(fullfile(dir,'isolated_tip_resolution_study.mat'));
        sourceInputs(end+1)=entry(fullfile(trajectory,'path_run_state.mat'));
        for scale=[2,.5]
            if scale==2,label='core_2_exterior_1';else,label='core_0p5_exterior_1';end
            for k=[17,21,22,23]
                fd=fullfile(dir,'fixed_geometry',label);
                cf=fullfile(fd,sprintf('step_%03d_candidate.mat',k));rf=fullfile(fd,sprintf('step_%03d_physical_small.mat',k));
                cp=fullfile(fd,sprintf('step_%03d_physical_solved.mat',k));
                z=load(cf,'candidate');c=z.candidate;z=load(rf,'R');R=z.R;
                assert(norm(c.path-E.vertices_m(1:k+1,:),'fro')<=2e-12,'isolatedpub:FixedGeometry','Wrong fixed reference geometry.');
                validate_isolated_fixed_result_cache(R,c,R0,cp);
                validate_candidate(c,scale);
                row=record(R,c,k,R.thetaCurrentDeg,R.deltaThetaNextDeg,R.thetaNextDeg);
                assert(row.pass,'isolatedpub:Acceptance','Unaccepted fixed state.');if isempty(F),F=row;else,F(end+1)=row;end
                sourceInputs=[sourceInputs,entry(cf),entry(rf),entry(cp)]; %#ok<AGROW>
            end
        end
        for k=1:23
            if k==1
                fd=fullfile(dir,'independent_core_2_exterior_1','p1_seed');
                cf=fullfile(fd,'P1_candidate.mat');rf=fullfile(fd,'P1_physical_small.mat');cp=fullfile(fd,'P1_physical_solved.mat');
                z=load(cf,'candidate');c=z.candidate;z=load(rf,'R');R=z.R;
                validate_increment_study_p1_cache(R,c,R0,cp);
                assert(isequal([R.EDI.KI_unit,R.EDI.KII_unit],S.rowsThroughCompleted(1,5:6)), ...
                    'isolatedpub:Seed','State is not seeded by its own accepted P1 field.');
                row=record(R,c,k,0,S.rowsThroughCompleted(1,8),S.rowsThroughCompleted(1,9));
            else
                cf=fullfile(trajectory,sprintf('step_%03d_candidate.mat',k));
                rf=fullfile(trajectory,sprintf('step_%03d_physical_small.mat',k));cp=fullfile(trajectory,sprintf('step_%03d_physical_solved.mat',k));
                z=load(cf,'candidate');c=z.candidate;z=load(rf,'R');R=z.R;
                validate_isolated_fixed_result_cache(R,c,R0,cp);
                assert(norm(R.pathFixed-S.vertices(1:k+1,:),'fro')<=2e-12,'isolatedpub:Path','Path/physical record disagreement.');
                row=record(R,c,k,R.thetaCurrentDeg,R.deltaThetaNextDeg,R.thetaNextDeg);
            end
            validate_candidate(c,2);
            assert(norm([row.tip_x_m,row.tip_y_m]-S.vertices(k+1,:))<=2e-12,'isolatedpub:Tip','Saved tip differs.');
            if isempty(I),I=row;else,I(end+1)=row;end
            sourceInputs=[sourceInputs,entry(cf),entry(rf),entry(cp)]; %#ok<AGROW>
        end
    end

    function validate_candidate(c,scale)
        controls=struct('coreScale',scale,'exteriorScale',1,'farCapOverIncrement',.625, ...
            'transitionOverIncrement',1,'calibrationOverride',struct());
        assert(crack_candidate_mesh_controls_match(c,controls)&& ...
            all(structfun(@logical,c.gates))&&all(structfun(@logical,c.syntheticGates)), ...
            'isolatedpub:Qualification','Wrong scale/exterior or failed qualified candidate.');
    end

    function row=record(R,c,k,theta,turn,next)
        assert(R.pass&&all(structfun(@logical,R.gates)),'isolatedpub:Acceptance','Physical gates failed.');
        tip=c.crack.Pmid(end,:);
        row=struct('segment',k,'core_scale',c.coreMeshControls.scale,'exterior_scale',1, ...
            'crack_length_mm',1e3*R0.summary.a0_reserved_m*k,'tip_x_m',tip(1),'tip_y_m',tip(2), ...
            'KI_unit',R.EDI.KI_unit,'KII_unit',R.EDI.KII_unit,'KII_over_KI',R.EDI.KII_over_KI, ...
            'theta_deg',theta,'delta_theta_next_deg',turn,'theta_next_deg',next, ...
            'hTip_m',c.coreMeshControls.hTip_m,'T3_elements',R.summary.T3_elements, ...
            'EDI_elements',R.EDI.EDI_elements,'PCG_iterations',R.solverInfo.iter, ...
            'true_rel_residual',R.solverInfo.trueRelResidual,'pcg_flag',R.solverInfo.flag, ...
            'pass',logical(R.pass),'qualificationPassed',all(structfun(@logical,c.gates)), ...
            'syntheticPassed',all(structfun(@logical,c.syntheticGates)), ...
            'nativePoints',c.coreMeshControls.expectedNativeSamples(:).','physicalGates',R.gates);
    end

    function item=entry(file)
        path=strrep(file,[root filesep],'');path=strrep(path,'\','/');
        item=struct('path',path,'sha256',publication_sha256(file));
    end
end

function row=numeric_row(s)
names={'segment','core_scale','exterior_scale','tip_x_m','tip_y_m','KI_unit','KII_unit','KII_over_KI', ...
    'theta_deg','delta_theta_next_deg','theta_next_deg','T3_elements','EDI_elements'};
row=cellfun(@(n)s.(n),names);
end
function write_json(file,data)
fid=fopen(file,'wb');assert(fid>=0);cleanup=onCleanup(@()fclose(fid));
fwrite(fid,unicode2native([jsonencode(data),char(10)],'UTF-8'),'uint8');
end
function normalize_lines(file)
text=strrep(fileread(file),sprintf('\r\n'),sprintf('\n'));
fid=fopen(file,'wb');assert(fid>=0);cleanup=onCleanup(@()fclose(fid));
fwrite(fid,unicode2native(text,'UTF-8'),'uint8');
end
function path=resolve(root,given,default)
path=char(given);if isempty(path),path=default;
elseif isempty(regexp(path,'^([A-Za-z]:[\\/]|[/\\])','once')),path=fullfile(root,path);end
end
function commit=git_commit(root)
[status,commit]=system(sprintf('git -C "%s" rev-parse HEAD',root));assert(status==0);
end
