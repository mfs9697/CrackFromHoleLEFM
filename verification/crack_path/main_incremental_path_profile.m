function Report=main_incremental_path_profile(varargin)
%MAIN_INCREMENTAL_PATH_PROFILE Phase A only: sequential strict benchmarks.
% Empty fresh directories are required. Existing investigator results are
% read as full-precision anchors and never overwritten. No optimization,
% changed gate, angle sweep, parallelism or new Stage-I solve is performed.
ip=inputParser;
addParameter(ip,'FrozenState',[],@(x)isstruct(x)&&isscalar(x));
addParameter(ip,'OutputDir','',@(x)ischar(x)||isstring(x));
addParameter(ip,'ReferenceDir','',@(x)ischar(x)||isstring(x));
addParameter(ip,'AllowPhysicalSolves',false,@(x)islogical(x)&&isscalar(x));
parse(ip,varargin{:});opt=ip.Results;
assert(~isempty(opt.FrozenState),'pathprofile:FrozenState','Pass the exact accepted R0.');
assert(opt.AllowPhysicalSolves,'pathprofile:SolveGuard','Fresh benchmarks require explicit authorization.');
root=fileparts(fileparts(fileparts(mfilename('fullpath'))));addpath(genpath(root));
out=char(opt.OutputDir);assert(~isempty(out),'pathprofile:OutputDir','Use an isolated benchmark directory.');
ref=char(opt.ReferenceDir);
if isempty(ref),ref=fullfile(root,'verification','crack_path','incremental_run');end
old=load(fullfile(ref,'incremental_path_result.mat'),'Path');
assert(old.Path.nSegments==5&&all(old.Path.stepTable.pass),'pathprofile:Reference','Accepted five-step reference required.');
for name={'fresh3','fresh5'}
    dirName=fullfile(out,name{1});
    assert(exist(dirName,'dir')~=7,'pathprofile:FreshDirectory','Fresh benchmark folder already exists: %s',dirName);
end
if exist(out,'dir')~=7,mkdir(out);end
cleanup=onCleanup(@()incremental_profile_clock('reset',false));
Report=struct('matlabVersion',version,'baselineCommit',git_head(root), ...
    'scientificBaseCommit','10c39ac1e38a76520817d44f4a2d31ef28cb25e2', ...
    'runOrder',{{'fresh3','fresh5','resumed5'}},'optimizationPerformed',false);
names={'fresh3','fresh5','resumed5'};
for j=1:3
    name=names{j};incremental_profile_clock('reset',true);wall=tic;
    if j==1
        Path=main_incremental_path_regression('FrozenState',opt.FrozenState, ...
            'AllowSolve',true,'OutputDir',fullfile(out,'fresh3'));
    else
        Path=run_incremental_crack_path('FrozenState',opt.FrozenState, ...
            'MaxSegments',5,'AllowPhysicalSolves',j==2,'RunSynthetic',true, ...
            'OutputDir',fullfile(out,'fresh5'));
    end
    wallSeconds=toc(wall);phases=incremental_profile_clock('snapshot');
    actual=Path.stepTable;reference=old.Path.stepTable(1:height(actual),:);
    fields={'KI_unit','KII_unit','delta_theta_next_deg','theta_next_deg','theta_deg'};
    tolerance=[5e-10,5e-10,5e-9,5e-9,5e-9];maxDelta=zeros(1,5);
    for k=1:5
        maxDelta(k)=max(abs(actual.(fields{k})-reference.(fields{k})));
        assert(maxDelta(k)<=tolerance(k),'pathprofile:Equivalence', ...
            '%s: full-precision %s difference %.15g exceeds %.15g.',name,fields{k},maxDelta(k),tolerance(k));
    end
    assert(all(actual.EDI_elements(2:end)==11316));
    assert(all(actual.true_rel_residual(2:end)<=5e-10));
    assert(all(actual.newSolve(2:end)==(j~=3)));
    tips=tip_table(Path,phases);
    writetable(phases,fullfile(out,[name '_phases.csv']));
    writetable(tips,fullfile(out,[name '_tips.csv']));
    result=struct('wallSeconds',wallSeconds,'phases',phases,'tips',tips, ...
        'maxFullPrecisionDifference',maxDelta,'tolerances',tolerance, ...
        'passed',true,'newPhysicalSolves',sum(actual.newSolve));
    Report.(name)=result;
    save(fullfile(out,'profile_small.mat'),'Report','-v7');
    fprintf('PROFILE %s PASS: wall %.6f s; %d new physical solves.\n',name,wallSeconds,result.newPhysicalSolves);
end
Report.complete=true;save(fullfile(out,'profile_small.mat'),'Report','-v7');
end
function t=tip_table(P,T)
rows=zeros(P.nSegments-1,25);
for k=2:P.nSegments
    r=P.stepResults{k};s=r.summary;f=r.fitTable;
    assert(all(structfun(@logical,r.gates))&&isequal(f.n_native,[38;38;55;55;44;44;34;34]));
    qual=seconds(T,'qualification',k,'TOTAL','total');
    post=sum(T.seconds(strcmp(T.scope,'physical')&T.segment==k&strcmp(T.kind,'phase')& ...
        ismember(T.phase,{'checkpoint_reload','native_cod','cod_fits','physical_edi','postprocessing_gates','compact_save'})));
    rows(k-1,:)=[k,s.T3_elements,s.T6_nodes,qual, ...
        seconds(T,'physical',k,'stiffness_assembly','phase'),s.PCG_iterations, ...
        seconds(T,'physical',k,'pcg','phase'),post,seconds(T,'step',k,'TOTAL','total'), ...
        s.KI_unit,s.KII_unit,s.KII_over_KI,s.theta_current_deg, ...
        s.delta_theta_next_MTS_deg,s.theta_next_local_deg, ...
        seconds(T,'physical',k,'free_dof_extraction','phase'), ...
        seconds(T,'physical',k,'symamd_and_permutation','phase'), ...
        seconds(T,'physical',k,'sgs_construction','phase'), ...
        seconds(T,'physical',k,'native_cod','phase'), ...
        seconds(T,'physical',k,'physical_edi','phase'), ...
        seconds(T,'physical',k,'checkpoint_save','phase'), ...
        seconds(T,'step',k,'qualification_or_candidate_reuse','phase')-qual, ...
        s.true_rel_residual,s.EDI_elements,r.newSolve];
end
t=array2table(rows,'VariableNames',{'segment','T3_elements','T6_nodes', ...
    'qualification_s','assembly_s','PCG_iterations','PCG_s','postprocessing_s','total_step_s', ...
    'KI','KII','KII_over_KI','theta_current_deg','delta_theta_next_deg','theta_next_deg', ...
    'free_dof_s','symamd_and_permutation_s','SGS_s','native_COD_s','physical_EDI_s', ...
    'checkpoint_write_s','candidate_lookup_and_load_s','true_relative_residual','EDI_elements','newSolve'});
end
function v=seconds(T,scope,k,phase,kind)
v=sum(T.seconds(strcmp(T.scope,scope)&T.segment==k&strcmp(T.phase,phase)&strcmp(T.kind,kind)));
end
function h=git_head(root)
[status,h]=system(sprintf('git -C "%s" rev-parse HEAD',root));assert(status==0);h=strtrim(h);
end
