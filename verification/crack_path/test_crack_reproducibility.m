function Report=test_crack_reproducibility(varargin)
% Lightweight parameter/cache/publication checks. No physical solve or PDE mesh.
ip=inputParser;addParameter(ip,'WorkDir',tempname,@(x)ischar(x)||isstring(x));
parse(ip,varargin{:});work=char(ip.Results.WorkDir);mkdir(work);
root=fileparts(fileparts(fileparts(mfilename('fullpath'))));addpath(genpath(root));
set(groot,'defaultFigureVisible','off');
z=load(fullfile(root,'paper','data','accepted_stage1_source.mat'),'R0');original=z.R0;
count=0;coreRows=[];
% Exact preserved Stage-I record and the dimensionless coarse-core family.
for mm=[1,2,4]
    R0=original;R0.summary.a0_reserved_m=1e-3*mm;
    assert(isequaln(R0.C,original.C));
    back=R0;back.summary.a0_reserved_m=original.summary.a0_reserved_m;
    assert(isequaln(back,original));
    da=R0.summary.a0_reserved_m;
    core=build_stage2_scaled_audited_core([0,0],[1,0],da,'Scale',2);
    assert(size(core.local.connect3,1)==3318);
    assert(abs(core.hTip-2*.00675308135*da)<=1e-14&&abs(core.rCore-.75*da)<=1e-14);
    cr=struct('Pmid',[-da,0;0,0],'tipNode',1, ...
        'upperNodes',core.crack.upperT3,'lowerNodes',core.crack.lowerT3);
    mesh=struct('coord',core.local.coord,'connect',core.local.connect);
    mat=struct('E',original.C.E,'nu',original.C.nu,'ps',original.C.ps);
    [r,~,~]=native_COD_polyline_audit(mesh,zeros(2*size(mesh.coord,1),1),mat,cr,8);
    win=[.04,.20;.04,.30;.08,.30;.12,.30];n=zeros(4,1);
    for j=1:4,n(j)=nnz(r/da>=win(j,1)&r/da<=win(j,2));end
    assert(isequal(n,[19;28;23;18]));
    coreRows=[coreRows;mm,core.hTip*1e3,.10*mm,.65*mm,.75*mm,mm,1.25*mm]; %#ok<AGROW>
    count=count+1;
end
carrier=fileread(fullfile(root,'main_stage2_embed_scaled_core_full_domain_theta0.m'));
assert(contains(carrier,'''A0Override'',a0'));
adapter=fileread(fullfile(root,'build_stage2_cracked_mesh_for_theta.m'));
assert(contains(adapter,'''a0'', ip.Results.A0Override'));
count=count+1;
% Physical signatures deliberately exclude only Stage-II a0.
C=original.C;s=struct('C',C,'mat',struct('E',C.E,'nu',C.nu,'ps',C.ps,'D',eye(3),'Dmat',eye(3)));
assert_crack_checkpoint_physics(s,s.mat,C,'test:Physics');count=count+1;
other=C;other.a0=.001;assert(isequaln(crack_physics_signature(C),crack_physics_signature(other)));
for field={'E','nu','ps','A','B'}
    other=C;other.(field{1})=other.(field{1})+1;
    expect(@()assert_crack_checkpoint_physics(s,s.mat,other,'test:Physics'),'test:Physics');count=count+1;
end
other=C;other.load.sig0=2;
expect(@()assert_crack_checkpoint_physics(s,s.mat,other,'test:Physics'),'test:Physics');count=count+1;
expect(@()main_increment_sensitivity_coarse('IncrementMM',4,'TargetLengthMM',93), ...
    'incstudy:TargetNotDivisible');count=count+1;
expect(@()main_increment_sensitivity_coarse('IncrementMM',1,'TargetLengthMM',1), ...
    'incstudy:TargetTooShort');count=count+1;
base=fullfile(root,'verification','crack_path','increment_sensitivity_coarse');
required=fullfile(base,'da_4mm','p1_seed','P1_physical_small.mat');
archiveChecks=0;
if exist(required,'file')==2
    z=load(required,'R');seed=z.R;
    z=load(fullfile(base,'da_4mm','p1_seed','P1_candidate.mat'),'candidate');candidate=z.candidate;
    cp=fullfile(base,'da_4mm','p1_seed','P1_physical_solved.mat');
    validate_increment_study_p1_cache(seed,candidate,original,cp);archiveChecks=archiveChecks+1;
    z=load(fullfile(base,'da_2mm','p1_seed','P1_physical_small.mat'),'R');bad=z.R;
    expect(@()validate_increment_study_p1_cache(bad,candidate,original,cp), ...
        'incstudy:P1ReuseMismatch');archiveChecks=archiveChecks+1;
    bad=seed;bad.exteriorMeshControls.farCapOverA0=.625;
    expect(@()validate_increment_study_p1_cache(bad,candidate,original,cp), ...
        'incstudy:P1ReuseMismatch');archiveChecks=archiveChecks+1;
    for mm=[2,4]
        R0=original;R0.summary.a0_reserved_m=1e-3*mm;
        z=load(fullfile(base,sprintf('da_%dmm',mm),'p1_seed','P1_physical_small.mat'),'R');p1=z.R;
        dirPath=fullfile(base,sprintf('da_%dmm',mm),'trajectory');
        z=load(fullfile(dirPath,'path_run_state.mat'),'State');S=z.State;
        z=load(fullfile(dirPath,'step_003_physical_small.mat'),'R');R=z.R;
        S.vertices=S.vertices(1:4,:);S.thetaDeg=S.thetaDeg(1:3);
        S.completedPhysicalSegments=2;S.nextThetaDeg=S.thetaDeg(3);
        S.rowsThroughCompleted=S.rowsThroughCompleted(1:2,:);
        opts={'FrozenState',R0,'MaxSegments',3,'AllowPhysicalSolves',false, ...
            'ResumeState',S,'ResumeAcceptedResult',R,'CoreScale',2,'RegressionGates',false, ...
            'SeedKI',p1.EDI.KI_unit,'SeedKII',p1.EDI.KII_unit, ...
            'ExteriorFarCapOverIncrement',1.25,'ExteriorTransitionOverIncrement',1, ...
            'ExteriorCalibration',struct('farSlope',.15,'boundaryMetricGrowth',.35)};
        P=run_incremental_crack_path(opts{:},'OutputDir',fullfile(work,sprintf('promote_%d',mm)));
        assert(P.nSegments==3&&P.promotedAcceptedSegment==3);
        assert(isequaln(P.stepTable.KI_unit(3),R.EDI.KI_unit)&& ...
            isequaln(P.stepTable.KII_unit(3),R.EDI.KII_unit));archiveChecks=archiveChecks+1;
        other=R0;other.summary.a0_reserved_m=.001;
        expect(@()run_incremental_crack_path(opts{:},'FrozenState',other, ...
            'OutputDir',fullfile(work,sprintf('wrong_increment_%d',mm))), ...
            'pathrun:ResumeIncrement');archiveChecks=archiveChecks+1;
        corrupt=S;corrupt.exteriorMeshControls.coreScale=1;
        expect(@()run_incremental_crack_path(opts{:},'ResumeState',corrupt, ...
            'OutputDir',fullfile(work,sprintf('wrong_core_%d',mm))), ...
            'pathrun:ResumeMeshFamily');archiveChecks=archiveChecks+1;
        corrupt=S;corrupt.exteriorMeshControls.farCapOverIncrement=NaN;
        expect(@()run_incremental_crack_path(opts{:},'ResumeState',corrupt, ...
            'OutputDir',fullfile(work,sprintf('nan_family_%d',mm))), ...
            'pathrun:ResumeMeshFamily');archiveChecks=archiveChecks+1;
        corrupt=S;corrupt.schemaVersion=99;
        expect(@()run_incremental_crack_path(opts{:},'ResumeState',corrupt, ...
            'OutputDir',fullfile(work,sprintf('wrong_schema_%d',mm))), ...
            'pathrun:ResumeSchema');archiveChecks=archiveChecks+1;
        corrupt=S;corrupt.rowVariableNames([5,6])=corrupt.rowVariableNames([6,5]);
        expect(@()run_incremental_crack_path(opts{:},'ResumeState',corrupt, ...
            'OutputDir',fullfile(work,sprintf('wrong_rows_%d',mm))), ...
            'pathrun:ResumeSchema');archiveChecks=archiveChecks+1;
        [T,V]=load_increment_comparison_run(fullfile(base,sprintf('da_%dmm',mm)),mm);
        assert(max(abs(T.tip_x_m-V.x_m(2:end)))<=2e-12);archiveChecks=archiveChecks+1;
    end
    % Full checkpoint defects reject before EDI/postprocessing, with solves disabled.
    dirPath=fullfile(base,'da_4mm','trajectory');
    z=load(fullfile(dirPath,'step_003_candidate.mat'),'candidate');candidate=z.candidate;
    source=load(fullfile(dirPath,'step_003_physical_solved.mat'));
    bad=source;bad.mat.E=bad.mat.E+1;
    badfile=fullfile(work,'material_mismatch.mat');save(badfile,'-struct','bad','-v7');
    expect(@()solve_incremental_crack_tip(candidate,'FrozenState',original, ...
        'AllowSolve',false,'CheckpointFile',badfile,'SaveFile',fullfile(work,'unused.mat')), ...
        'pathsolve:CheckpointMismatch');archiveChecks=archiveChecks+1;
    bad=source;bad.mesh.coord(end,1)=bad.mesh.coord(end,1)+1e-6;
    badfile=fullfile(work,'t6_mismatch.mat');save(badfile,'-struct','bad','-v7');
    expect(@()solve_incremental_crack_tip(candidate,'FrozenState',original, ...
        'AllowSolve',false,'CheckpointFile',badfile,'SaveFile',fullfile(work,'unused.mat')), ...
        'pathsolve:CheckpointMismatch');archiveChecks=archiveChecks+1;
    % Publication CSV identity/uniqueness mutants use only isolated copies.
    pub=fullfile(work,'publication_fixture');mkdir(pub);
    CopySource=fullfile(base,'da_4mm');
    for name={'states.csv','vertices.csv','increment_study_summary.mat'}
        copyfile(fullfile(CopySource,name{1}),fullfile(pub,name{1}));
    end
    z=load(fullfile(pub,'increment_study_summary.mat'),'M');M=z.M;M.controls.coreScale=1;
    save(fullfile(pub,'increment_study_summary.mat'),'M');
    expect(@()load_increment_comparison_run(pub,4),'incfig:FamilyMismatch');archiveChecks=archiveChecks+1;
    copyfile(fullfile(CopySource,'increment_study_summary.mat'),fullfile(pub,'increment_study_summary.mat'));
    T=readtable(fullfile(pub,'states.csv'));T=[T;T(3,:)];writetable(T,fullfile(pub,'states.csv'));
    expect(@()load_increment_comparison_run(pub,4),'incfig:GridMismatch');archiveChecks=archiveChecks+1;
    copyfile(fullfile(CopySource,'states.csv'),fullfile(pub,'states.csv'));
    T=readtable(fullfile(pub,'states.csv'));T.KI_unit(3)=T.KI_unit(3)+.01;
    writetable(T,fullfile(pub,'states.csv'));
    expect(@()load_increment_comparison_run(pub,4),'incfig:ArchiveMismatch');archiveChecks=archiveChecks+1;
else
    fprintf('Investigator-local archive absent: archive-specific checks skipped.\n');
end
% Synthetic CSV format fixtures exercise both availability branches and
% exact common-length matching; these numbers are not physical results.
dirs=cell(1,3);levels=[4,2,1];
for il=1:3
    mm=levels(il);dirName=fullfile(work,sprintf('synthetic_%dmm',mm));mkdir(dirName);dirs{il}=dirName;
    segment=(1:12/mm)';a=mm*segment;
    T=table(segment,.2+1e-3*a,zeros(size(a)),zeros(size(a)),ones(size(a)), ...
        1e-3*sin(a),1e-3*sin(a),.02*mm*ones(size(a)),zeros(size(a)),a, ...
        'VariableNames',{'segment','tip_x_m','tip_y_m','theta_deg','KI_unit','KII_unit', ...
        'KII_over_KI','delta_theta_next_deg','theta_next_deg','crack_length_mm'});
    vertex=(0:height(T))';V=table(vertex,mm*vertex,.2+1e-3*mm*vertex,zeros(size(vertex)), ...
        'VariableNames',{'vertex','crack_length_mm','x_m','y_m'});
    writetable(T,fullfile(dirName,'states.csv'));writetable(V,fullfile(dirName,'vertices.csv'));
    M=struct('study','crack_increment_sensitivity_coarse_family','increment_mm',mm, ...
        'controls',struct('coreScale',2,'farCapOverIncrement',1.25,'transitionOverIncrement',1, ...
        'farSlope',.15,'boundaryMetricGrowth',.35));
    T.pass=true(height(T),1);M.Path=struct('stepTable',T,'vertices',[V.x_m,V.y_m]);
    save(fullfile(dirName,'increment_study_summary.mat'),'M');
end
S=plot_increment_sensitivity_coarse_publication('Run4Dir',dirs{1},'Run2Dir',dirs{2}, ...
    'Run1Dir',fullfile(work,'absent_1mm'),'OutputDir',fullfile(work,'plot_two'),'Export',false);
assert(~S.has1mm&&isequal(S.metrics.crack_length_mm,[4;8;12]));count=count+1;
S=plot_increment_sensitivity_coarse_publication('Run4Dir',dirs{1},'Run2Dir',dirs{2}, ...
    'Run1Dir',dirs{3},'OutputDir',fullfile(work,'plot_three'),'Export',false);
assert(S.has1mm&&isequal(S.metrics.crack_length_mm,[4;8;12]));
assert(max(abs(S.native.turn_density_deg_per_mm-.02))<1e-14);count=count+1;
close all;
Report=struct('portableChecks',count,'archiveChecks',archiveChecks,'passed',true, ...
    'noPhysicalSolve',true,'coreFlow',array2table(coreRows,'VariableNames', ...
    {'increment_mm','hTip_mm','rInner_mm','rOuter_mm','rCore_mm','transition_mm','farCap_mm'}));
save(fullfile(work,'checks.mat'),'Report');
fid=fopen(fullfile(work,'checks.json'),'w');fprintf(fid,'%s\n',jsonencode(Report));fclose(fid);
disp(Report.coreFlow);fprintf('REPRODUCIBILITY CHECKS PASS: %d portable, %d archive checks; no solve.\n',count,archiveChecks);
end
function expect(fn,id)
try,fn();catch err,assert(strcmp(err.identifier,id),'Expected %s; got %s: %s',id,err.identifier,err.message);return,end
error('reprotest:ExpectedFailure','Expected %s.',id);
end
