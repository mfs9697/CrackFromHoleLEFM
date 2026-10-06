function E=export_pilot_evidence(varargin)
% Manuscript-only reader/auditor. No mesh generation, EDI replay, or FE solve.
ip=inputParser;
addParameter(ip,'RunDir','',@(x)ischar(x)||isstring(x));
addParameter(ip,'FrozenStateFile','',@(x)ischar(x)||isstring(x));
addParameter(ip,'OutputDir','',@(x)ischar(x)||isstring(x));
parse(ip,varargin{:});o=ip.Results;
root=fileparts(fileparts(mfilename('fullpath')));addpath(genpath(root));
run=char(o.RunDir);if isempty(run),run=fullfile(root,'verification','crack_path','final_clean_run');end
out=char(o.OutputDir);if isempty(out),out=fullfile(root,'paper','data');end
d=load(fullfile(run,'path_run_state.mat'),'State');S=d.State;
assert(S.completedPhysicalSegments==23&&size(S.vertices,1)==25&&S.rowHistoryComplete, ...
    'pilot:ArchiveShape','Stop: authoritative archive classification differs from the request.');
T=array2table(S.rowsThroughCompleted,'VariableNames',cellstr(string(S.rowVariableNames)));
assert(height(T)==23&&isequal(T.segment,(1:23)'));
assert(S.regression.theta2_pass&&S.regression.step2_pass, ...
    'pilot:Regression','Stop: stored initial/P2 regression did not pass.');
seg=diff(S.vertices,1,1);da=median(vecnorm(seg,2,2));
assert(abs(da-.004)<2e-12&&max(abs(vecnorm(seg,2,2)-da))<2e-12);
T.crack_length_mm=T.segment*da*1e3;
cp=load(fullfile(run,'step_002_physical_solved.mat'),'C','mat');
assert(cp.C.A==.30&&cp.C.B==.10&&cp.mat.E==210000&&cp.mat.nu==.30&&cp.mat.ps==1);
physical=dir(fullfile(run,'step_*_physical_small.mat'));
qualification=dir(fullfile(run,'step_*_qualification_small.mat'));
assert(numel(physical)==22&&numel(qualification)==23);
assert(exist(fullfile(run,'step_024_physical_small.mat'),'file')~=2&& ...
    exist(fullfile(run,'step_024_physical_solved.mat'),'file')~=2);
COD=table();Q=table();maxDiff=zeros(1,6);synthetic=[];
stateFields={'KI_unit','KII_unit','KII_over_KI','theta_deg','delta_theta_next_deg','theta_next_deg'};
resultFields={'KI_unit','KII_unit','KII_over_KI','theta_current_deg','delta_theta_next_MTS_deg','theta_next_local_deg'};
for k=2:24
    d=load(fullfile(run,sprintf('step_%03d_qualification_small.mat',k)),'Small');s=d.Small;
    assert(s.pass&&all(structfun(@logical,s.gates))&&all(structfun(@logical,s.syntheticGates)));
    assert(height(s.synthetic)==3&&all(s.synthetic.nElem_used==11316));
    assert(isequal(s.sampleCounts.nativePoints,[38;55;44;34]));
    row=s.summary;row.segment=k;Q=[Q;row]; %#ok<AGROW>
    synthetic=[synthetic;[k,max(abs(s.synthetic.KI_error)),max(abs(s.synthetic.KII_error)), ...
        abs(s.synthetic.KII_recovered(3)/1e-4-1)]]; %#ok<AGROW>
    if k==24,continue,end
    d=load(fullfile(run,sprintf('step_%03d_physical_small.mat',k)),'R');r=d.R;
    assert(r.pass&&all(structfun(@logical,r.gates))&&r.solverInfo.flag==0);
    assert(r.solverInfo.pcgTol==1e-10&&r.solverInfo.pcgMaxIt==5000&& ...
        r.solverInfo.relres<=1e-10&&r.solverInfo.trueRelResidual<=5e-10);
    assert(r.EDI.EDI_elements==11316&&isequal(r.fitTable.n_native,[38;38;55;55;44;44;34;34]));
    assert(norm(r.pathFixed-S.vertices(1:k+1,:),'fro')<=2e-12);
    for j=1:6
        delta=abs(T.(stateFields{j})(k)-r.summary.(resultFields{j}));
        maxDiff(j)=max(maxDiff(j),delta);
        tol=5e-9;if j<=3,tol=5e-10;end
        assert(delta<=tol,'pilot:InconsistentHistory','Stop: state/result disagreement P%d %s.',k,stateFields{j});
    end
    [~,turn]=kink_angle_LEFM_MTS(r.EDI.KI_unit,r.EDI.KII_unit);
    assert(abs(turn-T.delta_theta_next_deg(k))<=5e-9);
    f=r.fitTable;f.segment=repmat(k,height(f),1);
    f.crack_length_mm=repmat(T.crack_length_mm(k),height(f),1);
    f.EDI_ratio=repmat(T.KII_over_KI(k),height(f),1);
    f.EDI_turn_deg=repmat(T.delta_theta_next_deg(k),height(f),1);
    f.ratio_error=f.ratio_COD-f.EDI_ratio;
    f.turn_error_deg=f.delta_theta_next_MTS_deg-f.EDI_turn_deg;
    COD=[COD;f]; %#ok<AGROW>
end
assert(height(COD)==176&&all(diff(T.KI_unit)>0));
[qPeak,iQ]=max(T.KII_over_KI);[KIIpeak,iK]=max(T.KII_unit);
assert(T.segment(iQ)==17&&T.segment(iK)==17);
assert(all(T.KII_unit(1:21)>0)&&all(T.KII_unit(22:23)<0));
assert(all(diff(T.KII_unit(17:23))<0));
fraction=-T.KII_over_KI(21)/(T.KII_over_KI(22)-T.KII_over_KI(21));
aLS=T.crack_length_mm(21)+fraction*(T.crack_length_mm(22)-T.crack_length_mm(21));
assert(abs(aLS-84.146369911941)<1e-9);
LS=S.vertices(22,:)+fraction*(S.vertices(23,:)-S.vertices(22,:));
E=struct('sourceCommit','cfdf0f110010d688a4e3c48f6d88a00fd17dc698', ...
    'C',cp.C,'material',cp.mat,'increment_m',da,'stateRows',table2struct(T), ...
    'vertices_m',S.vertices,'theta_deg',S.thetaDeg,'fastEDI',S.fastEDI, ...
    'qualification',table2struct(Q),'codFits',table2struct(COD), ...
    'modeMixityPeakSegment',T.segment(iQ),'modeMixityPeak',qPeak, ...
    'KIIPeakSegment',T.segment(iK),'KIIPeak',KIIpeak,'interpolatedZeroLength_mm',aLS, ...
    'interpolatedZeroPoint_m',LS,'interpolationFraction',fraction);
E.audit=struct('acceptedPhysicalCount',22,'qualificationCount',23,'CODFitCount',176, ...
    'completedPhysicalSegments',23,'qualifiedUnsolvedSegment',24, ...
    'comparisonFields',{stateFields},'maxStateResultDifferences',maxDiff, ...
    'maxPCGReportedResidual',max(T.PCG_relres(2:end)), ...
    'maxTrueResidual',max(T.true_rel_residual(2:end)), ...
    'minPCGIterations',min(T.PCG_iterations(2:end)),'maxPCGIterations',max(T.PCG_iterations(2:end)), ...
    'maxCODTurnDifference_deg',max(abs(COD.turn_error_deg)), ...
    'maxCODRatioDifference',max(abs(COD.ratio_error)), ...
    'syntheticErrors',synthetic,'noPhysicalSolvePerformed',true);
E.regression=S.regression;
frozen=char(o.FrozenStateFile);
if ~isempty(frozen)
    d=load(frozen,'R0');r=d.R0;
    assert(r.stage1Pass&&isequaln(r.C,cp.C));
    p0=[r.summary.x_star_m,r.summary.y_star_m];
    n=[r.summary.nmat_x,r.summary.nmat_y];
    assert(norm(p0-S.vertices(1,:))<=2e-12&&norm(seg(1,:)/da-n)<=2e-12);
    E.stage1Summary=table2struct(r.summary);E.stage1Gates=r.gates;E.stage1Method=r.method;
    E.stage1Source='Investigator-saved accepted_R0.mat, saved for the preceding profiling task';
end
if exist(out,'dir')~=7,mkdir(out);end
writetable(T,fullfile(out,'accepted_states.csv'));writetable(COD,fullfile(out,'cod_fits.csv'));
writetable(Q,fullfile(out,'qualification.csv'));
fid=fopen(fullfile(out,'evidence.json'),'w');fprintf(fid,'%s\n',jsonencode(E));fclose(fid);
save(fullfile(out,'evidence_exact.mat'),'E','T','COD','Q','-v7');
fprintf('PILOT AUDIT PASS: 22 accepted physical files, 23 qualifications, 176 COD fits; no solve.\n');
disp(E.audit);fprintf('Linear zero estimate %.17g mm; not a solved state.\n',aLS);
end
