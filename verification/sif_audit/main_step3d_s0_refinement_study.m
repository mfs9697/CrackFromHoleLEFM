function Out=main_step3d_s0_refinement_study(varargin)
% Step 3D: refine canonical S0 and print/plot total mesh for every case.
ip=inputParser; addParameter(ip,'PlotMeshes',true,@(x)islogical(x)||isnumeric(x));
addParameter(ip,'ShowT6Nodes',true,@(x)islogical(x)||isnumeric(x)); parse(ip,varargin{:}); o=ip.Results;
here=fileparts(mfilename('fullpath')); addpath(genpath(fileparts(fileparts(here))));
fprintf('\n============================================================\n');
fprintf('SIF AUDIT STEP 3D: CANONICAL S0 MESH REFINEMENT\n');
fprintf('============================================================\n');
NrL=[8 16 32]; NtL=[64 128 256]; n=numel(NrL);
M=nan(n,10); B=nan(n,10); C=nan(n,8); SW=cell(n,1); Mesh=cell(n,1); Base=cell(n,1);
for k=1:n
  Nr=NrL(k); Nt=NtL(k);
  R=validate_EDI_Williams_fields('NrList',Nr,'NthList',Nt, ...
    'rMeshInner',0.005,'rMeshOuter',0.20,'rInner',0.024,'rOuter',0.12, ...
    'Verbose',false,'AssertFine',false,'MeshTopology','mirror_reflected');
  m=R.meshes{1}; a=R.meshAudit{1}; Mesh{k}=m; Base{k}=R;
  M(k,:)=[k,Nr,Nt,size(m.coord3,1),size(m.connect3,1),size(m.coord,1), ...
    size(m.connect,1),2*size(m.coord,1),a.maxMirrorCoordError,double(a.crackFacesDistinct)];
  T=R.table; pI=T(T.caseName=="pure_I",:); pII=T(T.caseName=="pure_II",:); pm=T(T.caseName=="mixed_I_II",:);
  B(k,:)=[k,Nr,Nt,pI.KI_recovered_over_input,abs(pI.KII_recovered), ...
    pII.KII_recovered_over_input,abs(pII.KI_recovered),pm.KI_recovered_over_input, ...
    pm.KII_recovered_over_input,hypot(pm.KI_recovered-1,pm.KII_recovered-0.35)/hypot(1,0.35)];
  SW{k}=run_s0_annulus_sweeps_for_mesh(m,R.displacements{1,3},R.material);
  q=SW{k}.range; C(k,:)=[k,Nr,Nt,q,max(q)];
  if logical(o.PlotMeshes), plot_s0_full_mesh(m,a,Nr,Nt,logical(o.ShowT6Nodes)); end
end
Tm=array2table(M,'VariableNames',{'level','Nr','Ntheta','nT3_vertices','nT3_elements','nT6_nodes','nT6_elements','nDOF','max_mirror_error','crack_faces_distinct'});
Tb=array2table(B,'VariableNames',{'level','Nr','Ntheta','pureI_KI_ratio','pureI_abs_crossKII','pureII_KII_ratio','pureII_abs_crossKI','mixed_KI_ratio','mixed_KII_ratio','mixed_vector_error_rel'});
Tc=array2table(C,'VariableNames',{'level','Nr','Ntheta','outer_KI_range','outer_KII_range','inner_KI_range','inner_KII_range','max_annulus_range'});
fprintf('\nTOTAL MESH FOR ALL S0 CASES\n'); disp(Tm);
fprintf('\nBASELINE EXACT-FIELD RECOVERY\n'); disp(Tb);
fprintf('\nANNULUS SENSITIVITY VS REFINEMENT\n'); disp(Tc);
for k=1:n
 fprintf('\nLevel %d outer sweep\n',k); disp(SW{k}.outer);
 fprintf('\nLevel %d inner sweep\n',k); disp(SW{k}.inner);
end
Out=struct('meshSummary',Tm,'baseline',Tb,'refinement',Tc,'sweeps',{SW},'meshes',{Mesh},'baselineRuns',{Base});
fprintf('\nSTEP 3D completed.\n');
end
