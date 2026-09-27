function Out=main_step3e_compare_edi_weight_functions()
% Step 3E: compare legacy clipped analytic q with FE-consistent nodal q.
here=fileparts(mfilename('fullpath')); addpath(genpath(fileparts(fileparts(here))));
fprintf('\n============================================================\n');
fprintf('SIF AUDIT STEP 3E: EDI WEIGHT-FUNCTION / QUADRATURE CHECK\n');
fprintf('============================================================\n');
NrL=[8 16 32]; NtL=[64 128 256]; n=numel(NrL);
names={"analytic_radial","fe_nodal"};
rows=[]; detail=cell(n,2);
for k=1:n
  R=validate_EDI_Williams_fields('NrList',NrL(k),'NthList',NtL(k), ...
    'rMeshInner',0.005,'rMeshOuter',0.20,'rInner',0.024,'rOuter',0.12, ...
    'Verbose',false,'AssertFine',false,'MeshTopology','mirror_reflected');
  mesh=R.meshes{1}; U=R.displacements{1,3}; mat=R.material;
  for j=1:2
    W=run_s0_annulus_sweeps_for_mesh(mesh,U,mat,names{j}); detail{k,j}=W;
    q=W.range;
    rows=[rows; k,NrL(k),NtL(k),j,q,max(q)]; %#ok<AGROW>
  end
end
T=array2table(rows,'VariableNames',{'level','Nr','Ntheta','weightID', ...
 'outer_KI_range','outer_KII_range','inner_KI_range','inner_KII_range','max_annulus_range'});
WName=strings(height(T),1); for i=1:height(T), WName(i)=names{T.weightID(i)}; end
T.weightFunction=WName; T=movevars(T,'weightFunction','After','weightID');
fprintf('\nANNULUS-SENSITIVITY COMPARISON\n'); disp(T);
fprintf('\nDetailed finest-mesh sweeps (Nr=32, Ntheta=256)\n');
for j=1:2
 fprintf('\n%s: outer sweep\n',names{j}); disp(detail{end,j}.outer);
 fprintf('\n%s: inner sweep\n',names{j}); disp(detail{end,j}.inner);
end
Out=struct('summary',T,'detail',{detail},'NrList',NrL,'NthList',NtL);
fprintf('\nSTEP 3E completed.\n');
fprintf('Interpretation: reduced and smoother ranges for fe_nodal indicate annulus-boundary quadrature aliasing in the clipped analytic-q implementation.\n');
end
