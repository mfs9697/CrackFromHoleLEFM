function R66=main_step66_solver_memory_preflight(varargin)
%MAIN_STEP66_SOLVER_MEMORY_PREFLIGHT
% Symbolic memory preflight for Level-0 and Level-1 C03 meshes.
%
% STRICTLY NO PHYSICAL SOLVE:
%   - no physical U is loaded;
%   - no stiffness VALUES are assembled;
%   - no load vector is assembled;
%   - no backslash/PCG/iterative solve is called.
%
% The routine reconstructs the deterministic C03 meshes in memory, upgrades
% their T3 connectivity to T6, builds only a NODE-ADJACENCY sparsity graph,
% and uses symamd + symbfact on that graph to estimate direct-factor fill.
%
% Important solver observation:
% current stif_assem.m imposes homogeneous Dirichlet conditions by clamping
% rows only. This makes the stored K nonsymmetric, so K\F is not guaranteed
% to use the memory-favorable SPD Cholesky path even though the underlying
% elasticity problem is SPD. Step66 therefore reports:
%   (a) exact assembly-triplet memory implied by current stif_assem.m;
%   (b) exact block-structural K sparsity implied by the mesh;
%   (c) a block-symbolic SPD Cholesky fill estimate;
%   (d) planning envelopes, explicitly heuristic, for deciding whether a
%       Level-0 solver qualification step is needed.
%
% Step66 never changes the production solver.
%
root=fileparts(fileparts(fileparts(mfilename('fullpath'))));
vdir=fullfile(root,'verification');

ip=inputParser;
addParameter(ip,'SourceCandidateFile', ...
    fullfile(vdir,'step62_structured_graded_mesh_candidate_T3.mat'), ...
    @(s)ischar(s)||(isstring(s)&&isscalar(s)));
addParameter(ip,'AvailableMemoryGiB',NaN, ...
    @(x)isnumeric(x)&&isscalar(x)&&(isnan(x)||(isfinite(x)&&x>0)));
addParameter(ip,'SaveFile', ...
    fullfile(vdir,'step66_solver_memory_preflight_small_data.mat'), ...
    @(s)ischar(s)||(isstring(s)&&isscalar(s)));
parse(ip,varargin{:});
opt=ip.Results;

addpath(genpath(root));
assert_step66_branch(root);

sourceFile=char(opt.SourceCandidateFile);
saveFile=char(opt.SaveFile);
if exist(sourceFile,'file')~=2
    error('step66:MissingSource', ...
        'Missing committed archived Step62 candidate: %s',sourceFile);
end
if ~matlab_callable('symamd') || ~matlab_callable('symbfact')
    error('step66:MissingSymbolicTools', ...
        ['Step66 requires callable MATLAB symamd and symbfact. ', ...
         'Availability accepts M-files, MEX/P-code, and built-ins.']);
end

cal=c03_calibration();

fprintf('\n============================================================\n');
fprintf('STEP 66: SOLVER-MEMORY PREFLIGHT — SYMBOLIC ONLY\n');
fprintf('============================================================\n');
fprintf('  NO physical FEM solve. NO stiffness-value assembly.\n');
fprintf('  NO load assembly. NO backslash. NO iterative solve.\n');
fprintf('  Reconstructing deterministic C03 Level 0 and Level 1 in memory.\n');

% -------------------------------------------------------------------------
% Level 0: exact accepted C03 baseline.
% -------------------------------------------------------------------------
O0=main_step62_structured_graded_mesh( ...
    'SourceCandidateFile',sourceFile, ...
    'SavePrefix',fullfile(vdir,'step66_internal_L0'), ...
    'Visible','off', ...
    'RunSynthetic',false, ...
    'ExteriorCalibration',cal, ...
    'WriteArtifacts',false, ...
    'ReturnCandidate',true, ...
    'Verbose',false, ...
    'Level',0);

if ~O0.gates.structuralPass || O0.summary.newT3~=32980 || ...
        O0.summary.newT6Nodes~=66854 || ...
        abs(O0.maxNeighborSizeRatio-1.79678451)>5e-7
    error('step66:Level0Mismatch', ...
        'Current family definition does not reproduce accepted C03 Level 0.');
end
fprintf('  Level 0 reproduced: T3=%d, T6=%d.\n', ...
    O0.summary.newT3,O0.summary.newT6Nodes);

M0=symbolic_mesh_memory(O0.candidate.p,O0.candidate.t,0);
L0summary=compact_level_summary(O0);

% Release Level-0 geometry before Level-1 symbolic work.
O0candidate=O0.candidate; %#ok<NASGU>
O0=rmfield(O0,'candidate');
clear O0candidate

% -------------------------------------------------------------------------
% Level 1: same qualified family, half structured scale.
% -------------------------------------------------------------------------
O1=main_step62_structured_graded_mesh( ...
    'SourceCandidateFile',sourceFile, ...
    'SavePrefix',fullfile(vdir,'step66_internal_L1'), ...
    'Visible','off', ...
    'RunSynthetic',false, ...
    'ExteriorCalibration',cal, ...
    'WriteArtifacts',false, ...
    'ReturnCandidate',true, ...
    'Verbose',false, ...
    'Level',1);

if ~O1.gates.structuralPass || abs(O1.design.scale-.5)>1e-14 || ...
        abs(O1.summary.pairedRadius_mm-6)>1e-12 || ...
        O1.maxNeighborSizeRatio>1.8+5e-12 || ...
        ~O1.gates.completeT3Pairing || ~O1.gates.completeT6Pairing
    error('step66:Level1Mismatch', ...
        'Level 1 no longer matches the qualified C03 family definition.');
end
fprintf('  Level 1 reproduced: T3=%d, T6=%d, max ratio=%.9g.\n', ...
    O1.summary.newT3,O1.summary.newT6Nodes,O1.maxNeighborSizeRatio);

M1=symbolic_mesh_memory(O1.candidate.p,O1.candidate.t,1);
L1summary=compact_level_summary(O1);
O1=rmfield(O1,'candidate');

% -------------------------------------------------------------------------
% Machine-memory observation. Override wins when supplied.
% -------------------------------------------------------------------------
[availGiB,totalGiB,memSource]=available_memory_gib(opt.AvailableMemoryGiB);

MemoryTable=struct2table([M0;M1]);
scaleDOF=M1.ndof/M0.ndof;
scaleTriplets=M1.assemblyTripletGiB/M0.assemblyTripletGiB;
scaleFactor=M1.spdFactorGiB/M0.spdFactorGiB;
scaleHighEnvelope=M1.spdPlanningHighGiB/M0.spdPlanningHighGiB;

ScalingTable=table(scaleDOF,scaleTriplets,scaleFactor,scaleHighEnvelope, ...
    'VariableNames',{'dofRatio_L1_over_L0','tripletMemoryRatio', ...
    'symbolicFactorMemoryRatio','planningHighRatio'});

% Current production matrix is row-clamped only and hence structurally
% nonsymmetric after BC application. This is a source-code fact, not a
% measurement from a physical solve.
currentRowOnlyClamp=true;
currentBackslashSPDPathGuaranteed=false;

if isnan(availGiB)
    riskClass="UNKNOWN_NO_MACHINE_MEMORY";
    if M1.spdPlanningHighGiB<=8
        recommendation="QUALIFY_SYMMETRIC_DIRECT_ON_LEVEL0";
    else
        recommendation="MEASURE_MACHINE_MEMORY_THEN_CHOOSE_SOLVER";
    end
else
    assemblyRatio=M1.assemblyPeakProxyGiB/availGiB;
    spdHighRatio=M1.spdPlanningHighGiB/availGiB;
    if assemblyRatio>=0.70
        riskClass="HIGH_ASSEMBLY_MEMORY_RISK";
        recommendation="QUALIFY_MEMORY_SAVING_ASSEMBLY_AND_ITERATIVE_SOLVER_ON_LEVEL0";
    elseif spdHighRatio>=0.70
        riskClass="HIGH_DIRECT_FACTORIZATION_RISK";
        recommendation="QUALIFY_ITERATIVE_SOLVER_ON_LEVEL0";
    elseif spdHighRatio>=0.45
        riskClass="CAUTION_DIRECT_SOLVE";
        recommendation="QUALIFY_SYMMETRIC_DIRECT_ON_LEVEL0_BEFORE_LEVEL1";
    else
        riskClass="SPD_DIRECT_MEMORY_APPEARS_PLAUSIBLE";
        recommendation="QUALIFY_SYMMETRIC_DIRECT_ON_LEVEL0_BEFORE_LEVEL1";
    end
end

MachineMemory=table(availGiB,totalGiB,string(memSource), ...
    'VariableNames',{'availableGiB','totalPhysicalGiB','source'});

Decision=table(string(riskClass),string(recommendation), ...
    currentRowOnlyClamp,currentBackslashSPDPathGuaranteed, ...
    'VariableNames',{'riskClass','recommendedNextQualification', ...
    'currentStifAssemRowOnlyClamp','currentBackslashSPDPathGuaranteed'});

fprintf('\nSYMBOLIC MEMORY PREFLIGHT\n');
disp(MemoryTable);
fprintf('\nLEVEL-1 / LEVEL-0 SCALING\n');
disp(ScalingTable);
fprintf('\nMACHINE MEMORY OBSERVATION\n');
disp(MachineMemory);
fprintf('\nDECISION\n');
disp(Decision);

fprintf('\nINTERPRETATION NOTES\n');
fprintf('  assemblyTripletGiB: exact storage of rw/cl/st triplet arrays only.\n');
fprintf('  KpatternGiB: approximate MATLAB sparse-double storage from structural nnz.\n');
fprintf('  spdFactorGiB: block-symbolic Cholesky factor-storage estimate after symamd.\n');
fprintf('  assemblyPeakProxyGiB = 2*triplets + K pattern (planning proxy).\n');
fprintf('  spdPlanningLowGiB  = triplets + K + 2*factor (heuristic).\n');
fprintf('  spdPlanningHighGiB = 2*triplets + K + 4*factor (conservative heuristic).\n');
fprintf('  Current row-only BC clamping makes K nonsymmetric; actual K\\F memory may\n');
fprintf('  exceed the SPD estimates. Step66 therefore NEVER authorizes Level-1 solve.\n');

R66=struct( ...
    'MemoryTable',MemoryTable, ...
    'ScalingTable',ScalingTable, ...
    'MachineMemory',MachineMemory, ...
    'Decision',Decision, ...
    'level0Summary',L0summary, ...
    'level1Summary',L1summary, ...
    'calibration',cal, ...
    'currentRowOnlyClamp',currentRowOnlyClamp, ...
    'currentBackslashSPDPathGuaranteed',currentBackslashSPDPathGuaranteed, ...
    'noPhysicalFEM',true, ...
    'noStiffnessValueAssembly',true, ...
    'noLoadAssembly',true, ...
    'noLinearSolve',true, ...
    'interpretation',['Symbolic memory preflight only. Any solver change ', ...
      'must be separately qualified on Level 0 before a Level-1 physical solve.']);

save(saveFile,'R66','-v7');
fprintf('  Compact Step66 result saved: %s\n',saveFile);
fprintf('STEP66 complete: symbolic preflight only; zero physical solves.\n');
end

% =========================================================================
function M=symbolic_mesh_memory(P,T3,level)
[P6,T6]=T3toT6_fast(P,T3);
nNode=size(P6,1);
nElem=size(T6,1);
ndof=2*nNode;

fprintf('\n  Level %d symbolic graph: T6 nodes=%d, elements=%d, DOF=%d.\n', ...
    level,nNode,nElem,ndof);

% Each T6 element contains 15 unordered distinct-node pairs. Build the
% node-block adjacency only; no element stiffness values are evaluated.
pairs=nchoosek(1:6,2);
np=size(pairs,1);
I=zeros(np*nElem,1);
J=zeros(np*nElem,1);
for k=1:np
    q=(k-1)*nElem+(1:nElem);
    I(q)=T6(:,pairs(k,1));
    J(q)=T6(:,pairs(k,2));
end
V=ones(size(I));
A=sparse(I,J,V,nNode,nNode);
clear I J V
A=spones(A+A'+speye(nNode));

nnzNodePattern=nnz(A);
nnzKPattern=4*nnzNodePattern; % dense 2x2 block per structural node pair
nnzKLower=2*nnzNodePattern+nNode;

fprintf('    node-block pattern nnz=%d; running symamd/symbfact...\n', ...
    nnzNodePattern);
p=symamd(A);
counts=symbfact(A(p,p));
nnzLNode=sum(double(counts));

% Every off-diagonal lower node block contributes at most four scalar
% factor entries; every diagonal 2x2 lower block contributes three.
nnzLScalar=4*nnzLNode-nNode;
fillRatioLower=nnzLScalar/nnzKLower;

GiB=1024^3;
bytesPerSparseNZ=16; % approximate: 8-byte value + 8-byte row index
bytesPerColPtr=8;
assemblyTripletEntries=144*nElem;
assemblyTripletGiB=(3*8*assemblyTripletEntries)/GiB;
KpatternGiB=(bytesPerSparseNZ*nnzKPattern+bytesPerColPtr*(ndof+1))/GiB;
spdFactorGiB=(bytesPerSparseNZ*nnzLScalar+bytesPerColPtr*(ndof+1))/GiB;

% Planning proxies: deliberately labeled heuristics, not MATLAB guarantees.
assemblyPeakProxyGiB=2*assemblyTripletGiB+KpatternGiB;
spdPlanningLowGiB=assemblyTripletGiB+KpatternGiB+2*spdFactorGiB;
spdPlanningHighGiB=2*assemblyTripletGiB+KpatternGiB+4*spdFactorGiB;

M=struct( ...
    'level',level, ...
    'nT3',size(T3,1), ...
    'nT6Nodes',nNode, ...
    'ndof',ndof, ...
    'assemblyTripletEntries',assemblyTripletEntries, ...
    'assemblyTripletGiB',assemblyTripletGiB, ...
    'nnzNodePattern',nnzNodePattern, ...
    'nnzKPattern',nnzKPattern, ...
    'KpatternGiB',KpatternGiB, ...
    'nnzLNode',nnzLNode, ...
    'nnzLScalarEstimate',nnzLScalar, ...
    'fillRatioLower',fillRatioLower, ...
    'spdFactorGiB',spdFactorGiB, ...
    'assemblyPeakProxyGiB',assemblyPeakProxyGiB, ...
    'spdPlanningLowGiB',spdPlanningLowGiB, ...
    'spdPlanningHighGiB',spdPlanningHighGiB);
end

function S=compact_level_summary(O)
S=struct( ...
    'nT3',O.summary.newT3, ...
    'nT6',O.summary.newT6Nodes, ...
    'tipMedian_mm',O.summary.newTipMedian_mm, ...
    'pairedRadius_mm',O.summary.pairedRadius_mm, ...
    'maxNeighborRatio',O.maxNeighborSizeRatio, ...
    'patchMinAngle_deg',O.summary.patchMinAngle_deg, ...
    'exteriorMinAngle_deg',O.summary.exteriorMinAngle_deg, ...
    'structuralPass',logical(O.gates.structuralPass));
end

function [availGiB,totalGiB,source]=available_memory_gib(override)
GiB=1024^3;
availGiB=NaN;totalGiB=NaN;source='unavailable';
if isfinite(override)
    availGiB=override;source='user override';return
end
if ~ispc,return,end
try
    [u,s]=memory; %#ok<ASGLU>
    if isfield(s,'PhysicalMemory')
        pm=s.PhysicalMemory;
        if isfield(pm,'Available')&&isfinite(pm.Available)
            availGiB=double(pm.Available)/GiB;
            source='MATLAB memory: physical available';
        end
        if isfield(pm,'Total')&&isfinite(pm.Total)
            totalGiB=double(pm.Total)/GiB;
        end
    end
    if isnan(availGiB)&&isfield(u,'MemAvailableAllArrays')&& ...
            isfinite(u.MemAvailableAllArrays)
        availGiB=double(u.MemAvailableAllArrays)/GiB;
        source='MATLAB memory: MemAvailableAllArrays';
    end
catch
    % Leave NaN: machine-memory reporting is optional.
end
end

function cal=c03_calibration()
cal=struct( ...
    'transitionLength_m',.008, ...
    'farSlope',.10, ...
    'boundaryMetricGrowth',.25, ...
    'smoothingSteps',6, ...
    'refinementMaxPasses',140, ...
    'refinementMinAngle_deg',25, ...
    'refinementLongestFactor',1.65, ...
    'neighborRatioTarget',1.8, ...
    'verbose',false);
end

function tf=matlab_callable(name)
tf=exist(name,'file')~=0 || exist(name,'builtin')~=0;
end

function assert_step66_branch(root)
[status,b]=system(sprintf('git -C "%s" branch --show-current',root));
assert(status==0&&strcmp(strtrim(b),'audit/step66-solver-memory-preflight'), ...
    'step66:Branch', ...
    'Run Step66 only on audit/step66-solver-memory-preflight.');
end
