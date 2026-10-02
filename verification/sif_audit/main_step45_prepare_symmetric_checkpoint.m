function P45=main_step45_prepare_symmetric_checkpoint(O18,varargin)
%MAIN_STEP45_PREPARE_SYMMETRIC_CHECKPOINT
% Reuse ONE previously solved Step-18 theta=0 half-domain case, if present.
% Otherwise, a SINGLE zero-angle Stage-II solve is possible ONLY when
% explicitly requested with 'AllowSolve',true. No EDI in this function.
% The immediately saved checkpoint intentionally contains no K, F, R or
% stress arrays. The original Step-18 data, when available, are not changed.
%
% Examples:
%   P45=main_step45_prepare_symmetric_checkpoint(O18);
%   P45=main_step45_prepare_symmetric_checkpoint(O18,'Npoly',240);
%   % Only if no O18 exists and a new single solve is approved:
%   P45=main_step45_prepare_symmetric_checkpoint([], ...
%       'AllowSolve',true,'Npoly',240);

if nargin<1,O18=[];end
ip=inputParser;
addParameter(ip,'AllowSolve',false,@(x)islogical(x)&&isscalar(x));
addParameter(ip,'Npoly',[],@(x)isempty(x) || ...
    (isnumeric(x)&&isscalar(x)&&isfinite(x)&&x>=32&&x==round(x)));
addParameter(ip,'CheckpointFile', ...
    fullfile(pwd,'step45_symmetric_theta0_solved.mat'), ...
    @(x)ischar(x)||(isstring(x)&&isscalar(x)));
addParameter(ip,'Overwrite',false,@(x)islogical(x)&&isscalar(x));
parse(ip,varargin{:});
opt=ip.Results;
cp=char(opt.CheckpointFile);
if exist(cp,'file')==2 && ~opt.Overwrite
    stored=load(cp,'meta');
    if ~isfield(stored,'meta') || ...
            ~isfield(stored.meta,'caseType') || ...
            ~strcmp(stored.meta.caseType,'step45_centered_half_theta0')
        error('step45:UnrecognizedCheckpoint', ...
            'Existing file is not an identified Step-45 checkpoint.');
    end
    if ~isempty(opt.Npoly) && stored.meta.Npoly~=opt.Npoly
        error('step45:ExistingNpoly', ...
            'Existing Npoly=%d differs from requested Npoly=%d.', ...
            stored.meta.Npoly,opt.Npoly);
    end
    P45=struct('checkpointPath',cp,'reusedFile',true, ...
        'newSolve',false,'meta',stored.meta);
    fprintf('STEP 45: reused existing identified checkpoint: %s\n',cp);
    return
end

legacy=struct('available',false);
if ~isempty(O18)
    for f={'Cases','NpolyList','thetaDeg','KI','KII','rOuterOverA0'}
        need(O18,f{1});
    end
    angleIds=find(abs(O18.thetaDeg(:))<=1e-12);
    if numel(angleIds)~=1
        error('step45:MissingZeroAngle', ...
            'Step18 output must contain exactly one theta=0 case.');
    end
    if isempty(opt.Npoly)
        [~,iMesh]=max(O18.NpolyList(:));
    else
        iMesh=find(O18.NpolyList(:)==opt.Npoly,1);
        if isempty(iMesh)
            error('step45:NpolyAbsent','Requested Npoly absent from O18.');
        end
    end
    i0=angleIds(1);
    caseData=O18.Cases{iMesh,i0};
    if isempty(caseData)
        error('step45:MissingSolvedCase', ...
            'O18.Cases{%d,%d} contains no saved FEM solution.',iMesh,i0);
    end
    for f={'C','Mc','S2'}
        need(caseData,f{1});
    end
    C=caseData.C;
    S2=caseData.S2;
    crack=caseData.Mc.crack;
    source='Step18 saved zero-angle FEM result';
    legacy.available=true;
    legacy.outerRatio=O18.rOuterOverA0(:).';
    legacy.KI=reshape(O18.KI(iMesh,i0,:),1,[]);
    legacy.KII=reshape(O18.KII(iMesh,i0,:),1,[]);
    legacy.ratio=legacy.KII./legacy.KI;
    legacy.note=['Original Step18 EDI used FE-nodal weight, ', ...
        'default SEVEN-point quadrature and variable inner radii. ', ...
        'NOT directly comparable to 16-point Step45 results.'];
    nPoly=C.hole.npoly;
    fprintf('STEP 45: extracting already solved Npoly=%d theta=0 case.\n',nPoly);
else
    if ~opt.AllowSolve
        error('step45:NeedsPriorField', ...
            ['No Step18 output was supplied. Load O18, or explicitly ', ...
             'authorize exactly ONE new zero-angle solve using ', ...
             '''AllowSolve'',true. Nothing was calculated.']);
    end
    if isempty(opt.Npoly),nPoly=240;else,nPoly=opt.Npoly;end
    C=cfg_centered_half_domain();
    C.hole.npoly=nPoly;C.holes={C.hole};
    hArc=2*pi*C.hole.r/nPoly;
    C.mesh1.hmin=hArc;
    C.mesh1.hhole=hArc;
    C.mesh1.hmax=20*hArc;
    C.mesh2.hmax=C.mesh1.hmax;
    C.mesh2.hhole=C.mesh1.hmin;
    C.mesh2.hcrack=C.mesh1.hmin;
    C.solver.verbose=0;
    % This builder overwrites the initiation position with the exact
    % rightmost hole point. A Stage-I solve is NOT required for this
    % independently prescribed zero-angle control geometry.
    [~,~,~,Mc]=build_stage2_centered_half_cracked_mesh_for_theta( ...
        C,struct(),0,'PlotGeom',false,'PlotMesh',false, ...
        'PlotCollapsed',false);
    crack=Mc.crack;
    fprintf('STEP 45: explicitly authorized ONE symmetric Stage-II solve.\n');
    S2=solve_cracked_LEFM(C,Mc,'lambda',1.0);
    source='one explicitly authorized symmetric Stage-II FEM solve';
end

for f={'mesh','U','mat'}
    need(S2,f{1});
end
need(crack,'Pmid');
need(C,'hole');
if ~isfield(C,'domain')|| ...
        ~strcmp(C.domain.mode,'centered_right_half') || ...
        ~strcmp(C.load.type,'remote_tension_y') || ...
        ~strcmp(C.bc.anchor_mode,'symmetry_half_x')
    error('step45:NotSymmetricProblem', ...
        'Requires the centered right-half, symmetry-x, tension-y setup.');
end
a0=norm(diff(crack.Pmid,1,1));
center=C.hole.center(:).';
mouth=[center(1)+C.hole.r,center(2)];
if size(crack.Pmid,1)~=2 || abs(a0-C.a0)>1e-10 || ...
        norm(crack.Pmid(1,:)-mouth)>1e-10 || ...
        abs(crack.Pmid(2,2)-crack.Pmid(1,2))>1e-12 || ...
        crack.Pmid(2,1)<=crack.Pmid(1,1) || ...
        abs(center(2))>1e-12
    error('step45:NotZeroModeIIControl', ...
        'The physical crack/hole geometry is not the symmetric theta=0 control.');
end
mesh=S2.mesh;
U=S2.U;
mat=S2.mat;
if size(mesh.connect,2)~=6 || numel(U)~=2*size(mesh.coord,1) || ...
        any(~isfinite(U))
    error('step45:InvalidField','Expected finite T6 displacement field.');
end
% Avoid checkpointing the original large FEM stiffness and stress arrays.
nT6=size(mesh.coord,1);
nT3=size(mesh.connect3,1);
meta=struct('caseType','step45_centered_half_theta0', ...
    'source',source,'Npoly',nPoly,'a0',a0,'nT6',nT6,'nT3',nT3, ...
    'zeroAngleDeg',0, ...
    'symmetryExpectation','KI nonzero and KII exactly zero in continuum', ...
    'previousStep18EDI',legacy);
clear S2 C Mc O18
[folder,~,~]=fileparts(cp);
if ~isempty(folder)&&exist(folder,'dir')~=7,mkdir(folder);end
tmp=[cp '.incomplete.mat'];
if exist(tmp,'file')==2,delete(tmp);end
save(tmp,'mesh','U','mat','crack','a0','meta','-v7.3');
[ok,msg]=movefile(tmp,cp,'f');
if ~ok,error('step45:SaveFailed','%s',msg);end
P45=struct('checkpointPath',cp,'reusedFile',false, ...
    'newSolve',strcmp(source,'one explicitly authorized symmetric Stage-II FEM solve'), ...
    'meta',meta);
fprintf('  checkpoint: %s\n',cp);
fprintf('  T6 nodes=%d, T3 triangles=%d, crack length=%.8g m\n', ...
    nT6,nT3,a0);
fprintf('  No EDI integration. Saved without FEM stiffness or stress arrays.\n');
end

function need(s,f)
if ~isstruct(s)||~isfield(s,f)||isempty(s.(f))
    error('step45:MissingField','Required field %s is missing.',f);
end
end
