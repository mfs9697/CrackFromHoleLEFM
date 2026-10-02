function R52=main_step52_interpolated_reflection(varargin)
%MAIN_STEP52_INTERPOLATED_REFLECTION
% Evaluate ACTUAL saved symmetric FEM fields at physically mirrored spatial
% points via their own T6 shape functions. The previous Step49 found ZERO
% coincident off-face mirror nodes, so direct nodal parity was unevaluable.
%
% Controls on each SAME existing mesh:
%   (A) computed FEM displacement, upper node vs interpolated mirror point;
%   (B) pure-I exact Williams sampled at T6 NODES, same mirror interpolation;
%   (C) exactly affine parity-compatible field, which T6 interpolation MUST
%       reproduce at arbitrary interior points to roundoff;
%   (D) pure-I exact analytical value directly at mirrored sample point.
%
% Only (A) measures actual solved-field parity. (B) quantifies the numerical
% interpolation effect for one prescribed leading-order analytical field;
% it is NOT a correction for (A). (C)/(D) catch sampling or branch errors.
%
% Two fixed physical regions:
%   COD tip neighborhood: r/a0 in [.04,.30];
%   Step48 matched EDI annulus: [r_inner,r_outer].
% A fixed nonzero y-offset excludes points almost on the slit/midline.
%
% NO mesh generation, FEM solve, EDI or SIF polynomial fits. All existing
% solved checkpoints and caches are read-only.
%
% Usage:
%   addpath(genpath(pwd));
%   R52=main_step52_interpolated_reflection();
%   disp(R52.summary);
root=fileparts(fileparts(fileparts(mfilename('fullpath'))));
outDir=fullfile(root,'verification');
ip=inputParser;
addParameter(ip,'Step48File',fullfile(outDir, ...
    'step48_refined_matched_edi_comparison_small_data.mat'), ...
    @(s)ischar(s)||(isstring(s)&&isscalar(s)));
addParameter(ip,'MinAbsYOffset_m',8e-5, ...
    @(s)isnumeric(s)&&isscalar(s)&&isfinite(s)&&s>0);
addParameter(ip,'SaveFile',fullfile(outDir, ...
    'step52_interpolated_reflection_small_data.mat'), ...
    @(s)ischar(s)||(isstring(s)&&isscalar(s)));
parse(ip,varargin{:});
opt=ip.Results;
addpath(genpath(root));
inputFile=char(opt.Step48File);
if exist(inputFile,'file')~=2
    error('step52:Step48Absent', ...
        'Missing existing matched Step48 result: %s',inputFile);
end
a=load(inputFile,'R48');
if ~isfield(a,'R48'),error('step52:Step48Format','Missing R48.');end
R48=a.R48;
if ~isfield(R48,'comparison')||~isfield(R48,'refined')|| ...
        ~isfield(R48,'baselineSource')
    error('step52:Step48Format','Incomplete saved Step48 result.');
end
T=R48.comparison;
if ~istable(T)||height(T)~=2 || ...
        ~all(ismember({'r_inner','r_outer'}, ...
                     T.Properties.VariableNames))
    error('step52:Step48Domains','Step48 comparison table missing.');
end
ri=T.r_inner(1);ro=T.r_outer(1);
if ~isfinite(ri)||~isfinite(ro)||ri<=0||ri>=ro|| ...
        max(abs(T.r_inner-ri))>1e-12|| ...
        max(abs(T.r_outer-ro))>1e-12
    error('step52:UnmatchedDomain','Original/refined annuli differ.');
end
if exist(char(R48.baselineSource),'file')~=2
    error('step52:BaselineResultAbsent', ...
        'Original Step45 compact EDI file missing.');
end
o=load(char(R48.baselineSource),'O45');
if ~isfield(o,'O45')|| ...
        ~isfield(o.O45,'checkpointPath')|| ...
        ~isfield(R48.refined,'checkpointPath')
    error('step52:CheckpointReference','Missing saved FEM checkpoint path.');
end
paths={char(o.O45.checkpointPath), ...
       char(R48.refined.checkpointPath)};
names={'Original Step45','Refined Step47'};
regionNames={'COD tip disk','EDI annulus'};
rows=cell(4,1);
k=0;
savedPhysical=struct();
fprintf('\n============================================================\n');
fprintf('STEP 52: T6-INTERPOLATED REFLECTED DISPLACEMENT PARITY\n');
fprintf('============================================================\n');
fprintf('  Two ALREADY SOLVED control fields, no FEM solve, mesh, EDI or fit.\n');
fprintf('  Matched EDI annulus: r_inner=%.12g, r_outer=%.12g m.\n',ri,ro);
fprintf('  Fixed |y| clearance: %.9g m.\n',opt.MinAbsYOffset_m);
for imesh=1:2
    cp=paths{imesh};
    if exist(cp,'file')~=2
        error('step52:MissingSolvedCheckpoint', ...
            'Previously saved FEM field is missing: %s',cp);
    end
    s=load(cp,'mesh','U','mat','crack','a0','meta');
    for field={'mesh','U','mat','crack','a0','meta'}
        need(s,field{1});
    end
    if ~isfield(s.meta,'caseType') || ...
            ~strcmp(s.meta.caseType,'step45_centered_half_theta0') || ...
            s.meta.Npoly~=240 || ...
            abs(s.a0-0.004)>1e-12 || ...
            abs(ro/s.a0-.65)>1e-12 || ...
            size(s.mesh.connect,2)~=6 || ...
            size(s.mesh.connect3,2)~=3 || ...
            size(s.mesh.connect,1)~=size(s.mesh.connect3,1) || ...
            ~isequal(s.mesh.connect(:,1:3),s.mesh.connect3) || ...
            any(~isfinite(s.U)) || ...
            numel(s.U)~=2*size(s.mesh.coord,1)
        error('step52:WrongFEMCheckpoint', ...
            'Expected unchanged saved horizontal Step45/47 T6 fields.');
    end
    if imesh==2
        if ~isfield(s.meta,'stage') || ...
                ~strcmp(s.meta.stage,'step47_refined_control') || ...
                s.meta.nT6~=5054 || ...
                norm(s.crack.Pmid-savedPhysical.Pmid,'fro')>1e-12 || ...
                abs(s.mat.E-savedPhysical.E)>1e-10 || ...
                abs(s.mat.nu-savedPhysical.nu)>1e-12 || ...
                s.mat.ps~=savedPhysical.ps
            error('step52:NotSamePhysicalProblem', ...
                'Refined checkpoint does not match original physical control.');
        end
    else
        savedPhysical=struct('Pmid',s.crack.Pmid, ...
            'E',s.mat.E,'nu',s.mat.nu,'ps',s.mat.ps);
    end
    tip=s.crack.Pmid(end,:);
    if norm(diff(s.crack.Pmid,1,1)-[s.a0 0])>1e-12 || ...
            abs(tip(2))>1e-12
        error('step52:CrackNotHorizontal', ...
            'Reflection test is restricted to the known horizontal crack.');
    end

    % Reuse the previously tested face labeling and exact-field generator.
    [faceR,apparent,face]=native_COD_audit( ...
        s.mesh,s.U,s.mat,s.crack,8,true);
    X=s.mesh.coord;
    Xlocal=X-tip;
    ids=find(face.faceSide~=0);
    Xlocal(ids,2)=0; % crack-face evaluation only, no mesh modification
    upperFace=find(face.faceSide==1);
    lowerFace=find(face.faceSide==-1);
    Usyn=exact_williams_displacement_audit( ...
        Xlocal,1,0,s.mat.E,s.mat.nu,s.mat.ps, ...
        'UpperFaceIDs',upperFace,'LowerFaceIDs',lowerFace);
    [rSyn,appSyn]=native_COD_audit( ...
        s.mesh,Usyn,s.mat,s.crack,8);
    if numel(faceR)~=numel(rSyn) || ...
            max(abs(faceR-rSyn))>1e-12 || ...
            max(abs(appSyn(:,1)-1))>1e-8 || ...
            max(abs(appSyn(:,2)))>1e-8
        error('step52:ExactCODSelfCheck', ...
            'Known pure-I nodal field did not pass native COD validation.');
    end
    mu=s.mat.E/(2*(1+s.mat.nu));
    if s.mat.ps==1
        kappa=3-4*s.mat.nu;
    else
        kappa=(3-s.mat.nu)/(1+s.mat.nu);
    end
    scale=mu/(kappa+1)*sqrt(2*pi./faceR);
    ref=faceR/s.a0>=.12 & faceR/s.a0<=.30;
    if nnz(ref)<2
        error('step52:MissingOpeningReference', ...
            'Insufficient native crack-face opening reference points.');
    end
    normActual=sqrt(mean((apparent(ref,1)./scale(ref)).^2));
    normSyn=sqrt(mean((appSyn(ref,1)./scale(ref)).^2));
    if ~all(isfinite([normActual,normSyn])) || ...
            min([normActual,normSyn])<=0
        error('step52:InvalidNormalOpening','Bad displacement normalization.');
    end
    Uactual=reshape(s.U,2,[]).';
    Uexact=reshape(Usyn,2,[]).';
    Usmooth=[X(:,1)-tip(1),X(:,2)-tip(2)]; % affine even/odd
    TR=triangulation(s.mesh.connect3,s.mesh.coord3);
    r=hypot(X(:,1)-tip(1),X(:,2)-tip(2));
    for region=1:2
        if region==1
            inRegion=r>=.04*s.a0 & r<=.30*s.a0;
        else
            inRegion=r>=ri & r<=ro;
        end
        chosen=find(inRegion & ...
            X(:,2)>tip(2)+opt.MinAbsYOffset_m & ...
            face.faceSide==0);
        if isempty(chosen)
            Q=nan(0,2);
        else
            Q=[X(chosen,1),2*tip(2)-X(chosen,2)];
        end
        if isempty(Q)
            elem=zeros(0,1);
            B=zeros(0,3);
        else
            [elem,B]=pointLocation(TR,Q);
        end
        hit=isfinite(elem);
        nHit=nnz(hit);
        if nHit>=2
            idx=chosen(hit);
            elem=elem(hit);
            bary=B(hit,:);
            if any(~isfinite(bary(:))) || ...
                    any(abs(sum(bary,2)-1)>1e-7)
                error('step52:BadPointLocation', ...
                    'Barycentric values inconsistent with T3 geometry.');
            end
            v1=bary(:,1);v2=bary(:,2);v3=bary(:,3);
            N=[v1.*(2*v1-1),v2.*(2*v2-1),v3.*(2*v3-1), ...
                4*v1.*v2,4*v2.*v3,4*v3.*v1];
            C=s.mesh.connect(elem,:);
            actualMirror=zeros(nHit,2);
            exactMirror=zeros(nHit,2);
            smoothMirror=zeros(nHit,2);
            for j=1:nHit
                actualMirror(j,:)=N(j,:)*Uactual(C(j,:),:);
                exactMirror(j,:)=N(j,:)*Uexact(C(j,:),:);
                smoothMirror(j,:)=N(j,:)*Usmooth(C(j,:),:);
            end
            % Affine parity of this mesh's T6 spatial evaluator MUST
            % hold to roundoff irrespective of unequal mesh geometry.
            aff=[Usmooth(idx,1)-smoothMirror(:,1), ...
                 Usmooth(idx,2)+smoothMirror(:,2)];
            maxAffine=max(abs(aff(:)));
            if maxAffine>1e-10*s.a0
                error('step52:ShapeFunctionValidation', ...
                    'Affine parity failed; check T6 interpolation ordering.');
            end

            % Pure-I true values directly at reflected spatial locations:
            % distinguishes interpolation artifacts from branch-cut errors.
            direct=exact_williams_displacement_audit( ...
                Q(hit,:)-tip,1,0,s.mat.E,s.mat.nu,s.mat.ps);
            direct=reshape(direct,2,[]).';
            exactDirectErr=[Uexact(idx,1)-direct(:,1), ...
                            Uexact(idx,2)+direct(:,2)];
            maxExactDirect=max(abs(exactDirectErr(:)));
            if maxExactDirect>1e-10*normSyn
                error('step52:ExactMirrorMismatch', ...
                    'Direct analytical field not mirror symmetric.');
            end

            ae=Uactual(idx,1)-actualMirror(:,1);
            ao=Uactual(idx,2)+actualMirror(:,2);
            ao=ao-median(ao); % rigid translation y-gauge removal
            ee=Uexact(idx,1)-exactMirror(:,1);
            eo=Uexact(idx,2)+exactMirror(:,2);
            eo=eo-median(eo); % identical gauge convention
            row=struct( ...
                'mesh',names{imesh},'region',regionNames{region}, ...
                'nUpperCandidates',numel(chosen), ...
                'nMirroredQueriesLocated',nHit, ...
                'queryCoverage',nHit/max(1,numel(chosen)), ...
                'maxAffineParityError_m',maxAffine, ...
                'maxDirectPureIParityError',maxExactDirect, ...
                'rmsActualEvenUx_over_openY', ...
                     sqrt(mean(ae.^2))/normActual, ...
                'rmsActualGaugeOddUy_over_openY', ...
                     sqrt(mean(ao.^2))/normActual, ...
                'rmsNodalPureIEvenUx_over_openY', ...
                     sqrt(mean(ee.^2))/normSyn, ...
                'rmsNodalPureIGaugeOddUy_over_openY', ...
                     sqrt(mean(eo.^2))/normSyn);
        else
            row=struct( ...
                'mesh',names{imesh},'region',regionNames{region}, ...
                'nUpperCandidates',numel(chosen), ...
                'nMirroredQueriesLocated',nHit, ...
                'queryCoverage',nHit/max(1,numel(chosen)), ...
                'maxAffineParityError_m',NaN, ...
                'maxDirectPureIParityError',NaN, ...
                'rmsActualEvenUx_over_openY',NaN, ...
                'rmsActualGaugeOddUy_over_openY',NaN, ...
                'rmsNodalPureIEvenUx_over_openY',NaN, ...
                'rmsNodalPureIGaugeOddUy_over_openY',NaN);
        end
        k=k+1;
        rows{k}=row;
    end
    clear s Usyn
end
R52=struct();
R52.summary=struct2table(vertcat(rows{:}));
R52.source=paths;
R52.rInner=ri;R52.rOuter=ro;
R52.minAbsYOffset_m=opt.MinAbsYOffset_m;
R52.numericalScope=['A point-location/T6 interpolation test on two ', ...
    'previously solved FEM fields, not a physical KII estimate.'];
R52.noNewFEM=true;
R52.noEDI=true;
fprintf('\nSTEP 52: REFLECTED-POINT T6 PARITY VS PRESCRIBED PURE-I\n');
disp(R52.summary);
fprintf(['  Genuine FEM and prescribed leading Williams parity are ', ...
    'different estimands. Do not subtract them as an error correction.\n']);
file=char(opt.SaveFile);
[folder,~,~]=fileparts(file);
if ~isempty(folder)&&exist(folder,'dir')~=7,mkdir(folder);end
save(file,'R52'); % compact summary only, no FEM mesh/displacements
fprintf('  Compact output: %s\n',file);
fprintf('  NO FEM solve, mesh or EDI calculation performed.\n');
end

function need(s,f)
if ~isstruct(s)||~isfield(s,f)||isempty(s.(f))
    error('step52:MissingSavedField', ...
        'Previously saved checkpoint is missing %s.',f);
end
end
