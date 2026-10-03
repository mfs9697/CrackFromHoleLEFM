function Out=main_step45_coarse_exact_pureI_replay(O45,varargin)
%MAIN_STEP45_COARSE_EXACT_PUREI_REPLAY
% Replay one prescribed EXACT pure Mode-I leading Williams displacement
% field on the SAVED, coarse, genuinely symmetric Step-45 FEM T6 mesh.
% Compare exact-field EDI pure-I->II leakage against the already measured
% actual-FEM EDI leakage, at EXACTLY the same 16-point FE-nodal-q annulus.
%
% This is Step 45b: NO FEM solve, NO remeshing, NO additional physical
% field EDI. Default: ONLY the previously computed rOuter/a0=0.65 domain.
% Requires O45 from main_step45_symmetric_field_leakage(...,'RunEDI',true).
% Use saved compact file if O45 was cleared from MATLAB workspace:
%   load('verification/step45_symmetric_field_leakage_small_data.mat','O45')
%
% Usage:
%   Out=main_step45_coarse_exact_pureI_replay(O45);
%   disp(Out.table)
% Exact displacement generation uses the SAME crack-face branch-cut labels
% and coordinate rotations as Step44. COD exact-field self-check runs first.
% This checks exact-field EDI interpolation/integration on THIS coarse mesh;
% it is not an independently solved physical equilibrium benchmark.
ip=inputParser;
addParameter(ip,'ROuterOverA0',0.65, ...
    @(x)isnumeric(x)&&isvector(x)&&~isempty(x) && ...
        all(isfinite(x))&&all(x>0)&&all(x<1));
addParameter(ip,'SavePrefix','', ...
    @(x)ischar(x)||(isstring(x)&&isscalar(x)));
parse(ip,varargin{:});
opt=ip.Results;
must(O45,'checkpointPath');
must(O45,'EDI');
if ~isfield(O45.EDI,'performed')||~O45.EDI.performed|| ...
        ~isfield(O45.EDI,'table')
    error('step45b:PhysicalEDIRequired', ...
        'Run actual FEM EDI first, using the saved Step45 field.');
end
cp=char(O45.checkpointPath);
if exist(cp,'file')~=2
    error('step45b:CheckpointMissing','Missing FEM checkpoint %s.',cp);
end
if ~isfield(O45,'meta') || ...
        ~strcmp(O45.meta.caseType,'step45_centered_half_theta0')
    error('step45b:NotSymmetricControl', ...
        'O45 must refer to the centered theta=0 FEM control.');
end
f=dir(cp);
key=sprintf('%s|%d|%.16g',cp,f.bytes,f.datenum);
if isempty(opt.SavePrefix)
    prefix=fullfile(fileparts(cp),'step45_coarse_exact_pureI');
else
    prefix=char(opt.SavePrefix);
end
[folder,~,~]=fileparts(prefix);
if ~isempty(folder)&&exist(folder,'dir')~=7,mkdir(folder);end
rat=sort(unique(opt.ROuterOverA0(:).'));
Tactual=O45.EDI.table;
rows=nan(numel(rat),7);
for k=1:numel(rat)
    ia=find(abs(Tactual.r_outer_over_a0-rat(k))<1e-12);
    if numel(ia)~=1
        error('step45b:MissingActualDomain', ...
            'The actual-FEM O45.EDI.table needs domain %.2f.',rat(k));
    end
    rows(k,1:4)=[rat(k),Tactual.r_inner(ia), ...
        Tactual.KI(ia),Tactual.KII(ia)];
    if ~isfinite(rows(k,2))||~isfinite(rows(k,3))|| ...
            rows(k,3)<=0 || ~isfinite(rows(k,4))
        error('step45b:InvalidActualData','Invalid physical EDI data.');
    end
end
s=load(cp,'mesh','mat','crack','a0','meta');
for fName={'mesh','mat','crack','a0','meta'}
    must(s,fName{1});
end
mesh=s.mesh;mat=s.mat;crack=s.crack;a0=s.a0;
if ~strcmp(s.meta.caseType,O45.meta.caseType) || ...
        s.meta.Npoly~=O45.meta.Npoly || ...
        abs(s.meta.a0-a0)>1e-12
    error('step45b:CheckpointMetadataMismatch', ...
        'Reported O45 and saved symmetric checkpoint differ.');
end
if size(mesh.connect,2)~=6 || ...
        abs(norm(diff(crack.Pmid,1,1))-a0)>1e-12 || ...
        abs(crack.Pmid(end,2)-crack.Pmid(1,2))>1e-12
    error('step45b:GeometryMismatch','Not a straight zero-angle T6 mesh.');
end

fprintf('\n============================================================\n');
fprintf('STEP 45b: EXACT PURE-I ON SAME COARSE SYMMETRIC T6 MESH\n');
fprintf('============================================================\n');
n=size(mesh.coord,1);
fprintf('  T6 nodes=%d; a0=%.7g m; no FEM solve.\n',n,a0);
tip=crack.Pmid(end,:);
ex=(crack.Pmid(end,:)-crack.Pmid(1,:))/a0;
ey=[-ex(2),ex(1)];
R=[ex(:),ey(:)];
xl=(mesh.coord-tip)*R;
% Reuse the exact Step44 face-topology classification rather than
% infer upper/lower identities from coincident node coordinates.
[~,~,fd]=native_COD_audit(mesh,zeros(2*n,1),mat,crack,8,true);
face=fd.faceSide;
upper=find(face==+1);lower=find(face==-1);
faceIDs=find(face~=0);
faceTol=max(1e-12,1e-8*a0);
if isempty(faceIDs)||...
        any(abs(xl(faceIDs,2))>faceTol) || ...
        any(xl(faceIDs,1)>=-faceTol)
    error('step45b:FaceGeometry', ...
        'Crack-face labels are inconsistent with local crack geometry.');
end
xEval=xl;
xEval(faceIDs,2)=0; % evaluation-only snap of tiny face roundoff
uLocal=exact_williams_displacement_audit(xEval,1,0, ...
    mat.E,mat.nu,mat.ps, ...
    'UpperFaceIDs',upper,'LowerFaceIDs',lower);
clear xEval
uGlobal=reshape(uLocal,2,[]).'*R.';
clear uLocal
Uexact=reshape(uGlobal.',[],1);
clear uGlobal
[r,p]=native_COD_audit(mesh,Uexact,mat,crack,8);
codKI=max(abs(p(:,1)-1));
codKII=max(abs(p(:,2)));
fprintf(['  EXACT COD self-check n=%d, max |KI-1|=%.3e, ', ...
    'max |spurious KII|=%.3e.\n'],numel(r),codKI,codKII);
if codKI>1e-8 || codKII>1e-8
    error('step45b:CODSelfCheck', ...
        'Exact-field face assignment/normalization did not pass COD.');
end
clear p r
if ~isfield(mat,'Dmat') && isfield(mat,'D')
    mat.Dmat=mat.D;
end

progressFile=[prefix '_edi_progress.mat'];
progress=struct('key',key,'rRatio',[], ...
    'rInner',[],'rule',16,'method','fe_nodal', ...
    'exactKI',[],'exactKII',[],'done',[]);
if exist(progressFile,'file')==2
    prev=load(progressFile,'progress');
    saved=prev.progress;
    if ~strcmp(saved.key,key)||saved.rule~=16||...
            ~strcmp(saved.method,'fe_nodal')
        error('step45b:StaleCache', ...
            'Exact EDI progress belongs to a different checkpoint/method.');
    end
    progress=saved;
end
for k=1:numel(rat)
    ro=rat(k)*a0;
    ri=rows(k,2); % EXACT SAME DOMAIN as already measured physical EDI
    if ~(ri>0 && ro>ri)
        error('step45b:BadDomain','Invalid previously measured EDI annulus.');
    end
    j=find(abs(progress.rRatio-rat(k))<1e-12 & ...
        abs(progress.rInner-ri)<1e-12,1);
    if isempty(j)
        progress.rRatio(end+1)=rat(k);
        progress.rInner(end+1)=ri;
        progress.exactKI(end+1)=NaN;
        progress.exactKII(end+1)=NaN;
        progress.done(end+1)=false;
        j=numel(progress.rRatio);
    end
    if ~progress.done(j)
        [ki,kii]=SIF_LEFM_interaction_EDI( ...
            mesh,Uexact,crack.Pmid,mat, ...
            struct('r_inner',ri,'r_outer',ro), ...
            'UsePlaneStrain',mat.ps==1,'Verbose',false, ...
            'WeightFunction','fe_nodal','QuadratureRule',16, ...
            'StoreGPDiagnostics',false);
        if ~isfinite(ki)||~isfinite(kii)||ki<=0
            error('step45b:InvalidExactEDI','Exact pure-I EDI invalid.');
        end
        progress.exactKI(j)=ki;
        progress.exactKII(j)=kii;
        progress.done(j)=true;
        tmp=[progressFile '.incomplete.mat'];
        save(tmp,'progress');
        [ok,msg]=movefile(tmp,progressFile,'f');
        if ~ok,error('step45b:SaveCache','%s',msg);end
        fprintf('  exact EDI %.2f: KI=%.10g, KII=%+.10g\n', ...
            rat(k),ki,kii);
    else
        fprintf('  exact EDI %.2f reused from matching cache.\n',rat(k));
    end
    rows(k,5:7)=[progress.exactKI(j),progress.exactKII(j), ...
        progress.exactKII(j)/progress.exactKI(j)];
end
T=table(rows(:,1),rows(:,2),rows(:,3),rows(:,4), ...
    rows(:,4)./rows(:,3),rows(:,5),rows(:,6),rows(:,7), ...
    'VariableNames',{'r_outer_over_a0','r_inner', ...
    'actual_KI','actual_KII','actual_ratio', ...
    'exact_KI','exact_KII_leakage','exact_ratio_leakage'});
fprintf('\nACTUAL FEM VERSUS EXACT PURE-I ON IDENTICAL T6 AND ANNULUS\n');
disp(T);
fprintf(['  Compare SIGNED exact leakage with ACTUAL-FEM symmetry residual. ', ...
    'Do not treat their difference as a validated numerical correction.\n']);
Out=struct('table',T,'checkpointPath',cp, ...
    'CODSelfCheck',struct('maxAbsKIError',codKI, ...
        'maxAbsKIILeakage',codKII), ...
    'cachePath',progressFile, ...
    'exactField','prescribed leading Williams, KI=1, KII=0', ...
    'physicalFieldUnchanged',true);
save([prefix '_small_data.mat'],'Out');
fprintf('STEP 45b completed: no FEM solution or geometry changes.\n');
end

function must(s,f)
if ~isstruct(s)||~isfield(s,f)||isempty(s.(f))
    error('step45b:MissingField','Required field %s is missing.',f);
end
end
