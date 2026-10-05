function R45c=main_step45c_exact_gauss_isolation(R45b,varargin)
%MAIN_STEP45C_EXACT_GAUSS_ISOLATION
% On the EXACT SAME coarse Step45 symmetric control T6 mesh and already
% tested annulus, compare two prescribed pure-I representations:
%   A: exact Williams displacement sampled at NODES then T6-interpolated
%      inside EDI (previously measured; imported from R45b, not rerun).
%   B: exact Williams stress/displacement gradients evaluated DIRECTLY
%      at integration points via SIF_LEFM_interaction_EDI('AnalyticActualK',[1 0]).
%
% Fixed in A and B: mesh, local axes, FE-nodal q, integration domain,
% 16-point quadrature, elasticity and auxiliary Williams convention.
% Thus A-vs-B isolates errors tied to T6 interpolation of the singular
% nodal displacement (plus potential interactions with shared routines).
% Remaining B leakage would implicate domain/q/quadrature implementation.
%
% This test NEVER assembles/solves FEM or modifies the saved checkpoint.
% The result is NOT a correction to physical asymmetric KII.
%
% Usage after the saved Step45b run:
%   R45c=main_step45c_exact_gauss_isolation(R45b);
%   disp(R45c.table);
%
% Start with only the previously tested r_outer/a0=0.65. Further domains
% require explicit opt-in and existing matching R45b rows.
ip=inputParser;
addParameter(ip,'ROuterOverA0',0.65, ...
    @(x)isnumeric(x)&&isvector(x)&&~isempty(x)&& ...
        all(isfinite(x(:)))&&all(x(:)>0)&&all(x(:)<1));
addParameter(ip,'SavePrefix','', ...
    @(x)ischar(x)||(isstring(x)&&isscalar(x)));
parse(ip,varargin{:});
opt=ip.Results;
must(R45b,'table');
must(R45b,'checkpointPath');
T=R45b.table;
required={'r_outer_over_a0','r_inner','actual_KI','actual_KII', ...
    'actual_ratio','exact_KI','exact_KII_leakage', ...
    'exact_ratio_leakage'};
if ~istable(T) || ~all(ismember(required,T.Properties.VariableNames))
    error('step45c:BadStep45bTable', ...
        'R45b must contain measured actual and nodal-exact EDI results.');
end
cp=char(R45b.checkpointPath);
if exist(cp,'file')~=2
    error('step45c:MissingCheckpoint', ...
        'Symmetric solved-field checkpoint not found: %s',cp);
end
info=dir(cp);
key=sprintf('%s|%d|%.16g',cp,info.bytes,info.datenum);
if isempty(opt.SavePrefix)
    prefix=fullfile(fileparts(cp),'step45c_exact_gauss');
else
    prefix=char(opt.SavePrefix);
end
[parent,~,~]=fileparts(prefix);
if ~isempty(parent)&&exist(parent,'dir')~=7,mkdir(parent);end
rad=sort(unique(opt.ROuterOverA0(:).'));
n=numel(rad);
out=nan(n,11);
for i=1:n
    hit=find(abs(T.r_outer_over_a0-rad(i))<1e-12);
    if numel(hit)~=1
        error('step45c:AbsentMatchedDomain', ...
            'Step45b must first measure the EXACT matching domain %.2f.', ...
            rad(i));
    end
    out(i,1:7)=[rad(i),T.r_inner(hit),T.actual_KI(hit), ...
        T.actual_KII(hit),T.actual_ratio(hit), ...
        T.exact_KI(hit),T.exact_KII_leakage(hit)];
    if any(~isfinite(out(i,2:7))) || out(i,2)<=0 || ...
            out(i,3)<=0 || out(i,6)<=0 || ...
            abs(T.exact_ratio_leakage(hit)-out(i,7)/out(i,6))>1e-12
        error('step45c:InvalidMatchedDomain', ...
            'Nonfinite/inconsistent earlier EDI data for radius %.2f.', ...
            rad(i));
    end
end
s=load(cp,'mesh','mat','crack','a0','meta');
for field={'mesh','mat','crack','a0','meta'}
    must(s,field{1});
end
mesh=s.mesh;mat=s.mat;crack=s.crack;a0=s.a0;
if ~isfield(s.meta,'caseType') || ...
        ~strcmp(s.meta.caseType,'step45_centered_half_theta0') || ...
        size(mesh.connect,2)~=6 || ...
        abs(norm(diff(crack.Pmid,1,1))-a0)>1e-12 || ...
        abs(crack.Pmid(end,2)-crack.Pmid(1,2))>1e-12
    error('step45c:BadCheckpoint', ...
        'Expected the exact previously solved horizontal symmetric T6 field.');
end
if ~isfield(mat,'Dmat') && isfield(mat,'D')
    mat.Dmat=mat.D;
end
% The production EDI function checks U size even when AnalyticActualK
% bypasses nodal displacement gradients. This zero vector is a DUMMY:
% it is not a physical FEM field and is never used in actual-field
% evaluation for this diagnostic option.
zeroU=zeros(2*size(mesh.coord,1),1);
cacheFile=[prefix '_progress.mat'];
cache=struct('key',key,'rule',16,'method','fe_nodal', ...
    'rRatio',[],'rInner',[],'gaussKI',[],'gaussKII',[],'done',[]);
if exist(cacheFile,'file')==2
    d=load(cacheFile,'cache');
    prior=d.cache;
    if ~strcmp(prior.key,key) || prior.rule~=16 || ...
            ~strcmp(prior.method,'fe_nodal')
        error('step45c:CacheMismatch', ...
            'Existing cache uses another checkpoint or EDI method.');
    end
    if any(~isfinite(prior.gaussKI(prior.done))) || ...
            any(~isfinite(prior.gaussKII(prior.done)))
        error('step45c:CorruptCache','Cached completed results are invalid.');
    end
    cache=prior;
end
fprintf('\n============================================================\n');
fprintf('STEP 45c: EXACT GAUSS vs NODAL EXACT on SAME COARSE T6\n');
fprintf('============================================================\n');
fprintf('  T6 nodes=%d; a0=%.8g m. No FEM solve, no remesh.\n', ...
    size(mesh.coord,1),a0);
for i=1:n
    ri=out(i,2);
    ro=rad(i)*a0;
    if ~(ro>ri && ri>0)
        error('step45c:BadAnnulus', ...
            'Expected r_outer>r_inner>0 for existing domain.');
    end
    j=find(abs(cache.rRatio-rad(i))<1e-12 & ...
        abs(cache.rInner-ri)<1e-12,1);
    if isempty(j)
        cache.rRatio(end+1)=rad(i);
        cache.rInner(end+1)=ri;
        cache.gaussKI(end+1)=NaN;
        cache.gaussKII(end+1)=NaN;
        cache.done(end+1)=false;
        j=numel(cache.rRatio);
    end
    if ~cache.done(j)
        % Already verified previous Step45b at precisely this domain:
        % DO NOT repeat nodal-exact or actual-FEM EDI calculations.
        [ki,kii,aux]=SIF_LEFM_interaction_EDI( ...
            mesh,zeroU,crack.Pmid,mat, ...
            struct('r_inner',ri,'r_outer',ro), ...
            'UsePlaneStrain',mat.ps==1,'Verbose',false, ...
            'WeightFunction','fe_nodal', ...
            'QuadratureRule',16, ...
            'StoreGPDiagnostics',false, ...
            'AnalyticActualK',[1,0]);
        if ~strcmp(aux.actualFieldSource,'exact_Williams_Gauss_diagnostic') ...
                || ~isfinite(ki) || ~isfinite(kii) || ki<=0
            error('step45c:InvalidGaussReplay', ...
                'Exact Gauss override was not activated or produced bad SIFs.');
        end
        cache.gaussKI(j)=ki;
        cache.gaussKII(j)=kii;
        cache.done(j)=true;
        tmp=[cacheFile '.incomplete.mat'];
        save(tmp,'cache');
        [ok,msg]=movefile(tmp,cacheFile,'f');
        if ~ok,error('step45c:CacheSave','%s',msg);end
        fprintf(['  r_o/a0=%.2f: exact Gauss KI=%.10g, ', ...
            'KII=%+.10g (saved)\n'],rad(i),ki,kii);
    else
        fprintf('  r_o/a0=%.2f exact Gauss reused from cache.\n',rad(i));
    end
    out(i,8)=out(i,7)/out(i,6); % exact nodal ratio
    out(i,9)=cache.gaussKI(j);
    out(i,10)=cache.gaussKII(j);
    out(i,11)=cache.gaussKII(j)/cache.gaussKI(j);
end
report=table(out(:,1),out(:,2),out(:,5),out(:,8),out(:,11), ...
    out(:,3),out(:,6),out(:,9),out(:,4),out(:,7),out(:,10), ...
    'VariableNames',{'r_outer_over_a0','r_inner', ...
    'actual_FEM_ratio','exact_nodal_ratio','exact_Gauss_ratio', ...
    'actual_FEM_KI','exact_nodal_KI','exact_Gauss_KI', ...
    'actual_FEM_KII','exact_nodal_KII','exact_Gauss_KII'});
fprintf('\nMATCHED EDI RECOVERY: SAME MESH, q, ANNULUS, 16-POINT RULE\n');
disp(report);
fprintf(['  Gauss minus nodal exact leakage changes ONLY the way the ', ...
    'prescribed singular displacement is represented.\n']);
fprintf(['  The results cannot be subtracted as a validated correction ', ...
    'to the unknown physical KII.\n']);
R45c=struct('table',report,'checkpointPath',cp, ...
    'method','pure-I prescribed exact GP versus exact sampled T6 nodal', ...
    'usesAnalyticActualK',true, ...
    'noNewFEM',true,'cachePath',cacheFile);
save([prefix '_small_data.mat'],'R45c');
fprintf('STEP 45c complete. No new FEM solution.\n');
end

function must(s,name)
if ~isstruct(s)||~isfield(s,name)||isempty(s.(name))
    error('step45c:MissingInput','Required field %s is missing.',name);
end
end
