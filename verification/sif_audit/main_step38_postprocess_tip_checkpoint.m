function O38=main_step38_postprocess_tip_checkpoint(checkpointPath,varargin)
%MAIN_STEP38_POSTPROCESS_TIP_CHECKPOINT
% POSTPROCESS only an existing solved-field Step-38 checkpoint.
% NEVER calls the FEM solver or re-generates meshes.
%
% Uses same three 16-point FE-nodal-q EDI annuli as completed Step 34,
% the EXACT tested native COD implementation, and identical COD
% fit windows/degrees. Saves per-radius, tiny EDI progress after EACH
% integration: if postprocessing is interrupted the completed EDI
% radii are automatically reused on the next call.
%
% The compact result file stores no mesh, U, or stiffness matrix.
%
% O38=main_step38_postprocess_tip_checkpoint(C38.checkpointPath);
% O38=main_step38_postprocess_tip_checkpoint(C38.checkpointPath,...
%      'Plot',true);

ip=inputParser;
addParameter(ip,'Plot',false,@(x)islogical(x)&&isscalar(x));
addParameter(ip,'MinFitPoints',8, ...
    @(x)isnumeric(x)&&isscalar(x)&&isfinite(x)&&x>=6&&x==round(x));
addParameter(ip,'FitWindows',[.04 .20;.04 .30;.08 .30], ...
    @(x)isnumeric(x)&&size(x,2)==2&&all(isfinite(x(:))) && ...
    all(x(:,1)>0)&all(x(:,2)<0.8)&all(x(:,1)<x(:,2)));
addParameter(ip,'FitDegrees',[1 2], ...
    @(x)isnumeric(x)&&isvector(x)&&all(ismember(x,[1 2])));
parse(ip,varargin{:});
opt=ip.Results;
checkpointPath=char(checkpointPath);
if exist(checkpointPath,'file')~=2
    error('step38:CheckpointAbsent','Checkpoint not found: %s', ...
        checkpointPath);
end
infocp=dir(checkpointPath);
key=sprintf('%s|%d|%.16g',checkpointPath, ...
    infocp.bytes,infocp.datenum);
[folder,stem,~]=fileparts(checkpointPath);
progressPath=fullfile(folder,[stem '_edi_progress.mat']);
resultsPath=fullfile(folder,[stem '_results.mat']);
s=load(checkpointPath,'mesh','U','mat','crack','baseline', ...
    'actualTip','a0');
required={'mesh','U','mat','crack','baseline','actualTip','a0'};
for j=1:numel(required)
    if ~isfield(s,required{j})
        error('step38:IncompleteCheckpoint', ...
            'Solved checkpoint is missing %s.',required{j});
    end
end
mesh=s.mesh;U=s.U;mat=s.mat;crack=s.crack;
baseline=s.baseline;
actualTip=s.actualTip;a0=s.a0;
clear s
if numel(U)~=2*size(mesh.coord,1) || ...
        any(~isfinite(U)) || ...
        norm(baseline.crack-crack.Pmid,'fro')>1e-11 || ...
        abs(norm(diff(crack.Pmid,1,1))/a0-1)>1e-10
    error('step38:CheckpointInconsistent', ...
        'Stored mesh/U or physical crack geometry is inconsistent.');
end
if any(~isfinite(baseline.KI)) || ...
        any(~isfinite(baseline.KII)) || baseline.hTip<=0 || ...
        baseline.rInner<=0
    error('step38:InvalidBaseline', ...
        'Reference Step34 SIFs or tip/EDI radii are invalid.');
end
rRat=baseline.rOuterOverA0(:).';
nR=numel(rRat);
if nR~=numel(baseline.KI) || nR~=numel(baseline.KII)
    error('step38:RadiusMismatch','Baseline EDI arrays mismatch.');
end
ri=baseline.rInner;
fprintf('\n============================================================\n');
fprintf('STEP 38 / PHASE 2: CHECKPOINTED EDI AND NATIVE COD\n');
fprintf('============================================================\n');
fprintf('  saved T6=%d; tip size %.10e m (Step34: %.10e m)\n', ...
    size(mesh.coord,1),actualTip,baseline.hTip);
fprintf('  fixed EDI inner radius %.7g m, radii=%s\n', ...
    ri,mat2str(rRat));
fprintf('  No FEM solve, no stiffness matrix, no global remeshing.\n');

KI=nan(1,nR); KII=nan(1,nR);completed=false(1,nR);
if exist(progressPath,'file')==2
    saved=load(progressPath,'progress');
    if ~isfield(saved,'progress') || ...
            ~isstruct(saved.progress) || ...
            ~isfield(saved.progress,'key') || ...
            ~strcmp(saved.progress.key,key) || ...
            ~isequal(saved.progress.rRat,rRat) || ...
            ~isequal(size(saved.progress.KI),size(KI))
        error('step38:StaleProgress', ...
            ['Existing EDI progress belongs to another checkpoint: %s. ', ...
             'Rename it before starting this analysis.'],progressPath);
    end
    KI=saved.progress.KI;
    KII=saved.progress.KII;
    completed=saved.progress.completed;
    clear saved
    fprintf('  REUSING completed EDI domains: %s\n', ...
        mat2str(rRat(completed)));
end

if ~isfield(mat,'D') && ~isfield(mat,'Dmat')
    error('step38:Material','Saved material has no elasticity tensor.');
end
for k=1:nR
    if completed(k)
        fprintf('  EDI radius %.2f: cached, no repeated integration.\n', ...
            rRat(k));
        continue
    end
    [KI(k),KII(k)]=SIF_LEFM_interaction_EDI( ...
        mesh,U,crack.Pmid,mat, ...
        struct('r_inner',ri,'r_outer',rRat(k)*a0), ...
        'UsePlaneStrain',mat.ps==1,'Verbose',false, ...
        'WeightFunction','fe_nodal','QuadratureRule',16, ...
        'StoreGPDiagnostics',false);
    completed(k)=true;
    progress=struct('key',key,'rRat',rRat,'KI',KI,'KII',KII, ...
        'completed',completed);
    % Progress contains just a handful of doubles (not FEM fields).
    tmpProgress=[progressPath '.incomplete.mat'];
    save(tmpProgress,'progress');
    [ok,msg]=movefile(tmpProgress,progressPath,'f');
    if ~ok,error('step38:ProgressSave','%s',msg);end
    fprintf(['  EDI radius %.2f DONE: KI=%.10e KII=%+.10e ', ...
        'KII/KI=%+.10e; progress saved\n'], ...
        rRat(k),KI(k),KII(k),KII(k)/KI(k));
end

[r,app,face]=native_COD_step38(mesh,U,mat,crack,8);
% All expensive FE arrays can now be released BEFORE reporting/plotting.
clear mesh U mat
rr=r/a0;
windows=opt.FitWindows;
degrees=sort(unique(opt.FitDegrees(:).'));
fitRows=nan(size(windows,1)*numel(degrees),14);
n=0;
old=baseline.CODTable;
for iw=1:size(windows,1)
    ids=find(rr>=windows(iw,1) & rr<=windows(iw,2));
    for d=degrees
        if numel(ids)<max(opt.MinFitPoints,2*(d+1))
            fprintf('  COD window [%.2f,%.2f] degree=%d skipped, n=%d\n', ...
                windows(iw,1),windows(iw,2),d,numel(ids));
            continue
        end
        j=find(abs(old.lower_r_over_a0-windows(iw,1))<1e-12 & ...
            abs(old.upper_r_over_a0-windows(iw,2))<1e-12 & ...
            old.degree==d,1);
        if isempty(j)
            error('step38:MissingMatchedBaseline', ...
                'No Step-34 COD fit for [%.2f,%.2f] degree=%d.', ...
                windows(iw,1),windows(iw,2),d);
        end
        pI=polyfit(rr(ids),app(ids,1),d);
        pII=polyfit(rr(ids),app(ids,2),d);
        k1=pI(end);
        k2=pII(end);
        newRatio=k2/k1;
        iRef=find(abs(rRat-.65)==min(abs(rRat-.65)),1);
        ediRatio=KII(iRef)/KI(iRef);
        n=n+1;
        fitRows(n,:)=[windows(iw,:),d, ...
            old.n_native_face(j),numel(ids), ...
            old.refined_KI(j),k1, ...
            old.refined_KII(j),k2, ...
            old.refined_ratio(j),newRatio, ...
            newRatio-old.refined_ratio(j), ...
            100*(ediRatio-newRatio)/ediRatio, ...
            100*(baseline.KII(iRef)/baseline.KI(iRef)- ...
                old.refined_ratio(j))/(baseline.KII(iRef)/baseline.KI(iRef))];
    end
end
CODTable=array2table(fitRows(1:n,:), ...
    'VariableNames',{'lower_r_over_a0','upper_r_over_a0', ...
    'degree','old_n_native','new_n_native', ...
    'old_KI_COD','new_KI_COD','old_KII_COD','new_KII_COD', ...
    'old_ratio_COD','new_ratio_COD','delta_ratio_COD', ...
    'new_EDI_COD_gap_pct','old_EDI_COD_gap_pct'});
qOld=baseline.KII(:)./baseline.KI(:);
qNew=KII(:)./KI(:);
EDIStudy=array2table([rRat(:), ...
    baseline.KI(:),baseline.KII(:),qOld(:), ...
    KI(:),KII(:),qNew(:),qNew(:)-qOld(:)], ...
    'VariableNames',{'r_outer_over_a0','old_KI','old_KII', ...
    'old_ratio','new_KI','new_KII','new_ratio','ratio_change'});

fprintf('\nTIP-REFINEMENT EDI COMPARISON\n');disp(EDIStudy);
fprintf('\nTIP-REFINEMENT NATIVE COD COMPARISON\n');disp(CODTable);
fprintf('  Old EDI ratio domain spread %.10e\n',max(qOld)-min(qOld));
fprintf('  New EDI ratio domain spread %.10e\n',max(qNew)-min(qNew));
fprintf(['  Tip refinement affects the common EDI inner neighborhood. ', ...
    'Outer protected T3 topology remains fixed, but this is not ', ...
    'a proof of physical tiny-kink sign or magnitude.\n']);

O38=struct('checkpointPath',checkpointPath, ...
    'resultsPath',resultsPath, ...
    'baselineTip',baseline.hTip,'refinedTip',actualTip, ...
    'nT6Baseline',baseline.nT6, ...
    'KI',KI,'KII',KII,'baselineKI',baseline.KI, ...
    'baselineKII',baseline.KII,'rOuterOverA0',rRat, ...
    'rInner',ri,'nativeR',r,'nativeApparent',app, ...
    'baselineNativeR',baseline.nativeR, ...
    'baselineNativeApparent',baseline.nativeApparent, ...
    'face',face,'EDIStudy',EDIStudy,'CODTable',CODTable, ...
    'EDIProgressPath',progressPath);

if opt.Plot
    fig=figure('Name','Step38: targeted tip-refinement effect', ...
        'Color','w','Position',[100 90 1250 560]);
    tl=tiledlayout(fig,1,2,'Padding','compact','TileSpacing','compact');
    ax=nexttile(tl);
    hold(ax,'on');grid(ax,'on');
    plot(ax,rRat,qOld,'-o','DisplayName','Step 34 EDI');
    plot(ax,rRat,qNew,'-s','DisplayName','Tip-refined EDI');
    xlabel(ax,'r_{outer}/a_0');
    ylabel(ax,'signed K_{II}/K_I');
    legend(ax,'Location','best');
    ax=nexttile(tl);
    hold(ax,'on');grid(ax,'on');
    xb=baseline.nativeR(:)/a0;
    x=r/a0;
    yb=baseline.nativeApparent;
    plot(ax,xb,yb(:,2)./yb(:,1),'-o','MarkerSize',3, ...
        'DisplayName','Step 34 native COD');
    plot(ax,x,app(:,2)./app(:,1),'-s','MarkerSize',3, ...
        'DisplayName','Tip-refined native COD');
    xlim(ax,[0 .32]);
    yline(ax,qOld(iRef),':','DisplayName','Step 34 EDI');
    yline(ax,qNew(iRef),'--','DisplayName','Tip-refined EDI');
    xlabel(ax,'distance behind tip r/a_0');
    ylabel(ax,'pointwise apparent COD ratio');
    legend(ax,'Location','best');
    pngPath=fullfile(folder,[stem '_comparison.png']);
    exportgraphics(fig,pngPath, ...
        'Resolution',220,'BackgroundColor','white');
    O38.figurePath=pngPath;
end

save(resultsPath,'O38'); % compact: no mesh or FEM U in O38
fprintf('  Compact comparison saved: %s\n',resultsPath);
fprintf('STEP 38 PHASE 2 complete.\n');
end

function [r,app,diag]=native_COD_step38(mesh,U,mat,crack,minPts)
X=mesh.coord;T=mesh.connect;
n=size(X,1);
tip=crack.Pmid(end,:);
vec=crack.Pmid(end,:)-crack.Pmid(1,:);
a0=norm(vec);vec=vec/a0;
perp=[-vec(2),vec(1)];
R=[vec(:),perp(:)];
xl=(X-tip)*R;
face=zeros(n,1);
tipID=crack.tipNode;
up=unique(crack.upperNodes(:));
lo=unique(crack.lowerNodes(:));
face(setdiff(up,[lo;tipID]))=1;
face(setdiff(lo,[up;tipID]))=-1;
emap=[1 2 4;2 3 5;3 1 6];
for j=1:3
    edge=T(:,emap(j,:));
    v1=edge(:,1);v2=edge(:,2);
    mU=(face(v1)==1&(face(v2)==1|v2==tipID)) | ...
       (face(v2)==1&(face(v1)==1|v1==tipID));
    mL=(face(v1)==-1&(face(v2)==-1|v2==tipID)) | ...
       (face(v2)==-1&(face(v1)==-1|v1==tipID));
    idsU=unique(edge(mU,3));idsL=unique(edge(mL,3));
    if any(face(idsU)==-1) || any(face(idsL)==1)
        error('step38:FaceConflict','Crack-face midside sets conflict.');
    end
    face(idsU)=+1;face(idsL)=-1;
end
tol=max(1e-12,1e-8*a0);
onFace=xl(:,1)<-tol & -xl(:,1)<=a0+tol & abs(xl(:,2))<tol;
up=find(onFace&face==1);lo=find(onFace&face==-1);
if numel(up)<minPts||numel(lo)<minPts
    error('step38:FaceNodes','Too few classified crack-face T6 nodes.');
end
[rU,iu]=sort(-xl(up,1));[rL,il]=sort(-xl(lo,1));
up=up(iu);lo=lo(il);
[rU,~,gU]=unique(rU);[rL,~,gL]=unique(rL);
u=reshape(U,2,[]).'*R;
Uu=zeros(numel(rU),2);Ul=zeros(numel(rL),2);
for k=1:2
    Uu(:,k)=accumarray(gU,u(up,k),[],@mean);
    Ul(:,k)=accumarray(gL,u(lo,k),[],@mean);
end
mask=rU>=min(rL)&rU<=max(rL);
r=rU(mask);
if isempty(r),error('step38:FaceOverlap','No face abscissa overlap.');end
jump=Uu(mask,:) - interp1(rL,Ul,r,'pchip');
if any(~isfinite(jump(:)))
    error('step38:NonfiniteCOD','Native COD interpolation failed.');
end
mu=mat.E/(2*(1+mat.nu));
if mat.ps==1
    kappa=3-4*mat.nu;
else
    kappa=(3-mat.nu)/(1+mat.nu);
end
scale=mu/(kappa+1)*sqrt(2*pi./r);
app=bsxfun(@times,[jump(:,2),jump(:,1)],scale);
mismatch=NaN;
if numel(rU)==numel(rL)
    mismatch=max(abs(rU-rL));
end
diag=struct('nUpper',numel(rU),'nLower',numel(rL), ...
    'gridMismatch',mismatch);
fprintf('  COD new mesh: native nodes upper/lower=%d/%d; ', ...
    diag.nUpper,diag.nLower);
fprintf('abscissa mismatch %.5e m\n',mismatch);
end
