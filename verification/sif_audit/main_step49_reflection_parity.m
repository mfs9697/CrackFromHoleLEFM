function R49=main_step49_reflection_parity(varargin)
%MAIN_STEP49_REFLECTION_PARITY
% Inspect geometrical reflection pairing and actual displacement parity on
% the already-SOLVED Step45/Step47 zero-angle symmetric control meshes.
%
% NO FEM solve, NO mesh generation, NO EDI, NO polynomial SIF extrapolation.
% Uses the same local physical disk/annulus and baseline r_inner/r_outer
% recorded by Step48. It diagnoses the observable *source* of asymmetric
% crack-face tangential opening, not an error correction to any SIF.
%
% Physically, under reflection about the horizontal crack y=y0:
%   u_x(x,+y) = u_x(x,-y)  (even), and
%   u_y(x,+y) =-u_y(x,-y)  (odd, up to a rigid-body y translation).
% At coincident crack-face points, an exact symmetric solution therefore
% has tangential jump u_x^upper-u_x^lower = 0.
%
% Distinguish: (i) exactly matched native opposite-face nodes,
% (ii) nearest reflection mates of OFF-FACE T6 nodes, and
% (iii) reflection pairing of complete tip-adjacent T3 triangles.
% Equal upper/lower triangle *counts* are NOT complete reflection pairing.
%
% Default:
%   addpath(genpath(pwd));
%   R49=main_step49_reflection_parity();
%   disp(R49.faceWindows);
%   disp(R49.offFaceMirror);
%   disp(R49.tipFanMirror);
root=fileparts(fileparts(fileparts(mfilename('fullpath'))));
vdir=fullfile(root,'verification');
ip=inputParser;
addParameter(ip,'Step48File',fullfile(vdir, ...
    'step48_refined_matched_edi_comparison_small_data.mat'), ...
    @(s)ischar(s)||(isstring(s)&&isscalar(s)));
addParameter(ip,'SaveFile',fullfile(vdir, ...
    'step49_reflection_parity_small_data.mat'), ...
    @(s)ischar(s)||(isstring(s)&&isscalar(s)));
parse(ip,varargin{:});
opt=ip.Results;
addpath(genpath(root));
if exist(char(opt.Step48File),'file')~=2
    error('step49:MissingComparison', ...
        'Step48 comparison MAT is required for exact matched annulus.');
end
tmp=load(char(opt.Step48File),'R48');
if ~isfield(tmp,'R48'),error('step49:BadComparison','R48 missing.');end
R48=tmp.R48;
need(R48,'comparison');need(R48,'refined');need(R48,'baselineSource');
T=R48.comparison;
if ~istable(T)||height(T)~=2 || ...
        ~all(ismember({'r_inner','r_outer','KI','KII'}, ...
             T.Properties.VariableNames))
    error('step49:MissingDomain','Step48 matched two-field table required.');
end
ri=T.r_inner(1);ro=T.r_outer(1);
if ~all(isfinite([ri ro]))||ri<=0||ro<=ri|| ...
        max(abs(T.r_inner-ri))>1e-12|| ...
        max(abs(T.r_outer-ro))>1e-12
    error('step49:UnmatchedDomains', ...
        'Step48 must have EXACTLY matched baseline/refined EDI domains.');
end
if exist(char(R48.baselineSource),'file')~=2
    error('step49:MissingOriginalResults', ...
        'Original Step45 compact baseline result no longer exists.');
end
original=load(char(R48.baselineSource),'O45');
if ~isfield(original,'O45')
    error('step49:BadOriginalResults','Original O45 variable absent.');
end
o=original.O45;
need(o,'checkpointPath');need(R48.refined,'checkpointPath');
paths={char(o.checkpointPath),char(R48.refined.checkpointPath)};
labels={'Original Step45';'Refined Step47'};
results=cell(2,1);
for i=1:2
    if exist(paths{i},'file')~=2
        error('step49:MissingSolvedCheckpoint', ...
            'Saved FEM checkpoint absent: %s',paths{i});
    end
    s=load(paths{i},'mesh','U','mat','crack','a0','meta');
    for field={'mesh','U','mat','crack','a0','meta'}
        need(s,field{1});
    end
    if ~strcmp(s.meta.caseType,'step45_centered_half_theta0') || ...
            s.meta.Npoly~=240 ||abs(s.a0-0.004)>1e-12 || ...
            size(s.mesh.connect,2)~=6 || ...
            size(s.mesh.connect3,2)~=3 || ...
            numel(s.U)~=2*size(s.mesh.coord,1)
        error('step49:NotMatchedControl', ...
            'Both inputs must be previously solved zero-angle T6 controls.');
    end
    if i==2
        if ~isfield(s.meta,'stage') || ...
                ~strcmp(s.meta.stage,'step47_refined_control') || ...
                s.meta.nT6~=5054
            error('step49:NotSelectedRefinedField', ...
                'The refined checkpoint must be the approved Step47 field.');
        end
    end
    if abs(ro/s.a0-0.65)>1e-12 || ...
            any(~isfinite(s.U)) || ...
            norm(s.crack.Pmid(end,:)-s.crack.Pmid(1,:)-[s.a0,0])>1e-12
        error('step49:ChangedSetup','Not the matched horizontal crack.');
    end
    if i==2
        if norm(s.crack.Pmid-results{1}.Pmid,'fro')>1e-12 || ...
                abs(s.mat.E-results{1}.E)>1e-10 || ...
                abs(s.mat.nu-results{1}.nu)>1e-12 || ...
                s.mat.ps~=results{1}.ps
            error('step49:PhysicalProblemMismatch', ...
                'Solved control crack/material differs between meshes.');
        end
    end
    pairTol=max(1e-12,1e-8*s.a0);
    [fw,native,rmsOpening]=local_face_parity( ...
        s.mesh,s.U,s.mat,s.crack,s.a0,pairTol);
    fan=local_tip_fan(s.mesh,s.crack,pairTol);
    off=cell(2,1);
    off{1}=local_off_face_mates(s.mesh,s.U, ...
        s.crack.Pmid(end,:),0,ro,pairTol,fan.hTip,rmsOpening);
    off{2}=local_off_face_mates(s.mesh,s.U, ...
        s.crack.Pmid(end,:),ri,ro,pairTol,fan.hTip,rmsOpening);
    results{i}=struct('Pmid',s.crack.Pmid,'E',s.mat.E,'nu', ...
        s.mat.nu,'ps',s.mat.ps,'face',fw, ...
        'fan',fan,'off',{off},'nT6',size(s.mesh.coord,1));
    fprintf('\n%s: T6=%d, tip T3 upper/lower=%d/%d\n', ...
        labels{i},size(s.mesh.coord,1),fan.nAbove,fan.nBelow);
    fprintf('  Face symmetry: native upper/lower=%d/%d, mismatch=%.3e m\n', ...
        fw.nUpper,fw.nLower,fw.mismatch);
    fprintf('  Fan mirror worst complete-triangle mismatch=%.6e m\n', ...
        fan.bestMaxMirrorDistance);
    fprintf('  Off-face exact reflection mates in disk: %d/%d upper nodes\n', ...
        off{1}.matched,off{1}.nUpper);
    clear s
end

W=[.04 .30;.08 .30;.12 .30];
fMesh=cell(6,1);winStart=nan(6,1);winEnd=nan(6,1);
nNative=nan(6,1);medianSigned=nan(6,1);medianAbsolute=nan(6,1);
rmsTangential=nan(6,1);rmsOpening=nan(6,1);
residualOddUy=nan(6,1);
row=0;
for i=1:2
    fw=results{i}.face;
    for j=1:3
        row=row+1;
        fMesh{row}=labels{i};
        winStart(row)=W(j,1);winEnd(row)=W(j,2);
        idx=fw.ratioR>=W(j,1)&fw.ratioR<=W(j,2);
        nNative(row)=nnz(idx);
        if nnz(idx)==0,continue;end
        medianSigned(row)=median(fw.rawRatio(idx));
        medianAbsolute(row)=median(abs(fw.rawRatio(idx)));
        rmsTangential(row)=sqrt(mean(fw.jumpX(idx).^2));
        rmsOpening(row)=sqrt(mean(fw.openY(idx).^2));
        % Vertical rigid-body gauge can add a constant to u_y^+ + u_y^-;
        % remove that constant before assessing the ODD component.
        yy=fw.sumY(idx)-median(fw.sumY(idx));
        residualOddUy(row)=sqrt(mean(yy.^2)) /rmsOpening(row);
    end
end
faceWindows=table(fMesh,winStart,winEnd,nNative, ...
    medianSigned,medianAbsolute, ...
    rmsTangential./rmsOpening,residualOddUy, ...
    'VariableNames',{'mesh','r_low_over_a0','r_high_over_a0', ...
    'n_native','median_signed_jumpX_over_openY', ...
    'median_abs_jumpX_over_openY', ...
    'rms_jumpX_over_rms_openY', ...
    'rms_gauge_removed_oddUy_over_openY'});
mirrorMesh=cell(4,1);mirrorRegion=cell(4,1);
nUpper=zeros(4,1);nLower=zeros(4,1);nMatched=zeros(4,1);
coverage=zeros(4,1);medianNearest_over_hTip=nan(4,1);
rms_evenUx_over_openY=nan(4,1);rms_oddUy_over_openY=nan(4,1);
row=0;
for i=1:2
    for region=1:2
        row=row+1;
        v=results{i}.off{region};
        mirrorMesh{row}=labels{i};
        if region==1,mirrorRegion{row}='Tip disk'; ...
        else,mirrorRegion{row}='EDI annulus';end
        nUpper(row)=v.nUpper;
        nLower(row)=v.nLower;
        nMatched(row)=v.matched;
        coverage(row)=v.coverage;
        medianNearest_over_hTip(row)=v.medianNearest_over_hTip;
        rms_evenUx_over_openY(row)=v.rmsEvenUx;
        rms_oddUy_over_openY(row)=v.rmsOddUy;
    end
end
offFaceMirror=table(mirrorMesh,mirrorRegion,nUpper,nLower, ...
    nMatched,coverage,medianNearest_over_hTip, ...
    rms_evenUx_over_openY,rms_oddUy_over_openY, ...
    'VariableNames',{'mesh','region','upper_offFace_T6', ...
    'lower_offFace_T6','exact_reflection_matches', ...
    'fraction_upper_exactly_matched', ...
    'median_nearest_dist_over_hTip', ...
    'rms_evenUx_residual_over_face_openY', ...
    'rms_gauge_removed_oddUy_over_face_openY'});
fanMesh=labels;
fanAbove=[results{1}.fan.nAbove;results{2}.fan.nAbove];
fanBelow=[results{1}.fan.nBelow;results{2}.fan.nBelow];
fanDistance=[results{1}.fan.bestMaxMirrorDistance; ...
    results{2}.fan.bestMaxMirrorDistance];
fanRelative=[results{1}.fan.mirrorMismatchOverHTip; ...
    results{2}.fan.mirrorMismatchOverHTip];
fanExact=[results{1}.fan.exactReflectionPaired; ...
    results{2}.fan.exactReflectionPaired];
tipFanMirror=table(fanMesh,fanAbove,fanBelow,fanDistance, ...
    fanRelative,fanExact, ...
    'VariableNames',{'mesh','tipT3_above','tipT3_below', ...
    'best_worst_reflected_T3_vertex_error_m', ...
    'error_over_hTip','complete_reflection_paired'});
fprintf('\nSTEP 49: FACE PARITY BY FIXED FIT WINDOW (no COD extrapolation)\n');
disp(faceWindows);
fprintf('\nSTEP 49: OFF-FACE REFLECTION PAIRING (unmatched FEM nodes excluded)\n');
disp(offFaceMirror);
fprintf('\nSTEP 49: WHOLE TIP TRIANGLE REFLECTION PAIRING\n');
disp(tipFanMirror);
fprintf(['  Low off-face pairing coverage makes displacement parity ', ...
    'INCONCLUSIVE, not zero.\n']);
fprintf(['  Crack-face jumpX/openY remains observable even on ', ...
    'unpaired off-face meshes.\n']);
R49=struct('faceWindows',faceWindows, ...
    'offFaceMirror',offFaceMirror, ...
    'tipFanMirror',tipFanMirror, ...
    'rInner',ri,'rOuter',ro, ...
    'noNewSolve',true,'noEDI',true, ...
    'note',['Node/triangle pairing and raw native displacements ', ...
       'are diagnostics, not physical KII or its error estimate.']);
saveFile=char(opt.SaveFile);
[folder,~,~]=fileparts(saveFile);
if ~isempty(folder)&&exist(folder,'dir')~=7,mkdir(folder);end
save(saveFile,'R49'); % compact tables ONLY, no FEM displacement checkpoint
fprintf('  Compact Step49 diagnostic saved at %s\n',saveFile);
end

function [F,~,openRMS]=local_face_parity(mesh,U,mat,crack,a0,tol)
% Reuse exact audited native face labeling, then bypass PCHIP since in both
% measured meshes the two native radial grids coincide to roundoff.
[~,~,diag]=native_COD_audit(mesh,U,mat,crack,8,true);
P=mesh.coord;
tip=crack.Pmid(end,:);
xDist=tip(1)-P(:,1);
onFace=xDist>tol & xDist<=a0+tol & ...
    abs(P(:,2)-tip(2))<tol;
up=find(onFace&diag.faceSide==1);
lo=find(onFace&diag.faceSide==-1);
[rU,ju]=sort(xDist(up));[rL,jl]=sort(xDist(lo));
up=up(ju);lo=lo(jl);
if numel(up)~=numel(lo) || ...
        numel(up)~=diag.nUpper || numel(lo)~=diag.nLower || ...
        isempty(up) || max(abs(rU-rL))>tol
    error('step49:UnpairedFaceNodes', ...
        'Exact face-by-face comparison requires both native grids to match.');
end
uu=reshape(U,2,[]).';
jumpX=uu(up,1)-uu(lo,1);
openY=uu(up,2)-uu(lo,2);
sumY=uu(up,2)+uu(lo,2);
if any(~isfinite([jumpX;openY;sumY])) || ...
        any(abs(openY)<1e-14)
    error('step49:InvalidFaceOpenings', ...
        'Cannot normalize native face tangential jump by zero opening.');
end
ratioR=(rU+rL)/(2*a0);
ref=ratioR>=.12 & ratioR<=.30;
if nnz(ref)<2, error('step49:InsufficientReferenceOpening', ...
        'Selected symmetric field has too few reference face nodes.');end
openRMS=sqrt(mean(openY(ref).^2));
F=struct('r_m',0.5*(rU+rL),'ratioR',ratioR, ...
    'nUpper',numel(up),'nLower',numel(lo), ...
    'mismatch',max(abs(rU-rL)),'jumpX',jumpX, ...
    'openY',openY,'sumY',sumY,'rawRatio',jumpX./openY);
end

function G=local_off_face_mates(mesh,U,tip,rin,rout,tol,htip,openRMS)
P=mesh.coord;
r=hypot(P(:,1)-tip(1),P(:,2)-tip(2));
use=r>=rin & r<=rout;
up=find(use & P(:,2)>tip(2)+tol);
lo=find(use & P(:,2)<tip(2)-tol);
d=zeros(numel(up),1);
id=zeros(numel(up),1);
for k=1:numel(up)
    if isempty(lo),d(k)=NaN;continue;end
    mirrored=[P(up(k),1),2*tip(2)-P(up(k),2)];
    dd=hypot(P(lo,1)-mirrored(1),P(lo,2)-mirrored(2));
    [d(k),j]=min(dd);
    id(k)=lo(j);
end
match=isfinite(d)&d<=tol;
n=nnz(match);
if isempty(up),coverage=NaN;else,coverage=n/numel(up);end
nearest=d(isfinite(d));
if isempty(nearest),medN=NaN;else,medN=median(nearest)/htip;end
even=NaN;odd=NaN;
if n>=2
    uu=reshape(U,2,[]).';
    ij=up(match);kl=id(match);
    xerr=uu(ij,1)-uu(kl,1);
    yerr=uu(ij,2)+uu(kl,2);
    yerr=yerr-median(yerr); % remove unknown rigid vertical translation
    even=sqrt(mean(xerr.^2))/openRMS;
    odd=sqrt(mean(yerr.^2))/openRMS;
end
G=struct('nUpper',numel(up),'nLower',numel(lo), ...
    'matched',n,'coverage',coverage, ...
    'medianNearest_over_hTip',medN, ...
    'rmsEvenUx',even,'rmsOddUy',odd);
end

function F=local_tip_fan(mesh,crack,tol)
P=mesh.coord3;T=mesh.connect3;
tip=crack.Pmid(end,:);
r=hypot(P(:,1)-tip(1),P(:,2)-tip(2));
tipIDs=find(r<=min(r)+tol);
near=T(any(ismember(T,tipIDs),2),:);
if isempty(near),error('step49:NoTipFan','No tip-adjacent T3.');end
p1=P(near(:,1),:);p2=P(near(:,2),:);p3=P(near(:,3),:);
L=[hypot(p1(:,1)-p2(:,1),p1(:,2)-p2(:,2)); ...
   hypot(p2(:,1)-p3(:,1),p2(:,2)-p3(:,2)); ...
   hypot(p3(:,1)-p1(:,1),p3(:,2)-p1(:,2))];
hTip=median(L(L>tol));
cent=(p1+p2+p3)/3;
above=near(cent(:,2)>tip(2)+tol,:);
below=near(cent(:,2)<tip(2)-tol,:);
nA=size(above,1);nB=size(below,1);
best=NaN;
if nA>0&&nB>0&&nB<=7 && nA<=nB
    cost=inf(nA,nB);
    vperms=perms(1:3);
    for i=1:nA
        A=P(above(i,:),:);
        A(:,2)=2*tip(2)-A(:,2);
        for j=1:nB
            B=P(below(j,:),:);
            for k=1:size(vperms,1)
                C=B(vperms(k,:),:);
                cost(i,j)=min(cost(i,j), ...
                    max(hypot(A(:,1)-C(:,1),A(:,2)-C(:,2))));
            end
        end
    end
    % Select the complete triangle assignment minimizing the worst
    % reflected-vertex difference; equal counts alone cannot pass.
    assignments=perms(1:nB);
    for i=1:size(assignments,1)
        j=assignments(i,1:nA);
        cur=max(cost(sub2ind(size(cost),1:nA,j)));
        if ~isfinite(best)||cur<best,best=cur;end
    end
end
F=struct('nAbove',nA,'nBelow',nB, ...
    'hTip',hTip,'bestMaxMirrorDistance',best, ...
    'mirrorMismatchOverHTip',best/hTip, ...
    'exactReflectionPaired',nA==nB && ...
          isfinite(best) && best<=tol);
end

function need(s,field)
if ~isstruct(s)||~isfield(s,field)||isempty(s.(field))
    error('step49:MissingField','Missing %s in saved control data.',field);
end
end
