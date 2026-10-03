function O55=main_step55_upper_half_boundary_preflight(varargin)
%MAIN_STEP55_UPPER_HALF_BOUNDARY_PREFLIGHT
% Step55: EXACT prescribed polygon split for a future reflection-paired
% T3 mesh. PURE GEOMETRY ONLY: NO PDE MODEL, generateMesh, FEM solve, EDI.
%
% Starting with the Step54-verified 127-vertex centered right-half polygon:
%   - trace its original UPPER exterior from (xSym,+B) through the
%     original UPPER circular arc and original upper pencil face to tip;
%   - add just ONE artificial internal cut along y=0 from the crack tip
%     to the original right plate boundary (A,0);
%   - close the upper domain along the unchanged original right and top
%     exterior boundaries.
%
% The cut tip -> (A,0) will later be SHARED (intact ligament), while the
% original physical upper/lower pencil faces will be DISTINCT (open crack
% faces) after mesh mirroring and pencil collapse. Those topological
% checks are for the LATER mesh-only stage, not this polygon-only script.
%
% The test reflects all NON-CUT boundary segments of the upper polygon
% and requires a one-to-one undirected full-edge match to the previously
% prescribed entire polygon. Only the original right plate segment is
% allowed to be split by the artificial midpoint (A,0).
%
% Usage:
%   addpath(genpath(pwd));
%   O55=main_step55_upper_half_boundary_preflight();
%   disp(O55.summary);
%
% Produces a COMPACT upper-half polygon and a geometry audit to support
% a separate explicitly inspected mesh-only prototype later.
root=fileparts(fileparts(fileparts(mfilename('fullpath'))));
vdir=fullfile(root,'verification');
ip=inputParser;
addParameter(ip,'Step54File', ...
    fullfile(vdir,'step54_geometry_symmetry_small_data.mat'), ...
    @(x)ischar(x)||(isstring(x)&&isscalar(x)));
addParameter(ip,'Tolerance_m',1e-12, ...
    @(x)isnumeric(x)&&isscalar(x)&&isfinite(x)&&x>0);
addParameter(ip,'SaveFile', ...
    fullfile(vdir,'step55_upper_half_boundary_small_data.mat'), ...
    @(x)ischar(x)||(isstring(x)&&isscalar(x)));
parse(ip,varargin{:});
opt=ip.Results;
addpath(genpath(root));
if exist(char(opt.Step54File),'file')~=2
    error('step55:MissingStep54', ...
        'Run the verified Step54 geometry-only boundary audit first.');
end
a=load(char(opt.Step54File),'O54');
if ~isfield(a,'O54')||~isfield(a.O54,'summary')|| ...
        ~isfield(a.O54,'crackEndpoints')|| ...
        ~isfield(a.O54,'checks')
    error('step55:BadStep54File','Step54 compact result lacks provenance.');
end
S=a.O54;
if height(S.summary)~=1 || ...
        ~logical(S.summary.prescribedPolygonReflectionPass(1)) || ...
        ~logical(S.checks.fullPass) || ...
        S.summary.nOuterPolygonVertices(1)~=127 || ...
        S.summary.nUpperArcVertices(1)~=61 || ...
        S.summary.nLowerArcVertices(1)~=61
    error('step55:UnapprovedPolygon', ...
        'Only the investigator-verified Step54 127-vertex geometry is accepted.');
end
tol=opt.Tolerance_m;
C=cfg_centered_half_domain();
C.hole.npoly=240;
C.holes={C.hole};
C.a0=0.004;
if ~strcmp(C.domain.mode,'centered_right_half') || ...
        abs(C.hole.center(2))>tol || ...
        ~strcmp(C.load.type,'remote_tension_y')
    error('step55:ChangedConfiguration', ...
        'The current prescribed control configuration has changed.');
end
mouth=[C.hole.center(1)+C.hole.r,C.hole.center(2)];
Pmid=[mouth;mouth+[C.a0,0]];
if norm(Pmid-S.crackEndpoints,'fro')>tol
    error('step55:ChangedCrack', ...
        'The crack differs from the Step54-approved endpoints.');
end
D=build_domain_centered_half_pencil(Pmid,C,C.mesh2.chw);
P=D.outerPoly;
if size(P,1)~=127 || size(P,2)~=2 || ...
        any(~isfinite(P(:))) || ...
        norm(D.Pmid-Pmid,'fro')>tol
    error('step55:ChangedPolygon', ...
        'Rebuilt original polygon does not reproduce Step54 geometry.');
end
xSym=C.hole.center(1);
y0=C.hole.center(2);
A=C.A;
B=C.B;
pLeftTop=[xSym,B];
pRightTop=[A,B];
pRightBottom=[A,-B];
pRightMid=[A,y0];
pTip=Pmid(end,:);
iLT=one_vertex(P,pLeftTop,tol);
iRT=one_vertex(P,pRightTop,tol);
iRB=one_vertex(P,pRightBottom,tol);
iTip=one_vertex(P,pTip,tol);
n=size(P,1);
if iLT>=iTip || mod(iRB,n)+1~=iRT || ...
        any(P(iLT:iTip,2)<y0-tol)
    error('step55:UnexpectedPolygonTraversal', ...
        ['Original upper exterior path and right plate segment ', ...
         'are not where this checked polygon-split algorithm expects.']);
end
% Preserve ALL actual upper boundary vertices: no interpolation of the
% circular hole or of the upper sharp-pencil flank takes place.
upperExterior=P(iLT:iTip,:);
Q=[upperExterior; pRightMid; pRightTop];
nQ=size(Q,1);
nNext=[2:nQ,1];
seamIndex=nQ-2;
if norm(Q(seamIndex,:)-pTip)>tol || ...
        norm(Q(seamIndex+1,:)-pRightMid)>tol || ...
        norm(Q(end,:)-pRightTop)>tol || ...
        norm(Q(1,:)-pLeftTop)>tol || ...
        any(Q(:,2)<y0-tol) || ...
        nnz(abs(Q(:,2)-y0)<=tol)~=2
    error('step55:ArtificialSeamInvalid', ...
        'Exactly two upper-half polygon vertices must lie on y=0.');
end
% Boundary corner tip is shared with both regions. Its two adjacent
% upper-half sides are the EXISTING upper pencil face and the NEW seam.
G=D.channelGeom.append;
if norm(Q(seamIndex-1,:)-G.Mup)>tol || ...
        norm(Q(seamIndex,:)-G.xtip)>tol || ...
        norm(G.Mup-[G.Mlo(1),2*y0-G.Mlo(2)])>tol
    error('step55:PencilFlankChanged', ...
        'Upper boundary must retain the original upper pencil face.');
end
signedArea=0.5*sum(Q(:,1).*Q(nNext,2)- ...
                   Q(nNext,1).*Q(:,2));
lengths=hypot(Q(:,1)-Q(nNext,1),Q(:,2)-Q(nNext,2));
if signedArea<=tol^2 || min(lengths)<=tol || ...
        abs(Q(seamIndex,2)-y0)>tol || ...
        abs(Q(seamIndex+1,2)-y0)>tol
    error('step55:DegenerateUpperDomain', ...
        'Split upper domain has incorrect orientation or degenerate edges.');
end
if abs(norm(diff(Q(seamIndex:seamIndex+1,:),1,1)) - ...
       (A-pTip(1)))>tol
    error('step55:IncorrectSeam', ...
        'The artificial cut must be the original intact horizontal ligament.');
end

% Split ONLY the originally unsplit right plate boundary at (A,0).
% This is a geometric edge subdivision, NOT a changed physical boundary.
Psplit=[P(1:iRB,:);pRightMid;P(iRB+1:end,:)];
nFull=size(Psplit,1);
nextFull=[2:nFull,1];
originalEdges=[Psplit,Psplit(nextFull,:)];
keep=(1:nQ)~=seamIndex;
upperEdges=[Q(keep,:),Q(nNext(keep),:)];
lowerEdges=[mirror(upperEdges(:,1:2),y0), ...
            mirror(upperEdges(:,3:4),y0)];
reconstructedEdges=[upperEdges;lowerEdges];
if size(reconstructedEdges,1)~=size(originalEdges,1)
    error('step55:EdgeCountNotReconstructible', ...
        ['Reflecting non-cut upper edges does not reconstruct ', ...
         'the original complete polygon's split boundary edge count.']);
end
% One-to-one UNDIRECTED segment matching, not merely vertex set matching.
used=false(size(reconstructedEdges,1),1);
edgeError=zeros(size(originalEdges,1),1);
for j=1:size(originalEdges,1)
    oldA=originalEdges(j,1:2);
    oldB=originalEdges(j,3:4);
    candA=reconstructedEdges(:,1:2);
    candB=reconstructedEdges(:,3:4);
    forward=max(hypot(candA(:,1)-oldA(1),candA(:,2)-oldA(2)), ...
                hypot(candB(:,1)-oldB(1),candB(:,2)-oldB(2)));
    reverse=max(hypot(candA(:,1)-oldB(1),candA(:,2)-oldB(2)), ...
                hypot(candB(:,1)-oldA(1),candB(:,2)-oldA(2)));
    d=min(forward,reverse);
    d(used)=Inf;
    [edgeError(j),picked]=min(d);
    used(picked)=true;
end
maxEdge=max(edgeError);
unmatched=nnz(edgeError>tol);
boundaryExact=all(used) && unmatched==0 && maxEdge<=tol;
tipMirror=hypot(G.xtip(1)-pTip(1),G.xtip(2)-y0);
seamLen=A-pTip(1);
nAxisUpper=nnz(abs(Q(:,2)-y0)<=tol);
allPass=boundaryExact && tipMirror<=tol && ...
        nAxisUpper==2 && ...
        norm(upperExterior-P(iLT:iTip,:),'fro')<=tol && ...
        norm(Q(seamIndex-1,:)-G.Mup)<=tol;
summary=table(n,nQ,nFull,size(reconstructedEdges,1), ...
    signedArea,min(lengths),seamLen,nAxisUpper, ...
    maxEdge,unmatched,tipMirror,allPass, ...
    'VariableNames',{'nOriginalPolygonVertices', ...
    'nUpperPolygonVertices','nOriginalSplitBoundaryEdges', ...
    'nReconstructedExternalEdges', ...
    'upperDomainArea_m2','smallestUpperEdge_m', ...
    'intactLigamentCutLength_m','nUpperVerticesOnSymmetryAxis', ...
    'maxExteriorEdgeReconstructionError_m', ...
    'unmatchedExteriorEdges','tipOffsetFromAxis_m', ...
    'upperSplitAndReconstructionPass'});
fprintf('\n============================================================\n');
fprintf('STEP 55: EXACT UPPER-HALF POLYGON AND ARTIFICIAL SEAM\n');
fprintf('============================================================\n');
fprintf('  Original source: investigator-verified Step54 geometry.\n');
fprintf('  The original upper circular arc and upper pencil face are UNCHANGED.\n');
fprintf('  The ONLY artificial upper boundary is the intact ligament [tip,A] on y=0.\n');
disp(summary);
if ~allPass
    warning('step55:UpperSplitGateFailed', ...
        'Do not proceed to mesh-only assembly: source polygon was not restored.');
else
    fprintf(['  Reflecting exterior edges reproduces ALL original ', ...
        'boundary segments after splitting the right plate edge.\n']);
end
O55=struct('summary',summary, ...
    'originalOuterPolygon',P, ...
    'upperHalfPolygon',Q, ...
    'artificialSeam',[pTip;pRightMid], ...
    'originalUpperPencilFace',G.face_upper, ...
    'originalLowerPencilFace',G.face_lower, ...
    'originalUpperExterior',upperExterior, ...
    'maxExternalEdgeError_m',maxEdge, ...
    'step54Source',char(opt.Step54File), ...
    'sourceCrack',Pmid, ...
    'note',['Geometry-only split. NO mesh exists yet. On future ', ...
        'mesh assembly, mirror the upper T3 mesh, share axis ligament ', ...
        'nodes from tip to plate right edge, keep two crack faces ', ...
        'distinct and share only the crack tip.'], ...
    'noNewMesh',true,'noNewFEM',true,'noEDI',true);
saveFile=char(opt.SaveFile);
[folder,~,~]=fileparts(saveFile);
if ~isempty(folder)&&exist(folder,'dir')~=7,mkdir(folder);end
save(saveFile,'O55');
fprintf('  Compact upper polygon description: %s\n',saveFile);
fprintf('  No PDE mesh, FEM solve or EDI has been performed.\n');
end

function i=one_vertex(P,q,tol)
d=hypot(P(:,1)-q(1),P(:,2)-q(2));
ids=find(d<=tol);
if numel(ids)~=1
    error('step55:VertexNotUnique', ...
        'Required exact original boundary vertex is missing or duplicated.');
end
i=ids(1);
end

function p=mirror(p,y0)
p(:,2)=2*y0-p(:,2);
end