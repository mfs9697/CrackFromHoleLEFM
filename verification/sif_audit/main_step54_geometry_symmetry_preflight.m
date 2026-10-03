function O54=main_step54_geometry_symmetry_preflight(varargin)
%MAIN_STEP54_GEOMETRY_SYMMETRY_PREFLIGHT
% Rebuild ONLY the Step45/47 prescribed centered, horizontal-crack polygon
% and audit its y-reflection at complete VERTEX and BOUNDARY EDGE levels.
% NO generateMesh, PDE model creation, FEM solve, EDI or COD operation.
%
% This is a prerequisite to a future deliberately mirror-paired FEM mesh:
% an asymmetric input polygon must NOT silently be replaced with a
% physically different symmetric polygon during mesh construction.
%
% Reconstructing this script's prescribed polygon does NOT demonstrate
% that every hidden PDE preprocessor preserves geometric reflection.
%
% Usage:
%  addpath(genpath(pwd));
%  O54=main_step54_geometry_symmetry_preflight();
%  disp(O54.summary);
root=fileparts(fileparts(fileparts(mfilename('fullpath'))));
vdir=fullfile(root,'verification');
ip=inputParser;
addParameter(ip,'OriginalCheckpoint', ...
    fullfile(vdir,'step45_symmetric_theta0_solved.mat'), ...
    @(x)ischar(x)||(isstring(x)&&isscalar(x)));
addParameter(ip,'RefinedCheckpoint', ...
    fullfile(vdir,'step47_refined_symmetric_theta0_solved.mat'), ...
    @(x)ischar(x)||(isstring(x)&&isscalar(x)));
addParameter(ip,'Tolerance_m',1e-12, ...
    @(x)isnumeric(x)&&isscalar(x)&&isfinite(x)&&x>0);
addParameter(ip,'SaveFile', ...
    fullfile(vdir,'step54_geometry_symmetry_small_data.mat'), ...
    @(x)ischar(x)||(isstring(x)&&isscalar(x)));
parse(ip,varargin{:});
opt=ip.Results;
addpath(genpath(root));
oldFile=char(opt.OriginalCheckpoint);
newFile=char(opt.RefinedCheckpoint);
if exist(oldFile,'file')~=2||exist(newFile,'file')~=2
    error('step54:MissingSavedCheckpoint', ...
        'Both previously saved Step45/Step47 checkpoints are required.');
end
a=load(oldFile,'crack','a0','meta');
b=load(newFile,'crack','a0','meta');
for x={a,b}
    z=x{1};
    if ~isfield(z,'crack')||~isfield(z,'a0')|| ...
            ~isfield(z,'meta')|| ...
            ~strcmp(z.meta.caseType,'step45_centered_half_theta0') || ...
            z.meta.Npoly~=240 || ...
            abs(z.a0-.004)>1e-12
        error('step54:WrongSavedCheckpoint', ...
            'Expected two previously solved Npoly=240, a0=0.004 controls.');
    end
end
if ~isfield(b.meta,'stage') || ...
        ~strcmp(b.meta.stage,'step47_refined_control') || ...
        norm(a.crack.Pmid-b.crack.Pmid,'fro')>1e-12
    error('step54:DifferentPhysicalControls', ...
        'Solved control checkpoints have different crack coordinates.');
end
C=cfg_centered_half_domain();
C.hole.npoly=240;
C.holes={C.hole};
C.a0=a.a0;
if abs(C.hole.center(2))>opt.Tolerance_m || ...
        ~strcmp(C.domain.mode,'centered_right_half') || ...
        ~strcmp(C.load.type,'remote_tension_y') || ...
        ~strcmp(C.bc.anchor_mode,'symmetry_half_x')
    error('step54:NotCenteredControl', ...
        'Current configuration is not centered-hole half-domain control.');
end
mouth=[C.hole.center(1)+C.hole.r,C.hole.center(2)];
Pmid=[mouth;mouth+[C.a0,0]];
if norm(Pmid-a.crack.Pmid,'fro')>opt.Tolerance_m
    error('step54:RebuiltCrackDiffers', ...
        'Reconstructed crack endpoints do not match saved FEM controls.');
end
% Direct construction: this builder returns a polygon description ONLY.
% It does NOT call generateMesh or solve any FEM equations.
D=build_domain_centered_half_pencil(Pmid,C,C.mesh2.chw);
P=D.outerPoly;
if size(P,2)~=2 || size(P,1)<12 || ...
        any(~isfinite(P(:))) || ...
        norm(D.Pmid-Pmid,'fro')>opt.Tolerance_m
    error('step54:InvalidPolygon','Prescribed polygon was not recovered.');
end
tol=opt.Tolerance_m;
n=size(P,1);
jNext=[2:n,1];
edgeLength=hypot(P(:,1)-P(jNext,1),P(:,2)-P(jNext,2));
if min(edgeLength)<=tol
    error('step54:DegeneratePolygon', ...
        'Prescribed polygon contains a degenerate boundary edge.');
end
% Mirror each boundary vertex (x,y) -> (x, 2*y0-y) and compare to the
% actual polygon's vertex list WITHOUT changing vertex coordinates.
y0=mouth(2);
distVertex=zeros(n,1);
for i=1:n
    q=[P(i,1),2*y0-P(i,2)];
    distVertex(i)=min(hypot(P(:,1)-q(1),P(:,2)-q(2)));
end
% A mirrored boundary vertex alone is insufficient: a different polygon
% can have the same vertex set but different edges. Check the FULL
% undirected segment set with the same 1:1 prescribed boundary.
distEdge=zeros(n,1);
starts=P;
ends=P(jNext,:);
for i=1:n
    qa=[starts(i,1),2*y0-starts(i,2)];
    qb=[ends(i,1),2*y0-ends(i,2)];
    % Both orientations, since reflection reverses polygon traversal.
    ds=hypot(starts(:,1)-qa(1),starts(:,2)-qa(2));
    de=hypot(ends(:,1)-qb(1),ends(:,2)-qb(2));
    forward=max(ds,de);
    dsr=hypot(ends(:,1)-qa(1),ends(:,2)-qa(2));
    der=hypot(starts(:,1)-qb(1),starts(:,2)-qb(2));
    reverse=max(dsr,der);
    distEdge(i)=min([forward;reverse]);
end
G=D.channelGeom.append;
required={'Mup','Mlo','xtip','face_upper','face_lower'};
for k=1:numel(required)
    if ~isfield(G,required{k})
        error('step54:IncompletePencilMetadata', ...
            'Missing sharp-pencil geometry metadata %s.',required{k});
    end
end
mouthPair=hypot(G.Mup(1)-G.Mlo(1), ...
    G.Mup(2)+G.Mlo(2)-2*y0);
tipY=abs(G.xtip(2)-y0);
faceUpper=G.face_upper;
faceLower=flipud(G.face_lower);
if ~isequal(size(faceUpper),size(faceLower))
    error('step54:FaceGeometry','Pencil faces have different vertex counts.');
end
faceLower(:,2)=2*y0-faceLower(:,2);
faceMirror=max(hypot(faceUpper(:,1)-faceLower(:,1), ...
                     faceUpper(:,2)-faceLower(:,2)));
% Check the upper/lower quarter-hole arc vertex counts separately (the
% overall vertex/edge tests above are the primary complete test).
isUpperArc= P(:,1)>=C.hole.center(1)-tol & ...
    P(:,1)<=C.hole.center(1)+C.hole.r+tol & ...
    P(:,2)>y0+tol & ...
    abs(hypot(P(:,1)-C.hole.center(1),P(:,2)-y0)-C.hole.r)<=tol;
isLowerArc= P(:,1)>=C.hole.center(1)-tol & ...
    P(:,1)<=C.hole.center(1)+C.hole.r+tol & ...
    P(:,2)<y0-tol & ...
    abs(hypot(P(:,1)-C.hole.center(1),P(:,2)-y0)-C.hole.r)<=tol;
nUpperArc=nnz(isUpperArc);
nLowerArc=nnz(isLowerArc);
vmax=max(distVertex);
emax=max(distEdge);
allPass=vmax<=tol && emax<=tol && ...
    mouthPair<=tol && tipY<=tol && ...
    faceMirror<=tol && nUpperArc==nLowerArc;
summary=table(n,nUpperArc,nLowerArc,vmax,emax, ...
    mouthPair,faceMirror,tipY, ...
    nnz(distVertex>tol),nnz(distEdge>tol),allPass, ...
    'VariableNames',{'nOuterPolygonVertices', ...
    'nUpperArcVertices','nLowerArcVertices', ...
    'maxReflectedVertexError_m','maxReflectedEdgeError_m', ...
    'mouthReflectionError_m','pencilFaceReflectionError_m', ...
    'tipOffsetFromMirrorAxis_m', ...
    'unpairedVertices','unpairedBoundaryEdges', ...
    'prescribedPolygonReflectionPass'});
fprintf('\n============================================================\n');
fprintf('STEP 54: PRESCRIBED BOUNDARY GEOMETRY MIRROR AUDIT\n');
fprintf('============================================================\n');
fprintf('  Same saved Step45/47 Npoly=240, a0=%.7g m endpoints.\n',a.a0);
fprintf('  Reflection across y=%.8g m; tolerance %.3e m.\n',y0,tol);
disp(summary);
if ~allPass
    warning('step54:PrescribedGeometryNotReflectionSymmetric', ...
        ['Do not introduce a reflection-paired mesh by silently ', ...
         'changing this failing prescribed polygon.']);
else
    fprintf(['  Prescribed polygon is reflection-paired at both ', ...
        'vertex and complete-edge levels to the stated tolerance.\n']);
end
O54=struct('summary',summary,'checks',struct( ...
    'polygonVertexDist_m',distVertex, ...
    'polygonEdgeDist_m',distEdge, ...
    'tolerance_m',tol, ...
    'fullPass',allPass), ...
    'crackEndpoints',Pmid,'checkpointPaths',{{oldFile,newFile}}, ...
    'note',['This reconstructs the prescribed boundary only; ', ...
    'it does NOT verify the unsaved internal PDE preprocessor ', ...
    'or a new reflection-paired volume mesh.'], ...
    'noNewMesh',true,'noNewFEM',true,'noEDI',true);
path=char(opt.SaveFile);
[folder,~,~]=fileparts(path);
if ~isempty(folder)&&exist(folder,'dir')~=7,mkdir(folder);end
save(path,'O54'); % compact geometry diagnostics, not FEM meshes
fprintf('  Compact Step54 geometry report saved: %s\n',path);
end