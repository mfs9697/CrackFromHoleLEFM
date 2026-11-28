% IAM-UA8-F02: read-only verification of the regenerated meshes and raw export of the P21 mesh for the figure.
% 1 regenerated P5 vs the preserved historical P5 candidate (profiling run of the same accepted path): bitwise T6 nodes/connectivity
%   and crack-face node lists; 2 regenerated P5/P21 fingerprints vs paper/data/qualification.csv (production records);
% 3 raw little-endian float64 exports of the P21 arrays (global mm), the crack faces, the paired core, the EDI support and the mirror map.
out='<scratch>\f02';
wt='<repo>';
pres='<investigator_codex_work>\2026-10-04\files-pasted-by-the-user-you\work\incremental-profile\phaseA_baseline\fresh5\step_005_candidate.mat';
V=struct();
P=load(pres);c=P.candidate;A=load(fullfile(out,'F02_regen_P05.mat'));F5=A.F;
V.preservedFile=strrep(pres,'<home>','<home>');V.preservedSHA256=lower(sha(pres));
V.preservedSizes=struct('p',size(c.p),'t',size(c.t));V.regenSizes=struct('coord',size(F5.mesh.coord),'connect',size(F5.mesh.connect),'coord3',size(F5.mesh.coord3),'connect3',size(F5.mesh.connect3));
V.preservedCrackFields=fieldnames(c.crack);
V.P5_T6nodesBitwise=isequal(c.p,F5.mesh.coord);V.P5_T6connectExact=isequal(double(c.t),double(F5.mesh.connect));
up=fieldOr(c.crack,'upperNodes');lo=fieldOr(c.crack,'lowerNodes');
V.P5_crackFacesExact=isequal(up(:),F5.candidate.crack.upperNodes(:))&&isequal(lo(:),F5.candidate.crack.lowerNodes(:));
V.P5_pairedElementsExact=isequal(c.pairedElementIDs(:),F5.candidate.pairedElementIDs(:));
V.P5_mirrorMapExact=isequal(c.pairedMirrorLocal,F5.candidate.pairedMirrorLocal);
V.P5_maxNodeDiff_m=NaN;if isequal(size(c.p),size(F5.mesh.coord)),V.P5_maxNodeDiff_m=max(abs(c.p-F5.mesh.coord),[],'all');end
% fingerprints vs production qualification records
Ev=load(fullfile(wt,'paper','data','evidence_exact.mat'),'E');Q=Ev.E.qualification;if isstruct(Q),Q=struct2table(Q);end   % full precision (the CSV export has 15 significant digits)
cols={'T3_nodes','T3_elements','T6_nodes','core_T3_elements','exterior_T3_elements','EDI_elements','physical_clearance_m','exterior_crack_m', ...
 'area_rel_error','min_angle_deg','max_neighbor_ratio','min_T6_detJ_m2','core_coord_error_m','pair_error_T3_m','pair_error_T6_m'};
for k=[5 21]
 B=load(fullfile(out,sprintf('F02_regen_P%02d.mat',k)));S=B.F.summary;row=Q(Q.segment==k,:);
 r=struct();for j=1:numel(cols),a=S.(cols{j});b=row.(cols{j});r.(cols{j})=struct('regenerated',a,'production',b,'equal',isequal(a,b),'relDiff',abs(a-b)/max(abs(b),realmin));end
 V.(sprintf('P%d_fingerprint',k))=r;
 fprintf('P%d fingerprint equal: %s\n',k,mat2str(cellfun(@(f)r.(f).equal,cols)));
end
% raw export of P21 (global coordinates in mm)
B=load(fullfile(out,'F02_regen_P21.mat'));F=B.F;cd2=F.candidate;raw=fullfile(out,'raw');if ~isfolder(raw),mkdir(raw);end;IX=struct('file',{},'shape',{});
IX=ex(IX,raw,'X3_mm',1e3*F.mesh.coord3);IX=ex(IX,raw,'T3',double(F.mesh.connect3));IX=ex(IX,raw,'X6_mm',1e3*F.mesh.coord);IX=ex(IX,raw,'T6',double(F.mesh.connect));
IX=ex(IX,raw,'XL3_mm',1e3*F.localMesh.coord3);IX=ex(IX,raw,'XL6_mm',1e3*F.localMesh.coord);
IX=ex(IX,raw,'upperNodes',double(F.crackLocal.upperNodes(:)));IX=ex(IX,raw,'lowerNodes',double(F.crackLocal.lowerNodes(:)));
IX=ex(IX,raw,'upperS_m',F.crackLocal.upperS(:));IX=ex(IX,raw,'lowerS_m',F.crackLocal.lowerS(:));IX=ex(IX,raw,'tipNode',double(F.crackLocal.tipNode));
IX=ex(IX,raw,'pathGlobal_mm',1e3*cd2.path);IX=ex(IX,raw,'pairedElementIDs',double(cd2.pairedElementIDs(:)));IX=ex(IX,raw,'exteriorElementIDs',double(cd2.exteriorElementIDs(:)));
IX=ex(IX,raw,'EDIsupportIDs',double(cd2.primarySupportElementIDs(:)));IX=ex(IX,raw,'pairedNodeIDs',double(cd2.pairedNodeIDs(:)));
mm=cd2.pairedMirrorLocal;V.mirrorMapClass=class(mm);V.mirrorMapSize=size(mm);if isnumeric(mm),IX=ex(IX,raw,'pairedMirrorLocal',double(mm));end
if isfield(F.crackLocal,'lowerMatchForUpper'),IX=ex(IX,raw,'lowerMatchForUpper',double(F.crackLocal.lowerMatchForUpper(:)));end
V.P21=struct('tip_mm',1e3*cd2.path(end,:),'directionLast',cd2.directionLast,'thetaCurrentDeg',cd2.thetaCurrentDeg,'design',small(cd2.structuredDesign),'candidateFields',{fieldnames(cd2)});
V.R0hole=struct('A_m',0.3,'B_m',0.1);
fid=fopen(fullfile(raw,'RAW_INDEX.json'),'w');fwrite(fid,jsonencode(IX),'char');fclose(fid);
fid=fopen(fullfile(out,'F02_VERIFY.json'),'w');fwrite(fid,jsonencode(V,'PrettyPrint',true),'char');fclose(fid);
fprintf('P5 bitwise nodes=%d connect=%d faces=%d paired=%d mirror=%d maxdiff=%g\n',V.P5_T6nodesBitwise,V.P5_T6connectExact,V.P5_crackFacesExact,V.P5_pairedElementsExact,V.P5_mirrorMapExact,V.P5_maxNodeDiff_m);
fprintf('VERIFY_EXPORT_DONE raw=%d mirror class=%s size=%s\n',numel(IX),V.mirrorMapClass,mat2str(V.mirrorMapSize));
function IX=ex(IX,raw,name,x),f=fullfile(raw,[name '.f64']);fid=fopen(f,'w','l');fwrite(fid,double(x(:)),'float64');fclose(fid);IX(end+1)=struct('file',[name '.f64'],'shape',size(x));end
function v=fieldOr(S,f),if isstruct(S)&&isfield(S,f),v=S.(f);else,v=[];end,end
function h=sha(f),md=java.security.MessageDigest.getInstance('SHA-256');fid=fopen(f,'r');d=fread(fid,inf,'*uint8');fclose(fid);md.update(d);h=sprintf('%02X',typecast(md.digest(),'uint8'));end
function S=small(S)
if ~isstruct(S)||~isscalar(S),return,end
f=fieldnames(S);for k=1:numel(f),x=S.(f{k});if isstruct(x)&&isscalar(x),S.(f{k})=small(x);elseif (isnumeric(x)||islogical(x))&&numel(x)>64,S=rmfield(S,f{k});elseif iscell(x)||isa(x,'function_handle'),S=rmfield(S,f{k});end,end
end
