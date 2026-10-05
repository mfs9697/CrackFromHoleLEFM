function out=test_step62_saved_candidate(path)
%TEST_STEP62_SAVED_CANDIDATE Reload integrity, reflection and support closure.
% Reads only the accepted T3 artifact; reconstructs T6 without physical U.
folder=fileparts(mfilename('fullpath'));root=fileparts(fileparts(folder));
if nargin<1,path=fullfile(root,'verification','step62_structured_graded_mesh_candidate_T3.mat');end
inventory=whos('-file',path);assert(numel(inventory)==1&&strcmp(inventory.name,'candidate'));
s=load(path,'candidate');c=s.candidate;
assert(~any(isfield(c,{'U','u','K','stiffness','forces'})));
assert(c.gates.structuralPass&&all(isfinite(c.p(:)))&&all(c.t(:)>0));
names={'main_step62_structured_graded_mesh','build_step62_structured_patch', ...
    'build_step62_graded_exterior'};
hashes={c.provenance.driverSHA256,c.provenance.patchBuilderSHA256, ...
    c.provenance.exteriorBuilderSHA256};
for k=1:3
    assert(strcmp(file_hash(fullfile(folder,[names{k} '.m'])),hashes{k}), ...
        'step62:SourceIntegrity','Source bytes differ from candidate provenance.');
end
tip=c.crack.Pmid(end,:);e=diff(c.crack.Pmid)/norm(diff(c.crack.Pmid));
R=[e(:),[-e(2);e(1)]];Z=(c.p-tip)*R;
ids=c.pairedNodeMap;nU=numel(c.localMirrorMap);
A=Z(ids(1:nU),:);B=Z(ids(c.localMirrorMap),:);
reflection=max(vecnorm(B-A.*[1 -1],2,2));assert(reflection<1e-12);
map=zeros(size(c.p,1),1);map(ids(1:nU))=ids(c.localMirrorMap);
T=c.t(c.pairedElementIDs,:);half=size(T,1)/2;
assert(isequal(map(T(1:half,[1 3 2])),T(half+1:end,:)));
assert(nnz(any(c.t==c.crack.tipNode,2))==6);
assert(numel(setdiff(unique(c.t(any(c.t==c.crack.tipNode,2),:)),c.crack.tipNode))==7);
[P6,T6]=T3toT6_fast(c.p,c.t);Z6=(P6-tip)*R;
assert(all(ismember(c.literalPrimarySupportElementIDs,c.pairedElementIDs)));
nodes=unique(T6(c.literalPrimarySupportElementIDs,:));
assert(all(inpolygon(Z6(nodes,1),Z6(nodes,2), ...
    c.pairedBoundaryLocal(:,1),c.pairedBoundaryLocal(:,2))));
assert(isempty(c.retainedOriginalElementIDs));
if c.scientificallyReadyForFEMProposal,assert(c.synthetic.passed);end
out=table(c.structuredDesign.level,size(c.p,1),size(c.t,1),size(P6,1), ...
    reflection,max(vecnorm(Z6(nodes,:),2,2))*1e3,c.scientificallyReadyForFEMProposal, ...
    'VariableNames',{'level','T3Nodes','T3Triangles','T6Nodes','reflectionError_m', ...
    'literalSupportMaxNodeRadius_mm','prescribedFieldsQualified'});
disp(out);fprintf('PASS: saved candidate source hashes, T6 rebuild, fan, reflection, full support closure.\n');
end
function h=file_hash(path)
f=fopen(path,'rb');assert(f>=0);guard=onCleanup(@()fclose(f));
md=java.security.MessageDigest.getInstance('SHA-256');
while ~feof(f),b=fread(f,1024*1024,'*uint8');md.update(typecast(b,'int8'));end
h=lower(reshape(dec2hex(typecast(md.digest(),'uint8'),2).',1,[]));clear guard
end
