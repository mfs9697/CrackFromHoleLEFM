function O62=main_step62_structured_graded_mesh(varargin)
%MAIN_STEP62_STRUCTURED_GRADED_MESH One mesh-only asymmetric Step38 candidate.
% NEVER loads physical U, assembles stiffness, or solves equilibrium.
% New explicit sector/ring topology and a completely remeshed exterior.
% Audits reuse Step61 conventions; no Step61 candidate geometry is reused.
% Prescribed-field EDI is permitted ONLY after all structural gates pass.
root=fileparts(fileparts(fileparts(mfilename('fullpath'))));
ip=inputParser;
addParameter(ip,'CheckpointFile',fullfile(root,'verification', ...
    'step38_tip_refined_solved.mat'),@(x)ischar(x)||isstring(x));
addParameter(ip,'SourceCandidateFile','',@(x)ischar(x)||isstring(x));
addParameter(ip,'SavePrefix',fullfile(root,'verification', ...
    'step62_structured_graded_mesh'),@(x)ischar(x)||isstring(x));
addParameter(ip,'Visible','off',@(x)ischar(x)||isstring(x));
addParameter(ip,'RunSynthetic',true,@(x)islogical(x)&&isscalar(x));
addParameter(ip,'ExteriorCalibration',struct(),@(x)isstruct(x)&&isscalar(x));
addParameter(ip,'WriteArtifacts',true,@(x)islogical(x)&&isscalar(x));
addParameter(ip,'ReturnCandidate',false,@(x)islogical(x)&&isscalar(x));
addParameter(ip,'Verbose',true,@(x)islogical(x)&&isscalar(x));
addParameter(ip,'Level',0,@(x)isnumeric(x)&&isscalar(x)&& ...
    x==fix(x)&&x>=0&&x<=4);
addParameter(ip,'DebugFile','',@(x)ischar(x)||isstring(x));
parse(ip,varargin{:}); opt=ip.Results;
addpath(genpath(root));
assert_audit_branch(root);
cp=char(opt.CheckpointFile); sourceCandidateFile=char(opt.SourceCandidateFile);
prefix=char(opt.SavePrefix);
if opt.Level>0,prefix=sprintf('%s_L%d',prefix,opt.Level);end

% Two source modes are supported. Historical mode selectively loads the
% Step38 checkpoint without U. Recovery mode uses the exact archived Step62
% candidate when the untracked Step38 checkpoint is unavailable.
if ~isempty(sourceCandidateFile)
    d=load(sourceCandidateFile,'candidate');
    assert(isfield(d,'candidate'),'step62:SourceCandidate', ...
        'SourceCandidateFile does not contain candidate.');
    src=d.candidate;
    for ff={'p','t','crack','mat','structuredDesign','provenance'}
        assert(isfield(src,ff{1}),'step62:SourceCandidateField', ...
            'Source candidate missing field %s.',ff{1});
    end
    P=src.p;T=src.t;cr=src.crack;mat=src.mat;
    a0=norm(diff(cr.Pmid));hTip=src.structuredDesign.hBase_m;
    ri=src.provenance.rInner;ro=src.provenance.rOuter;
    [Psrc6,Tsrc6]=T3toT6_fast(P,T);
    sourceMesh=struct('coord3',P,'connect3',T,'coord',Psrc6,'connect',Tsrc6);
    if isfield(src.provenance,'baselineOuterRatios')
        sourceOuterRatios=src.provenance.baselineOuterRatios;
    else
        sourceOuterRatios=[.50 .65 .80];
    end
    sourceMode='archived_step62_candidate';
    sourcePath=sourceCandidateFile;sourceSHA=sha256(sourceCandidateFile);
    checkpointPath='';
    expectedSupport=src.primarySupportElementIDs(:);
    if isfield(src,'literalPrimarySupportElementIDs')
        expectedLiteralSupport=src.literalPrimarySupportElementIDs(:);
    else
        expectedLiteralSupport=[];
    end
else
    % Intentionally selective: U, forces and stiffness are NEVER loaded.
    d=load(cp,'mesh','crack','baseline','actualTip','a0','mat');
    P=d.mesh.coord3;T=d.mesh.connect3;cr=d.crack;mat=d.mat;
    a0=d.a0;hTip=d.actualTip;ri=d.baseline.rInner;ro=.65*a0;
    sourceMesh=d.mesh;sourceOuterRatios=d.baseline.rOuterOverA0;
    sourceMode='step38_checkpoint';sourcePath=cp;sourceSHA=sha256(cp);
    checkpointPath=cp;
    expectedSupport=[];
    expectedLiteralSupport=[];
    assert(size(P,1)==20164&&size(T,1)==39441&& ...
        size(d.mesh.coord,1)==79769,'step62:SourceCount','Wrong Step38 mesh.');
end
assert(abs(a0-.008)<1e-12&&abs(ri-.0008)<1e-12&& ...
    abs(hTip-5.4024650785e-5)<1e-12,'step62:SourceScale','Wrong source scales.');
assert(abs(ro-.65*a0)<1e-12,'step62:SourceOuterRadius','Wrong primary EDI radius.');
assert(norm(cr.Pmid(1,:)-[.199988872196,-.0208170339049])<2e-12&& ...
    norm(cr.Pmid(end,:)-[.207985904782,-.0210349096129])<2e-12, ...
    'step62:SourceCrack','Wrong saved crack endpoints.');
sourceMat=mat;
tip=cr.Pmid(end,:);e1=diff(cr.Pmid)/a0;
R=[e1(:),[-e1(2);e1(1)]];X=(P-tip)*R;
if opt.Verbose
    fprintf('STEP62: source=%s, no physical U loaded.\n',sourceMode);
end
% Step60 convention and literal production participation are both audited.
oldSupport=q_support(sourceMesh,cr,ri,ro,true);
oldLiteralSupport=q_support(sourceMesh,cr,ri,ro,false);
if strcmp(sourceMode,'step38_checkpoint')
    assert(numel(oldSupport)==14215,'step62:SourceSupport','Step60 support changed.');
else
    assert(isequal(oldSupport(:),expectedSupport), ...
        'step62:SourceSupport','Archived Step62 support changed on reload.');
    if ~isempty(expectedLiteralSupport)
        assert(isequal(oldLiteralSupport(:),expectedLiteralSupport), ...
            'step62:SourceLiteralSupport', ...
            'Archived Step62 literal support changed on reload.');
    end
end
supportNodes=unique(sourceMesh.connect(oldSupport,:));
maxSupport=max(vecnorm((sourceMesh.coord(supportNodes,:)-tip)*R,2,2));
oldQuality=quality(X,T);
[oldFree,~]=edge_inventory(T);
oldFree=oldFree(:,1:2);
upOld=unique(cr.upperNodes(:)); loOld=unique(cr.lowerNodes(:));
isFace=face_edge_mask(oldFree,upOld,loOld,cr.tipNode);
physicalOther=oldFree(~isFace,:);
physicalClearance=min_segment_distance([0 0],X,physicalOther);

% Predeclared 6 mm structured extent continues grading beyond the fixed
% 5.2 mm annulus. User-authorized full exterior remeshing retains every
% exact physical boundary geometry, but replaces all original triangles.
rp=.006;removed=true(size(T,1),1);outside=~removed;
assert(rp-maxSupport>.15e-3&&physicalClearance-rp>.4e-3);
extentTable=table(rp,maxSupport,rp-maxSupport,physicalClearance, ...
    'VariableNames',{'pairedRadius_m','oldSupportMaxRadius_m', ...
    'supportClearance_m','nearestPhysicalBoundary_m'});
[Zpatch,Tpatch,mirrorMap,nUpper,ringRadii,axisUp,design]= ...
    build_step62_structured_patch(rp,hTip,opt.Level);
design.exteriorCalibration=opt.ExteriorCalibration;
[Zexterior,Texterior,exteriorMeta,physicalIDs]= ...
    build_step62_graded_exterior(X,T,cr,Zpatch,Tpatch,rp,design);
outerIDs=exteriorMeta.originalPhysicalIDs;outerXY=X(physicalIDs,:);
% Repeat the entire construction to verify exact deterministic output.
[Zagain,Tagain]=build_step62_structured_patch(rp,hTip,opt.Level);
[Eagain,Fagain]=build_step62_graded_exterior(X,T,cr,Zagain,Tagain,rp,design);
deterministic=isequal(Zpatch,Zagain)&&isequal(Tpatch,Tagain)&& ...
    isequal(Zexterior,Eagain)&&isequal(Texterior,Fagain);

% Map the new geometry onto existing interface IDs. Retained old vertices
% stay bitwise unchanged; no coordinate welding across the crack.
Q=P; pairedIDs=zeros(size(Zpatch,1),1);
for j=1:size(Zpatch,1)
    Q(end+1,:)=tip+Zpatch(j,:)*R.'; pairedIDs(j)=size(Q,1); %#ok<AGROW>
end
patchTriangles=pairedIDs(Tpatch);
exteriorIDs=zeros(size(Zexterior,1),1);
for j=1:size(Zexterior,1)
    z=Zexterior(j,:);
    side=exteriorMeta.nodeSide(j);
    candidates=outerIDs;
    if z(1)<0&&abs(z(2))<1e-12
        if side>=0,candidates=setdiff(candidates,loOld); ...
        else,candidates=setdiff(candidates,upOld);end
    end
    d=vecnorm(X(candidates,:)-z,2,2); [dd,ii]=min(d);
    if dd<1e-12
        exteriorIDs(j)=candidates(ii); continue
    end
    d=vecnorm(Zpatch-z,2,2);
    if z(1)<0&&abs(z(2))<1e-12
        upperAxis=axisUp(Zpatch(axisUp,1)<0);
        lowerAxis=mirrorMap(upperAxis);
        if side>=0,d(lowerAxis)=Inf;else,d(upperAxis)=Inf;end
    end
    [dd,ii]=min(d);
    if dd<1e-12,exteriorIDs(j)=pairedIDs(ii);continue,end
    Q(end+1,:)=tip+z*R.'; exteriorIDs(j)=size(Q,1); %#ok<AGROW>
end
Tc=[T(outside,:);patchTriangles;exteriorIDs(Texterior)];
% Compact unused old interior nodes; retain exact exterior coordinates/IDs
% through a saved original-to-candidate map.
used=unique(Tc(:)); oldToNew=zeros(size(Q,1),1);
oldToNew(used)=(1:numel(used))'; Pc=Q(used,:); Tc=oldToNew(Tc);
pairedIDs=oldToNew(pairedIDs); exteriorIDs=oldToNew(exteriorIDs);
nExterior=nnz(outside); nPatch=size(Tpatch,1);
patchRows=nExterior+(1:nPatch)';
exteriorRows=(nExterior+nPatch+1:size(Tc,1))';
Xc=(Pc-tip)*R;
% Derive crack IDs from topology, never from coincident coordinates alone.
crNew=cr; crNew.tipNode=pairedIDs(1);
newUpper=pairedIDs(axisUp(Zpatch(axisUp,1)<-1e-12));
newLower=pairedIDs(mirrorMap(axisUp(Zpatch(axisUp,1)<-1e-12)));
onExterior=abs(Zexterior(:,2))<1e-12 & Zexterior(:,1)<0 & ...
    Zexterior(:,1)>=-a0-1e-12;
newUpper=[newUpper;exteriorIDs(onExterior&exteriorMeta.nodeSide>=0)];
newLower=[newLower;exteriorIDs(onExterior&exteriorMeta.nodeSide<0)];
keptUpper=oldToNew(upOld); keptLower=oldToNew(loOld);
crNew.upperNodes=unique([keptUpper(keptUpper>0);newUpper;crNew.tipNode]);
crNew.lowerNodes=unique([keptLower(keptLower>0);newLower;crNew.tipNode]);
crNew=refresh_crack_metadata(crNew,Pc,e1,a0);
[P6,T6]=T3toT6_fast(Pc,Tc);
candidateMesh=struct('coord3',Pc,'connect3',Tc,'coord',P6,'connect',T6);
if strlength(string(opt.DebugFile))>0
    save(char(opt.DebugFile),'candidateMesh','crNew','P','T','Xc','design','exteriorMeta','-v7');
end
newSupport=q_support(candidateMesh,crNew,ri,ro,true);
newLiteralSupport=q_support(candidateMesh,crNew,ri,ro,false);
newQuality=quality(Xc,Tc);
gates=struct();
gates.physicalGeometryUnchanged=boundary_same(P,T,Pc,Tc,cr,crNew) && ...
    isequal(cr.Pmid,crNew.Pmid);
gates.physicalBoundaryNodesUnchanged= ...
    isequal(Pc(oldToNew(physicalIDs),:),P(physicalIDs,:));
gates.materialUnchanged=isequal(mat,sourceMat);
gates.deterministicWholeMesh=deterministic;
[pair3,pair6]=pair_errors(candidateMesh,pairedIDs,mirrorMap,nUpper, ...
    patchRows,tip,R);
gates.completeT3Pairing=pair3<1e-12;
gates.completeT6Pairing=pair6<1e-12;
[topology,~,nativeR]=topology_audit(candidateMesh,crNew, ...
    pairedIDs,axisUp,mirrorMap,Zpatch);
gates.crackFacesDistinct=topology.crackDistinct;
gates.intactLigamentShared=topology.intactShared;
gates.crackT6MidsidesDistinct=topology.crackMidsDistinct;
gates.intactT6MidsidesShared=topology.intactMidsShared;
gates.positiveT3Areas=all(newQuality.area>0);
[minJ,midsError]=jacobian_audit(P6,T6);
gates.positiveT6Jacobians=minJ>0&&midsError<1e-12;
gates.noHangingGapsOverlaps=global_planar_audit(P,T,Pc,Tc, ...
    Xc,patchRows,exteriorRows,exteriorMeta.innerPolygon,crNew);
gates.originalSupportInsidePairedRegion= ...
    all(inpolygon((sourceMesh.coord(supportNodes,:)-tip)*e1.', ...
    (sourceMesh.coord(supportNodes,:)-tip)*R(:,2), ...
    exteriorMeta.innerPolygon(:,1),exteriorMeta.innerPolygon(:,2)));
gates.candidateSupportInsidePairedRegion=all(ismember(newSupport,patchRows))&& ...
    all(ismember(newLiteralSupport,patchRows));
gates.exteriorOutsideSupport=~any(ismember(newSupport,exteriorRows))&& ...
    ~any(ismember(newLiteralSupport,exteriorRows));
windows=[.04 .20;.04 .30;.08 .30;.12 .30];
sampleN=zeros(4,1);
for k=1:4,sampleN(k)=nnz(nativeR/a0>=windows(k,1)& ...
    nativeR/a0<=windows(k,2));end
gates.nativeSamplingAdequate=all(sampleN>=12);
newTip=tip_edges(Pc,Tc,crNew.tipNode);
oldTip=tip_edges(P,T,cr.tipNode);
sourceLabel=char(string(sourceMode));
tipEdgeTable=table({sourceLabel;'Candidate'},[numel(oldTip);numel(newTip)], ...
    [min(oldTip);min(newTip)]*1e3,[median(oldTip);median(newTip)]*1e3, ...
    [prctile(oldTip,90);prctile(newTip,90)]*1e3, ...
    [max(oldTip);max(newTip)]*1e3,'VariableNames', ...
    {'mesh','incidentEdges','min_mm','median_mm','p90_mm','max_mm'});
gates.tipResolutionMatched=abs(median(newTip)/(hTip*design.scale)-1)<=.05;
gates.tipFanSixTriangles=nnz(any(Tc==crNew.tipNode,2))==6&&numel(newTip)==7;
gates.monotoneStructuredLaw=all(diff(design.ringTable.targetH_mm)>=0)&& ...
    all(diff(design.ringTable.bandWidth_mm(2:end))>=0);
gates.monotoneExteriorLaw=all(diff(exteriorMeta.ringTable.targetH_m)>=0)&& ...
    all(diff(exteriorMeta.ringTable.bandWidth_m)>=0);
% Engineering gates, not physical error bars.
affected=[patchRows;exteriorRows];
gates.minimumAngle20deg=min(newQuality.minAngle(affected))>=20;
neighborTarget=2.5;
if isfield(opt.ExteriorCalibration,'neighborRatioTarget')
    neighborTarget=opt.ExteriorCalibration.neighborRatioTarget;
end
gates.noAbruptGrading=max_neighbor_size_ratio(Tc,newQuality.longest,affected)<=neighborTarget;
gates.structuralPass=all(structfun(@(x)logical(x),gates));
summary=table(size(P,1),size(T,1),size(sourceMesh.coord,1), ...
    size(Pc,1),size(Tc,1),size(P6,1),nnz(removed),nnz(outside), ...
    rp*1e3,physicalClearance*1e3, ...
    max(vecnorm(outerXY,2,2))*1e3,median(oldTip)*1e3,median(newTip)*1e3, ...
    numel(oldSupport),numel(newSupport), ...
    min(newQuality.minAngle(patchRows)),min(newQuality.minAngle(exteriorRows)), ...
    pair3,pair6,'VariableNames',{'oldT3Nodes','oldT3','oldT6Nodes', ...
    'newT3Nodes','newT3','newT6Nodes','replacedOldT3','retainedOldT3', ...
    'pairedRadius_mm','nearestPhysicalBoundary_mm','maxPhysicalRadius_mm', ...
    'oldTipMedian_mm','newTipMedian_mm','oldPrimarySupport','newPrimarySupport', ...
    'patchMinAngle_deg','exteriorMinAngle_deg','pairErrorT3_m','pairErrorT6_m'});
qualityTable=quality_comparison(oldQuality,removed,newQuality,patchRows,exteriorRows);
radialTable=radial_comparison(X,T,Xc,Tc,oldQuality,newQuality,rp);
samplingTable=table(windows(:,1),windows(:,2),sampleN, ...
    'VariableNames',{'lower_r_over_a0','upper_r_over_a0','nativePoints'});
if opt.Verbose
    disp(summary);disp(qualityTable);disp(samplingTable);disp(gates);
end
if opt.WriteArtifacts
    files=plot_candidate(prefix,opt.Visible,X,T,Xc,Tc,removed,patchRows, ...
        exteriorRows,oldSupport,newSupport,outerXY,exteriorMeta.innerPolygon,rp,ri,ro, ...
        cr,crNew,oldQuality,newQuality,design,exteriorMeta);
else
    files=struct();
end
provenance=struct('sourceMode',sourceMode,'sourcePath',sourcePath, ...
    'sourceSHA256',sourceSHA,'checkpointPath',checkpointPath, ...
    'branch',current_branch(root),'sourceCommit',git_head(root), ...
    'driverSHA256',sha256(mfilename('fullpath')), ...
    'physicalUWasLoaded',false,'sourceCrack',cr.Pmid, ...
    'sourceMaterial',mat,'baselineOuterRatios',sourceOuterRatios, ...
    'rInner',ri,'rOuter',ro,'upperGeneratedOnce',true, ...
    'outsideTrianglesExact',false,'allOriginalTrianglesReplaced',true, ...
    'physicalBCsAndLoadsUnmodified',true, ...
    'patchBuilderSHA256',sha256(fullfile(fileparts(mfilename('fullpath')), ...
        'build_step62_structured_patch.m')), ...
    'exteriorBuilderSHA256',sha256(fullfile(fileparts(mfilename('fullpath')), ...
        'build_step62_graded_exterior.m')), ...
    'matlabVersion',version,'supportDefinition','Step60 plus literal EDI audit', ...
    'exteriorCalibration',opt.ExteriorCalibration);
O62=struct('feasible',true,'summary',summary,'extentTable',extentTable, ...
    'qualityTable',qualityTable,'radialTable',radialTable, ...
    'samplingTable',samplingTable,'gates',gates,'files',files, ...
    'tipEdgeTable',tipEdgeTable,'originalTipEdgeDistribution_m',oldTip, ...
    'targetTipEdge_m',hTip*design.scale,'design',design,'exteriorDesign',exteriorMeta, ...
    'ringTable',design.ringTable,'exteriorRingTable',exteriorMeta.ringTable, ...
    'exteriorRefinementTable',exteriorMeta.refinementTable, ...
    'topology',topology,'minT6Jacobian',minJ, ...
    'maxNeighborSizeRatio',max_neighbor_size_ratio(Tc,newQuality.longest,affected), ...
    'originalLiteralSupportCount',numel(oldLiteralSupport), ...
    'candidateLiteralSupportCount',numel(newLiteralSupport), ...
    'minimumNoncrackBoundaryDistance_m',physicalClearance, ...
    'excludedZeroInteriorHullSlivers', ...
        exteriorMeta.excludedZeroInteriorHullSlivers, ...
    'tipEdgeDistribution_m',newTip,'provenance',provenance, ...
    'synthetic',struct('performed',false,'passed',false), ...
    'readyForOneAsymmetricFEMProposal',false, ...
    'noPhysicalU',true,'noFEM',true,'noStiffness',true);
% Save diagnostic reports even on rejection; only structurally accepted
% candidates receive the candidate MAT file.
if opt.WriteArtifacts,save_report(prefix,O62);end
if ~gates.structuralPass
    if opt.Verbose,fprintf('STOP: structural/design gate failed. No synthetic EDI, no FEM.\n');end
    return
end
candidate=struct('p',Pc,'t',Tc,'crack',crNew,'mat',mat, ...
    'pairedElementIDs',patchRows,'exteriorElementIDs',exteriorRows, ...
    'primarySupportElementIDs',newSupport, ...
    'literalPrimarySupportElementIDs',newLiteralSupport, ...
    'pairedNodeMap',pairedIDs,'localMirrorMap',mirrorMap, ...
    'pairedBoundaryLocal',exteriorMeta.innerPolygon,'physicalBoundaryLocal',outerXY, ...
    'replacedOriginalElementIDs',find(removed), ...
    'retainedOriginalElementIDs',find(outside), ...
    'originalNodeToCandidate',oldToNew(1:size(P,1)), ...
    'ringRadii_m',ringRadii,'structuredDesign',design, ...
    'exteriorDesign',exteriorMeta,'provenance',provenance,'gates',gates);
candidateFile=[prefix '_candidate_T3.mat'];
candidate.scientificallyReadyForFEMProposal=false;
if opt.ReturnCandidate,O62.candidate=candidate;end
if opt.WriteArtifacts
    save(candidateFile,'candidate','-v7'); O62.files.candidateMAT=candidateFile;
end
if opt.RunSynthetic
    O62.synthetic=synthetic_controls(candidateMesh,crNew,mat,ri,ro);
    candidate.scientificallyReadyForFEMProposal=O62.synthetic.passed;
    candidate.synthetic=O62.synthetic;
    if opt.ReturnCandidate,O62.candidate=candidate;end
    if opt.WriteArtifacts,save(candidateFile,'candidate','-v7');end
end
O62.readyForOneAsymmetricFEMProposal=gates.structuralPass&&O62.synthetic.passed;
if opt.WriteArtifacts,save_report(prefix,O62);end
if opt.Verbose
    fprintf('Step62 finished: proposal ready=%d; zero FEM solves.\n', ...
        O62.readyForOneAsymmetricFEMProposal);
end
end

function p=conditional_checkpoint(mode,cp)
if strcmp(mode,'step38_checkpoint'),p=cp;else,p='';end
end
function [free,all]=edge_inventory(T)
E=sort([T(:,[1 2]);T(:,[2 3]);T(:,[3 1])],2);
[U,~,ic]=unique(E,'rows');count=accumarray(ic,1);
all=[U,count];free=all(count==1,:);
end
function mask=face_edge_mask(E,up,lo,tip)
mask=(all(ismember(E,[up(:);tip]),2)|all(ismember(E,[lo(:);tip]),2));
end
function a=signed_area(P,T)
b=P(T(:,2),:)-P(T(:,1),:);c=P(T(:,3),:)-P(T(:,1),:);
a=.5*(b(:,1).*c(:,2)-b(:,2).*c(:,1));
end
function q=quality(P,T)
a=P(T(:,1),:);b=P(T(:,2),:);c=P(T(:,3),:);
L=[vecnorm(b-c,2,2),vecnorm(a-c,2,2),vecnorm(a-b,2,2)];
ang=zeros(size(L));
for k=1:3
    j=mod(k,3)+1;h=mod(k+1,3)+1;
    ang(:,k)=acosd(max(-1,min(1,(L(:,j).^2+L(:,h).^2-L(:,k).^2)./ ...
        (2*L(:,j).*L(:,h)))));
end
area=signed_area(P,T);
q=struct('area',area,'longest',max(L,[],2),'minAngle',min(ang,[],2), ...
    'shape',4*sqrt(3)*area./sum(L.^2,2));
end
function d=min_segment_distance(z,P,E)
a=P(E(:,1),:);b=P(E(:,2),:);v=b-a;
t=sum((z-a).*v,2)./sum(v.^2,2);t=max(0,min(1,t));
d=min(vecnorm(a+t.*v-z,2,2));
end
function cr=refresh_crack_metadata(cr,P,e,a0)
cr.upperS=(P(cr.upperNodes,:)-cr.Pmid(1,:))*e.'/a0;
cr.lowerS=(P(cr.lowerNodes,:)-cr.Pmid(1,:))*e.'/a0;
[cr.upperS,ix]=sort(cr.upperS);cr.upperNodes=cr.upperNodes(ix);
[cr.lowerS,ix]=sort(cr.lowerS);cr.lowerNodes=cr.lowerNodes(ix);
cr.nUpper=numel(cr.upperNodes);cr.nLower=numel(cr.lowerNodes);
cr.sameCount=cr.nUpper==cr.nLower;
cr.upperTarget=P(cr.upperNodes,:);cr.lowerTarget=P(cr.lowerNodes,:);
cr.lowerMatchForUpper=(1:cr.nUpper)';
end
function ids=q_support(mesh,cr,ri,ro,skipConstant)
P=mesh.coord;T=mesh.connect;tip=cr.Pmid(end,:);
e=diff(cr.Pmid)/norm(diff(cr.Pmid));R=[e(:),[-e(2);e(1)]];
r=vecnorm((P-tip)*R,2,2);q=ones(size(r));q(r>=ro)=0;
mid=r>ri&r<ro;q(mid)=(ro-r(mid))/(ro-ri);
xip=rule16();use=false(size(T,1),1);
for k=1:size(T,1)
    X=P(T(k,:),:);qe=q(T(k,:));
    if skipConstant&&all(qe==qe(1)),continue,end
    if all(qe==0),continue,end
    c=mean(X(1:3,:),1);
    if norm((c-tip)*R)>ro+max(vecnorm(X-c,2,2)),continue,end
    for j=1:16
        [N,d]=shape(xip(:,j));J=d*X;detJ=det(J);
        assert(detJ>0,'step62:InvalidT6','Nonpositive T6 Jacobian.');
        if norm(N.'*X-tip)<1e-12,continue,end
        grad=R.'*((J\d)*qe);
        if norm(grad)>1e-14,use(k)=true;break,end
    end
end
ids=find(use);
end
function xi=rule16()
xi=zeros(2,16);xi(:,1)=[1/3;1/3];
a=.170569307751760;b=.658861384496480;xi(:,2:4)=[a a b;a b a];
a=.050547228317031;b=.898905543365938;xi(:,5:7)=[a a b;a b a];
a=.459292588292723;b=.081414823414554;xi(:,8:10)=[a a b;a b a];
a=.263112829634638;b=.728492392955404;c=.008394777409958;
xi(:,11:16)=[a a b b c c;b c a c a b];
end
function [N,d]=shape(xi)
a=xi(1);b=xi(2);c=1-a-b;
N=[a*(2*a-1);b*(2*b-1);c*(2*c-1);4*a*b;4*b*c;4*c*a];
d=[4*a-1,0,-(4*c-1),4*b,-4*b,4*(c-a); ...
   0,4*b-1,-(4*c-1),4*a,4*(c-b),-4*a];
end
function [err3,err6]=pair_errors(mesh,ids,map,nU,rows,tip,R)
n=size(rows,1)/2;U=mesh.connect3(rows(1:n),:);
L=mesh.connect3(rows(n+1:end),[1 3 2]);
globalMap=zeros(size(mesh.coord3,1),1);globalMap(ids(1:nU))=ids(map);
assert(isequal(globalMap(U),L),'step62:PairConnectivity','T3 connectivity not mirrored.');
x=(mesh.coord3-tip)*R;
err3=max(vecnorm(x(L(:),:)-[x(U(:),1),-x(U(:),2)],2,2));
U6=mesh.connect(rows(1:n),:);L6=mesh.connect(rows(n+1:end),[1 3 2 6 5 4]);
x=(mesh.coord-tip)*R;
err6=max(vecnorm(x(L6(:),:)-[x(U6(:),1),-x(U6(:),2)],2,2));
% Verify a consistent involution for shared midsides across ALL paired rows.
pairs=[U6(:),L6(:)];[u,~,g]=unique(pairs(:,1));
mins=accumarray(g,pairs(:,2),[],@min);maxs=accumarray(g,pairs(:,2),[],@max);
assert(all(mins==maxs)&&numel(unique(mins))==numel(u), ...
    'step62:PairT6Map','Inconsistent/non-bijective T6 reflection map.');
end
function [out,counts,r]=topology_audit(mesh,cr,ids,axisIDs,map,Z)
up=cr.upperNodes;lo=cr.lowerNodes;shared=intersect(up,lo);
out=struct('crackDistinct',isequal(shared,cr.tipNode), ...
    'intactShared',true,'crackMidsDistinct',true,'intactMidsShared',true);
T=mesh.connect;maps=[1 2 4;2 3 5;3 1 6];
intact=axisIDs(Z(axisIDs,1)>=0);[~,ix]=sort(Z(intact,1));intact=intact(ix);
negative=axisIDs(Z(axisIDs,1)<=0);[~,ix]=sort(Z(negative,1));negative=negative(ix);
for k=1:numel(intact)-1
    E=ids(intact(k:k+1));m=find_mid(T,maps,E);
    out.intactShared=out.intactShared&&numel(m)==1&& ...
        ids(intact(k))==ids(map(intact(k)));
    out.intactMidsShared=out.intactMidsShared&&numel(m)==1;
end
for k=1:numel(negative)-1
    u=ids(negative(k:k+1));l=ids(map(negative(k:k+1)));
    mu=find_mid(T,maps,u);ml=find_mid(T,maps,l);
    out.crackMidsDistinct=out.crackMidsDistinct&& ...
        numel(mu)==1&&numel(ml)==1&&mu~=ml;
end
% Geometry-only native sampling: zero field has NO physical interpretation.
matDummy=struct('E',1,'nu',.3,'ps',1);
[r,~,diag]=native_COD_audit(mesh,zeros(2*size(mesh.coord,1),1), ...
    matDummy,cr,2);
counts=[diag.nUpper,diag.nLower];
out.faceGridMismatch=diag.gridMismatch;
out.crackDistinct=out.crackDistinct&&counts(1)==counts(2)&&diag.gridMismatch<1e-12;
end
function m=find_mid(T,maps,E)
m=[];
for k=1:3
    hit=all(sort(T(:,maps(k,1:2)),2)==sort(E(:).'),2);
    m=[m;T(hit,maps(k,3))]; %#ok<AGROW>
end
m=unique(m);
end
function [minJ,err]=jacobian_audit(P,T)
err=0;
for edge=[1 2 4;2 3 5;3 1 6].'
    err=max(err,max(vecnorm(P(T(:,edge(3)),:)- ...
        .5*(P(T(:,edge(1)),:)+P(T(:,edge(2)),:)),2,2)));
end
minJ=Inf;xi=rule16();
for j=1:16
    [~,d]=shape(xi(:,j));
    j11=sum(reshape(P(T(:),1),size(T)).*d(1,:),2);
    j12=sum(reshape(P(T(:),2),size(T)).*d(1,:),2);
    j21=sum(reshape(P(T(:),1),size(T)).*d(2,:),2);
    j22=sum(reshape(P(T(:),2),size(T)).*d(2,:),2);
    minJ=min(minJ,min(j11.*j22-j12.*j21));
end
end
function pass=boundary_same(P,T,Q,S,cr,crNew)
[E,~]=edge_inventory(T);[F,~]=edge_inventory(S);
E=E(~face_edge_mask(E(:,1:2),cr.upperNodes,cr.lowerNodes,cr.tipNode),1:2);
F=F(~face_edge_mask(F(:,1:2),crNew.upperNodes,crNew.lowerNodes,crNew.tipNode),1:2);
A=[P(E(:,1),:),P(E(:,2),:)];B=[Q(F(:,1),:),Q(F(:,2),:)];
A=canonical_segments(A);B=canonical_segments(B);
% Exact saved polygon coverage permits subdivision of its straight segments.
% Original vertices are separately required to remain bitwise unchanged.
assigned=zeros(size(B,1),1);intervals=zeros(size(B,1),2);pass=true;
for j=1:size(B,1)
    v=A(:,3:4)-A(:,1:2);ell=sum(v.^2,2);
    u=(B(j,1:2)-A(:,1:2));w=(B(j,3:4)-A(:,1:2));
    t=sum(u.*v,2)./ell;s=sum(w.*v,2)./ell;
    match=vecnorm(u-t.*v,2,2)<1e-12&vecnorm(w-s.*v,2,2)<1e-12& ...
        t>=-1e-9&t<=1+1e-9&s>=-1e-9&s<=1+1e-9;
    k=find(match,1);if isempty(k),pass=false;return,end
    assigned(j)=k;intervals(j,:)=sort([t(k),s(k)]);
end
for j=1:size(A,1)
    spans=sortrows(intervals(assigned==j,:));
    if isempty(spans)||abs(spans(1,1))>1e-9||abs(spans(end,2)-1)>1e-9|| ...
            any(abs(spans(2:end,1)-spans(1:end-1,2))>1e-9)
        pass=false;return
    end
end
end
function A=canonical_segments(A)
flip=A(:,1)>A(:,3)|(A(:,1)==A(:,3)&A(:,2)>A(:,4));
A(flip,:)=A(flip,[3 4 1 2]);
end
function pass=global_planar_audit(P,T,Q,S,X,patchRows,extRows,inner,cr)
[~,inventory]=edge_inventory(S);
pass=all(inventory(:,3)<=2)&&size(unique(sort(S,2),'rows'),1)==size(S,1)&& ...
    all(vecnorm(Q(inventory(:,1),:)-Q(inventory(:,2),:),2,2)>1e-14)&& ...
    numel(unique(S(:)))==size(Q,1);
areaOld=sum(signed_area(P,T));areaNew=sum(signed_area(Q,S));
pass=pass&&abs(areaOld-areaNew)<1e-12*areaOld;
% Every structured outer edge is shared by patch and exterior, except
% intentional coincident-but-distinct slit topology.
[free,~]=edge_inventory(S(patchRows,:));E=free(:,1:2);
rr=reshape(vecnorm(X(E(:),:),2,2),size(E));rp=max(vecnorm(inner,2,2));
E=E(all(rr>rp-1e-12,2),:);
[has,where]=ismember(sort(E,2),inventory(:,1:2),'rows');
pass=pass&&all(has)&&all(inventory(where(has),3)==2);
cent=(X(S(extRows,1),:)+X(S(extRows,2),:)+X(S(extRows,3),:))/3;
pass=pass&&~any(inpolygon(cent(:,1),cent(:,2),inner(:,1),inner(:,2)));
% Boundary equality is separately required. All other free edges must be
% geometrically valid classified crack edges, with no unclassified seam.
[free,~]=edge_inventory(S);E=free(:,1:2);
face=face_edge_mask(E,cr.upperNodes,cr.lowerNodes,cr.tipNode);
xx=reshape(X(E(face,:)',1),2,[])';yy=reshape(X(E(face,:)',2),2,[])';
pass=pass&&all(xx(:)<=1e-12)&&all(xx(:)>=-norm(diff(cr.Pmid))-1e-12)&& ...
    all(abs(yy(:))<1e-12);
end

function m=max_neighbor_size_ratio(T,L,affected)
n=size(T,1);E=sort([T(:,[1 2]);T(:,[2 3]);T(:,[3 1])],2);
which=repmat((1:n)',3,1);[~,~,g]=unique(E,'rows');
lo=accumarray(g,L(which),[],@min);hi=accumarray(g,L(which),[],@max);
touch=accumarray(g,ismember(which,affected),[],@max)>0;
m=max(hi(touch)./lo(touch));
end
function L=tip_edges(P,T,tip)
tt=T(any(T==tip,2),:);
neighbors=setdiff(unique(tt(:)),tip);
L=vecnorm(P(neighbors,:)-P(tip,:),2,2);
end
function tableOut=quality_comparison(old,mask,new,patch,exterior)
names={'Source affected';'Candidate paired';'Candidate exterior';'Candidate affected'};
q={old,new,new,new};ii={find(mask),patch,exterior,[patch;exterior]};rows=zeros(4,6);
for k=1:4
    v=q{k};i=ii{k};rows(k,:)=[numel(i),min(v.minAngle(i)), ...
    min(v.shape(i)),median(v.shape(i)),median(v.longest(i))*1e3, ...
    prctile(v.longest(i),90)*1e3];
end
tableOut=array2table(rows,'VariableNames',{'nT3','minAngle_deg','minShape', ...
    'medianShape','longestMedian_mm','longestP90_mm'});
tableOut.region=names;tableOut=movevars(tableOut,'region','Before',1);
end
function out=radial_comparison(X,T,Y,S,old,new,rp)
edges=[0 .2 .4 .8 1.2 1.6 2.4 3.2 4 4.8 5.2 rp*1e3 6.4 7];
edges=unique(sort(edges));rows=zeros(numel(edges)-1,10);
r=vecnorm((X(T(:,1),:)+X(T(:,2),:)+X(T(:,3),:))/3,2,2)*1e3;
s=vecnorm((Y(S(:,1),:)+Y(S(:,2),:)+Y(S(:,3),:))/3,2,2)*1e3;
for k=1:size(rows,1)
    a=r>=edges(k)&r<edges(k+1);b=s>=edges(k)&s<edges(k+1);
    rows(k,:)=[edges(k:k+1),nnz(a),nnz(b),median(old.longest(a))*1e3, ...
    median(new.longest(b))*1e3,prctile(old.longest(a),90)*1e3, ...
    prctile(new.longest(b),90)*1e3,min(old.minAngle(a)),min(new.minAngle(b))];
end
out=array2table(rows,'VariableNames',{'rMin_mm','rMax_mm','oldCount','newCount', ...
    'oldMedian_mm','newMedian_mm','oldP90_mm','newP90_mm', ...
    'oldMinAngle_deg','newMinAngle_deg'});
end
function syn=synthetic_controls(mesh,cr,mat,ri,ro)
fprintf('STRUCTURAL PASS: starting prescribed fields ONLY.\n');
tip=cr.Pmid(end,:);a0=norm(diff(cr.Pmid));e=diff(cr.Pmid)/a0;
R=[e(:),[-e(2);e(1)]];Z=(mesh.coord-tip)*R;
% Affine sanity at all sixteen points in every affected element, including
% both displacement components and analytical gradient reproduction.
u=[1+2*Z(:,1)-3*Z(:,2),-2+.5*Z(:,1)+4*Z(:,2)];
xi=rule16();affineError=0;
for j=1:16
    [N,d]=shape(xi(:,j));
    for k=1:size(mesh.connect,1)
        ids=mesh.connect(k,:);x=Z(ids,:);v=u(ids,:);J=d*x;
        affineError=max(affineError,max(abs((N.'*v)- ...
            [1+2*(N.'*x(:,1))-3*(N.'*x(:,2)), ...
             -2+.5*(N.'*x(:,1))+4*(N.'*x(:,2))])));
        affineError=max(affineError,max(abs((J\d)*v-[2 .5;-3 4]),[],'all'));
    end
end
syn=struct('performed',true,'passed',false,'affineError',affineError);
if affineError>1e-8,fprintf('STOP: affine sanity failed.\n');return,end
[~,~,face]=native_COD_audit(mesh,zeros(2*size(mesh.coord,1),1),mat,cr,2,true);
up=find(face.faceSide==1);lo=find(face.faceSide==-1);
Z([up;lo],2)=0;
% Independent tiny mixed input: never use the physical Step38 K values.
cases=[1 0;0 1;1 1e-4];
rows=nan(3,7);domain=struct('r_inner',ri,'r_outer',ro);
for k=1:3
    Ulocal=exact_williams_displacement_audit(Z,cases(k,1),cases(k,2), ...
        mat.E,mat.nu,mat.ps,'UpperFaceIDs',up,'LowerFaceIDs',lo);
    v=reshape(Ulocal,2,[]).'*R.';U=zeros(2*size(Z,1),1);
    U(1:2:end)=v(:,1);U(2:2:end)=v(:,2);
    [ki,kii,aux]=SIF_LEFM_interaction_EDI(mesh,U,cr.Pmid,mat,domain, ...
        'UsePlaneStrain',mat.ps==1,'WeightFunction','fe_nodal', ...
        'QuadratureRule',16,'StoreGPDiagnostics',false);
    rows(k,:)=[cases(k,:),ki,kii,ki-cases(k,1),kii-cases(k,2),aux.nElem_used];
    fprintf('prescribed case %d: KI=%.12g KII=%+.12g\n',k,ki,kii);
    if k==1&&abs(kii/ki)>1e-10
        fprintf('STOP: prescribed pure-I leakage gate failed.\n');break
    end
    if k==2&&(kii<=0||abs(kii-1)>2e-4||abs(ki)>1e-10)
        fprintf('STOP: prescribed pure-II gate failed.\n');break
    end
end
syn.table=array2table(rows,'VariableNames',{'inputKI','inputKII','recoveredKI', ...
    'recoveredKII','deltaKI','deltaKII','literalEDIElementsUsed'});
M=rows(1:2,3:4).';syn.recoveryMatrix=M;
syn.matrixError=norm(M-eye(2),'fro');
syn.mixedRelativeKIIError=abs(rows(3,4)/cases(3,2)-1);
syn.superpositionError=norm(rows(3,3:4).'-M*cases(3,:).');
% Prospective verification criteria; never physical uncertainty estimates.
syn.passed=all(isfinite(rows(:)))&&abs(rows(1,4)/rows(1,3))<=1e-10&& ...
    abs(rows(2,3))<=1e-10&&syn.matrixError<=2e-4&& ...
    syn.mixedRelativeKIIError<=2e-4&&syn.superpositionError<=1e-10;
disp(syn.table);fprintf('Synthetic qualification passed=%d\n',syn.passed);
end
function files=plot_candidate(prefix,vis,X,T,Y,S,~,patchIDs,extIDs, ...
    oldSupport,newSupport,~,inner,~,ri,ro,cr,crNew,oq,nq,design,exterior)
f=figure('Visible',char(vis),'Color','w','Position',[30 30 800 600], ...
    'WindowStyle','normal','Renderer','painters', ...
    'DefaultAxesFontSize',7,'DefaultTextFontSize',8);
% Explicit subplot axes retain equal scales during offscreen export.
for j=1:4
    ax=subplot(3,2,j,'Parent',f);hold(ax,'on');
    if mod(j,2)==1,Z=X;C=T;support=oldSupport;else,Z=Y;C=S;support=newSupport;end
    patch(ax,'Faces',C,'Vertices',Z*1e3,'FaceColor','none', ...
        'EdgeColor',[.38 .44 .49],'LineWidth',.25);
    if mod(j,2)==0
        patch(ax,'Faces',S(patchIDs,:),'Vertices',Y*1e3, ...
            'FaceColor',[.64 .85 .96],'FaceAlpha',.25,'EdgeColor','none');
    end
    if j>2
        patch(ax,'Faces',C(support,:),'Vertices',Z*1e3, ...
            'FaceColor',[.97 .72 .25],'FaceAlpha',.4,'EdgeColor','none');
        plot(ax,1e3*inner([1:end,1],1),1e3*inner([1:end,1],2),'b-','LineWidth',1.5);
        th=linspace(0,2*pi,300);
        for r=[ri ro],plot(ax,r*1e3*cos(th),r*1e3*sin(th),'k--');end
        xlim(ax,[-9 7]);ylim(ax,[-7 7]);
    end
    axis(ax,'equal');
    if j>2,xlim(ax,[-9 7]);ylim(ax,[-7 7]);end
    pbaspect(ax,[diff(xlim(ax)),diff(ylim(ax)),1]);
    xlabel(ax,'local x_1 (mm)');ylabel(ax,'local x_2 (mm)');
    if mod(j,2)==1,title(ax,sprintf('Source baseline: %d triangles',size(T,1))); ...
    else,title(ax,sprintf('Candidate L%d: full graded remesh, %d triangles',design.level,size(S,1)));end
end
ax=subplot(3,2,5,'Parent',f);hold(ax,'on');
patch(ax,'Faces',S,'Vertices',Y*1e3,'FaceColor','none','EdgeColor',[.2 .35 .5],'LineWidth',.6);
plot(ax,1e3*Y(crNew.upperNodes,1),1e3*Y(crNew.upperNodes,2),'bo','MarkerSize',4);
plot(ax,1e3*Y(crNew.lowerNodes,1),1e3*Y(crNew.lowerNodes,2),'rx','MarkerSize',5);
axis(ax,'equal');xlim(ax,[-.25 .25]);ylim(ax,[-.25 .25]);
pbaspect(ax,[1 1 1]);
title(ax,'Six equilateral tip triangles; seven incident edges');
xlabel(ax,'local x_1 (mm)');ylabel(ax,'local x_2 (mm)');
ax=subplot(3,2,6,'Parent',f);hold(ax,'on');
r=vecnorm((X(T(:,1),:)+X(T(:,2),:)+X(T(:,3),:))/3,2,2)*1e3;
s=vecnorm((Y(S(:,1),:)+Y(S(:,2),:)+Y(S(:,3),:))/3,2,2)*1e3;
scatter(ax,r,oq.longest*1e3,2,[.65 .65 .65],'.');
scatter(ax,s(patchIDs),nq.longest(patchIDs)*1e3,4,[0 .4 .8],'.');
scatter(ax,s(extIDs),nq.longest(extIDs)*1e3,3,[.55 .2 .65],'.');
plot(ax,design.ringTable.radius_mm,design.ringTable.targetH_mm,'k-','LineWidth',1.5);
plot(ax,exterior.ringTable.radius_m*1e3,exterior.ringTable.targetH_m*1e3,'k-','LineWidth',1.5);
xlim(ax,[0 40]);ylim(ax,[0 4]);grid(ax,'on');
title(ax,'Gray source; blue paired patch; purple exterior; black target h');
xlabel(ax,'centroid radius (mm)');ylabel(ax,'edge length / target spacing (mm)');
sgtitle(f,sprintf('Structured graded family L%d; exact physical boundary; no physical FEM',design.level), ...
    'FontSize',10);
files=struct('overviewPNG',[prefix '_overview.png'],'tipPNG',[prefix '_tip.png']);
drawnow;exportgraphics(f,files.overviewPNG,'Resolution',300);close(f);
f=figure('Visible',char(vis),'Color','w','Position',[30 30 800 450], ...
    'WindowStyle','normal','Renderer','painters');
% Separate equal-scale tip panels.
for j=1:2
    ax=subplot(1,2,j,'Parent',f);hold(ax,'on');
    if j==1,Z=X;C=T;cc=cr;else,Z=Y;C=S;cc=crNew;end
    patch(ax,'Faces',C,'Vertices',Z*1e3,'FaceColor','none', ...
        'EdgeColor',[.25 .4 .55],'LineWidth',.6);
    plot(ax,Z(cc.upperNodes,1)*1e3,Z(cc.upperNodes,2)*1e3,'bo','MarkerSize',4);
    plot(ax,Z(cc.lowerNodes,1)*1e3,Z(cc.lowerNodes,2)*1e3,'rx','MarkerSize',5);
    axis(ax,'equal');xlim(ax,[-.22 .22]);ylim(ax,[-.22 .22]);grid(ax,'on');
    pbaspect(ax,[1 1 1]);
    xlabel(ax,'local x_1 (mm)');ylabel(ax,'local x_2 (mm)');
    if j==1,title(ax,'Source tip fan');else,title(ax,sprintf('Candidate L%d six-triangle fan',design.level));end
end
drawnow;exportgraphics(f,files.tipPNG,'Resolution',220);close(f);
end

function save_report(prefix,O62)
folder=fileparts(prefix);if ~isempty(folder)&&exist(folder,'dir')~=7,mkdir(folder);end
save([prefix '_small_data.mat'],'O62','-v7');
names={'summary','extentTable','qualityTable','radialTable','samplingTable', ...
    'tipEdgeTable','ringTable','exteriorRingTable','exteriorRefinementTable'};
for k=1:numel(names)
    if isfield(O62,names{k}),writetable(O62.(names{k}), ...
        [prefix '_' names{k} '.csv']);end
end
if isfield(O62,'synthetic')&&isfield(O62.synthetic,'table')
    writetable(O62.synthetic.table,[prefix '_syntheticTable.csv']);
end
end
function assert_audit_branch(root)
branch=current_branch(root);
ok=strcmp(branch,'sif-asymmetric-mesh-audit')||startsWith(branch,'audit/step62')|| ...
    strcmp(branch,'audit/step63-calibrated-physical-solve')|| ...
    strcmp(branch,'audit/step63r-recover-lost-physical-field')|| ...
    strcmp(branch,'audit/step65-level1-mesh-qualification')|| ...
    strcmp(branch,'audit/step66-solver-memory-preflight')|| ...
    strcmp(branch,'audit/step67-level0-iterative-solver-qualification')|| ...
    strcmp(branch,'audit/step67a-level0-sgs-solver-qualification');
assert(ok,'step62:Branch', ...
    'Step62 mesh-only construction is not authorized on the current branch.');
end
function branch=current_branch(root)
[status,branch]=system(sprintf('git -C "%s" branch --show-current',root));
assert(status==0,'step62:GitBranch','Cannot resolve current Git branch.');
branch=strtrim(branch);
end
function s=git_head(root)
[status,s]=system(sprintf('git -C "%s" rev-parse HEAD',root));
assert(status==0);s=strtrim(s);
end
function h=sha256(path)
if ~endsWith(path,'.m')&&exist([path '.m'],'file')==2,path=[path '.m'];end
f=fopen(path,'rb');assert(f>=0);guard=onCleanup(@()fclose(f));
md=java.security.MessageDigest.getInstance('SHA-256');
while ~feof(f),b=fread(f,1024*1024,'*uint8');md.update(typecast(b,'int8'));end
h=lower(reshape(dec2hex(typecast(md.digest(),'uint8'),2).',1,[]));
clear guard
end
