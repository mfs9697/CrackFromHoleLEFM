function O61=main_step61_local_paired_patch(varargin)
%MAIN_STEP61_LOCAL_PAIRED_PATCH One mesh-only asymmetric Step38 candidate.
% NEVER loads physical U, assembles stiffness, or solves equilibrium.
% Uses the saved collapsed geometry. One upper graded triangulation is
% reflected; only a bounded collar connects it to the unchanged exterior.
% Prescribed-field EDI is permitted ONLY after all structural gates pass.
root=fileparts(fileparts(fileparts(mfilename('fullpath'))));
ip=inputParser;
addParameter(ip,'CheckpointFile',fullfile(root,'verification', ...
    'step38_tip_refined_solved.mat'),@(x)ischar(x)||isstring(x));
addParameter(ip,'SavePrefix',fullfile(root,'verification', ...
    'step61_local_paired_patch'),@(x)ischar(x)||isstring(x));
addParameter(ip,'Visible','off',@(x)ischar(x)||isstring(x));
addParameter(ip,'RunSynthetic',true,@(x)islogical(x)&&isscalar(x));
parse(ip,varargin{:}); opt=ip.Results;
addpath(genpath(root));
assert_audit_branch(root);
cp=char(opt.CheckpointFile); prefix=char(opt.SavePrefix);
% Intentionally selective: U, forces and stiffness are NEVER loaded.
s=load(cp,'mesh','crack','baseline','actualTip','a0','mat');
P=s.mesh.coord3; T=s.mesh.connect3; cr=s.crack; mat=s.mat;
a0=s.a0; hTip=s.actualTip; ri=s.baseline.rInner; ro=.65*a0;
assert(size(P,1)==20164&&size(T,1)==39441&& ...
    size(s.mesh.coord,1)==79769,'step61:SourceCount','Wrong Step38 mesh.');
assert(abs(a0-.008)<1e-12&&abs(ri-.0008)<1e-12&& ...
    abs(hTip-5.4024650785e-5)<1e-12,'step61:SourceScale','Wrong source scales.');
assert(norm(cr.Pmid(1,:)-[.199988872196,-.0208170339049])<2e-12&& ...
    norm(cr.Pmid(end,:)-[.207985904782,-.0210349096129])<2e-12, ...
    'step61:SourceCrack','Wrong saved crack endpoints.');
tip=cr.Pmid(end,:); e1=diff(cr.Pmid)/a0;
R=[e1(:),[-e1(2);e1(1)]]; X=(P-tip)*R;
fprintf('STEP61: saved asymmetric mesh, no physical U loaded.\n');
% Step60 convention and literal production participation are both audited.
oldSupport=q_support(s.mesh,cr,ri,ro,true);
oldLiteralSupport=q_support(s.mesh,cr,ri,ro,false);
assert(numel(oldSupport)==14215,'step61:SourceSupport','Step60 support changed.');
supportNodes=unique(s.mesh.connect(oldSupport,:));
maxSupport=max(vecnorm((s.mesh.coord(supportNodes,:)-tip)*R,2,2));
oldQuality=quality(X,T);
[oldFree,~]=edge_inventory(T);
oldFree=oldFree(:,1:2);
upOld=unique(cr.upperNodes(:)); loOld=unique(cr.lowerNodes(:));
isFace=face_edge_mask(oldFree,upOld,loOld,cr.tipNode);
physicalOther=oldFree(~isFace,:);
physicalClearance=min_segment_distance([0 0],X,physicalOther);

% First choose the BOUNDARY from measured geometry, before triangulation.
% 0.15 mm support clearance; collar >=0.45 mm; no noncrack boundary touched.
radii=(5.4:.1:6.0)'*1e-3;
extentRows=zeros(numel(radii),7); masks=cell(numel(radii),1);
for k=1:numel(radii)
    rp=radii(k); cut=rp+.65e-3;
    mask=any(reshape(vecnorm(X(T(:),:),2,2),size(T))<cut,2);
    [boundary,~,~]=cavity_boundary(T,mask,cr);
    rb=vecnorm(X(unique(boundary),:),2,2);
    minCollar=min_segment_distance([0 0],X,boundary)-rp;
    clearPhysical=~any(ismember(sort(boundary,2), ...
        sort(physicalOther,2),'rows'));
    eligible=rp-maxSupport>=.15e-3 && minCollar>=.45e-3 && ...
        clearPhysical && max(rb)<physicalClearance-.4e-3;
    extentRows(k,:)=[rp,maxSupport,rp-maxSupport,minCollar, ...
        max(rb),nnz(mask),eligible]; masks{k}=mask;
end
extentTable=array2table(extentRows,'VariableNames', ...
    {'pairedRadius_m','oldSupportMaxRadius_m','supportClearance_m', ...
     'minimumCollar_m','cavityMaxRadius_m','replacedT3','eligible'});
disp(extentTable);
pick=find(extentRows(:,7)>0,1);
if isempty(pick)
    O61=struct('feasible',false,'extentTable',extentTable, ...
        'reason','No support-enclosing boundary with a safe conforming collar.', ...
        'noPhysicalU',true,'noFEM',true);
    save_report(prefix,O61); return
end
rp=radii(pick); removed=masks{pick}; outside=~removed;
[outerEdges,outerIDs,outerXY]=cavity_boundary(T,removed,cr,X);
% Order the geometric outer boundary after merging ONLY coincident slit ends.
fprintf('Chosen paired radius %.4f mm, replacing %d old triangles.\n', ...
    rp*1e3,nnz(removed));
[Zpatch,Tpatch,mirrorMap,nUpper,ringRadii,axisUp]=graded_patch(rp,hTip);
[Zcollar,Tcollar,collarMeta]=transition_collar( ...
    X,T,removed,outerXY,Zpatch,Tpatch,rp);

% Map the new geometry onto existing interface IDs. Retained old vertices
% stay bitwise unchanged; no coordinate welding across the crack.
Q=P; pairedIDs=zeros(size(Zpatch,1),1);
for j=1:size(Zpatch,1)
    Q(end+1,:)=tip+Zpatch(j,:)*R.'; pairedIDs(j)=size(Q,1); %#ok<AGROW>
end
patchTriangles=pairedIDs(Tpatch);
collarIDs=zeros(size(Zcollar,1),1);
for j=1:size(Zcollar,1)
    z=Zcollar(j,:);
    side=collarMeta.nodeSide(j);
    candidates=outerIDs;
    if z(1)<0&&abs(z(2))<1e-12
        if side>=0,candidates=setdiff(candidates,loOld); ...
        else,candidates=setdiff(candidates,upOld);end
    end
    d=vecnorm(X(candidates,:)-z,2,2); [dd,ii]=min(d);
    if dd<1e-12
        collarIDs(j)=candidates(ii); continue
    end
    d=vecnorm(Zpatch-z,2,2);
    if z(1)<0&&abs(z(2))<1e-12
        upperAxis=axisUp(Zpatch(axisUp,1)<0);
        lowerAxis=mirrorMap(upperAxis);
        if side>=0,d(lowerAxis)=Inf;else,d(upperAxis)=Inf;end
    end
    [dd,ii]=min(d);
    if dd<1e-12,collarIDs(j)=pairedIDs(ii);continue,end
    Q(end+1,:)=tip+z*R.'; collarIDs(j)=size(Q,1); %#ok<AGROW>
end
Tc=[T(outside,:);patchTriangles;collarIDs(Tcollar)];
% Compact unused old interior nodes; retain exact exterior coordinates/IDs
% through a saved original-to-candidate map.
used=unique(Tc(:)); oldToNew=zeros(size(Q,1),1);
oldToNew(used)=(1:numel(used))'; Pc=Q(used,:); Tc=oldToNew(Tc);
pairedIDs=oldToNew(pairedIDs); collarIDs=oldToNew(collarIDs);
nExterior=nnz(outside); nPatch=size(Tpatch,1);
patchRows=nExterior+(1:nPatch)';
collarRows=(nExterior+nPatch+1:size(Tc,1))';
Xc=(Pc-tip)*R;
% Derive crack IDs from topology, never from coincident coordinates alone.
crNew=cr; crNew.tipNode=pairedIDs(1);
newUpper=pairedIDs(axisUp(Zpatch(axisUp,1)<-1e-12));
newLower=pairedIDs(mirrorMap(axisUp(Zpatch(axisUp,1)<-1e-12)));
onCollar=abs(Zcollar(:,2))<1e-12 & Zcollar(:,1)<0;
newUpper=[newUpper;collarIDs(onCollar&collarMeta.nodeSide>=0)];
newLower=[newLower;collarIDs(onCollar&collarMeta.nodeSide<0)];
keptUpper=oldToNew(upOld); keptLower=oldToNew(loOld);
crNew.upperNodes=unique([keptUpper(keptUpper>0);newUpper;crNew.tipNode]);
crNew.lowerNodes=unique([keptLower(keptLower>0);newLower;crNew.tipNode]);
crNew=refresh_crack_metadata(crNew,Pc,e1,a0);
[P6,T6]=T3toT6_fast(Pc,Tc);
candidateMesh=struct('coord3',Pc,'connect3',Tc,'coord',P6,'connect',T6);
newSupport=q_support(candidateMesh,crNew,ri,ro,true);
newLiteralSupport=q_support(candidateMesh,crNew,ri,ro,false);
newQuality=quality(Xc,Tc);
gates=struct();
gates.physicalGeometryUnchanged=boundary_same(P,T,Pc,Tc,cr,crNew) && ...
    isequal(cr.Pmid,crNew.Pmid);
gates.exteriorUnchanged=isequal(Tc(1:nExterior,:),oldToNew(T(outside,:)))&& ...
    isequal(Pc(oldToNew(unique(T(outside,:))),:),P(unique(T(outside,:)),:));
[pair3,pair6]=pair_errors(candidateMesh,pairedIDs,mirrorMap,nUpper, ...
    patchRows,tip,R);
gates.completeT3Pairing=pair3<1e-12;
gates.completeT6Pairing=pair6<1e-12;
[topology,faceCounts,nativeR]=topology_audit(candidateMesh,crNew, ...
    pairedIDs,axisUp,mirrorMap,Zpatch);
gates.crackFacesDistinct=topology.crackDistinct;
gates.intactLigamentShared=topology.intactShared;
gates.crackT6MidsidesDistinct=topology.crackMidsDistinct;
gates.intactT6MidsidesShared=topology.intactMidsShared;
gates.positiveT3Areas=all(newQuality.area>0);
[minJ,midsError]=jacobian_audit(P6,T6);
gates.positiveT6Jacobians=minJ>0&&midsError<1e-12;
gates.noHangingGapsOverlaps=planar_audit(P,T,Pc,Tc,removed, ...
    oldToNew,outerEdges,Xc,patchRows,collarRows,outerXY,Zpatch,Tpatch);
gates.originalSupportInsidePairedRegion= ...
    all(inpolygon((s.mesh.coord(supportNodes,:)-tip)*e1.', ...
    (s.mesh.coord(supportNodes,:)-tip)*R(:,2), ...
    collarMeta.innerPolygon(:,1),collarMeta.innerPolygon(:,2)));
gates.candidateSupportInsidePairedRegion=all(ismember(newSupport,patchRows))&& ...
    all(ismember(newLiteralSupport,patchRows));
gates.collarOutsideSupport=~any(ismember(newSupport,collarRows))&& ...
    ~any(ismember(newLiteralSupport,collarRows));
windows=[.04 .20;.04 .30;.08 .30;.12 .30];
sampleN=zeros(4,1);
for k=1:4,sampleN(k)=nnz(nativeR/a0>=windows(k,1)& ...
    nativeR/a0<=windows(k,2));end
gates.nativeSamplingAdequate=all(sampleN>=12);
newTip=tip_edges(Pc,Tc,crNew.tipNode);
oldTip=tip_edges(P,T,cr.tipNode);
tipEdgeTable=table({'Step38';'Candidate'},[numel(oldTip);numel(newTip)], ...
    [min(oldTip);min(newTip)]*1e3,[median(oldTip);median(newTip)]*1e3, ...
    [prctile(oldTip,90);prctile(newTip,90)]*1e3, ...
    [max(oldTip);max(newTip)]*1e3,'VariableNames', ...
    {'mesh','incidentEdges','min_mm','median_mm','p90_mm','max_mm'});
gates.tipResolutionMatched=abs(median(newTip)/hTip-1)<=.05;
% Engineering gates, not physical error bars.
affected=[patchRows;collarRows];
gates.minimumAngle20deg=min(newQuality.minAngle(affected))>=20;
gates.noAbruptGrading=max_neighbor_size_ratio(Tc,newQuality.longest,affected)<=2.5;
gates.structuralPass=all(structfun(@(x)logical(x),gates));
summary=table(size(P,1),size(T,1),size(s.mesh.coord,1), ...
    size(Pc,1),size(Tc,1),size(P6,1),nnz(removed),nnz(outside), ...
    rp*1e3,min(vecnorm(outerXY,2,2))*1e3, ...
    max(vecnorm(outerXY,2,2))*1e3,median(oldTip)*1e3,median(newTip)*1e3, ...
    numel(oldSupport),numel(newSupport), ...
    min(newQuality.minAngle(patchRows)),min(newQuality.minAngle(collarRows)), ...
    pair3,pair6,'VariableNames',{'oldT3Nodes','oldT3','oldT6Nodes', ...
    'newT3Nodes','newT3','newT6Nodes','replacedOldT3','retainedOldT3', ...
    'pairedRadius_mm','cavityMinVertexRadius_mm','cavityMaxRadius_mm', ...
    'oldTipMedian_mm','newTipMedian_mm','oldPrimarySupport','newPrimarySupport', ...
    'patchMinAngle_deg','collarMinAngle_deg','pairErrorT3_m','pairErrorT6_m'});
qualityTable=quality_comparison(oldQuality,removed,newQuality,patchRows,collarRows);
radialTable=radial_comparison(X,T,Xc,Tc,oldQuality,newQuality,rp);
samplingTable=table(windows(:,1),windows(:,2),sampleN, ...
    'VariableNames',{'lower_r_over_a0','upper_r_over_a0','nativePoints'});
disp(summary);disp(qualityTable);disp(samplingTable);disp(gates);
files=plot_candidate(prefix,opt.Visible,X,T,Xc,Tc,removed,patchRows, ...
    collarRows,oldSupport,newSupport,outerXY,collarMeta.innerPolygon,rp,ri,ro, ...
    cr,crNew,oldQuality,newQuality);
provenance=struct('checkpointPath',cp,'checkpointSHA256',sha256(cp), ...
    'branch','sif-asymmetric-mesh-audit','sourceCommit',git_head(root), ...
    'driverSHA256',sha256(mfilename('fullpath')), ...
    'physicalUWasLoaded',false,'sourceCrack',cr.Pmid, ...
    'sourceMaterial',mat,'baselineOuterRatios',s.baseline.rOuterOverA0, ...
    'rInner',ri,'rOuter',ro,'upperGeneratedOnce',true, ...
    'outsideTrianglesExact',true,'physicalBCsAndLoadsUnmodified',true, ...
    'matlabVersion',version,'supportDefinition','Step60 plus literal EDI audit');
O61=struct('feasible',true,'summary',summary,'extentTable',extentTable, ...
    'qualityTable',qualityTable,'radialTable',radialTable, ...
    'samplingTable',samplingTable,'gates',gates,'files',files, ...
    'tipEdgeTable',tipEdgeTable,'originalTipEdgeDistribution_m',oldTip, ...
    'targetTipEdge_m',hTip, ...
    'topology',topology,'minT6Jacobian',minJ, ...
    'maxNeighborSizeRatio',max_neighbor_size_ratio(Tc,newQuality.longest,affected), ...
    'originalLiteralSupportCount',numel(oldLiteralSupport), ...
    'candidateLiteralSupportCount',numel(newLiteralSupport), ...
    'minimumNoncrackBoundaryDistance_m',physicalClearance, ...
    'minimumCollarWidth_m',collarMeta.minimumCollar_m, ...
    'excludedZeroInteriorHullSlivers', ...
        collarMeta.excludedZeroInteriorHullSlivers, ...
    'tipEdgeDistribution_m',newTip,'provenance',provenance, ...
    'synthetic',struct('performed',false,'passed',false), ...
    'readyForOneAsymmetricFEMProposal',false, ...
    'noPhysicalU',true,'noFEM',true,'noStiffness',true);
% Save diagnostic reports even on rejection; only structurally accepted
% candidates receive the candidate MAT file.
save_report(prefix,O61);
if ~gates.structuralPass
    fprintf('STOP: structural/design gate failed. No synthetic EDI, no FEM.\n');
    return
end
candidate=struct('p',Pc,'t',Tc,'crack',crNew,'mat',mat, ...
    'pairedElementIDs',patchRows,'transitionElementIDs',collarRows, ...
    'primarySupportElementIDs',newSupport, ...
    'literalPrimarySupportElementIDs',newLiteralSupport, ...
    'pairedNodeMap',pairedIDs,'localMirrorMap',mirrorMap, ...
    'pairedBoundaryLocal',collarMeta.innerPolygon,'cavityBoundaryLocal',outerXY, ...
    'replacedOriginalElementIDs',find(removed), ...
    'retainedOriginalElementIDs',find(outside), ...
    'originalNodeToCandidate',oldToNew(1:size(P,1)), ...
    'ringRadii_m',ringRadii,'provenance',provenance,'gates',gates);
candidateFile=[prefix '_candidate_T3.mat'];
candidate.scientificallyReadyForFEMProposal=false;
save(candidateFile,'candidate','-v7'); O61.files.candidateMAT=candidateFile;
if opt.RunSynthetic
    O61.synthetic=synthetic_controls(candidateMesh,crNew,mat,ri,ro);
    candidate.scientificallyReadyForFEMProposal=O61.synthetic.passed;
    candidate.synthetic=O61.synthetic;
    save(candidateFile,'candidate','-v7');
end
O61.readyForOneAsymmetricFEMProposal=gates.structuralPass&&O61.synthetic.passed;
save_report(prefix,O61);
fprintf('Step61 finished: proposal ready=%d; zero FEM solves.\n', ...
    O61.readyForOneAsymmetricFEMProposal);
end

function [Z,T,map,nU,rings,axisIDs]=graded_patch(rp,hTip)
% Smooth radial sizing; no attempt to copy nested red/green complexity.
rings=hTip; r=hTip;
while r<rp
    h=target_spacing(r);
    dr=.86*h;
    if rp-r<1.5*dr,r=rp;else,r=r+dr;end
    rings(end+1)=r; %#ok<AGROW>
end
U=[0 0];
for k=1:numel(rings)
    r=rings(k);
    if k==1,n=6;else
        h=target_spacing(r);
        n=ceil(pi*r/h);
    end
    th=(0:n)'*pi/n;
    z=r*[cos(th),sin(th)];z([1 end],2)=0;
    U=[U;z]; %#ok<AGROW>
end
axisIDs=find(U(:,2)==0);
[~,order]=sort(U(axisIDs,1)); ax=axisIDs(order);
outer=find(abs(vecnorm(U,2,2)-rp)<1e-12);
[~,order]=sort(atan2(U(outer,2),U(outer,1)));outer=outer(order);
C=[outer(1:end-1),outer(2:end);ax(1:end-1),ax(2:end)];
dt=delaunayTriangulation(U,C); Tu=dt.ConnectivityList;
cent=(U(Tu(:,1),:)+U(Tu(:,2),:)+U(Tu(:,3),:))/3;
Tu=Tu(inpolygon(cent(:,1),cent(:,2),U([outer;ax(2:end-1)],1), ...
    U([outer;ax(2:end-1)],2)),:);
Tu=orient(U,Tu);
nU=size(U,1);shared=find(U(:,2)==0&U(:,1)>=0);
nonShared=setdiff((1:nU)',shared);map=(1:nU)';
map(nonShared)=nU+(1:numel(nonShared))';
Z=[U;[U(nonShared,1),-U(nonShared,2)]];
T=[Tu;map(Tu(:,[1 3 2]))];
end

function h=target_spacing(r)
% Smooth monotone sizing chosen from Step38 edge statistics, not K values.
% Delay growth to retain ~0.07--0.08 mm edges at the 0.8 mm inner annulus.
h=.056e-3+.079e-3*(1-exp(-(max(0,r-.0005)/.0013)^2));
end

function [Z,T,info]=transition_collar(X,Told,removed,outer,Zp,Tp,rp)
% One constrained planar triangulation, between exact interfaces.
% Coincident crack coordinates are split by side AFTER triangulation.
[free,~]=edge_inventory(Tp);
% Only the outer circumference is a collar interface; remove slit edges.
edges=free(:,1:2);
rr=reshape(vecnorm(Zp(edges(:),:),2,2),size(edges));
edges=edges(all(rr>rp-1e-12,2),:);
xy=Zp(unique(edges(:)),:);
[~,ix]=unique(round(xy/1e-12),'rows','stable');xy=xy(ix,:);
[~,ix]=sort(atan2(xy(:,2),xy(:,1)));inner=xy(ix,:);
nO=size(outer,1);nI=size(inner,1);
Q=[outer;inner];
C=[(1:nO)',[2:nO 1]';nO+(1:nI)',nO+[2:nI 1]'];
minOuter=min_segment_distance([0 0],outer,[(1:nO)',[2:nO 1]']);
% Graded ring seeds kept away from the fixed irregular outer boundary.
for r=rp+.12e-3:.12e-3:max(vecnorm(outer,2,2))-.06e-3
    n=ceil(2*pi*r/.14e-3);th=(0:n-1)'*2*pi/n;
    seeds=r*[cos(th),sin(th)];
    seeds=seeds(inpolygon(seeds(:,1),seeds(:,2),outer(:,1),outer(:,2)),:);
    keep=false(size(seeds,1),1);
    for j=1:size(seeds,1)
        keep(j)=min_segment_distance(seeds(j,:),outer, ...
            [(1:nO)',[2:nO 1]'])>.075e-3 && ...
            abs(seeds(j,2))>.05e-3;
    end
    Q=[Q;seeds(keep,:)]; %#ok<AGROW>
end
% The crack is a prescribed straight constraint in the collar.
cutOuter=find(abs(outer(:,2))<1e-12&outer(:,1)<0);
assert(numel(cutOuter)==1,'step61:CavityCrack','Cavity must have one slit endpoint.');
cutInner=find(abs(inner(:,2))<1e-12&inner(:,1)<0)+nO;
x0=outer(cutOuter,1);
% Equal subdivision avoids a short last crack segment at the paired edge.
nCut=max(2,round((-rp-x0)/.12e-3));
xx=linspace(x0,-rp,nCut+1).';xx=xx(2:end-1);
cutIDs=[cutOuter;(size(Q,1)+(1:numel(xx)))';cutInner];
Q=[Q;xx,zeros(size(xx))];
C=[C;cutIDs(1:end-1),cutIDs(2:end)];
dt=delaunayTriangulation(Q,C); Q=dt.Points;
Tc=dt.ConnectivityList;
cent=(Q(Tc(:,1),:)+Q(Tc(:,2),:)+Q(Tc(:,3),:))/3;
use=inpolygon(cent(:,1),cent(:,2),outer(:,1),outer(:,2)) & ...
    ~inpolygon(cent(:,1),cent(:,2),inner(:,1),inner(:,2));
% The original cavity has collinear boundary subdivisions. Roundoff in
% the crack-frame transform can make Delaunay emit isolated hull slivers
% using three such vertices. They have zero geometric interior and are
% not collar cells. Exclude them using a dimensionless area test; retain
% every original interface vertex/segment and verify exact incidence below.
scale=max([sum((Q(Tc(:,1),:)-Q(Tc(:,2),:)).^2,2), ...
    sum((Q(Tc(:,2),:)-Q(Tc(:,3),:)).^2,2), ...
    sum((Q(Tc(:,3),:)-Q(Tc(:,1),:)).^2,2)],[],2);
degenerate=abs(signed_area(Q,Tc))<=1e-12*scale;
excludedHullSlivers=nnz(use&degenerate);
use=use&~degenerate;
Tc=orient(Q,Tc(use,:));
% Eight bounded Laplacian smoothing steps, constrained interfaces fixed.
fixed=unique(C(:));
for it=1:8
    E=unique(sort([Tc(:,[1 2]);Tc(:,[2 3]);Tc(:,[3 1])],2),'rows');
    A=sparse([E(:,1);E(:,2)],[E(:,2);E(:,1)],1,size(Q,1),size(Q,1));
    trial=Q; movable=setdiff(unique(Tc(:)),fixed);
    avg=(A*Q)./max(1,sum(A,2));
    trial(movable,:)=.5*Q(movable,:)+.5*avg(movable,:);
    areas=signed_area(trial,Tc);
    if all(areas>0),Q=trial;end
    % Restore Delaunay connectivity after moving nodes; keeping the initial
    % diagonals can leave thin cells even when every area stays positive.
    if mod(it,2)==0
        dt=delaunayTriangulation(Q,C);Q=dt.Points;Tc=dt.ConnectivityList;
        cent=(Q(Tc(:,1),:)+Q(Tc(:,2),:)+Q(Tc(:,3),:))/3;
        use=inpolygon(cent(:,1),cent(:,2),outer(:,1),outer(:,2)) & ...
            ~inpolygon(cent(:,1),cent(:,2),inner(:,1),inner(:,2));
        scale=max([sum((Q(Tc(:,1),:)-Q(Tc(:,2),:)).^2,2), ...
            sum((Q(Tc(:,2),:)-Q(Tc(:,3),:)).^2,2), ...
            sum((Q(Tc(:,3),:)-Q(Tc(:,1),:)).^2,2)],[],2);
        use=use&abs(signed_area(Q,Tc))>1e-12*scale;
        Tc=orient(Q,Tc(use,:));
    end
end
% Duplicate each negative-axis collar node for the lower crack face.
cut=find(abs(Q(:,2))<1e-12&Q(:,1)<0);
lower=(1:size(Q,1))';lower(cut)=size(Q,1)+(1:numel(cut))';
cent=(Q(Tc(:,1),:)+Q(Tc(:,2),:)+Q(Tc(:,3),:))/3;
bottom=cent(:,2)<0; Tc(bottom,:)=lower(Tc(bottom,:));
Z=[Q;Q(cut,:)]; nodeSide=ones(size(Z,1),1);
nodeSide(lower(cut))=-1;
T=Tc;info=struct('innerPolygon',inner,'nodeSide',nodeSide, ...
    'excludedZeroInteriorHullSlivers',excludedHullSlivers, ...
    'minimumCollar_m',minOuter-rp,'originalRemovedArea', ...
    sum(signed_area(X,Told(removed,:))));
end

function [edges,ids,poly]=cavity_boundary(T,mask,cr,X)
[free,~]=edge_inventory(T(mask,:));edges=free(:,1:2);
edges=edges(~face_edge_mask(edges,cr.upperNodes,cr.lowerNodes,cr.tipNode),:);
ids=unique(edges(:));poly=[];
if nargin<4,return,end
xy=X(ids,:);[~,first,grp]=unique(round(xy/1e-12),'rows','stable');
xy=xy(first,:);local=zeros(max(ids),1);local(ids)=grp;
E=local(edges); deg=accumarray(E(:),1,[size(xy,1),1]);
assert(all(deg==2),'step61:CavityBoundary','Cavity boundary is not one closed loop.');
order=zeros(size(xy,1),1);order(1)=1;prev=0;
for k=2:numel(order)
    adjacent=E(any(E==order(k-1),2),:);adjacent=adjacent(:);
    next=setdiff(adjacent,[order(k-1);prev]);
    assert(~isempty(next),'step61:CavityLoop','Broken cavity loop.');
    order(k)=next(1);prev=order(k-1);
end
assert(numel(unique(order))==numel(order),'step61:CavityHoles','Cavity has extra loops.');
poly=xy(order,:);
if sum(poly(:,1).*circshift(poly(:,2),-1)- ...
       poly(:,2).*circshift(poly(:,1),-1))<0,poly=flipud(poly);end
end

function [free,all]=edge_inventory(T)
E=sort([T(:,[1 2]);T(:,[2 3]);T(:,[3 1])],2);
[U,~,ic]=unique(E,'rows');count=accumarray(ic,1);
all=[U,count];free=all(count==1,:);
end
function mask=face_edge_mask(E,up,lo,tip)
mask=(all(ismember(E,[up(:);tip]),2)|all(ismember(E,[lo(:);tip]),2));
end
function T=orient(P,T)
a=signed_area(P,T);bad=a<0;T(bad,[2 3])=T(bad,[3 2]);
assert(all(signed_area(P,T)>0),'step61:Degenerate','Degenerate triangle.');
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
        assert(detJ>0,'step61:InvalidT6','Nonpositive T6 Jacobian.');
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
assert(isequal(globalMap(U),L),'step61:PairConnectivity','T3 connectivity not mirrored.');
x=(mesh.coord3-tip)*R;
err3=max(vecnorm(x(L(:),:)-[x(U(:),1),-x(U(:),2)],2,2));
U6=mesh.connect(rows(1:n),:);L6=mesh.connect(rows(n+1:end),[1 3 2 6 5 4]);
x=(mesh.coord-tip)*R;
err6=max(vecnorm(x(L6(:),:)-[x(U6(:),1),-x(U6(:),2)],2,2));
% Verify a consistent involution for shared midsides across ALL paired rows.
pairs=[U6(:),L6(:)];[u,~,g]=unique(pairs(:,1));
mins=accumarray(g,pairs(:,2),[],@min);maxs=accumarray(g,pairs(:,2),[],@max);
assert(all(mins==maxs)&&numel(unique(mins))==numel(u), ...
    'step61:PairT6Map','Inconsistent/non-bijective T6 reflection map.');
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
pass=isequal(sortrows(A),sortrows(B));
end
function A=canonical_segments(A)
flip=A(:,1)>A(:,3)|(A(:,1)==A(:,3)&A(:,2)>A(:,4));
A(flip,:)=A(flip,[3 4 1 2]);
end
function pass=planar_audit(P,T,Q,S,mask,map,E,X,patchRows,collarRows,outer,Zp,Tp)
[free,inventory]=edge_inventory(S);
% All old interface edges survive as a shared two-triangle edge.
mapped=sort(map(E),2);[present,where]=ismember(mapped,inventory(:,1:2),'rows');
pass=all(present)&&all(inventory(where(present),3)==2)&&all(inventory(:,3)<=2);
areaOld=sum(signed_area(P,T));areaNew=sum(signed_area(Q,S));
pass=pass&&abs(areaOld-areaNew)<1e-12*areaOld;
% Prove the new triangulations cover exactly the declared cavity with
% disjoint interiors: constrained collar excludes the constrained patch;
% each is a valid Delaunay embedding, and their areas exhaust old cavity.
cent=(X(S(collarRows,1),:)+X(S(collarRows,2),:)+X(S(collarRows,3),:))/3;
pass=pass&&all(inpolygon(cent(:,1),cent(:,2),outer(:,1),outer(:,2)));
newArea=sum(signed_area(Q,S([patchRows;collarRows],:)));
oldArea=sum(signed_area(P,T(mask,:)));
pass=pass&&abs(oldArea-newArea)<1e-12*oldArea;
% Nonmanifold duplicate triangles, zero-length edges and dangling nodes.
pass=pass&&size(unique(sort(S,2),'rows'),1)==size(S,1)&& ...
    all(vecnorm(Q(inventory(:,1),:)-Q(inventory(:,2),:),2,2)>1e-14)&& ...
    numel(unique(S(:)))==size(Q,1);
clear free Zp Tp
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
function tableOut=quality_comparison(old,mask,new,patch,collar)
names={'Step38 affected';'Candidate paired';'Candidate collar';'Candidate affected'};
q={old,new,new,new};ii={find(mask),patch,collar,[patch;collar]};rows=zeros(4,6);
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
function files=plot_candidate(prefix,vis,X,T,Y,S,mask,patchIDs,collarIDs, ...
    oldSupport,newSupport,outer,inner,rp,ri,ro,cr,crNew,oq,nq)
folder=fileparts(prefix);if ~isempty(folder)&&exist(folder,'dir')~=7,mkdir(folder);end
f=figure('Visible',char(vis),'Color','w','Position',[50 50 1500 980]);
tl=tiledlayout(f,2,2,'Padding','compact','TileSpacing','compact');
for j=1:2
    ax=nexttile(tl);hold(ax,'on');
    if j==1,Z=X;C=T;support=oldSupport;else,Z=Y;C=S;support=newSupport;end
    if j==2
        patch(ax,'Faces',S(patchIDs,:),'Vertices',1e3*Y, ...
            'FaceColor',[.8 .91 .97],'EdgeColor','none');
        patch(ax,'Faces',S(collarIDs,:),'Vertices',1e3*Y, ...
            'FaceColor',[.93 .84 .96],'EdgeColor','none');
    end
    patch(ax,'Faces',C(support,:),'Vertices',1e3*Z, ...
        'FaceColor',[.97 .78 .35],'FaceAlpha',.65,'EdgeColor','none');
    patch(ax,'Faces',C,'Vertices',1e3*Z,'FaceColor','none', ...
        'EdgeColor',[.45 .48 .5],'LineWidth',.18);
    plot(ax,1e3*[inner(:,1);inner(1,1)],1e3*[inner(:,2);inner(1,2)], ...
        'b-','LineWidth',1.5);
    plot(ax,1e3*[outer(:,1);outer(1,1)],1e3*[outer(:,2);outer(1,2)], ...
        'm-','LineWidth',1.5);circles(ax,[ri ro]);
    axis(ax,'equal');xlim(ax,[-8.5 7.3]);ylim(ax,[-7.3 7.3]);
    xlabel(ax,'local x_1 (mm)');ylabel(ax,'local x_2 (mm)');
    if j==1,title(ax,'Original Step38; gold = primary q support'); ...
    else,title(ax,'Candidate; blue = paired, purple = collar, gold = q support');end
end
ax=nexttile(tl);hold(ax,'on');
patch(ax,'Faces',S,'Vertices',1e3*Y,'FaceColor','none', ...
    'EdgeColor',[.35 .4 .45],'LineWidth',.45);
plot(ax,1e3*Y(crNew.upperNodes,1),1e3*Y(crNew.upperNodes,2),'bo','MarkerSize',4);
plot(ax,1e3*Y(crNew.lowerNodes,1),1e3*Y(crNew.lowerNodes,2),'rx','MarkerSize',5);
axis(ax,'equal');xlim(ax,[-.35 .35]);ylim(ax,[-.35 .35]);
title(ax,'Paired tip: coincident distinct faces; shared ligament');
xlabel(ax,'local x_1 (mm)');ylabel(ax,'local x_2 (mm)');
ax=nexttile(tl);hold(ax,'on');
r=vecnorm((X(T(:,1),:)+X(T(:,2),:)+X(T(:,3),:))/3,2,2)*1e3;
s=vecnorm((Y(S(:,1),:)+Y(S(:,2),:)+Y(S(:,3),:))/3,2,2)*1e3;
scatter(ax,r(mask),oq.longest(mask)*1e3,3,[.65 .65 .65],'.');
scatter(ax,s(patchIDs),nq.longest(patchIDs)*1e3,4,[0 .4 .7],'.');
scatter(ax,s(collarIDs),nq.longest(collarIDs)*1e3,4,[.65 .15 .65],'.');
xlim(ax,[0 7]);ylim(ax,[0 .5]);grid(ax,'on');
xlabel(ax,'centroid radius (mm)');ylabel(ax,'longest T3 edge (mm)');
title(ax,'Measured grading: gray original, blue paired, purple collar');
title(tl,sprintf('Step61 asymmetric mesh-only candidate; paired radius %.3f mm; no physical U',rp*1e3));
files=struct('overviewPNG',[prefix '_overview.png'],'tipPNG',[prefix '_tip.png']);
exportgraphics(f,files.overviewPNG,'Resolution',200);close(f);
f=figure('Visible',char(vis),'Color','w','Position',[50 50 1100 650]);
tl=tiledlayout(f,1,2,'Padding','compact');
for j=1:2
    ax=nexttile(tl);hold(ax,'on');
    if j==1,Z=X;C=T;cc=cr;else,Z=Y;C=S;cc=crNew;end
    patch(ax,'Faces',C,'Vertices',1e3*Z,'FaceColor','none', ...
        'EdgeColor',[.3 .4 .5],'LineWidth',.55);
    plot(ax,1e3*Z(cc.upperNodes,1),1e3*Z(cc.upperNodes,2),'bo','MarkerSize',4);
    plot(ax,1e3*Z(cc.lowerNodes,1),1e3*Z(cc.lowerNodes,2),'rx','MarkerSize',5);
    axis(ax,'equal');xlim(ax,[-.22 .22]);ylim(ax,[-.22 .22]);grid(ax,'on');
    xlabel(ax,'local x_1 (mm)');ylabel(ax,'local x_2 (mm)');
    if j==1,title(ax,'Step38 actual topology');else,title(ax,'Step61 paired topology');end
end
exportgraphics(f,files.tipPNG,'Resolution',240);close(f);
end
function circles(ax,rr)
t=linspace(0,2*pi,300);
for r=rr,plot(ax,r*1e3*cos(t),r*1e3*sin(t),'k--','LineWidth',.9);end
end
function save_report(prefix,O61)
folder=fileparts(prefix);if ~isempty(folder)&&exist(folder,'dir')~=7,mkdir(folder);end
save([prefix '_small_data.mat'],'O61','-v7');
names={'summary','extentTable','qualityTable','radialTable','samplingTable','tipEdgeTable'};
for k=1:numel(names)
    if isfield(O61,names{k}),writetable(O61.(names{k}), ...
        [prefix '_' names{k} '.csv']);end
end
if isfield(O61,'synthetic')&&isfield(O61.synthetic,'table')
    writetable(O61.synthetic.table,[prefix '_syntheticTable.csv']);
end
end
function assert_audit_branch(root)
[status,branch]=system(sprintf('git -C "%s" branch --show-current',root));
assert(status==0&&strcmp(strtrim(branch),'sif-asymmetric-mesh-audit'), ...
    'step61:Branch','Step61 must run on sif-asymmetric-mesh-audit.');
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
