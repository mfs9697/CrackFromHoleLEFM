function F = qualify_incremental_crack_candidate(pathGlobal,varargin)
%QUALIFY_INCREMENTAL_CRACK_CANDIDATE
% Build and qualify an arbitrary N-segment incremental crack candidate.
%
% pathGlobal is (N+1)-by-2, ordered mouth -> current tip, with N>=2.
% Production use assumes the frozen 4-mm increment on every segment.
%
% This function performs NO physical FEM solve. It generalizes the accepted
% Stage III-C true-polyline qualification while preserving the same audited
% new-tip core, exterior grading, topology, native COD, and synthetic EDI.
%
% The current-tip asymptotic frame is always the LAST crack segment.

    validateattributes(pathGlobal,{'numeric'},{'2d','ncols',2,'finite'});
    assert(size(pathGlobal,1)>=3,'pathqual:NeedTwoSegments', ...
        'At least two crack segments are required for incremental qualification.');

    ip=inputParser;
    addParameter(ip,'FrozenState',[],@(x)isempty(x)||isstruct(x));
    addParameter(ip,'StateFile','',@(x)ischar(x)||isstring(x));
    addParameter(ip,'NArc',480,@(x)isnumeric(x)&&isscalar(x)&&isfinite(x)&&x>=32&&x==round(x));
    addParameter(ip,'ExteriorVerbose',true,@(x)islogical(x)&&isscalar(x));
    addParameter(ip,'RunSynthetic',true,@(x)islogical(x)&&isscalar(x));
    addParameter(ip,'SaveCandidate',true,@(x)islogical(x)&&isscalar(x));
    addParameter(ip,'CandidateFile','',@(x)ischar(x)||isstring(x));
    addParameter(ip,'SaveCompact',true,@(x)islogical(x)&&isscalar(x));
    addParameter(ip,'CompactFile','',@(x)ischar(x)||isstring(x));
    addParameter(ip,'Plot',false,@(x)islogical(x)&&isscalar(x));
    parse(ip,varargin{:});
    opt=ip.Results;

    incremental_profile_clock('begin','qualification',size(pathGlobal,1)-1);
    profileCleanup=onCleanup(@()incremental_profile_clock('end','qualification'));
    root=fileparts(mfilename('fullpath'));
    [R0,sourceLabel]=local_load_frozen_state(root,opt.FrozenState,char(opt.StateFile));

    assert(isfield(R0,'summary')&&istable(R0.summary)&&height(R0.summary)==1, ...
        'pathqual:FrozenSummary','Frozen R0.summary is required.');
    assert(isfield(R0,'C')&&isstruct(R0.C), ...
        'pathqual:FrozenConfig','Frozen R0.C is required.');

    S0=R0.summary;
    req={'a0_reserved_m','x_star_m','y_star_m', ...
        'nmat_x','nmat_y','that_x','that_y','stage1_pass','hole_npoly'};
    local_require_table_variables(S0,req);
    row=S0(1,:);
    assert(logical(row.stage1_pass),'pathqual:Stage1NotPassed', ...
        'Frozen Stage-I state did not pass.');
    assert(opt.NArc==row.hole_npoly,'pathqual:NArcMismatch', ...
        'Full-domain carrier must use frozen Stage-I hole_npoly=%d.',row.hole_npoly);

    C=R0.C;
    increment=row.a0_reserved_m;
    mouth=[row.x_star_m,row.y_star_m];
    nMat=[row.nmat_x,row.nmat_y]; nMat=nMat/norm(nMat);
    tHat=[row.that_x,row.that_y]; tHat=tHat/norm(tHat);
    if abs(dot(nMat,tHat))>1e-10 || det([nMat(:),tHat(:)])<=0
        error('pathqual:FrozenFrame','Frozen Stage-I frame is inconsistent.');
    end

    p0=pathGlobal(1,:);
    assert(norm(p0-mouth)<=2e-12,'pathqual:MouthChanged', ...
        'Supplied path does not start at the frozen initiation point.');

    seg=diff(pathGlobal,1,1);
    segLength=vecnorm(seg,2,2);
    nSegments=numel(segLength);
    assert(all(segLength>1e-12),'pathqual:DegenerateSegment', ...
        'Crack path contains a degenerate segment.');
    assert(max(abs(segLength-increment))<=2e-12,'pathqual:IncrementLength', ...
        'Production path must use the frozen %.12g-m increment on every segment.',increment);

    directions=seg./segLength;
    assert(norm(directions(1,:)-nMat)<=1e-12,'pathqual:FirstDirection', ...
        'First crack segment must retain theta_1=0 along the frozen material normal.');

    thetaSegments=atan2(directions*tHat(:),directions*nMat(:));
    thetaSegmentsDeg=rad2deg(thetaSegments);
    deltaThetaDeg=[thetaSegmentsDeg(1);diff(thetaSegmentsDeg)];

    eFirst=directions(1,:);
    eLast=directions(end,:);
    thetaLast=thetaSegments(end);
    thetaLastDeg=thetaSegmentsDeg(end);
    thetaPrevDeg=thetaSegmentsDeg(end-1);

    priorTip=pathGlobal(end-1,:);
    tip=pathGlobal(end,:);
    pathLength=sum(segLength);

    e1=eLast;
    e2=[-e1(2),e1(1)];
    R=[e1(:),e2(:)];
    pathLocal=(pathGlobal-tip)*R;

    fprintf('\n============================================================\n');
    fprintf('INCREMENTAL PATH: QUALIFY %d-SEGMENT CANDIDATE\n',nSegments);
    fprintf('============================================================\n');
    fprintf('  Frozen source   : %s\n',sourceLabel);
    fprintf('  NO physical FEM solve is performed.\n');
    fprintf('  increment       : %.9f mm\n',1e3*increment);
    fprintf('  path length     : %.9f mm\n',1e3*pathLength);
    fprintf('  theta_current   : %+.12g deg\n',thetaLastDeg);
    fprintf('  Delta theta     : %+.12g deg from prior leg\n',thetaLastDeg-thetaPrevDeg);
    fprintf('  current tip     : [%.12g, %.12g] m\n',tip(1),tip(2));
    fprintf('  prior tip local : [%.12g, %.12g] mm\n', ...
        1e3*pathLocal(end-1,1),1e3*pathLocal(end-1,2));

    % ------------------------------------------------------------------
    % True polyline source carrier. No source triangle is retained later.
    % ------------------------------------------------------------------
    incremental_profile_clock('phase','qualification','source_carrier');
    D=build_domain_hole_true_polyline( ...
        pathGlobal,C.A,C.B,C.holes,C.mesh2.chw, ...
        'corner_tol',1e-10,'epsMode','arclength', ...
        'nArc',opt.NArc,'orientation','cw','miter_limit',6);

    M=mesh_hole_pencil_domain(D, ...
        'Hmin',C.mesh1.hmin,'Hmax',C.mesh1.hmax,'Hgrad',C.mesh1.hgrad, ...
        'PlotGeom',false,'PlotMesh',false);

    incremental_profile_clock('phase','qualification','geometry_id_recovery');
    polyIDs=identify_polyline_pencil_edge_sets(M,D,'Verbose',true);
    incremental_profile_clock('phase','qualification','crack_face_collapse');
    Mc=collapse_polyline_pencil_faces_to_midline(M,D, ...
        'UpperEdgeIDs',polyIDs.upperEdges, ...
        'LowerEdgeIDs',polyIDs.lowerEdges, ...
        'TipVertexID',polyIDs.v_tip);

    incremental_profile_clock('phase','qualification','carrier_checks');
    assert(norm(Mc.crack.x0-p0)<=1e-12*max(1,pathLength), ...
        'pathqual:CarrierMouth','Carrier mouth changed.');
    assert(norm(Mc.crack.xtip-tip)<=1e-12*max(1,pathLength), ...
        'pathqual:CarrierTip','Carrier tip changed.');
    assert(norm(Mc.crack.Pmid-pathGlobal,'fro')<=1e-12, ...
        'pathqual:CarrierPolyline','Collapsed carrier did not retain the supplied polyline.');

    X=(Mc.p-tip)*R;
    Told=Mc.t;

    crLocal=Mc.crack;
    crLocal.Pmid=pathLocal;
    crLocal.x0=pathLocal(1,:);
    crLocal.xtip=[0,0];

    sourceArea=sum(local_signed_area(X,Told));
    assert(all(local_signed_area(X,Told)>0), ...
        'pathqual:CarrierOrientation','Carrier contains nonpositive T3 area.');

    interiorVertices=2:size(pathGlobal,1)-1;
    carrierVertexPairs=true;
    for kk=interiorVertices
        uu=Mc.crack.pathVertexUpperNodes{kk};
        ll=Mc.crack.pathVertexLowerNodes{kk};
        carrierVertexPairs=carrierVertexPairs && numel(uu)==1 && ...
            numel(ll)==1 && uu~=ll;
    end
    carrierMultiEdge=numel(polyIDs.upperEdges)>=nSegments && ...
        numel(polyIDs.lowerEdges)>=nSegments;
    chordDeviation=0;
    for kk=interiorVertices
        chordDeviation=max(chordDeviation, ...
            local_point_segment_distance(pathGlobal(kk,:),p0,tip));
    end

    % ------------------------------------------------------------------
    % Exact already-qualified structured core at scale s=1.
    % ------------------------------------------------------------------
    incremental_profile_clock('phase','qualification','structured_core');
    Core=build_stage2_scaled_audited_core([0 0],[1 0],increment,'Scale',1);
    Zp=Core.local.coord3;
    Tp=Core.local.connect3;
    rp=Core.rCore;
    design=Core.design;
    design.rCore_m=rp;
    design.transitionLength_m=1.0*increment;
    design.farCap_m=0.625*increment;
    design.exteriorCalibration=struct('verbose',opt.ExteriorVerbose);

    ri=.10*increment;
    ro=.65*increment;

    assert(abs(Core.hTip/increment-0.00675308135)<1e-14, ...
        'pathqual:HTipScale','Unexpected core hTip/a0.');
    assert(abs(rp/increment-.75)<1e-14,'pathqual:CoreScale', ...
        'Unexpected core radius/a0.');

    % Physical-boundary clearance in the actual carrier geometry.
    incremental_profile_clock('phase','qualification','carrier_clearance');
    [freeOld,~]=local_edge_inventory(Told);
    freeOld=freeOld(:,1:2);
    faceOld=local_face_edge_mask(freeOld,crLocal.upperNodes, ...
        crLocal.lowerNodes,crLocal.tipNode);
    physicalOld=freeOld(~faceOld,:);
    physicalClearance=local_min_segment_distance([0 0],X,physicalOld);
    assert(physicalClearance>rp, ...
        'pathqual:CoreHitsPhysicalBoundary', ...
        'Paired core reaches a noncrack physical boundary.');

    % ------------------------------------------------------------------
    % Audited C03-style full exterior, scaled by a0.
    % ------------------------------------------------------------------
    incremental_profile_clock('phase','qualification','exterior');
    [Ze,Te,ext,physicalIDs]=build_stage3c_polyline_exterior( ...
        X,Told,crLocal,Zp,Tp,rp,design);

    incremental_profile_clock('phase','qualification','geometry_mapping');
    outerIDs=ext.originalPhysicalIDs;
    upOld=unique(crLocal.upperNodes(:));
    loOld=unique(crLocal.lowerNodes(:));
    axisUp=Core.design.axisIDs;
    mirrorMap=Core.mirror.T3;
    nUpper=Core.mirror.nUpperT3Nodes;

    % Assemble in LOCAL coordinates. Preserve source physical vertices
    % exactly; weld the exterior to the core boundary by exact coordinates,
    % while never welding opposite crack-face IDs.
    Q=X;
    pairedIDs=zeros(size(Zp,1),1);
    for j=1:size(Zp,1)
        Q(end+1,:)=Zp(j,:); %#ok<AGROW>
        pairedIDs(j)=size(Q,1);
    end
    patchTriangles=pairedIDs(Tp);

    exteriorIDs=zeros(size(Ze,1),1);
    tol=1e-12;
    nRetainedCrackExteriorNodes=0;
    nRetainedCrackSourceMatches=0;
    upperAxis=axisUp(Zp(axisUp,1)<0);
    lowerAxis=mirrorMap(upperAxis);

    for j=1:size(Ze,1)
        z=Ze(j,:);
        side=ext.nodeSide(j);

        candidates=outerIDs;

        % For the retained crack OUTSIDE the paired core, use topology rather
        % than the old straight-crack test x_2=0. In a kinked current-tip
        % frame the earlier leg and crack mouth are generally off-axis.
        % ext.crackNodeMask identifies both copies of every retained crack
        % vertex after exterior duplication; nodeSide selects the source face.
        onRetainedCrack=logical(ext.crackNodeMask(j));
        if onRetainedCrack
            nRetainedCrackExteriorNodes=nRetainedCrackExteriorNodes+1;
            if side>=0
                candidates=setdiff(candidates,loOld);
            else
                candidates=setdiff(candidates,upOld);
            end
        end

        d=vecnorm(X(candidates,:)-z,2,2);
        [dd,ii]=min(d);
        if dd<tol
            exteriorIDs(j)=candidates(ii);
            if onRetainedCrack
                nRetainedCrackSourceMatches=nRetainedCrackSourceMatches+1;
            end
            continue
        end

        d=vecnorm(Zp-z,2,2);
        if z(1)<0&&abs(z(2))<tol
            if side>=0,d(lowerAxis)=Inf;else,d(upperAxis)=Inf;end
        end
        [dd,ii]=min(d);
        if dd<tol
            exteriorIDs(j)=pairedIDs(ii);
            continue
        end

        Q(end+1,:)=z; %#ok<AGROW>
        exteriorIDs(j)=size(Q,1);
    end

    Tc=[patchTriangles;exteriorIDs(Te)];
    nPatch=size(Tp,1);
    patchRows=(1:nPatch)';
    exteriorRows=(nPatch+1:size(Tc,1))';

    % Compact unused carrier/intermediate nodes.
    used=unique(Tc(:));
    oldToNew=zeros(size(Q,1),1);
    oldToNew(used)=(1:numel(used))';
    Pc=Q(used,:);
    Tc=oldToNew(Tc);
    pairedIDs=oldToNew(pairedIDs);
    exteriorIDs=oldToNew(exteriorIDs);

    % Final crack metadata from topology.
    cr=struct();
    cr.Pmid=pathLocal;
    cr.x0=pathLocal(1,:);
    cr.xtip=[0,0];
    cr.tipNode=pairedIDs(Core.crack.tipNode);

    newUpper=pairedIDs(axisUp(Zp(axisUp,1)<-tol));
    newLower=pairedIDs(mirrorMap(axisUp(Zp(axisUp,1)<-tol)));

    extCrack=ext.crackNodeMask;
    newUpper=[newUpper;exteriorIDs(extCrack & ext.nodeSide>=0)];
    newLower=[newLower;exteriorIDs(extCrack & ext.nodeSide<0)];

    cr.upperNodes=unique([newUpper;cr.tipNode]);
    cr.lowerNodes=unique([newLower;cr.tipNode]);
    cr=local_refresh_crack_metadata_polyline(cr,Pc);

    incremental_profile_clock('phase','qualification','t3_t6_local');
    [P6,T6]=T3toT6_fast(Pc,Tc);
    meshLocal=struct('coord3',Pc,'connect3',Tc,'coord',P6,'connect',T6);

    % Rotate final candidate back to global coordinates. Then restore every
    % original carrier physical-boundary vertex BITWISE, following the
    % closed-audit preservation rule. Rebuild T6 midsides from this final
    % T3 geometry; connectivity ordering must remain identical.
    incremental_profile_clock('phase','qualification','global_boundary_restore');
    Pg=tip+Pc*R.';
    assert(all(oldToNew(physicalIDs)>0), ...
        'pathqual:LostPhysicalVertex','A source physical vertex was lost.');
    Pg(oldToNew(physicalIDs),:)=Mc.p(physicalIDs,:);

    incremental_profile_clock('phase','qualification','t3_t6_global');
    [P6g,T6g]=T3toT6_fast(Pg,Tc);
    assert(isequal(T6g,T6), ...
        'pathqual:T6OrderingChanged','Global T6 connectivity ordering changed.');

    crGlobal=cr;
    crGlobal.Pmid=pathGlobal;
    crGlobal.x0=p0;
    crGlobal.xtip=tip;
    crGlobal.upperTarget=Pg(cr.upperNodes,:);
    crGlobal.lowerTarget=Pg(cr.lowerNodes,:);

    meshGlobal=struct('coord3',Pg,'connect3',Tc,'coord',P6g,'connect',T6g);

    % ------------------------------------------------------------------
    % Structural qualification.
    % ------------------------------------------------------------------
    incremental_profile_clock('phase','qualification','structural_qualification');
    area=local_signed_area(Pc,Tc);
    q=local_quality(Pc,Tc);
    [minJ,midErr]=local_jacobian_audit(P6,T6);

    finalArea=sum(area);
    areaRel=abs(finalArea-sourceArea)/max(abs(sourceArea),eps);

    [freeNew,allNew]=local_edge_inventory(Tc);
    maxIncidence=max(allNew(:,3));
    duplicateTri=local_has_duplicate_triangles(Tc);

    % Every geometric outer edge of the structured core must be shared by
    % exactly one exterior triangle after assembly. Exterior centroids must
    % remain outside the paired-core polygon.
    [seamShared,exteriorOutsideCore]=local_seam_audit( ...
        Pc,Tc,patchRows,exteriorRows,ext.innerPolygon);

    % Exact preservation of the structured core after assembly.
    coreCoordErr=max(vecnorm(Pc(pairedIDs,:)-Zp,2,2));
    coreConnExact=isequal(Tc(patchRows,:),pairedIDs(Tp));

    [pair3,pair6]=local_pair_errors(meshLocal,pairedIDs,mirrorMap,nUpper,patchRows);
    [topology,nativeR,faceSide]=local_topology_audit(meshLocal,cr, ...
        pairedIDs,axisUp,mirrorMap,Zp,C);

    % Original physical boundary geometry, permitting exact straight-segment
    % subdivision but no displacement of the carrier boundary.
    physicalBoundarySame=local_boundary_same( ...
        X,Told,Pc,Tc,crLocal,cr);
    physicalVerticesSame=all(oldToNew(physicalIDs)>0) && ...
        isequal(Pc(oldToNew(physicalIDs),:),X(physicalIDs,:));
    globalPhysicalVerticesBitwise= ...
        isequal(Pg(oldToNew(physicalIDs),:),Mc.p(physicalIDs,:));

    % EDI support must be wholly within the untouched core rows.
    support=local_q_support(meshLocal,cr,ri,ro,false);
    supportSkipConstant=local_q_support(meshLocal,cr,ri,ro,true);
    supportInside=all(ismember(support,patchRows)) && ...
        all(ismember(supportSkipConstant,patchRows));
    exteriorOutside=~any(ismember(support,exteriorRows)) && ...
        ~any(ismember(supportSkipConstant,exteriorRows));

    % The retained polyline outside the new-tip core spans the crack mouth
    % through every historical kink to the rear core intersection.
    exteriorCrackLength=pathLength-rp;
    sU=cr.upperS;
    sL=cr.lowerS;
    mouthUpper=min(sU);
    mouthLower=min(sL);

    % Native COD windows are entirely inside the paired core.
    windows=[.04 .20;.04 .30;.08 .30;.12 .30];
    sampleN=zeros(4,1);
    for k=1:4
        sampleN(k)=nnz(nativeR/increment>=windows(k,1) & nativeR/increment<=windows(k,2));
    end

    affected=[patchRows;exteriorRows];
    maxNeighbor=local_max_neighbor_size_ratio(Tc,q.longest,affected);

    % Previous Stage-II tip must survive as an exact duplicated crack-face
    % vertex one increment behind the new tip.
    cornerLocal=pathLocal(end-1,:);
    upCorner=cr.upperNodes(vecnorm(Pc(cr.upperNodes,:)-cornerLocal,2,2)<=1e-12);
    loCorner=cr.lowerNodes(vecnorm(Pc(cr.lowerNodes,:)-cornerLocal,2,2)<=1e-12);
    cornerDistance=norm(cornerLocal);
    mouthDistance=norm(pathLocal(1,:));

    gates=struct();
    gates.frozenStage1Pass=logical(row.stage1_pass);
    gates.carrierUsesFrozenHoleResolution=opt.NArc==row.hole_npoly;
    gates.fullPathRetained=norm(crGlobal.Pmid-pathGlobal,'fro')<=1e-12;
    gates.mouthUnchanged=norm(crGlobal.Pmid(1,:)-p0)<=1e-12;
    gates.priorTipUnchanged=norm(crGlobal.Pmid(end-1,:)-priorTip)<=1e-12;
    gates.newTipUnchanged=norm(crGlobal.Pmid(end,:)-tip)<=1e-12;
    gates.uniformIncrement=max(abs(segLength-increment))<=2e-12;
    gates.firstLegDirection=norm(directions(1,:)-nMat)<=1e-12;
    gates.lastLegDirection=norm(directions(end,:)-eLast)<=1e-13;
    gates.segmentAnglesFinite=all(isfinite(thetaSegmentsDeg));
    gates.carrierUsesAllPolylineEdges=carrierMultiEdge;
    gates.carrierInteriorVerticesDistinct=carrierVertexPairs;
    gates.retainedExteriorPathExact=norm(ext.retainedCrackPath- ...
        [pathLocal(1:end-1,:);[-rp,0]],'fro')<=1e-12;
    gates.sideAwareRetainedCrackMapping=nRetainedCrackExteriorNodes>0 && ...
        nRetainedCrackSourceMatches>=2;
    gates.pathNonStraight=chordDeviation>1e-12;
    gates.priorTipUpperNode=numel(upCorner)==1;
    gates.priorTipLowerNode=numel(loCorner)==1;
    gates.priorTipFacesDistinct=numel(upCorner)==1&&numel(loCorner)==1&&upCorner~=loCorner;
    gates.priorTipDistance=abs(cornerDistance-increment)<=1e-12;
    gates.priorTipOutsideCore=cornerDistance>rp+1e-12;
    gates.priorTipOutsideEDI=cornerDistance>ro+1e-12;
    gates.mouthOutsideEDI=mouthDistance>ro+1e-12;
    gates.physicalBoundaryOutsideCore=physicalClearance>rp+1e-12;
    gates.physicalBoundaryOutsideEDI=physicalClearance>ro+1e-12;
    gates.coreBitwiseCoordinates=coreCoordErr==0;
    gates.coreConnectivityExact=coreConnExact;
    gates.completeT3CorePairing=pair3<=1e-13;
    gates.completeT6CorePairing=pair6<=1e-13;
    gates.crackFacesDistinct=topology.crackDistinct;
    gates.intactLigamentShared=topology.intactShared;
    gates.crackT6MidsidesDistinct=topology.crackMidsDistinct;
    gates.intactT6MidsidesShared=topology.intactMidsShared;
    gates.T6FaceAbscissaeMatched=topology.faceGridMismatch<=1e-13;
    gates.lastSegmentCODFrame=topology.usesLastSegmentCODFrame;
    gates.positiveT3Areas=all(area>0);
    gates.positiveT6Jacobians=minJ>0 && midErr<=1e-13;
    gates.edgeIncidenceValid=maxIncidence<=2;
    gates.noDuplicateTriangles=~duplicateTri;
    gates.structuredCoreSeamShared=seamShared;
    gates.exteriorOutsideCore=exteriorOutsideCore;
    gates.domainAreaPreserved=areaRel<=5e-10;
    gates.physicalBoundaryGeometryPreserved=physicalBoundarySame;
    gates.originalPhysicalVerticesPreserved=physicalVerticesSame;
    gates.globalPhysicalVerticesBitwise=globalPhysicalVerticesBitwise;
    gates.EDIInsideUntouchedCore=supportInside;
    gates.exteriorExcludedFromEDI=exteriorOutside;
    gates.nativeSamplingAdequate=all(sampleN>=12);
    gates.nativeSamplingExact=isequal(sampleN,[38;55;44;34]);
    gates.exteriorCrackLengthCorrect=abs(exteriorCrackLength-(pathLength-rp))<=1e-14;
    gates.upperMouthAtZero=mouthUpper<=1e-12;
    gates.lowerMouthAtZero=mouthLower<=1e-12;
    gates.minimumAngle20deg=min(q.minAngle(affected))>=20;
    gates.neighborRatioAtMost18=maxNeighbor<=1.8+1e-10;

    structuralPass=all(structfun(@logical,gates));

    fprintf('\nFULL-DOMAIN MESH SUMMARY\n');
    fprintf('  carrier T3           = %d nodes / %d elements\n',size(X,1),size(Told,1));
    fprintf('  final T3             = %d nodes / %d elements\n',size(Pc,1),size(Tc,1));
    fprintf('  final T6             = %d nodes / %d elements\n',size(P6,1),size(T6,1));
    fprintf('  paired core T3       = %d elements\n',numel(patchRows));
    fprintf('  exterior T3          = %d elements\n',numel(exteriorRows));
    fprintf('  physical clearance   = %.6f mm\n',1e3*physicalClearance);
    fprintf('  exterior crack       = %.6f mm (mouth to rear core)\n',1e3*exteriorCrackLength);
    fprintf('  prior-tip distance   = %.6f mm from new tip\n',1e3*cornerDistance);
    fprintf('  kink/chord deviation = %.6e mm\n',1e3*chordDeviation);
    fprintf('  carrier face edges   = %d upper / %d lower\n',numel(polyIDs.upperEdges),numel(polyIDs.lowerEdges));
    fprintf('  retained crack nodes = %d exterior copies; %d source-side matches\n', ...
        nRetainedCrackExteriorNodes,nRetainedCrackSourceMatches);
    fprintf('  primary EDI elements = %d (literal), %d (skip-constant)\n', ...
        numel(support),numel(supportSkipConstant));
    fprintf('  core coord error     = %.3e m\n',coreCoordErr);
    fprintf('  core pair error T3/6 = %.3e / %.3e m\n',pair3,pair6);
    fprintf('  area relative error  = %.3e\n',areaRel);
    fprintf('  min angle            = %.6f deg\n',min(q.minAngle));
    fprintf('  max neighbor ratio   = %.9f\n',maxNeighbor);
    fprintf('  min T6 detJ          = %.6e m^2\n',minJ);

    fprintf('\nFULL-DOMAIN STRUCTURAL GATES\n');
    local_print_gates(gates);

    if ~structuralPass
        F=local_partial_output();
        error('pathqual:StructuralFailed', ...
            'Full-domain embedding failed at least one structural gate; no synthetic EDI was run.');
    end

    % ------------------------------------------------------------------
    incremental_profile_clock('phase','qualification','synthetic_qualification');
    % Full-mesh prescribed Williams qualification.
    % ------------------------------------------------------------------
    mat=local_material(C);
    if opt.RunSynthetic
    cases=[1 0;0 1;1 1e-4];
    names=["pure_I";"pure_II";"tiny_mixed"];
    rows=nan(3,8);

    up=find(faceSide==1);
    lo=find(faceSide==-1);
    Z=P6;
    Z([up;lo],2)=0;

    for k=1:3
        incremental_profile_clock('begin','williams_replay',nSegments);
        Uloc=exact_williams_displacement_audit( ...
            Z,cases(k,1),cases(k,2),mat.E,mat.nu,mat.ps, ...
            'UpperFaceIDs',up,'LowerFaceIDs',lo);
        ug=reshape(Uloc,2,[]).'*R.';
        U=zeros(2*size(P6,1),1);
        U(1:2:end)=ug(:,1);
        U(2:2:end)=ug(:,2);

        incremental_profile_clock('end','williams_replay');
        incremental_profile_clock('begin','synthetic_edi',nSegments);
        [ki,kii,aux]=SIF_LEFM_interaction_EDI( ...
            meshGlobal,U,crGlobal.Pmid,mat, ...
            struct('r_inner',ri,'r_outer',ro), ...
            'UsePlaneStrain',mat.ps==1, ...
            'WeightFunction','fe_nodal', ...
            'QuadratureRule',16, ...
            'StoreGPDiagnostics',false, ...
            'Verbose',false);

        incremental_profile_clock('end','synthetic_edi');
        rows(k,:)=[cases(k,:),ki,kii, ...
            ki-cases(k,1),kii-cases(k,2),aux.nElem_used,aux.nGP_used];
    end

    Synthetic=array2table(rows,'VariableNames',{ ...
        'KI_input','KII_input','KI_recovered','KII_recovered', ...
        'KI_error','KII_error','nElem_used','nGP_used'});
    Synthetic.caseName=names;
    Synthetic=movevars(Synthetic,'caseName','Before','KI_input');

    M=rows(1:2,3:4).';
    matrixError=norm(M-eye(2),'fro');
    mixedRel=abs(rows(3,4)/1e-4-1);
    superErr=norm(rows(3,3:4).'-M*[1;1e-4]);

    syntheticGates=struct();
    syntheticGates.sameEDISupportCountAsCore=all(Synthetic.nElem_used==11316);
    syntheticGates.pureIRecovery=abs(rows(1,3)-1)<=2e-4;
    syntheticGates.pureICrossLeakage=abs(rows(1,4))<=1e-10;
    syntheticGates.pureIIRecovery=abs(rows(2,4)-1)<=2e-4;
    syntheticGates.pureIICrossLeakage=abs(rows(2,3))<=1e-10;
    syntheticGates.recoveryMatrix=matrixError<=2e-4;
    syntheticGates.tinyMixedKII=mixedRel<=2e-4;
    syntheticGates.superposition=superErr<=1e-10;
    syntheticPass=all(structfun(@logical,syntheticGates));

    fprintf('\nFULL-DOMAIN PRESCRIBED WILLIAMS QUALIFICATION\n');
    disp(Synthetic);
    fprintf('  ||M-I||_F                = %.6e\n',matrixError);
    fprintf('  tiny-mixed KII rel error = %.6e\n',mixedRel);
    fprintf('  superposition residual   = %.6e\n',superErr);
    fprintf('\nSYNTHETIC GATES\n');
    local_print_gates(syntheticGates);

    else
        Synthetic=table();
        matrixError=NaN; mixedRel=NaN; superErr=NaN;
        syntheticGates=struct('skippedByExplicitOption',true);
        syntheticPass=true;
        fprintf('\nSYNTHETIC REPLAY SKIPPED BY EXPLICIT OPTION.\n');
    end

    incremental_profile_clock('phase','qualification','result_packaging');
    pass=structuralPass&&syntheticPass;

    Summary=table( ...
        nSegments,increment,pathLength,thetaPrevDeg,thetaLastDeg,Core.hTip,ri,ro,rp, ...
        design.transitionLength_m,design.farCap_m, ...
        size(Pc,1),size(Tc,1),size(P6,1),numel(patchRows),numel(exteriorRows), ...
        physicalClearance,cornerDistance,exteriorCrackLength,numel(support), ...
        coreCoordErr,pair3,pair6,areaRel,min(q.minAngle),maxNeighbor,minJ, ...
        matrixError,mixedRel,superErr,pass, ...
        'VariableNames',{ ...
        'n_segments','increment_m','path_length_m','theta_prev_deg','theta_current_deg', ...
        'hTip_m','rInner_m','rOuter_m','rCore_m','transition_m','farCap_m', ...
        'T3_nodes','T3_elements','T6_nodes','core_T3_elements','exterior_T3_elements', ...
        'physical_clearance_m','prior_tip_distance_m','exterior_crack_m','EDI_elements', ...
        'core_coord_error_m','pair_error_T3_m','pair_error_T6_m','area_rel_error', ...
        'min_angle_deg','max_neighbor_ratio','min_T6_detJ_m2', ...
        'recovery_matrix_error','tiny_mixed_KII_rel_error', ...
        'superposition_residual','pass'});

    candidate=struct();
    candidate.p=Pg;
    candidate.t=Tc;
    candidate.crack=crGlobal;
    candidate.mat=mat;
    candidate.nSegments=nSegments;
    candidate.segmentLengths=segLength;
    candidate.currentIncrementLength=increment;
    candidate.totalCrackLength=pathLength;
    candidate.thetaSegments=thetaSegments;
    candidate.thetaSegmentsDeg=thetaSegmentsDeg;
    candidate.deltaThetaSegmentsDeg=deltaThetaDeg;
    candidate.thetaCurrent=thetaLast;
    candidate.thetaCurrentDeg=thetaLastDeg;
    candidate.thetaPreviousDeg=thetaPrevDeg;
    candidate.deltaThetaCurrentDeg=thetaLastDeg-thetaPrevDeg;
    candidate.directionFirst=eFirst;
    candidate.directionLast=eLast;
    candidate.path=pathGlobal;
    candidate.priorTip=priorTip;
    candidate.nMatFrozen=nMat;
    candidate.tHatFrozen=tHat;
    candidate.pairedElementIDs=patchRows;
    candidate.exteriorElementIDs=exteriorRows;
    candidate.primarySupportElementIDs=support;
    candidate.skipConstantSupportElementIDs=supportSkipConstant;
    candidate.pairedNodeIDs=pairedIDs;
    candidate.pairedMirrorLocal=mirrorMap;
    candidate.structuredDesign=design;
    candidate.exteriorDesign=ext;
    candidate.polylineCarrierEdgeIDs=polyIDs;
    candidate.maxChordDeviation=chordDeviation;
    candidate.gates=gates;
    candidate.synthetic=Synthetic;
    candidate.syntheticGates=syntheticGates;
    candidate.syntheticReplayPerformed=logical(opt.RunSynthetic);
    candidate.scientificallyReadyForIncrementalPhysicalSolve=pass;
    candidate.source=sourceLabel;

    F=struct();
    F.summary=Summary;
    F.gates=gates;
    F.synthetic=Synthetic;
    F.syntheticGates=syntheticGates;
    F.candidate=candidate;
    F.mesh=meshGlobal;
    F.localMesh=meshLocal;
    F.crackLocal=cr;
    F.Core=Core;
    F.exteriorDesign=ext;
    F.sampleCounts=table(windows(:,1),windows(:,2),sampleN, ...
        'VariableNames',{'lower_r_over_DeltaA','upper_r_over_DeltaA','nativePoints'});
    F.pass=pass;

    if opt.Plot
        local_plot_full_mesh(Pg,Tc,crGlobal,patchRows,exteriorRows,ri,ro,tip,R);
    end

    outDir=fullfile(root,'verification','crack_path');
    if exist(outDir,'dir')~=7,mkdir(outDir);end

    incremental_profile_clock('phase','qualification','candidate_save');
    if opt.SaveCandidate
        file=char(opt.CandidateFile);
        if isempty(file)
            file=fullfile(outDir, ...
                sprintf('incremental_step_%03d_candidate_T3.mat',nSegments));
        end
        save(file,'candidate','-v7');
        F.candidateFile=file;
        fprintf('  Candidate MAT: %s\n',file);
    else
        F.candidateFile='';
    end

    incremental_profile_clock('phase','qualification','compact_save');
    if opt.SaveCompact
        file=char(opt.CompactFile);
        if isempty(file)
            file=fullfile(outDir, ...
                sprintf('incremental_step_%03d_qualification_small.mat',nSegments));
        end
        Small=rmfield(F,{'candidate','mesh','localMesh','Core'});
        save(file,'Small','-v7');
        F.compactFile=file;
        fprintf('  Compact MAT: %s\n',file);
    else
        F.compactFile='';
    end

    incremental_profile_clock('phase','qualification','finish');
    if pass
        fprintf('\nINCREMENTAL CANDIDATE QUALIFICATION PASS.\n');
        fprintf('  Full assembled T3/T6 mesh and prescribed EDI are qualified.\n');
        fprintf('  Candidate is ready for one guarded incremental physical solve.\n');
        fprintf('  No physical solve has been performed here.\n');
    end

    function O=local_partial_output()
        O=struct('gates',gates,'pass',false,'Core',Core, ...
            'localMesh',meshLocal,'mesh',meshGlobal);
    end
end


% =========================================================================
function [R0,label]=local_load_frozen_state(root,Rin,stateFile)
    if ~isempty(Rin)
        R0=Rin;label='<in-memory FrozenState>';return
    end
    if isempty(strtrim(stateFile))
        stateFile=fullfile(root,'verification','crack_path','stage1_starting_state.mat');
    elseif exist(stateFile,'file')~=2 && ~local_is_absolute_path(stateFile)
        q=fullfile(root,stateFile);if exist(q,'file')==2,stateFile=q;end
    end
    if exist(stateFile,'file')~=2
        error('pathqual:MissingFrozenState', ...
            'Frozen Stage-I MAT not found; pass ''FrozenState'',R0.');
    end
    d=load(stateFile);assert(isfield(d,'R0')&&isstruct(d.R0));
    R0=d.R0;label=stateFile;
end

function local_require_table_variables(T,names)
    miss=names(~ismember(names,T.Properties.VariableNames));
    if ~isempty(miss),error('pathqual:FrozenFields', ...
            'Missing frozen summary fields: %s',strjoin(miss,', '));end
end

function cr=local_refresh_crack_metadata_polyline(cr,P)
    cr.upperS=local_polyline_parameter(P(cr.upperNodes,:),cr.Pmid);
    cr.lowerS=local_polyline_parameter(P(cr.lowerNodes,:),cr.Pmid);
    [cr.upperS,ix]=sort(cr.upperS);cr.upperNodes=cr.upperNodes(ix);
    [cr.lowerS,ix]=sort(cr.lowerS);cr.lowerNodes=cr.lowerNodes(ix);
    cr.nUpper=numel(cr.upperNodes);cr.nLower=numel(cr.lowerNodes);
    cr.sameCount=cr.nUpper==cr.nLower;
    cr.upperTarget=P(cr.upperNodes,:);cr.lowerTarget=P(cr.lowerNodes,:);
    cr.lowerMatchForUpper=zeros(cr.nUpper,1);
    for k=1:cr.nUpper
        [~,cr.lowerMatchForUpper(k)]=min(abs(cr.lowerS-cr.upperS(k)));
    end
end

function s=local_polyline_parameter(X,P)
    seg=diff(P,1,1);L=vecnorm(seg,2,2);cum=[0;cumsum(L)];Lt=cum(end);
    s=zeros(size(X,1),1);
    for i=1:size(X,1)
        best=inf;bs=0;
        for k=1:numel(L)
            v=seg(k,:);tt=dot(X(i,:)-P(k,:),v)/dot(v,v);
            tt=max(0,min(1,tt));q=P(k,:)+tt*v;d=norm(X(i,:)-q);
            if d<best,best=d;bs=(cum(k)+tt*L(k))/Lt;end
        end
        s(i)=bs;
    end
end

function d=local_point_segment_distance(x,A,B)
    v=B-A;L2=max(dot(v,v),1e-30);
    t=dot(x-A,v)/L2;t=max(0,min(1,t));
    d=norm(x-(A+t*v));
end

function a=local_signed_area(P,T)
    b=P(T(:,2),:)-P(T(:,1),:);c=P(T(:,3),:)-P(T(:,1),:);
    a=.5*(b(:,1).*c(:,2)-b(:,2).*c(:,1));
end

function q=local_quality(P,T)
    a=P(T(:,1),:);b=P(T(:,2),:);c=P(T(:,3),:);
    L=[vecnorm(b-c,2,2),vecnorm(a-c,2,2),vecnorm(a-b,2,2)];
    ang=zeros(size(L));
    for k=1:3
        j=mod(k,3)+1;h=mod(k+1,3)+1;
        ang(:,k)=acosd(max(-1,min(1,(L(:,j).^2+L(:,h).^2-L(:,k).^2)./ ...
            (2*L(:,j).*L(:,h)))));
    end
    area=local_signed_area(P,T);
    q=struct('area',area,'longest',max(L,[],2),'minAngle',min(ang,[],2), ...
        'shape',4*sqrt(3)*area./sum(L.^2,2));
end

function [free,all]=local_edge_inventory(T)
    E=sort([T(:,[1 2]);T(:,[2 3]);T(:,[3 1])],2);
    [U,~,g]=unique(E,'rows');count=accumarray(g,1);
    all=[U,count];free=all(count==1,:);
end

function mask=local_face_edge_mask(E,up,lo,tip)
    mask=all(ismember(E,[up(:);tip]),2)|all(ismember(E,[lo(:);tip]),2);
end

function d=local_min_segment_distance(z,P,E)
    a=P(E(:,1),:);b=P(E(:,2),:);v=b-a;
    t=sum((z-a).*v,2)./sum(v.^2,2);t=max(0,min(1,t));
    d=min(vecnorm(a+t.*v-z,2,2));
end

function [err3,err6]=local_pair_errors(mesh,ids,map,nU,rows)
    n=size(rows,1)/2;
    U=mesh.connect3(rows(1:n),:);
    L=mesh.connect3(rows(n+1:end),[1 3 2]);
    globalMap=zeros(size(mesh.coord3,1),1);
    globalMap(ids(1:nU))=ids(map);
    assert(isequal(globalMap(U),L), ...
        'pathqual:PairT3Connectivity','Core T3 connectivity is not reflected.');

    x=mesh.coord3;
    err3=max(vecnorm(x(L(:),:)-[x(U(:),1),-x(U(:),2)],2,2));

    U6=mesh.connect(rows(1:n),:);
    L6=mesh.connect(rows(n+1:end),[1 3 2 6 5 4]);
    x=mesh.coord;
    err6=max(vecnorm(x(L6(:),:)-[x(U6(:),1),-x(U6(:),2)],2,2));
end

function [out,r,face]=local_topology_audit(mesh,cr,ids,axisIDs,map,Z,C)
    up=cr.upperNodes;lo=cr.lowerNodes;shared=intersect(up,lo);
    out=struct('crackDistinct',isequal(shared,cr.tipNode), ...
        'intactShared',true,'crackMidsDistinct',true,'intactMidsShared',true);

    T=mesh.connect;maps=[1 2 4;2 3 5;3 1 6];
    intact=axisIDs(Z(axisIDs,1)>=0);[~,ix]=sort(Z(intact,1));intact=intact(ix);
    negative=axisIDs(Z(axisIDs,1)<=0);[~,ix]=sort(Z(negative,1));negative=negative(ix);

    for k=1:numel(intact)-1
        E=ids(intact(k:k+1));m=local_find_mid(T,maps,E);
        out.intactShared=out.intactShared&&numel(m)==1&& ...
            ids(intact(k))==ids(map(intact(k)));
        out.intactMidsShared=out.intactMidsShared&&numel(m)==1;
    end

    for k=1:numel(negative)-1
        u=ids(negative(k:k+1));l=ids(map(negative(k:k+1)));
        mu=local_find_mid(T,maps,u);ml=local_find_mid(T,maps,l);
        out.crackMidsDistinct=out.crackMidsDistinct&& ...
            numel(mu)==1&&numel(ml)==1&&mu~=ml;
    end

    mat=local_material(C);
    [r,~,diag]=native_COD_polyline_audit(mesh,zeros(2*size(mesh.coord,1),1), ...
        mat,cr,2,true);
    face=diag.faceSide;
    out.faceGridMismatch=diag.gridMismatch;
    out.usesLastSegmentCODFrame=logical(diag.usesLastSegmentFrame);
    out.crackDistinct=out.crackDistinct&& ...
        diag.nUpper==diag.nLower&&diag.gridMismatch<1e-12;
end

function m=local_find_mid(T,maps,E)
    m=[];
    for k=1:3
        hit=all(sort(T(:,maps(k,1:2)),2)==sort(E(:).'),2);
        m=[m;T(hit,maps(k,3))]; %#ok<AGROW>
    end
    m=unique(m);
end

function [minJ,err]=local_jacobian_audit(P,T)
    err=0;
    for edge=[1 2 4;2 3 5;3 1 6].'
        err=max(err,max(vecnorm(P(T(:,edge(3)),:)- ...
            .5*(P(T(:,edge(1)),:)+P(T(:,edge(2)),:)),2,2)));
    end
    minJ=Inf;xi=local_rule16();
    for j=1:16
        [~,d]=local_shape(xi(:,j));
        j11=sum(reshape(P(T(:),1),size(T)).*d(1,:),2);
        j12=sum(reshape(P(T(:),2),size(T)).*d(1,:),2);
        j21=sum(reshape(P(T(:),1),size(T)).*d(2,:),2);
        j22=sum(reshape(P(T(:),2),size(T)).*d(2,:),2);
        minJ=min(minJ,min(j11.*j22-j12.*j21));
    end
end

function ids=local_q_support(mesh,cr,ri,ro,skipConstant)
    P=mesh.coord;T=mesh.connect;tip=cr.Pmid(end,:);
    e=cr.Pmid(end,:)-cr.Pmid(end-1,:);e=e/norm(e);R=[e(:),[-e(2);e(1)]];
    r=vecnorm((P-tip)*R,2,2);q=ones(size(r));q(r>=ro)=0;
    mid=r>ri&r<ro;q(mid)=(ro-r(mid))/(ro-ri);
    xip=local_rule16();use=false(size(T,1),1);

    for k=1:size(T,1)
        X=P(T(k,:),:);qe=q(T(k,:));
        if skipConstant&&all(qe==qe(1)),continue,end
        if all(qe==0),continue,end
        c=mean(X(1:3,:),1);
        if norm((c-tip)*R)>ro+max(vecnorm(X-c,2,2)),continue,end
        for j=1:16
            [N,d]=local_shape(xip(:,j));J=d*X;
            assert(det(J)>0,'pathqual:InvalidT6','Nonpositive T6 Jacobian.');
            if norm(N.'*X-tip)<1e-12,continue,end
            grad=R.'*((J\d)*qe);
            if norm(grad)>1e-14,use(k)=true;break,end
        end
    end
    ids=find(use);
end

function pass=local_boundary_same(P,T,Q,S,cr,crNew)
    [E,~]=local_edge_inventory(T);[F,~]=local_edge_inventory(S);
    E=E(~local_face_edge_mask(E(:,1:2),cr.upperNodes,cr.lowerNodes,cr.tipNode),1:2);
    F=F(~local_face_edge_mask(F(:,1:2),crNew.upperNodes,crNew.lowerNodes,crNew.tipNode),1:2);

    A=[P(E(:,1),:),P(E(:,2),:)];
    B=[Q(F(:,1),:),Q(F(:,2),:)];
    A=local_canonical_segments(A);
    B=local_canonical_segments(B);

    % Exact carrier-polygon coverage permits subdivision of its straight
    % segments. Original physical vertices are checked separately.
    assigned=zeros(size(B,1),1);
    intervals=zeros(size(B,1),2);
    pass=true;

    for j=1:size(B,1)
        v=A(:,3:4)-A(:,1:2);
        ell=sum(v.^2,2);
        u=B(j,1:2)-A(:,1:2);
        w=B(j,3:4)-A(:,1:2);
        t=sum(u.*v,2)./ell;
        s=sum(w.*v,2)./ell;

        match=vecnorm(u-t.*v,2,2)<1e-12 & ...
              vecnorm(w-s.*v,2,2)<1e-12 & ...
              t>=-1e-9 & t<=1+1e-9 & s>=-1e-9 & s<=1+1e-9;

        k=find(match,1);
        if isempty(k),pass=false;return,end
        assigned(j)=k;
        intervals(j,:)=sort([t(k),s(k)]);
    end

    for j=1:size(A,1)
        spans=sortrows(intervals(assigned==j,:));
        if isempty(spans) || abs(spans(1,1))>1e-9 || ...
                abs(spans(end,2)-1)>1e-9 || ...
                any(abs(spans(2:end,1)-spans(1:end-1,2))>1e-9)
            pass=false;
            return
        end
    end
end

function A=local_canonical_segments(A)
    flip=A(:,1)>A(:,3) | (A(:,1)==A(:,3) & A(:,2)>A(:,4));
    A(flip,:)=A(flip,[3 4 1 2]);
end

function tf=local_has_duplicate_triangles(T)
    C=sort(T,2);
    tf=size(unique(C,'rows'),1)~=size(C,1);
end

function [seamShared,exteriorOutside]=local_seam_audit(P,T,patchRows,extRows,inner)
    [~,inventory]=local_edge_inventory(T);

    [freePatch,~]=local_edge_inventory(T(patchRows,:));
    E=freePatch(:,1:2);
    rr=reshape(vecnorm(P(E(:),:),2,2),size(E));
    rp=max(vecnorm(inner,2,2));
    E=E(all(rr>rp-1e-12,2),:);

    [has,where]=ismember(sort(E,2),inventory(:,1:2),'rows');
    seamShared=all(has) && all(inventory(where(has),3)==2);

    cent=(P(T(extRows,1),:)+P(T(extRows,2),:)+P(T(extRows,3),:))/3;
    exteriorOutside=~any(inpolygon(cent(:,1),cent(:,2),inner(:,1),inner(:,2)));
end

function ratio=local_max_neighbor_size_ratio(T,L,affected)
    n=size(T,1);
    E=sort([T(:,[1 2]);T(:,[2 3]);T(:,[3 1])],2);
    which=repmat((1:n)',3,1);
    [~,~,g]=unique(E,'rows');
    lo=accumarray(g,L(which),[],@min);
    hi=accumarray(g,L(which),[],@max);
    touch=accumarray(g,ismember(which,affected),[],@max)>0;
    ratio=max(hi(touch)./lo(touch));
end

function xi=local_rule16()
    xi=zeros(2,16);xi(:,1)=[1/3;1/3];
    a=.170569307751760;b=.658861384496480;xi(:,2:4)=[a a b;a b a];
    a=.050547228317031;b=.898905543365938;xi(:,5:7)=[a a b;a b a];
    a=.459292588292723;b=.081414823414554;xi(:,8:10)=[a a b;a b a];
    a=.263112829634638;b=.728492392955404;c=.008394777409958;
    xi(:,11:16)=[a a b b c c;b c a c a b];
end

function [N,d]=local_shape(xi)
    a=xi(1);b=xi(2);c=1-a-b;
    N=[a*(2*a-1);b*(2*b-1);c*(2*c-1);4*a*b;4*b*c;4*c*a];
    d=[4*a-1,0,-(4*c-1),4*b,-4*b,4*(c-a); ...
       0,4*b-1,-(4*c-1),4*a,4*(c-b),-4*a];
end

function mat=local_material(C)
    E=C.E;nu=C.nu;ps=C.ps;
    if ps==1
        coef=E/((1+nu)*(1-2*nu));
        D=coef*[1-nu,nu,0;nu,1-nu,0;0,0,(1-2*nu)/2];
    else
        coef=E/(1-nu^2);
        D=coef*[1,nu,0;nu,1,0;0,0,(1-nu)/2];
    end
    mat=struct('E',E,'nu',nu,'ps',ps,'D',D,'Dmat',D);
end

function local_print_gates(S)
    fn=fieldnames(S);
    for k=1:numel(fn)
        fprintf('  %-38s : %d\n',fn{k},logical(S.(fn{k})));
    end
end

function local_plot_full_mesh(P,T,cr,patchRows,extRows,ri,ro,tip,R)
    X=(P-tip)*R;
    figure('Name','Incremental crack candidate','Color','w');
    tiledlayout(1,3,'Padding','compact','TileSpacing','compact');

    nexttile;hold on;axis equal;box on
    triplot(T,P(:,1),P(:,2));
    plot(cr.Pmid(:,1),cr.Pmid(:,2),'k-','LineWidth',1.5);
    title('Full candidate');xlabel('x');ylabel('y');

    nexttile;hold on;axis equal;box on
    patch('Faces',T(extRows,:),'Vertices',X*1e3,'FaceColor','none', ...
        'EdgeColor',[.65 .65 .65]);
    patch('Faces',T(patchRows,:),'Vertices',X*1e3,'FaceColor','none', ...
        'EdgeColor',[.15 .15 .15]);
    th=linspace(0,2*pi,300);
    plot(1e3*ri*cos(th),1e3*ri*sin(th),'--');
    plot(1e3*ro*cos(th),1e3*ro*sin(th),'--');
    xlim([-5 5]);ylim([-5 5]);title('Core and EDI support');
    xlabel('x_1 [mm]');ylabel('x_2 [mm]');

    nexttile;hold on;axis equal;box on
    patch('Faces',T(patchRows,:),'Vertices',X*1e3,'FaceColor','none');
    plot(1e3*X(cr.upperNodes,1),1e3*X(cr.upperNodes,2),'o');
    plot(1e3*X(cr.lowerNodes,1),1e3*X(cr.lowerNodes,2),'x');
    xlim([-8.2 .3]);ylim([-1 1]);title('Polyline path and new-tip 3-mm core');
    xlabel('x_1 [mm]');ylabel('x_2 [mm]');
end

function tf=local_is_absolute_path(p)
    p=char(p);
    tf=startsWith(p,filesep)|| ...
        ~isempty(regexp(p,'^[A-Za-z]:[\\/]','once'))||startsWith(p,'\\');
end

