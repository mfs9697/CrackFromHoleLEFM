function Core = build_stage2_scaled_audited_core(xTip, e1, a0, varargin)
%BUILD_STAGE2_SCALED_AUDITED_CORE
% Construct the reflection-paired structured crack-tip/core mesh inherited
% from the closed Step62 SIF audit, rescaled by the current trial-crack
% length a0 and rotated/translated into the current global crack frame.
%
% This function builds the CORE ONLY. It does not mesh the outer plate,
% assemble stiffness, solve equilibrium, or evaluate SIFs.
%
% Audit-derived nondimensional defaults (historical a0 = 8 mm):
%   hTip/a0      = 0.0540246508 / 8 = 0.00675308135
%   rCore/a0     = 6 / 8            = 0.75
%   h(r)/a0      = scale * [hTip/a0 + 0.028*(r/a0)]
%
% The topology is the exact Step62 three-sector radial zipper:
%   - 6 T3 triangles incident on the tip;
%   - 7 topological incident edges (duplicated crack edge);
%   - upper half explicitly built then reflected;
%   - positive crack-axis ligament shared;
%   - negative crack-axis faces topologically distinct.
%
% Inputs
% ------
%   xTip  [1x2] global crack-tip coordinate
%   e1    [1x2] global crack tangent, mouth -> tip
%   a0    trial-crack length
%
% Name-value
% ----------
%   'Scale'        family scale s (default 1)
%   'HTipOverA0'   audited base hTip/a0 (default 0.00675308135)
%   'RCoreOverA0'  audited paired-core radius/a0 (default 0.75)
%
% Output
% ------
%   Core.local.coord3/connect3   local T3 mesh
%   Core.local.coord/connect     local T6 mesh
%   Core.global.coord3/connect3  global T3 mesh
%   Core.global.coord/connect    global T6 mesh
%   Core.crack                   crack-face/topology metadata
%   Core.design                  audited nondimensional design metadata
%   Core.mirror                  T3/T6 reflection maps for the upper half

    ip = inputParser;
    addParameter(ip,'Scale',1,@(x)isnumeric(x)&&isscalar(x)&&isfinite(x)&&x>0&&x<=1);
    addParameter(ip,'HTipOverA0',0.00675308135,@(x)isnumeric(x)&&isscalar(x)&&isfinite(x)&&x>0);
    addParameter(ip,'RCoreOverA0',0.75,@(x)isnumeric(x)&&isscalar(x)&&isfinite(x)&&x>0);
    parse(ip,varargin{:});
    opt = ip.Results;

    validateattributes(xTip,{'numeric'},{'vector','numel',2,'finite','real'});
    validateattributes(e1,{'numeric'},{'vector','numel',2,'finite','real'});
    validateattributes(a0,{'numeric'},{'scalar','positive','finite','real'});

    xTip = reshape(double(xTip),1,2);
    e1 = reshape(double(e1),1,2);
    e1 = e1 / norm(e1);
    e2 = [-e1(2), e1(1)];

    hBase = opt.HTipOverA0 * a0;
    rCore = opt.RCoreOverA0 * a0;
    scale = opt.Scale;

    [Z,T,mirrorT3,nUpper,rings,axisIDs,design] = ...
        local_build_step62_patch(rCore,hBase,scale);

    [coord6,connect6] = T3toT6_fast(Z,T);

    nElemUpper = size(T,1)/2;
    assert(nElemUpper==round(nElemUpper), ...
        'stage2core:UpperLowerElementCount','Expected equal reflected halves.');

    mirrorT6 = local_build_t6_mirror(T,connect6,mirrorT3,nUpper,nElemUpper);

    % Global transform from local crack frame.
    Rgl = [e1(:),e2(:)];
    coord3Global = xTip + Z*Rgl.';
    coord6Global = xTip + coord6*Rgl.';

    % Crack faces are the duplicated negative local x1 axis.
    tipNode = 1;
    up3 = find((1:size(Z,1))'<=nUpper & Z(:,1)<=0 & abs(Z(:,2))<=64*eps(max(1,rCore)));
    up3 = unique([tipNode;up3(:)]);
    lo3 = mirrorT3(up3);
    [~,iu] = sort(vecnorm(Z(up3,:),2,2));
    [~,il] = sort(vecnorm(Z(lo3,:),2,2));
    up3 = up3(iu);
    lo3 = lo3(il);

    [up6,upMid] = local_face_t6_nodes(T,connect6,up3);
    [lo6,loMid] = local_face_t6_nodes(T,connect6,lo3);

    [~,iu6] = sort(vecnorm(coord6(up6,:),2,2));
    [~,il6] = sort(vecnorm(coord6(lo6,:),2,2));
    up6 = up6(iu6);
    lo6 = lo6(il6);

    Core = struct();
    Core.a0 = a0;
    Core.xTip = xTip;
    Core.e1 = e1;
    Core.e2 = e2;
    Core.Rgl = Rgl;
    Core.hBase = hBase;
    Core.hTip = scale*hBase;
    Core.rCore = rCore;
    Core.scale = scale;

    Core.local = struct( ...
        'coord3',Z, ...
        'connect3',T, ...
        'coord',coord6, ...
        'connect',connect6);

    Core.global = struct( ...
        'coord3',coord3Global, ...
        'connect3',T, ...
        'coord',coord6Global, ...
        'connect',connect6);

    Core.crack = struct();
    Core.crack.tipNode = tipNode;
    Core.crack.upperT3 = up3;
    Core.crack.lowerT3 = lo3;
    Core.crack.upperT6 = up6;
    Core.crack.lowerT6 = lo6;
    Core.crack.upperMidsideT6 = upMid;
    Core.crack.lowerMidsideT6 = loMid;

    Core.mirror = struct();
    Core.mirror.T3 = mirrorT3;
    Core.mirror.T6 = mirrorT6;
    Core.mirror.nUpperT3Nodes = nUpper;
    Core.mirror.nUpperT3Elements = nElemUpper;

    Core.design = design;
    Core.design.a0_m = a0;
    Core.design.hTipOverA0 = opt.HTipOverA0;
    Core.design.rCoreOverA0 = opt.RCoreOverA0;
    Core.design.hTip_m = scale*hBase;
    Core.design.rCore_m = rCore;
    Core.design.rings = rings;
    Core.design.axisIDs = axisIDs;
    Core.design.source = ...
        'closed Step62 SIF-audit three-sector reflection-paired core';
end


% =========================================================================
% Exact Step62 structured-patch construction, expressed with explicit scale.
% =========================================================================
function [Z,T,mirror,nUpper,rings,axisIDs,design] = ...
        local_build_step62_patch(rp,hBase,scale)

    rFirst = scale*hBase;
    assert(rp>4*rFirst,'stage2core:PatchSize', ...
        'Insufficient room for the audited graded rings.');

    slope = 0.028;
    radialFactor = sqrt(3)/2;

    metricSpan = log((hBase+slope*rp)/(hBase+slope*rFirst))/slope;
    nBands = ceil(metricSpan/(radialFactor*scale));
    ds = metricSpan/nBands;

    k = (0:nBands)';
    rings = ((hBase+slope*rFirst)*exp(slope*ds*k)-hBase)/slope;
    rings(1) = rFirst;
    rings(end) = rp;

    h = scale*(hBase+slope*rings);
    nTheta = 3*ceil(pi*rings./h/3);
    nTheta(1) = 3;

    assert(all(diff(nTheta)>=0) && all(mod(nTheta,3)==0), ...
        'stage2core:AngularLaw','Invalid audited angular-count progression.');
    assert(all(diff(diff(rings))>=-1e-15), ...
        'stage2core:RadialLaw','Audited ring widths are not monotone.');

    U = [0 0];
    ringIDs = cell(numel(rings),1);
    nodeRing = 0;
    nodeAngle = 0;

    for j = 1:numel(rings)
        th = (0:nTheta(j))'*pi/nTheta(j);
        xy = rings(j)*[cos(th),sin(th)];
        xy([1 end],2) = 0;

        ringIDs{j} = (size(U,1)+(1:size(xy,1)))';
        U = [U;xy]; %#ok<AGROW>
        nodeRing = [nodeRing;repmat(j,size(xy,1),1)]; %#ok<AGROW>
        nodeAngle = [nodeAngle;th]; %#ok<AGROW>
    end

    first = ringIDs{1};
    Tu = [ones(3,1),first(1:3),first(2:4)];
    bandRange = zeros(nBands+1,2);
    bandRange(1,:) = [1 3];

    for j = 2:numel(rings)
        start = size(Tu,1)+1;
        A = ringIDs{j-1};
        B = ringIDs{j};
        na = nTheta(j-1)/3;
        nb = nTheta(j)/3;

        for sector = 0:2
            aa = A(sector*na+(1:na+1));
            bb = B(sector*nb+(1:nb+1));
            i = 1;
            l = 1;

            while i<=na || l<=nb
                if i>na
                    advanceInner = false;
                elseif l>nb
                    advanceInner = true;
                else
                    di = sum((U(aa(i+1),:)-U(bb(l),:)).^2);
                    dout = sum((U(aa(i),:)-U(bb(l+1),:)).^2);
                    tol = 64*eps(max(di,dout));

                    if abs(di-dout)<=tol
                        advanceInner = mod(j+sector+i+l,2)==0;
                    else
                        advanceInner = di<dout;
                    end
                end

                if advanceInner
                    Tu(end+1,:) = [aa(i),bb(l),aa(i+1)]; %#ok<AGROW>
                    i = i+1;
                else
                    Tu(end+1,:) = [aa(i),bb(l),bb(l+1)]; %#ok<AGROW>
                    l = l+1;
                end
            end
        end

        bandRange(j,:) = [start,size(Tu,1)];
    end

    a = U(Tu(:,2),:)-U(Tu(:,1),:);
    b = U(Tu(:,3),:)-U(Tu(:,1),:);
    assert(all(a(:,1).*b(:,2)-a(:,2).*b(:,1)>0), ...
        'stage2core:UpperOrientation','Invalid upper structured strip.');
    assert(size(Tu,1)==3+sum(nTheta(1:end-1)+nTheta(2:end)), ...
        'stage2core:StripCount','Unexpected audited strip count.');

    nUpper = size(U,1);
    axisIDs = find(U(:,2)==0);
    shared = axisIDs(U(axisIDs,1)>=0);
    nonShared = setdiff((1:nUpper)',shared);

    mirror = (1:nUpper)';
    mirror(nonShared) = nUpper+(1:numel(nonShared))';

    Z = [U;U(nonShared,1),-U(nonShared,2)];
    T = [Tu;mirror(Tu(:,[1 3 2]))];

    width = [NaN;diff(rings)];
    design = struct();
    design.family = 'explicit three-sector radial zipper';
    design.scale = scale;
    design.hBase_m = hBase;
    design.slope = slope;
    design.radialFactor = radialFactor;
    design.metricStep = ds;
    design.upperRingIDs = ringIDs;
    design.upperNodeRing = nodeRing;
    design.upperNodeAngle = nodeAngle;
    design.upperBandElementRange = bandRange;
    design.angularIntervalsUpper = nTheta;
    design.tipTriangles = 6;
    design.tipIncidentTopologicalEdges = 7;
    design.hLaw = 'h(r)=scale*(hBase+0.028*r)';
    design.connectivityRule = ...
        'three sector strips; shortest new bridge; fixed parity ties';
    design.usesDelaunayInPairedRegion = false;
    design.usesRandomness = false;
    design.usesSmoothingInPairedRegion = false;
    design.allBandWidthsIncreasing = true;
    design.nRings = numel(rings);
    design.ringTable = table( ...
        (1:numel(rings))',rings*1e3,h*1e3,width*1e3,nTheta, ...
        pi*1e3*rings./nTheta, ...
        'VariableNames',{'ring','radius_mm','targetH_mm','bandWidth_mm', ...
        'upperAngularIntervals','arcSpacing_mm'});
end


function mirrorT6 = local_build_t6_mirror(T3,T6,mirrorT3,nUpper,nElemUpper)
% Build the upper-to-lower T6 reflection map from the T3 vertex map and the
% global T3->T6 edge mapping. Shared positive-axis nodes map to themselves.

    n6 = max(T6,[],'all');
    mirrorT6 = zeros(n6,1);
    mirrorT6(1:nUpper) = mirrorT3;

    E = [sort(T3(:,[1 2]),2),T6(:,4); ...
         sort(T3(:,[2 3]),2),T6(:,5); ...
         sort(T3(:,[3 1]),2),T6(:,6)];
    [edgePairs,ia] = unique(E(:,1:2),'rows');
    edgeMids = E(ia,3);

    upperElems = (1:nElemUpper)';
    upperEdges = [sort(T3(upperElems,[1 2]),2),T6(upperElems,4); ...
                  sort(T3(upperElems,[2 3]),2),T6(upperElems,5); ...
                  sort(T3(upperElems,[3 1]),2),T6(upperElems,6)];

    upperEdges = unique(upperEdges,'rows');

    for k = 1:size(upperEdges,1)
        ab = upperEdges(k,1:2);
        midU = upperEdges(k,3);
        % mirrorT3 is a column vector. MATLAB preserves that vector
        % orientation for vector indexing, so force the reflected edge back
        % to a 1x2 row before the row-wise lookup.
        abM = sort(mirrorT3(ab(:))).';
        [tf,j] = ismember(abM,edgePairs,'rows');
        assert(tf,'stage2core:MirrorEdge','Reflected T3 edge not found.');
        mirrorT6(midU) = edgeMids(j);
    end

    upperT6 = unique(T6(upperElems,:));
    assert(all(mirrorT6(upperT6)>0), ...
        'stage2core:MirrorT6Incomplete','Incomplete T6 mirror map.');

    % Reflected element connectivity including the T6 edge-node permutation.
    lowerPattern = T6(upperElems,[1 3 2 6 5 4]);
    expectedLower = reshape(mirrorT6(lowerPattern(:)),size(lowerPattern));
    actualLower = T6(nElemUpper+(1:nElemUpper),:);

    assert(isequal(expectedLower,actualLower), ...
        'stage2core:MirrorT6Connectivity', ...
        'T6 reflected connectivity is not exactly paired.');
end


function [nodes6,mids] = local_face_t6_nodes(T3,T6,face3)
    face3 = face3(:);

    [~,ord] = sort(face3); %#ok<ASGLU>
    faceSet = false(max(T3,[],'all'),1);
    faceSet(face3) = true;

    rows = [];
    edgeEnds = [1 2;2 3;3 1];
    midCols = [4 5 6];

    for j = 1:3
        a = T3(:,edgeEnds(j,1));
        b = T3(:,edgeEnds(j,2));
        hit = faceSet(a) & faceSet(b);
        rows = [rows;[a(hit),b(hit),T6(hit,midCols(j))]]; %#ok<AGROW>
    end

    assert(~isempty(rows),'stage2core:FaceEdges','No T6 crack-face edges found.');
    mids = unique(rows(:,3));
    nodes6 = unique([face3;mids]);
end
