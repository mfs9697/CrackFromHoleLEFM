function ids=identify_polyline_pencil_edge_sets(M,D,varargin)
%IDENTIFY_POLYLINE_PENCIL_EDGE_SETS
% Recover every PDE geometry edge belonging to the upper/lower finite-width
% polyline pencil faces after meshing. Unlike the historical helper, this
% supports multiple face edges on each side.

    ip=inputParser;
    addParameter(ip,'Tol',[],@(x)isempty(x)||(isnumeric(x)&&isscalar(x)&&x>0));
    addParameter(ip,'Verbose',false,@(x)islogical(x)&&isscalar(x));
    parse(ip,varargin{:});
    opt=ip.Results;

    assert(isfield(M,'meshobj')&&isfield(M,'geom')&&isfield(M,'p'), ...
        'stage3c:MeshFields','Mesh object/geometry/coordinates are required.');
    assert(isfield(D,'channelGeom')&&isfield(D.channelGeom,'append'), ...
        'stage3c:AppendMeta','Polyline append metadata are required.');

    G=D.channelGeom.append;
    U=G.face_upper_chain;
    L=G.face_lower_chain;
    tip=G.xtip(:).';

    scale=max([1;vecnorm(D.Pmid-tip,2,2)]);
    if isempty(opt.Tol),tol=1e-9*scale;else,tol=opt.Tol;end

    geom=M.geom.Geometry;
    nEdges=geom.NumEdges;
    p=M.p;

    upper=[];lower=[];
    for eid=1:nEdges
        nd=unique(findNodes(M.meshobj,'region','Edge',eid));
        if isempty(nd),continue,end
        X=p(nd,:);
        dU=point_polyline_distance(X,U);
        dL=point_polyline_distance(X,L);
        if max(dU)<=50*tol && max(dU)<=max(dL)
            upper(end+1)=eid; %#ok<AGROW>
        elseif max(dL)<=50*tol
            lower(end+1)=eid; %#ok<AGROW>
        end
    end

    upper=unique(upper);lower=unique(lower);
    upper=setdiff(upper,lower);
    lower=setdiff(lower,upper);

    assert(numel(upper)>=size(U,1)-1,'stage3c:UpperEdges', ...
        'Did not recover all upper polyline face edges.');
    assert(numel(lower)>=size(L,1)-1,'stage3c:LowerEdges', ...
        'Did not recover all lower polyline face edges.');

    V=geom.Vertices;
    if size(V,1)==2,V=V.';end
    [dTip,vTip]=min(vecnorm(V-tip,2,2));
    assert(dTip<=100*tol,'stage3c:TipVertex','Could not recover sharp-tip vertex.');

    ids=struct('upperEdges',upper,'lowerEdges',lower,'v_tip',vTip, ...
        'xtip',tip,'tol',tol);

    if opt.Verbose
        fprintf('Stage III-C polyline face IDs: upper=%s lower=%s tipVertex=%d\n', ...
            mat2str(upper),mat2str(lower),vTip);
    end
end

function d=point_polyline_distance(X,P)
    d=inf(size(X,1),1);
    for k=1:size(P,1)-1
        A=P(k,:);B=P(k+1,:);v=B-A;L2=max(dot(v,v),1e-30);
        t=((X-A)*v.')/L2;t=max(0,min(1,t));
        Q=A+t.*v;
        d=min(d,vecnorm(X-Q,2,2));
    end
end
