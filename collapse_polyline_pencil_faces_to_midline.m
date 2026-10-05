function Mc=collapse_polyline_pencil_faces_to_midline(M,D,varargin)
%COLLAPSE_POLYLINE_PENCIL_FACES_TO_MIDLINE
% Collapse all upper/lower geometry edges of a finite-width polyline pencil
% onto the complete crack midline while retaining separate node IDs/topology.

    ip=inputParser;
    addParameter(ip,'UpperEdgeIDs',[],@(x)isnumeric(x)&&~isempty(x));
    addParameter(ip,'LowerEdgeIDs',[],@(x)isnumeric(x)&&~isempty(x));
    addParameter(ip,'TipVertexID',[],@(x)isnumeric(x)&&isscalar(x));
    addParameter(ip,'Tol',1e-12,@(x)isnumeric(x)&&isscalar(x)&&x>0);
    parse(ip,varargin{:});
    opt=ip.Results;

    Pmid=D.Pmid;
    assert(size(Pmid,1)>=3,'stage3c:NeedPolyline','Polyline requires >=3 vertices.');

    up=[];
    for eid=opt.UpperEdgeIDs(:).'
        up=[up;findNodes(M.meshobj,'region','Edge',eid).']; %#ok<AGROW>
    end
    lo=[];
    for eid=opt.LowerEdgeIDs(:).'
        lo=[lo;findNodes(M.meshobj,'region','Edge',eid).']; %#ok<AGROW>
    end
    up=unique(up);lo=unique(lo);
    assert(~isempty(up)&&~isempty(lo),'stage3c:EmptyFaces','Empty face node set.');

    tip=unique(findNodes(M.meshobj,'region','Vertex',opt.TipVertexID));
    if isempty(tip)
        [~,tip]=min(vecnorm(M.p-Pmid(end,:),2,2));
    else
        tip=tip(1);
    end

    sU=polyline_parameter(M.p(up,:),Pmid,opt.Tol);
    sL=polyline_parameter(M.p(lo,:),Pmid,opt.Tol);

    p=M.p;
    p(up,:)=polyline_points(Pmid,sU);
    p(lo,:)=polyline_points(Pmid,sL);
    p(tip,:)=Pmid(end,:);

    % Orthogonal projection of a finite-width miter vertex onto a kinked
    % midline is not guaranteed to land on the path vertex itself. The
    % carrier metadata gives the exact corresponding upper/lower geometry
    % vertex for every Pmid vertex, so snap those specific mesh nodes to the
    % matching crack-path vertex after the general projection.
    G=D.channelGeom.append;
    assert(isfield(G,'face_upper_chain')&&isfield(G,'face_lower_chain'), ...
        'stage3c:MissingFaceChains','Polyline carrier face chains are required.');
    Uchain=G.face_upper_chain;
    Lchain=G.face_lower_chain;
    nVert=size(Pmid,1);
    assert(size(Uchain,1)==nVert&&size(Lchain,1)==nVert, ...
        'stage3c:FaceChainSize','Carrier face-chain/path-vertex counts differ.');

    snapUp=zeros(nVert,1);
    snapLo=zeros(nVert,1);
    snapTol=1e-8*max(1,norm(Pmid(end,:)-Pmid(1,:)));
    for k=1:nVert
        [du,iu]=min(vecnorm(M.p(up,:)-Uchain(k,:),2,2));
        [dl,il]=min(vecnorm(M.p(lo,:)-Lchain(k,:),2,2));
        assert(du<=snapTol&&dl<=snapTol, ...
            'stage3c:CarrierVertexNodeMissing', ...
            'Could not recover a mesh node for carrier path vertex %d.',k);
        snapUp(k)=up(iu);
        snapLo(k)=lo(il);
        p(snapUp(k),:)=Pmid(k,:);
        p(snapLo(k),:)=Pmid(k,:);
    end
    p(tip,:)=Pmid(end,:);

    % Recompute arc-length parameters after the exact vertex snaps, then
    % sort both faces consistently from mouth to current tip.
    sU=polyline_parameter(p(up,:),Pmid,opt.Tol);
    sL=polyline_parameter(p(lo,:),Pmid,opt.Tol);
    [sU,iu]=sort(sU);up=up(iu);
    [sL,il]=sort(sL);lo=lo(il);

    % Every historical path vertex before the current tip must now exist on
    % both faces after collapse. Interior upper/lower IDs remain distinct.
    vUp=cell(nVert,1);vLo=cell(nVert,1);
    for k=1:nVert
        vUp{k}=up(vecnorm(p(up,:)-Pmid(k,:),2,2)<=50*opt.Tol);
        vLo{k}=lo(vecnorm(p(lo,:)-Pmid(k,:),2,2)<=50*opt.Tol);
    end

    crack=struct();
    crack.x0=Pmid(1,:);crack.xtip=Pmid(end,:);crack.Pmid=Pmid;
    crack.upperEdgeIDs=opt.UpperEdgeIDs(:).';
    crack.lowerEdgeIDs=opt.LowerEdgeIDs(:).';
    crack.tipVertexID=opt.TipVertexID;crack.tipNode=tip;
    crack.upperNodes=up;crack.lowerNodes=lo;
    crack.upperS=sU;crack.lowerS=sL;
    crack.upperTarget=p(up,:);crack.lowerTarget=p(lo,:);
    crack.nUpper=numel(up);crack.nLower=numel(lo);
    crack.sameCount=numel(up)==numel(lo);
    crack.pathVertexUpperNodes=vUp;
    crack.pathVertexLowerNodes=vLo;
    crack.snappedUpperPathVertexNodes=snapUp;
    crack.snappedLowerPathVertexNodes=snapLo;

    Mc=M;
    Mc.p0=M.p;Mc.p=p;Mc.t=M.t;
    Mc.geom0=M.geom;Mc.meshobj0=M.meshobj;
    Mc.crack=crack;
    if ~isfield(Mc,'edgeSets'),Mc.edgeSets=struct();end
    Mc.edgeSets.crackUpper=up;
    Mc.edgeSets.crackLower=lo;
    Mc.edgeSets.crackTip=tip;
    if ~isfield(Mc,'region'),Mc.region=struct();end
    Mc.region.mode='collapsed_true_polyline_crack';
end

function s=polyline_parameter(X,P,tol)
    seg=diff(P,1,1);L=vecnorm(seg,2,2);cum=[0;cumsum(L)];Lt=cum(end);
    assert(Lt>tol&&all(L>tol),'stage3c:DegeneratePolyline','Degenerate polyline.');
    s=zeros(size(X,1),1);
    for i=1:size(X,1)
        best=inf;bs=0;
        for k=1:numel(L)
            v=seg(k,:);tt=dot(X(i,:)-P(k,:),v)/dot(v,v);
            tt=max(0,min(1,tt));q=P(k,:)+tt*v;
            d=norm(X(i,:)-q);
            if d<best,best=d;bs=(cum(k)+tt*L(k))/Lt;end
        end
        s(i)=bs;
    end
end

function X=polyline_points(P,s)
    seg=diff(P,1,1);L=vecnorm(seg,2,2);cum=[0;cumsum(L)];Lt=cum(end);
    X=zeros(numel(s),2);
    for i=1:numel(s)
        z=max(0,min(1,s(i)))*Lt;
        if abs(z-Lt)<=1e-14*max(1,Lt),X(i,:)=P(end,:);continue,end
        k=find(cum<=z,1,'last');if k>=numel(cum),k=numel(cum)-1;end
        t=(z-cum(k))/L(k);
        X(i,:)=P(k,:)+t*seg(k,:);
    end
end
