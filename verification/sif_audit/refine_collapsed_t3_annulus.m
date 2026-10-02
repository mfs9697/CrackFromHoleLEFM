function [Pnew,Tnew,crackNew,info]=refine_collapsed_t3_annulus( ...
    P,T,crack,rInner,rOuter,targetEdge,varargin)
%REFINE_COLLAPSED_T3_ANNULUS
% Conforming local edge-bisection refinement on an ALREADY COLLAPSED T3
% crack mesh. Coordinates of all ORIGINAL vertices are untouched.
% Duplicate crack faces have different IDs and therefore refine as
% separate topological edges, retaining the crack discontinuity.
%
% Refinement is driven by physical triangle locations and measured edge
% lengths, not by global PDE Toolbox Hmax. Split only triangle edges
% longer than targetEdge in a circle covering the entire EDI annulus,
% including a buffer. Neighbors with marked edges are handled by
% 1-/2-/3-edge conforming templates (no hanging nodes).
%
% Every new point is an EXACT edge midpoint in the already-collapsed
% polygon mesh. Physical boundary geometry and crack path are thus
% unchanged as piecewise-linear curves, although the T3 triangulation
% changes. Existing PDE mesh/geometry objects are archival only.
%
% Output face lists (crackNew.upperNodes/lowerNodes) include any newly
% split crack-boundary T3 vertices. New T6 midsides are classified later
% using the refined T3 topology.
%
% Intended only for the current 2D straight appended-hole LEFM audit.
%
% Example:
%  [p,t,cr,info]=refine_collapsed_t3_annulus( ...
%      mesh.coord3,mesh.connect3,crack,0.0008,0.0064,0.0004, ...
%      'OuterBuffer',0.0006,'MaxPasses',2);

ip=inputParser;
addParameter(ip,'OuterBuffer',0.0006, ...
    @(x)isnumeric(x)&&isscalar(x)&&isfinite(x)&&x>=0);
addParameter(ip,'InnerBuffer',0.0008, ...
    @(x)isnumeric(x)&&isscalar(x)&&isfinite(x)&&x>=0);
addParameter(ip,'MaxPasses',2, ...
    @(x)isnumeric(x)&&isscalar(x)&&isfinite(x)&&x>=1&&x<=5&&x==round(x));
addParameter(ip,'Verbose',true,@(x)islogical(x)||isnumeric(x));
parse(ip,varargin{:});
O=ip.Results;

if size(P,2)~=2||size(T,2)~=3|| ...
        any(~isfinite(P(:)))|| ...
        any(T(:)<1)||any(T(:)>size(P,1))||any(T(:)~=round(T(:)))
    error('annulus_refine:InvalidMesh','Requires finite 2D T3 data.');
end
fields={'Pmid','tipNode','upperNodes','lowerNodes'};
for k=1:numel(fields)
    if ~isfield(crack,fields{k})
        error('annulus_refine:MissingCrack', ...
            'crack.%s is required.',fields{k});
    end
end
if size(crack.Pmid,1)~=2 || ...
        ~(rInner>0&&rOuter>rInner&&targetEdge>0)
    error('annulus_refine:Geometry', ...
        'Requires straight crack and 0<rInner<rOuter, targetEdge>0.');
end
tip=crack.Pmid(end,:);
a0=norm(diff(crack.Pmid,1,1));
e1=diff(crack.Pmid,1,1)/a0;
startArea=total_t3_area(P,T);
if ~all(triangle_signed_twice_area(P,T)>0)
    error('annulus_refine:OriginalOrientation', ...
        'Original collapsed T3 elements must be positively oriented.');
end

Pnew=P;
Tnew=T;
crackNew=crack;
passes=nan(O.MaxPasses,7);
for pass=1:O.MaxPasses
    nt=size(Tnew,1);
    neP=size(Pnew,1);
    E=[Tnew(:,[1 2]);Tnew(:,[2 3]);Tnew(:,[3 1])];
    [ue,~,which]=unique(sort(E,2),'rows');
    ie=reshape(which,nt,3);
    L=hypot(Pnew(ue(:,1),1)-Pnew(ue(:,2),1), ...
            Pnew(ue(:,1),2)-Pnew(ue(:,2),2));
    r=hypot(Pnew(:,1)-tip(1),Pnew(:,2)-tip(2));
    rv=r(Tnew);
    C=(Pnew(Tnew(:,1),:)+Pnew(Tnew(:,2),:)+ ...
        Pnew(Tnew(:,3),:))/3;
    rc=hypot(C(:,1)-tip(1),C(:,2)-tip(2));
    inRegion=(rc<=rOuter+O.OuterBuffer) & ...
        (max(rv,[],2)>=max(0,rInner-O.InnerBuffer));

    % Include any triangle overlapping the annulus, not only centroid-
    % selected triangles. This prevents an oversized crossing element
    % from escaping refinement because its centroid lies just outside.
    overlap=(min(rv,[],2)<=rOuter+O.OuterBuffer) & ...
        (max(rv,[],2)>=max(0,rInner-O.InnerBuffer));
    useTri=inRegion | overlap;

    long=(L(ie)>targetEdge);
    markTri=useTri & any(long,2);
    if ~any(markTri)
        passes(pass,:)=[pass,neP,nt,0,0,0,0];
        if logical(O.Verbose)
            fprintf(['  annulus edge refinement pass %d: ', ...
                'no marked edges; target reached.\n'],pass);
        end
        break;
    end

    % Mark all THREE edges of each selected triangle. This produces
    % well-shaped red elements there; adjacent triangles are split by
    % conforming green 1-/2-edge templates if their edges are shared.
    selected=ie(markTri,:);
    marked=false(size(ue,1),1);
    marked(unique(selected(:)))=true;
    mids=zeros(size(ue,1),1);
    edgeSelected=find(marked);
    mids(edgeSelected)=neP+(1:numel(edgeSelected)).';
    pmid=0.5*(Pnew(ue(edgeSelected,1),:)+ ...
              Pnew(ue(edgeSelected,2),:));
    Pnext=[Pnew;pmid];

    % Update original crack-boundary T3 face nodes. Each face's identity
    % is topological; no coordinate-based merging across the two faces.
    oldUp=unique(crackNew.upperNodes(:));
    oldLo=unique(crackNew.lowerNodes(:));
    commonTip=crackNew.tipNode;
    eA=ue(edgeSelected,1);
    eB=ue(edgeSelected,2);
    upMark=ismember(eA,[oldUp;commonTip]) & ...
           ismember(eB,[oldUp;commonTip]);
    loMark=ismember(eA,[oldLo;commonTip]) & ...
           ismember(eB,[oldLo;commonTip]);
    % Never label an interior chord as a crack face: midpoint must lie
    % on the negative local crack axis and between mouth and tip.
    tMid=pmid-tip;
    sMid=tMid*e1.';
    tCross=tMid(:,1)*(-e1(2))+tMid(:,2)*e1(1);
    faceGeom=abs(tCross)<max(1e-12,1e-8*a0) & ...
        sMid<=1e-10 & sMid>=-a0-1e-10;
    upMark=upMark & faceGeom;
    loMark=loMark & faceGeom;
    if any(upMark&loMark)
        error('annulus_refine:CrossFace', ...
            'New edge midpoint belongs simultaneously to both crack faces.');
    end
    crackNew.upperNodes=unique([oldUp;mids(edgeSelected(upMark))]);
    crackNew.lowerNodes=unique([oldLo;mids(edgeSelected(loMark))]);

    Tnext=zeros(4*nt,3);
    outN=0;
    for k=1:nt
        v=Tnew(k,:);
        m=mids(ie(k,:)); % m12,m23,m31 in original CCW order
        ns=nnz(m);
        switch ns
            case 0
                outN=outN+1;
                Tnext(outN,:)=v;
            case 1
                s=find(m,1);
                rr=mod((s-1)+(0:2),3)+1;
                a=v(rr(1));b=v(rr(2));d=v(rr(3));
                mm=m(s);
                Tnext(outN+(1:2),:)=[a mm d;mm b d];
                outN=outN+2;
            case 2
                s=NaN;
                for ii=1:3
                    j=mod(ii,3)+1;
                    if m(ii)>0&&m(j)>0
                        s=ii;
                        break;
                    end
                end
                if ~isfinite(s)
                    error('annulus_refine:TwoEdgeLogic', ...
                        'Could not identify the shared vertex.');
                end
                rr=mod((s-1)+(0:2),3)+1;
                a=v(rr(1));b=v(rr(2));d=v(rr(3));
                mab=m(s);
                mbd=m(mod(s,3)+1);
                % Common corner at b, remaining quadrilateral
                % [a,mab,mbd,d]. Choose its shorter diagonal.
                corner=[b mbd mab];
                diag1=sum((Pnext(mab,:)-Pnext(d,:)).^2);
                diag2=sum((Pnext(a,:)-Pnext(mbd,:)).^2);
                if diag1<=diag2
                    fill=[a mab d;mab mbd d];
                else
                    fill=[a mab mbd;a mbd d];
                end
                Tnext(outN+(1:3),:)=[corner;fill];
                outN=outN+3;
            case 3
                a=v(1);b=v(2);d=v(3);
                m12=m(1);m23=m(2);m31=m(3);
                Tnext(outN+(1:4),:)=[ ...
                    a m12 m31;
                    m12 b m23;
                    m31 m23 d;
                    m12 m23 m31];
                outN=outN+4;
        end
    end
    Tnext=Tnext(1:outN,:);
    signed=triangle_signed_twice_area(Pnext,Tnext);
    if any(~isfinite(signed)) || any(signed<=1e-22)
        error('annulus_refine:InvertedTriangle', ...
            'Refinement produced an inverted/degenerate T3 triangle.');
    end
    newArea=sum(signed)/2;
    if abs(newArea-startArea)>1e-11*max(1,startArea)
        error('annulus_refine:AreaMismatch', ...
            'Refinement changed area: old %.15g, new %.15g.', ...
            startArea,newArea);
    end
    passes(pass,:)=[pass,neP,nt,numel(edgeSelected),nnz(markTri), ...
        size(Pnext,1),size(Tnext,1)];
    if logical(O.Verbose)
        fprintf(['  annulus edge refinement pass %d: ', ...
            '%d selected T3, %d split edges; ', ...
            'T3 %d -> %d, vertices %d -> %d\n'], ...
            pass,nnz(markTri),numel(edgeSelected), ...
            nt,size(Tnext,1),neP,size(Pnext,1));
    end
    Pnew=Pnext;
    Tnew=Tnext;
end

% Metadata is consistently recomputed for the final T3 face boundary.
for side=1:2
    if side==1
        nodes=crackNew.upperNodes(:);
        key='upper';
    else
        nodes=crackNew.lowerNodes(:);
        key='lower';
    end
    proj=(Pnew(nodes,:)-crack.Pmid(1,:))*e1.';
    [val,idx]=sort(proj);
    nodes=nodes(idx);
    if side==1
        crackNew.upperNodes=nodes;
        crackNew.upperS=val/a0;
        crackNew.upperTarget=Pnew(nodes,:);
        crackNew.nUpper=numel(nodes);
    else
        crackNew.lowerNodes=nodes;
        crackNew.lowerS=val/a0;
        crackNew.lowerTarget=Pnew(nodes,:);
        crackNew.nLower=numel(nodes);
    end
end
shared=setdiff(intersect(crackNew.upperNodes,crackNew.lowerNodes), ...
    crackNew.tipNode);
if ~isempty(shared)
    error('annulus_refine:SharedCrackFaces', ...
        'Refinement merged %d upper/lower crack-face vertices.', ...
        numel(shared));
end
crackNew.sameCount=(crackNew.nUpper==crackNew.nLower);
crackNew.lowerMatchForUpper=zeros(crackNew.nUpper,1);
for j=1:crackNew.nUpper
    [~,crackNew.lowerMatchForUpper(j)]=min( ...
        abs(crackNew.lowerS-crackNew.upperS(j)));
end
info=struct('passLog',passes(all(isfinite(passes),2),:), ...
    'initialVertices',size(P,1),'finalVertices',size(Pnew,1), ...
    'initialTriangles',size(T,1),'finalTriangles',size(Tnew,1), ...
    'originalArea',startArea,'finalArea',total_t3_area(Pnew,Tnew), ...
    'rInner',rInner,'rOuter',rOuter,'targetEdge',targetEdge, ...
    'outerBuffer',O.OuterBuffer,'innerBuffer',O.InnerBuffer);
end

function a=triangle_signed_twice_area(P,T)
p1=P(T(:,1),:);p2=P(T(:,2),:);p3=P(T(:,3),:);
u=p2-p1;v=p3-p1;
a=u(:,1).*v(:,2)-u(:,2).*v(:,1);
end

function a=total_t3_area(P,T)
a=sum(triangle_signed_twice_area(P,T))/2;
end
