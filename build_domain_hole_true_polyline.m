function D=build_domain_hole_true_polyline(Pmid,A,B,holes,w,varargin)
%BUILD_DOMAIN_HOLE_TRUE_POLYLINE
% Geometry description for a plate with a circular hole and a genuinely
% polyline sharp appended slit. This is the Stage III-C replacement for the
% historical straight-only appended-hole route.

    validateattributes(Pmid,{'numeric'},{'2d','ncols',2,'finite'});
    validateattributes(A,{'numeric'},{'scalar','positive','finite'});
    validateattributes(B,{'numeric'},{'scalar','positive','finite'});
    validateattributes(w,{'numeric'},{'scalar','positive','finite'});
    assert(size(Pmid,1)>=3,'stage3c:NeedPolyline', ...
        'Stage III-C requires mouth, prior tip, and new tip.');

    ip=inputParser;
    addParameter(ip,'corner_tol',1e-10,@(x)isnumeric(x)&&isscalar(x)&&x>0);
    addParameter(ip,'mouth_eps',[],@(x)isempty(x)||(isnumeric(x)&&isscalar(x)&&x>0));
    addParameter(ip,'epsMode','arclength',@(s)ischar(s)||isstring(s));
    addParameter(ip,'nArc',160,@(x)isnumeric(x)&&isscalar(x)&&x>=8);
    addParameter(ip,'orientation','cw',@(s)ischar(s)||isstring(s));
    addParameter(ip,'miter_limit',6,@(x)isnumeric(x)&&isscalar(x)&&x>1);
    parse(ip,varargin{:});
    opt=ip.Results;
    if isempty(opt.mouth_eps),mouth_eps=w;else,mouth_eps=opt.mouth_eps;end

    holes=normalize_holes(holes);
    assert(~isempty(holes),'stage3c:NoHole','At least one hole is required.');

    [itouch,dtouch]=detect_touch(Pmid(1,:),holes);
    assert(~isempty(itouch),'stage3c:NoTouchedHole', ...
        'The first path point must lie on a circular hole.');

    hk=holes{itouch};
    assert(strcmpi(strtrim(hk.type),'circle'),'stage3c:HoleType', ...
        'Stage III-C currently supports circular touched holes only.');

    outer=[0,-B;A,-B;A,B;0,B];
    loops=holes_to_loops_local(holes);

    [Vapp,Gapp]=build_appended_hole_polyline_loop( ...
        hk,Pmid,mouth_eps, ...
        'epsMode',opt.epsMode, ...
        'nArc',opt.nArc, ...
        'orientation',opt.orientation, ...
        'miter_limit',opt.miter_limit, ...
        'corner_tol',opt.corner_tol);
    loops{itouch}=Vapp;

    D=struct();
    D.outerPoly=outer;
    D.holeLoops=loops;
    D.channelPoly=[];
    D.Pmid=Pmid;
    D.A=A;D.B=B;D.w=w;D.holes=holes;

    D.channelGeom=struct();
    D.channelGeom.mode='merged_appended_hole_true_polyline';
    D.channelGeom.append=Gapp;

    D.topology=struct();
    D.topology.mode='merged_appended_hole_true_polyline';
    D.topology.channelTouchesHole=true;
    D.topology.touchingHoleIndex=itouch;
    D.topology.touchDistance=dtouch;
    D.topology.appendedHoleIndex=itouch;
    D.topology.appendedHoleArea=abs(signed_area(Vapp));
    D.topology.flags=struct( ...
        'tipInsidePlate',in_box(Pmid(end,:),A,B), ...
        'midlineInsidePlate',all(Pmid(:,1)>=-1e-12&Pmid(:,1)<=A+1e-12& ...
                                  Pmid(:,2)>=-B-1e-12&Pmid(:,2)<=B+1e-12), ...
        'startOnHole',true);

    D.options=struct('corner_tol',opt.corner_tol,'mouth_eps',mouth_eps, ...
        'epsMode',char(opt.epsMode),'nArc',opt.nArc, ...
        'orientation',char(opt.orientation),'miter_limit',opt.miter_limit);
end

function holes=normalize_holes(holes)
    if isstruct(holes),holes=num2cell(holes);end
    assert(iscell(holes),'stage3c:HoleInput','holes must be a cell/struct array.');
    for k=1:numel(holes)
        assert(isstruct(holes{k})&&isfield(holes{k},'type'), ...
            'stage3c:HoleSpec','Each hole needs a type field.');
    end
end

function [idx,dmin]=detect_touch(x,holes)
    idx=[];dmin=inf;
    for k=1:numel(holes)
        h=holes{k};
        if strcmpi(strtrim(h.type),'circle')
            d=abs(norm(x-h.center(:).')-h.r);
        else
            continue
        end
        if d<dmin,dmin=d;idx=k;end
    end
    tol=1e-8*max(1,norm(x));
    if dmin>tol,idx=[];end
end

function loops=holes_to_loops_local(holes)
    loops=cell(size(holes));
    for k=1:numel(holes)
        h=holes{k};
        switch lower(strtrim(h.type))
            case 'circle'
                n=160;
                if isfield(h,'npoly')&&~isempty(h.npoly),n=h.npoly;end
                th=linspace(0,2*pi,n+1).';th(end)=[];
                V=h.center(:).'+h.r*[cos(th),sin(th)];
                if signed_area(V)>0,V=flipud(V);end
                loops{k}=V;
            otherwise
                error('stage3c:UnsupportedHole','Unsupported hole type %s.',h.type);
        end
    end
end

function A=signed_area(P)
    x=P(:,1);y=P(:,2);
    A=.5*sum(x.*y([2:end 1])-x([2:end 1]).*y);
end

function tf=in_box(x,A,B)
    tf=x(1)>=-1e-12&&x(1)<=A+1e-12&&x(2)>=-B-1e-12&&x(2)<=B+1e-12;
end
