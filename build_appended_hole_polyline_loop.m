function [Vapp,G]=build_appended_hole_polyline_loop(hole,Pmid,eps,varargin)
%BUILD_APPENDED_HOLE_POLYLINE_LOOP
% Build one circular-hole inner loop with a genuinely polyline sharp pencil.
%
% Pmid is ordered mouth -> ... -> tip. The two finite-width pencil faces
% follow the complete polyline through mitered offset vertices and taper
% linearly from the circular-hole mouth shift to zero at the sharp tip.
% The geometry is only a meshing carrier: the two faces are later collapsed
% onto Pmid while retaining distinct upper/lower topology.

    ip=inputParser;
    addParameter(ip,'epsMode','arclength',@(s)ischar(s)||isstring(s));
    addParameter(ip,'nArc',160,@(x)isnumeric(x)&&isscalar(x)&&x>=8);
    addParameter(ip,'orientation','cw',@(s)ischar(s)||isstring(s));
    addParameter(ip,'miter_limit',6,@(x)isnumeric(x)&&isscalar(x)&&x>1);
    addParameter(ip,'corner_tol',1e-10,@(x)isnumeric(x)&&isscalar(x)&&x>0);
    parse(ip,varargin{:});
    opt=ip.Results;

    validateattributes(Pmid,{'numeric'},{'2d','ncols',2,'finite'});
    assert(size(Pmid,1)>=3, ...
        'stage3c:NeedPolyline','Stage III-C requires at least three midline vertices.');

    must_field(hole,'type');must_field(hole,'center');must_field(hole,'r');
    assert(strcmpi(strtrim(hole.type),'circle'), ...
        'stage3c:HoleType','Only circular holes are supported.');
    assert(isnumeric(eps)&&isscalar(eps)&&isfinite(eps)&&eps>0, ...
        'stage3c:BadMouthEps','eps must be positive.');

    c=hole.center(:).';
    r=hole.r;
    A=Pmid(1,:);

    rc=A-c;
    assert(norm(rc)>0,'stage3c:BadMouth','Crack mouth coincides with hole center.');
    A0=c+r*rc/norm(rc);

    % Snap only the mathematical mouth to the exact circle; the supplied
    % path must already agree with the frozen initiation point.
    P=Pmid;
    P(1,:)=A0;

    seg=diff(P,1,1);
    L=vecnorm(seg,2,2);
    assert(all(L>1e-12),'stage3c:DegenerateLeg','Polyline contains a degenerate leg.');
    T=seg./L;
    N=[-T(:,2),T(:,1)];

    % Mouth points are defined on the actual circle, as in the accepted
    % straight appended-hole construction.
    phiA=atan2(A0(2)-c(2),A0(1)-c(1));
    tcirc=[-sin(phiA),cos(phiA)];
    firstNormal=N(1,:);
    if dot(tcirc,firstNormal)>=0
        sgn=1;
    else
        sgn=-1;
    end

    switch lower(char(opt.epsMode))
        case 'arclength'
            dphi=eps/r;
        case 'angle'
            dphi=eps;
        otherwise
            error('stage3c:BadEpsMode','epsMode must be arclength or angle.');
    end
    assert(dphi>0&&dphi<pi/2,'stage3c:BadDphi','Invalid mouth half-angle.');

    M1=c+r*[cos(phiA+sgn*dphi),sin(phiA+sgn*dphi)];
    M2=c+r*[cos(phiA-sgn*dphi),sin(phiA-sgn*dphi)];

    if dot(M1-A0,firstNormal)>=dot(M2-A0,firstNormal)
        Mup=M1;Mlo=M2;
        phiUp=phiA+sgn*dphi;phiLo=phiA-sgn*dphi;
    else
        Mup=M2;Mlo=M1;
        phiUp=phiA-sgn*dphi;phiLo=phiA+sgn*dphi;
    end

    n=size(P,1);
    cum=[0;cumsum(L)];
    Ltot=cum(end);
    width=eps*(1-cum/Ltot);
    width(end)=0;

    Up=zeros(n,2);Lo=zeros(n,2);
    Up(1,:)=Mup;Lo(1,:)=Mlo;
    Up(end,:)=P(end,:);Lo(end,:)=P(end,:);

    for k=2:n-1
        n1=N(k-1,:);n2=N(k,:);
        b=n1+n2;
        if norm(b)<opt.corner_tol
            % Near 180-degree reversal is not an admissible crack increment.
            error('stage3c:Hairpin','Polyline contains a near-180-degree turn.');
        end
        b=b/norm(b);
        den=dot(b,n1);
        if abs(den)<opt.corner_tol
            error('stage3c:BadMiter','Ill-conditioned crack-face miter.');
        end
        m=width(k)/den;
        if abs(m)>opt.miter_limit*max(width(k),eps*1e-6)
            error('stage3c:MiterLimit','Crack-face miter exceeds the allowed limit.');
        end
        Up(k,:)=P(k,:)+m*b;
        Lo(k,:)=P(k,:)-m*b;
    end

    % Long retained hole arc excludes the removed mouth neighborhood.
    phiArc=long_arc(phiUp,phiLo,phiA,opt.nArc);
    arcMain=c+r*[cos(phiArc),sin(phiArc)];

    % Inner-loop path:
    % Mup -> retained hole -> Mlo -> lower face -> tip -> upper face -> Mup.
    Vapp=[arcMain;Lo(2:end);flipud(Up(1:end-1))];
    Vapp=remove_consecutive(Vapp,1e-12);
    if norm(Vapp(end,:)-Vapp(1,:))<1e-12
        Vapp(end,:)=[];
    end

    Aloop=signed_area(Vapp);
    switch lower(char(opt.orientation))
        case 'cw'
            if Aloop>0,Vapp=flipud(Vapp);end
        case 'ccw'
            if Aloop<0,Vapp=flipud(Vapp);end
        otherwise
            error('stage3c:Orientation','orientation must be cw or ccw.');
    end

    G=struct();
    G.mode='merged_appended_hole_true_polyline';
    G.center=c;G.r=r;G.A=A0;G.A_input=A;
    G.Pmid=P;
    G.xtip=P(end,:);
    G.Mup=Mup;G.Mlo=Mlo;
    G.phiA=phiA;G.dphi=dphi;
    G.arcMain=arcMain;
    G.face_upper_chain=Up;
    G.face_lower_chain=Lo;
    % Last face segments preserve compatibility with the existing automatic
    % tip-edge finder; all face edges are identified separately after mesh.
    G.face_upper=Up(end-1:end,:);
    G.face_lower=flipud(Lo(end-1:end,:));
    G.segmentTangents=T;
    G.segmentNormals=N;
    G.widthProfile=width;
    G.miterLimit=opt.miter_limit;
    G.orientation=char(opt.orientation);
end

function phi=long_arc(a,b,x,n)
    a=wrap2(a);b=wrap2(b);x=wrap2(x);
    bi=b;while bi<a,bi=bi+2*pi;end
    xi=x;while xi<a,xi=xi+2*pi;end
    if xi>=a&&xi<=bi
        bd=b;while bd>a,bd=bd-2*pi;end
        phi=linspace(a,bd,n).';
    else
        phi=linspace(a,bi,n).';
    end
end

function x=wrap2(x)
    x=mod(x,2*pi);if x<0,x=x+2*pi;end
end

function P=remove_consecutive(P,tol)
    keep=[true;vecnorm(diff(P,1,1),2,2)>tol];
    P=P(keep,:);
end

function A=signed_area(P)
    x=P(:,1);y=P(:,2);
    A=.5*sum(x.*y([2:end 1])-x([2:end 1]).*y);
end

function must_field(S,n)
    if ~isfield(S,n)||isempty(S.(n))
        error('stage3c:MissingField','Missing field %s.',n);
    end
end
