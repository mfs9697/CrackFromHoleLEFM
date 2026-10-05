function [r,app,diag]=native_COD_polyline_audit(mesh,U,mat,crack,minPts,includeFace)
%NATIVE_COD_POLYLINE_AUDIT
% Polyline-safe wrapper around the qualified historical native_COD_audit.
%
% The historical extractor defines its local frame from crack.Pmid(1,:) to
% crack.Pmid(end,:), which is correct only for a straight crack. For an
% incremental polyline path, the asymptotic crack-tip frame must instead be
% defined by the LAST segment. This wrapper constructs an extraction-only
% crack descriptor containing [P_{k-1};P_k] while preserving the full
% upper/lower face topology and tip node.
%
% No geometry, connectivity, displacement, or COD formula is changed.

    if nargin<5||isempty(minPts),minPts=8;end
    if nargin<6||isempty(includeFace),includeFace=false;end

    assert(isstruct(crack)&&isfield(crack,'Pmid')&&size(crack.Pmid,1)>=2, ...
        'nativeCODpoly:BadCrack','crack.Pmid must contain at least two points.');

    P=crack.Pmid;
    last=P(end,:)-P(end-1,:);
    assert(norm(last)>1e-14,'nativeCODpoly:DegenerateLastSegment', ...
        'Last crack segment is degenerate.');

    tipCrack=crack;
    tipCrack.Pmid=P(end-1:end,:);
    tipCrack.x0=P(end-1,:);
    tipCrack.xtip=P(end,:);

    [r,app,diag]=native_COD_audit(mesh,U,mat,tipCrack,minPts,includeFace);
    diag.fullPathVertexCount=size(P,1);
    diag.usesLastSegmentFrame=true;
    diag.lastSegmentLength=norm(last);
    diag.lastSegmentDirection=last/norm(last);
end
