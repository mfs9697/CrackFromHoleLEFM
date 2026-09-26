function M = SIF_mesh_mirror_metrics(mesh, V, rI, varargin)
%SIF_MESH_MIRROR_METRICS  Quantify local mesh mismatch at mirrored P/Q points.
%
% This diagnostic is deliberately independent of the displacement field.
% It samples the same mirror-point pairs used by SIF_LEFM_circle2_debug and
% compares the parent T3 elements containing P and Q.
%
% Primary asymmetry metric
%   A_h = median(abs(log(h_P/h_Q)))
% where h=sqrt(2*A_tri) is a characteristic element size. A_h=0 means
% identical local element size at every pair; for small mismatch it is
% approximately the relative size difference.
%
% Secondary diagnostics
%   centroid_mismatch = distance between the P-element centroid and the
%     mirror of the Q-element centroid, normalized by mean(h_P,h_Q)
%   bary_margin_mismatch = |min(baryP)-min(baryQ)|
%   quality_mismatch = |q_P-q_Q|, q=4*sqrt(3)*A/sum(edge_length^2)
%
% Inputs
%   mesh  .coord3, .connect3
%   V     crack polyline; last segment defines the local tip frame
%   rI    contour radius
%
% Name-value options
%   'nthet'     number of upper-half contour points, default 200
%   'eps_th'    angular exclusion from crack faces, default 1e-3
%   'theta'     optional explicit theta vector in radians
%   'edge_tol'  barycentric near-edge threshold, default 1e-4

    ip = inputParser;
    addParameter(ip,'nthet',200,@(x)isnumeric(x)&&isscalar(x)&&x>=10);
    addParameter(ip,'eps_th',1e-3,@(x)isnumeric(x)&&isscalar(x)&&x>0);
    addParameter(ip,'theta',[],@(x)isnumeric(x));
    addParameter(ip,'edge_tol',1e-4,@(x)isnumeric(x)&&isscalar(x)&&x>0);
    parse(ip,varargin{:});

    must_have(mesh,'coord3');
    must_have(mesh,'connect3');

    coord3 = mesh.coord3;
    conn3 = mesh.connect3;

    if size(V,1)<2 || size(V,2)~=2
        error('SIF_mesh_mirror_metrics:BadV','V must be [n x 2], n>=2.');
    end

    tip = V(end,:).';
    e1 = tip - V(end-1,:).';
    if norm(e1)<=eps
        error('SIF_mesh_mirror_metrics:DegenerateLastSegment', ...
            'Final crack segment has zero length.');
    end
    e1 = e1/norm(e1);
    e2 = [-e1(2);e1(1)];
    Rgl = [e1,e2];
    Rloc = Rgl.';

    if isempty(ip.Results.theta)
        th = linspace(ip.Results.eps_th,pi-ip.Results.eps_th, ...
            round(ip.Results.nthet)).';
    else
        th = ip.Results.theta(:);
    end

    PLoc = [rI*cos(th), rI*sin(th)];
    QLoc = [rI*cos(th),-rI*sin(th)];
    P = (Rgl*PLoc.').'+tip.';
    Q = (Rgl*QLoc.').'+tip.';

    TR = triangulation(conn3,coord3);
    [eP,bP] = pointLocation(TR,P(:,1),P(:,2));
    [eQ,bQ] = pointLocation(TR,Q(:,1),Q(:,2));

    if any(isnan(eP)) || any(isnan(eQ))
        error('SIF_mesh_mirror_metrics:ContourOutsideMesh', ...
            'One or more mirrored contour points lie outside the mesh.');
    end

    n = numel(th);
    hP = zeros(n,1); hQ = zeros(n,1);
    qP = zeros(n,1); qQ = zeros(n,1);
    cP = zeros(n,2); cQ = zeros(n,2);

    for k=1:n
        [hP(k),qP(k),cP(k,:)] = tri_metrics(coord3(conn3(eP(k),:),:));
        [hQ(k),qQ(k),cQ(k,:)] = tri_metrics(coord3(conn3(eQ(k),:),:));
    end

    cPloc = (Rloc*(cP-tip.').').';
    cQloc = (Rloc*(cQ-tip.').').';
    cQmir = [cQloc(:,1),-cQloc(:,2)];

    hbar = 0.5*(hP+hQ);
    hLog = abs(log(hP./hQ));
    cMis = sqrt(sum((cPloc-cQmir).^2,2))./max(hbar,eps);
    bMinP = min(bP,[],2);
    bMinQ = min(bQ,[],2);
    bMis = abs(bMinP-bMinQ);
    qMis = abs(qP-qQ);

    M = struct();
    M.rI = rI;
    M.theta = th;
    M.thetaDeg = rad2deg(th);
    M.tip = tip.';
    M.e1 = e1.';
    M.e2 = e2.';
    M.elemP = eP;
    M.elemQ = eQ;
    M.baryP = bP;
    M.baryQ = bQ;
    M.hP = hP;
    M.hQ = hQ;
    M.qualityP = qP;
    M.qualityQ = qQ;
    M.centroidP_local = cPloc;
    M.centroidQ_local = cQloc;

    M.h_log_mismatch = hLog;
    M.centroid_mismatch = cMis;
    M.bary_margin_mismatch = bMis;
    M.quality_mismatch = qMis;

    M.A_h_median = median(hLog);
    M.A_h_rms = sqrt(mean(hLog.^2));
    M.A_h_max = max(hLog);
    M.A_centroid_median = median(cMis);
    M.A_bary_median = median(bMis);
    M.A_quality_median = median(qMis);

    etol = ip.Results.edge_tol;
    M.edge_tol = etol;
    M.frac_near_edge_P = mean(bMinP<etol);
    M.frac_near_edge_Q = mean(bMinQ<etol);
    M.frac_size_mismatch_gt_10pct = mean(hLog>log(1.10));
    M.frac_size_mismatch_gt_25pct = mean(hLog>log(1.25));

    M.PairTable = table(M.thetaDeg,eP,eQ,hP,hQ,hLog,bMinP,bMinQ,bMis, ...
        qP,qQ,qMis,cMis, ...
        'VariableNames',{'thetaDeg','elemP','elemQ','hP','hQ', ...
        'hLogMismatch','baryMinP','baryMinQ','baryMarginMismatch', ...
        'qualityP','qualityQ','qualityMismatch','centroidMismatch'});
end


function [h,q,c] = tri_metrics(X)
    e12 = norm(X(2,:)-X(1,:));
    e23 = norm(X(3,:)-X(2,:));
    e31 = norm(X(1,:)-X(3,:));
    A = 0.5*abs(det([X(2,:)-X(1,:);X(3,:)-X(1,:)]));
    h = sqrt(2*A);
    den = e12^2+e23^2+e31^2;
    if den<=eps
        q=0;
    else
        q=4*sqrt(3)*A/den;
    end
    c=mean(X,1);
end


function must_have(S,f)
    if ~isstruct(S) || ~isfield(S,f) || isempty(S.(f))
        error('SIF_mesh_mirror_metrics:MissingField', ...
            'Required field "%s" missing.',f);
    end
end
