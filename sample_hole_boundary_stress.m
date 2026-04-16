function B = sample_hole_boundary_stress(C, G, S1)
%SAMPLE_HOLE_BOUNDARY_STRESS  Sample tangential stress along a circular hole boundary.
%
%   B = sample_hole_boundary_stress(C, G, S1)
%
% Input:
%   C   configuration struct
%   G   geometry struct from geom_hole_only
%   S1  stage-I solution struct from solve_hole_only
%
% Output:
%   B   boundary-stress struct with fields:
%       .phi         [nphi x 1] sampling angles
%       .x           [nphi x 2] sampled boundary points on the exact circle
%       .xq          [nphi x 2] shifted query points inside the body
%       .n_out       [nphi x 2] outward normals
%       .t_hat       [nphi x 2] tangents (CCW orientation)
%       .sig_tt      [nphi x 1] tangential stress
%       .sig_nn      [nphi x 1] normal stress
%       .sig_nt      [nphi x 1] shear stress
%       .sig_tt_eff  [nphi x 1] effective tangential stress used in criterion
%       .hole        hole spec used
%
% Notes:
%   - This revised version uses scattered interpolation of recovered nodal
%     stresses instead of element-by-element point location.
%   - This is much more robust for Stage I and is good enough for detecting
%     the boundary stress maximum.
%   - A mesh-sized inward shift is used so the query points lie safely in
%     the body.

    % ------------------------------------------------------------
    % checks
    % ------------------------------------------------------------
    must(C,  'stage1');
    must(S1, 'mesh');
    must(S1, 'stress');

    hole = get_single_circular_hole(C, G);

    nphi     = getf(C.stage1, 'nphi', 720);
    avg_mode = lower(strtrim(getf(C.stage1, 'avg_mode', 'none')));
    avg_rho  = getf(C.stage1, 'avg_rho', 0.0);

    coord  = S1.mesh.coord;      % [nnod x 2]
    sNodal = S1.stress;          % [nnod x 3] = [sxx syy sxy]

    if size(sNodal,2) ~= 3
        error('sample_hole_boundary_stress:BadStress', ...
            'Expected S1.stress to have 3 columns: [sxx syy sxy].');
    end

    % ------------------------------------------------------------
    % angular sampling on the hole boundary
    % ------------------------------------------------------------
    c = hole.center(:).';
    R = hole.r;

    phi = linspace(0, 2*pi, nphi+1).';
    phi(end) = [];  % remove duplicate endpoint

    n_out = [cos(phi), sin(phi)];
    t_hat = [-sin(phi), cos(phi)];

    xb = c + R * n_out;   % exact boundary points [nphi x 2]

    % ------------------------------------------------------------
    % inward shift based on mesh size, not on radius
    % ------------------------------------------------------------
    hhole = [];
    if isfield(C, 'mesh1') && isfield(C.mesh1, 'hhole') && ~isempty(C.mesh1.hhole)
        hhole = C.mesh1.hhole;
    end
    hmax = [];
    if isfield(C, 'mesh1') && isfield(C.mesh1, 'hmax') && ~isempty(C.mesh1.hmax)
        hmax = C.mesh1.hmax;
    end

    if ~isempty(hhole)
        eps_shift = 0.25 * hhole;
    elseif ~isempty(hmax)
        eps_shift = 0.10 * hmax;
    else
        eps_shift = 1e-3 * R;
    end

    % cap to avoid over-shifting
    eps_shift = min(eps_shift, 0.10 * R);
    eps_shift = max(eps_shift, 1e-8 * max(R,1));

    xq = xb - eps_shift * n_out;

    % ------------------------------------------------------------
    % scattered interpolation of recovered nodal stresses
    % ------------------------------------------------------------
    Fx  = scatteredInterpolant(coord(:,1), coord(:,2), sNodal(:,1), 'linear', 'nearest');
    Fy  = scatteredInterpolant(coord(:,1), coord(:,2), sNodal(:,2), 'linear', 'nearest');
    Fxy = scatteredInterpolant(coord(:,1), coord(:,2), sNodal(:,3), 'linear', 'nearest');

    sig_xx = Fx(xq(:,1), xq(:,2));
    sig_yy = Fy(xq(:,1), xq(:,2));
    sig_xy = Fxy(xq(:,1), xq(:,2));

    % ------------------------------------------------------------
    % local stress components on the hole boundary
    % ------------------------------------------------------------
    sig_tt = nan(nphi,1);
    sig_nn = nan(nphi,1);
    sig_nt = nan(nphi,1);

    for k = 1:nphi
        S = [sig_xx(k), sig_xy(k);
             sig_xy(k), sig_yy(k)];

        n = n_out(k,:).';
        t = t_hat(k,:).';

        sig_nn(k) = n.' * S * n;
        sig_tt(k) = t.' * S * t;
        sig_nt(k) = n.' * S * t;
    end

    % ------------------------------------------------------------
    % optional averaging
    % ------------------------------------------------------------
    switch avg_mode
        case 'none'
            sig_tt_eff = sig_tt;

        case 'arc_average'
            if avg_rho <= 0
                warning('sample_hole_boundary_stress:BadAvgRho', ...
                    'avg_rho <= 0 with avg_mode=''arc_average''; using no averaging.');
                sig_tt_eff = sig_tt;
            else
                sig_tt_eff = circular_arc_average(phi, sig_tt, R, avg_rho);
            end

        otherwise
            error('sample_hole_boundary_stress:UnknownAvgMode', ...
                'Unsupported C.stage1.avg_mode = "%s".', avg_mode);
    end

    % ------------------------------------------------------------
    % output
    % ------------------------------------------------------------
    B = struct();
    B.phi        = phi;
    B.x          = xb;
    B.xq         = xq;
    B.n_out      = n_out;
    B.t_hat      = t_hat;

    B.sig_tt     = sig_tt;
    B.sig_nn     = sig_nn;
    B.sig_nt     = sig_nt;
    B.sig_tt_eff = sig_tt_eff;

    B.hole       = hole;
    B.eps_shift  = eps_shift;
end


% =========================================================================
% helpers
% =========================================================================

function hole = get_single_circular_hole(C, G)
%GET_SINGLE_CIRCULAR_HOLE Return the single circular hole spec.

    if isfield(G, 'hole') && ~isempty(G.hole)
        hole = G.hole;
    elseif isfield(C, 'hole') && ~isempty(C.hole)
        hole = C.hole;
    elseif isfield(C, 'holes') && numel(C.holes) == 1
        hole = C.holes{1};
    else
        error('sample_hole_boundary_stress:HoleSpec', ...
            'This first draft supports exactly one circular hole.');
    end

    if ~isstruct(hole) || ~isfield(hole, 'type')
        error('sample_hole_boundary_stress:BadHoleSpec', ...
            'Hole specification must be a struct with field "type".');
    end

    if ~strcmpi(strtrim(hole.type), 'circle')
        error('sample_hole_boundary_stress:UnsupportedHoleType', ...
            'Only circular holes are supported in this first draft.');
    end

    must(hole, 'center');
    must(hole, 'r');
end


function yavg = circular_arc_average(phi, y, R, rho)
%CIRCULAR_ARC_AVERAGE Circular moving average over arc half-length rho.

    n = numel(phi);
    dphi = rho / R;

    phi_ext = [phi - 2*pi; phi; phi + 2*pi];
    y_ext   = [y; y; y];

    yavg = zeros(n,1);

    for k = 1:n
        pk = phi(k);
        mask = (phi_ext >= pk - dphi) & (phi_ext <= pk + dphi);

        phw = phi_ext(mask);
        yw  = y_ext(mask);

        if numel(yw) < 2
            yavg(k) = y(k);
        else
            yavg(k) = trapz(phw, yw) / max(phw(end) - phw(1), eps);
        end
    end
end


function must(S, field)
%MUST Error if field does not exist or is empty.

    if ~isfield(S, field) || isempty(S.(field))
        error('sample_hole_boundary_stress:MissingField', ...
            'Required field "%s" is missing or empty.', field);
    end
end


function v = getf(S, field, default)
%GETF Get struct field or default.

    if isstruct(S) && isfield(S, field) && ~isempty(S.(field))
        v = S.(field);
    else
        v = default;
    end
end