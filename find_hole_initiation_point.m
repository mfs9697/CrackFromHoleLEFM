function I = find_hole_initiation_point(C, B)
%FIND_HOLE_INITIATION_POINT  Determine initiation point and initiation load.
%
%   I = find_hole_initiation_point(C, B)
%
% Input:
%   C   configuration struct from cfg_hole_initiation
%   B   boundary-stress struct from sample_hole_boundary_stress
%
% Output:
%   I   initiation struct with fields:
%       .idx_star         index of the selected sampling point
%       .phi_star         angular position of initiation point
%       .x_star           initiation point coordinates [1 x 2]
%       .n_out_star       outward normal at initiation point [1 x 2]
%       .t_hat_star       tangent at initiation point [1 x 2]
%       .sig_tt_unit      effective tangential stress at unit load
%       .sig_tt_pos_unit  positive part of effective tangential stress
%       .lambda_ini       initiation load factor
%       .sig_applied_ini  applied nominal stress/load at initiation
%       .all_max_idx      indices whose values are tied (within tolerance)
%       .selection_rule   text label of the rule used
%
% Notes:
%   The initiation criterion is
%
%       lambda_ini = sig_c / max(sig_tt_eff^+)
%
%   where sig_tt_eff^+ = max(sig_tt_eff, 0).
%
%   The selected initiation point is the first maximizer by default.
%   If several points are tied within tolerance, the first one is used
%   and all such indices are returned in I.all_max_idx.

    % ------------------------------------------------------------
    % checks
    % ------------------------------------------------------------
    must(C, 'sig_c');
    must(B, 'phi');
    must(B, 'x');
    must(B, 'n_out');
    must(B, 't_hat');
    must(B, 'sig_tt_eff');

    sig_c = C.sig_c;

    phi   = B.phi(:);
    x     = B.x;
    n_out = B.n_out;
    t_hat = B.t_hat;
    sig   = B.sig_tt_eff(:);

    if numel(phi) ~= size(x,1) || numel(phi) ~= size(n_out,1) || numel(phi) ~= size(t_hat,1)
        error('find_hole_initiation_point:SizeMismatch', ...
            'Inconsistent sizes in boundary-stress struct B.');
    end

    % ------------------------------------------------------------
    % tensile part of tangential stress at unit load
    % ------------------------------------------------------------
    sig_pos = max(sig, 0);

    sigmax = max(sig_pos);

    if ~(isfinite(sigmax) && sigmax > 0)
        error('find_hole_initiation_point:NoTension', ...
            ['The maximum positive tangential stress is non-positive or non-finite. ', ...
             'The initiation criterion cannot be applied.']);
    end

    % ------------------------------------------------------------
    % identify maximizer(s)
    % ------------------------------------------------------------
    % tolerance for ties
    tol = 1e-10 * max(1, abs(sigmax));

    all_max_idx = find(abs(sig_pos - sigmax) <= tol);

    idx_star = all_max_idx(1);   % default choice: first one found
    selection_rule = 'first_maximizer';

    % ------------------------------------------------------------
    % local frame and initiation load
    % ------------------------------------------------------------
    x_star     = x(idx_star, :);

    % For a circular hole, the sampled radial normal points from the hole
    % center into the solid. That is the physically admissible crack-normal
    % direction at initiation.
    n_mat_star = n_out(idx_star, :);   % material-side normal
    n_hole_star = -n_mat_star;         % points into the hole cavity

    t_hat_star = t_hat(idx_star, :);
    phi_star   = phi(idx_star);

    lambda_ini = sig_c / sigmax;

    sig0 = 1.0;
    if isfield(C, 'load') && isfield(C.load, 'sig0') && ~isempty(C.load.sig0)
        sig0 = C.load.sig0;
    end
    sig_applied_ini = lambda_ini * sig0;

    % ------------------------------------------------------------
    % optional console message
    % ------------------------------------------------------------
    verbose = 0;
    if isfield(C, 'solver') && isfield(C.solver, 'verbose') && ~isempty(C.solver.verbose)
        verbose = C.solver.verbose;
    end

    if verbose
        fprintf(['find_hole_initiation_point: idx=%d, phi=%.6f rad (%.3f deg), ', ...
                 'sig_tt^+_max(unit)=%.6e, lambda_ini=%.6e, applied=%.6e\n'], ...
                idx_star, phi_star, phi_star*180/pi, sigmax, lambda_ini, sig_applied_ini);

        if numel(all_max_idx) > 1
            fprintf(['  note: %d tied maximizers detected within tolerance; ', ...
                     'using the first one.\n'], numel(all_max_idx));
        end
    end

    % ------------------------------------------------------------
    % output
    % ------------------------------------------------------------
    I = struct();

    I.idx_star         = idx_star;
    I.phi_star         = phi_star;
    I.x_star           = x_star;

    I.n_mat_star   = n_mat_star;
    I.n_hole_star  = n_hole_star;
    I.t_hat_star   = t_hat_star;

    I.sig_tt_unit      = sig(idx_star);
    I.sig_tt_pos_unit  = sig_pos(idx_star);

    I.lambda_ini       = lambda_ini;
    I.sig_applied_ini  = sig_applied_ini;

    I.all_max_idx      = all_max_idx;
    I.selection_rule   = selection_rule;
end


% =========================================================================
% helpers
% =========================================================================

function must(S, field)
%MUST Error if field does not exist or is empty.

    if ~isfield(S, field) || isempty(S.(field))
        error('find_hole_initiation_point:MissingField', ...
            'Required field "%s" is missing or empty.', field);
    end
end