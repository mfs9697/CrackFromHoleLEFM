function I = stage1_force_initiation_phi(C, B, phi0)
%STAGE1_FORCE_INITIATION_PHI
% Build an initiation struct at a prescribed hole angle.
%
% This is intended for validation benchmarks where symmetry should determine
% the initiation point exactly.

    hole = B.hole;
    c = hole.center(:).';
    R = hole.r;

    phi0 = mod(phi0, 2*pi);

    n_mat = [cos(phi0), sin(phi0)];
    t_hat = [-sin(phi0), cos(phi0)];
    x_star = c + R*n_mat;

    dphi = angle(exp(1i*(B.phi(:) - phi0)));
    [~, idx] = min(abs(dphi));

    sig_tt_unit = B.sig_tt_eff(idx);
    sig_tt_pos_unit = max(sig_tt_unit, 0);

    sigmax = max(max(B.sig_tt_eff, 0));

    sig0 = 1.0;
    if isfield(C, 'load') && isfield(C.load, 'sig0') && ~isempty(C.load.sig0)
        sig0 = C.load.sig0;
    end

    I = struct();

    I.idx_star = idx;
    I.phi_star = phi0;
    I.x_star = x_star;

    I.n_mat_star = n_mat;
    I.n_hole_star = -n_mat;
    I.t_hat_star = t_hat;

    I.sig_tt_unit = sig_tt_unit;
    I.sig_tt_pos_unit = sig_tt_pos_unit;

    I.lambda_ini = C.sig_c / sigmax;
    I.sig_applied_ini = I.lambda_ini * sig0;

    I.all_max_idx = idx;
    I.selection_rule = 'forced_phi';
end