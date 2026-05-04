function stage1_plot_boundary_stress(B, Iauto, Iforced)
%STAGE1_PLOT_BOUNDARY_STRESS
% Plot tangential stress along the hole boundary.

    figure('Name', 'Stage I: hole-boundary tangential stress', 'Color', 'w');
    clf;
    hold on;
    box on;
    grid on;

    plot(rad2deg(B.phi), B.sig_tt_eff, 'LineWidth', 1.2, ...
        'DisplayName', '\sigma_{\theta\theta}^{eff}');

    plot(rad2deg(Iauto.phi_star), Iauto.sig_tt_unit, ...
        'ro', 'MarkerSize', 7, 'LineWidth', 1.5, ...
        'DisplayName', 'automatic maximum');

    plot(rad2deg(Iforced.phi_star), Iforced.sig_tt_unit, ...
        'ks', 'MarkerSize', 7, 'LineWidth', 1.5, ...
        'DisplayName', 'forced validation point');

    xlabel('\phi, deg');
    ylabel('\sigma_{\theta\theta}^{eff} at unit load');
    title('Stage I: tangential stress along circular hole');
    xlim([0, 360]);
    legend('Location', 'best');
end