C  = cfg_hole_initiation();
G  = geom_hole_only(C);
S1 = solve_hole_only(C, G, 'lambda', 1.0);
B  = sample_hole_boundary_stress(C, G, S1);
I  = find_hole_initiation_point(C, B);

%{
figure; plot(B.phi*180/pi, B.sig_tt_eff, 'LineWidth', 1.2); hold on
plot(I.phi_star*180/pi, I.sig_tt_unit, 'ro', 'MarkerSize', 8, 'LineWidth', 1.5)
xlabel('\phi, deg'); ylabel('\sigma_{tt}');
grid on
%}

theta = 0;
G2 = geom_hole_shortcrack(C, I, theta);

Pmid