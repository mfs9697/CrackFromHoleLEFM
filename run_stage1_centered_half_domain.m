function R=run_stage1_centered_half_domain(C)
%RUN_STAGE1_CENTERED_HALF_DOMAIN  Stage-I symmetry-reduced benchmark.
%
% Uses the right-half plate with a semicircular cutout, ux=0 on the
% vertical symmetry boundary, recovered-T6 material-side stress sampling,
% eps->0 extrapolation, and the configured angular peak refinement.

if nargin<1 || isempty(C)
    C=cfg_centered_half_domain();
end

if ~isfield(C,'domain') || ~isfield(C.domain,'mode') || ...
        ~strcmpi(C.domain.mode,'centered_right_half')
    error('run_stage1_centered_half_domain:WrongConfig', ...
        'C.domain.mode must be "centered_right_half".');
end

G=geom_centered_half_hole(C);
S1=solve_hole_only(C,G,'lambda',1.0);
B=sample_hole_boundary_stress_v2(C,G,S1);
I=find_hole_initiation_point_v2(C,B);

R=struct();
R.C=C;
R.G=G;
R.S1=S1;
R.B=B;
R.I=I;
R.method='centered_right_half_boundary_extrapolated_t6';
end
