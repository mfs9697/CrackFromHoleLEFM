function R=run_stage1_hole_initiation(C)
%RUN_STAGE1_HOLE_INITIATION  Production Stage-I orchestration.
%
% Centralizes the hole-only solve, boundary-stress postprocessing, and
% initiation-point selection so production drivers cannot silently diverge.
%
% C.stage1.method:
%   'boundary_extrapolated_t6'  redesigned production method
%   'legacy_scattered'          historical method retained for reproduction

if ~isfield(C,'stage1')||isempty(C.stage1)
    error('run_stage1_hole_initiation:MissingStage1', ...
        'C.stage1 configuration is required.');
end

method=lower(strtrim(getf(C.stage1,'method','boundary_extrapolated_t6')));

G=geom_hole_only(C);
S1=solve_hole_only(C,G,'lambda',1.0);

switch method
    case 'boundary_extrapolated_t6'
        B=sample_hole_boundary_stress_v2(C,G,S1);
        I=find_hole_initiation_point_v2(C,B);

    case 'legacy_scattered'
        B=sample_hole_boundary_stress(C,G,S1);
        I=find_hole_initiation_point(C,B);

    otherwise
        error('run_stage1_hole_initiation:UnknownMethod', ...
            'Unknown C.stage1.method="%s".',method);
end

R=struct();
R.C=C;
R.G=G;
R.S1=S1;
R.B=B;
R.I=I;
R.method=method;
end


function v=getf(S,f,d)
if isstruct(S)&&isfield(S,f)&&~isempty(S.(f))
    v=S.(f);
else
    v=d;
end
end
