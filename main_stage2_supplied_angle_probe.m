function P = main_stage2_supplied_angle_probe(varargin)
%MAIN_STAGE2_SUPPLIED_ANGLE_PROBE
% Qualify and, only with explicit permission, solve ONE supplied Stage-II
% first-segment angle. This routine never performs an angle sweep.
%
% Recommended next experiment:
%
%   P = main_stage2_supplied_angle_probe( ...
%       'FrozenState',R0, ...
%       'Theta1Deg',-0.01, ...
%       'ReferenceResult',Rphys, ... % prior theta_1=0 result, optional
%       'AllowSolve',true);
%
% Sequence:
%   1) build the full-domain supplied-angle carrier;
%   2) embed the unchanged a0-scaled audited core;
%   3) pass all T3/T6/exterior/prescribed-Williams gates;
%   4) only then permit one physical SGS-PCG solve;
%   5) checkpoint before COD/EDI;
%   6) if ReferenceResult is supplied, report sign bracket/secant estimate.
%
% No subsequent angle is selected or solved automatically.

    ip=inputParser;
    addParameter(ip,'FrozenState',[],@(x)isempty(x)||isstruct(x));
    addParameter(ip,'Theta1Deg',[],@(x)isempty(x)||(isnumeric(x)&&isscalar(x)&&isfinite(x)));
    addParameter(ip,'ReferenceResult',[],@(x)isempty(x)||isstruct(x));
    addParameter(ip,'AllowSolve',false,@(x)islogical(x)&&isscalar(x));
    addParameter(ip,'NArc',480,@(x)isnumeric(x)&&isscalar(x)&&isfinite(x)&&x>=32&&x==round(x));
    addParameter(ip,'ExteriorVerbose',true,@(x)islogical(x)&&isscalar(x));
    addParameter(ip,'Plot',true,@(x)islogical(x)&&isscalar(x));
    parse(ip,varargin{:});
    opt=ip.Results;

    if isempty(opt.Theta1Deg)
        error('stage2probe:MissingTheta1Deg', ...
            'Supply one probe angle explicitly via ''Theta1Deg''.');
    end

    theta1Deg=double(opt.Theta1Deg);

    fprintf('\n============================================================\n');
    fprintf('STAGE II: ONE SUPPLIED-ANGLE PROBE\n');
    fprintf('============================================================\n');
    fprintf('  requested theta_1 = %+.12g deg\n',theta1Deg);
    fprintf('  physical solve    = %d\n',logical(opt.AllowSolve));
    fprintf('  angle sweep       = 0\n\n');

    Fangle=main_stage2_embed_scaled_core_full_domain( ...
        'FrozenState',opt.FrozenState, ...
        'Theta1Deg',theta1Deg, ...
        'NArc',opt.NArc, ...
        'ExteriorVerbose',opt.ExteriorVerbose, ...
        'SaveCandidate',true, ...
        'SaveCompact',true, ...
        'Plot',opt.Plot);

    if ~Fangle.pass || ...
            ~all(structfun(@logical,Fangle.gates)) || ...
            ~all(structfun(@logical,Fangle.syntheticGates))
        error('stage2probe:MeshQualificationFailed', ...
            'Supplied-angle mesh/synthetic qualification did not pass.');
    end

    fprintf('\nSUPPLIED-ANGLE MESH QUALIFIED.\n');
    fprintf('  Reaching physical solver gate.\n');

    Rangle=main_stage2_supplied_angle_physical_solve( ...
        'FrozenState',opt.FrozenState, ...
        'Candidate',Fangle.candidate, ...
        'Theta1Deg',theta1Deg, ...
        'AllowSolve',opt.AllowSolve);

    if ~Rangle.pass
        error('stage2probe:PhysicalQualificationFailed', ...
            'Supplied-angle physical solve/postprocessing did not pass.');
    end

    Bracket=local_bracket(opt.ReferenceResult,Rangle);

    P=struct();
    P.theta1Deg=theta1Deg;
    P.meshQualification=Fangle;
    P.physical=Rangle;
    P.bracket=Bracket;
    P.pass=Fangle.pass&&Rangle.pass;
    P.noAngleSweep=true;
    P.noAutomaticNextProbe=true;

    fprintf('\nSUPPLIED-ANGLE PROBE COMPLETE.\n');
    fprintf('  theta_1 = %+.12g deg\n',theta1Deg);
    fprintf('  KII/KI  = %+.12g\n',Rangle.EDI.KII_over_KI);
    if Bracket.referenceAvailable
        fprintf('  reference theta = %+.12g deg, KII/KI=%+.12g\n', ...
            Bracket.referenceThetaDeg,Bracket.referenceRatio);
        if Bracket.signBracket
            fprintf('  SIGN BRACKET ESTABLISHED: [%+.12g, %+.12g] deg\n', ...
                Bracket.lowerThetaDeg,Bracket.upperThetaDeg);
            fprintf('  secant theta(KII=0) = %+.12g deg\n',Bracket.secantThetaDeg);
        else
            fprintf('  No sign bracket with supplied reference.\n');
        end
    end
    fprintf('  No additional physical angle was solved.\n');
end


function B=local_bracket(ref,cur)
    B=struct( ...
        'referenceAvailable',false, ...
        'referenceThetaDeg',NaN, ...
        'referenceKII',NaN, ...
        'referenceRatio',NaN, ...
        'currentThetaDeg',cur.theta1Deg, ...
        'currentKII',cur.EDI.KII_unit, ...
        'currentRatio',cur.EDI.KII_over_KI, ...
        'signBracket',false, ...
        'lowerThetaDeg',NaN, ...
        'upperThetaDeg',NaN, ...
        'secantThetaDeg',NaN);

    if isempty(ref)
        return
    end

    if ~isfield(ref,'pass') || ~logical(ref.pass) || ...
            ~isfield(ref,'EDI') || ~istable(ref.EDI) || height(ref.EDI)~=1
        error('stage2probe:BadReference', ...
            'ReferenceResult must be one passed physical Stage-II result.');
    end

    if isfield(ref,'theta1Deg')
        th0=ref.theta1Deg;
    elseif isfield(ref,'summary') && ...
            ismember('theta1_deg',ref.summary.Properties.VariableNames)
        th0=ref.summary.theta1_deg;
    else
        % The accepted legacy Stage-II-D result is the known theta=0 case.
        th0=0;
    end

    if ~ismember('KII_unit',ref.EDI.Properties.VariableNames) || ...
       ~ismember('KII_over_KI',ref.EDI.Properties.VariableNames)
        error('stage2probe:BadReferenceEDI', ...
            'ReferenceResult.EDI lacks KII_unit/KII_over_KI.');
    end

    f0=ref.EDI.KII_unit;
    q0=ref.EDI.KII_over_KI;
    f1=cur.EDI.KII_unit;
    th1=cur.theta1Deg;

    B.referenceAvailable=true;
    B.referenceThetaDeg=th0;
    B.referenceKII=f0;
    B.referenceRatio=q0;

    if isfinite(f0)&&isfinite(f1)&&f0*f1<=0 && th0~=th1
        B.signBracket=true;
        B.lowerThetaDeg=min(th0,th1);
        B.upperThetaDeg=max(th0,th1);

        den=f1-f0;
        if abs(den)>eps(max(abs([f0 f1])))
            B.secantThetaDeg=th0-f0*(th1-th0)/den;
        end
    end
end
