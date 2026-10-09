function tf=crack_exterior_controls_match(a,b,coreA,coreB)
% Labels and reference flags are descriptive; requested numeric controls decide reuse.
tf=false;
try
    % A partially recorded family is not a metadata-free historical default.
    for controls={a,b}
        c=controls{1};
        if ~isstruct(c)||~isscalar(c)|| ...
                ~(isfield(c,'farCapOverIncrement')||isfield(c,'farCapOverA0'))|| ...
                ~(isfield(c,'transitionOverIncrement')||isfield(c,'transitionOverA0'))
            return
        end
    end
    a=crack_exterior_identity(a,coreA);b=crack_exterior_identity(b,coreB);
    for name={'exteriorScale','farCapOverIncrement','transitionOverIncrement'}
        x=a.(name{1});y=b.(name{1});
        if ~isnumeric(x)||~isnumeric(y)||~isreal(x)||~isreal(y)|| ...
                ~isscalar(x)||~isscalar(y)||~isfinite(x)||~isfinite(y)|| ...
                abs(x-y)>1e-14
            return
        end
    end
    tf=isequaln(a.calibrationOverride,b.calibrationOverride);
catch
    % Missing/malformed identity never permits reuse.
end
end
