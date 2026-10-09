function identity=crack_exterior_identity(controls,coreScale)
% Canonical exterior-only identity. Legacy records used their parent core scale.
if nargin<2,coreScale=1;end
es=[];
if isfield(controls,'exteriorScale'),es=controls.exteriorScale;end
identity=struct('exteriorScale',resolve_crack_exterior_scale(coreScale,es), ...
    'farCapOverIncrement',.625,'transitionOverIncrement',1, ...
    'calibrationOverride',struct());
for pair={'farCapOverIncrement','farCapOverA0'; ...
          'transitionOverIncrement','transitionOverA0'}.'
    key=pair{1};alias=pair{2};
    if isfield(controls,key),identity.(key)=controls.(key);
    elseif isfield(controls,alias),identity.(key)=controls.(alias);end
end
if isfield(controls,'calibrationOverride')
    identity.calibrationOverride=controls.calibrationOverride;
end
end
