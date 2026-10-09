function exteriorScale=resolve_crack_exterior_scale(coreScale,exteriorScale)
% Omitted/empty ExteriorScale keeps the historical CoreScale coupling.
if nargin<2||isempty(exteriorScale),exteriorScale=coreScale;end
assert(isnumeric(exteriorScale)&&isreal(exteriorScale)&& ...
    isscalar(exteriorScale)&&isfinite(exteriorScale)&&exteriorScale>0, ...
    'crackmesh:ExteriorScale','ExteriorScale must be a finite positive scalar.');
end
