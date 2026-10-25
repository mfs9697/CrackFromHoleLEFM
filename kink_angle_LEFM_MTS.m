function [theta_k, theta_deg, sigma_tt_val] = kink_angle_LEFM_MTS(KI, KII)
%KINK_ANGLE_LEFM_MTS Signed tensile maximum tangential stress direction.
% [theta_k,theta_deg,sigma_tt_val]=kink_angle_LEFM_MTS(KI,KII)
%
% theta_k is measured from the last crack tangent (radians).
% Uses true stationary roots of
%     KI*sin(theta)+KII*(3*cos(theta)-1) = 0
% and selects the largest *positive local tensile maximum* in the original
% [-80,80] degree forward search sector. The two stationary roots are not
% in general related by reflection theta -> -theta.
%
% No physically admissible tensile interior root -> returns NaN triplet
% (warns), rather than supplying a spurious propagation direction.
% In particular negative pure Mode I is NOT assigned theta=0.
% No FEM, SIF evaluation, mesh change or geometry solve performed here.
%
% This helper does NOT determine whether a material growth threshold
% has been satisfied; the caller must separately apply that condition.

if ~(isscalar(KI)&&isscalar(KII)&&isreal(KI)&&isreal(KII)&& ...
        isfinite(KI)&&isfinite(KII))
    error('kink_angle_LEFM_MTS:NonfiniteOrNonscalarSIF', ...
        'KI and KII must be finite, real scalar signed SIFs.');
end
theta_k=NaN;theta_deg=NaN;sigma_tt_val=NaN;
thetaBound=80*pi/180;
scale=max(1,hypot(KI,KII));
if abs(KII)<=1e-14*scale
    if KI>0
        theta_k=0;theta_deg=0;sigma_tt_val=KI;
    else
        warning('kink_angle_LEFM_MTS:NoTensileStationaryRoot', ...
            'No tensile forward MTS direction for pure compressive or zero Mode I.');
    end
    return;
end
D=hypot(KI,sqrt(8)*KII);
den=KI+D;
% Numerically stable branch for KI>=0 plus the complementary root.
% For extreme KI<0 with loss of significance in den, the finite
% complementary root remains available. Nonfinite candidates are ignored.
if den>0
    t1=-2*KII/den;
else
    t1=sign(-KII)*Inf;
end
t2=den/(4*KII);
candidates=2*atan([t1 t2]);
rootAccepted=false(size(candidates));
peak=-Inf(size(candidates));
for k=1:numel(candidates)
    theta=candidates(k);
    if ~isfinite(theta)||abs(theta)>thetaBound
        continue
    end
    c=cos(theta/2); s=sin(theta/2);
    stationaryResidual=KI*sin(theta)+KII*(3*cos(theta)-1);
    stress=KI*c^3-3*KII*s*c^2;
    % Exact curvature of leading sigma_tt at a stationary root:
    curvature=-.75*c*(KI*cos(theta)-3*KII*sin(theta));
    if abs(stationaryResidual)<=1e-11*scale&& ...
            stress>1e-12*scale&&curvature<0
        rootAccepted(k)=true;
        peak(k)=stress;
    end
end
if ~any(rootAccepted)
    warning('kink_angle_LEFM_MTS:NoTensileStationaryRoot', ...
        'No admissible tensile stationary maximum in +/-80 deg; MTS is undefined.');
    return
end
[~,idx]=max(peak);
theta_k=candidates(idx);
theta_deg=theta_k*180/pi;
sigma_tt_val=peak(idx);
end
