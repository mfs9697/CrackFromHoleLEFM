function C=cfg_centered_half_domain()
%CFG_CENTERED_HALF_DOMAIN  Symmetry-reduced benchmark for the centered hole.
%
% Right-half model of the full centered-hole plate:
%   x in [A/2,A], y in [-B,B],
% with the right semicircle of the circular hole represented as a
% traction-free notch on the left symmetry boundary.
%
% This is a verification configuration.  It removes the physically
% equivalent left-hand initiation site without imposing horizontal
% symmetry about y=0, so Stage-II up/down kinking remains free.

C=cfg_hole_initiation();

xc=C.hole.center(1);

C.domain=struct();
C.domain.mode='centered_right_half';
C.domain.xmin=xc;
C.domain.xmax=C.A;
C.domain.symmetry_x=xc;

% ux=0 along the vertical symmetry boundary; one uy gauge constraint is
% added by solve_hole_only to eliminate rigid vertical translation.
C.bc.anchor_mode='symmetry_half_x';

% Keep the same angular resolution as the full 1440-point circle:
% 180 deg / 0.25 deg = 720 intervals -> 721 samples including endpoints.
C.stage1.nphi=721;
C.stage1.phi_range=[-pi/2,pi/2];
C.stage1.periodic=false;

% Candidate mesh-scaled angular regression from Step 16.  It is enabled
% here for the centered half-domain verification only; the general full
% production configuration remains unchanged until asymmetric validation.
C.stage1.angular_fit_halfwidth_factor=3.0;
end
