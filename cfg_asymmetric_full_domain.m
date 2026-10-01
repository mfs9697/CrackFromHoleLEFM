function C=cfg_asymmetric_full_domain()
%CFG_ASYMMETRIC_FULL_DOMAIN  Full-domain asymmetric-hole benchmark.
%
% Purpose:
%   Production-style numerical example after the centered symmetry audit.
%   Both x- and y-reflection symmetries are broken so Stage I should select
%   one genuinely preferred initiation site.
%
% Geometry candidate:
%   plate x in [0,A], y in [-B,B]
%   circular hole center = [0.17,-0.02] m, R = 0.03 m
%
% This leaves minimum ligaments of 0.05 m to the bottom boundary and
% 0.10 m to the right boundary, so the case is asymmetric without being
% close to a near-contact singular geometry.

C=cfg_hole_initiation();

C.domain=struct();
C.domain.mode='asymmetric_full';

C.hole.center=[0.17,-0.02];
C.holes={C.hole};

% Recompute mesh scales explicitly after changing the benchmark geometry.
C.mesh1.hmin=2*pi*C.hole.r/C.hole.npoly;
C.mesh1.hhole=C.mesh1.hmin;
C.mesh1.hmax=20*C.mesh1.hmin;

C.mesh2.hmax=C.mesh1.hmax;
C.mesh2.hhole=C.mesh1.hmin;
C.mesh2.hcrack=C.mesh1.hmin;

% Full circular sampling.
C.stage1.nphi=1440;
C.stage1.periodic=true;

% Candidate from the centered-hole Step-16 audit.  Step 19 tests whether
% this mesh-scaled window preserves a genuine asymmetric peak.
C.stage1.angular_fit_halfwidth_factor=3.0;
end
