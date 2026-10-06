function R0 = make_symmetric_crack_path_state()
%MAKE_SYMMETRIC_CRACK_PATH_STATE Build the prescribed symmetry benchmark state.
%
% This is NOT a Stage-I numerical solution. The initiation point and local
% frame are fixed exactly by symmetry at the rightmost point of the centered
% circular hole. The struct intentionally follows the FrozenState interface
% used by the qualified incremental crack machinery.
%
% No FEM solve is performed.

    C=cfg_symmetric_crack_path();

    P0=C.symmetricBenchmark.P0;
    nMat=C.symmetricBenchmark.nMat;
    tHat=C.symmetricBenchmark.tHat;

    Summary=table( ...
        C.a0,P0(1),P0(2),nMat(1),nMat(2),tHat(1),tHat(2), ...
        true,C.hole.npoly, ...
        'VariableNames',{ ...
        'a0_reserved_m','x_star_m','y_star_m', ...
        'nmat_x','nmat_y','that_x','that_y', ...
        'stage1_pass','hole_npoly'});

    gates=struct();
    gates.centeredHole=norm(C.hole.center-[0.5*C.A,0])<=1e-14;
    gates.rightmostInitiation=norm(P0-[C.hole.center(1)+C.hole.r,0])<=1e-14;
    gates.exactFrame=isequal(nMat,[1 0]) && isequal(tHat,[0 1]);
    gates.fourMillimeterIncrement=abs(C.a0-0.004)<=1e-14;
    gates.symmetricExternalGeometry=abs(C.hole.center(2))<=1e-14;
    gates.unitRemoteYTension=strcmp(C.load.type,'remote_tension_y') && ...
        abs(C.load.sig0-1)<=1e-14;

    assert(all(structfun(@logical,gates)), ...
        'symstate:Fingerprint','Symmetric benchmark fingerprint failed.');

    R0=struct();
    R0.C=C;
    R0.summary=Summary;
    R0.gates=gates;
    R0.stage1Pass=true;
    R0.isPrescribedSymmetryBenchmark=true;
    R0.method='symmetry_prescribed_centered_hole_benchmark';
    R0.interpretation=[ ...
        'Centered-hole full-domain benchmark. Initiation point and local ', ...
        'frame are prescribed exactly by symmetry; no Stage-I solve was ', ...
        'used to define them.'];
end
