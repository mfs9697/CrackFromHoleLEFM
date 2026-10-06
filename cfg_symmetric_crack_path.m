function C = cfg_symmetric_crack_path()
%CFG_SYMMETRIC_CRACK_PATH Full-domain centered-hole benchmark.
%
% Exact external symmetry:
%   plate  [0,A] x [-B,B]
%   hole   center [A/2,0], radius R
%   load   remote vertical tension
%
% Crack-path benchmark prescription:
%   P0 = [A/2+R,0]  (rightmost hole point)
%   n_mat = [1,0]
%   t_hat = [0,1]
%   Delta a = 4 mm
%
% The configuration deliberately inherits the already-qualified material,
% mesh scales, loading, anchoring and cohesive-carrier width from the
% asymmetric crack-path line. Only the hole position is changed.

    C=cfg_first_segment_asymmetric();

    C.hole.center=[0.5*C.A,0.0];
    C.holes={C.hole};

    % Explicit benchmark metadata.
    C.symmetricBenchmark=struct();
    C.symmetricBenchmark.axisY=0.0;
    C.symmetricBenchmark.side='right';
    C.symmetricBenchmark.P0=[C.hole.center(1)+C.hole.r,0.0];
    C.symmetricBenchmark.nMat=[1.0,0.0];
    C.symmetricBenchmark.tHat=[0.0,1.0];
    C.symmetricBenchmark.increment=C.a0;

    local_assert_fingerprint(C);
end


function local_assert_fingerprint(C)
    assert(abs(C.A-0.30)<=1e-14 && abs(C.B-0.10)<=1e-14, ...
        'symcfg:Plate','Unexpected plate dimensions.');
    assert(norm(C.hole.center-[0.15,0])<=1e-14 && ...
        abs(C.hole.r-0.030)<=1e-14, ...
        'symcfg:Hole','Unexpected centered-hole geometry.');
    assert(C.hole.npoly==480,'symcfg:Npoly','Expected Npoly=480.');
    assert(abs(C.a0-0.004)<=1e-14,'symcfg:Increment','Expected Delta a=4 mm.');
    assert(strcmp(C.load.type,'remote_tension_y') && ...
        abs(C.load.sig0-1)<=1e-14, ...
        'symcfg:Load','Expected unit remote-y traction.');
    assert(strcmp(C.bc.anchor_mode,'minimal'), ...
        'symcfg:Anchoring','Expected minimal anchoring.');
end
