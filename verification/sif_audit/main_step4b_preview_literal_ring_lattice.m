function Out=main_step4b_preview_literal_ring_lattice()
% Geometry-only preview of the literal first-picture ring lattice.

here=fileparts(mfilename('fullpath'));
addpath(genpath(fileparts(fileparts(here))));

[mesh,info]=build_literal_ring_lattice('r0',0.005,'r1',0.20,'Ntheta',64);

fprintf('\n============================================================\n');
fprintf('STEP 4B: LITERAL FIRST-PICTURE RING LATTICE PREVIEW\n');
fprintf('============================================================\n');
fprintf('Nr=%d, Ntheta=%d\n',info.Nr,info.Ntheta);
fprintf('q_eq=%.8f, q_actual=%.8f\n',info.qEquilateral,info.qActual);
fprintf('ring segment relative spread = %.3e\n',info.ringSegmentRelSpreadMax);
fprintf('quality min/p05/median = %.5f / %.5f / %.5f\n', ...
    info.qualityMin,info.qualityP05,info.qualityMedian);
fprintf('T3 vertices/elements = %d / %d\n', ...
    info.nT3Vertices,info.nT3Elements);

figure('Name','Literal first-picture ring lattice','NumberTitle','off');
triplot(mesh.connect3,mesh.coord3(:,1),mesh.coord3(:,2));
axis equal; grid on; xlabel('x_1'); ylabel('x_2');
title('Literal first-picture ring lattice before crack cut');

Out=struct('mesh',mesh,'info',info);
end
