function Out=main_step4a_preview_graded_ring_mesh()
%MAIN_STEP4A_PREVIEW_GRADED_RING_MESH
% Geometry-only preview/gate for the corrected graded-ring S0 mesh.
% No SIF extraction is performed.

here=fileparts(mfilename('fullpath'));
addpath(genpath(fileparts(fileparts(here))));

[mesh,info]=build_graded_ring_crack_mesh( ...
    'r0',0.005,'r1',0.20,'Ntheta',64, ...
    'Variant','mirror_reflected');

fprintf('\n============================================================\n');
fprintf('STEP 4A: CORRECTED GRADED-RING S0 MESH PREVIEW\n');
fprintf('============================================================\n');
fprintf('Nr = %d, nominal Ntheta = %d\n',info.Nr,info.Ntheta);
fprintf('q_target = %.8f, q_actual = %.8f\n', ...
    info.qTargetNearEquilateral,info.qActual);
fprintf('T3 vertices/elements = %d / %d\n', ...
    info.nT3Vertices,info.nT3Elements);
fprintf('T6 nodes/elements = %d / %d\n', ...
    info.nT6Nodes,info.nT6Elements);
fprintf('upper ring segment spread = %.3e\n', ...
    info.upperRingSegmentRelSpreadMax);
fprintf('lower ring segment spread = %.3e\n', ...
    info.lowerRingSegmentRelSpreadMax);
fprintf('quality min/p05/median = %.5f / %.5f / %.5f\n', ...
    info.qualityMin,info.qualityP05,info.qualityMedian);
fprintf('quality at crack seam min = %.5f\n',info.qualityCrackMin);
fprintf('triangles with quality < 0.8 = %d\n',info.nQualityBelow08);
fprintf('max crack-face |x2| = %.3e\n',info.maxCrackFaceAbsY);
fprintf('crack faces distinct = %d\n',info.crackFacesDistinct);
fprintf('max mirror-coordinate error = %.3e\n',info.maxMirrorCoordError);

figure('Name','Step 4A corrected S0 full mesh','NumberTitle','off');
triplot(mesh.connect3,mesh.coord3(:,1),mesh.coord3(:,2));
axis equal; grid on; xlabel('x_1'); ylabel('x_2');
title('Corrected graded-ring S0 mesh');

figure('Name','Step 4A corrected S0 crack-seam zoom','NumberTitle','off');
triplot(mesh.connect3,mesh.coord3(:,1),mesh.coord3(:,2));
axis equal; grid on; xlabel('x_1'); ylabel('x_2');
xlim([-0.08,-0.005]);
ylim([-0.015,0.015]);
title('Corrected S0: negative-x crack-seam zoom');

Out=struct('mesh',mesh,'info',info);
end
