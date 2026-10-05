function plot_s0_full_mesh(mesh,audit,Nr,Nth,showT6)
if nargin<5, showT6=true; end
figure;
triplot(mesh.connect3,mesh.coord3(:,1),mesh.coord3(:,2));
hold on;
if showT6
    plot(mesh.coord(:,1),mesh.coord(:,2),'.','MarkerSize',2);
end
axis equal; grid on;
xlabel('x_1'); ylabel('x_2');
title(sprintf('S0 Nr=%d Ntheta=%d T3=%d T6nodes=%d err=%.1e',Nr,Nth,size(mesh.connect3,1),size(mesh.coord,1),audit.maxMirrorCoordError));
hold off;
end
