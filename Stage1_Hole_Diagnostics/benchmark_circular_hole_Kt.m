function P = benchmark_circular_hole_Kt()
%BENCHMARK_CIRCULAR_HOLE_KT
% Stress concentration near a circular hole in a finite rectangular plate.

%% Geometry / material / load
A = 0.30;          % plate length, x in [0,A]
B = 0.10;          % half-height, y in [-B,B]
c = [0.5*A, 0];    % hole center
R = 0.030;         % hole radius

E = 210e3;
nu = 0.30;
sig0 = 1.0;        % remote nominal traction

nphi = 1440;
hmax = 2*pi*R/360 * 20;
hmin = 2*pi*R/360;

%% Structural PDE model
model = createpde('structural','static-planestrain');

% Rectangle minus exact circle using decsg
R1 = [3 4 0 A A 0 -B -B B B]';
C1 = [1 c(1) c(2) R 0 0 0 0 0 0]';

gd = [R1 C1];
ns = char('R1','C1')';
sf = 'R1-C1';

dl = decsg(gd, sf, ns);
geometryFromEdges(model, dl);

structuralProperties(model, ...
    'YoungsModulus', E, ...
    'PoissonsRatio', nu);

%% Boundary conditions and loading
geom = model.Geometry;

eTop = nearestEdge(geom, [0.5*A, B]);
eBot = nearestEdge(geom, [0.5*A, -B]);

structuralBoundaryLoad(model, 'Edge', eTop, ...
    'SurfaceTraction', [0; sig0]);

structuralBoundaryLoad(model, 'Edge', eBot, ...
    'SurfaceTraction', [0; -sig0]);

% Minimal rigid-body constraints
vLB = nearest_vertex(geom.Vertices, [0, -B]);
vRB = nearest_vertex(geom.Vertices, [A, -B]);

structuralBC(model, 'Vertex', vLB, ...
    'XDisplacement', 0, 'YDisplacement', 0);

structuralBC(model, 'Vertex', vRB, ...
    'YDisplacement', 0);

%% Mesh and solve
generateMesh(model, ...
    'Hmax', hmax, ...
    'Hmin', hmin, ...
    'Hgrad', 1.3, ...
    'GeometricOrder', 'quadratic');

figure('Color','w');
pdemesh(model);

Rsol = solve(model);

%% Sample tangential stress near hole contour
phi = linspace(0, 2*pi, nphi+1).';
phi(end) = [];

n = [cos(phi), sin(phi)];       % radial direction into material
t = [-sin(phi), cos(phi)];      % tangent direction

epsShift = 0.1*hmin;
xb = c + R*n;
xq = xb + epsShift*n;           % slightly inside material, outside hole

S = interpolateStress(Rsol, xq(:,1), xq(:,2));

sig_tt = zeros(nphi,1);
for k = 1:nphi
    Smat = [S.sxx(k), S.sxy(k);
            S.sxy(k), S.syy(k)];

    tk = t(k,:).';
    sig_tt(k) = tk.' * Smat * tk;
end

[sig_tt_peak, idx] = max(max(sig_tt, 0));
Kt = sig_tt_peak / sig0;

%% Output
P = struct();
P.model = model;
P.result = Rsol;
P.phi = phi;
P.phi_deg = rad2deg(phi);
P.x_boundary = xb;
P.x_query = xq;
P.sig_tt = sig_tt;
P.sig_tt_peak = sig_tt_peak;
P.Kt = Kt;
P.phi_peak = phi(idx);
P.phi_peak_deg = rad2deg(phi(idx));
P.epsShift = epsShift;

fprintf('Kt = %.8f\n', Kt);
fprintf('Peak phi = %.5f deg\n', P.phi_peak_deg);

%% Plot
figure('Color','w');
plot(P.phi_deg, sig_tt, 'LineWidth', 1.2); hold on
plot(P.phi_peak_deg, sig_tt_peak, 'ro', 'LineWidth', 1.5)
grid on; box on
xlabel('\phi, deg')
ylabel('\sigma_{\theta\theta}')
title(sprintf('Circular-hole stress concentration: Kt = %.4f', Kt))
xlim([0 360])
end


function idx = nearest_vertex(V, xq)
%NEAREST_VERTEX Return nearest PDE geometry vertex id.

if size(V,1) == 2
    P = V.';
else
    P = V;
end

[~, idx] = min(sum((P - xq).^2, 2));
end