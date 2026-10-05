function [G2,D,M,Mc]=build_stage2_centered_half_cracked_mesh_for_theta(C,I,theta,varargin)
ip=inputParser;
addParameter(ip,'PlotGeom',false,@(x)islogical(x)||isnumeric(x));
addParameter(ip,'PlotMesh',false,@(x)islogical(x)||isnumeric(x));
addParameter(ip,'PlotCollapsed',false,@(x)islogical(x)||isnumeric(x));
addParameter(ip,'a0',C.a0,@(x)isnumeric(x)&&isscalar(x)&&x>0);
parse(ip,varargin{:});
O=ip.Results;

hole=C.hole;
c=hole.center(:).';
R=hole.r;
A0=[c(1)+R,c(2)];

if nargin<2 || isempty(I)
    I=struct();
end

Iexact=I;
Iexact.phi_star=0;
Iexact.x_star=A0;
Iexact.n_mat_star=[1,0];
Iexact.n_hole_star=[-1,0];
Iexact.t_hat_star=[0,1];

edir=[cos(theta),sin(theta)];
xtip=A0+O.a0*edir;
Pmid=[A0;xtip];

G2=struct();
G2.hole=hole;
G2.crack=struct('x0',A0,'xtip',xtip,'a0',O.a0, ...
    'theta',theta,'thetaDeg',rad2deg(theta), ...
    'e_dir',edir,'n_mat',[1,0],'t_hat',[0,1],'polyline',Pmid);
G2.tip=struct('x',xtip,'tangent',edir,'normal',[-edir(2),edir(1)], ...
    'radiusJ',C.mesh2.tip_radius);
G2.initiation=Iexact;
G2.meta=struct('domain_mode','centered_right_half');

D=build_domain_centered_half_pencil(Pmid,C,C.mesh2.chw);

if logical(O.PlotGeom)
    figure('Name','Centered half-domain Stage II geometry','Color','w'); clf
    hold on; axis equal; box on
    P=D.outerPoly;
    plot([P(:,1);P(1,1)],[P(:,2);P(1,2)],'k-','LineWidth',1.2);
    plot(Pmid(:,1),Pmid(:,2),'r-','LineWidth',1.8);
    plot(xtip(1),xtip(2),'kp','MarkerSize',10,'LineWidth',1.4);
    xlabel('x'); ylabel('y');
    title(sprintf('Centered half-domain Stage II, theta=%+.3f deg',rad2deg(theta)));
end

M=mesh_hole_pencil_domain(D, ...
    'Hmin',C.mesh1.hmin, ...
    'Hmax',C.mesh1.hmax, ...
    'Hgrad',C.mesh1.hgrad, ...
    'PlotGeom',logical(O.PlotGeom), ...
    'PlotMesh',logical(O.PlotMesh));

if ~isfield(M,'region') || ~isfield(M.region,'geomIDs') || isempty(M.region.geomIDs)
    error('build_stage2_centered_half_cracked_mesh_for_theta:NoGeomIDs', ...
        'Could not recover the two appendix-face geometry IDs.');
end

ids=M.region.geomIDs;
Mc=collapse_pencil_faces_to_midline(M,D, ...
    'EdgeIDs',ids.e_tip, ...
    'TipVertexID',ids.v_tip);

Mc.region.domain_mode='centered_right_half';

if logical(O.PlotCollapsed)
    plot_collapsed_pencil_mesh(Mc,'ShowOriginalFaces',true);
end
end
