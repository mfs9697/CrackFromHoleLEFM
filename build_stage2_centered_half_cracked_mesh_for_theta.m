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
