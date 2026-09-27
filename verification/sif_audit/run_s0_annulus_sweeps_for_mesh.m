function R = run_s0_annulus_sweeps_for_mesh(mesh,U,mat,weightFunction)
% Run Step-3C annulus sweeps on one fixed S0 mesh and one mixed field.
if nargin<4 || isempty(weightFunction), weightFunction='analytic_radial'; end
V=[-1,0;0,0];

outer=[0.06 0.08 0.10 0.12 0.16];
fout=0.20;
A=nan(numel(outer),6);
for i=1:numel(outer)
    ro=outer(i); ri=fout*ro;
    d=struct('r_inner',ri,'r_outer',ro);
    [KI,KII]=SIF_LEFM_interaction_EDI(mesh,U,V,mat,d, ...
        'UsePlaneStrain',true,'Verbose',false,'WeightFunction',weightFunction);
    A(i,:)=[ro,ri,KI,KII,KI,KII/0.35];
end
R.outer=array2table(A,'VariableNames', ...
    {'r_outer','r_inner','KI','KII','KI_ratio','KII_ratio'});

ro=0.12;
fin=[0.10 0.20 0.30 0.40 0.50];
B=nan(numel(fin),6);
for i=1:numel(fin)
    f=fin(i); ri=f*ro;
    d=struct('r_inner',ri,'r_outer',ro);
    [KI,KII]=SIF_LEFM_interaction_EDI(mesh,U,V,mat,d, ...
        'UsePlaneStrain',true,'Verbose',false,'WeightFunction',weightFunction);
    B(i,:)=[f,ri,KI,KII,KI,KII/0.35];
end
R.inner=array2table(B,'VariableNames', ...
    {'inner_over_outer','r_inner','KI','KII','KI_ratio','KII_ratio'});

R.range=[ ...
    max(R.outer.KI_ratio)-min(R.outer.KI_ratio), ...
    max(R.outer.KII_ratio)-min(R.outer.KII_ratio), ...
    max(R.inner.KI_ratio)-min(R.inner.KI_ratio), ...
    max(R.inner.KII_ratio)-min(R.inner.KII_ratio)];
R.weightFunction=weightFunction;
end
