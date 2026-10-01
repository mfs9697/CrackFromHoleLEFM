function Out=main_step17_centered_half_domain_verification(varargin)
%MAIN_STEP17_CENTERED_HALF_DOMAIN_VERIFICATION
% Verify the centered right-half Stage-I model against the full-domain
% centered-hole solution at matching nominal boundary resolution.
%
% The half-domain is the preferred symmetry benchmark because it removes
% the physically equivalent left-hand initiation site while retaining both
% y>0 and y<0 material, so later Stage-II up/down kinking is not prescribed.
%
% Defaults:
%   NpolyList = [240 480]
%   Plot = true
%
% Comparisons:
%   - half-domain boundary-limit sigma_tt(phi) on [-90,90] deg
%     versus the right half of the full-domain boundary-limit field;
%   - sigma_tt(0), local symmetry defect, and fitted right-peak angle;
%   - relative Linf/L2 curve differences;
%   - half-domain symmetry BC diagnostics.

ip=inputParser;
addParameter(ip,'NpolyList',[240 480], ...
    @(x)isnumeric(x)&&isvector(x)&&all(isfinite(x))&&all(x>=32));
addParameter(ip,'Plot',true,@(x)islogical(x)||isnumeric(x));
parse(ip,varargin{:});
O=ip.Results;

addpath(genpath(pwd));

fprintf('\n============================================================\n');
fprintf('SIF AUDIT STEP 17: CENTERED HALF-DOMAIN VERIFICATION\n');
fprintf('============================================================\n');

Nlist=round(O.NpolyList(:));
nL=numel(Nlist);
rows=nan(nL,25);
Half=cell(nL,1);
Full=cell(nL,1);

for il=1:nL
    Np=Nlist(il);

    %% Half-domain
    Ch=cfg_centered_half_domain();
    Ch.hole.npoly=Np;
    Ch.holes={Ch.hole};
    hArc=2*pi*Ch.hole.r/Np;
    Ch.mesh1.hmin=hArc;
    Ch.mesh1.hhole=hArc;
    Ch.mesh1.hmax=20*hArc;
    Ch.solver.verbose=0;
    Ch.plot.show_mesh1=false;

    Rh=run_stage1_centered_half_domain(Ch);
    Half{il}=Rh;

    %% Full-domain reference on an independently generated mesh
    Cf=cfg_hole_initiation();
    Cf.hole.npoly=Np;
    Cf.holes={Cf.hole};
    Cf.mesh1.hmin=hArc;
    Cf.mesh1.hhole=hArc;
    Cf.mesh1.hmax=20*hArc;
    Cf.solver.verbose=0;
    Cf.plot.show_mesh1=false;

    Gf=geom_hole_only(Cf);
    Sf=solve_hole_only(Cf,Gf,'lambda',1.0);
    Bf=sample_hole_boundary_stress_v2(Cf,Gf,Sf);

    % Right-hand full-domain local peak, using the same mesh-scaled
    % half-width factor as the half-domain verification configuration.
    factor=Ch.stage1.angular_fit_halfwidth_factor;
    hw=factor*hArc/Ch.hole.r;
    Pf=local_fit_right_peak(Bf,hw);

    Rf=struct('C',Cf,'G',Gf,'S1',Sf,'B',Bf,'rightPeak',Pf);
    Full{il}=Rf;

    %% Common right-semicircle grid
    phiH=Rh.B.phi(:);
    sigH=Rh.B.sig_tt_eff(:);

    phiF=local_wrap(Bf.phi(:));
    [phiFs,ord]=sort(phiF);
    sigFs=Bf.sig_tt_eff(ord);

    sigF=interp1(phiFs,sigFs,phiH,'linear');

    if any(~isfinite(sigF))
        error('step17:InterpolationFailure', ...
            'Could not map the full-domain right-half stress to the half-domain grid.');
    end

    scale=max([abs(sigH);abs(sigF);eps]);
    d=sigH-sigF;
    relLinf=max(abs(d))/scale;
    relL2=norm(d)/max(norm([sigH;sigF])/sqrt(2),eps);

    [~,i0h]=min(abs(phiH));
    [~,i0f]=min(abs(phiFs));

    Mh=local_symmetry(phiH,sigH);
    Mf=local_symmetry(phiH,sigF);

    % Half-domain traction residuals after radial boundary extrapolation.
    rightWin=abs(phiH)<=deg2rad(12);
    hscale=max(abs(sigH(rightWin)));
    nnH=max(abs(Rh.B.sig_nn(rightWin)))/max(hscale,eps);
    ntH=max(abs(Rh.B.sig_nt(rightWin)))/max(hscale,eps);

    % Symmetry-BC quality directly from T6 displacement solution.
    symNodes=Rh.S1.bc.symmetry_nodes;
    uxSym=max(abs(Rh.S1.U(2*symNodes-1)));

    halfPhi=rad2deg(Rh.I.phi_star);
    fullPhi=rad2deg(local_wrap(Pf.phiStar));

    rows(il,:)=[ ...
        Np,hArc,hArc/Ch.hole.r, ...
        size(Rh.G.p,1),size(Rh.G.t,1),size(Rh.S1.mesh.coord,1), ...
        size(Gf.p,1),size(Gf.t,1),size(Sf.mesh.coord,1), ...
        sigH(i0h),sigF(i0f), ...
        (sigH(i0h)-sigF(i0f))/max([abs(sigH(i0h)),abs(sigF(i0f)),eps]), ...
        relLinf,relL2, ...
        halfPhi,fullPhi,halfPhi-fullPhi, ...
        Mh.symRel,Mf.symRel, ...
        Rh.I.sig_tt_pos_unit,Pf.sigStar, ...
        Rh.I.lambda_ini,Ch.sig_c/Pf.sigStar, ...
        nnH,ntH];

    fprintf('\n--- Npoly=%d | h/R=%.8g rad = %.6f deg ---\n', ...
        Np,hArc/Ch.hole.r,rad2deg(hArc/Ch.hole.r));
    fprintf('  half mesh T3/T6 nodes = %d / %d | elements=%d\n', ...
        size(Rh.G.p,1),size(Rh.S1.mesh.coord,1),size(Rh.G.t,1));
    fprintf('  full mesh T3/T6 nodes = %d / %d | elements=%d\n', ...
        size(Gf.p,1),size(Sf.mesh.coord,1),size(Gf.t,1));
    fprintf('  max |ux| on symmetry boundary = %.6e\n',uxSym);
    fprintf('  sigma_tt(0): half=%.10e | full=%.10e | rel diff=% .3e\n', ...
        sigH(i0h),sigF(i0f),rows(il,12));
    fprintf('  right-semicircle curve: Linf=% .3e | L2=% .3e\n', ...
        relLinf,relL2);
    fprintf('  symmetry defect: half=% .3e | full-right=% .3e\n', ...
        Mh.symRel,Mf.symRel);
    fprintf('  fitted right peak: half=%+.8f deg | full=%+.8f deg | diff=%+.3e deg\n', ...
        halfPhi,fullPhi,halfPhi-fullPhi);
    fprintf('  peak stress: half=%.10e | full=%.10e\n', ...
        Rh.I.sig_tt_pos_unit,Pf.sigStar);
    fprintf('  lambda_ini: half=%.10e | full-right=%.10e\n', ...
        Rh.I.lambda_ini,Ch.sig_c/Pf.sigStar);
    fprintf('  half boundary residuals near right peak: nn=% .3e | nt=% .3e\n', ...
        nnH,ntH);
end

T=array2table(rows,'VariableNames',{ ...
    'Npoly','h_arc','h_arc_over_R', ...
    'half_T3_nodes','half_T3_elements','half_T6_nodes', ...
    'full_T3_nodes','full_T3_elements','full_T6_nodes', ...
    'half_sig0','full_sig0','sig0_rel_diff', ...
    'curve_rel_Linf','curve_rel_L2', ...
    'half_phi_fit_deg','full_right_phi_fit_deg','phi_fit_diff_deg', ...
    'half_sym_rel','full_right_sym_rel', ...
    'half_peak_sig','full_right_peak_sig', ...
    'half_lambda_ini','full_right_lambda_ini', ...
    'half_sig_nn_relmax','half_sig_nt_relmax'});

fprintf('\nCENTERED HALF-DOMAIN VERIFICATION TABLE\n');
disp(T);

if logical(O.Plot)
    local_plot_overlay(Half,Full,Nlist);
    local_plot_difference(Half,Full,Nlist);
    local_plot_convergence(T);
    local_plot_half_geometry(Half{1});
end

Out=struct();
Out.table=T;
Out.half=Half;
Out.full=Full;
Out.settings=O;

fprintf('\nSTEP 17 completed.\n');
fprintf(['Gate: half-domain and full-domain right-semicircle boundary stresses ', ...
    'should converge toward one another, while the half-domain fitted peak ', ...
    'should remain near phi=0 without left/right degeneracy.\n']);
end


function P=local_fit_right_peak(B,halfWidth)
phi=local_wrap(B.phi(:));
sig=max(B.sig_tt_eff(:),0);

right=abs(phi)<=deg2rad(12);
ids=find(right);
[~,j]=max(sig(right));
idx0=ids(j);
phi0=phi(idx0);

d=local_wrap(phi-phi0);
keep=abs(d)<=halfWidth+100*eps;
x=d(keep);
y=sig(keep);

[x,ord]=sort(x);
y=y(ord);

if numel(x)<5
    error('step17:TooFewFitPoints','Too few samples for full-domain right-peak fit.');
end

p=polyfit(x,y,2);
if p(1)>=0
    error('step17:BadCurvature','Full-domain right-peak quadratic is not concave.');
end

dv=-p(2)/(2*p(1));
if abs(dv)>halfWidth
    error('step17:VertexOutside','Full-domain fitted vertex left the fit window.');
end

P=struct();
P.phiStar=local_wrap(phi0+dv);
P.sigStar=polyval(p,dv);
P.rmse=sqrt(mean((y-polyval(p,x)).^2));
P.curvature=p(1);
P.halfWidth=halfWidth;
end


function M=local_symmetry(phi,sig)
pos=find(phi>0);
d=[];
for j=1:numel(pos)
    [dm,im]=min(abs(phi+phi(pos(j))));
    if dm<1e-10
        d(end+1,1)=abs(sig(pos(j))-sig(im)); %#ok<AGROW>
    end
end
scale=max(abs(sig));
if isempty(d)
    s=NaN;
else
    s=max(d)/max(scale,eps);
end
M=struct('symRel',s);
end


function a=local_wrap(a)
a=mod(a+pi,2*pi)-pi;
end


function local_plot_overlay(Half,Full,Nlist)
figure('Name','Step 17: half vs full boundary stress','Color','w');
clf; tiledlayout(numel(Nlist),1,'Padding','compact','TileSpacing','compact');
for il=1:numel(Nlist)
    nexttile; hold on; box on; grid on
    Bh=Half{il}.B;
    Bf=Full{il}.B;

    ph=rad2deg(Bh.phi);
    pf=rad2deg(local_wrap(Bf.phi));
    keep=abs(pf)<=90+1e-9;
    [x,ord]=sort(pf(keep));
    yf=Bf.sig_tt_eff(keep); yf=yf(ord);

    plot(ph,Bh.sig_tt_eff,'-','LineWidth',1.3,'DisplayName','half-domain');
    plot(x,yf,'--','LineWidth',1.1,'DisplayName','full-domain right half');
    xline(0,'k:','HandleVisibility','off');
    xlabel('\phi [deg]'); ylabel('\sigma_{tt}^{boundary}');
    title(sprintf('Npoly=%d',Nlist(il)));
    legend('Location','best');
end
end


function local_plot_difference(Half,Full,Nlist)
figure('Name','Step 17: half-full stress difference','Color','w');
clf; tiledlayout(numel(Nlist),1,'Padding','compact','TileSpacing','compact');
for il=1:numel(Nlist)
    nexttile; hold on; box on; grid on
    Bh=Half{il}.B;
    Bf=Full{il}.B;

    ph=Bh.phi(:);
    pf=local_wrap(Bf.phi(:));
    [pfs,ord]=sort(pf);
    sfs=Bf.sig_tt_eff(ord);
    sf=interp1(pfs,sfs,ph,'linear');

    plot(rad2deg(ph),Bh.sig_tt_eff-sf,'LineWidth',1.1);
    xline(0,'k:');
    xlabel('\phi [deg]');
    ylabel('\Delta\sigma_{tt}');
    title(sprintf('half - full, Npoly=%d',Nlist(il)));
end
end


function local_plot_convergence(T)
figure('Name','Step 17: symmetry-reduction convergence','Color','w');
clf; hold on; box on; grid on
semilogy(T.h_arc_over_R,abs(T.sig0_rel_diff),'-o','LineWidth',1.1, ...
    'DisplayName','|\Delta\sigma_{tt}(0)|');
semilogy(T.h_arc_over_R,T.curve_rel_Linf,'-s','LineWidth',1.1, ...
    'DisplayName','curve L_\infty');
semilogy(T.h_arc_over_R,T.curve_rel_L2,'-^','LineWidth',1.1, ...
    'DisplayName','curve L_2');
set(gca,'XDir','reverse');
xlabel('h_{arc}/R (refinement -> right)');
ylabel('relative difference');
legend('Location','best');
title('Half-domain versus full-domain convergence');
end


function local_plot_half_geometry(Rh)
G=Rh.G; S=Rh.S1;
figure('Name','Step 17: centered half-domain geometry','Color','w'); clf
hold on; axis equal; box on
triplot(G.t,G.p(:,1),G.p(:,2),'Color',[0.82 0.82 0.82]);
plot(G.holeArc(:,1),G.holeArc(:,2),'LineWidth',1.5,'DisplayName','semicircular free boundary');

sym=S.bc.symmetry_nodes;
plot(S.mesh.coord(sym,1),S.mesh.coord(sym,2),'o','MarkerSize',3, ...
    'DisplayName','u_x=0 symmetry nodes');
plot(Rh.I.x_star(1),Rh.I.x_star(2),'p','MarkerSize',10,'LineWidth',1.5, ...
    'DisplayName','fitted initiation point');

xlabel('x'); ylabel('y');
legend('Location','best');
title('Centered-hole right half-domain verification model');
end
