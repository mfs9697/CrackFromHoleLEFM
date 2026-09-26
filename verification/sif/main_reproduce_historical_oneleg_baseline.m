function Results = main_reproduce_historical_oneleg_baseline(varargin)
%MAIN_REPRODUCE_HISTORICAL_ONELEG_BASELINE  Reproduce the published LEFM path.
%
% This is the historical baseline actually supported by the IJF source:
% the SIFs are evaluated at the tip of Leg 1 (the physical crack), before
% the cohesive Leg 2 is introduced into the LEFM predictor.
%
% The old production call used rI = 0.5*a.  We also sweep nearby radii as a
% diagnostic, but Historical is the rI/a=0.5 row.
%
% Name-value options
%   'Sigma0'       nominal remote stress in MPa (default 1)
%   'rOverA'       default [0.2 0.3 0.4 0.5 0.6]
%   'SaveOutputs'  true
%   'OutputDir'    verification/sif/outputs

    ip=inputParser;
    addParameter(ip,'Sigma0',1.0,@(x)isnumeric(x)&&isscalar(x)&&isfinite(x));
    addParameter(ip,'rOverA',[0.2 0.3 0.4 0.5 0.6],@isnumeric);
    addParameter(ip,'SaveOutputs',true,@(x)islogical(x)||isnumeric(x));
    addParameter(ip,'OutputDir',fullfile('verification','sif','outputs'),@(x)ischar(x)||isstring(x));
    parse(ip,varargin{:});

    cm=.01;
    C=struct();
    C.A=10*cm;
    C.B=10*cm;
    C.a=2*cm;
    C.eps1=1e-8;
    C.hmax=C.B/40;
    C.ncoh=40;
    C.hgrad=1.1;
    C.da=0.06*C.a;
    C.chw=(1/8)*C.da/C.ncoh;
    C.theta_deg=20;
    C.alf1=-20*pi/180;
    C.alf2=C.theta_deg*pi/180;
    C.V0=[0,0];
    C.V1=C.a*[cos(C.alf1),sin(C.alf1)];
    C.V2=C.V1+C.da*[cos(C.alf2),sin(C.alf2)];
    C.nu=.3;
    C.E2=4e3;
    C.Dmat=C.E2*[1-C.nu,C.nu,0; C.nu,1-C.nu,0; 0,0,(1-2*C.nu)/2] ...
        /(1+C.nu)/(1-2*C.nu);
    C.G12=C.E2/2/(1+C.nu);
    C.plotMesh=false;

    % Historical one-leg meshing call uses arrow=0.005.
    [mesh,mat,~,elod,K,fixvar]=historical_geom_pencil_1leg(C,0.005);
    mat.Dmat=C.Dmat;
    mat.ps=1;

    F=zeros(2*size(mesh.coord,1),1);
    ids=elod(:,1);
    F(2*ids)=F(2*ids)+ip.Results.Sigma0*elod(:,2);
    F(fixvar)=0;
    U=K\F;

    V=[C.V0;C.V1];
    rr=ip.Results.rOverA(:);
    rows=zeros(numel(rr),10);

    for k=1:numel(rr)
        rI=rr(k)*C.a;
        [KI,KII,Dbg]=SIF_LEFM_circle2_debug(mesh,U,V,mat,rI, ...
            'nthet',240,'eps_th',1e-3,'plot',false,'verbose',false);
        [th,thdeg]=kink_angle_LEFM_MTS(KI,KII);
        M=SIF_mesh_mirror_metrics(mesh,V,rI,'theta',Dbg.theta);

        rows(k,:)=[rr(k),rI,KI,KII,KII/KI,th,thdeg, ...
            M.A_h_median,Dbg.diagnostics.min_baryP,Dbg.diagnostics.min_baryQ];
    end

    T=array2table(rows,'VariableNames',{ ...
        'rI_over_a','rI','KI_old','KII_old','KII_over_KI', ...
        'thetaMTS_rad','thetaMTS_deg','A_h_median','minBaryP','minBaryQ'});

    [~,ih]=min(abs(T.rI_over_a-0.5));
    Historical=T(ih,:);

    Results=struct();
    Results.C=C;
    Results.Sigma0=ip.Results.Sigma0;
    Results.Table=T;
    Results.Historical=Historical;
    Results.Mesh=mesh;
    Results.U=U;

    fprintf('\nHistorical one-leg LEFM baseline (old mirror-J):\n');
    fprintf('  sigma0 = %.8g MPa\n',ip.Results.Sigma0);
    fprintf('  rI/a   = %.3f\n',Historical.rI_over_a);
    fprintf('  KI      = %.10e MPa*sqrt(m)\n',Historical.KI_old);
    fprintf('  KII     = %+.10e MPa*sqrt(m)\n',Historical.KII_old);
    fprintf('  MTS     = %+.8f deg\n',Historical.thetaMTS_deg);

    if logical(ip.Results.SaveOutputs)
        outdir=char(ip.Results.OutputDir);
        if ~exist(outdir,'dir'),mkdir(outdir);end
        writetable(T,fullfile(outdir,'historical_oneleg_old_baseline_radius_sweep.csv'));
        writetable(Historical,fullfile(outdir,'historical_oneleg_old_baseline.csv'));
        save(fullfile(outdir,'historical_oneleg_old_baseline.mat'),'Results','-v7.3');
    end
end
