function Results = main_twoleg_mirror_sign_check(varargin)
%MAIN_TWOLEG_MIRROR_SIGN_CHECK  Reflect the full two-leg benchmark in y=0.
%
% Under symmetric remote vertical tension, reflection should preserve KI
% and reverse KII when the local crack frame is defined consistently.

    ip=inputParser;
    addParameter(ip,'MeshScale',1.0,@(x)isnumeric(x)&&isscalar(x)&&x>0);
    addParameter(ip,'rOverDelta',[0.25 0.35 0.45 0.55],@isnumeric);
    addParameter(ip,'InnerFraction',0.20,@(x)isnumeric(x)&&isscalar(x)&&x>=0&&x<1);
    addParameter(ip,'SaveOutputs',true,@(x)islogical(x)||isnumeric(x));
    addParameter(ip,'OutputDir',fullfile('verification','sif','outputs'),@(x)ischar(x)||isstring(x));
    parse(ip,varargin{:});

    Cplus=cfg_twoleg_historical();
    Cminus=mirror_cfg(Cplus);

    SP=build_twoleg_sharp_benchmark(Cplus,'MeshScale',ip.Results.MeshScale, ...
        'AsymmetryLevel',0,'Plot',false);
    SM=build_twoleg_sharp_benchmark(Cminus,'MeshScale',ip.Results.MeshScale, ...
        'AsymmetryLevel',0,'Plot',false);

    rows=[];
    for rr=ip.Results.rOverDelta(:).'
        r=rr*Cplus.delta;

        [KIoP,KIIoP]=SIF_LEFM_circle2_debug(SP.mesh,SP.U,SP.V,SP.mat,r, ...
            'nthet',Cplus.nthet,'eps_th',Cplus.eps_th);
        [KIoM,KIIoM]=SIF_LEFM_circle2_debug(SM.mesh,SM.U,SM.V,SM.mat,r, ...
            'nthet',Cminus.nthet,'eps_th',Cminus.eps_th);

        D=struct('r_inner',ip.Results.InnerFraction*r,'r_outer',r);
        [KIeP,KIIeP]=SIF_LEFM_interaction_EDI(SP.mesh,SP.U,SP.V,SP.mat,D);
        [KIeM,KIIeM]=SIF_LEFM_interaction_EDI(SM.mesh,SM.U,SM.V,SM.mat,D);

        rows=[rows;rr,r, ...
            KIoP,KIIoP,KIoM,KIIoM, ...
            abs(KIoP-KIoM)/max(abs(0.5*(KIoP+KIoM)),eps), ...
            abs(KIIoP+KIIoM)/max(abs(0.5*(KIIoP-KIIoM)),eps), ...
            KIeP,KIIeP,KIeM,KIIeM, ...
            abs(KIeP-KIeM)/max(abs(0.5*(KIeP+KIeM)),eps), ...
            abs(KIIeP+KIIeM)/max(abs(0.5*(KIIeP-KIIeM)),eps)]; %#ok<AGROW>
    end

    T=array2table(rows,'VariableNames',{ ...
        'rOverDelta','r', ...
        'KIold_original','KIIold_original','KIold_mirror','KIIold_mirror', ...
        'KIold_even_error','KIIold_odd_error', ...
        'KIEDI_original','KIIEDI_original','KIEDI_mirror','KIIEDI_mirror', ...
        'KIEDI_even_error','KIIEDI_odd_error'});

    Results=struct('OriginalConfig',Cplus,'MirrorConfig',Cminus,'Table',T);

    if logical(ip.Results.SaveOutputs)
        outdir=char(ip.Results.OutputDir);
        if ~exist(outdir,'dir'),mkdir(outdir);end
        writetable(T,fullfile(outdir,'twoleg_full_geometry_mirror_sign_check.csv'));
        save(fullfile(outdir,'twoleg_full_geometry_mirror_sign_check.mat'),'Results','-v7.3');
    end
end


function C=mirror_cfg(C)
    C.theta1=-C.theta1;
    C.theta2=-C.theta2;
    C.theta1_deg=rad2deg(C.theta1);
    C.theta2_deg=rad2deg(C.theta2);

    C.V0=[0,0];
    C.V1=C.V0+C.a*[cos(C.theta1),sin(C.theta1)];
    C.V2=C.V1+C.delta*[cos(C.theta2),sin(C.theta2)];
    C.Pmid=[C.V0;C.V1;C.V2];

    C.e1=[cos(C.theta2);sin(C.theta2)];
    C.e2=[-sin(C.theta2);cos(C.theta2)];
    C.R_gl=[C.e1,C.e2];
    C.R_loc=C.R_gl.';
end
