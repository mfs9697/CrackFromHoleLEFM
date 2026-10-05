function Probe=test_incremental_fast_edi(varargin)
% Protect optional unused-work removal, including diagnostic fallback.
ip=inputParser;
addParameter(ip,'CheckpointFile','',@(x)ischar(x)||isstring(x));
parse(ip,varargin{:});
[p,t]=deal([-2,-1;2,-1;2,1;-2,1],[1,2,3;1,3,4]);
[p6,t6]=T3toT6_fast(p,t);mesh=struct('coord',p6,'connect',t6);
U=reshape([sin(p6(:,1)),cos(p6(:,2))].',[],1);
E=100;nu=.25;D=E/((1+nu)*(1-2*nu))*[1-nu,nu,0;nu,1-nu,0;0,0,(1-2*nu)/2];
mat=struct('E',E,'nu',nu,'D',D,'ps',1);
V=[-1,0;0,0];domain=struct('r_inner',.2,'r_outer',1.8);
cases=0;
for rule=[7,12,16]
    for weight={'analytic_radial','fe_nodal'}
        for diagnostics=[false,true]
            for exact=[false,true]
                analytic=[];if exact,analytic=[.7,-.2];end
                args={'QuadratureRule',rule,'WeightFunction',weight{1}, ...
                    'StoreGPDiagnostics',diagnostics,'AnalyticActualK',analytic,'AuxK',.37};
                [ki,kii,a]=SIF_LEFM_interaction_EDI(mesh,U,V,mat,domain,args{:});
                [kf,kff,b]=SIF_LEFM_interaction_EDI(mesh,U,V,mat,domain,args{:},'SkipUnusedAuxWork',true);
                assert(isequaln([ki,kii],[kf,kff])&&isequaln(a,b), ...
                    'fastedi:Bitwise','Integral or diagnostics changed.');
                cases=cases+1;
            end
        end
    end
end
Probe=struct('fixtureCases',cases,'bitwiseIdentical',true);
file=char(ip.Results.CheckpointFile);
if ~isempty(file)
    cp=load(file,'mesh','U','mat','crack','currentIncrement');
    domain=struct('r_inner',.10*cp.currentIncrement,'r_outer',.65*cp.currentIncrement);
    args={'UsePlaneStrain',cp.mat.ps==1,'WeightFunction','fe_nodal', ...
        'QuadratureRule',16,'StoreGPDiagnostics',false};
    % Warm both paths outside the paired timing samples.
    [ki,kii,a]=SIF_LEFM_interaction_EDI(cp.mesh,cp.U,cp.crack.Pmid,cp.mat,domain,args{:});
    [kf,kff,b]=SIF_LEFM_interaction_EDI(cp.mesh,cp.U,cp.crack.Pmid,cp.mat,domain,args{:},'SkipUnusedAuxWork',true);
    assert(isequaln([ki,kii],[kf,kff])&&isequaln(a,b)&&a.nElem_used==11316);
    timings=zeros(3,2);
    for j=1:3
        order=[1,2];if mod(j,2)==0,order=[2,1];end
        for mode=order
            timer=tic;
            [x,y,z]=SIF_LEFM_interaction_EDI(cp.mesh,cp.U,cp.crack.Pmid,cp.mat,domain,args{:},'SkipUnusedAuxWork',mode==2);
            timings(j,mode)=toc(timer);
            assert(isequaln([x,y],[ki,kii])&&isequaln(z,a));
        end
    end
    Probe.pairedSeconds=timings;Probe.medianSeconds=median(timings,1);
    Probe.speedup=Probe.medianSeconds(1)/Probe.medianSeconds(2);
    Probe.KI=ki;Probe.KII=kii;Probe.support=a.nElem_used;
    disp(array2table(timings,'VariableNames',{'strictSeconds','fastSeconds'}));
    fprintf('Stored P5 bitwise PASS; median EDI speedup %.6fx.\n',Probe.speedup);
end
fprintf('PASS: %d EDI equivalence cases, including diagnostics-enabled fallback.\n',cases);
end
