function Out=main_step4_graded_ring_asymmetry_experiment(varargin)
%MAIN_STEP4_GRADED_RING_ASYMMETRY_EXPERIMENT
% Controlled old-mirror/J versus FE-nodal EDI experiment on the graded
% concentric-ring mesh requested for the audit.
%
% The baseline mesh is made of geometrically growing concentric rings. Every
% ring is uniformly subdivided; adjacent rings alternate between M and M+1
% half-ring intervals to interlace nodes without half-size sectors at either
% x-axis. The radial growth ratio is chosen from a near-equilateral shape
% target and adjusted only enough to land exactly at r1.
%
% Cases:
%   S0       exact mirror-reflected lower half
%   A_shift  same nominal density; lower interior angular nodes shifted
%   A_1p5    lower-half angular spacing about 1.5 times larger
%   A_2      lower-half angular spacing 2 times larger
%
% Exact fields:
%   pure I        KI=1, KII=0
%   pure II       KI=0, KII=1
%   small mixed   KI=1, KII=0.05
%
% The old method is swept over circular contour radius. Canonical EDI uses
% FE-nodal q on a fixed annulus. Because the exact KI/KII are known, neither
% extractor is treated as ground truth.
%
% Usage:
%   O4 = main_step4_graded_ring_asymmetry_experiment();

ip=inputParser;
addParameter(ip,'PlotMeshes',true,@(x)islogical(x)||isnumeric(x));
addParameter(ip,'Ntheta',64,@(x)isnumeric(x)&&isscalar(x)&&x>=16);
addParameter(ip,'OldNtheta',240,@(x)isnumeric(x)&&isscalar(x)&&x>=40);
parse(ip,varargin{:}); opt=ip.Results;

here=fileparts(mfilename('fullpath'));
addpath(genpath(fileparts(fileparts(here))));

fprintf('\n============================================================\n');
fprintf('SIF AUDIT STEP 4: GRADED-RING ASYMMETRY EXPERIMENT\n');
fprintf('============================================================\n');

r0=0.005; r1=0.20;
oldR=[0.04 0.06 0.08 0.10 0.12];
rOldRef=0.08;
ediDomain=struct('r_inner',0.02,'r_outer',0.12);

E=4e3; nu=0.30; ps=1;
coef=E/((1+nu)*(1-2*nu));
D=coef*[1-nu,nu,0;nu,1-nu,0;0,0,(1-2*nu)/2];
mat=struct('E',E,'nu',nu,'ps',ps,'D',D,'Dmat',D);
V=[-1,0;0,0];

Kcases=[1,0;0,1;1,0.05];
caseNames=["pure_I";"pure_II";"small_mixed"];

spec={ ...
 'S0','mirror_reflected',0,1.0; ...
 'A_shift','lower_shift',0.20,1.0; ...
 'A_1p5','lower_coarse',0,1.5; ...
 'A_2','lower_coarse',0,2.0};

nM=size(spec,1); nK=size(Kcases,1);
Meshes=cell(nM,1); Info=cell(nM,1);
OldSweeps=cell(nM,nK);
rows=[]; meshRows=[];

for im=1:nM
    label=spec{im,1}; variant=spec{im,2};
    shift=spec{im,3}; coarse=spec{im,4};

    [mesh,info]=build_graded_ring_crack_mesh( ...
        'r0',r0,'r1',r1,'Ntheta',opt.Ntheta, ...
        'Variant',variant, ...
        'LowerShiftFraction',shift, ...
        'LowerCoarsenFactor',coarse);

    Meshes{im}=mesh; Info{im}=info;

    meshRows=[meshRows; im,info.Nr,info.Ntheta, ...
        info.nT3Vertices,info.nT3Elements,info.nT6Nodes,info.nDOF, ...
        info.qTargetNearEquilateral,info.qActual,info.qualityMin, ...
        info.qualityP05,info.qualityMedian,info.lowerAngularFactor, ...
        info.upperRingSegmentRelSpreadMax,info.lowerRingSegmentRelSpreadMax, ...
        info.nQualityBelow08]; %#ok<AGROW>

    fprintf('\n------------------------------------------------------------\n');
    fprintf('%s: %s\n',label,variant);
    fprintf('Nr=%d, Ntheta=%d, q_target=%.8f, q=%.8f\n', ...
        info.Nr,info.Ntheta,info.qTargetNearEquilateral,info.qActual);
    fprintf('T3=%d, T6 nodes=%d, DOF=%d\n', ...
        info.nT3Elements,info.nT6Nodes,info.nDOF);
    fprintf('triangle quality min/p05/median = %.4f / %.4f / %.4f\n', ...
        info.qualityMin,info.qualityP05,info.qualityMedian);
    fprintf('ring segment spread upper/lower = %.3e / %.3e\n', ...
        info.upperRingSegmentRelSpreadMax,info.lowerRingSegmentRelSpreadMax);
    fprintf('triangles with quality < 0.8 = %d\n',info.nQualityBelow08);
    fprintf('lower angular factor = %.5f\n',info.lowerAngularFactor);

    if logical(opt.PlotMeshes)
        figure('Name',['graded ring ',label],'NumberTitle','off');
        triplot(mesh.connect3,mesh.coord3(:,1),mesh.coord3(:,2));
        axis equal; grid on;
        xlabel('x_1'); ylabel('x_2');
        title(sprintf('%s: graded concentric-ring mesh',label));
    end

    for ic=1:nK
        KIin=Kcases(ic,1); KIIin=Kcases(ic,2);
        U=exact_williams_displacement_audit(mesh.coord,KIin,KIIin,E,nu,ps);

        [KIe,KIIe,AE]=SIF_LEFM_interaction_EDI( ...
            mesh,U,V,mat,ediDomain, ...
            'UsePlaneStrain',true,'Verbose',false, ...
            'WeightFunction','fe_nodal');

        O=nan(numel(oldR),5);
        dbg=cell(numel(oldR),1);
        for ir=1:numel(oldR)
            [KIo,KIIo,DO]=SIF_LEFM_circle2_debug( ...
                mesh,U,V,mat,oldR(ir), ...
                'nthet',opt.OldNtheta,'plot',false,'verbose',false);
            O(ir,:)=[oldR(ir),KIo,KIIo,KIo-KIin,KIIo-KIIin];
            dbg{ir}=DO;
        end

        Told=array2table(O,'VariableNames', ...
            {'rI','KI_old','KII_old','KI_error','KII_error'});
        OldSweeps{im,ic}=struct('table',Told,'debug',{dbg});

        [~,iref]=min(abs(oldR-rOldRef));
        KIo=Told.KI_old(iref); KIIo=Told.KII_old(iref);

        oldKIRange=max(Told.KI_old)-min(Told.KI_old);
        oldKIIRange=max(Told.KII_old)-min(Told.KII_old);

        trueNorm=hypot(KIin,KIIin);
        oldVecErr=hypot(KIo-KIin,KIIo-KIIin)/max(trueNorm,1);
        ediVecErr=hypot(KIe-KIin,KIIe-KIIin)/max(trueNorm,1);

        rows=[rows; im,ic,KIin,KIIin, ...
            KIo,KIIo,KIe,KIIe, ...
            KIo-KIin,KIIo-KIIin,KIe-KIin,KIIe-KIIin, ...
            oldVecErr,ediVecErr,oldKIRange,oldKIIRange, ...
            AE.nElem_used,AE.nGP_used]; %#ok<AGROW>
    end
end

Tmesh=array2table(meshRows,'VariableNames',{ ...
 'meshID','Nr','Ntheta','nT3_vertices','nT3_elements','nT6_nodes','nDOF', ...
 'q_target','q_actual','quality_min','quality_p05','quality_median', ...
 'lower_angular_factor','upper_ring_segment_spread', ...
 'lower_ring_segment_spread','n_quality_below_08'});
meshName=strings(height(Tmesh),1);
for i=1:height(Tmesh), meshName(i)=string(spec{Tmesh.meshID(i),1}); end
Tmesh.meshName=meshName; Tmesh=movevars(Tmesh,'meshName','After','meshID');

T=array2table(rows,'VariableNames',{ ...
 'meshID','caseID','KI_true','KII_true', ...
 'KI_old','KII_old','KI_EDI','KII_EDI', ...
 'KI_old_error','KII_old_error','KI_EDI_error','KII_EDI_error', ...
 'old_vector_error','EDI_vector_error','old_KI_radius_range', ...
 'old_KII_radius_range','EDI_nElem','EDI_nGP'});
meshName=strings(height(T),1); caseName=strings(height(T),1);
for i=1:height(T)
    meshName(i)=string(spec{T.meshID(i),1});
    caseName(i)=caseNames(T.caseID(i));
end
T.meshName=meshName; T.caseName=caseName;
T=movevars(T,'meshName','After','meshID');
T=movevars(T,'caseName','After','caseID');

fprintf('\n============================================================\n');
fprintf('TOTAL MESH SUMMARY\n');
fprintf('============================================================\n');
disp(Tmesh);

fprintf('\n============================================================\n');
fprintf('EXACT-FIELD EXTRACTION SUMMARY\n');
fprintf('old method shown at rI=%.4f; EDI annulus=[%.4f,%.4f]\n', ...
    rOldRef,ediDomain.r_inner,ediDomain.r_outer);
fprintf('============================================================\n');
disp(T);

for im=1:nM
    for ic=1:nK
        fprintf('\n%s / %s : old contour sweep\n', ...
            spec{im,1},caseNames(ic));
        disp(OldSweeps{im,ic}.table);
    end
end

Out=struct();
Out.meshSummary=Tmesh;
Out.results=T;
Out.meshes=Meshes;
Out.info=Info;
Out.oldSweeps=OldSweeps;
Out.settings=struct('r0',r0,'r1',r1,'oldR',oldR, ...
    'rOldRef',rOldRef,'ediDomain',ediDomain, ...
    'Kcases',Kcases,'caseNames',caseNames);

fprintf('\nSTEP 4 completed.\n');
fprintf(['Interpret pure-mode false cross terms and small-mixed errors against ', ...
    'the known exact SIFs; do not treat either extractor as truth.\n']);
end
