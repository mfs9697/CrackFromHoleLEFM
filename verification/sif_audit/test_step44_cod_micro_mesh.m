function report=test_step44_cod_micro_mesh()
%TEST_STEP44_COD_MICRO_MESH Cheap exact COD extractor regression.
% Self-contained synthetic T6 crack-face topology; no project checkpoint,
% stiffness assembly, EDI integration, FEM solve or mesh generator.
% Tests pure I, pure II and very small mixed II independently against
% imposed leading-order Williams displacements at actual crack-face nodes.
%
% This verifies the COD formula and face bookkeeping on a CONTROLLED
% mesh; it does NOT establish accuracy for a computed physical FEM field.
a0=.008;
nseg=8;
r=linspace(0,a0,nseg+1);
P=[0 0]; % shared physical mathematical tip, one T3 node ID
up=zeros(nseg+1,1);lo=zeros(nseg+1,1);
up(1)=1;lo(1)=1;
for j=2:nseg+1
    up(j)=size(P,1)+1;
    P(end+1,:)=[-r(j) 0]; %#ok<AGROW>
    lo(j)=size(P,1)+1;
    P(end+1,:)=[-r(j) 0]; %#ok<AGROW>
end
T=zeros(2*nseg,3);
for j=1:nseg
    pMid=-0.5*(r(j)+r(j+1));
    apexU=size(P,1)+1;
    P(end+1,:)=[pMid .20*a0]; %#ok<AGROW>
    apexL=size(P,1)+1;
    P(end+1,:)=[pMid -.20*a0]; %#ok<AGROW>
    T(2*j-1,:)=[up(j),up(j+1),apexU];
    T(2*j,:)=[lo(j),lo(j+1),apexL];
end
% Dedicated T6 midside IDs per triangle. Their physical coordinates
% are exact arithmetic edge midpoints; faces remain topologically split.
T6=zeros(size(T,1),6);
T6(:,1:3)=T;
for j=1:size(T,1)
    for k=1:3
        e=[1 2;2 3;3 1];
        p1=P(T(j,e(k,1)),:);
        p2=P(T(j,e(k,2)),:);
        T6(j,k+3)=size(P,1)+1;
        P(end+1,:)=(p1+p2)/2; %#ok<AGROW>
    end
end
mesh=struct('coord',P,'connect',T6);
mat=struct('E',210e9,'nu',.30,'ps',1);
crack=struct('Pmid',[-a0 0;0 0], ...
    'tipNode',1,'upperNodes',up,'lowerNodes',lo);
emptyU=zeros(2*size(P,1),1);
[~,~,f]=native_COD_audit(mesh,emptyU,mat,crack,8,true);
face=f.faceSide;
ids=find(face~=0);
assert(numel(ids)==4*nseg,'Expected duplicated crack-face T6 samples.');
testK=[1 0;0 1;.43784 4.7269e-5];
maxErr=zeros(3,2);
nNative=zeros(3,1);
for j=1:3
    uf=exact_williams_displacement_audit(P(ids,:), ...
        testK(j,1),testK(j,2),mat.E,mat.nu,mat.ps, ...
        'UpperFaceIDs',find(face(ids)==+1), ...
        'LowerFaceIDs',find(face(ids)==-1));
    u=zeros(2*size(P,1),1);
    u(2*ids-1)=uf(1:2:end);
    u(2*ids)=uf(2:2:end);
    [rr,app]=native_COD_audit(mesh,u,mat,crack,8);
    nNative(j)=numel(rr);
    maxErr(j,:)=max(abs(bsxfun(@minus,app,testK(j,:))),[],1);
end
pass=all(maxErr(1,:)<1e-9) && ...
    all(maxErr(2,:)<1e-9) && ...
    maxErr(3,1)<1e-9 && ...
    maxErr(3,2)/abs(testK(3,2))<1e-4 && ...
    all(nNative==2*nseg);
report=struct('passed',pass,'testK',testK, ...
    'nNative',nNative,'maxAbsError',maxErr);
fprintf('\nSTEP 44 SYNTHETIC SMALL-MESH COD TEST\n');
disp(report);
assert(pass,'Step44:MicroMeshTestFailed', ...
    'Exact COD recovery failed on the controlled small T6 mesh.');
end
