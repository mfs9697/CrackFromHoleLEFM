function O29=main_step29_crack_face_COD_audit(O25,O26,O27,varargin)
%MAIN_STEP29_CRACK_FACE_COD_AUDIT
% Independent displacement-jump estimate of crack-tip SIFs on the SAME
% solved T6 mesh used in Steps 25-28. NO new mesh or FEM solve.
%
% Compare actual FEM displacement jump to the exact unit Mode-I and
% Mode-II Williams fields already sampled on this mesh in Step 26.
% Upper/lower T6 crack-face node classes from O26 are essential:
% duplicated crack faces may have coincident coordinates but DIFFERENT U.
% For a straight traction-free crack with local tip coordinate -r:
%   jump(u_local_2) = (kappa+1)/mu * sqrt(r/(2*pi)) * KI + O(r^(3/2))
%   jump(u_local_1) = (kappa+1)/mu * sqrt(r/(2*pi)) * KII + O(r^(3/2))
% so Kapp(r) = mu/(kappa+1)*sqrt(2*pi/r) * jump(r).
%
% Sample each face separately with pchip interpolation at common physical
% distances behind the tip, then fit Kapp(r) = K0 + c*(r/a0) over
% several windows. Raw K0 is the asymptotic COD estimate, NOT corrected by
% the exact-Williams same-mesh calibration matrix. The synthetic matrix
% is only a numerical check for signs/interpolation of the COD extractor.
% A stable K0 requires both sampling and fit-window convergence, and
% it need not agree with finite-radius EDI until EDI path dependence
% is resolved.
%
% Usage:
%   O29=main_step29_crack_face_COD_audit(O25,O26,O27);

p=inputParser;
addParameter(p,'MinRoverA0',0.04, ...
    @(v)isnumeric(v)&&isscalar(v)&&isfinite(v)&&v>0&&v<0.2);
addParameter(p,'MaxRoverA0',0.40, ...
    @(v)isnumeric(v)&&isscalar(v)&&isfinite(v)&&v>=0.2&&v<0.8);
addParameter(p,'FitUpperBounds',[0.12 0.20 0.30 0.40], ...
    @(v)isnumeric(v)&&isvector(v)&&all(isfinite(v))&&all(v>0 & v<1));
addParameter(p,'NQuery',61, ...
    @(v)isnumeric(v)&&isscalar(v)&&isfinite(v)&&v==round(v)&&v>=15);
addParameter(p,'Plot',true,@(v)islogical(v)||isnumeric(v));
parse(p,varargin{:});
O=p.Results;

requiredO25={'mesh','U','mat','crack','hTip'};
requiredO26={'UI_exact','UII_exact','crack_face_side'};
for i=1:numel(requiredO25),must(O25,requiredO25{i},'O25');end
for i=1:numel(requiredO26),must(O26,requiredO26{i},'O26');end

mesh=O25.mesh;
V=O25.crack.Pmid;
n=size(mesh.coord,1);
if numel(O25.U)~=2*n || numel(O26.UI_exact)~=2*n || ...
        numel(O26.UII_exact)~=2*n || ...
        numel(O26.crack_face_side)~=n
    error('step29:MeshMismatch', ...
        'O25 and O26 must contain the same T6 mesh/DOF ordering.');
end
tip=V(end,:);
ev=tip-V(end-1,:);
ev=ev/norm(ev);
en=[-ev(2),ev(1)];
Rgl=[ev(:),en(:)];
xl=(mesh.coord-tip)*Rgl;
a0=sum(sqrt(sum(diff(V,1,1).^2,2)));
if size(V,1)~=2
    error('step29:NeedStraightCrack', ...
        'This COD audit currently assumes a straight 2-point crack.');
end

face=O26.crack_face_side(:);
faceTol=max(1e-12,1e-8*a0);
onFace=(abs(xl(:,2))<faceTol) & xl(:,1)<-faceTol ...
       & -xl(:,1)<=a0+faceTol;
uIdx=find(onFace & face==+1);
lIdx=find(onFace & face==-1);
if numel(uIdx)<6||numel(lIdx)<6
    error('step29:InsufficientFaceNodes', ...
        'Insufficient classified T6 crack-face nodes (up=%d,lo=%d).', ...
        numel(uIdx),numel(lIdx));
end

% Distances r are POSITIVE behind the crack tip.
rUpper=-xl(uIdx,1);
rLower=-xl(lIdx,1);
[rUpper,ju]=sort(rUpper);
[rLower,jl]=sort(rLower);
uIdx=uIdx(ju);
lIdx=lIdx(jl);

% Duplicated abscissae cannot be passed to pchip; their exact U values
% should coincide on the same face. Average only within that SAME face.
[rUpper,~,jU]=unique(rUpper);
[rLower,~,jL]=unique(rLower);
queryFrac=linspace(O.MinRoverA0,O.MaxRoverA0,O.NQuery).';
rQuery=a0*queryFrac;
if rQuery(1)<max(min(rUpper),min(rLower)) || ...
        rQuery(end)>min(max(rUpper),max(rLower))
    error('step29:OutsideFaceSupport', ...
        ['Requested r interval [%.4g,%.4g] m is outside overlap ', ...
        'of upper/lower crack-face node supports [%.4g,%.4g] m.'], ...
        rQuery(1),rQuery(end), ...
        max(min(rUpper),min(rLower)),min(max(rUpper),max(rLower)));
end
if rQuery(1)<2*O25.hTip
    warning('step29:NearTipSample', ...
        'Smallest r=%.4g is below two median tip edges %.4g.', ...
        rQuery(1),2*O25.hTip);
end

mu=O25.mat.E/(2*(1+O25.mat.nu));
if O25.mat.ps==1
    kappa=3-4*O25.mat.nu;
else
    kappa=(3-O25.mat.nu)/(1+O25.mat.nu);
end
factor=mu/(kappa+1)*sqrt(2*pi./rQuery);

Uset={O25.U,O26.UI_exact,O26.UII_exact};
names={'Actual FEM','Exact unit I','Exact unit II'};
profiles=nan(O.NQuery,2,numel(Uset));
jumps=nan(O.NQuery,2,numel(Uset));
for ifield=1:numel(Uset)
    u=reshape(Uset{ifield},2,[]).';
    ul=u*Rgl;
    % Both U values and interpolation abscissae are local crack-tip
    % quantities; positive e2 is geometrically above the crack.
    uUpper=nan(numel(rUpper),2);
    uLower=nan(numel(rLower),2);
    for ic=1:2
        valU=ul(uIdx,ic);
        valL=ul(lIdx,ic);
        uUpper(:,ic)=accumarray(jU,valU,[],@mean);
        uLower(:,ic)=accumarray(jL,valL,[],@mean);
    end
    jmp=interp1(rUpper,uUpper,rQuery,'pchip') ...
       -interp1(rLower,uLower,rQuery,'pchip');
    if any(~isfinite(jmp),'all')
        error('step29:BadInterpolation', ...
            'Nonfinite displacement jump for field %s.',names{ifield});
    end
    % Column 1 is KI from normal jump; column 2 is signed KII from
    % tangential jump. Unit exact modes validate these signs.
    jumps(:,:,ifield)=[jmp(:,2),jmp(:,1)];
    profiles(:,:,ifield)=bsxfun(@times,jumps(:,:,ifield),factor);
end

upperBounds=O.FitUpperBounds(:).';
upperBounds=sort(unique(upperBounds));
if upperBounds(1)<=O.MinRoverA0 || upperBounds(end)>O.MaxRoverA0
    error('step29:BadFitWindows', ...
        'Every fit upper bound must be > MinRoverA0 and <= MaxRoverA0.');
end
nW=numel(upperBounds);
KFit=nan(2,numel(Uset),nW);
slopes=nan(2,numel(Uset),nW);
fitRows=nan(nW,12);
fprintf('\n============================================================\n');
fprintf('STEP 29: SAME-MESH CRACK-FACE DISPLACEMENT-JUMP AUDIT\n');
fprintf('============================================================\n');
fprintf('  a0=%.8e; phi_tip=%+.8f deg; upper/lower T6 nodes=%d/%d\n', ...
    a0,atan2d(ev(2),ev(1)),numel(uIdx),numel(lIdx));
fprintf('  r/a0 sampled=[%.3f,%.3f], n=%d; actual h_tip/a0=%.5f\n', ...
    queryFrac(1),queryFrac(end),numel(queryFrac),O25.hTip/a0);
fprintf('  Asymptotic COD: KI from normal jump, signed KII from tangential jump.\n');

for iw=1:nW
    j=find(queryFrac<=upperBounds(iw)+1e-12);
    if numel(j)<6
        error('step29:TooFewWindowPoints', ...
            'At least six samples needed in fit window %.3f.',upperBounds(iw));
    end
    for ifield=1:numel(Uset)
        for mode=1:2
            poly=polyfit(queryFrac(j),profiles(j,mode,ifield),1);
            KFit(mode,ifield,iw)=poly(2);
            slopes(mode,ifield,iw)=poly(1);
        end
    end
    % Synthetic recovery matrix from both independently sampled unit modes:
    % columns are exact unit I and exact unit II.
    synthetic=[KFit(:,2,iw),KFit(:,3,iw)];
    raw=KFit(:,1,iw);
    [~,kRef]=min(abs(O27.rOuterOverA0(:)-0.65));
    ir16=find(O27.rules==16,1);
    if isempty(ir16)
        edi16=[NaN;NaN];
    else
        edi16=O27.Kactual(:,ir16,kRef);
    end
    fitRows(iw,:)=[upperBounds(iw),numel(j), ...
        raw(1),raw(2),raw(2)/raw(1), ...
        norm(synthetic-eye(2),'fro'), ...
        synthetic(2,1),synthetic(1,2), ...
        synthetic(1,1),synthetic(2,2), ...
        edi16(1),edi16(2)];

    fprintf([' r/a0=[%.3f,%.3f] | COD intercept KI=%.8e, ', ...
        'KII=%+.8e, KII/KI=%+.6e | synthetic ||M-I||=%.3e\n'], ...
        O.MinRoverA0,upperBounds(iw),raw(1),raw(2), ...
        raw(2)/raw(1),norm(synthetic-eye(2),'fro'));
end

T=array2table(fitRows,'VariableNames',{ ...
    'fit_upper_r_over_a0','n_points', ...
    'KI_COD_raw','KII_COD_raw','KII_over_KI_COD_raw', ...
    'synthetic_recovery_error','synthetic_KII_leak_from_I', ...
    'synthetic_KI_leak_from_II','synthetic_I_recovery', ...
    'synthetic_II_recovery', ...
    'EDI16_KI_ref','EDI16_KII_ref'});

fprintf('\nCOD FIT-WINDOW RESULTS\n');disp(T);
fprintf(['Interpretation: COD extraction is independent of EDI and ', ...
    'assesses the SAME FEM displacement field. Synthetic unit-mode ', ...
    'COD recovery checks its own interpolation and fitting error. ', ...
    'Window convergence is required before interpreting its ', ...
    'near-zero signed KII as physical.\n']);

if logical(O.Plot)
    figure('Name','Step29 independent crack-face COD SIF profiles','Color','w');
    t=tiledlayout(1,2,'Padding','compact','TileSpacing','compact');
    nexttile;
    plot(queryFrac,profiles(:,1,1),'-o','MarkerSize',3);hold on;
    plot(queryFrac,profiles(:,1,2),'--');
    grid on;box on;xlabel('r/a_0');ylabel('apparent K_I');
    title('Crack-face normal displacement jump');
    legend('Actual FEM','Exact mode I','Location','best');
    nexttile;
    plot(queryFrac,profiles(:,2,1),'-o','MarkerSize',3);hold on;
    plot(queryFrac,profiles(:,2,2),'--');
    plot(queryFrac,profiles(:,2,3),'-.');
    grid on;box on;xlabel('r/a_0');ylabel('apparent signed K_{II}');
    title('Crack-face tangential displacement jump');
    legend('Actual FEM','Exact mode I leakage','Exact mode II','Location','best');
    title(t,'Independent local displacement-jump verification');
end

O29=struct();
O29.settings=O;
O29.rOverA0=queryFrac;
O29.tip=tip;
O29.axis=ev;
O29.nUpper=numel(uIdx);
O29.nLower=numel(lIdx);
O29.profiles=profiles;
O29.jumps=jumps;
O29.fieldNames=names;
O29.KFit=KFit;
O29.KFitSlope=slopes;
O29.table=T;
fprintf('STEP 29 completed; no new FEM calculations.\n');
end

function must(S,key,name)
if ~isstruct(S)||~isfield(S,key)||isempty(S.(key))
    error('step29:MissingInput','Missing %s.%s.',name,key);
end
end
