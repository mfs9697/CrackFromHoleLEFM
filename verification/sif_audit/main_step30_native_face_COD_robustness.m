function O30=main_step30_native_face_COD_robustness(O25,O26,O27,O29,varargin)
%MAIN_STEP30_NATIVE_FACE_COD_ROBUSTNESS
% Audit the independent crack-face COD SIF estimates without generating
% 61 correlated points from a much smaller number of native face nodes.
%
% Use each face's ACTUAL T6 nodal distances behind the tip as sampling
% abscissae. Interpolate only the OPPOSITE face to those abscissae, then
% compare upper-native vs lower-native COD fits over varied physical
% windows and linear vs quadratic terms in r/a0.
%
% All fields are from the EXISTING Step-25 FEM solution and Step-26 exact
% unit-I/II synthetic fields; no new FEM calculation or EDI evaluation.
% The fitted zero-radius intercept is an extrapolation. Stability across
% low/high radius limits, polynomial orders, and native-face sampling
% directions is a NECESSARY but not sufficient condition for physical
% resolution of a tiny signed KII.
%
% Usage:
%   O30=main_step30_native_face_COD_robustness(O25,O26,O27,O29);

p=inputParser;
addParameter(p,'LowerRoverA0',[0.04 0.08], ...
    @(v)isnumeric(v)&&isvector(v)&&all(isfinite(v))&&all(v>0)&&all(v<0.2));
addParameter(p,'UpperRoverA0',[0.20 0.30 0.40], ...
    @(v)isnumeric(v)&&isvector(v)&&all(isfinite(v))&&all(v>0.15)&&all(v<0.8));
addParameter(p,'PolynomialDegrees',[1 2], ...
    @(v)isnumeric(v)&&isvector(v)&&all(ismember(v,[1 2])));
addParameter(p,'MinNativePoints',8, ...
    @(v)isnumeric(v)&&isscalar(v)&&isfinite(v)&&v>=5&&v==round(v));
addParameter(p,'Plot',true,@(v)islogical(v)||isnumeric(v));
parse(p,varargin{:});
P=p.Results;

req25={'mesh','U','mat','crack','hTip'};
req26={'UI_exact','UII_exact','crack_face_side'};
for z=1:numel(req25),must(O25,req25{z},'O25');end
for z=1:numel(req26),must(O26,req26{z},'O26');end

V=O25.crack.Pmid;
if size(V,1)~=2
    error('step30:NotStraight','Native-face COD currently requires a straight 2-point crack.');
end
coord=O25.mesh.coord;
n=size(coord,1);
assert(numel(O25.U)==2*n && numel(O26.UI_exact)==2*n && ...
       numel(O26.UII_exact)==2*n && numel(O26.crack_face_side)==n, ...
       'O25 and O26 must have the same T6 nodal order.');
tip=V(end,:);
a0=norm(V(2,:)-V(1,:));
e1=(V(2,:)-V(1,:))/a0;
e2=[-e1(2),e1(1)];
R=[e1(:),e2(:)];
xl=(coord-tip)*R;
face=O26.crack_face_side(:);
tol=max(1e-12,1e-8*a0);
onFace=(abs(xl(:,2))<tol & xl(:,1)<-tol & -xl(:,1)<=a0+tol);
up=find(onFace & face==1);
lo=find(onFace & face==-1);
if numel(up)<P.MinNativePoints || numel(lo)<P.MinNativePoints
    error('step30:TooFewFaces','Need more classified T6 crack-face nodes.');
end
[rUp,iu]=sort(-xl(up,1));
[rLo,il]=sort(-xl(lo,1));
up=up(iu);lo=lo(il);
[rUp,~,gUp]=unique(rUp);
[rLo,~,gLo]=unique(rLo);
lower=sort(unique(P.LowerRoverA0(:).'));
upper=sort(unique(P.UpperRoverA0(:).'));
degrees=sort(unique(P.PolynomialDegrees(:).'));
overlap=[max(min(rUp),min(rLo)),min(max(rUp),max(rLo))];
if overlap(2)<=overlap(1)
    error('step30:NoOverlap','Upper and lower crack-face abscissae do not overlap.');
end

mat=O25.mat;
mu=mat.E/(2*(1+mat.nu));
if mat.ps==1
    kappa=3-4*mat.nu;
else
    kappa=(3-mat.nu)/(1+mat.nu);
end
Usets={O25.U,O26.UI_exact,O26.UII_exact};
fields={'Actual FEM','Exact I','Exact II'};
nativeR={rUp,rLo};
nativeProfiles=cell(2,1);
for side=1:2
    r=nativeR{side};
    % Evaluation at native positions of the selected face; interpolate
    % the opposite face ONLY. This is independent of Step-29's regular grid.
    valid=(r>=overlap(1) & r<=overlap(2) & r>tol);
    rq=r(valid);
    pr=nan(numel(rq),2,numel(Usets));
    for field=1:numel(Usets)
        uf=reshape(Usets{field},2,[]).'*R;
        upLocal=nan(numel(rUp),2);
        loLocal=nan(numel(rLo),2);
        for comp=1:2
            upLocal(:,comp)=accumarray(gUp,uf(up,comp),[],@mean);
            loLocal(:,comp)=accumarray(gLo,uf(lo,comp),[],@mean);
        end
        if side==1
            jump=upLocal(valid,:) - interp1(rLo,loLocal,rq,'pchip');
        else
            jump=interp1(rUp,upLocal,rq,'pchip') - loLocal(valid,:);
        end
        fac=mu/(kappa+1)*sqrt(2*pi./rq);
        pr(:,:,field)=bsxfun(@times,[jump(:,2),jump(:,1)],fac);
    end
    nativeR{side}=rq;
    nativeProfiles{side}=pr;
end

% Regression is performed against ACTUAL native face positions. Its sample
% count represents real data points, not Step 29's 61 virtual query points.
nMax=2*numel(lower)*numel(upper)*numel(degrees);
data=nan(nMax,16);
k=0;
EDI=[NaN;NaN];
if nargin>=3 && isstruct(O27) && isfield(O27,'rules') && ...
        isfield(O27,'Kactual')
    idxQ=find(O27.rules==16,1);
    [~,idxR]=min(abs(O27.rOuterOverA0-0.65));
    if ~isempty(idxQ),EDI=O27.Kactual(:,idxQ,idxR);end
end

fprintf('\n============================================================\n');
fprintf('STEP 30: NATIVE CRACK-FACE COD WINDOW/ORDER AUDIT\n');
fprintf('============================================================\n');
fprintf('  a0=%.7g m; h_tip/a0=%.5g, unique face nodes (upper/lower)=%d/%d\n', ...
    a0,O25.hTip/a0,numel(rUp),numel(rLo));
fprintf('  independent native-face sampling, no resampled virtual 61-point grid\n');
fprintf('  EDI16 reference: KI=%.8e, signed KII=%+.8e\n',EDI(1),EDI(2));

for side=1:2
    rr=nativeR{side}/a0;
    pr=nativeProfiles{side};
    sideNames={'upper-native','lower-native'};
    for il=1:numel(lower)
        for ih=1:numel(upper)
            if lower(il)>=upper(ih),continue;end
            idx=find(rr>=lower(il) & rr<=upper(ih));
            for jd=1:numel(degrees)
                degree=degrees(jd);
                k=k+1;
                if numel(idx)<max(P.MinNativePoints,2*(degree+1))
                    fprintf('  SKIP %-12s r/a0=[%.2f,%.2f], degree=%d: only %d real nodes\n', ...
                        sideNames{side},lower(il),upper(ih),degree,numel(idx));
                    continue;
                end
                fit=nan(2,numel(Usets));
                for field=1:numel(Usets)
                    for mode=1:2
                        coeff=polyfit(rr(idx),pr(idx,mode,field),degree);
                        fit(mode,field)=coeff(end);
                    end
                end
                M=[fit(:,2),fit(:,3)];
                ka=fit(:,1);
                data(k,:)=[side,lower(il),upper(ih),degree,numel(idx), ...
                    min(rr(idx)),max(rr(idx)),ka(1),ka(2), ...
                    ka(2)/ka(1),norm(M-eye(2),'fro'), ...
                    M(2,1),M(1,2),M(1,1),M(2,2), ...
                    max(diff(rr(idx)))];
                fprintf(['  %-12s [%4.2f,%4.2f] degree=%d n=%2d | ', ...
                    'KI=%+.8e KII=%+.8e ratio=%+.5e | exact ||M-I||=%.2e\n'], ...
                    sideNames{side},lower(il),upper(ih),degree,numel(idx), ...
                    ka(1),ka(2),ka(2)/ka(1),norm(M-eye(2),'fro'));
            end
        end
    end
end
validRows=isfinite(data(:,5));
T=array2table(data(validRows,:), ...
    'VariableNames',{'side_1upper_2lower','lower_r_over_a0', ...
    'upper_r_over_a0','polynomial_degree','n_actual_face_samples', ...
    'actual_min_r_over_a0','actual_max_r_over_a0', ...
    'KI_COD','KII_COD','KII_over_KI_COD', ...
    'synthetic_matrix_error','synthetic_KII_leak_from_I', ...
    'synthetic_KI_leak_from_II','synthetic_I_recovery', ...
    'synthetic_II_recovery','max_native_spacing_over_a0'});
fprintf('\nNATIVE-FACE COD FIT MATRIX\n');disp(T);

% Pair the two independently sampled faces when windows/degrees coincide.
% No automatic selection of a 'best' fit; only display the disagreement.
rPaired=nan(0,8);
for il=1:numel(lower)
    for ih=1:numel(upper)
        for jd=1:numel(degrees)
            if lower(il)>=upper(ih),continue;end
            A=find(T.side_1upper_2lower==1 & ...
                   abs(T.lower_r_over_a0-lower(il))<1e-12 & ...
                   abs(T.upper_r_over_a0-upper(ih))<1e-12 & ...
                   T.polynomial_degree==degrees(jd),1);
            B=find(T.side_1upper_2lower==2 & ...
                   abs(T.lower_r_over_a0-lower(il))<1e-12 & ...
                   abs(T.upper_r_over_a0-upper(ih))<1e-12 & ...
                   T.polynomial_degree==degrees(jd),1);
            if isempty(A)||isempty(B),continue;end
            rPaired(end+1,:)=[lower(il),upper(ih),degrees(jd), ...
                T.KII_over_KI_COD(A),T.KII_over_KI_COD(B), ...
                T.KII_over_KI_COD(A)-T.KII_over_KI_COD(B), ...
                T.KI_COD(A)-T.KI_COD(B), ...
                max(T.synthetic_matrix_error([A B]))]; %#ok<AGROW>
        end
    end
end
Pairs=array2table(rPaired,'VariableNames',{ ...
    'lower_r_over_a0','upper_r_over_a0','degree', ...
    'ratio_upper_native','ratio_lower_native','ratio_difference', ...
    'KI_difference','max_synthetic_matrix_error'});
fprintf('\nMATCHED FACE-SAMPLING CROSS-CHECKS\n');disp(Pairs);

if nargin>=4 && ~isempty(O29) && isstruct(O29) && ...
        isfield(O29,'table')
    fprintf('\nStep-29 regular-grid COD baseline for reference:\n');
    disp(O29.table(:,{'fit_upper_r_over_a0','KI_COD_raw','KII_over_KI_COD_raw'}));
end

if logical(P.Plot)
    figure('Name','Step 30: native face COD profiles','Color','w');
    t=tiledlayout(1,2,'Padding','compact','TileSpacing','compact');
    labels={'KI apparent','KII apparent'};
    for mode=1:2
        nexttile;hold on;grid on;box on
        for side=1:2
            plot(nativeR{side}/a0,nativeProfiles{side}(:,mode,1),'-o', ...
                'MarkerSize',3,'DisplayName',sideNames{side});
        end
        xlabel('r/a_0');ylabel(labels{mode});
        legend('Location','best');
    end
    title(t,'Actual FEM COD on native upper/lower face abscissae');
end

O30=struct('settings',P,'table',T,'paired',Pairs, ...
    'nativeR', {nativeR},'nativeProfiles',{nativeProfiles}, ...
    'EDI16_ref',EDI,'a0',a0,'hTip',O25.hTip);
fprintf(['STEP 30 completed. No new FEM solves; synthetic accuracy and ', ...
    'window/face-sampling stability must be considered together.\n']);
end

function must(S,f,who)
if ~isstruct(S)||~isfield(S,f)||isempty(S.(f))
    error('step30:MissingInput','Missing %s.%s.',who,f);
end
end
