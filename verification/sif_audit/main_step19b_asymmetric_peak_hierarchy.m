function Out=main_step19b_asymmetric_peak_hierarchy(O19,varargin)
%MAIN_STEP19B_ASYMMETRIC_PEAK_HIERARCHY
% Re-rank physically independent Stage-I peaks using an existing Step-19 run.
%
% Step 19 showed that the two highest "local peaks" can be duplicate
% quadratic refinements of the same broad physical maximum: their fitted
% angles differ by only O(1e-2 deg), far below the mesh angular scale.
%
% This postprocessor performs non-maximum suppression in angle after local
% refinement.  Two fitted maxima are considered the same physical peak if
% their angular separation is less than
%
%     clusterFactor * h_hole/R.
%
% No FEM solve is repeated; O19.runs{...}.B is reused directly.

ip=inputParser;
addParameter(ip,'WindowFactor',3.0,@(x)isnumeric(x)&&isscalar(x)&&x>0);
addParameter(ip,'ClusterFactor',3.0,@(x)isnumeric(x)&&isscalar(x)&&x>0);
addParameter(ip,'MaxPeaks',6,@(x)isnumeric(x)&&isscalar(x)&&x>=2);
addParameter(ip,'Plot',true,@(x)islogical(x)||isnumeric(x));
parse(ip,varargin{:});
O=ip.Results;

if nargin<1 || ~isstruct(O19) || ~isfield(O19,'runs')
    error('step19b:NeedO19', ...
        'Pass the existing O19 returned by main_step19_asymmetric_stage1_validation.');
end

nL=numel(O19.runs);
Tables=cell(nL,1);

fprintf('\n============================================================\n');
fprintf('SIF AUDIT STEP 19B: INDEPENDENT ASYMMETRIC PEAK HIERARCHY\n');
fprintf('============================================================\n');

for il=1:nL
    R=O19.runs{il};
    B=R.B;
    C=R.C;

    T=local_rank_independent_peaks(B,O.WindowFactor,O.ClusterFactor,O.MaxPeaks);
    Tables{il}=T;

    meshAngle=B.offset.hhole/B.hole.r;

    fprintf('\n--- Npoly=%d | h/R=%.8g rad = %.6f deg ---\n', ...
        C.hole.npoly,meshAngle,rad2deg(meshAngle));
    fprintf('  regression half-width = %.6f deg\n', ...
        rad2deg(O.WindowFactor*meshAngle));
    fprintf('  clustering separation = %.6f deg\n', ...
        rad2deg(O.ClusterFactor*meshAngle));

    disp(T);

    if height(T)>=2
        fprintf(['  dominant independent peak: phi=%+.8f deg, sigma=%.10e\n', ...
                 '  second independent peak  : phi=%+.8f deg, sigma=%.10e\n', ...
                 '  relative stress gap      : %.6e\n', ...
                 '  angular separation       : %.6f deg\n'], ...
            T.phi_fit_deg(1),T.sig_fit(1), ...
            T.phi_fit_deg(2),T.sig_fit(2), ...
            T.top_minus_this_rel(2),T.separation_from_top_deg(2));
    end
end

% Cross-mesh comparison of the first two independent peaks.
rows=nan(nL,8);
for il=1:nL
    T=Tables{il};
    R=O19.runs{il};

    if height(T)>=2
        rows(il,:)=[ ...
            R.C.hole.npoly, ...
            T.phi_fit_deg(1),T.sig_fit(1), ...
            T.phi_fit_deg(2),T.sig_fit(2), ...
            T.top_minus_this_rel(2), ...
            T.separation_from_top_deg(2), ...
            T.sig_fit(2)/T.sig_fit(1)];
    else
        rows(il,:)=[R.C.hole.npoly,nan(1,7)];
    end
end

S=array2table(rows,'VariableNames',{ ...
    'Npoly','phi1_deg','sig1','phi2_deg','sig2', ...
    'top_second_gap_rel','peak_separation_deg','sig2_over_sig1'});

fprintf('\nINDEPENDENT-PEAK CROSS-MESH SUMMARY\n');
disp(S);

if logical(O.Plot)
    figure('Name','Step 19B: independent peak hierarchy','Color','w');
    clf; hold on; box on; grid on

    for il=1:nL
        T=Tables{il};
        if isempty(T), continue; end
        k=(1:height(T)).';
        plot(k,T.sig_fit/T.sig_fit(1),'-o','LineWidth',1.1, ...
            'DisplayName',sprintf('Npoly=%d',O19.runs{il}.C.hole.npoly));
    end

    xlabel('independent peak rank');
    ylabel('sigma_{peak}/sigma_{peak,1}');
    title('Physically separated Stage-I tensile maxima');
    legend('Location','best');
end

Out=struct();
Out.settings=O;
Out.tables=Tables;
Out.summary=S;

fprintf('\nSTEP 19B completed without any new FEM solves.\n');
end


function T=local_rank_independent_peaks(B,windowFactor,clusterFactor,nKeep)
phi=B.phi(:);
sig=max(B.sig_tt_eff(:),0);
n=numel(sig);

prev=sig(1+mod((0:n-1)-1,n));
next=sig(1+mod((0:n-1)+1,n));
idx=find(sig>=prev(:) & sig>=next(:) & sig>0);

meshAngle=B.offset.hhole/B.hole.r;
halfWidth=windowFactor*meshAngle;
clusterSep=clusterFactor*meshAngle;

cand=[];
for k=1:numel(idx)
    F=local_refine(phi,sig,idx(k),halfWidth);
    if F.accepted
        ph=F.phiStar;
        sf=F.sigStar;
    else
        ph=phi(idx(k));
        sf=sig(idx(k));
    end

    cand(end+1,:)=[ ...
        idx(k), ...
        rad2deg(local_wrap(phi(idx(k)))), ...
        rad2deg(local_wrap(ph)), ...
        sf,double(F.accepted),F.rmse]; %#ok<AGROW>
end

if isempty(cand)
    T=array2table(zeros(0,9),'VariableNames',{ ...
        'idx','phi_discrete_deg','phi_fit_deg','sig_fit','fit_accepted', ...
        'fit_rmse','cluster_size','top_minus_this_rel','separation_from_top_deg'});
    return;
end

[~,ord]=sort(cand(:,4),'descend');
cand=cand(ord,:);

accepted=zeros(0,size(cand,2)+1);

for k=1:size(cand,1)
    ph=deg2rad(cand(k,3));

    if isempty(accepted)
        accepted=[cand(k,:),1]; %#ok<AGROW>
        continue;
    end

    phAcc=deg2rad(accepted(:,3));
    sep=abs(local_wrap(ph-phAcc));

    [dmin,j]=min(sep);

    if dmin<clusterSep
        accepted(j,end)=accepted(j,end)+1;
    else
        accepted=[accepted;cand(k,:),1]; %#ok<AGROW>
    end

    if size(accepted,1)>=nKeep
        % Continue only to accumulate duplicates into already accepted
        % clusters; no need to add further low-ranked clusters.
    end
end

% Sort clusters again by representative peak stress and retain nKeep.
[~,ord]=sort(accepted(:,4),'descend');
accepted=accepted(ord,:);
accepted=accepted(1:min(nKeep,size(accepted,1)),:);

topSig=accepted(1,4);
topPhi=deg2rad(accepted(1,3));

gap=(topSig-accepted(:,4))/max(abs(topSig),eps);
sep=abs(rad2deg(local_wrap(deg2rad(accepted(:,3))-topPhi)));

T=array2table([accepted,gap,sep],'VariableNames',{ ...
    'idx','phi_discrete_deg','phi_fit_deg','sig_fit','fit_accepted', ...
    'fit_rmse','cluster_size','top_minus_this_rel','separation_from_top_deg'});
end


function F=local_refine(phi,sig,idx0,halfWidth)
phi0=phi(idx0);
d=atan2(sin(phi-phi0),cos(phi-phi0));
keep=abs(d)<=halfWidth+100*eps;

x=d(keep);
y=sig(keep);
[x,ord]=sort(x);
y=y(ord); %#ok<NASGU>

F=struct('accepted',false,'phiStar',phi0,'sigStar',sig(idx0), ...
    'rmse',NaN,'curvature',NaN,'vertexOffset',NaN);

if numel(x)<5
    return;
end

p=polyfit(x,y,2);
F.rmse=sqrt(mean((y-polyval(p,x)).^2));
F.curvature=p(1);

if any(~isfinite(p)) || p(1)>=0
    return;
end

dv=-p(2)/(2*p(1));
F.vertexOffset=dv;

if ~isfinite(dv) || abs(dv)>halfWidth
    return;
end

sf=polyval(p,dv);
if ~isfinite(sf) || sf<=0
    return;
end

F.accepted=true;
F.phiStar=mod(phi0+dv,2*pi);
F.sigStar=sf;
end


function a=local_wrap(a)
a=mod(a+pi,2*pi)-pi;
end
