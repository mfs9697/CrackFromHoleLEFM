function O36=main_step36_cod_tip_window_audit(O33,O34,varargin)
%MAIN_STEP36_COD_TIP_WINDOW_AUDIT
% Diagnose whether near-tip conventional T6 displacement-jump profiles
% bias polynomial COD extrapolation on two ALREADY SOLVED nested meshes.
%
% NO new mesh, FEM solve, EDI integration or loading of large FE arrays.
% Only the stored native face r and apparent KI/KII vectors are accessed.
%
% Uses the SAME physical windows and polynomial orders for both meshes.
% Intercepts at r=0 are extrapolations, NOT directly measured values.
% Regression RMSE only describes fit quality; do not use as a physical
% error bar or as proof that either EDI or COD is exact.
%
% Usage:
%   O36=main_step36_cod_tip_window_audit(O33,O34);
%
% Outputs:
% - pointwise apparent SIF ratio by physical distance band
% - complete cutoff/order COD comparison and percent difference from EDI
% - charts exposing near-tip behavior and cutoff sensitivity
% - small numeric MAT and PNG; NO FE meshes or U fields are duplicated.

ip=inputParser;
addParameter(ip,'LowerCutoffs',[.02 .04 .06 .08 .10 .12], ...
    @(x)isnumeric(x)&&isvector(x)&&all(isfinite(x))&&all(x>0)&&all(x<.3));
addParameter(ip,'UpperCutoffs',[.16 .20 .30], ...
    @(x)isnumeric(x)&&isvector(x)&&all(isfinite(x))&&all(x>.10)&&all(x<.65));
addParameter(ip,'Degrees',[1 2], ...
    @(x)isnumeric(x)&&isvector(x)&&all(ismember(x,[1 2])));
addParameter(ip,'MinNodes',10, ...
    @(x)isnumeric(x)&&isscalar(x)&&isfinite(x)&&x>=8&&x==round(x));
addParameter(ip,'Plot',true,@(x)islogical(x)&&isscalar(x));
addParameter(ip,'Save',true,@(x)islogical(x)&&isscalar(x));
addParameter(ip,'Prefix',fullfile(pwd,'step36_cod_tip_cutoff'), ...
    @(x)ischar(x)||(isstring(x)&&isscalar(x)));
parse(ip,varargin{:});
opt=ip.Results;

for field={'nativeR','nativeApparent','crack','newKI','newKII', ...
        'rOuterOverA0','hTipNew'}
    required(O33,field{1},'O33');
    required(O34,field{1},'O34');
end
if norm(O33.crack.Pmid-O34.crack.Pmid,'fro')>1e-11 || ...
        abs(O33.hTipNew/O34.hTipNew-1)>1e-9
    error('step36:GeometryMismatch', ...
        'Nested levels must have the same crack and measured tip scale.');
end
a0=norm(diff(O34.crack.Pmid,1,1));
[~,ref33]=min(abs(O33.rOuterOverA0(:)-.65));
[~,ref34]=min(abs(O34.rOuterOverA0(:)-.65));
qEDI=[O33.newKII(ref33)/O33.newKI(ref33), ...
      O34.newKII(ref34)/O34.newKI(ref34)];
data={O33,O34};
lows=sort(unique(opt.LowerCutoffs(:).'));
highs=sort(unique(opt.UpperCutoffs(:).'));
degrees=sort(unique(opt.Degrees(:).'));
Nmax=2*numel(lows)*numel(highs)*numel(degrees);
rows=nan(Nmax,13);
k=0;
for iMesh=1:2
    Q=data{iMesh};
    x=Q.nativeR(:)/a0;
    y=Q.nativeApparent;
    if size(y,1)~=numel(x) || size(y,2)~=2 || ...
            any(~isfinite(x)) || any(~isfinite(y(:))) || any(diff(x)<=0)
        error('step36:InvalidNativeData', ...
            'Stored native face COD data must be sorted and finite.');
    end
    for low=lows
        for high=highs
            if low>=high,continue;end
            idx=find(x>=low & x<=high);
            for deg=degrees
                if numel(idx)<max(opt.MinNodes,4*(deg+1)),continue;end
                xx=x(idx);
                % Check that the basis is adequately conditioned. This
                % number does not quantify FEM or SIF extraction error.
                V=zeros(numel(xx),deg+1);
                for j=0:deg
                    V(:,deg+1-j)=xx.^j;
                end
                basisCond=cond(V);
                if ~isfinite(basisCond)||basisCond>1e6,continue;end
                pI=polyfit(xx,y(idx,1),deg);
                pII=polyfit(xx,y(idx,2),deg);
                ki=pI(end);
                kii=pII(end);
                rmse=sqrt(mean((polyval(pII,xx)-y(idx,2)).^2));
                k=k+1;
                rows(k,:)=[iMesh+32,low,high,deg,numel(idx), ...
                    min(xx)*a0/Q.hTipNew,ki,kii,kii/ki, ...
                    qEDI(iMesh),100*(qEDI(iMesh)-kii/ki)/qEDI(iMesh), ...
                    rmse,basisCond];
            end
        end
    end
end
rows=rows(1:k,:);
T=array2table(rows,'VariableNames',{ ...
    'step','lower_r_over_a0','upper_r_over_a0','polynomial_degree', ...
    'n_native_nodes','closest_fitted_node_r_over_htip', ...
    'KI_intercept','KII_intercept','COD_ratio','EDI_ratio', ...
    'EDI_minus_COD_pct_of_EDI','KII_fit_RMSE','polynomial_basis_condition'});
% Steps are 33 and 34, not arbitrary mesh indices.
fprintf('\n============================================================\n');
fprintf('STEP 36: COD NEAR-TIP CUTOFF SENSITIVITY; NO NEW FEM\n');
fprintf('============================================================\n');
fprintf('  tip edge=%.8g m, a0=%.8g m\n',O34.hTipNew,a0);
fprintf('  existing EDI ratios: Step33=%+.9e; Step34=%+.9e\n',qEDI);
fprintf('  COD fit RMSE is descriptive, NOT a physical uncertainty.\n');
fprintf('\nCOD INTERCEPT MATRIX (NATIVE MEASUREMENTS ONLY)\n');
disp(T);

% Describe RAW near-tip data before fitting: changes with r may reveal
% a non-Williams conventional-T6 inner zone or higher-order effects.
bands=[.0 .04;.04 .08;.08 .12;.12 .20;.20 .30];
B=nan(size(bands,1)*2,7);
row=0;
for iMesh=1:2
    Q=data{iMesh};
    x=Q.nativeR(:)/a0;
    Y=Q.nativeApparent;
    for ib=1:size(bands,1)
        row=row+1;
        idx=(x>=bands(ib,1)) & (x<bands(ib,2));
        if nnz(idx)<2,continue;end
        q=Y(idx,2)./Y(idx,1);
        B(row,:)=[iMesh+32,bands(ib,:),nnz(idx), ...
            median(q),min(q),max(q)];
    end
end
BandTable=array2table(B(isfinite(B(:,4)),:), ...
    'VariableNames',{'step','lower_r_over_a0','upper_r_over_a0', ...
    'n_native_nodes','median_pointwise_COD_ratio', ...
    'min_pointwise_COD_ratio','max_pointwise_COD_ratio'});
fprintf('\nRAW POINTWISE COD RATIOS BY PHYSICAL DISTANCE BAND\n');
disp(BandTable);

% There is deliberately no 'best window' selection. A credible COD
% intercept should be stable across a scientifically justified range,
% not cherry-picked to agree with EDI.
fig=[];
if opt.Plot
    fig=figure('Name','Step36 native COD cutoff stability','Color','w', ...
        'Position',[90 90 1250 570]);
    tl=tiledlayout(fig,1,2,'Padding','compact','TileSpacing','compact');
    ax=nexttile(tl);
    hold(ax,'on');grid(ax,'on');box(ax,'on');
    for iMesh=1:2
        Q=data{iMesh};
        x=Q.nativeR(:)/a0;
        y=Q.nativeApparent;
        plot(ax,x,y(:,2)./y(:,1),'-o','MarkerSize',3, ...
            'DisplayName',sprintf('Step %d native points',iMesh+32));
    end
    yline(ax,qEDI(1),':','DisplayName','Step 33 EDI');
    yline(ax,qEDI(2),'--','DisplayName','Step 34 EDI');
    xlim(ax,[0 .32]);
    xline(ax,O34.rInner/a0,':', ...
        'DisplayName','EDI inner radius');
    xlabel(ax,'r/a_0');ylabel(ax,'Pointwise apparent COD K_{II}/K_I');
    title(ax,'Native near-tip response (not extrapolated)');
    legend(ax,'Location','best','FontSize',8);

    ax=nexttile(tl);
    hold(ax,'on');grid(ax,'on');box(ax,'on');
    colors=lines(2);
    for iMesh=1:2
        for d=degrees
            sel=(T.step==iMesh+32) & ...
                T.polynomial_degree==d & ...
                abs(T.upper_r_over_a0-.30)<1e-12;
            t=T(sel,:);
            if isempty(t),continue;end
            if d==1,style='-o';else,style='--s';end
            plot(ax,t.lower_r_over_a0,t.COD_ratio,style, ...
                'MarkerSize',5,'Color',colors(iMesh,:), ...
                'DisplayName',sprintf('Step %d, degree %d',iMesh+32,d));
        end
    end
    yline(ax,qEDI(2),':','DisplayName','Step 34 EDI');
    xlabel(ax,'Excluded near-tip cutoff, r_{min}/a_0');
    ylabel(ax,'COD intercept K_{II}/K_I');
    title(ax,'Same upper bound r_{max}/a_0=0.30');
    legend(ax,'Location','best','FontSize',8);
    title(tl,'COD method dependence after two nested mesh levels');
end

O36=struct('settings',opt,'table',T,'bands',BandTable, ...
    'EDI_ratios',qEDI,'figure',fig,'files', ...
    struct('png','','fig','','mat',''));
if opt.Save
    prefix=char(opt.Prefix);
    [folder,~,~]=fileparts(prefix);
    if ~isempty(folder)&&exist(folder,'dir')~=7,mkdir(folder);end
    O36.files.mat=[prefix '_small_data.mat'];
    FitTable=T; %#ok<NASGU>
    RawBands=BandTable; %#ok<NASGU>
    EDIRatios=qEDI; %#ok<NASGU>
    save(O36.files.mat,'FitTable','RawBands','EDIRatios');
    if opt.Plot
        O36.files.png=[prefix '.png'];
        O36.files.fig=[prefix '.fig'];
        exportgraphics(fig,O36.files.png, ...
            'Resolution',220,'BackgroundColor','white');
        savefig(fig,O36.files.fig);
    end
    fprintf('\n  Compact numeric results: %s\n',O36.files.mat);
    if opt.Plot
        fprintf('  Figure: %s\n',O36.files.png);
    end
end
fprintf('STEP 36 completed: existing numerical data only.\n');
end

function required(S,field,name)
if ~isstruct(S)||~isfield(S,field)||isempty(S.(field))
    error('step36:MissingField','Missing %s.%s.',name,field);
end
end
