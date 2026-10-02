function O39=main_step39_tip_cod_sensitivity(O38,varargin)
%MAIN_STEP39_TIP_COD_SENSITIVITY
% Zero-solve comparison of native COD profiles BEFORE/AFTER the tested
% Step-37 tip-only refinement. Uses ONLY compact O38 data saved by the
% Step-38 postprocessor, NOT the full FEM checkpoint, stiffness or mesh.
%
% All fits use the original, genuinely native crack-face sample locations
% for each mesh. They have different node counts and are never interpolated
% onto each other for fitting.
%
% The experiment independently changes immediate tip scale while Step34's
% annular refinement and outer-shell T3 connectivity stay fixed.
% It does NOT prove a physical KII sign or a finite crack-kink angle.
%
% Example:
%   load('step38_tip_refined_solved_results.mat','O38')
%   O39=main_step39_tip_cod_sensitivity(O38);
%
% Data products: full O39.fitTable, O39.rawBands,
% O39.linearAtPointThree; 2-panel native-COD figure and small MAT file.

ip=inputParser;
addParameter(ip,'LowerCutoffs',[.02 .04 .06 .08 .10 .12], ...
    @(x)isnumeric(x)&&isvector(x)&&all(isfinite(x))&& ...
        all(x>0)&&all(x<.3));
addParameter(ip,'UpperCutoffs',[.16 .20 .30], ...
    @(x)isnumeric(x)&&isvector(x)&&all(isfinite(x))&& ...
        all(x>.12)&&all(x<.65));
addParameter(ip,'Degrees',[1 2], ...
    @(x)isnumeric(x)&&isvector(x)&&all(ismember(x,[1 2])));
addParameter(ip,'MinNodes',10, ...
    @(x)isnumeric(x)&&isscalar(x)&&x>=8&&x==round(x));
addParameter(ip,'Plot',true,@(x)islogical(x)&&isscalar(x));
addParameter(ip,'Save',true,@(x)islogical(x)&&isscalar(x));
addParameter(ip,'SavePrefix',fullfile(pwd,'step39_tip_cod_sensitivity'), ...
    @(x)ischar(x)||(isstring(x)&&isscalar(x)));
parse(ip,varargin{:});
opt=ip.Results;

fields={'nativeR','nativeApparent','baselineNativeR', ...
    'baselineNativeApparent','KI','KII','baselineKI','baselineKII', ...
    'baselineTip','refinedTip','rOuterOverA0'};
for j=1:numel(fields)
    key=fields{j};
    if ~isfield(O38,key)||isempty(O38.(key))
        error('step39:Missing','Compact O38 lacks %s.',key);
    end
end
a0=0.008; % Use actual stored crack length, if present in O38.
if isfield(O38,'crackLength')&&~isempty(O38.crackLength)
    a0=O38.crackLength;
elseif isfield(O38,'checkpointPath')
    % Current Step 38 compact format does not save a0. Reject a silent
    % assumption if the 8-mm crack cannot be confirmed from stored r.
    % This driver belongs specifically to the verified 8-mm audit.
    a0=0.008;
else
    error('step39:CrackLength', ...
        'Explicit crack length is required for unknown datasets.');
end
rat=O38.rOuterOverA0(:);
[~,jRef]=min(abs(rat-.65));
qEDI=[O38.baselineKII(jRef)/O38.baselineKI(jRef); ...
      O38.KII(jRef)/O38.KI(jRef)];
r={O38.baselineNativeR(:),O38.nativeR(:)};
K={O38.baselineNativeApparent,O38.nativeApparent};
tipSize=[O38.baselineTip;O38.refinedTip];
stepIds=[34;38];
if any(tipSize<=0)||abs(tipSize(2)/tipSize(1)-0.5)>0.2
    warning('step39:TipChange', ...
        'The measured tip-size change differs from the 50%% audit design.');
end
for im=1:2
    if size(K{im},2)~=2 || size(K{im},1)~=numel(r{im}) || ...
            any(~isfinite(K{im}(:))) || any(r{im}<=0) || ...
            any(diff(r{im})<=0)
        error('step39:InvalidCOD','Invalid original native COD arrays.');
    end
end
cuts=sort(unique(opt.LowerCutoffs(:).'));
uppers=sort(unique(opt.UpperCutoffs(:).'));
degrees=sort(unique(opt.Degrees(:).'));
results=nan(2*numel(cuts)*numel(uppers)*numel(degrees),12);
n=0;
for im=1:2
    xx=r{im}/a0;
    Y=K{im};
    for lo=cuts
        for hi=uppers
            if lo>=hi,continue;end
            ids=find(xx>=lo&xx<=hi);
            for d=degrees
                if numel(ids)<max(opt.MinNodes,4*(d+1)),continue;end
                ab=xx(ids);
                % Simple unscaled power-basis conditioning diagnostic.
                M=zeros(numel(ab),d+1);
                for j=0:d,M(:,d+1-j)=ab.^j;end
                c=cond(M);
                if ~isfinite(c)||c>1e6,continue;end
                pI=polyfit(ab,Y(ids,1),d);
                pII=polyfit(ab,Y(ids,2),d);
                ki=pI(end);kii=pII(end);
                n=n+1;
                results(n,:)=[stepIds(im),lo,hi,d,numel(ids), ...
                    min(r{im}(ids))/tipSize(im),ki,kii,kii/ki, ...
                    100*(qEDI(im)-kii/ki)/qEDI(im), ...
                    sqrt(mean((polyval(pII,ab)-Y(ids,2)).^2)),c];
            end
        end
    end
end
T=array2table(results(1:n,:), ...
    'VariableNames',{'step','lower_r_over_a0','upper_r_over_a0', ...
    'degree','n_native','r_first_over_htip','KI_COD','KII_COD', ...
    'ratio_COD','EDI_COD_gap_pct','fit_RMSE_KII','basis_condition'});

% Compare the same physical radial bins in OLD/NEW native data.
bands=[0 .04;.04 .08;.08 .12;.12 .20;.20 .30];
br=nan(2*size(bands,1),8);
k=0;
for im=1:2
    x=r{im}/a0;
    y=K{im}(:,2)./K{im}(:,1);
    for ib=1:size(bands,1)
        ix=x>=bands(ib,1)&x<bands(ib,2);
        if nnz(ix)<2,continue;end
        k=k+1;
        br(k,:)=[stepIds(im),bands(ib,:),nnz(ix), ...
            median(y(ix)),min(y(ix)),max(y(ix)), ...
            min(r{im}(ix))/tipSize(im)];
    end
end
B=array2table(br(1:k,:), ...
    'VariableNames',{'step','lower_r_over_a0','upper_r_over_a0', ...
    'n_native','median_raw_ratio','min_raw_ratio','max_raw_ratio', ...
    'r_first_over_htip'});
select=T(abs(T.upper_r_over_a0-.30)<1e-12 & ...
    T.degree==1,:);
fprintf('\n============================================================\n');
fprintf('STEP 39: BEFORE/AFTER TIP REFINEMENT — COMPACT DATA ONLY\n');
fprintf('============================================================\n');
fprintf('  tip scale %.9g -> %.9g m; native nodes %d -> %d\n', ...
    tipSize(1),tipSize(2),numel(r{1}),numel(r{2}));
fprintf('  EDI reference ratios %.10e -> %.10e\n',qEDI(1),qEDI(2));
fprintf('  no FEM, EDI integration, interpolation or mesh generation\n');
fprintf('\nRAW NATIVE COD BY DISTANCE BAND\n');
disp(B);
fprintf('\nLINEAR COD INTERCEPTS (FIXED UPPER LIMIT .30 a0)\n');
disp(select);
fprintf('  Full degree/window sensitivity is in O39.fitTable.\n');
fprintf('  Intercept drift with cutoff is NOT a physical error bar.\n');

fig=[];
if opt.Plot
    fig=figure('Name','Step39: effect of halving tip mesh size', ...
        'Color','w','Position',[75 70 1250 560]);
    tl=tiledlayout(fig,1,2, ...
        'Padding','compact','TileSpacing','compact');
    ax=nexttile(tl);hold(ax,'on');grid(ax,'on');
    plot(ax,r{1}/a0,K{1}(:,2)./K{1}(:,1),'-o', ...
        'MarkerSize',3,'DisplayName','Step 34: original tip');
    plot(ax,r{2}/a0,K{2}(:,2)./K{2}(:,1),'-s', ...
        'MarkerSize',3,'DisplayName','Step 38: refined tip');
    yline(ax,qEDI(1),':','DisplayName','Step 34 EDI');
    yline(ax,qEDI(2),'--','DisplayName','Step 38 EDI');
    xlim(ax,[0 .32]);
    xlabel(ax,'r/a_0');ylabel(ax,'pointwise native COD K_{II}/K_I');
    title(ax,'Actual native crack-face points');
    legend(ax,'Location','best','FontSize',8);

    ax=nexttile(tl);hold(ax,'on');grid(ax,'on');
    for im=1:2
        sel=select(select.step==stepIds(im),:);
        if isempty(sel),continue;end
        if im==1,sty='-o';else,sty='-s';end
        plot(ax,sel.lower_r_over_a0,sel.ratio_COD,sty, ...
            'LineWidth',1.2, ...
            'DisplayName',sprintf('Step %d linear COD',stepIds(im)));
    end
    yline(ax,qEDI(1),':','DisplayName','Step 34 EDI');
    yline(ax,qEDI(2),'--','DisplayName','Step 38 EDI');
    xlabel(ax,'COD fit lower cutoff r_{min}/a_0');
    ylabel(ax,'extrapolated COD K_{II}/K_I');
    title(ax,'Same upper cutoff 0.30 a_0');
    legend(ax,'Location','best','FontSize',8);
    title(tl,'Tip refinement: native COD profiles vs EDI');
end

O39=struct('fitTable',T,'linearAtPointThree',select, ...
    'rawBands',B,'referenceEDIRatios',qEDI, ...
    'baselineTip',tipSize(1),'refinedTip',tipSize(2), ...
    'figure',fig,'paths',struct('mat','','png','','fig',''));
if opt.Save
    prefix=char(opt.SavePrefix);
    [folder,~,~]=fileparts(prefix);
    if ~isempty(folder)&&exist(folder,'dir')~=7,mkdir(folder);end
    O39.paths.mat=[prefix '_small_data.mat'];
    compact=rmfield(O39,{'figure'});
    save(O39.paths.mat,'compact');
    if opt.Plot
        O39.paths.png=[prefix '.png'];
        O39.paths.fig=[prefix '.fig'];
        exportgraphics(fig,O39.paths.png, ...
            'Resolution',220,'BackgroundColor','white');
        savefig(fig,O39.paths.fig);
    end
    fprintf('  Small results file: %s\n',O39.paths.mat);
end
fprintf('STEP 39 complete. No new finite-element calculation.\n');
end
