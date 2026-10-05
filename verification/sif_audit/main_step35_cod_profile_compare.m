function R=main_step35_cod_profile_compare(O33,O34,varargin)
%MAIN_STEP35_COD_PROFILE_COMPARE
% Small-memory POSTPROCESSING ONLY of the two ALREADY SOLVED, genuinely
% nested Step-33/34 meshes. No FEM solves, new EDI integrations, T3 mesh
% copies, or interpolation-generated data points in any regression.
%
% Displays physical native crack-face COD apparent SIF profiles for both
% meshes, the existing 16-point EDI references, and six native-node COD
% extrapolation comparisons. Unlike a new fit, this only VISUALIZES native
% data previously recovered in the two completed steps.
%
% Examples:
% R=main_step35_cod_profile_compare(O33,O34);
% R=main_step35_cod_profile_compare(O33,O34,'Save',false);
% R=main_step35_cod_profile_compare(O33,O34,'Plot',false);
%
% If saving, produces a compact MAT file WITHOUT either full FEM mesh,
% plus a PNG and an editable MATLAB figure.

p=inputParser;
addParameter(p,'XLimits',[0.035 0.35], ...
    @(x)isnumeric(x)&&numel(x)==2&&all(isfinite(x))&&x(1)>=0&&x(2)>x(1));
addParameter(p,'Plot',true,@(x)islogical(x)&&isscalar(x));
addParameter(p,'Save',true,@(x)islogical(x)&&isscalar(x));
addParameter(p,'SavePrefix',fullfile(pwd,'step35_cod_profiles'), ...
    @(x)ischar(x)||(isstring(x)&&isscalar(x)));
parse(p,varargin{:});
opt=p.Results;

fields={'nativeR','nativeApparent','crack','newKI','newKII', ...
        'rOuterOverA0','CODTable'};
for j=1:numel(fields)
    required(O33,fields{j},'O33');
    required(O34,fields{j},'O34');
end
if norm(O33.crack.Pmid-O34.crack.Pmid,'fro')>1e-11
    error('step35:DifferentCracks', ...
        'Compare only solutions using exactly the same crack geometry.');
end
a0=norm(diff(O33.crack.Pmid,1,1));
x33=O33.nativeR(:)/a0;
x34=O34.nativeR(:)/a0;
p33=O33.nativeApparent;
p34=O34.nativeApparent;
if size(p33,1)~=numel(x33)||size(p34,1)~=numel(x34) || ...
        size(p33,2)~=2 || size(p34,2)~=2 || ...
        any(~isfinite(p33(:))) || any(~isfinite(p34(:)))
    error('step35:CODProfiles', ...
        'Native COD profiles must have N-by-2 finite apparent SIFs.');
end
if any(x33<=0)||any(x34<=0)|| ...
        any(diff(x33)<=0)||any(diff(x34)<=0)
    error('step35:CODAbscissae', ...
        'Native crack-face distances must be positive and increasing.');
end
rat=O33.rOuterOverA0(:).';
if numel(rat)~=numel(O34.rOuterOverA0)|| ...
        any(abs(rat-O34.rOuterOverA0(:).')>1e-12)
    error('step35:DifferentEDIDomains','EDI radii do not match.');
end
[~,jRef]=min(abs(rat-0.65));
q33=O33.newKII(jRef)/O33.newKI(jRef);
q34=O34.newKII(jRef)/O34.newKI(jRef);
EDI=table([33;34], ...
    [O33.newKI(jRef);O34.newKI(jRef)], ...
    [O33.newKII(jRef);O34.newKII(jRef)], ...
    [q33;q34], ...
    'VariableNames',{'step','KI_EDI16','KII_EDI16','KII_over_KI_EDI16'});

A=O33.CODTable;
B=O34.CODTable;
fitRows=nan(height(B),12);
for j=1:height(B)
    lower=B.lower_r_over_a0(j);
    upper=B.upper_r_over_a0(j);
    deg=B.degree(j);
    ix=find(abs(A.lower_r_over_a0-lower)<1e-12 & ...
        abs(A.upper_r_over_a0-upper)<1e-12 & ...
        A.degree==deg,1);
    if isempty(ix)
        error('step35:MissingPairedCOD', ...
            'The two solved meshes have incompatible COD fit windows.');
    end
    fitRows(j,:)=[lower upper deg ...
        A.n_native_face(ix) B.n_native_face(j) ...
        A.refined_KI(ix) B.refined_KI(j) ...
        A.refined_KII(ix) B.refined_KII(j) ...
        A.refined_ratio(ix) B.refined_ratio(j) ...
        100*(q34-B.refined_ratio(j))/q34];
end
Fits=array2table(fitRows,'VariableNames',{ ...
    'lower_r_over_a0','upper_r_over_a0','degree', ...
    'n_native_step33','n_native_step34', ...
    'KI_step33','KI_step34','KII_step33','KII_step34', ...
    'ratio_step33','ratio_step34','EDI34_vs_COD34_gap_pct_of_EDI'});

fprintf('\n============================================================\n');
fprintf('STEP 35: EXISTING-FIELD COD PROFILES AND EDI COMPARISON\n');
fprintf('============================================================\n');
fprintf('  a0=%.8g; original COD samples=%d, nested=%d\n', ...
    a0,numel(x33),numel(x34));
fprintf('  No new FEM, EDI integration, or fitted virtual sample points.\n');
fprintf('\nFIXED REFERENCE EDI\n');disp(EDI);
fprintf('\nNATIVE COD FITS COMPARED (different real node counts)\n');
disp(Fits);
fprintf(['  A consistent ~20%% EDI/COD gap despite nested mesh ', ...
    'convergence indicates residual method dependence; the plot ', ...
    'tests whether a near-tip COD plateau is visually supported.\n']);

fig=[];
if logical(opt.Plot)
    fig=figure('Name','Step 35: already solved native COD convergence', ...
        'NumberTitle','off','Color','w', ...
        'Position',[95 100 1260 580]);
    tl=tiledlayout(fig,1,3, ...
        'TileSpacing','compact','Padding','compact');
    names={'Apparent K_I','Apparent signed K_{II}', ...
        'Apparent signed K_{II}/K_I'};
    for j=1:3
        ax=nexttile(tl);
        hold(ax,'on');
        if j==1
            yy33=p33(:,1);
            yy34=p34(:,1);
            ref33=O33.newKI(jRef);
            ref34=O34.newKI(jRef);
        elseif j==2
            yy33=p33(:,2);
            yy34=p34(:,2);
            ref33=O33.newKII(jRef);
            ref34=O34.newKII(jRef);
        else
            yy33=p33(:,2)./p33(:,1);
            yy34=p34(:,2)./p34(:,1);
            ref33=q33;
            ref34=q34;
        end
        plot(ax,x33,yy33,'-o','MarkerSize',3, ...
            'LineWidth',1.0,'DisplayName','Step 33 native COD');
        plot(ax,x34,yy34,'-s','MarkerSize',3, ...
            'LineWidth',1.0,'DisplayName','Step 34 native COD');
        yline(ax,ref33,':','LineWidth',1.1, ...
            'DisplayName','Step 33 EDI');
        yline(ax,ref34,'--','LineWidth',1.1, ...
            'DisplayName','Step 34 EDI');
        xlim(ax,opt.XLimits);
        xlabel(ax,'Distance behind crack tip, r/a_0');
        ylabel(ax,names{j});
        grid(ax,'on');box(ax,'on');
        if j==3
            legend(ax,'Location','best','FontSize',8);
        end
    end
    title(tl,'Independent native crack-face COD vs 16-point EDI');
end

R=struct();
R.EDI=EDI;
R.CODFits=Fits;
R.xNative33=x33;
R.xNative34=x34;
R.apparent33=p33;
R.apparent34=p34;
R.ratioEDI33=q33;
R.ratioEDI34=q34;
R.fig=fig;
R.paths=struct('png','','fig','','mat','');
if logical(opt.Save)
    prefix=char(opt.SavePrefix);
    [dirName,~,~]=fileparts(prefix);
    if ~isempty(dirName)&&exist(dirName,'dir')~=7
        mkdir(dirName);
    end
    R.paths.mat=[prefix '_small_data.mat'];
    % No coordinates/connectivity or FEM field in the output file.
    EDIResult=EDI; %#ok<NASGU>
    CODFits=Fits; %#ok<NASGU>
    xNative33=x33; xNative34=x34; %#ok<NASGU>
    apparent33=p33; apparent34=p34; %#ok<NASGU>
    save(R.paths.mat,'EDIResult','CODFits', ...
        'xNative33','xNative34','apparent33','apparent34');
    fprintf('  Saved compact COD/EDI comparison: %s\n',R.paths.mat);
    if logical(opt.Plot)
        R.paths.png=[prefix '.png'];
        R.paths.fig=[prefix '.fig'];
        exportgraphics(fig,R.paths.png, ...
            'Resolution',230,'BackgroundColor','white');
        savefig(fig,R.paths.fig);
        fprintf('  Plots: %s and %s\n',R.paths.png,R.paths.fig);
    end
end
fprintf('STEP 35 completed without a new FEM solve.\n');
end

function required(S,f,name)
if ~isstruct(S)||~isfield(S,f)||isempty(S.(f))
    error('step35:MissingField','Missing %s.%s.',name,f);
end
end
